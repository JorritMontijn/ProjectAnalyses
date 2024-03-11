%% aim
%{
Decode using sets of 20 spikes; is duration correlated with confidence of correct stim/accuracy?
I.e., are codes equally efficient during high and low rate periods?
%}

%% define qualifying areas
clear all;
runHeaderPopTimeCoding;
cellUseAreas = {...
	'Primary visual area',...
	...'posteromedial visual area',...
	};
strArea = cellUseAreas{1};

%% select all neurons in LP and drifting grating stimuli
if ~exist('sAggStim','var') || isempty(sAggStim) || isempty(sAggNeuron)
	[sAggStim,sAggNeuron]=loadDataNpx(strArea,'driftinggrating',strDataPath);
end
indRemDBA = strcmpi({sAggNeuron.SubjectType},'DBA');
fprintf('Removing %d cells of DBA animals; %d remaining [%s]\n',sum(indRemDBA),sum(~indRemDBA),getTime);
sAggNeuron(indRemDBA) = [];

%% pre-allocate matrices
strRunType = 'Npx';
strRunStim = 'DG';
cellTypes = {'Real','ShuffTid'};%, 'UniformTrial', 'ShuffTid'};%, 'PoissGain'}; %PoissGain not done
intNumTypes = numel(cellTypes);
boolFixSpikeGroupSize = true;

%% go through recordings
tic
for intRec=1:numel(sAggStim) %19 || weird: 11
	% get matching recording data
	strRec = sAggStim(intRec).Exp;
	sThisRec = sAggStim(strcmpi(strRec,{sAggStim(:).Exp}));
	
	%prep grating data
	[sUseNeuron,vecStimOnTime,vecStimOffTime,vecOrientation] = NpxPrepGrating(sAggNeuron,sThisRec,cellUseAreas);
	if isempty(sUseNeuron),continue;end
	[vecOriIdx,vecUniqueOris,vecRepNum,cellSelect,vecTrialRepetition] = val2idx(vecOrientation);
	intTrialNum = numel(vecStimOnTime);
	intStimNum = numel(unique(vecOrientation));
	intRepNum = intTrialNum/intStimNum;
	
	%% select area
	indArea1Neurons = contains({sUseNeuron.Area},strArea,'IgnoreCase',true);
	if sum(indArea1Neurons) == 0, continue;end
	
	%% get orientation responses & single-trial population noise
	sArea1Neurons = sUseNeuron(indArea1Neurons);
	
	%% prep data
	%get data matrix
	cellSpikeTimesRaw = {sArea1Neurons.SpikeTimes};
	[matData,indTuned,cellSpikeTimes,sOut,cellSpikeTimesPerCellPerTrial] = ...
		NpxPrepData(cellSpikeTimesRaw,vecStimOnTime,vecStimOffTime,vecOrientation);
	intTunedN = sum(indTuned);
	intNumN = size(matData,1);
	
	%get ori vars
	vecOri180 = mod(vecOrientation,180)*2;
	[vecOriIdx,vecUnique,vecPriorDistribution,cellSelect,vecRepetition] = val2idx(vecOri180);
	intOriNum = numel(vecUnique);
	intRepNum = min(vecPriorDistribution);
	dblStimDur = roundi(median(vecStimOffTime - vecStimOnTime),1,'ceil');
	dblBinWidth = 5/1000;%/32
	dblPreTime = 0.3;%10*dblBinWidth;
	dblPostTime = 0;%30*dblBinWidth;
	dblMaxDur = dblStimDur+dblPreTime+dblPostTime;
	if mean(sum(matData)) < 90
		fprintf('Avg # of spikes per trial was %.1f for %s; skipping...\n',mean(sum(matData)),strRec);
		continue;
	end
	if boolFixSpikeGroupSize
		intSpikeGroupSize = 10;
		strSGS = ['Fixed' num2str(intSpikeGroupSize)];
	else
		intSpikeGroupSize = ceil(mean(sum(matData))/30); %20
		strSGS = ['Var' num2str(intSpikeGroupSize)];
	end
	
	for intType=1%:numel(cellRunTypes)
		%% load prepro T0 data
		% pool spikes from all neurons, but save the time+id per spike, then calculate IFR over all
		% spikes at pop level, detect peaks, and split into trials.
		%get spikes per trial per neuron
		strType = cellTypes{intType};
		sSource = load(fullpath(strTargetDataPath,sprintf('T0Data_%s%s%s%s',strRec,strRunType,strRunStim,strType)));

		%% detect peaks
		% filter
		%real
		intLag = round((numel(sSource.vecIFR)/numel(vecStimOnTime))/2);
		if (intLag/2) == round(intLag/2)
			intLag = intLag - 1;
		end
		dblThreshZ = 1;
		dblInfluence = 0.5;
		
		%% transform time indices
		%events
		dblStartEpoch = vecStimOnTime(1)-dblStimDur;
		dblEpochDur = vecStimOnTime(end)-vecStimOnTime(1)+dblStimDur;
		dblStopEpoch = dblStartEpoch + dblEpochDur;
		
		%% plot
		vecTime = sSource.vecTime;
		vecIFR = sSource.vecIFR;
		vecAllSpikeTime = sSource.vecAllSpikeTime;
		vecAllSpikeNeuron = sSource.vecAllSpikeNeuron;
		intNumN = max(vecAllSpikeNeuron);
		
		%% replace IFR with binned rates?
		if 0
			dblBinDur = (5/1000);
			vecBins=(dblStartEpoch:dblBinDur:dblStopEpoch)';
			vecIFR = flat(histcounts(vecTime,vecBins)./dblBinDur);
			vecTime = vecBins(2:end)-dblBinDur/2;
		end
		
		%% run analyses
		%remove spikes outside epoch & build shuffled array
		cellUseSpikeTimes = cell(1,intNumN);
		for intN=1:intNumN
			vecT = vecAllSpikeTime(vecAllSpikeNeuron==intN);
			indRem = (vecT > dblStopEpoch) | (vecT < dblStartEpoch);
			cellUseSpikeTimes{intN} = vecT(~indRem);
		end
		%add data
		sAllSpike = struct;
		sAllSpike.vecAllSpikeTime = vecAllSpikeTime;
		sAllSpike.vecAllSpikeNeuron = vecAllSpikeNeuron;
		sAllSpike.dblStartEpoch = dblStartEpoch;
		sAllSpike.dblEpochDur = dblEpochDur;
		
		%set params and run pop event detection
		dblCutOff = 0.65;
		sPeakOpts = struct;
		sPeakOpts.intLag = intLag;
		sPeakOpts.dblThreshZ = dblThreshZ;
		sPeakOpts.dblInfluence = dblInfluence;
		[sPopEvents,sMergedPopEvents,vecStimTime,vecStimIFR,vecStimIFR_Raw] = getPopEvents(vecIFR,vecTime,vecStimOnTime,sPeakOpts,sAllSpike,dblCutOff);
		
		%extract
		vecPopEventTimes = [sPopEvents.Time];
		vecPopEventLocs = [sPopEvents.Loc];
		vecMergedPopEventTimes = [sMergedPopEvents.Time];
		vecMergedPopEventLocs = [sMergedPopEvents.Loc];
		vecPeakHeight = vecStimIFR(vecMergedPopEventLocs);
		
		%% are spikes in first half of peaks more tuned than spikes in second half of peak?
		dblAvgRate = sum(cellfun(@numel,cellUseSpikeTimes))/dblEpochDur;
		if boolFixSpikeGroupSize
			intSpikeGroupSize = 10;
			strSGS = ['Fixed' num2str(intSpikeGroupSize)];
		else
			intSpikeGroupSize = ceil(dblAvgRate/30); %20
			strSGS = ['Var' num2str(intSpikeGroupSize)];
		end

		%% randomize spikes (or not)
		if strcmp(strType,'ShuffTid')
			%shuffle across repetitions per neuron
			cellUseSpikeTimesPerCellPerTrial = cell(size(cellSpikeTimesPerCellPerTrial));
			for intStimType=1:intOriNum
				vecOrigTrials = find(vecOriIdx==intStimType);
				for intN=1:intNumN
					vecShuffTrials = vecOrigTrials(randperm(numel(vecOrigTrials)));
					cellUseSpikeTimesPerCellPerTrial(intN,vecShuffTrials) = cellSpikeTimesPerCellPerTrial(intN,vecOrigTrials);
				end
			end
		elseif strcmp(strType,'Real')
			%do nothing
			cellUseSpikeTimesPerCellPerTrial = cellSpikeTimesPerCellPerTrial;
		else
			error no type
		end
		
		%% go through pop spikes and group into sets of 20
		%throw away spikes not in stimulus period and assign all spikes
		intTotalSpikeNum = sum(sum(cellfun(@numel,cellUseSpikeTimesPerCellPerTrial)));
		matSpikeGroupData = nan(floor(intTotalSpikeNum/intSpikeGroupSize),intNumN);
		vecSpikeGroupDuration = nan(floor(intTotalSpikeNum/intSpikeGroupSize),1);
		vecSpikeGroupLatency = nan(floor(intTotalSpikeNum/intSpikeGroupSize),1);
		vecSpikeGroupTrialType = nan(floor(intTotalSpikeNum/intSpikeGroupSize),1);
		vecSpikeGroupTrialNumber = nan(floor(intTotalSpikeNum/intSpikeGroupSize),1);
		vecTrialType = nan(floor(intTotalSpikeNum/intSpikeGroupSize),1);
		intSpikeGroupCounter = 0;
		for intTrial=1:size(cellUseSpikeTimesPerCellPerTrial,2)
			cellSpikesInTrial = cellUseSpikeTimesPerCellPerTrial(:,intTrial);
			intSpikesInTrial = sum(sum(cellfun(@numel,cellSpikesInTrial)));
			vecSpikeTimes = nan(1,intSpikesInTrial);
			vecSpikeNeuron= nan(1,intSpikesInTrial);
			intSpikeC = 1;
			for intN=1:intNumN
				vecSpikeT = cellSpikesInTrial{intN};
				intNumS = numel(vecSpikeT);
				vecAssign = intSpikeC:(intSpikeC+intNumS-1);
				intSpikeC = intSpikeC + intNumS;
				
				vecSpikeTimes(vecAssign) = vecSpikeT;
				vecSpikeNeuron(vecAssign) = intN;
			end
			
			%sort spikes
			[vecSpikeTimes,vecSort]=sort(vecSpikeTimes);
			vecSpikeNeuron = vecSpikeNeuron(vecSort);
			
			%assign to spike group
			intAssignGroups = floor(numel(vecSpikeTimes)/intSpikeGroupSize);
			for intGroup=1:intAssignGroups
				%get data
				intEndSpike = intGroup*intSpikeGroupSize;
				vecUseSpikes = (intEndSpike-intSpikeGroupSize+1):intEndSpike;
				vecT = vecSpikeTimes(vecUseSpikes);
				vecN = vecSpikeNeuron(vecUseSpikes);
				vecCounts = accumarray(vecN(:),1,[intNumN 1]);
				
				%assign
				intAssignTo = intGroup + intSpikeGroupCounter;
				matSpikeGroupData(intAssignTo,:) = vecCounts;
				vecSpikeGroupDuration(intAssignTo) = range(vecT);
				vecSpikeGroupLatency(intAssignTo) = mean(vecT);
				vecSpikeGroupTrialType(intAssignTo) = vecOriIdx(intTrial);
				vecSpikeGroupTrialNumber(intAssignTo) = intTrial;
			end
			intSpikeGroupCounter = intSpikeGroupCounter + intAssignGroups;
		end
		
		%remove trailing entries
		intSpikeGroupNum = intSpikeGroupCounter;
		matSpikeGroupData = matSpikeGroupData(1:intSpikeGroupNum,:);
		vecSpikeGroupDuration = vecSpikeGroupDuration(1:intSpikeGroupNum);
		vecSpikeGroupLatency = vecSpikeGroupLatency(1:intSpikeGroupNum);
		vecSpikeGroupTrialType = vecSpikeGroupTrialType(1:intSpikeGroupNum);
		vecSpikeGroupTrialNumber = vecSpikeGroupTrialNumber(1:intSpikeGroupNum);
		
		%% decode all spike groups
		dblLambda = 1;
		intTypeCV = 0;
		
		%do old logistic regression
		vecPriorDistribution = accumarray(vecSpikeGroupTrialType(:),1,[intOriNum 1]);
		intVerbose = 1;
		[dblPerformanceCV,vecSpikeGroupDecodedTrial,matPosteriorProbability,dblMeanErrorDegs,matConfusion,matWeights] = ...
			doCrossValidatedDecodingLR(matSpikeGroupData,vecSpikeGroupTrialType,intTypeCV,[],dblLambda,intVerbose);
		
		vecSpikeGroupCorrect = vecSpikeGroupDecodedTrial(:) == vecSpikeGroupTrialType;
		vecSpikeGroupConfidence = nan(intSpikeGroupNum,1);
		for intStimType=1:intOriNum
			vecSpikeGroupConfidence(vecSpikeGroupTrialType==intStimType) = matPosteriorProbability(intStimType,vecSpikeGroupTrialType==intStimType);
		end
		
		%put in struct
		sSpikeGroup = struct;
		for intSpikeGroup=intSpikeGroupNum:-1:1
			sSpikeGroup(intSpikeGroup).Duration = vecSpikeGroupDuration(intSpikeGroup);
			sSpikeGroup(intSpikeGroup).Latency = vecSpikeGroupLatency(intSpikeGroup);
			sSpikeGroup(intSpikeGroup).TrialType = vecSpikeGroupTrialType(intSpikeGroup);
			sSpikeGroup(intSpikeGroup).TrialNumber = vecSpikeGroupTrialNumber(intSpikeGroup);
			sSpikeGroup(intSpikeGroup).DecodedType = vecSpikeGroupDecodedTrial(intSpikeGroup);
			sSpikeGroup(intSpikeGroup).Correct = vecSpikeGroupCorrect(intSpikeGroup);
			sSpikeGroup(intSpikeGroup).Confidence = vecSpikeGroupConfidence(intSpikeGroup);
		end
		return
		%% plot
		figure;maxfig;
		
		%make example
		intPlotTrial = 10;
		cellSpikesInTrial = cellUseSpikeTimesPerCellPerTrial(:,intPlotTrial);
		intSpikesInTrial = sum(sum(cellfun(@numel,cellSpikesInTrial)));
		vecSpikeTimes = nan(1,intSpikesInTrial);
		vecSpikeNeuron= nan(1,intSpikesInTrial);
		vecSpikeGroup = nan(1,intSpikesInTrial);
		intSpikeC = 1;
		for intN=1:numel(cellSpikesInTrial)
			vecSpikeT = cellSpikesInTrial{intN};
			intNumS = numel(vecSpikeT);
			vecAssign = intSpikeC:(intSpikeC+intNumS-1);
			intSpikeC = intSpikeC + intNumS;
			
			vecSpikeTimes(vecAssign) = vecSpikeT;
			vecSpikeNeuron(vecAssign) = intN;
		end
		
		%sort spikes
		[vecSpikeTimes,vecSort]=sort(vecSpikeTimes);
		vecSpikeNeuron = vecSpikeNeuron(vecSort);
		
		%assign to spike group
		intAssignGroups = floor(numel(vecSpikeTimes)/intSpikeGroupSize);
		vecGroupStart = nan(1,intAssignGroups);
		vecGroupStop = nan(1,intAssignGroups);
		for intGroup=1:intAssignGroups
			%get data
			intEndSpike = intGroup*intSpikeGroupSize;
			vecUseSpikes = (intEndSpike-intSpikeGroupSize+1):intEndSpike;
			vecSpikeGroup(vecUseSpikes) = intGroup;
			vecGroupStart(intGroup) = vecSpikeTimes(vecUseSpikes(1));
			vecGroupStop(intGroup) = vecSpikeTimes(vecUseSpikes(end));
		end
		vecSpikeGroup(isnan(vecSpikeGroup))=0;
		
		%match local groups to global groups
		vecGlobalGroups = find(vecSpikeGroupTrialNumber==intPlotTrial);
		
		% plot example trial
		%make plot
		h1=subplot(3,4,1);
		[vecT,vecR]=getIFR(vecSpikeTimes,0,1);
		plot(vecT,vecR,'color',[0.5 0.5 0.5]);
		hold on
		title(sprintf('%s, %s, trial %d',strRec,strType,intPlotTrial),'interpreter','none');
		xlabel('Time after stim onset (s)');
		ylabel('Pop IFR');		
		ylim([0 max(get(gca,'ylim'))]);
		
		matCol = [0 0 0; lines(intAssignGroups)];
		colormap(matCol);
		
		h2=subplot(3,4,5);
		hold on
		plot([vecGroupStart(1) vecGroupStop(end)],[1 1]*(1/intOriNum),'--','color',[0.5 0.5 0.5]);
		for intGroup=1:intAssignGroups
			plot(h2,[vecGroupStart(intGroup) vecGroupStop(intGroup)],[1 1]*vecSpikeGroupConfidence(vecGlobalGroups(intGroup)),'color',matCol(intGroup+1,:));
			intSp1 = find(vecT>=vecGroupStart(intGroup),1);
			intSp2 = find(vecT>=vecGroupStop(intGroup),1);
			plot(h1,vecT(intSp1:intSp2),vecR(intSp1:intSp2),'color',matCol(intGroup+1,:));
		end
		hold(h1,'off');
		hold(h2,'off');
		xlabel('Time after stim onset (s)');
		ylabel('Decoder confidence');		

		subplot(3,4,9);
		scatter(vecSpikeTimes,vecSpikeNeuron,[],vecSpikeGroup,'o','filled')
		xlabel('Time after stim onset (s)');
		ylabel('Neuron ID');		
		
		% plot groups of whole recording
		%sort
		[vecSortedDur,vecSort]=sort(vecSpikeGroupDuration);
		vecSortedCorr = vecSpikeGroupCorrect(vecSort);
		vecSortedConf = vecSpikeGroupConfidence(vecSort);
		vecSortedLat = vecSpikeGroupLatency(vecSort);
		%remove long durs
		%indRem = vecSortedDur > 0.05;
		%indRem = vecSortedLat < 0.15 | vecSortedDur > 0.05;
		indRem = [];
		vecSortedDur(indRem) = [];
		vecSortedCorr(indRem) = [];
		vecSortedConf(indRem) = [];
		vecSortedLat(indRem) = [];
		
		%calculate fraction correct and confidence per bin of equal size
		intBins = 10;
		intSperBin = floor(numel(vecSortedDur)/intBins);
		vecMeanDur = nan(1,intBins);
		vecSemDur = nan(1,intBins);
		vecMeanCorr = nan(1,intBins);
		matCiCorr = nan(2,intBins);
		cellValsCorr = cell(1,intBins);
		vecMeanConf = nan(1,intBins);
		vecSemConf = nan(1,intBins);
		cellValsConf = cell(1,intBins);
		vecSampleGroup = zeros(size(vecSortedDur));
		for intBin=1:intBins
			intEndS = intSperBin*intBin;
			vecSamples = (intEndS-intSperBin+1):intEndS;
			
			vecMeanDur(intBin) = mean(vecSortedDur(vecSamples));
			vecSemDur(intBin) = std(vecSortedDur(vecSamples))./sqrt(intSperBin);
			
			[phat,pci] = binofit(sum(vecSortedCorr(vecSamples)),intSperBin);
			vecMeanCorr(intBin) = phat;
			matCiCorr(:,intBin) = pci;
			cellValsCorr{intBin} = vecSortedCorr(vecSamples);
			
			vecMeanConf(intBin) = mean(vecSortedConf(vecSamples));
			vecSemConf(intBin) = std(vecSortedConf(vecSamples))./sqrt(intSperBin);
			cellValsConf{intBin} = vecSortedConf(vecSamples);
			vecSampleGroup(vecSamples) = intBin;
		end
		[r,p]=corr(vecSortedDur,vecSortedConf);
		[r2,p2]=corr(vecSortedDur,vecSortedCorr);
		
		subplot(2,4,2);
		hold on
		h=scatter(vecSortedDur*1000,vecSortedConf,100,[0.2 0.2 0.2],'.');
		h.MarkerFaceAlpha = 0.1;
		h.MarkerEdgeAlpha = 0.1;
		h2=scatter(vecSpikeGroupDuration(vecGlobalGroups)*1000,vecSpikeGroupConfidence(vecGlobalGroups),200,matCol(2:end,:),'.');
		hold off
		ylim([0 1]);
		xlim([0 200]);
		xlabel(sprintf('%d-spike block duration (ms)',intSpikeGroupSize));
		ylabel('Decoder confidence');
		
		subplot(2,4,6)
		hold on
		h=scatter(vecSpikeGroupLatency,vecSpikeGroupConfidence,100,[0.2 0.2 0.2],'.');
		h.MarkerFaceAlpha = 0.1;
		h.MarkerEdgeAlpha = 0.1;
		h2=scatter(vecSpikeGroupLatency(vecGlobalGroups),vecSpikeGroupConfidence(vecGlobalGroups),200,matCol(2:end,:),'.');
		hold off
		xlim([0 1]);
		ylim([0 1]);
		xlabel('Time after stim onset (s)');
		ylabel('Decoder confidence');
		
		subplot(2,4,8)
		vecBinsE = 0:0.002:max(vecSortedDur);
		vecBinsC = vecBinsE(2:end)-diff(vecBinsE(1:2))/2;
		[vecCounts,vecMeans,vecSDs] = makeBins(vecSortedDur,vecSortedCorr,vecBinsE);
		plot(vecBinsC*1000,(ones(size(vecBinsC))/intOriNum)*100,'--','color',[0.5 0.5 0.5]);
		hold on
		errorbar(vecBinsC*1000,vecMeans*100,(vecSDs*100)./sqrt(intSperBin),'color',lines(1));
		hold off
		xlabel(sprintf('%d-spike block duration (ms)',intSpikeGroupSize));
		ylabel('Decoder accuracy (%)');
		title(sprintf('r=%.3f, p=%.3f',r2,p2));
		ylim([0 max(get(gca,'ylim'))]);
		
		subplot(2,4,4)
		vecBinsE = 0:0.002:max(vecSortedDur);
		vecBinsC = vecBinsE(2:end)-diff(vecBinsE(1:2))/2;
		[vecCounts,vecMeans,vecSDs] = makeBins(vecSortedDur,vecSortedConf,vecBinsE);
		plot(vecBinsC*1000,(ones(size(vecBinsC))/intOriNum),'--','color',[0.5 0.5 0.5]);
		hold on
		errorbar(vecBinsC*1000,vecMeans,vecSDs./sqrt(intSperBin),'color',lines(1));
		hold off
		xlabel(sprintf('%d-spike block duration (ms)',intSpikeGroupSize));
		ylabel('Decoder confidence');
		title(sprintf('r=%.3f, p=%.3f',r,p));
		ylim([0 max(get(gca,'ylim'))]);
		
		%anova corr
		indRemS = vecSampleGroup==0;
		vecSampleGroup(indRemS) = [];
		vecSampleCorr = vecSortedCorr(~indRemS);
		[pA,table,stats] = anova1(vecSampleCorr,vecSampleGroup,'off');
		[c,~,~,gnames] = multcompare(stats,'CType','bonferroni','display','off');
		subplot(2,4,7)
		plot(vecMeanDur*1000,(ones(size(vecMeanDur))/intOriNum)*100,'--','color',[0.5 0.5 0.5]);
		hold on
		errorbar(vecMeanDur*1000,vecMeanCorr*100,(vecMeanCorr-matCiCorr(1,:))*100,(vecMeanCorr-matCiCorr(2,:))*100,vecSemDur,vecSemDur,'color',lines(1));
		hold off
		xlabel(sprintf('%d-spike block duration (ms)',intSpikeGroupSize));
		ylabel('Decoder accuracy (%)');
		title(sprintf('Deciles, ANOVA, p=%.1e',pA));
		ylim([0 max(get(gca,'ylim'))]);
		
		
		%anova conf
		indRemS = vecSampleGroup==0;
		vecSampleGroup(indRemS) = [];
		vecSampleConf = vecSortedConf(~indRemS);
		[pA2,table2,stats2] = anova1(vecSampleConf,vecSampleGroup,'off');
		[c2,~,~,gnames2] = multcompare(stats,'CType','bonferroni','display','off');
		
		subplot(2,4,3)
		plot(vecMeanDur*1000,(ones(size(vecMeanDur))/intOriNum),'--','color',[0.5 0.5 0.5]);
		hold on
		errorbar(vecMeanDur*1000,vecMeanConf,vecSemConf,vecSemConf,vecSemDur,vecSemDur,'color',lines(1));
		hold off
		xlabel(sprintf('%d-spike block duration (ms)',intSpikeGroupSize));
		ylabel('Decoder confidence');
		title(sprintf('Deciles, ANOVA, p=%.1e',pA2));
		ylim([0 max(get(gca,'ylim'))]);
		
		fixfig;
		
		%%
		export_fig(fullpath(strFigurePathSR,sprintf('A1c_PopActDynamics_%s_%s_SGS%s.tif',strRec,strType,strSGS)));
		export_fig(fullpath(strFigurePathSR,sprintf('A1c_PopActDynamics_%s_%s_SGS%s.pdf',strRec,strType,strSGS)));
		
		%% save data
		save(fullpath(strTargetDataPath,sprintf('Q1cData_%s_%s_SGS%s',strRec,strType,strSGS)),...
			'vecStimOnTime',...
			'vecStimOffTime',...
			'vecOrientation',...
			'vecOri180',...
			'intSpikeGroupSize',...
			'strRec',...
			'strType',...
			'sSpikeGroup');
	end
	close all;
end
toc
