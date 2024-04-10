%% aim
%{
Decode using sets of 20 spikes; is duration correlated with confidence of correct stim/accuracy?
I.e., are codes equally efficient during high and low rate periods?
%}


%% set parameters
cellDataTypes = {'Npx','Sim','ABI','SWN'};%topo, model, allen, nora
intRunDataType = 1;
strRunStim = 'DG';%DG or NM? => superseded to WS by SWN
cellTypes = {'Real','ShuffTid','Uniform'};%, 'UniformTrial', 'ShuffTid'};%, 'PoissGain'}; %PoissGain not done
dblRemOnset = 0; %remove onset period in seconds; 0.125 for sim, 0.25 for npx
runHeaderPopTimeCoding;
cellSupraGranuInfra = {'Supragranular','Granular','Infragranular'};

%% go through recordings
tic
for boolFixSpikeGroupSize = false%[true false]
for dblRemOnset = 0.25%[0 0.25] %remove onset period in seconds
for intRec=1:numel(sAggStim) %19 || weird: 11
	%% prep ABI or Npx data
	if strcmp(strRunType,'ABI')
		error to be updated
		runRecPrepABI;
		strThisRec = strRec;
	elseif strcmp(strRunType,'Sim')
		%load
		runRecPrepSim;
		
		%edit vars
		strThisRec = strRec;
		strDataPathT0=strDataPathSimT0;
		vecOri180 = mod(vecOrientation,180)*2;
		vecStimIdx = vecOri180;
		
		%get cell props
		%vecNeuronPrefOri = pi-[sLoad.sNeuron.PrefOri];
		vecNeuronType = [sLoad.sNeuron.Types]; %1=pyr,2=interneuron
		
	elseif strcmp(strRunType,'Npx') || strcmp(strRunType,'SWN')
		%prep
		runRecPrepNpx;
		strThisRec = strRec;
		strDataPathT0 = strTargetDataPath;
		
		
		%get cell props
		vecNeuronType = ones(size(indTuned)); %1=pyr,2=interneuron
		%narrow vs broad not done
	else
		error impossible
	end
	
	%% get cortical depth per cell
	sArea1Neurons = sUseNeuron(indArea1Neurons);
	cellAreas = {sArea1Neurons.Area};
	vecCorticalLayer = cellfun(@(x) str2double(x(regexp(x,'layer.*')+6)),cellAreas);
	vecDepth = [sArea1Neurons.DepthBelowIntersect];
	vecSupraGranuInfra = double(vecCorticalLayer < 4) + 2*double(vecCorticalLayer == 4) + 3*double(vecCorticalLayer > 4);
	%scatter(vecSupraGranuInfra,vecDepth)
	%hold on
	%text(vecSupraGranuInfra,vecDepth,cellAreas)
	for intCortLayer = 4%1:3
		if intCortLayer == 4
			indUseNeurons = true(size(vecSupraGranuInfra));
			strLayer = '';
		else
			indUseNeurons = vecSupraGranuInfra(:)==intCortLayer;
			strLayer = cellSupraGranuInfra{intCortLayer};
		end
	sArea1Neurons = sUseNeuron(indArea1Neurons & indUseNeurons);
	
	%% move onset
	%remove first x ms
	vecStimOnTime = vecStimOnTime + dblRemOnset;
	
	%% load prepped data
	strTarget = fullpath(strDataPathT0,[sprintf('T0Data_%s%s%s%s%s',strThisRec,strRunType,strRunStim,'Real') '.mat']);
	if ~exist(strTarget,'file')
		fprintf('Prepped T0 file did not exist for %s; skipping...\n',strThisRec);
		continue;
	end
	
	%check fr
	sSource = load(strTarget);
	cellSpikeTimes = sSource.cellSpikeTimes;
	intNumN = numel(cellSpikeTimes);
	cellSpikeTimesPerCellPerTrial = cell(intNumN,intOrigTrialNum);
	for intN=1:intNumN
		%real
		[vecTrialPerSpike,vecTimePerSpike] = getSpikesInTrial(cellSpikeTimes{intN},vecStimOnTime,dblStimDur);
		for intTrial=1:numel(vecStimOnTime)
			vecSpikeT = vecTimePerSpike(vecTrialPerSpike==intTrial);
			cellSpikeTimesPerCellPerTrial{intN,intTrial} = vecSpikeT;
		end
	end
	matData = cellfun(@numel,cellSpikeTimesPerCellPerTrial)./dblStimDur;
	if mean(sum(matData)) < 90%90 / 50
		fprintf('Avg # of spikes per trial was %.1f for %s; skipping...\n',mean(sum(matData)),strThisRec);
		continue;
	end
	
	
	%% get ori vars
	intTrialNum = numel(vecStimOnTime);
	vecOri180 = mod(vecOrientation,180);
	[vecOriIdx,vecUnique,vecPriorDistribution,cellSelect,vecRepetition] = val2idx(vecOri180);
	intOriNum = numel(vecUnique);
	intRepNum = min(vecPriorDistribution);
	dblStimDur = median(vecStimOffTime - vecStimOnTime);
	
	for intType=1:numel(cellTypes)
		%% load prepro T0 data
		% pool spikes from all neurons, but save the time+id per spike, then calculate IFR over all
		% spikes at pop level, detect peaks, and split into trials.
		%get spikes per trial per neuron
		strType = cellTypes{intType};
		fprintf('Running rec %d/%d for %s, %s, %s, %s, Rem %.3f s; %s [%s]\n',intRec,numel(sAggStim),strRunType,strLayer,strRunStim,strType,dblRemOnset,strRec,getTime);
		
		%% load ifr
		sSource = load(fullpath(strTargetDataPath,sprintf('T0Data_%s%s%s%s%s',strRec,strRunType,strRunStim,strType,strLayer)));
		vecTime = sSource.vecTime;
		vecIFR = sSource.vecIFR;
		vecAllSpikeTime = sSource.vecAllSpikeTime;
		vecAllSpikeNeuron = sSource.vecAllSpikeNeuron;
		if isempty(vecTime),continue;end
		
		%check fr
		cellSpikeTimes = sSource.cellSpikeTimes;
		intNumN = numel(cellSpikeTimes);
		cellSpikeTimesPerCellPerTrial = cell(intNumN,intOrigTrialNum);
		for intN=1:intNumN
			%real
			[vecTrialPerSpike,vecTimePerSpike] = getSpikesInTrial(cellSpikeTimes{intN},vecStimOnTime,dblStimDur);
			for intTrial=1:numel(vecStimOnTime)
				vecSpikeT = vecTimePerSpike(vecTrialPerSpike==intTrial);
				cellSpikeTimesPerCellPerTrial{intN,intTrial} = vecSpikeT;
			end
		end
		matData = cellfun(@numel,cellSpikeTimesPerCellPerTrial)./dblStimDur;
		if mean(sum(matData)) < 90%90 / 50
			fprintf('Avg # of spikes per trial was %.1f for %s; skipping...\n',mean(sum(matData)),strThisRec);
			continue;
		end
	
		%% check rates
		if boolFixSpikeGroupSize
			intSpikeGroupSize = 10;
			strSGS = ['Fixed' num2str(intSpikeGroupSize)];
		else
			intSpikeGroupSize = ceil(mean(sum(matData))/30); %20
			strSGS = ['Var' num2str(intSpikeGroupSize)];
		end
		
		%% go through pop spikes and group into sets of 20
		fprintf('   Collecting n-spike groups and decoding [%s]\n',getTime);
		sSpikeGroup = getSpikeGroupData(cellSpikeTimesPerCellPerTrial,intSpikeGroupSize,vecOriIdx,vecStimOnTime,vecTime,vecIFR);
		if isempty(sSpikeGroup),continue;end
		
		%is confidence correlated with rate change?
		vecRateChanges = [sSpikeGroup.RateChange]';
		vecConfidence = [sSpikeGroup.Confidence]';
		vecCorrect = [sSpikeGroup.Correct]';
		vecAvgRate = [sSpikeGroup.AvgRate]';
		vecSpikeGroupTrialNumber = [sSpikeGroup.TrialNumber]';
		vecSpikeGroupConfidence = [sSpikeGroup.Confidence]';
		vecSpikeGroupDuration = [sSpikeGroup.Duration]';
		vecSpikeGroupCorrect = [sSpikeGroup.Correct]';
		vecSpikeGroupLatency = [sSpikeGroup.Latency]';
		
		
		vecEdgesX = linspace(min(vecRateChanges),max(vecRateChanges),100);
		vecEdgesY = linspace(min(vecAvgRate),max(vecAvgRate),100);
		vecBinsX = vecEdgesX(2:end) - diff(vecEdgesX(1:2))/2;
		vecBinsY = vecEdgesY(2:end) - diff(vecEdgesY(1:2))/2;
		[matCounts,matValMeans,matValSDs,cellVals,cellIDs] = makeBins2(vecRateChanges,vecAvgRate,vecConfidence,vecEdgesX,vecEdgesY);
		
		%% calc deciles
		fprintf('   Separating into quantiles [%s]\n',getTime);
		[vecSortedDeltaRate,vecSort] = sort(vecRateChanges);
		vecSortedConfidence = vecConfidence(vecSort);
		vecSortedCorrect = vecCorrect(vecSort);
		vecSortedNormConfidence = zscore(vecSortedConfidence);
		
		%calculate fraction correct and confidence per bin of equal size
		intBins = 10;
		intSperBin = floor(numel(vecSortedDeltaRate)/intBins);
		vecMeanDeltaR = nan(1,intBins);
		vecSemDeltaR = nan(1,intBins);
		vecMeanCorrect = nan(1,intBins);
		matCiCorrect = nan(2,intBins);
		cellValsCorrect = cell(1,intBins);
		vecMeanConf = nan(1,intBins);
		vecSemConf = nan(1,intBins);
		cellValsConf = cell(1,intBins);
		vecMeanNormConf = nan(1,intBins);
		vecSemNormConf = nan(1,intBins);
		cellValsNormConf = cell(1,intBins);
		vecSampleGroup = zeros(size(vecSortedDeltaRate));
		for intBin=1:intBins
			intEndS = intSperBin*intBin;
			vecSamples = (intEndS-intSperBin+1):intEndS;
			
			vecMeanDeltaR(intBin) = mean(vecSortedDeltaRate(vecSamples));
			vecSemDeltaR(intBin) = std(vecSortedDeltaRate(vecSamples))./sqrt(intSperBin);
			
			[phat,pci] = binofit(sum(vecSortedCorrect(vecSamples)),intSperBin);
			vecMeanCorrect(intBin) = phat;
			matCiCorrect(:,intBin) = pci;
			cellValsCorrect{intBin} = vecSortedCorrect(vecSamples);
			
			vecMeanConf(intBin) = mean(vecSortedConfidence(vecSamples));
			vecSemConf(intBin) = std(vecSortedConfidence(vecSamples))./sqrt(intSperBin);
			cellValsConf{intBin} = vecSortedConfidence(vecSamples);
			
			vecMeanNormConf(intBin) = mean(vecSortedNormConfidence(vecSamples));
			vecSemNormConf(intBin) = std(vecSortedNormConfidence(vecSamples))./sqrt(intSperBin);
			cellValsNormConf{intBin} = vecSortedNormConfidence(vecSamples);
			
			vecSampleGroup(vecSamples) = intBin;
		end
		
		%% plot
		
		%make example
		intPlotTrial = 10;
		vecSpPerT = sum(cellfun(@numel,cellSpikeTimesPerCellPerTrial),1);
		if vecSpPerT(intPlotTrial) > intSpikeGroupSize
		figure;maxfig;
		cellSpikesInTrial = cellSpikeTimesPerCellPerTrial(:,intPlotTrial);
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
		plot(vecT+dblRemOnset,vecR,'color',[0.5 0.5 0.5]);
		hold on
		title(sprintf('%s, %s, %s, trial %d',strRec,strType,strLayer,intPlotTrial),'interpreter','none');
		xlabel('Time after stim onset (s)');
		ylabel('Pop IFR');		
		ylim([0 max(get(gca,'ylim'))]);
		xlim([0 max(get(gca,'xlim'))]);
		
		matCol = [0 0 0; lines(intAssignGroups)];
		colormap(matCol);
		
		h2=subplot(3,4,5);
		hold on
		plot([vecGroupStart(1) vecGroupStop(end)],[1 1]*(1/intOriNum),'--','color',[0.5 0.5 0.5]);
		for intGroup=1:intAssignGroups
			plot(h2,dblRemOnset+[vecGroupStart(intGroup) vecGroupStop(intGroup)],[1 1]*vecSpikeGroupConfidence(vecGlobalGroups(intGroup)),'color',matCol(intGroup+1,:));
			intSp1 = find(vecT>=vecGroupStart(intGroup),1);
			intSp2 = find(vecT>=vecGroupStop(intGroup),1);
			plot(h1,dblRemOnset+vecT(intSp1:intSp2),vecR(intSp1:intSp2),'color',matCol(intGroup+1,:));
		end
		hold(h1,'off');
		hold(h2,'off');
		xlabel('Time after stim onset (s)');
		ylabel('Decoder confidence');		
		xlim([0 max(get(gca,'xlim'))]);
		
		subplot(3,4,9);
		scatter(dblRemOnset+vecSpikeTimes,vecSpikeNeuron,[],vecSpikeGroup,'o','filled')
		xlabel('Time after stim onset (s)');
		ylabel('Neuron ID');		
		xlim([0 max(get(gca,'xlim'))]);
		
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
		vecMeanCorrect = nan(1,intBins);
		matCiCorrect = nan(2,intBins);
		cellValsCorrect = cell(1,intBins);
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
			vecMeanCorrect(intBin) = phat;
			matCiCorrect(:,intBin) = pci;
			cellValsCorrect{intBin} = vecSortedCorr(vecSamples);
			
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
		h=scatter(dblRemOnset+vecSpikeGroupLatency,vecSpikeGroupConfidence,100,[0.2 0.2 0.2],'.');
		h.MarkerFaceAlpha = 0.1;
		h.MarkerEdgeAlpha = 0.1;
		h2=scatter(dblRemOnset+vecSpikeGroupLatency(vecGlobalGroups),vecSpikeGroupConfidence(vecGlobalGroups),200,matCol(2:end,:),'.');
		hold off
		xlim([0 max(get(gca,'xlim'))]);
		ylim([0 1]);
		xlabel('Time after stim onset (s)');
		ylabel('Decoder confidence');
		
		
		%dur/confidence
		indRemS = vecSampleGroup==0;
		vecSampleGroup(indRemS) = [];
		vecSampleConf = vecSortedConf(~indRemS);
		[pA2,table2,stats2] = anova1(vecSampleConf,vecSampleGroup,'off');
		[c2,~,~,gnames2] = multcompare(stats2,'CType','bonferroni','display','off');
		
		subplot(2,4,3)
		plot(vecMeanDur*1000,(ones(size(vecMeanDur))/intOriNum),'--','color',[0.5 0.5 0.5]);
		hold on
		errorbar(vecMeanDur*1000,vecMeanConf,vecSemConf,vecSemConf,vecSemDur,vecSemDur,'color',lines(1));
		hold off
		xlabel(sprintf('%d-spike block duration (ms)',intSpikeGroupSize));
		ylabel('Decoder confidence');
		title(sprintf('Deciles, ANOVA, p=%.1e',pA2));
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
		
		
		%delta-r/confidence
		subplot(2,4,7)
		plot(vecMeanDeltaR,(zeros(size(vecMeanDeltaR))/intOriNum),'--','color',[0.5 0.5 0.5]);
		hold on
		errorbar(vecMeanDeltaR,vecMeanNormConf,vecSemNormConf,vecSemNormConf,vecSemDeltaR,vecSemDeltaR,'color',lines(1));
		hold off
		xlabel(sprintf('Rate change during spike block (Hz)'));
		ylabel('Relative decoder confidence (z-score)');
		ylim([-1 1]*max(abs(get(gca,'ylim'))));
		
		subplot(2,4,8)
		vecEdgesX = -500:10:500;
		vecCenterX = vecEdgesX(2:end) - diff(vecEdgesX(1:2))/2;
		[vecCounts,vecMeans,vecSDs,cellVals,cellIDs] = makeBins(vecRateChanges,vecConfidence,vecEdgesX);
		
		[r,p]=corr(vecRateChanges,vecConfidence);
		plot(vecCenterX,(ones(size(vecCenterX))/intOriNum),'--','color',[0.5 0.5 0.5]);
		hold on
		errorbar(vecCenterX,vecMeans,vecSDs./sqrt(vecCounts),'color',lines(1));
		hold off
		xlabel(sprintf('Rate change during spike block (Hz)'));
		ylabel('Decoder confidence');
		ylim([0 max(get(gca,'ylim'))]);
		title(sprintf('Pearson r=%.3f, p=%.3f',r,p));
		
		fixfig;
		
		%%
		strRemDur = sprintf('%.0f',dblRemOnset*1000);
		export_fig(fullpath(strFigurePathSR,sprintf('A1d_PopActDynamics_%s_%s_SGS%s_%s_%s_%s_Onset%s.tif',strRec,strType,strSGS,strRunType,strRunStim,strLayer,strRemDur)));
		export_fig(fullpath(strFigurePathSR,sprintf('A1d_PopActDynamics_%s_%s_SGS%s_%s_%s_%s_Onset%s.pdf',strRec,strType,strSGS,strRunType,strRunStim,strLayer,strRemDur)));
		end
		%% save data
		strRemDur = sprintf('%.0f',dblRemOnset*1000);
		save(fullpath(strTargetDataPath,sprintf('Q1dData_%s_%s_SGS%s_%s_%s_%s_Onset%s',strRec,strType,strSGS,strRunType,strRunStim,strLayer,strRemDur)),...
			'vecStimOnTime',...
			'vecStimOffTime',...
			'vecOrientation',...
			'vecOri180',...
			'intSpikeGroupSize',...
			'strRec',...
			'strType',...
			'strSGS',...
			'strLayer',...
			'strRunType',...
			'dblRemOnset',...
			'sSpikeGroup');
	end
	close all;
end
end
end
end
toc
%{
%% detect peaks
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
		fprintf('   Detecting peaks [%s]\n',getTime);
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
		
		%% go through peaks and decode leading and lagging spikes
		intSpikeGroupSize=20
		intPeakNum = numel(vecMergedPopEventLocs);
		matPeakGroupData = nan(intPeakNum*2,intNumN);
		vecPeakGroupTrial = nan(intPeakNum*2,1);
		vecPeakGroupStimType = nan(intPeakNum*2,1);
		vecPeakGroupIsPost = nan(intPeakNum*2,1);
		sPeakGroup = struct;
		for intPeakGroup = 1:intPeakNum
			%get trial
			intGroupCenter = vecMergedPopEventLocs(intPeakGroup);
			dblGroupTime = vecMergedPopEventTimes(intPeakGroup);
			intTrial = sum(dblGroupTime>vecStimOnTime);
			if isempty(intTrial),intTrial=0;end
			if intTrial==0 || dblGroupTime>vecStimOffTime(intTrial)
				%iti
				intStimType = 0;
				dblTrialStart = dblStartEpoch;
			else
				%during
				intStimType = vecOriIdx(intTrial);
				dblTrialStart = vecStimOnTime(intTrial);
			end
			%pre
			vecSpikesPre = (intGroupCenter-intSpikeGroupSize):(intGroupCenter-1);
			vecPreNeurons = vecAllSpikeNeuron(vecSpikesPre);
			
			matPeakGroupData(intPeakGroup*2-1,:) = accumarray(vecPreNeurons(:),1,[intNumN 1]);
			vecPeakGroupTrial(intPeakGroup*2-1) = intTrial;
			vecPeakGroupStimType(intPeakGroup*2-1) = intStimType;
			vecPeakGroupIsPost(intPeakGroup*2-1) = 0;
			
			%post
			vecSpikesPost = (intGroupCenter+1):(intGroupCenter+intSpikeGroupSize);
			vecPostNeurons = vecAllSpikeNeuron(vecSpikesPost);
			
			matPeakGroupData(intPeakGroup*2,:) = accumarray(vecPostNeurons(:),1,[intNumN 1]);
			vecPeakGroupTrial(intPeakGroup*2) = intTrial;
			vecPeakGroupStimType(intPeakGroup*2) = intStimType;
			vecPeakGroupIsPost(intPeakGroup*2) = 1;
		end
		
		%reduce data to only events in trials
		indRem = vecPeakGroupStimType==0;
		matPeakGroupData(indRem,:)=[];
		vecPeakGroupTrial(indRem)=[];
		vecPeakGroupStimType(indRem)=[];
		vecPeakGroupIsPost(indRem)=[];
		intDoublePeakGroup = numel(vecPeakGroupIsPost);
		
		%do logistic regression
		vecPriorDistribution = accumarray(vecPeakGroupStimType(:),1,[intOriNum 1]);
		intVerbose = 1;
		[dblPeakPerformanceCV,vecPeakGroupDecodedTrial,matPeakPosteriorProbability,dblPeakMeanErrorDegs,matPeakConfusion,matPeakWeights] = ...
			doCrossValidatedDecodingLR(matPeakGroupData,vecPeakGroupStimType,intTypeCV,[],dblLambda,intVerbose);
		
		vecPeakGroupCorrect = vecPeakGroupDecodedTrial(:) == vecPeakGroupStimType;
		vecPeakGroupConfidence = nan(intDoublePeakGroup,1);
		for intStimType=1:intOriNum
			vecPeakGroupConfidence(vecPeakGroupStimType==intStimType) = matPeakPosteriorProbability(intStimType,vecPeakGroupStimType==intStimType);
		end
		
		%
		% plot histo
		dblStep = 0.01;
		dblMax = max(vecPeakGroupConfidence);
		vecBinsE = 0:dblStep:dblMax;
		vecConfPre = histcounts(vecPeakGroupConfidence(vecPeakGroupIsPost==0),vecBinsE);
		vecConfPost = histcounts(vecPeakGroupConfidence(vecPeakGroupIsPost==1),vecBinsE);
		vecBinsC = vecBinsE(2:end) - dblStep;
		
		figure
		subplot(2,3,1)
		plot(vecBinsC,vecConfPre,'b');
		hold on
		plot(vecBinsC,vecConfPost,'r');
		hold off
		
		%diff
		subplot(2,3,2)
		vecPeakGroupConfDiff = vecPeakGroupConfidence(2:2:end) - vecPeakGroupConfidence(1:2:end);
		[h,p]=ttest(vecPeakGroupConfDiff);
		dblStep = 0.01;
		dblMax = 0.3;%max(abs(vecPeakGroupConfDiff));
		vecBinsE = -dblMax:dblStep:dblMax;
		vecBinsC = vecBinsE(2:end) - dblStep;
		vecConfDiffCount = histcounts(vecPeakGroupConfDiff,vecBinsE);
		plot(vecBinsC,vecConfDiffCount,'k');
		title(sprintf('%d, mu=%.3f, p=%.3f',intSpikeGroupSize,mean(vecPeakGroupConfDiff),p));
		
		return
%}