%% aim
%{
Decode using sets of 20 spikes; is duration correlated with confidence of correct stim/accuracy?
I.e., are codes equally efficient during high and low rate periods?
%}

%% define qualifying areas
clear all;
runHeaderPopTimeCoding;

%% pre-allocate matrices
cellTypes = {'Real','ShuffTid','Uniform'};%, 'UniformTrial', 'ShuffTid'};%, 'PoissGain'}; %PoissGain not done
intNumTypes = numel(cellTypes);
boolFixSpikeGroupSize = false;
dblRemOnset = 0.125; %remove onset period in seconds

if strcmp(strRunType,'ABI')
	runLoadABI;
elseif strcmp(strRunType,'Sim')
	%get prepped sim files
	sSimRecs = dir(fullpath(strDataPathSim,'SimDG18_MatchedTo*.mat'));
	intRecNum = numel(sSimRecs);
else
	runLoadNpx;
end

%% go through recordings
if dblRemOnset == 0
	strOnset = '';
else
	strOnset = sprintf('%.2f',dblRemOnset);
end
tic

for intRec=1:intRecNum %19 || weird: 11
	%% prep ABI or Npx data
	if strcmp(strRunType,'ABI')
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
		
		%% move onset
		%remove first x ms
		vecStimOnTime = vecStimOnTime + dblRemOnset;
		
		%get data matrix
		[matData,indTuned,cellSpikeTimes,sOut,cellSpikeTimesPerCellPerTrial,indResp] = ...
			SimPrepData(cellSpikeTimesRaw,vecStimOnTime,vecStimOffTime,vecOrientation);
		
		%get cell props
		vecNeuronPrefOri = pi-[sLoad.sNeuron.PrefOri];
		vecNeuronType = [sLoad.sNeuron.Types]; %1=pyr,2=interneuron
		
	elseif strcmp(strRunType,'Npx')
		%prep
		runRecPrepNpx;
		strThisRec = strRec;
		strDataPathT0 = strTargetDataPath;
		
		%% move onset
		%remove first x ms
		vecStimOnTime = vecStimOnTime + dblRemOnset;
		
		%get data matrix
		[matData,indTuned,cellSpikeTimes,sTuning24,cellSpikeTimesPerCellPerTrial] = ...
			NpxPrepData(cellSpikeTimesRaw,vecStimOnTime,vecStimOffTime,vecOrientation);
		
		%get cell props
		vecNeuronPrefOri = [];
		vecNeuronType = []; %1=pyr,2=interneuron
		narrow vs broad not done
	end
	
	%get ori vars
	intTunedN = sum(indTuned);
	intRespN = size(matData,1);
	intNumN = numel(cellSpikeTimes);
	intTrialNum = numel(vecStimOnTime);
	vecOri180 = mod(vecOrientation,180);
	[vecOriIdx,vecUnique,vecPriorDistribution,cellSelect,vecRepetition] = val2idx(vecOri180);
	sTuning = getTuningCurves(matData,vecOri180,0);
	vecNeuronPrefOri = sTuning.matFittedParams(:,1);
	vecNeuronBandwidth = real(sTuning.matBandwidth);
	
	intOriNum = numel(vecUnique);
	intRepNum = min(vecPriorDistribution);
	dblStimDur = median(vecStimOffTime - vecStimOnTime);
	if mean(sum(matData)) < 90
		fprintf('Avg # of spikes per trial was %.1f for %s; skipping...\n',mean(sum(matData)),strThisRec);
		continue;
	end
	if boolFixSpikeGroupSize
		intSpikeGroupSize = 10;
		strSGS = ['Fixed' num2str(intSpikeGroupSize)];
	else
		intSpikeGroupSize = ceil(mean(sum(matData))/30); %20
		strSGS = ['Var' num2str(intSpikeGroupSize)];
	end
	
	%types: Real, UniformTrial, ShuffTid, PoissGain
	for intType=3%1:numel(cellTypes)
		strType = cellTypes{intType};
		
		
		%% load prepped data and ifr
		sSource = load(fullpath(strDataPathT0,sprintf('T0Data_%s%s%s%s%s',strThisRec,strRunType,strRunStim,strType)));
		vecTime = sSource.vecTime;
		vecIFR = sSource.vecIFR;
		vecAllSpikeTime = sSource.vecAllSpikeTime;
		vecAllSpikeNeuron = sSource.vecAllSpikeNeuron;
		cellSpikeTimes = sSource.cellSpikeTimes;
		
		if isempty(vecTime),continue;end
		
		%% replace IFR with binned rates?
		if 0
			%events
			dblStartEpoch = vecStimOnTime(1)-dblStimDur;
			dblEpochDur = vecStimOnTime(end)-vecStimOnTime(1)+dblStimDur;
			dblStopEpoch = dblStartEpoch + dblEpochDur;
			dblBinDur = (5/1000);
			vecBins=(dblStartEpoch:dblBinDur:dblStopEpoch)';
			vecIFR = flat(histcounts(vecTime,vecBins)./dblBinDur);
			vecTime = vecBins(2:end)-dblBinDur/2;
		end
		
		%% build trial-neuron cell matrix
		cellSpikeTimesPerCellPerTrial = cell(intNumN,intOrigTrialNum);
		for intN=1:intNumN
			%real
			[vecTrialPerSpike,vecTimePerSpike] = getSpikesInTrial(cellSpikeTimes{intN},vecStimOnTime,dblStimDur);
			for intTrial=1:intOrigTrialNum
				vecSpikeT = vecTimePerSpike(vecTrialPerSpike==intTrial);
				cellSpikeTimesPerCellPerTrial{intN,intTrial} = vecSpikeT;
			end
		end
		
		%% go through pop spikes and group into sets of 20
		fprintf('   Collecting n-spike groups and decoding [%s]\n',getTime);
		[sSpikeGroup,matSpikeGroupData] = getSpikeGroupData(cellUseSpikeTimesPerCellPerTrial,intSpikeGroupSize,vecOriIdx,vecStimOnTime,vecTime,vecIFR);
		intSpikeGroupNum = numel(sSpikeGroup);
		if intSpikeGroupNum < intTrialNum,continue;end
		
		%is confidence correlated with rate change?
		vecRateChanges = [sSpikeGroup.RateChange]';
		vecConfidence = [sSpikeGroup.Confidence]';
		vecCorrect = [sSpikeGroup.Correct]';
		vecSpikeGroupTrialNumber = [sSpikeGroup.TrialNumber]';
		vecSpikeGroupConfidence = [sSpikeGroup.Confidence]';
		vecSpikeGroupDuration = [sSpikeGroup.Duration]';
		vecSpikeGroupCorrect = [sSpikeGroup.Correct]';
		vecSpikeGroupLatency = [sSpikeGroup.Latency]';
		vecSpikeGroupAvgIFR = [sSpikeGroup.AvgRate]';
		
		vecEdgesX = linspace(min(vecRateChanges),max(vecRateChanges),100);
		vecEdgesY = linspace(min(vecSpikeGroupAvgIFR),max(vecSpikeGroupAvgIFR),100);
		vecBinsX = vecEdgesX(2:end) - diff(vecEdgesX(1:2))/2;
		vecBinsY = vecEdgesY(2:end) - diff(vecEdgesY(1:2))/2;
		[matCounts,matValMeans,matValSDs,cellVals,cellIDs] = makeBins2(vecRateChanges,vecSpikeGroupAvgIFR,vecConfidence,vecEdgesX,vecEdgesY);
		
		%% why is conf lower with high act? 1) untuned cells are also active; or 2) cells tuned to other oris are also active
		%gather information on tuning per cell; actually, can also run spike-by-spike correlation of IFR with tuning
		
		%result: during high ifr epochs, cells with low firing rate and low tuning are also active,
		%while during low ifr epochs, only cells with high firing rate or high tuning are active
		
		
		vecTuningPerCell = sTuning.vecFitT;
		vecRatePerCell = mean(sTuning.matMeanResp,2);
		vecMaxOriResp = max(sTuning.matMeanResp,[],2);
		matNormRespPerCellPerOri = sTuning.matMeanResp ./ vecMaxOriResp;
		vecSpikeGroupAvgTuningOfCells = nan(size(vecSpikeGroupLatency)); %t-statistic
		vecSpikeGroupAvgBandwidthOfCells = nan(size(vecSpikeGroupLatency)); %t-statistic
		vecSpikeGroupAvgRateOfCells = nan(size(vecSpikeGroupLatency)); %hz
		vecSpikeGroupAvgNormRespOfCellsToStim = nan(size(vecSpikeGroupLatency)); %normalized response to stim ori
		vecSpikeGroupNumOfCells = nan(size(vecSpikeGroupLatency)); %how many cells participate?
		vecSpikeGroupFractionInterneurons = nan(size(vecSpikeGroupLatency)); %fraction of interneurons
		vecSpikeGroupAvgPrefDistToStim = nan(size(vecSpikeGroupLatency)); %average distance to stim ori
		
		for intSpikeGroup=1:intSpikeGroupNum
			%get stim
			intTrial = sSpikeGroup(intSpikeGroup).TrialNumber;
			intStimIdx = vecOriIdx(intTrial);
			dblStimRad = deg2rad(vecUnique(intStimIdx));
			
			%get cells
			vecActiveCellSpikes = matSpikeGroupData(intSpikeGroup,:);
			vecSpikeGroupNumOfCells(intSpikeGroup) = sum(vecActiveCellSpikes>0);
			
			vecActiveCells = [];
			while max(vecActiveCellSpikes) > 0
				vecActiveCells = cat(2,vecActiveCells,find(vecActiveCellSpikes>0));
				vecActiveCellSpikes = vecActiveCellSpikes - 1;
			end
			vecSpikeGroupAvgRateOfCells(intSpikeGroup) = mean(vecRatePerCell(vecActiveCells));
			
			%get avg tuning
			vecSpikeGroupAvgTuningOfCells(intSpikeGroup) = mean(vecTuningPerCell(vecActiveCells));
			
			%get avg resp to stim
			vecSpikeGroupAvgNormRespOfCellsToStim(intSpikeGroup) = mean(matNormRespPerCellPerOri(vecActiveCells,intStimIdx));
			
			%get interneurons
			vecSpikeGroupFractionInterneurons(intSpikeGroup) = mean(vecNeuronType(vecActiveCells));
			
			%get avg pref dist to stim
			vecSpikeGroupAvgPrefDistToStim(intSpikeGroup) = mean(abs(circ_dist(2*vecNeuronPrefOri(vecActiveCells),dblStimRad*2)));
			vecSpikeGroupAvgBandwidthOfCells(intSpikeGroup) = mean(vecNeuronBandwidth(vecActiveCells));
			
			%assign to spike group
			sSpikeGroup(intSpikeGroup).AvgTuningOfCells = vecSpikeGroupAvgTuningOfCells(intSpikeGroup);
			sSpikeGroup(intSpikeGroup).AvgRateOfCells = vecSpikeGroupAvgRateOfCells(intSpikeGroup);
			sSpikeGroup(intSpikeGroup).AvgNormRespOfCellsToStim = vecSpikeGroupAvgNormRespOfCellsToStim(intSpikeGroup);
			sSpikeGroup(intSpikeGroup).NumOfCells = vecSpikeGroupNumOfCells(intSpikeGroup);
			
			sSpikeGroup(intSpikeGroup).FractionInterneurons = vecSpikeGroupFractionInterneurons(intSpikeGroup);
			sSpikeGroup(intSpikeGroup).AvgPrefDistToStim = vecSpikeGroupAvgPrefDistToStim(intSpikeGroup);
			sSpikeGroup(intSpikeGroup).AvgBandwidthOfCells = vecSpikeGroupAvgBandwidthOfCells(intSpikeGroup);
		end
		
		%{
		figure
		subplot(2,3,1)
		scatter(vecSpikeGroupAvgIFR,vecSpikeGroupNumOfCells);
		[r,p]=corr(vecSpikeGroupAvgIFR,vecSpikeGroupNumOfCells);
		title(sprintf('avg ifr vs # of cells, r=%.3f, p=%.2e',r,p));
		
		subplot(2,3,2)
		scatter(vecSpikeGroupAvgIFR,vecSpikeGroupAvgRateOfCells);
		[r,p]=corr(vecSpikeGroupAvgIFR,vecSpikeGroupAvgRateOfCells);
		title(sprintf('avg ifr vs avg rate of cells, r=%.3f, p=%.2e',r,p));
		
		subplot(2,3,3)
		scatter(vecSpikeGroupAvgIFR,vecSpikeGroupAvgTuningOfCells);
		[r,p]=corr(vecSpikeGroupAvgIFR,vecSpikeGroupAvgTuningOfCells);
		title(sprintf('avg ifr vs avg tuning of cells, r=%.3f, p=%.2e',r,p));
		
		subplot(2,3,4)
		scatter(vecSpikeGroupAvgIFR,vecSpikeGroupAvgNormRespOfCellsToStim);
		[r,p]=corr(vecSpikeGroupAvgIFR,vecSpikeGroupAvgNormRespOfCellsToStim);
		title(sprintf('avg ifr vs avg resp to stim of cells, r=%.3f, p=%.2e',r,p));
		%}
		%predict conf
		vecSpikeGroupAvgIFR_log=vecSpikeGroupAvgIFR;
		[r_Num_Conf,p_Num_Conf]=corr(vecSpikeGroupNumOfCells,vecSpikeGroupConfidence);
		[r_Rate_Conf,p_Rate_Conf]=corr(vecSpikeGroupAvgRateOfCells,vecSpikeGroupConfidence);
		[r_Tune_Conf,p_Tune_Conf]=corr(vecSpikeGroupAvgTuningOfCells,vecSpikeGroupConfidence);
		%[r_Resp_Conf,p_Resp_Conf]=corr(vecSpikeGroupAvgNormRespOfCellsToStim,vecSpikeGroupConfidence);
		[r_IFR_Conf,p_IFR_Conf]=corr(vecSpikeGroupAvgIFR_log,vecSpikeGroupConfidence);
		
		tbl = table(vecSpikeGroupNumOfCells,vecSpikeGroupAvgRateOfCells,vecSpikeGroupAvgTuningOfCells,vecSpikeGroupAvgIFR_log,vecSpikeGroupConfidence,...
			'VariableNames',{'CellNumber','AvgCellRate','AvgCellTuning','AvgPopIFR','Conf'});
		mdl_predconf = fitlm(tbl,'linear');
		
		%predict ifr
		[r_Num_IFR,p_Num_IFR]=corr(vecSpikeGroupNumOfCells,vecSpikeGroupAvgIFR_log);
		[r_Rate_IFR,p_Rate_IFR]=corr(vecSpikeGroupAvgRateOfCells,vecSpikeGroupAvgIFR_log);
		[r_Tune_IFR,p_Tune_IFR]=corr(vecSpikeGroupAvgTuningOfCells,vecSpikeGroupAvgIFR_log);
		tbl = table(vecSpikeGroupAvgRateOfCells,vecSpikeGroupAvgTuningOfCells,vecSpikeGroupAvgIFR_log,...
			'VariableNames',{'AvgCellRate','AvgCellTuning','AvgPopIFR'});
		mdl_predifr = fitlm(tbl,'linear');
		%scatter(vecSpikeGroupAvgIFR_log,mdl_predifr.Fitted)
		
		%rate vs tuning
		[rGroup_Rate_Tune,pGroup_Rate_Tune]=corr(vecTuningPerCell,vecRatePerCell);
		
		%save data: r_Num_IFR, r_Rate_IFR, r_Tune_IFR
		
		%result: during high ifr epochs, cells with low firing rate and low tuning are also active,
		%while during low ifr epochs, only cells with high firing rate or high tuning are active
		
		%% run single-spike analysis
		%vecSpikeCellTuning = vecTuningPerCell(vecAllSpikeNeuron);
		%vecSpikeCellRate = vecRatePerCell(vecAllSpikeNeuron);
		%[rSpike_IFR_Rate,pSpike_IFR_Rate]=corr(vecIFR,vecSpikeCellRate);
		%[rSpike_IFR_Tune,pSpike_IFR_Tune]=corr(vecIFR,vecSpikeCellTuning);
		
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
		
		%% plot example trial
		%make plot
		h1=subplot(3,4,1);
		[vecT,vecR]=getIFR(vecSpikeTimes,0,1);
		plot(dblRemOnset+vecT,vecR,'color',[0.5 0.5 0.5]);
		hold on
		title(sprintf('%s, %s, %s trial %d',strThisRec,strType,strOnset,intPlotTrial),'interpreter','none');
		xlabel('Time after stim onset (s)');
		ylabel('Pop IFR');
		ylim([0 max(get(gca,'ylim'))]);
		
		matCol = [0 0 0; lines(intAssignGroups)];
		colormap(matCol);
		
		h2=subplot(3,4,5);
		hold on
		plot(dblRemOnset+[vecGroupStart(1) vecGroupStop(end)],[1 1]*(1/intOriNum),'--','color',[0.5 0.5 0.5]);
		for intGroup=1:intAssignGroups
			plot(h2,dblRemOnset+[vecGroupStart(intGroup) vecGroupStop(intGroup)],[1 1]*vecSpikeGroupConfidence(vecGlobalGroups(intGroup)),'color',matCol(intGroup+1,:));
			intSp1 = find(vecT>=vecGroupStart(intGroup),1);
			intSp2 = find(vecT>=vecGroupStop(intGroup),1);
			plot(h1,dblRemOnset+vecT(intSp1:intSp2),vecR(intSp1:intSp2),'color',matCol(intGroup+1,:));
		end
		hold(h1,'off');
		hold(h2,'off');
		xlim(h1,[0 max(get(h1,'xlim'))]);
		xlim(h2,[0 max(get(h2,'xlim'))]);
		xlabel('Time after stim onset (s)');
		ylabel('Decoder confidence');
		
		subplot(3,4,9);
		scatter(dblRemOnset+vecSpikeTimes,vecSpikeNeuron,[],vecSpikeGroup,'o','filled')
		xlim([0 max(get(gca,'xlim'))]);
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
		
		%calculate confidence per bin of equal size
		intQuantileNum = 10;
		[vecMeanDur,vecSemDur,vecMeanConf,vecSemConf,vecQuantile]=getQuantiles(vecSortedDur,vecSortedConf,intQuantileNum);
		
		%calculate accuracy per bin of equal size
		intSperBin = floor(numel(vecSortedDur)/intQuantileNum);
		vecMeanCorr = nan(1,intQuantileNum);
		matCiCorr = nan(2,intQuantileNum);
		cellValsCorr = cell(1,intQuantileNum);
		vecSampleGroup = zeros(size(vecSortedDur));
		for intBin=1:intQuantileNum
			intEndS = intSperBin*intBin;
			vecSamples = (intEndS-intSperBin+1):intEndS;
			
			[phat,pci] = binofit(sum(vecSortedCorr(vecSamples)),intSperBin);
			vecMeanCorr(intBin) = phat;
			matCiCorr(:,intBin) = pci;
			cellValsCorr{intBin} = vecSortedCorr(vecSamples);
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
		xlim([0 dblStimDur+dblRemOnset]);
		ylim([0 1]);
		xlabel('Time after stim onset (s)');
		ylabel('Decoder confidence');
		
		subplot(2,4,8)
		vecBinsE = 0:0.002:max(vecSortedDur);
		vecBinsC = vecBinsE(2:end)-diff(vecBinsE(1:2))/2;
		[vecCounts,vecMeans,vecSDs] = makeBins(vecSortedDur,vecSortedCorr,vecBinsE);
		plot(vecBinsC*1000,(ones(size(vecBinsC))/intOriNum)*100,'--','color',[0.5 0.5 0.5]);
		hold on
		errorbar(vecBinsC*1000,vecMeans*100,(vecSDs*100)./sqrt(vecCounts),'color',lines(1));
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
		errorbar(vecBinsC*1000,vecMeans,vecSDs./sqrt(vecCounts),'color',lines(1));
		hold off
		xlabel(sprintf('%d-spike block duration (ms)',intSpikeGroupSize));
		ylabel('Decoder confidence');
		title(sprintf('r=%.3f, p=%.3f',r,p));
		ylim([0 max(get(gca,'ylim'))]);
		
		%anova corr
		indRemS = vecSampleGroup==0;
		vecSampleGroup(indRemS) = [];
		vecSampleCorr = vecSortedCorr(~indRemS);
		[pA,tbl,stats] = anova1(vecSampleCorr,vecSampleGroup,'off');
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
		export_fig(fullpath(strFigurePathSR,sprintf('A1c_PopActDynamics_%s_%s_SGS%s%s.tif',strThisRec,strType,strSGS,strOnset)));
		export_fig(fullpath(strFigurePathSR,sprintf('A1c_PopActDynamics_%s_%s_SGS%s%s.pdf',strThisRec,strType,strSGS,strOnset)));
		
		%% save data
		save(fullpath(strTargetDataPath,sprintf('Q1cData_%s_%s_SGS%s%s.mat',strThisRec,strType,strSGS,strOnset)),...
			'vecStimOnTime',...
			'vecStimOffTime',...
			'vecOrientation',...
			'vecOri180',...
			'intSpikeGroupSize',...
			'strRec',...
			'strType',...
			'vecTuningPerCell',...
			'vecRatePerCell',...
			...'rSpike_IFR_Rate',...
			...'pSpike_IFR_Rate',...
			...'rSpike_IFR_Tune',...
			...'pSpike_IFR_Tune',...
			'sSpikeGroup');
		
	end
	close all;
end
toc
