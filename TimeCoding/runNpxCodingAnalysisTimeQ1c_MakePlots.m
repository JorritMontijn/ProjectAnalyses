%% aim
%{
Decode using sets of 20 spikes; is duration correlated with confidence of correct stim/accuracy?
I.e., are codes equally efficient during high and low rate periods?
%}

%% set parameters
cellDataTypes = {'Npx','Sim','ABI','SWN'};%topo, model, allen, nora
intRunDataType = 3;
strRunStim = 'NS';%DG or NM or NS? => superseded to WS by SWN
cellTypes = {'Real','ShuffTid','Uniform'};%, 'UniformTrial', 'ShuffTid'};%, 'PoissGain'}; %PoissGain not done
boolFixSpikeGroupSize = false;
dblRemOnset = 0; %remove onset period in seconds; 0.125 for sim, 0.25 for npx
runHeaderPopTimeCoding;

%% go through recs
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
		
		%get cell props
		%vecNeuronPrefOri = pi-[sLoad.sNeuron.PrefOri];
		vecNeuronType = [sLoad.sNeuron.Types]; %1=pyr,2=interneuron
		%swap types to match npx
		vecSwap = [2 1];
		vecNeuronType = vecSwap(vecNeuronType);
		intUseMaxRep = 40;
		indRemTrials = vecTrialRepetition>intUseMaxRep;
		vecOri180(indRemTrials) = [];
		vecStimIdx(indRemTrials) = [];
		vecStimOnTime(indRemTrials) = [];
		vecStimOffTime(indRemTrials) = [];
		vecOrientation(indRemTrials) = [];
		dblMinHz = 0;
	elseif strcmp(strRunType,'Npx') || strcmp(strRunType,'SWN')
		%prep
		runRecPrepNpx;
		strThisRec = strRec;
		strDataPathT0 = strTargetDataPath;
		dblMinHz = 90;
		
		
		%% layers
		sArea1Neurons = sUseNeuron(indArea1Neurons);
		cellAreas = {sArea1Neurons.Area};
		vecCorticalLayer = cellfun(@(x) str2double(x(regexp(x,'layer.*')+6)),cellAreas);
		vecDepth = [sArea1Neurons.DepthBelowIntersect];
		vecSupraGranuInfra = double(vecCorticalLayer < 4) + 2*double(vecCorticalLayer == 4) + 3*double(vecCorticalLayer > 4);
		vecSupraGranuInfra = vecSupraGranuInfra(indResp);
	else
		error impossible
	end
	
	
	%% move onset
	%remove first x ms
	vecStimOnTime = vecStimOnTime + dblRemOnset;
	
	%% load prepped data
	strTarget = fullpath(strDataPathT0,[sprintf('T0Data_%s%s%s%s%s',strThisRec,strRunType,strRunStim,'Real') '.mat']);
	if ~exist(strTarget,'file')
		fprintf('Prepped T0 file did not exist for %s; skipping...\n',strThisRec);
		continue;
	else
		fprintf('Running %s... [%s]\n',strRec,getTime);
	end
	
	%check fr
	sSource = load(strTarget);
	cellSpikeTimes = sSource.cellSpikeTimes;
	intNumN = numel(cellSpikeTimes);
	cellSpikeTimesPerCellPerTrial = cell(intNumN,intTrialNum);
	for intN=1:intNumN
		%real
		[vecTrialPerSpike,vecTimePerSpike] = getSpikesInTrial(cellSpikeTimes{intN},vecStimOnTime,dblStimDur);
		for intTrial=1:numel(vecStimOnTime)
			vecSpikeT = vecTimePerSpike(vecTrialPerSpike==intTrial);
			cellSpikeTimesPerCellPerTrial{intN,intTrial} = vecSpikeT;
		end
	end
	matData = cellfun(@numel,cellSpikeTimesPerCellPerTrial)./dblStimDur;
	if mean(sum(matData)) < dblMinHz%90 / 50
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
	
	%types: Real, UniformTrial, ShuffTid, PoissGain
	for intType=1:numel(cellTypes)
		strType = cellTypes{intType};
		
		%% load prepped data and ifr
		sSource = load(fullpath(strDataPathT0,sprintf('T0Data_%s%s%s%s%s',strThisRec,strRunType,strRunStim,strType)));
		vecTime = sSource.vecTime;
		vecIFR = sSource.vecIFR;
		vecAllSpikeTime = sSource.vecAllSpikeTime;
		vecAllSpikeNeuron = sSource.vecAllSpikeNeuron;
		cellSpikeTimes = sSource.cellSpikeTimes;
		if numel(cellSpikeTimes) ~= intNumN
			error('neuron # mismatch!')
		end
		if isempty(vecTime),continue;end
		
		
		%% cell type (narrow/broad)
		%get cell props
		vecNeuronType = sSource.vecNeuronType(indResp); %1=interneuron/narrow,2=pyramidal/broad
		vecSupraGranuInfra = sSource.vecSupraGranuInfra;%1=supra,2=granu,3=infra
		
		%% build trial-neuron cell matrix
		cellSpikeTimesPerCellPerTrial = cell(intNumN,intTrialNum);
		for intN=1:intNumN
			%real
			cellSpikeTimes{intN}(cellSpikeTimes{intN} > vecStimOnTime(end)+dblStimDur*2)=[];
			[vecTrialPerSpike,vecTimePerSpike] = getSpikesInTrial(cellSpikeTimes{intN},vecStimOnTime,dblStimDur);
			for intTrial=1:intTrialNum
				vecSpikeT = vecTimePerSpike(vecTrialPerSpike==intTrial);
				cellSpikeTimesPerCellPerTrial{intN,intTrial} = vecSpikeT;
			end
		end
		matData = cellfun(@numel,cellSpikeTimesPerCellPerTrial)./dblStimDur;
		sTuning = getTuningCurves(matData,vecOri180,0);
		vecNeuronPrefOri = sTuning.matFittedParams(:,1);
		vecNeuronBandwidth = real(sTuning.matBandwidth);
		vecNeuronLayer = vecSupraGranuInfra;
		
		%% check rates
		if boolFixSpikeGroupSize
			intSpikeGroupSize = 10;
			strSGS = ['Fixed' num2str(intSpikeGroupSize)];
		else
			intSpikeGroupSize = ceil(mean(sum(matData))/30); %20
			strSGS = ['Var' num2str(intSpikeGroupSize)];
		end
		
		%% go through pop spikes and group into sets of 20
		fprintf('   Collecting n-spike groups and decoding for %s-%s [%s]\n',strRec,strType,getTime);
		[sSpikeGroup,matSpikeGroupData] = getSpikeGroupData(cellSpikeTimesPerCellPerTrial,intSpikeGroupSize,vecOriIdx,vecStimOnTime,vecTime,vecIFR);
		
		%delete outliers
		dblRemOutliers = 0;%0.05
		[vecSortedRate,vecSort]=sort(intSpikeGroupSize./[sSpikeGroup.Duration]');
		vecRem = vecSort([1:round(numel(vecSortedRate)*(dblRemOutliers/2)) round(numel(vecSortedRate)*(1-dblRemOutliers/2)):numel(vecSortedRate)]);
		sSpikeGroup(vecRem) = [];
		matSpikeGroupData(vecRem,:) = [];
		
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
		
		
		vecTuningPerCell = real(sTuning.vecFitT);
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
		vecSpikeGroupFractionSupragranular = nan(size(vecSpikeGroupLatency)); %fraction of supragranular cells
		vecSpikeGroupFractionGranular = nan(size(vecSpikeGroupLatency)); %fraction of granular cells
		vecSpikeGroupFractionInfragranular = nan(size(vecSpikeGroupLatency)); %fraction of infragranular cells
		
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
			
			%get layer
			vecSpikeGroupFractionSupragranular(intSpikeGroup) = mean(vecSupraGranuInfra(vecActiveCells)==1);
			vecSpikeGroupFractionGranular(intSpikeGroup) = mean(vecSupraGranuInfra(vecActiveCells)==2);
			vecSpikeGroupFractionInfragranular(intSpikeGroup) = mean(vecSupraGranuInfra(vecActiveCells)==3);
			
			%get avg pref dist to stim
			vecSpikeGroupAvgPrefDistToStim(intSpikeGroup) = mean(abs(circ_dist(2*vecNeuronPrefOri(vecActiveCells),dblStimRad*2)));
			vecSpikeGroupAvgBandwidthOfCells(intSpikeGroup) = mean(vecNeuronBandwidth(vecActiveCells));
			
			%assign to spike group
			sSpikeGroup(intSpikeGroup).AvgTuningOfCells = vecSpikeGroupAvgTuningOfCells(intSpikeGroup);
			sSpikeGroup(intSpikeGroup).AvgRateOfCells = vecSpikeGroupAvgRateOfCells(intSpikeGroup);
			sSpikeGroup(intSpikeGroup).AvgNormRespOfCellsToStim = vecSpikeGroupAvgNormRespOfCellsToStim(intSpikeGroup);
			sSpikeGroup(intSpikeGroup).NumOfCells = vecSpikeGroupNumOfCells(intSpikeGroup);
			
			sSpikeGroup(intSpikeGroup).FractionInterneurons = vecSpikeGroupFractionInterneurons(intSpikeGroup);
			sSpikeGroup(intSpikeGroup).FractionSupra = vecSpikeGroupFractionSupragranular(intSpikeGroup);
			sSpikeGroup(intSpikeGroup).FractionGranu = vecSpikeGroupFractionGranular(intSpikeGroup);
			sSpikeGroup(intSpikeGroup).FractionInfra = vecSpikeGroupFractionInfragranular(intSpikeGroup);
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
		
		%% plot cell tune/ifr
		%make example
		intPlotTrial = 10;
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
		
		%match local groups to global groups
		vecGlobalGroups = find(vecSpikeGroupTrialNumber==intPlotTrial);
		vecGlobalLatencies = vecSpikeGroupLatency(vecGlobalGroups);
		
		%assign to spike group
		intAssignGroups = floor(numel(vecSpikeTimes)/intSpikeGroupSize);
		vecGroupStart = nan(1,intAssignGroups);
		vecGroupStop = nan(1,intAssignGroups);
		vecGroupLatency = nan(1,intAssignGroups);
		for intGroup=1:intAssignGroups
			%get data
			intEndSpike = intGroup*intSpikeGroupSize;
			vecUseSpikes = (intEndSpike-intSpikeGroupSize+1):intEndSpike;
			vecSpikeGroup(vecUseSpikes) = intGroup;
			vecGroupStart(intGroup) = vecSpikeTimes(vecUseSpikes(1));
			vecGroupStop(intGroup) = vecSpikeTimes(vecUseSpikes(end));
			vecGroupLatency(intGroup) = mean(vecSpikeTimes(vecUseSpikes));
		end
		vecSpikeGroup(isnan(vecSpikeGroup))=0;
		
		%remove groups
		indRem = ~ismember(vecGroupLatency,vecGlobalLatencies);
		vecGroupStart(indRem) = [];
		vecGroupStop(indRem) = [];
		vecGroupLatency(indRem) = [];
		intAssignGroups = numel(vecGroupLatency);
		
		%calc local ifr
		[vecT,vecR]=getIFR(vecSpikeTimes,0,1);
		
		if 0
			%% make uncolored example plot
			%get entry in global time
			vecTime = sSource.vecTime;
			vecIFR = sSource.vecIFR;
			dblGlobalTimeStart = vecStimOnTime(intPlotTrial);
			dblGlobalTimeStop = dblGlobalTimeStart + dblStimDur;
			intStartEntry = [];
			intStopEntry = [];
			for intGlobalSpikeEntry=1:numel(vecTime)
				if isempty(intStartEntry) && vecTime(intGlobalSpikeEntry) >= dblGlobalTimeStart
					intStartEntry = intGlobalSpikeEntry;
				end
				if vecTime(intGlobalSpikeEntry) >= dblGlobalTimeStop
					intStopEntry = intGlobalSpikeEntry;
					break;
				end
			end
			vecLocalTime = vecTime(intStartEntry:intStopEntry);
			vecLocalIFR = vecIFR(intStartEntry:intStopEntry);
			figure
			
			subplot(2,3,1);
			vecRandLoc = randperm(intNumN);
			scatter(dblRemOnset+vecSpikeTimes,vecSpikeNeuron,[],'k','o','filled')
			xlim([0 max(get(gca,'xlim'))]);
			xlabel('Time after stim onset (s)');
			ylabel('Neuron ID');
			
			vecTime = sSource.vecTime;
			vecIFR = sSource.vecIFR;
			
			
% 			subplot(2,3,2);
% 			plot(dblRemOnset+vecLocalTime,vecLocalIFR,'color',[0.5 0.5 0.5]);
% 			hold on
% 			title(sprintf('%s, %s, %s trial %d',strThisRec,strType,strOnset,intPlotTrial),'interpreter','none');
% 			xlabel('Time after stim onset (s)');
% 			ylabel('Pop IFR');
% 			ylim([0 max(get(gca,'ylim'))]);
			
			subplot(2,3,4);
			plot(dblRemOnset+vecT,vecR,'color',0*[0.5 0.5 0.5]);
			hold on
			title(sprintf('%s, %s, %s trial %d',strThisRec,strType,strOnset,intPlotTrial),'interpreter','none');
			xlabel('Time after stim onset (s)');
			ylabel('Pop IFR');
			ylim([0 max(get(gca,'ylim'))]);
			fixfig
		end
		
		
		% plot example trial
		strLabelX = 'avg pop rate (Hz)';
		figure;maxfig;
		
		%make plot
		h1=subplot(3,4,1);
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
		%[vecSortedIfr,vecSort]=sort(vecSpikeGroupAvgIFR);
		[vecSortedIfr,vecSort]=sort(intSpikeGroupSize./vecSpikeGroupDuration);
		vecSortedCorr = vecSpikeGroupCorrect(vecSort);
		if 0 %tuning
			strProp = 'Avg cell tuning (t-stat)';
			vecLimY = [2 12];
			vecLimX = [0 4000];
			vecSortedProp = vecSpikeGroupAvgTuningOfCells(vecSort);
			vecProp = vecSpikeGroupAvgTuningOfCells;
		else % # cells
			strProp = '# of contrib cells';
			vecLimY = [14 17];
			vecLimX = [0 4000];
			vecSortedProp = vecSpikeGroupNumOfCells(vecSort);
			vecProp = vecSpikeGroupNumOfCells;
		end
		strLabelY = strProp;

		vecSortedLat = vecSpikeGroupLatency(vecSort);
		%remove long durs
		%indRem = vecSortedDur > 0.05;
		%indRem = vecSortedLat < 0.15 | vecSortedDur > 0.05;
		indRem = [];
		vecSortedIfr(indRem) = [];
		vecSortedCorr(indRem) = [];
		vecSortedProp(indRem) = [];
		vecSortedLat(indRem) = [];
		
		%calculate confidence per bin of equal size
		intQuantileNum = 10;
		[vecMeanIfr,vecSemIfr,vecMeanProp,vecSemProp,vecQuantile]=getQuantiles(vecSortedIfr,vecSortedProp,intQuantileNum);
		
		%calculate accuracy per bin of equal size
		intSperBin = floor(numel(vecSortedIfr)/intQuantileNum);
		vecMeanCorr = nan(1,intQuantileNum);
		matCiCorr = nan(2,intQuantileNum);
		cellValsCorr = cell(1,intQuantileNum);
		vecSampleGroup = zeros(size(vecSortedIfr));
		for intBin=1:intQuantileNum
			intEndS = intSperBin*intBin;
			vecSamples = (intEndS-intSperBin+1):intEndS;
			
			[phat,pci] = binofit(sum(vecSortedCorr(vecSamples)),intSperBin);
			vecMeanCorr(intBin) = phat;
			matCiCorr(:,intBin) = pci;
			cellValsCorr{intBin} = vecSortedCorr(vecSamples);
			vecSampleGroup(vecSamples) = intBin;
		end
		[r,p]=corr(vecSortedIfr,vecSortedProp);
		[r2,p2]=corr(vecSortedIfr,vecSortedCorr);
		
		
		[bandwidth,density,X,Y]=kde2d([vecSortedIfr vecSortedProp],512);
		h=subplot(2,4,2);
		imagesc(X(1,1:end),Y(1:end,1),density);colormap(h,flipud(gray(1024)));axis xy;
		xlabel(sprintf('%d-spike block %s',intSpikeGroupSize,strLabelX));
		ylabel(strLabelY);
		xlim([0 2000]);
		%ylim(vecLimY);
		
		subplot(2,4,6)
		hold on
		h=scatter(vecSortedIfr,vecSortedProp,100,[0.2 0.2 0.2],'.');
		h.MarkerFaceAlpha = 0.1;
		h.MarkerEdgeAlpha = 0.1;
		h2=scatter(vecSpikeGroupAvgIFR(vecGlobalGroups),vecProp(vecGlobalGroups),200,matCol(2:end,:),'.');
		hold off
		xlabel(sprintf('%d-spike block %s',intSpikeGroupSize,strLabelX));
		ylabel(strLabelY);
		xlim(vecLimX);
		%ylim(vecLimY);
		title(sprintf('r=%.3f, p=%.3e',r,p));
		
		subplot(2,4,8)
		vecBinsE = 0:25:max(vecSortedIfr);
		vecBinsC = vecBinsE(2:end)-diff(vecBinsE(1:2))/2;
		[vecCounts,vecMeans,vecSDs] = makeBins(vecSortedIfr,vecSortedCorr,vecBinsE);
		hold on
		plot(vecBinsC,vecMeans,'color',lines(1));
		plot(vecBinsC,vecMeans+vecSDs./sqrt(vecCounts),'color',lines(1));
		plot(vecBinsC,vecMeans-vecSDs./sqrt(vecCounts),'color',lines(1));
		hold off
		xlabel(sprintf('%d-spike block %s',intSpikeGroupSize,strLabelX));
		ylabel('Decoder accuracy');
		title(sprintf('r=%.3f, p=%.3f',r2,p2));
		
		subplot(2,4,4)
		[vecCounts,vecMeans,vecSDs] = makeBins(vecSortedIfr,vecSortedProp,vecBinsE);
		indPlot = vecCounts>25;
		hold on
		plot(vecBinsC(indPlot),vecMeans(indPlot),'color',lines(1));
		plot(vecBinsC(indPlot),vecMeans(indPlot)+vecSDs(indPlot)./sqrt(vecCounts(indPlot)),'color',lines(1));
		plot(vecBinsC(indPlot),vecMeans(indPlot)-vecSDs(indPlot)./sqrt(vecCounts(indPlot)),'color',lines(1));
		hold off
		xlabel(sprintf('%d-spike block %s',intSpikeGroupSize,strLabelX));
		ylabel(strLabelY);
		
		xlim([0 2500]);
		%ylim([5.5 6.5]);
		
		%anova corr
		indRemS = vecSampleGroup==0;
		vecSampleGroup(indRemS) = [];
		vecSampleCorr = vecSortedCorr(~indRemS);
		[pA,tbl,stats] = anova1(vecSampleCorr,vecSampleGroup,'off');
		[c,~,~,gnames] = multcompare(stats,'CType','bonferroni','display','off');
		subplot(2,4,7)
		hold on
		errorbar(vecMeanIfr,vecMeanCorr*100,(vecMeanCorr-matCiCorr(1,:))*100,(vecMeanCorr-matCiCorr(2,:))*100,vecSemIfr,vecSemIfr,'color',lines(1));
		hold off
		xlabel(sprintf('%d-spike block %s',intSpikeGroupSize,strLabelX));
		ylabel('Decoder accuracy');
		title(sprintf('Deciles, ANOVA, p=%.1e',pA));
		
		
		%anova conf
		indRemS = vecSampleGroup==0;
		vecSampleGroup(indRemS) = [];
		vecSampleConf = vecSortedProp(~indRemS);
		[pA2,table2,stats2] = anova1(vecSampleConf,vecSampleGroup,'off');
		[c2,~,~,gnames2] = multcompare(stats,'CType','bonferroni','display','off');
		
		subplot(2,4,3)
		hold on
		
		%errorbar(vecMeanIfr,vecMeanTune,vecSemTune,vecSemTune,vecSemIfr,vecSemIfr,'color',lines(1));
		errorbar(vecMeanIfr,vecMeanProp,vecSemProp,vecSemProp,vecSemIfr,vecSemIfr,'color',lines(1));
		hold off
		xlabel(sprintf('%d-spike block %s',intSpikeGroupSize,strLabelX));
		ylabel(strLabelY);
		title(sprintf('Deciles, ANOVA, p=%.1e',pA2));
		
		fixfig;
			
		%%
		export_fig(fullpath(strFigurePathSR,sprintf('A1c_Example_%s_%s_SGS%s%s.tif',strThisRec,strType,strSGS,strOnset)));
		export_fig(fullpath(strFigurePathSR,sprintf('A1c_Example_%s_%s_SGS%s%s.pdf',strThisRec,strType,strSGS,strOnset)));
		
		%% plot multiple timescales
		vecStimOnTime = sSource.vecStimOnTime;
		intTrial=317;
		dblStimT = vecStimOnTime(intTrial);
		dblCenterTime = dblStimT+0.4;
		
		vecAllSpikeT = sSource.vecAllSpikeTime;
		dblPlotWindow = 50;
		figure
		vecTimescales = 10.^(1:-1:-2);
		intPlotNum = numel(vecTimescales)+1;
		vecRateBins = 0:100:2000;
		for intScaleIdx=1:numel(vecTimescales)
			dblScale = vecTimescales(intScaleIdx);
			subplot(intPlotNum,2,intScaleIdx*2-1);
			
			dblStartT = dblCenterTime-dblScale*dblPlotWindow;
			dblStopT = dblCenterTime+dblScale*dblPlotWindow;
			vecBinE = dblStartT:dblScale:dblStopT;
			vecSpikeRate = histcounts(vecAllSpikeT,vecBinE)./dblScale;
			vecBinC = vecBinE(2:end)-dblScale/2;
			
			indPlot = vecT>dblStartT & vecT<dblStopT;
			stairs(vecBinC,vecSpikeRate);
			hold on
			%plot([dblStimT dblStimT],[0 2000])
			xlim([dblStartT dblStopT]);
			ylim([0 2000]);
			
			dblNextStartT = dblCenterTime-(dblScale/10)*dblPlotWindow;
			dblNextStopT = dblCenterTime+(dblScale/10)*dblPlotWindow;
			plot([dblNextStartT dblNextStartT],[0 2000],'color',[0.5 0.5 0.5])
			plot([dblNextStopT dblNextStopT],[0 2000],'color',[0.5 0.5 0.5])
			
			
			subplot(intPlotNum,4,intScaleIdx*4-1);
			vecCounts = histcounts(vecSpikeRate,vecRateBins);
			stairs(vecRateBins(2:end)-diff(vecRateBins(1:2))/2,vecCounts);
			xlim([0 2000]);
			ylim([0 40]);
			title(sprintf('mu=%.3f, sd=%.3f, cv=%.3f',mean(vecSpikeRate),std(vecSpikeRate),std(vecSpikeRate)/mean(vecSpikeRate)));
		end
		
		subplot(intPlotNum,2,intPlotNum*2-1);
		vecIFR_T = vecAllSpikeT(vecAllSpikeT>dblStartT & vecAllSpikeT<dblStopT);
		[vecT,vecIFR]=getIFR(vecIFR_T,dblStartT,dblStopT-dblStartT);
		plot(vecT+dblStartT,vecIFR);
		xlim([dblStartT dblStopT]);
		%ylim([0 2000]);
		
		%%
		export_fig(fullpath(strFigurePathSR,sprintf('A1c_Example2_%s_%s_%s.tif',strThisRec,strType,strOnset)));
		export_fig(fullpath(strFigurePathSR,sprintf('A1c_Example2_%s_%s_%s.pdf',strThisRec,strType,strOnset)));
		return
		
		%%
		vecT=10.^(-4:0.1:1);
		vecLinCV = 0.001*(1./vecT);
		
		plot(vecT,vecLinCV);
		set(gca,'xscale','log','yscale','log')
		
		%% plot quantiles
		vecPopRate = sum(matData,1);
		figure
		subplot(2,3,1)
		vecBinE = 0:50:1200;
		vecBinC = vecBinE(2:end)-diff(vecBinE(1:2))/2;
		vecPopR = histcounts(vecPopRate,vecBinE);
		plot(vecBinC,vecPopR);hold on
		%plot quantiles
		[vecSortedR,vecReorder]=sort(vecPopRate);
		dblQsize = 960/5;
		vecQedges = dblQsize:dblQsize:959;
		for i =1:numel(vecQedges)
			plot(vecSortedR(vecQedges(i))*[1 1],[0 150]);
		end
		
		subplot(2,3,2);
		colormap(redbluepurple);
		[coeff,score,latent,tsquared,explained,mu] = pca(matData);
		h=scatter(coeff(:,1),coeff(:,2),[],vecPopRate,'marker','o');
		h.MarkerFaceColor='flat';
		h.MarkerEdgeColor='none';
		colorbar
		
		%% 35 54
		matCol=redbluepurple(5);
		n1 = 83;%35
		n2 = 86;

			cla;hold on;
			vec1 = matData(n1,:);%+rand(1,960)-0.5;
			vec2 = matData(n2,:);%+rand(1,960)-0.5;
			vecTot =vec1+vec2;
			[vecMeanX,vecSemX,vecMeanY,vecSemY,vecQuantileAssignment] = getQuantiles(vecTot,vecTot,5);
			for j=1:5
				scatter(vec1(vecQuantileAssignment==j),vec2(vecQuantileAssignment==j),[],matCol(j,:),'filled')
			end
			title(sprintf('%d',i))
		%xlim([0 100]);
		%ylim([0 100]);
		%pause

	end
	close all;
end
toc
