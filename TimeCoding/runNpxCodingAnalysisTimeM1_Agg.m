%% overview
%show spiking transfer from v1 to v2 occurs primarily during peaks in model


%% set parameters
cellDataTypes = {'Npx','Sim','ABI','SWN'};%topo, model, allen, nora
intRunDataType = 2;
strRunStim = 'DG';%DG or NM? => superseded to WS by SWN
cellTypes = {'Real','ShuffTid','Uniform'};%, 'UniformTrial', 'ShuffTid'};%, 'PoissGain'}; %PoissGain not done
boolFixSpikeGroupSize = false;
dblRemOnset = 0; %remove onset period in seconds; 0.125 for sim, 0.25 for npx
runHeaderPopTimeCoding;

%% load V2
strFileV2 = 'SimDG18_V2pop.mat';
sV2 = load(fullpath(strTargetDataPath,strFileV2));

%% go through recs
tic
for intRec=1:intRecNum %19 || weird: 11
	%% prep sim
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
		intSpikeGroupNum = numel(sSpikeGroup);
		if intSpikeGroupNum < intTrialNum,continue;end
		
		%% go through spike groups
		%show spiking transfer from v1 to v2 occurs primarily during peaks in model and specifically
		%in connected V2 cells(?)
		
		% what is the hypothesis exactly?
		
		% how do we answer this?
		
		% > decode stimuli in V2 time bins (or spike blocks?) following spike blocks in V1
		% is decoding accuracy higher when rate was higher in V1?
		
		%what are they connected to?
		vecCellIdsV1 = sLoad.vecOrigIds;
		vecCellIdsV2 = 1201:2400;
		matV1ToV2 = sV2.matConn(vecCellIdsV1,vecCellIdsV2);
		vecInputV2 = sum(sV2.matConn(vecCellIdsV1,vecCellIdsV2));
		cellSpikesV2 = sV2.cellSpikeTimes;
		vecNeuronTypeV2 = [sV2.sNeuron.Types];
		
		%build trial matrix
		cellSpikeTimesPerCellPerTrialV2 = cell(numel(cellSpikesV2),intTrialNum);
		for intN=1:numel(cellSpikesV2)
			%real
			cellSpikesV2{intN}(cellSpikesV2{intN} > vecStimOnTime(end)+dblStimDur*2)=[];
			[vecTrialPerSpike,vecTimePerSpike] = getSpikesInTrial(cellSpikesV2{intN},vecStimOnTime,dblStimDur);
			for intTrial=1:intTrialNum
				vecSpikeT = vecTimePerSpike(vecTrialPerSpike==intTrial);
				cellSpikeTimesPerCellPerTrialV2{intN,intTrial} = vecSpikeT;
			end
		end
		
		%perform analysis in spike blocks and give it the property of average ifr in V1
		%also give it the property of average v1-v2 input strength per neuron
		
		%select x# of most connected cells
		[dummy,vecSortedV2]=sort(vecInputV2,'descend');
		for intConnType=1:2
			if intConnType == 1
				%most
				strConn = 'MostConn';
				vecUseCells = vecSortedV2(1:intNumN);
			else
				%least
				strConn = 'LeastConn';
				vecUseCells = vecSortedV2((end-intNumN+1):end);
			end
			cellST_V2 = cellSpikeTimesPerCellPerTrialV2(vecUseCells,:);
			matDataV2 = cellfun(@numel,cellST_V2)./dblStimDur;
			sTuningV2 = getTuningCurves(matDataV2,vecOri180,0);
			vecNeuronType = vecNeuronTypeV2(vecUseCells);
			
			%% check rates
			if boolFixSpikeGroupSize
				intSpikeGroupSizeV2 = 10;
				strSGS = ['Fixed' num2str(intSpikeGroupSizeV2)];
			else
				intSpikeGroupSizeV2 = ceil(mean(sum(matDataV2))/30); %20
				strSGS = ['Var' num2str(intSpikeGroupSizeV2)];
			end
			
			%get spike blocks
			[sSpikeGroup,matSpikeGroupData] = getSpikeGroupData(cellST_V2,intSpikeGroupSizeV2,vecOriIdx,vecStimOnTime,vecTime,vecIFR);
			%ifr is v1 IFR
			
			vecIFR_V1 = [sSpikeGroup.AvgRate]';
			vecDur_V2 = [sSpikeGroup.Duration]';
			vecRate_V2 = 1./(vecDur_V2*intSpikeGroupSizeV2);
			vecConf_V2 = [sSpikeGroup.Confidence]';
			vecCorr_V2 = [sSpikeGroup.Correct]';
			[rRC,pRC]=corr(vecIFR_V1,vecConf_V2);
			[rRD,pRD]=corr(vecIFR_V1,vecRate_V2);
			[rCD,pCD]=corr(vecConf_V2,vecRate_V2);
			
			%% add additional properties
			vecTuningPerCell = real(sTuningV2.vecFitT);
			vecRatePerCell = mean(sTuningV2.matMeanResp,2);
			vecMaxOriResp = max(sTuningV2.matMeanResp,[],2);
			matNormRespPerCellPerOri = sTuningV2.matMeanResp ./ vecMaxOriResp;
			vecSpikeGroupAvgTuningOfCells = nan(size(vecIFR_V1)); %t-statistic
			vecSpikeGroupAvgBandwidthOfCells = nan(size(vecIFR_V1)); %t-statistic
			vecSpikeGroupAvgRateOfCells = nan(size(vecIFR_V1)); %hz
			vecSpikeGroupAvgNormRespOfCellsToStim = nan(size(vecIFR_V1)); %normalized response to stim ori
			vecSpikeGroupNumOfCells = nan(size(vecIFR_V1)); %how many cells participate?
			vecSpikeGroupFractionInterneurons = nan(size(vecIFR_V1)); %fraction of interneurons
			vecSpikeGroupAvgPrefDistToStim = nan(size(vecIFR_V1)); %average distance to stim ori
			intSpikeGroupNum = numel(sSpikeGroup);
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
			
			%% save data
			save(fullpath(strTargetDataPath,sprintf('M1Data_%s_%s_SGS%s%s%s.mat',strThisRec,strType,strSGS,strConn,strOnset)),...
				'vecStimOnTime',...
				'vecStimOffTime',...
				'vecOrientation',...
				'vecOri180',...
				'intSpikeGroupSize',...
				'intSpikeGroupSizeV2',...
				'strRec',...
				'strType',...
				'strConn',...
				'sTuningV2',...'rSpike_IFR_Rate',...
				...'pSpike_IFR_Rate',...
				...'rSpike_IFR_Tune',...
				...'pSpike_IFR_Tune',...
				'sSpikeGroup');
			
		end
		close all;
	end
end
toc
