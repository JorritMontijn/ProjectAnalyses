%% overview
%show spiking transfer from v1 to v2 occurs primarily during peaks in model


%% set parameters
cellDataTypes = {'Npx','Sim','ABI','SWN'};%topo, model, allen, nora
intRunDataType = 2;
strRunStim = 'DG';%DG or NM? => superseded to WS by SWN
cellTypes = {'Real'};%{'Real','ShuffTid','Uniform'};%, 'UniformTrial', 'ShuffTid'};%, 'PoissGain'}; %PoissGain not done
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
		
		
		%% cell type (narrow/broad)
		%get cell props
		vecNeuronType = sSource.vecNeuronType; %1=interneuron/narrow,2=pyramidal/broad
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
		intSpikeGroupNum = numel(sSpikeGroup);
		if intSpikeGroupNum < intTrialNum,continue;end
		
		%% go through spike groups
		%show spiking transfer from v1 to v2 occurs primarily during peaks in model and specifically
		%in connected V2 cells(?)
		
		% what is the hypothesis exactly?
		
		% how do we answer this?
		
		error
		
	end
	close all;
end
toc
