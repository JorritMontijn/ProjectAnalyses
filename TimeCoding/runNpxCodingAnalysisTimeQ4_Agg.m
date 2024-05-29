%% overview
%what is better model for pop activity? multivariate gaussian noise or fixed sd/mu over a
%gain-scaling axis?

%% set parameters
cellDataTypes = {'Npx','Sim','ABI','SWN'};%topo, model, allen, nora
intRunDataType = 1;
strRunStim = 'DG';%DG or NM? => superseded to WS by SWN
dblRemOnset = 0; %remove onset period in seconds; 0.125 for sim, 0.25 for npx
cellTypes = {'Real','Poiss','ShuffTid','RandTid','RandTxClass','Uniform'};%{'Real','Poiss','ShuffTid','RandTid','RandTxClass','Uniform'};
runHeaderPopTimeCoding;
boolMakeFigs = false;
vecTimescales = [1e-2 5e-2 1e-1 5e-1 1e-0]; %timescales for rapid flucts, noise corrs, tuning
vecJitter = [0 1e3];%[0 (2.^(-9:10))];
intPopSize = inf; %24 (smallest pop of all recs) or inf (uses full pop for each rec)

%% go through recs
tic
for intRec=1:intRecNum %19 || weird: 11
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
		
		%% move onset
		%remove first x ms
		vecOrigStimOnTime = vecStimOnTime + dblRemOnset;
		
		%get data matrix
		[matData,indTuned,cellSpikeTimes,sOut,cellSpikeTimesPerCellPerTrial,indResp] = ...
			SimPrepData(cellSpikeTimesRaw,vecOrigStimOnTime,vecStimOffTime,vecOrientation);
		
		%get cell props
		%vecNeuronPrefOri = pi-[sLoad.sNeuron.PrefOri];
		vecNeuronType = [sLoad.sNeuron.Types]; %1=pyr,2=interneuron
		
	elseif strcmp(strRunType,'Npx') || strcmp(strRunType,'SWN')
		%prep
		runRecPrepNpx;
		strThisRec = strRec;
		strDataPathT0 = strTargetDataPath;
		
		%% move onset
		%remove first x ms
		vecOrigStimOnTime = vecStimOnTime + dblRemOnset;
		
		%get data matrix
		[matData,indTuned,cellSpikeTimes,sTuning24,cellSpikeTimesPerCellPerTrial] = ...
			NpxPrepData(cellSpikeTimesRaw,vecOrigStimOnTime,vecStimOffTime,vecOrientation);
		
		%get cell props
		vecNeuronType = ones(size(indTuned)); %1=pyr,2=interneuron
		%narrow vs broad not done
	else
		error impossible
	end
	
	%get ori vars
	intTunedN = sum(indTuned);
	intRespN = size(matData,1);
	intNumN = numel(cellSpikeTimes);
	intTrialNum = numel(vecOrigStimOnTime);
	vecOri180 = mod(vecOrientation,180);
	[vecOriIdx,vecUnique,vecPriorDistribution,cellSelect,vecRepetition] = val2idx(vecOri180);
	sTuning = getTuningCurves(matData,vecOri180,0);
	vecNeuronPrefOri = sTuning.matFittedParams(:,1);
	vecNeuronBandwidth = real(sTuning.matBandwidth);
	
	intOriNum = numel(vecUnique);
	intRepNum = min(vecPriorDistribution);
	dblStimDur = median(vecStimOffTime - vecOrigStimOnTime);
	if mean(sum(matData)) < 90%90 / 50
		fprintf('Avg # of spikes per trial was %.1f for %s; skipping...\n',mean(sum(matData)),strThisRec);
		continue;
	else
		fprintf('Running %s... [%s]\n',strThisRec,getTime);
	end
	%% go through types
	clear sAggData;
	vecRunTypes = 1:numel(cellTypes);
	for intType=vecRunTypes
		strType = cellTypes{intType};
		
		%% load prepped data and ifr
		sSource = load(fullpath(strDataPathT0,sprintf('T0Data_%s%s%s%s%s',strThisRec,strRunType,strRunStim,strType)));
		cellOrigSpikeTimes = sSource.cellSpikeTimes;
		vecPseudoStartT = vecOrigStimOnTime;
		
		vecTrialT = sSource.vecStimOnTime;
		vecTrialT_end = vecTrialT+1;
		vecEdges = sort(cat(2,vecTrialT,vecTrialT_end));
		vecSpikeCounts = histcounts(sSource.vecAllSpikeTime,vecEdges);
		vecSpikesPerTrial = vecSpikeCounts(1:2:end);
		
		[a,b,c,d,e,cellValsX]=getQuantiles(vecSpikesPerTrial,vecSpikesPerTrial,7);
		vecMu = cellfun(@mean,cellValsX);
		vecSd = cellfun(@std,cellValsX);
		plot(vecMu,vecSd)
		return
	end
	%% save agg data
	save(fullpath(strTargetDataPath,sprintf('Q3Data%s_%s.mat',strThisRec,strOnset)),...
		'sAggData');
end
toc
