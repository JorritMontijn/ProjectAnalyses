
%% get matching recording data
close all;
strRec = sAggStim(intRec).Exp;
strRecOrig = strRec;
sThisRec = sAggStim(strcmpi(strRec,{sAggStim(:).Exp}));
sThisSource = sAggSources(strcmpi(strRec,{sAggSources(:).Exp}));
if strcmp(strRunStim,'DG') || strcmp(strRunType,'SWN')
	if strcmp(strRunStim,'DG')
		%prep grating data
		[sUseNeuron,vecStimOnTime,vecStimOffTime,vecOrientation] = NpxPrepGrating(sAggNeuron,sThisRec,cellUseAreas);
	elseif strcmp(strRunType,'SWN')
		%pretend that whisk stim and no whisk stim are two different orientations
		[sUseNeuron,vecStimOnTime,vecStimOffTime,vecOrientation] = NpxPrepWhisking(sAggNeuron,sThisRec,cellUseAreas);
	end
	[vecStimIdx,vecUnique,vecRepNum,cellSelect,vecTrialRepetition] = val2idx(vecOrientation);
	intTrialNum = numel(vecStimOnTime);
	intOriNum = numel(unique(vecOrientation));
	intRepNum = intTrialNum/intOriNum;
	intNeuronsInArea = numel(sUseNeuron);
	intNeuronNum = intNeuronsInArea;
	if intNeuronsInArea==0,return;end
	intOrigTrialNum = intTrialNum;
	
	%% get neurons in this area
	indArea1Neurons = contains({sUseNeuron.Area},strArea,'IgnoreCase',true);
	intNeuronsInArea = sum(indArea1Neurons);
	intNeuronNum = intNeuronsInArea;
	if intNeuronsInArea==0,return;end

	%% prep data
	sArea1Neurons = sUseNeuron(indArea1Neurons);
	cellSpikeTimesRaw = {sArea1Neurons.SpikeTimes};
	vecOri180 = mod(vecOrientation,180)*2;
	vecStimIdx = vecOri180;
	[matData,indTuned,cellSpikeTimes,sOut,cellSpikeTimesPerCellPerTrial,vecStimOnStitched,vecNonStat,dblBC,dblMaxDevFrac,indResp] = ...
		NpxPrepData(cellSpikeTimesRaw,vecStimOnTime,vecStimOffTime,vecOrientation);
	intTunedN = sum(indTuned);
	intRespN = size(matData,1);
	dblStimDur = roundi(median(vecStimOffTime - vecStimOnTime),1,'ceil');
	dblPreTime = 0;%0.3;
	dblPostTime = 0;%0.3;
	dblMaxDur = dblStimDur+dblPreTime+dblPostTime;
	intTrialNum = numel(vecStimOnTime);
elseif strcmp(strRunStim,'NM')
	%to do
	%prep move data
	[sUseNeuron,vecStimOnTimeBinned,vecStimOffTimeBinned,vecStimIdx,structStim] = NpxPrepMovies(sAggNeuron,sThisRec,cellUseAreas);
	intNeuronsInArea = numel(sUseNeuron);
	intNeuronNum = intNeuronsInArea;
	if intNeuronsInArea==0,return;end
	vecStimOnTime = structStim.vecStimOnTime;
	vecOrigStimOnTime = vecStimOnTime;
	vecStimIdx = 1:numel(vecStimOnTime);
	dblDur = median(diff(vecStimOnTime));
	vecStimOffTime = vecStimOnTime+dblDur;
	vecOrigStimOffTime = vecStimOffTime;
	[vecStimIdx,vecUnique,vecRepNum,cellSelect,vecTrialRepetition] = val2idx(vecStimIdx);
	intTrialNum = numel(vecStimOnTime);
	intRepNum = numel(vecStimOnTime);
	
	%% get neurons in this area
	indArea1Neurons = contains({sUseNeuron.Area},strArea,'IgnoreCase',true);
	intNeuronsInArea = sum(indArea1Neurons);
	intNeuronNum = intNeuronsInArea;
	if intNeuronsInArea==0,return;end

	%% prep data
	sArea1Neurons = sUseNeuron(indArea1Neurons);
	cellSpikeTimesRaw = {sArea1Neurons.SpikeTimes};
	
	dblMinRate = 0.1;
	
	%get dur
	cellSpikeTimes = cellSpikeTimesRaw;
	dblDur = median(vecStimOffTime-vecStimOnTime);
	matRawData = getSpikeCounts(cellSpikeTimes,vecStimOnTime,dblDur)./dblDur;
	indResp = sum(matRawData,2)'>(size(matRawData,2)/dblDur)*dblMinRate;
	matData = matRawData(indResp,:);
	
	%% remove non-responsive cells
	cellSpikeTimes(~indResp)=[];
	indTuned = indResp;
	intTrialNum = numel(vecStimIdx);
	intRespN = size(matData,1);
	intStimNum = numel(vecUnique);
	intRepNum = min(vecRepNum);
	
	%% check non-stationarity
	%get spikes per trial per neuron
	cellSpikeTimesPerCellPerTrial = cell(intRespN,intTrialNum);
	vecNonStat = nan(1,intRespN);
	boolDiscardEdges = true;
	for intN=1:intRespN
		% build pseudo data, stitching stimulus periods
		[vecPseudoSpikeTimes,vecPseudoEventT] = getPseudoSpikeVectors(cellSpikeTimes{intN},vecStimOnTime,dblDur,boolDiscardEdges);
		
		%real
		[vecTrialPerSpike,vecTimePerSpike] = getSpikesInTrial(vecPseudoSpikeTimes,vecPseudoEventT,dblDur);
		for intTrial=1:intTrialNum
			vecSpikeT = vecTimePerSpike(vecTrialPerSpike==intTrial);
			cellSpikeTimesPerCellPerTrial{intN,intTrial} = vecSpikeT;
		end
		
		%calc non-stationarity
		vecSortedSpikeTimes = sort(vecPseudoSpikeTimes,'ascend') - min(vecPseudoSpikeTimes);
		dblAUC = sum(vecSortedSpikeTimes);
		dblLinAUC = (max(vecSortedSpikeTimes) * numel(vecSortedSpikeTimes) ) / 2;
		vecNonStat(intN) = (dblAUC - dblLinAUC) / dblLinAUC;
	end
	vecStimOnStitched = vecPseudoEventT;
	matDataZ = zscore(log(1+matData),[],2);
	vecMeanZ = mean(matDataZ,1);
	vecFilt = normpdf(-4:4,0,1)/sum(normpdf(-4:4,0,1));
	vecFiltM = imfilt(vecMeanZ,vecFilt);
	
	%calc metrics
	[h,pKS,ksstat,cv] = kstest(vecFiltM);
	[BF, dblBC] = bimodalitycoeff(vecFiltM);
	dblMaxDevFrac = max(abs(vecFiltM));
	
	intTunedN = sum(indTuned);
	intRespN = size(matData,1);
	dblStimDur = roundi(median(vecStimOffTime - vecStimOnTime),1,'ceil');
	dblPreTime = 0;%0.3;
	dblPostTime = 0;%0.3;
	dblMaxDur = dblStimDur+dblPreTime+dblPostTime;
end

%% define quantiles and remove zero-variance neurons
%constants
[vecStimIdx,vecUnique,vecPriorDistribution,cellSelect,vecRepetition] = val2idx(vecStimIdx);
intStimNum = numel(vecUnique);
dblLambda = 1;%1
intTypeCV = 2;
dblUseStartT = 0;
dblUseMaxDur = dblMaxDur-dblUseStartT;
intUseMax = inf;
intRepNum = min(vecPriorDistribution);

%remove zero-variance neurons
matSpikeCounts_pre = cellfun(@(x) sum(x>dblUseStartT),cellSpikeTimesPerCellPerTrial);
matMeanRate_pre = matSpikeCounts_pre./dblUseMaxDur;
vecPopRate_pre = sum(matMeanRate_pre,1);

intQuantiles = 5;
vecStartTrials = round(intRepNum*linspace(1/intRepNum,(1+1/intRepNum),intQuantiles+1));
vecStartTrials(end)=[];
intSplitTrialsPerOri = min(floor(cellfun(@sum,cellSelect)/intQuantiles));
indZeroVarNeurons = false(intRespN,1);
for intQ=1:intQuantiles
	vecUseTrialsTemp = vecStartTrials(intQ):(vecStartTrials(intQ)+intSplitTrialsPerOri-1);
	for intStim=1:intStimNum
		vecThisStim = find(cellSelect{intStim});
		[vecSorted,vecReorder]=sort(vecPopRate_pre(vecThisStim));
		vecQualifyingTrials = vecThisStim(vecReorder(vecUseTrialsTemp));
		indZeroVarNeurons = indZeroVarNeurons | (var(matMeanRate_pre(:,vecQualifyingTrials),[],2) == 0);
	end
end
indZeroVarNeurons = false(size(indZeroVarNeurons));
vecUseNeurons = find(~indZeroVarNeurons);
vecRemNeurons = find(indZeroVarNeurons);
cellSpikeTimesPerCellPerTrial(vecRemNeurons,:) = [];
intNeuronNum = numel(vecUseNeurons);
intTrialsPerQ = intSplitTrialsPerOri*intStimNum;

%simple "rate code"
matSpikeCounts = cellfun(@(x) sum(x>dblUseStartT),cellSpikeTimesPerCellPerTrial);
matMeanRate = matSpikeCounts./dblUseMaxDur;
