
%% get matching recording data
close all;
strRec = sAggStim(intRec).Exp;
strRecOrig = strRec;
sThisRec = sAggStim(strcmpi(strRec,{sAggStim(:).Exp}));
sThisSource = sAggSources(strcmpi(strRec,{sAggSources(:).Exp}));
if strcmp(strRunStim,'DG')
	%prep grating data
	[sUseNeuron,vecStimOnTime,vecStimOffTime,vecOrientation] = NpxPrepGrating(sAggNeuron,sThisRec,cellUseAreas);
	[vecStimIdx,vecUniqueOris,vecRepNum,cellSelect,vecTrialRepetition] = val2idx(vecOrientation);
	intTrialNum = numel(vecStimOnTime);
	intOriNum = numel(unique(vecOrientation));
	intRepNum = intTrialNum/intOriNum;
	intNeuronsInArea = numel(sUseNeuron);
	intNeuronNum = intNeuronsInArea;
	if intNeuronsInArea==0,return;end

	%% get neurons in this area
	indArea1Neurons = contains({sUseNeuron.Area},strArea,'IgnoreCase',true);
	intNeuronsInArea = sum(indArea1Neurons);
	intNeuronNum = intNeuronsInArea;
	if intNeuronsInArea==0,return;end

	%% prep data
	sArea1Neurons = sUseNeuron(indArea1Neurons);
	cellSpikeTimesRaw = {sArea1Neurons.SpikeTimes};
	vecOri180 = mod(vecOrientation,180)*2;
	[matMeanRate,indTuned,cellSpikeTimes,sOut,cellSpikeTimesPerCellPerTrial,vecStimOnStitched,vecNonStat,dblBC,dblMaxDevFrac] = ...
		NpxPrepData(cellSpikeTimesRaw,vecStimOnTime,vecStimOffTime,vecOrientation);
	intTunedN = sum(indTuned);
	intRespN = size(matMeanRate,1);
	dblStimDur = roundi(median(vecStimOffTime - vecStimOnTime),1,'ceil');
	dblPreTime = -dblStartT;%0.3;
	dblPostTime = 0;%0.3;
	dblMaxDur = dblStimDur+dblPreTime+dblPostTime;
elseif strcmp(strRunStim,'NM')
	%to do

end

%% define quantiles and remove zero-variance neurons
%constants
[vecStimIdx,vecUnique,vecPriorDistribution,cellSelect,vecRepetition] = val2idx(vecOri180);
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
