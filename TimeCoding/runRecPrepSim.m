
%% get matching recording data
strTarget = fullpath(sSimRecs(intRec).folder,sSimRecs(intRec).name);
close all;
strRec = sSimRecs(intRec).name(1:(end-4));
strRecOrig = strRec;

if strcmp(strRunStim,'DG')
	%load
	sLoad = load(strTarget);
	
	
	%% transform format
	cellUseAreas{1} = strArea;
	% concatenate stimulus structures
	vecStimOnTime = sLoad.vecStimStartSecs;
	vecStimOffTime = sLoad.vecStimStopSecs;
	dblStimDur = roundi(median(vecStimOffTime - vecStimOnTime),1,'ceil');
	
	vecOrientation = sLoad.vecTrialOris;
	vecTypeTFs = [sLoad.sStim.TFs];
	vecTempFreq = vecTypeTFs(sLoad.vecTrialStimType);
	vecPhase = mod(sLoad.vecTrialPhaseNoise,2*pi)./(2*pi);
	vecDelayTimeBy = vecPhase./vecTempFreq;
	[vecOriIdx,vecUniqueOris,vecRepNum,cellSelect,vecTrialRepetition] = val2idx(vecOrientation);
	indRem=vecTrialRepetition>min(vecRepNum);
	vecOrientation(indRem) = [];
	vecDelayTimeBy(indRem) = [];
	vecStimOnTime(indRem) = [];
	vecStimOffTime(indRem) = [];
	
	
	%remove neurons in incorrect areas
	indConsiderNeurons = [sLoad.sNeuron.Area] == 1;
	
	%subselect from total
	indUseNeurons = indConsiderNeurons(:);
	cellSpikeTimesRaw = sLoad.cellSpikeTimes(indUseNeurons);
	
	%% prep grating data
	[vecStimIdx,vecUnique,vecRepNum,cellSelect,vecTrialRepetition] = val2idx(vecOrientation);
	intTrialNum = numel(vecStimOnTime);
	intOriNum = numel(unique(vecOrientation));
	intRepNum = intTrialNum/intOriNum;
	intNeuronsInArea = numel(cellSpikeTimesRaw);
	intNeuronNum = intNeuronsInArea;
	if intNeuronsInArea==0,return;end
	
elseif strcmp(strRunStim,'NM')
	%to do
	
end
%{
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
%}