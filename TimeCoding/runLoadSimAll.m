%% set variables
strArea = 'SimV1';
strAreaAbbr = 'V1';
if strcmp(strRunStim,'DG')
	strLoadStim = 'driftinggrating';
	strLoadFile= 'Simulation_xAreaDistributed_SG18_2019-07-04.mat';
else
	strLoadFile='Simulation_xAreaDistributed_IndRetNoise0_0_2019-08-28.mat';
	strLoadStim = '';
end

%% load data
sLoad = load(fullpath(strDataPathSim,strLoadFile));

%%
cellUseAreas{1} = strArea;
% concatenate stimulus structures
vecStimOnTime = sLoad.vecStimStartSecs;
vecStimOffTime = sLoad.vecStimStopSecs;

vecOrientation = sLoad.vecTrialOris;
vecTempFreq = sLoad.vecTrialTFs;
vecPhase = mod(sLoad.vecTrialPhaseNoise,2*pi)./(2*pi);
vecDelayTimeBy = vecPhase./vecTempFreq;
[vecOriIdx,vecUniqueOris,vecRepNum,cellSelect,vecTrialRepetition] = val2idx(vecOrientation);
indRem=vecTrialRepetition>min(vecRepNum);
vecOrientation(indRem) = [];
vecDelayTimeBy(indRem) = [];
vecStimOnTime(indRem) = [];
vecStimOffTime(indRem) = [];


%remove neurons in incorrect areas
indConsiderNeurons = sLoad.vecCellArea == 1;

%remove neurons outside subselection
vecConsider = find(indConsiderNeurons);
if ~exist('vecSubSelectIdx','var')
vecSubSelect = randperm(numel(vecConsider),intSelectCells);
vecSubSelectIdx = sort(vecConsider(vecSubSelect));
end
indSelectCells = false(1,numel(indConsiderNeurons));
indSelectCells(vecSubSelectIdx) = true;
strRec = 'SimSG18';

%subselect from total
indUseNeurons = indConsiderNeurons(:) & indSelectCells(:);
cellSpikeTimesSamples = sLoad.cellSpikeTimesCortex(indUseNeurons);
cellSpikeTimesRaw = cell(size(cellSpikeTimesSamples));
for i=1:numel(cellSpikeTimesRaw)
	cellSpikeTimesRaw{i} = sLoad.vecOverallT(cellSpikeTimesSamples{i}-1);
end

%% prep grating data
[vecStimIdx,vecUnique,vecRepNum,cellSelect,vecTrialRepetition] = val2idx(vecOrientation);
intTrialNum = numel(vecStimOnTime);
intOriNum = numel(unique(vecOrientation));
intRepNum = intTrialNum/intOriNum;
intNeuronsInArea = numel(cellSpikeTimesRaw);
intNeuronNum = intNeuronsInArea;
if intNeuronsInArea==0,return;end
