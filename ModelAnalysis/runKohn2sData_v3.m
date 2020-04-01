%{
The variable 'output' contains both the V1 (output.nev) and V2 (output.plx) data.
The spike time (*.TimeStamps) are matrices with many rows and 4 columns:

col 1: the electrode number (1-96), electrode 2000 is the digital stimulus
	marker which indicates the onset of the stimulus;

col 2: sort code, where 0 and/or 255 are junk, but any other number is a
	valid SUA or MUA

col 3: the time at which the event occurred relative to stimulus onset

col 4: the trial number

The stimulus are presented for 1280 ms and separated by 1500 ms of blank screen (ISI).

The stimulus order is in the variable 'stim', where 1-8 are the 8 gratings
	(0:22.5:157.5); stimulus 0 is the blank (gray screen) between stimuli.
%}
%% set recording
%anesthetized V1 data with 18 stimuli, each shown 200 times, 1.28 s presentation. The stimuli are drifting gratings of 2 sets of 3 finely-spaced orientations, each shown at 3 finely-spaced spatial frequencies.
clearvars;
intLoadExp = -6;
if intLoadExp == -6
	strExperiment = '139l001p113_alex';
elseif intLoadExp == -7
	strExperiment = '139l001p115_alex';
end
vecOris = [-5 -5 -5 0 0 0 5 5 5 85 85 85 90 90 90 95 95 95];
vecSFs = repmat([1 2 0.5],[1 6]);
intStimTypes = numel(vecOris);
dblStimDur = 1.28;
intReps = 200;
%format: resp_times_stim %[stim x cell x trial]

%build indices
vecOriIdx = label2idx(vecOris);

%% set paths
strSourceDir = 'D:\Data\Processed\AdamKohnData\';
strTargetDir = 'A:\SimAggregates\';

%% load data file
fprintf('Starting transformation of <%s> [%s]\n',strExperiment,getTime);
sLoad = load([strSourceDir strExperiment '.mat']);

%% define transformation settings
vecOrientations =vecOris;
vecSpatialFrequencies = vecSFs; %cycle per degree
vecTemporalFrequencies = 1*ones(size(vecOrientations)); %3-6.25Hz
if ~exist('vecContrasts','var'),vecContrasts = ones(size(vecOrientations));end
vecLuminance = ones(size(vecOrientations));

dblTransformFromDeltaT = 1e1;
dblTransformToDeltaT = 0.5e-3;
dblTotalTrialDur = dblStimDur;
intTrials = intReps * numel(vecOris);

%% get stimulus timing and types
%stimulus properties:
%vecSFs = repmat([1 2 1/2],[1 6]);
%vecOris = reshape(repmat([-5 0 5 85 90 95],[3 1]),[1 18]);
vecTrialStartOneRep = (0:(intStimTypes-1))'*dblTotalTrialDur;
dblOneRepDur = intStimTypes * dblTotalTrialDur;
vecRepStart = reshape((0:(intReps-1))'*dblOneRepDur,[1 1 intReps]);
matTrialOffsets = bsxfun(@plus,vecTrialStartOneRep,vecRepStart);
matStimType = repmat((1:18)',[1 1 intReps]);
matStimRep = repmat(reshape((1:intReps),[1 1 intReps]),[intStimTypes 1 1]);

%% transform V1
fprintf('Transforming V1 data... [%s]\n',getTime);

intNeuronsV1 = size(sLoad.resp_times_stim,2);
cellOrigin0 = arrayfun(@(cell,mat) cell{1} + mat,sLoad.resp_times_stim,repmat(matTrialOffsets,[1 intNeuronsV1 1]),'UniformOutput',false);
cellSpikeTimesV1 = cell(1,intNeuronsV1);
for intNeuron=1:intNeuronsV1
	cellSpikeTimesV1{intNeuron} = cell2vec(cellOrigin0(:,intNeuron,:));
end
	
%% build stimulus vectors
%identifier
strConnFile = strcat('KohnExp3_',strExperiment);
vecStimStart = matTrialOffsets(:)';
[vecStimStartSecs,vecSortTrials] = sort(vecStimStart,'ascend');
vecTrialStartSecs = vecStimStartSecs - vecStimStartSecs(1);
vecStimStopSecs = vecStimStartSecs+dblStimDur;
vecTrialEndSecs = vecStimStopSecs;
vecTrialEndSecs(end) = vecTrialEndSecs(end)+0.1;
vecStimIdx = matStimType(:)';
vecTrialStimType = vecStimIdx(vecSortTrials);
vecTrialOriIdx=vecOriIdx(vecTrialStimType);
vecTrialOris = vecOris(vecTrialStimType);
vecStimRep = matStimRep(:)';
vecTrialStimRep = vecStimRep(vecSortTrials);
intStimTypes = numel(unique(vecTrialStimType));


%% build stimulus properties
cellParamTypes{1} = vecOrientations;
cellParamTypes{2} = vecSpatialFrequencies;
cellParamTypes{3} = vecTemporalFrequencies;
cellParamTypes{4} = vecContrasts;
cellParamTypes{5} = vecLuminance;
matStimTypeCombos = buildStimCombos(cellParamTypes);

%% build sData
sData = struct;

%identifier
sData.strConnFile = strConnFile;
sData.vecOverallT = dblTransformToDeltaT:dblTransformToDeltaT:vecTrialEndSecs(end);
sData.dblDeltaT = dblTransformToDeltaT;
sData.strStimType = sprintf('OriDrift%d',intStimTypes);
sData.cellSpikeTimesCortex = cellSpikeTimesV1;

sData.vecTrialStartSecs = vecTrialStartSecs(:)';
sData.vecStimStartSecs = vecStimStartSecs(:)';
sData.vecStimStopSecs = vecStimStopSecs(:)';
sData.vecTrialEndSecs = vecTrialEndSecs(:)';
sData.vecTrialStimRep = vecTrialStimRep(:)';
sData.vecTrialStimType = vecTrialStimType(:)';

%sConnParams
sData.intCellsV1 = intNeuronsV1;
sData.intCellsV2 = 0;
sData.intCortexCells = sData.intCellsV1;
sData.vecCellArea = ones(1,sData.intCortexCells);

%sStimInputs
sData.dblTrialDur=dblTotalTrialDur;
sData.vecTrialOris=vecTrialOris;
sData.vecTrialOriIdx=vecTrialOriIdx;
sData.vecTrialStimRep = vecTrialStimRep;
sData.vecStimTypeOris=vecOrientations;
sData.vecStimTypeSFs=vecSpatialFrequencies;
sData.vecStimTypeTFs=vecTemporalFrequencies;
sData.vecStimTypeContrasts = vecContrasts;
sData.vecStimTypeLuminance = vecLuminance;

sData.matStimTypeCombos=matStimTypeCombos;

%% save
strOutputFile = strcat('Simulation_xAreaExperiment_',strExperiment,'_',getDate,'.mat');
fprintf('Transformation complete, saving data to <%s%s> [%s]\n',strTargetDir,strOutputFile,getTime);
save([strTargetDir strOutputFile],'-struct','sData','-v7.3');
fprintf('Done! [%s]\n',getTime);

