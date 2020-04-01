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
clearvars;
intLoadExp = -2;
if intLoadExp == -1
	strExperiment = '106r001p26-t';
elseif intLoadExp == -2
	strExperiment = '107l003p143-t'
end

%% set paths
strSourceDir = 'D:\Data\Processed\AdamKohnData\';
strTargetDir = 'A:\SimAggregates\';

%% define transformation settings
vecOrientations = 0:22.5:157.5;
vecSpatialFrequencies = ones(size(vecOrientations)); %cycle per degree
vecTemporalFrequencies = 4*ones(size(vecOrientations)); %3-6.25Hz
vecContrasts = ones(size(vecOrientations));
vecLuminance = ones(size(vecOrientations));

vecStampOrientations = [nan vecOrientations];
dblTransformFromDeltaT = 1e-3;
dblTransformToDeltaT = 0.5e-3;
dblStimDurMS=1280;
dblITIDurMS=1500;
dblTotalTrialDurMS = dblStimDurMS + dblITIDurMS;
vecStampDurMS = isnan(vecStampOrientations)*dblITIDurMS + ~isnan(vecStampOrientations)*dblStimDurMS;

%% load data file
fprintf('Starting transformation of <%s> [%s]\n',strExperiment,getTime);
sLoad = load([strSourceDir strExperiment '.mat']);

%% get stimulus timing and types
vecTrialStimType = double(sLoad.stim); %vector [0:8] -> [blank 0:22.5:157.5]
vecTrialStimDur = vecStampDurMS(vecTrialStimType+1);
vecAddCumTimeStim = int32([0 cumsum(vecTrialStimDur(1:(end-1)))])';
vecStimStart = vecAddCumTimeStim(vecTrialStimType ~= 0);
vecStimStop = vecAddCumTimeStim(vecTrialStimType == 0);
vecStimType = vecTrialStimType(vecTrialStimType ~= 0);

%% transform V1
fprintf('Transforming V1 data... [%s]\n',getTime);
matTimeStampsV1 = sLoad.output.nev.TimeStamps;
%get events to be removed
indStimStampsV1 = matTimeStampsV1(:,1)==2000;
indJunkClustersV1 = matTimeStampsV1(:,2)==0 | matTimeStampsV1(:,2)==255;
%remove events
matTimeStampsV1 = matTimeStampsV1(~(indStimStampsV1 | indJunkClustersV1),:);

%get trial indices for all events
vecTrialIdxV1 = matTimeStampsV1(:,4);

%transform time stamps to cumulative time
vecTimeStampsV1 = int32(matTimeStampsV1(:,3));
vecTimeStampsV1_AddStimCumSum = vecAddCumTimeStim(vecTrialIdxV1);
vecTimeStampsV1CumSum = vecTimeStampsV1 + vecTimeStampsV1_AddStimCumSum;

%% transform event vector to cell-based spiking arrays for V1
intElectrodesV1 = max(matTimeStampsV1(:,1))+1;
vecUniqueUnitEventsV1 = matTimeStampsV1(:,1)+(matTimeStampsV1(:,2)*intElectrodesV1);
vecUniqueUnitsV1 = unique(vecUniqueUnitEventsV1);
intCellsV1 = numel(unique(vecUniqueUnitEventsV1));
cellSpikeTimesV1 = cell(1,intCellsV1);
for intCellV1=1:intCellsV1
	intUnit = vecUniqueUnitsV1(intCellV1);
	indEvents = vecUniqueUnitEventsV1==intUnit;
	vecSpikeTimes = vecTimeStampsV1CumSum(indEvents)';
	cellSpikeTimesV1{intCellV1} = ((double(vecSpikeTimes)-1)*dblTransformFromDeltaT);
end
fprintf('%d V1 cells detected and transformed [%s]\n',intCellsV1,getTime);

%% transform V2
fprintf('Transforming V2 data... [%s]\n',getTime);
matTimeStampsV2 = sLoad.output.plx.TimeStamps;
%get events to be removed
indStimStampsV2 = matTimeStampsV2(:,1)==2000;
indJunkClustersV2 = matTimeStampsV2(:,2)==0 | matTimeStampsV2(:,2)==255;
%remove events
matTimeStampsV2 = matTimeStampsV2(~(indStimStampsV2 | indJunkClustersV2),:);

%get trial indices for all events
vecTrialIdxV2 = matTimeStampsV2(:,4);

%transform time stamps to cumulative time
vecTimeStampsV2 = int32(matTimeStampsV2(:,3));
vecTimeStampsV2_AddStimCumSum = vecAddCumTimeStim(vecTrialIdxV2);
vecTimeStampsV2CumSum = vecTimeStampsV2 + vecTimeStampsV2_AddStimCumSum;


%% transform event vector to cell-based spiking arrays for V2
intElectrodesV2 = max(matTimeStampsV2(:,1))+1;
vecUniqueUnitEventsV2 = matTimeStampsV2(:,1)+(matTimeStampsV2(:,2)*intElectrodesV2);
vecUniqueUnitsV2 = unique(vecUniqueUnitEventsV2);
intCellsV2 = numel(unique(vecUniqueUnitEventsV2));
cellSpikeTimesV2 = cell(1,intCellsV2);
for intCellV2=1:intCellsV2
	intUnit = vecUniqueUnitsV2(intCellV2);
	indEvents = vecUniqueUnitEventsV2==intUnit;
	vecSpikeTimes = vecTimeStampsV2CumSum(indEvents)';
	cellSpikeTimesV2{intCellV2} = ((double(vecSpikeTimes)-1)*dblTransformFromDeltaT);
end
fprintf('%d V2 cells detected and transformed [%s]\n',intCellsV2,getTime);

%build V1/V2 combos
intCortexCells = intCellsV1 + intCellsV2;
vecCellArea = cat(2,ones(1,intCellsV1),2*ones(1,intCellsV2));

%% build stimulus vectors
%identifier
strConnFile = strcat('KohnExp',strExperiment);
vecStimStartSecs = double(vecStimStart)*dblTransformFromDeltaT;
vecTrialStartSecs = vecStimStartSecs;
vecStimStopSecs = double(vecStimStop)*dblTransformFromDeltaT;
vecTrialEndSecs = vecStimStopSecs + dblStimDurMS*dblTransformFromDeltaT;
vecTrialStimType = vecStimType;
vecTrialOriIdx=vecTrialStimType;
vecTrialOris = vecOrientations(vecTrialOriIdx);
vecTrialStimRep = zeros(size(vecTrialOris));
for intTrialType = unique(vecTrialStimType)
	indTheseTrials = vecTrialStimType==intTrialType;
	vecTrialStimRep(indTheseTrials) = 1:sum(indTheseTrials);
end

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
sData.strConnFile = strcat('KohnExp',strExperiment);
sData.vecOverallT = dblTransformToDeltaT:dblTransformToDeltaT:vecTrialEndSecs(end);
sData.dblDeltaT = dblTransformToDeltaT;
sData.strStimType = 'OriDrift8';
sData.cellSpikeTimesCortex = cat(2,cellSpikeTimesV1,cellSpikeTimesV2);

sData.vecTrialStartSecs = vecTrialStartSecs(:)';
sData.vecStimStartSecs = vecStimStartSecs(:)';
sData.vecStimStopSecs = vecStimStopSecs(:)';
sData.vecTrialEndSecs = vecTrialEndSecs(:)';
sData.vecTrialStimRep = vecTrialStimRep(:)';
sData.vecTrialStimType = vecTrialStimType(:)';

%sConnParams
sData.intCellsV1 = intCellsV1;
sData.intCellsV2 = intCellsV2;
sData.intCortexCells = intCortexCells;
sData.vecCellArea = vecCellArea;

%sStimInputs
sData.dblTrialDur=dblTotalTrialDurMS*dblTransformFromDeltaT;
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

