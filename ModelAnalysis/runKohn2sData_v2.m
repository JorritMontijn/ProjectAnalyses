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
intLoadExp = -8;
if intLoadExp == -8
	%monyetp31: 124 neurons (recorded on 144 channels) in awake V1 in a monkey doing passive fixation (1 deg window). 2 gratings with 5 deg offset in orientation. Animal performed 831 trials with 3 stim presentations per trial, 200 ms each  with 150 ms ISI. Effectively ~2500 stimulus presentations divided across the two stimuli.
	vecOris = [0 5];
	strExperiment = 'cadetp70_alex';
elseif intLoadExp == -9
	%monyetp33: 106 neurons in awake V1 with 4 gratings (0,5 and 90,95). 845 trials or~600 presentations of each stimulus.
	vecOris = [0 5 90 95];
	strExperiment = 'cadetp149_alex';
elseif intLoadExp == -10
	%cadetp13: 126 neurons in awake V1, 1348 trials or~1000 presentations of each stimulus, 4 gratings (**as in monyetp33 but 0 deg grating is 0.25 contrast whereas rest are contrast 1**).
	vecOris = [0 5 90 95];
	vecContrasts = [0.25 1 1 1];
	strExperiment = 'cadetp13_alex';
end



%% set paths
strSourceDir = 'D:\Data\Processed\AdamKohnData\';
strTargetDir = 'A:\SimAggregates\';

%% load data file
fprintf('Starting transformation of <%s> [%s]\n',strExperiment,getTime);
sLoad = load([strSourceDir strExperiment '.mat']);

%% define transformation settings
vecOrientations =vecOris;
vecSpatialFrequencies = 1*ones(size(vecOrientations)); %cycle per degree
vecTemporalFrequencies = 1*ones(size(vecOrientations)); %3-6.25Hz
if ~exist('vecContrasts','var'),vecContrasts = ones(size(vecOrientations));end
vecLuminance = ones(size(vecOrientations));

dblTransformFromDeltaT = 1e1;
dblTransformToDeltaT = 0.5e-3;
dblTotalTrialDur = mean(sum(sLoad.dur,2));
intTrials = size(sLoad.dur,1);
intStimsPerTrial = size(sLoad.dur,2);

%% get stimulus timing and types
%stimulus properties:
%vecSFs = repmat([1 2 1/2],[1 6]);
%vecOris = reshape(repmat([-5 0 5 85 90 95],[3 1]),[1 18]);
vecTrialStart = (0:(intTrials-1))'*dblTotalTrialDur;
matStimStop = cumsum(sLoad.dur,2)+repmat(vecTrialStart,[1 intStimsPerTrial]);
matStimStart = matStimStop-sLoad.dur;
vecStimPeriods = sLoad.stim(1,:)~=1;
matStimStop = matStimStop(:,vecStimPeriods);
matStimStart = matStimStart(:,vecStimPeriods);
matStimIdx = label2idx(sLoad.stim(:,vecStimPeriods));

%% transform V1A1
fprintf('Transforming V1A1 data... [%s]\n',getTime);

intNeuronsV1A1 = size(sLoad.resp_v1_times,1);
matOrigin0 = arrayfun(@(cell,mat) cell{1} + mat,sLoad.resp_v1_times',repmat(vecTrialStart,[1 intNeuronsV1A1]),'UniformOutput',false);
cellSpikeTimesV1A1 = cell(1,intNeuronsV1A1);
for intNeuron=1:intNeuronsV1A1
	cellSpikeTimesV1A1{intNeuron} = cell2mat(matOrigin0(:,intNeuron));
end
	
%% transform V1A2
fprintf('Transforming V1A2 data... [%s]\n',getTime);

intNeuronsV1A2 = size(sLoad.resp_v1b_times,1);
matOrigin0A2 = arrayfun(@(cell,mat) cell{1} + mat,sLoad.resp_v1b_times',repmat(vecTrialStart,[1 intNeuronsV1A2]),'UniformOutput',false);
cellSpikeTimesV1A2 = cell(1,intNeuronsV1A2);
for intNeuron=1:intNeuronsV1A2
	cellSpikeTimesV1A2{intNeuron} = cell2mat(matOrigin0A2(:,intNeuron));
end
	
%% build stimulus vectors
%identifier
strConnFile = strcat('KohnExp2_',strExperiment);
vecStimStart = matStimStart(:)';
[vecStimStartSecs,vecSortTrials] = sort(vecStimStart,'ascend');
vecTrialStartSecs = vecStimStartSecs - vecStimStartSecs(1);
vecStimStop = matStimStop(:)';
vecStimStopSecs = vecStimStop(vecSortTrials);
vecTrialEndSecs = vecStimStopSecs;
vecTrialEndSecs(end) = vecTrialEndSecs(end)+0.1;
vecStimIdx = matStimIdx(:)';
vecTrialStimType = vecStimIdx(vecSortTrials);
vecTrialOriIdx=vecTrialStimType;
vecTrialOris = vecOrientations(vecTrialOriIdx);
vecTrialStimRep = zeros(size(vecTrialOris));
intStimTypes = numel(unique(vecTrialStimType));
for intTrialType = 1:intStimTypes
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
sData.strConnFile = strConnFile;
sData.vecOverallT = dblTransformToDeltaT:dblTransformToDeltaT:vecTrialEndSecs(end);
sData.dblDeltaT = dblTransformToDeltaT;
sData.strStimType = sprintf('OriDrift%d',intStimTypes);
sData.cellSpikeTimesCortex = cat(2,cellSpikeTimesV1A1,cellSpikeTimesV1A2);

sData.vecTrialStartSecs = vecTrialStartSecs(:)';
sData.vecStimStartSecs = vecStimStartSecs(:)';
sData.vecStimStopSecs = vecStimStopSecs(:)';
sData.vecTrialEndSecs = vecTrialEndSecs(:)';
sData.vecTrialStimRep = vecTrialStimRep(:)';
sData.vecTrialStimType = vecTrialStimType(:)';

%sConnParams
sData.intCellsV1 = intNeuronsV1A1+intNeuronsV1A2;
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

