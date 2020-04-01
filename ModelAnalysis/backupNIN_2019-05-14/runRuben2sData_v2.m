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
intLoadExp = -13;
if intLoadExp == -12
	strExperiment = '12712hs';
	vecOrientations = [-7 0 14]+90;
	vecSpatialFrequencies = 1.5*ones(size(vecOrientations)); %cycle per degree
	vecContrasts = 0.25*ones(size(vecOrientations));
	strStimField = 'respt';
	strITIField = 'respt_blk';
elseif intLoadExp == -13
	strExperiment = '12842hs';
	vecOrientations = [-9 0 18]+45;
	vecSpatialFrequencies = 2*ones(size(vecOrientations)); %cycle per degree
	vecContrasts = 0.25*ones(size(vecOrientations));
	strStimField = 'respt_corr';
	strITIField = 'respt_corr_blk';
elseif intLoadExp == -14
	strExperiment = '129153hs';
	vecOrientations = [-7 0 14]+60;
	vecSpatialFrequencies = 2*ones(size(vecOrientations)); %cycle per degree
	vecContrasts = 0.1*ones(size(vecOrientations));
	strStimField = 'respt_corr';
	strITIField = 'respt_corr_blk';
end

%% set paths
strSourceDir = 'D:\Data\Processed\CoenCagliData\';
strTargetDir = 'A:\SimAggregates\';

%% define transformation settings
intOris = numel(vecOrientations);
vecTemporalFrequencies = 0*ones(size(vecOrientations)); %3-6.25Hz
vecLuminance = ones(size(vecOrientations));

%% load data file
fprintf('Starting transformation of <%s> [%s]\n',strExperiment,getTime);
sLoad = load([strSourceDir strExperiment '.mat']);
return
%[trials x neurons]
%sLoad.V1data_ORI0_NOISE0, sLoad.V1data_ORI1_NOISE0, sLoad.V1data_ORI2_NOISE0, sLoad.V1data_ORI0_NOISE1, sLoad.V1data_ORI1_NOISE1, sLoad.V1data_ORI2_NOISE1, sLoad.ORI0, sLoad.ORI1, sLoad.ORI2

%% meta data
%identifier
strConnFile = strcat('RubenExp',strExperiment);
dblTransformToDeltaT = 0.5/1000;

%% transform data
vecRemNeurons = isnan(sLoad.SNR) | sLoad.SNR < 2.5;
matStim = sLoad.(strStimField);
[intNeurons,intOris,intNoise,intReps,intStimMS] = size(matStim);
dblStimDur = sLoad.stimdur;

matITI = sLoad.(strITIField);
[intNeurons,intOris,intNoise,intReps,intITIMS] = size(matITI);
dblITIDur = sLoad.blkdur;
vecUseNeurons = find(~vecRemNeurons);
cellSpikeTimesV1 = cell(1,numel(vecUseNeurons));
dblTrialDur = dblStimDur + dblITIDur;


%% transform data
vecTrialStimRep = [];
vecTrialStimType = [];
vecTrialStartSecs = [];
dblCurT = 0;
vecRem = false(1,numel(vecUseNeurons));
vecUseOris = [1 2];
for intOri=vecUseOris
	vecTrialStartSecs = [vecTrialStartSecs (dblCurT+(dblTrialDur:dblTrialDur:(dblTrialDur*intReps))-dblTrialDur)];
	vecTrialStimType = [vecTrialStimType intOri*ones(1,intReps)];
	vecTrialStimRep = [vecTrialStimRep 1:intReps];
	for intNeuronIdx = 1:numel(vecUseNeurons)
		fprintf('Ori %d/%d, neuron %d/%d [%s]\n',intOri,2,intNeuronIdx,numel(vecUseNeurons),getTime);
		
		intNeuron = vecUseNeurons(intNeuronIdx);
		matTheseSpikes = squeeze(matStim(intNeuron,intOri,1,:,:));
		[row,col]=find(matTheseSpikes>0);
		vecSpikeTimes = dblTrialDur*(row-1) + col/1000;
		if numel(vecSpikeTimes) < 10
			vecRem(intNeuronIdx) = true;
		end
		matTheseITISpikes = squeeze(matITI(intNeuron,intOri,1,:,:));
		[row,col]=find(matTheseITISpikes>0);
		vecSpikeTimesITI = dblTrialDur*(row-1) + col/1000 + dblStimDur;
		
		cellSpikeTimesV1{intNeuronIdx} = [cellSpikeTimesV1{intNeuronIdx} dblCurT+sort(cat(1,vecSpikeTimes(:),vecSpikeTimesITI(:)),'ascend')'];
	end
	dblCurT = dblCurT + vecTrialStartSecs(end)+dblTrialDur;
end
cellSpikeTimesV1(vecRem) = [];
intCellsV1 = numel(cellSpikeTimesV1);

%% V2 placeholders
intCellsV2 = 0;

%build V1/V2 combos
cellSpikeTimesV2 = {};
intCortexCells = intCellsV1 + intCellsV2;
vecCellArea = ones(1,intCellsV1);

%% build sData
sData = struct;

%identifier
sData.strConnFile = strConnFile;
sData.vecOverallT = dblTransformToDeltaT:dblTransformToDeltaT:(vecTrialStartSecs(end)+dblTrialDur);
sData.dblDeltaT = dblTransformToDeltaT;
sData.strStimType = 'OriDrift2';
sData.cellSpikeTimesCortex = cat(2,cellSpikeTimesV1,cellSpikeTimesV2);

sData.vecTrialStartSecs = vecTrialStartSecs(:)';
sData.vecStimStartSecs = vecTrialStartSecs(:)';
sData.vecStimStopSecs = vecTrialStartSecs(:)'+dblStimDur;
sData.vecTrialEndSecs = vecTrialStartSecs(:)'+dblStimDur+dblITIDur;
sData.vecTrialStimRep = vecTrialStimRep(:)';
sData.vecTrialStimType = vecTrialStimType(:)';

%sConnParams
sData.intCellsV1 = intCellsV1;
sData.intCellsV2 = 0;
sData.intCortexCells = intCortexCells;
sData.vecCellArea = vecCellArea;

%sStimInputs
sData.dblTrialDur=dblTrialDur;
sData.vecTrialOris=vecOrientations(vecTrialStimType);
sData.vecTrialOriIdx=vecTrialStimType;
sData.vecTrialStimRep = vecTrialStimRep;
sData.vecStimTypeOris=vecOrientations(vecUseOris);
sData.vecStimTypeSFs=vecSpatialFrequencies(vecUseOris);
sData.vecStimTypeTFs=vecTemporalFrequencies(vecUseOris);
sData.vecStimTypeContrasts = vecContrasts(vecUseOris);
sData.vecStimTypeLuminance = vecLuminance(vecUseOris);

%% build stimulus properties
cellParamTypes{1} = sData.vecStimTypeOris;
cellParamTypes{2} = sData.vecStimTypeSFs;
cellParamTypes{3} = sData.vecStimTypeTFs;
cellParamTypes{4} = sData.vecStimTypeContrasts ;
cellParamTypes{5} = sData.vecStimTypeLuminance;
matStimTypeCombos = buildStimCombos(cellParamTypes);

sData.matStimTypeCombos=matStimTypeCombos;

%% save
strOutputFile = strcat('Simulation_xAreaExperiment_Ruben',strExperiment,'_',getDate,'.mat');
fprintf('%s transformation complete, saving data to <%s%s> [%s]\n',strExperiment,strTargetDir,strOutputFile,getTime);
save([strTargetDir strOutputFile],'-struct','sData','-v7.3');
fprintf('Done! [%s]\n',getTime);



