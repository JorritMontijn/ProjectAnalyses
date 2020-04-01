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
intLoadExp = -11;
if intLoadExp == -11
	strExperiment = 'V1data';
end

%% set paths
strSourceDir = 'D:\Data\Processed\CoenCagliData\';
strTargetDir = 'A:\SimAggregates\';

%% define transformation settings
vecOrientations = [-7 0 14];
intOris = numel(vecOrientations);
vecSpatialFrequencies = ones(size(vecOrientations)); %cycle per degree
vecTemporalFrequencies = 0*ones(size(vecOrientations)); %3-6.25Hz
vecContrasts = 0.25*ones(size(vecOrientations));
vecLuminance = ones(size(vecOrientations));

vecStampOrientations = [vecOrientations];
dblTransformFromDeltaT = 250e-3;
dblTransformToDeltaT = 0.5e-3;
dblStimDur=0.250;
dblITIDur=0.250;
dblTotalTrialDurMS = dblStimDur + dblITIDur;
vecStampDurMS = isnan(vecStampOrientations)*dblITIDur + ~isnan(vecStampOrientations)*dblStimDur;

%% load data file
fprintf('Starting transformation of <%s> [%s]\n',strExperiment,getTime);
sLoad = load([strSourceDir strExperiment '.mat']);
%[trials x neurons]
%sLoad.V1data_ORI0_NOISE0, sLoad.V1data_ORI1_NOISE0, sLoad.V1data_ORI2_NOISE0, sLoad.V1data_ORI0_NOISE1, sLoad.V1data_ORI1_NOISE1, sLoad.V1data_ORI2_NOISE1, sLoad.ORI0, sLoad.ORI1, sLoad.ORI2

%% meta data
%identifier
strConnFile = strcat('RubenExp',strExperiment);

%% build stimulus properties
cellParamTypes{1} = vecOrientations;
cellParamTypes{2} = vecSpatialFrequencies;
cellParamTypes{3} = vecTemporalFrequencies;
cellParamTypes{4} = vecContrasts;
cellParamTypes{5} = vecLuminance;
matStimTypeCombos = buildStimCombos(cellParamTypes);

%% transform data
cellOriSubFields = {'ORI0','ORI1','ORI2'};
cellNoise = {'NOISE0','NOISE1'};
intCellsV1 =size(sLoad.V1data_ORI0_NOISE0,2);
for intNoise = [1 2]
	%% V1 data & stimulus vectors
	dblTimeOffset = 0;
	cellSpikeTimesV1 = cell(1,intCellsV1);
	vecTrialStartSecs = [];
	vecTrialStimRep = [];
	vecTrialStimType = [];
	for intOri=1:intOris
		strField = sprintf('%s_%s_%s',strExperiment,cellOriSubFields{intOri},cellNoise{intNoise});
		intRepetitions = size(sLoad.(strField),1);
		vecTrialStartSecs = [vecTrialStartSecs dblTimeOffset+[0:dblTotalTrialDurMS:(dblTotalTrialDurMS*intRepetitions-dblTotalTrialDurMS/2)]];
		vecTrialStimRep = [vecTrialStimRep 1:intRepetitions];
		vecTrialStimType = [vecTrialStimType intOri*ones(1,intRepetitions)];
		for intCellV1=1:intCellsV1
			vecAct = sLoad.(strField)(:,intCellV1);
			vecSpikeTimes = nan(1,sum(vecAct));
			intC=1;
			for intSpikeTrial=find(vecAct>0)'
				intSpikes = vecAct(intSpikeTrial);
				vecSpikeTimes(intC:(intC+intSpikes-1)) = (intSpikeTrial*dblITIDur + intSpikeTrial*dblStimDur - dblStimDur/2)+dblTransformToDeltaT*[1:intSpikes]+dblTimeOffset;
				intC = intC + intSpikes;
			end
			cellSpikeTimesV1{intCellV1} = [cellSpikeTimesV1{intCellV1} vecSpikeTimes];
		end
		dblTimeOffset = dblTimeOffset + dblTotalTrialDurMS*intRepetitions;
	end
	fprintf('%s; %d V1 cells detected and transformed;  [%s]\n',cellNoise{intNoise},intCellsV1,getTime);
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
	sData.vecOverallT = dblTransformToDeltaT:dblTransformToDeltaT:dblTimeOffset;
	sData.dblDeltaT = dblTransformToDeltaT;
	sData.strStimType = 'OriDrift3';
	sData.cellSpikeTimesCortex = cat(2,cellSpikeTimesV1,cellSpikeTimesV2);
	
	sData.vecTrialStartSecs = vecTrialStartSecs(:)';
	sData.vecStimStartSecs = vecTrialStartSecs(:)'+dblITIDur;
	sData.vecStimStopSecs = vecTrialStartSecs(:)'+dblITIDur+dblStimDur;
	sData.vecTrialEndSecs = vecTrialStartSecs(:)'+dblITIDur+dblStimDur;
	sData.vecTrialStimRep = vecTrialStimRep(:)';
	sData.vecTrialStimType = vecTrialStimType(:)';
	
	%sConnParams
	sData.intCellsV1 = intCellsV1;
	sData.intCellsV2 = 0;
	sData.intCortexCells = intCortexCells;
	sData.vecCellArea = vecCellArea;
	
	%sStimInputs
	sData.dblTrialDur=dblTotalTrialDurMS;
	sData.vecTrialOris=vecOrientations(vecTrialStimType);
	sData.vecTrialOriIdx=vecTrialStimType;
	sData.vecTrialStimRep = vecTrialStimRep;
	sData.vecStimTypeOris=vecOrientations;
	sData.vecStimTypeSFs=vecSpatialFrequencies;
	sData.vecStimTypeTFs=vecTemporalFrequencies;
	sData.vecStimTypeContrasts = vecContrasts;
	sData.vecStimTypeLuminance = vecLuminance;
	
	sData.matStimTypeCombos=matStimTypeCombos;
	
	%% save
	strOutputFile = strcat('Simulation_xAreaExperiment_Ruben',strExperiment,'_',cellNoise{intNoise},'_',getDate,'.mat');
	fprintf('%s transformation complete, saving data to <%s%s> [%s]\n',cellNoise{intNoise},strTargetDir,strOutputFile,getTime);
	save([strTargetDir strOutputFile],'-struct','sData','-v7.3');
	fprintf('Done! [%s]\n',getTime);
	
	
end



