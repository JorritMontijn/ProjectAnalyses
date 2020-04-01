%{
Code	Description
1	Start trial
2	End trial
3	Fix acquired
4	Stim On
5	Target acquired
7	Juice
8	No Juice
11	Scan Code On
12	scan Code Off
13	Target1
14	Target2
40	Left
41	Right
79	Increment TrialRecord.CorrectTrialsWithinBlock
%}

%% set recording
clearvars;
intLoadExp = -21;
if intLoadExp == -21
	strFile = 'SaccadeTask-voltaire-10-06-2016-SpikeDataNoWaves';
	strExperiment = 'voltaire20160610';
end

%% set paths
strSourceDir = 'D:\Data\Processed\Averbeck\';
strTargetDir = 'A:\SimAggregates\';

%% load data file
fprintf('Starting transformation of <%s> [%s]\n',strFile,getTime);
sLoad = load([strSourceDir strFile '.mat']);
fprintf('Data loaded. Transforming... [%s]\n',getTime);

%% define transformation settings
dblTransformToDeltaT = 0.5e-3;
dblTransformFromDeltaT = 1e-3;
dblMinimumCellQuality = 1;

%% transform data
intTrials = size(sLoad.SPKstruc,2);
intElectrodes = numel(sLoad.SortInfo);
vecUnitsPerElectrode = nan(1,intElectrodes);
cellUnitRatings = cell(1,intElectrodes);
for intElectrode=1:intElectrodes
	intUnits = sLoad.SortInfo(intElectrode).numUnits;
	if isempty(intUnits),intUnits=0;end
	vecUnitsPerElectrode(intElectrode) = intUnits;
	cellUnitRatings{intElectrode} = sLoad.SortInfo(intElectrode).unitRating;
end
sRawData = sLoad.SPKstruc(vecUnitsPerElectrode>0,:);
vecErrorTrials = sLoad.errortrials(:,1);
cellUnitRatings(vecUnitsPerElectrode==0) = [];
vecUnitsPerElectrode(vecUnitsPerElectrode==0) = [];
vecCellNum = cellfun(@sum,cellfun(@ge,cellUnitRatings,cellfill(dblMinimumCellQuality,size(cellUnitRatings)),'UniformOutput',false));
clear sLoad;

%% build trial vectors & cell array for spiking
vecTrialStartSecs = nan(1,intTrials);
vecStimStartSecs  = nan(1,intTrials);
vecStimStopSecs = nan(1,intTrials);
vecTrialEndSecs = nan(1,intTrials);
vecTrialStimType = nan(1,intTrials);
	
%loop
intUnits = sum(vecCellNum);
cellSpikeTimesCortex = cell(1,intUnits);
intUnit = 0;
for intElectrode=1:size(sRawData,1)
	fprintf('Processing electrode %d/%d [%s]\n',intElectrode,size(sRawData,1),getTime);
	%get unit data
	intElecUnits = vecUnitsPerElectrode(intElectrode);
	vecQuality = cellUnitRatings{intElectrode};
	vecUseCells = flat(find(vecQuality >= dblMinimumCellQuality))';
	for intTrial=1:intTrials
		
		sElectrode = sRawData(intElectrode,intTrial);
		intOrigTrial = sElectrode.Trial;
		
		if intElectrode==1
			if ~ismember(intOrigTrial,vecErrorTrials)
				vecTrialStartSecs(intTrial) = sElectrode.BHV(find(sElectrode.BHV(:,1)==1,1,'last'),2);
				vecStimStartSecs(intTrial)  = sElectrode.BHV(find(sElectrode.BHV(:,1)==4,1,'last'),2);
				vecStimStopSecs(intTrial) = sElectrode.BHV(find(sElectrode.BHV(:,1)==5,1,'last'),2);
				vecTrialEndSecs(intTrial) = sElectrode.BHV(find(sElectrode.BHV(:,1)==5,1,'last'),2);
				vecTrialStimType(intTrial) = ismember(41,sElectrode.BHV(:,1))+1;
			end
		end
		
		%get spiking
		for intUseUnit=1:numel(vecUseCells)
			vecSpikes = vecUseCells(intUseUnit);
			cellSpikeTimesCortex{intUseUnit+intUnit} = [cellSpikeTimesCortex{intUseUnit+intUnit} sElectrode.TimeStamp(sElectrode.SpikeID==vecUseCells(intUseUnit))];
		end
		
	end
	
	%increment units with units on electrode
	intUnit = intUnit + numel(vecUseCells);
end
%clear data
clearvars sRawData sElectrode;

%transform to seconds
cellSpikeTimesCortex = cellfun(@mtimes,cellSpikeTimesCortex,cellfill(dblTransformFromDeltaT,size(cellSpikeTimesCortex)),'UniformOutput',false);

%remove error trials
vecTrialStartSecs(vecErrorTrials) = [];
vecStimStartSecs(vecErrorTrials) = [];
vecStimStopSecs(vecErrorTrials) = [];
vecTrialEndSecs(vecErrorTrials) = [];
vecTrialStimType(vecErrorTrials) = [];
intStimTypes = numel(unique(roundi(vecTrialStimType,3)));

%build stim rep vector
vecTrialStimRep = nan(size(vecTrialStimType));
for intStimType = 1:intStimTypes
	vecTheseTrials = vecTrialStimType==intStimType;
	vecRep = cumsum(vecTheseTrials);
	vecTrialStimRep(vecTheseTrials) = vecRep(vecTheseTrials);
end

%% build meta data
%identifier
strConnFile = strcat('AverExp1_',strExperiment);
vecTrialOriIdx=vecTrialStimType;
vecTrialOris = vecTrialStimType;

%% build stimulus properties
cellParamTypes = {1 2};
matStimTypeCombos = buildStimCombos(cellParamTypes);

%% build sData
sData = struct;

%identifier
sData.strConnFile = strConnFile;
sData.vecOverallT = dblTransformToDeltaT:dblTransformToDeltaT:(vecTrialEndSecs(end)*dblTransformFromDeltaT);
sData.dblDeltaT = dblTransformToDeltaT;
sData.strStimType = sprintf('Choice%d',intStimTypes);
sData.cellSpikeTimesCortex = cellSpikeTimesCortex;

sData.vecTrialStartSecs = dblTransformFromDeltaT*vecTrialStartSecs(:)';
sData.vecStimStartSecs = dblTransformFromDeltaT*vecStimStartSecs(:)';
sData.vecStimStopSecs = dblTransformFromDeltaT*vecStimStopSecs(:)';
sData.vecTrialEndSecs = dblTransformFromDeltaT*vecTrialEndSecs(:)';
sData.vecTrialStimRep = vecTrialStimRep(:)';
sData.vecTrialStimType = vecTrialStimType(:)';

%sConnParams
sData.intCellsV1 = intUnits;
sData.intCellsV2 = 0;
sData.intCortexCells = sData.intCellsV1;
sData.vecCellArea = ones(1,sData.intCortexCells);

%sStimInputs
sData.dblTrialDur=dblTransformFromDeltaT*min(vecStimStopSecs-vecStimStartSecs);
sData.vecTrialOris=vecTrialOris;
sData.vecTrialOriIdx=vecTrialOriIdx;
sData.vecTrialStimRep = vecTrialStimRep;
sData.vecStimTypeOris=1:intStimTypes;
sData.matStimTypeCombos=matStimTypeCombos;

%% save
strOutputFile = strcat('Simulation_xAreaExperiment_Averbeck',strExperiment,'_',getDate,'.mat');
fprintf('Transformation complete, saving data to <%s%s> [%s]\n',strTargetDir,strOutputFile,getTime);
save([strTargetDir strOutputFile],'-struct','sData','-v7.3');
fprintf('Done! [%s]\n',getTime);

