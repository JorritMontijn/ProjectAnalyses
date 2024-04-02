%define files
strDataPathSim = 'F:\Data\Processed\PopTimeCoding\';
if ~exist(strDataPathSim,'dir')
	strDataPathSim = 'D:\Data\Processed\PopTimeCoding\';
end
strLoadFile= 'Simulation_xAreaDistributed_SG18_2019-07-04.mat';
strSimName = 'SimDG18';

strNpxDir= 'C:\Drive\PopTimeCoding\data';
strNpxFiles = 'Q1cData_Rec*_Real_SGSFixed10.mat';
sNpx = dir(fullpath(strNpxDir,strNpxFiles));

%find neurons per rec
vecNeuronsPerRec=nan(size(sNpx));
cellRec = cell(size(sNpx));
for intRec=1:numel(vecNeuronsPerRec)
	fprintf('   Loading %d/%d - %s [%s]\n',intRec,numel(vecNeuronsPerRec),sNpx(intRec).name,getTime);
	sLoad = load(fullpath(sNpx(intRec).folder,sNpx(intRec).name));
	vecNeuronsPerRec(intRec) = numel(sLoad.vecTuningPerCell);
	cellRec{intRec} = getFlankedBy(sLoad.strRec,'','R01');
end
vecUseNeurons = randperm(1200,sum(vecNeuronsPerRec));
vecEndNeuronPerRec = cumsum(vecNeuronsPerRec);

%% reformat data
%match numbers in model and split data in random sets
sSim = load(fullpath(strDataPathSim,strLoadFile));

%remove synapse variables;
cellFields = fieldnames(sSim);
indSyn = contains(cellFields,'vecSyn');
sSim=rmfield(sSim,cellFields(indSyn));

%save some trial vars, delete rest
vecTrialOris = sSim.vecTrialOris;
vecTrialPhaseNoise = sSim.vecTrialPhaseNoise;
vecTrialEndSecs = sSim.vecTrialEndSecs;
vecTrialStartSecs = sSim.vecTrialStartSecs;
vecTrialStimType = sSim.vecTrialStimType;
indTrial = contains(cellFields,'vecTrial');
sSim=rmfield(sSim,cellFields(indTrial));

%remove exemplar matrices
cellFields = fieldnames(sSim);
indMat = contains(cellFields,'mat');
sSim=rmfield(sSim,cellFields(indMat));

%remove cells except spike times
cellSpikeTimesCortex = sSim.cellSpikeTimesCortex;
cellFields = fieldnames(sSim);
indCell = contains(cellFields,'cell');
sSim=rmfield(sSim,cellFields(indCell));

%move stim types
indStim = contains(cellFields,'vecStimType');
cellStimFields = cellFields(indStim);
sStim = struct;
for intStimVar=1:numel(cellStimFields)
	strField = cellStimFields{intStimVar};
	strAssignField = strField((numel('vecStimType')+1):end);
	vecData = sSim.(strField);
	for intType = 1:numel(vecData)
		sStim(intType).(strAssignField) = vecData(intType);
	end
end
sSim=rmfield(sSim,cellFields(indStim));

%move cell data
intNeuronNum = numel(cellSpikeTimesCortex);
cellFields = fieldnames(sSim);
vecEntriesPerField = structfun(@numel,sSim);
indNeuron = vecEntriesPerField == intNeuronNum;
cellNeuronFields = cellFields(indNeuron);
sNeuron = struct;
for intNeuronVar=1:numel(cellNeuronFields)
	strField = cellNeuronFields{intNeuronVar};
	if strcmp(strField(1:numel('vecCell')),'vecCell')
		strAssignField = strField((numel('vecCell')+1):end);
	elseif strcmp(strField(1:numel('vecPref')),'vecPref')
		strAssignField = strField((numel('vec')+1):end);
	elseif strcmp(strField,'cellSpikeTimesCortex')
		%do nothing
		continue;
	else
		error
	end
	vecData = sSim.(strField);
	for intType = 1:numel(vecData)
		sNeuron(intType).(strAssignField) = vecData(intType);
	end
end
sSim=rmfield(sSim,cellNeuronFields);

%remove scalars and short vectors
indLowFields = vecEntriesPerField < 10;
cellLowFields = cellFields(indLowFields);
sSim=rmfield(sSim,cellLowFields);

%remove misc
vecOverallT = sSim.vecOverallT;
sSim=rmfield(sSim,{'varDeltaSyn','vecCellTypesV2','vecDefinitionV1PrefOri','vecOrientationNoise','vecOverallT'});

%add data
sSim.sStim = sStim;
sSim.sNeuron = sNeuron;
sSim.vecTrialOris = vecTrialOris;
sSim.vecTrialPhaseNoise = vecTrialPhaseNoise;
sSim.vecTrialEndSecs = vecTrialEndSecs;
sSim.vecTrialStartSecs = vecTrialStartSecs;
sSim.vecTrialStimType = vecTrialStimType;


%% save per rec
for intRec=1:numel(vecNeuronsPerRec)
	strFile = [strSimName '_MatchedTo' cellRec{intRec} '.mat'];
	fprintf('   Saving %d/%d - %s [%s]\n',intRec,numel(vecNeuronsPerRec),strFile,getTime);
	%select
	intRecNeurons = vecNeuronsPerRec(intRec);
	vecRecNeurons = (vecEndNeuronPerRec(intRec)-intRecNeurons+1):vecEndNeuronPerRec(intRec);
	vecOrigIds = vecUseNeurons(vecRecNeurons);
	
	%transform spike times
	cellSpikeTimesSamples = cellSpikeTimesCortex(vecOrigIds);
	cellSpikeTimesSecs = cell(size(cellSpikeTimesSamples));
	for i=1:numel(cellSpikeTimesSecs)
		cellSpikeTimesSecs{i} = vecOverallT(cellSpikeTimesSamples{i}-1);
	end
	
	%build rec struct
	sRec = sSim;
	sRec.vecOrigIds = vecOrigIds;
	sRec.cellSpikeTimes = cellSpikeTimesSecs;
	sRec.sNeuron = sRec.sNeuron(vecOrigIds);
	
	%save sim rec file
	save(fullpath(strDataPathSim,strFile),'-struct','sRec');
end
fprintf('Done!\n');