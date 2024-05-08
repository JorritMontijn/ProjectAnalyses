%define files
strDataPathSim = 'F:\Data\Processed\PopTimeCoding\';
if ~exist(strDataPathSim,'dir')
	strDataPathSim = 'D:\Data\Processed\PopTimeCoding\';
end
strLoadFile= 'Simulation_xAreaDistributed_SG18_2019-07-04.mat';
strSimName = 'SimDG18';

strNpxDir= 'F:\Drive\PopTimeCoding\data';
strNpxFiles = 'Q1cData_Rec*_Real_SGSVar*.mat';
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

%% reformat data
%match numbers in model and split data in random sets
sSim = load(fullpath(strDataPathSim,strLoadFile));

%define neurons
vecCellArea = sSim.vecCellArea;
vecNeuronsV1 = find(vecCellArea==1);
vecNeuronsV2 = find(vecCellArea==2);
vecUseNeurons = vecNeuronsV1(randperm(numel(vecNeuronsV1),sum(vecNeuronsPerRec)));
vecEndNeuronPerRec = cumsum(vecNeuronsPerRec);
intTotNum = numel(vecCellArea);

%transform synapse data to connectivity matrix
matSynFromTo=sSim.matSynFromTo;
vecSynVal = sSim.vecSynWeight .* sSim.vecSynConductance .* ((sSim.vecSynExcInh==1)*2-1);
intSynNum = numel(vecSynVal);
matConn = zeros(intTotNum,intTotNum);
vecLim = max(abs([min(vecSynVal) max(vecSynVal)]))*[-1 1];
colormap(redblue)
for i=1:intSynNum
	intFrom=matSynFromTo(i,1);
	intTo=matSynFromTo(i,2);
	matConn(intFrom,intTo) = vecSynVal(i);
end
%{
%filter
matConnPlus = matConn.*double(matConn>0);
matConnMinus = matConn.*double(matConn<0);

matConnMinusNorm=imnorm(-matConnMinus);
vecV1=1:1200;
vecV2=1201:2400;
matConnPlusNorm=zeros(size(matConnPlus));
matConnPlusNorm(vecV1,vecV1)=imnorm(matConnPlus(vecV1,vecV1));
matConnPlusNorm(vecV1,vecV2)=imnorm(matConnPlus(vecV1,vecV2));
matConnPlusNorm(vecV2,vecV2)=imnorm(matConnPlus(vecV2,vecV2));

matConnBig = ones([size(matConnPlus) 3]);
matConnBig(:,:,1) = 1-matConnMinusNorm;
matConnBig(:,:,2) = 1-sign(matConnPlusNorm+matConnMinusNorm);
matConnBig(:,:,3) = 1-matConnPlusNorm;

figure;
imshow(matConnBig);
%}

%%
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

%remove cells except spike times and area
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

%% save V2
strFile = [strSimName '_V2pop.mat'];
fprintf('   Saving V2 data - %s [%s]\n',strFile,getTime);
%select
vecOrigIds = vecNeuronsV2;

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
sRec.matConn = matConn;

%save sim rec file
save(fullpath(strDataPathSim,strFile),'-struct','sRec');

%% save V1 per rec
for intRec=[]%1:numel(vecNeuronsPerRec)
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