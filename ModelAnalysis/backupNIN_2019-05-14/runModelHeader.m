%runModelHeader
%% load data
indRequiredPrePros = [true; ... %prepro
	false;... %prepro2
	true;...%preproLGN
	false;...%preproSync
	true]; %prepro0
indLoadedPrePros = false(size(indRequiredPrePros));
indRequiredPrePros(3) = boolLoadLGN;
indRequiredPrePros(4) = boolLoadSync;
if ~exist('boolLoadSpikeTimes','var')
	boolLoadSpikeTimes = false;
end
if (~exist('sData','var') && ~exist('sSimRun','var')) || ~exist('strLastSim','var') || ~strcmp(strLastSim,strSimulation)
	fprintf('Loading simulation <%s> [%s]\n',strSimulation,getTime);
	%clearvars -except intLoadSim strSimulation strFigDir bool*
	if ~exist('boolLoadLGN','var'),boolLoadLGN=false;end
	strLastSim = strSimulation;
	strPath = 'A:\SimAggregates\';
	if exist([strPath 'Simulation_' strSimulation '_prepro.mat'],'file')
		indLoadedPrePros(1) = true;
		sLoadPrePro = load([strPath 'Simulation_' strSimulation '_prepro.mat']); %matModelResp
	end
	if exist([strPath 'Simulation_' strSimulation '_prepro2.mat'],'file')
		indLoadedPrePros(2) = true;
		sLoadPrePro2 = load([strPath 'Simulation_' strSimulation '_prepro2.mat']);
	end
	if boolLoadLGN && exist([strPath 'Simulation_' strSimulation '_preproLGN.mat'],'file')
		indLoadedPrePros(3) = true;
		sLoadPreProLGN = load([strPath 'Simulation_' strSimulation '_preproLGN.mat']);%matModelRespLGN_ON
	end
	if boolLoadSync && exist([strPath 'Simulation_' strSimulation '_preproSyncCort.mat'],'file')
		indLoadedPrePros(4) = true;
		sLoadPreProSyncCort = load([strPath 'Simulation_' strSimulation '_preproSyncCort.mat']);%matModelRespLGN_ON
	end
	if exist([strPath 'Simulation_' strSimulation '_prepro0.mat'],'file')
		indLoadedPrePros(end) = true;
		sLoad = load([strPath 'Simulation_' strSimulation '_prepro0.mat']);%sData
	end
	clear sData;
	if ~all(indLoadedPrePros(indRequiredPrePros)) || boolLoadSpikeTimes
		sLoad = load([strPath 'Simulation_' strSimulation '.mat']);
	end
end

%% check if vars are saved as fields
if exist('sData','var') %do nothing
elseif isfield(sLoad,'sData')
	sData = sLoad.sData;
	clear sLoad;
elseif isfield(sLoad,'vecOverallT')
	sData = sLoad;
	clear sLoad;
end

%check if cellSpike arrays are uint32
if exist('sData','var') && isfield(sData,'cellSpikeTimesCortex') && ~isa(sData.cellSpikeTimesCortex{1},'double')
	dblStepSize = (sData.vecOverallT(end)-sData.vecOverallT(1))/(numel(sData.vecOverallT)-1);
	sData.cellSpikeTimesCortex = doSpikeIntToDouble(sData.cellSpikeTimesCortex,dblStepSize,sData.vecOverallT(1));
	fprintf('   Cortex spikes have been transformed... [%s]\n',getTime);
	if boolLoadLGN && isfield(sData,'cellSpikeTimesLGN_ON')
		sData.cellSpikeTimesLGN_ON = doSpikeIntToDouble(sData.cellSpikeTimesLGN_ON,dblStepSize,sData.vecOverallT(1));
		fprintf('   LGN ON spikes have been transformed... [%s]\n',getTime);
		sData.cellSpikeTimesLGN_OFF = doSpikeIntToDouble(sData.cellSpikeTimesLGN_OFF,dblStepSize,sData.vecOverallT(1));
		fprintf('   LGN OFF spikes have been transformed... [%s]\n',getTime);
	end
end

%% further processing
if isfield(sData,'matStimTypeCombos') && exist('vecOverallT','var')
	fprintf('Detected previous recording...\n')
	if strcmp(strConnFile,sData.strConnFile)
		fprintf('\b Connectivity structures match; starting concatenation [%s]\n',getTime)
	else
		fprintf('Error: connectivity identifiers do not match! <%s> <%s> [%s]\n',strConnFile,sData.strConnFile,getTime);
		error
	end
	%% change previous trial data to backup
	vecTrialOris1=vecTrialOris; 
	vecTrialOriIdx1=vecTrialOriIdx; 
	vecTrialSFs1=vecTrialSFs; 
	vecTrialSFIdx1=vecTrialSFIdx; 
	vecTrialTFs1=vecTrialTFs; 
	vecTrialTFIdx1=vecTrialTFIdx; 
	vecTrialContrasts1=vecTrialContrasts; 
	vecTrialContrastIdx1=vecTrialContrastIdx; 
	vecTrialLuminance1=vecTrialLuminance;
	vecTrialLuminanceIdx1=vecTrialLuminanceIdx;
	matStimTypeCombos1=matStimTypeCombos; 
	vecStimTypeOris1=vecStimTypeOris; 
	vecStimTypeSFs1=vecStimTypeSFs; 
	vecStimTypeTFs1=vecStimTypeTFs;
	vecStimTypeContrasts1 = vecStimTypeContrasts;
	vecStimTypeLuminance1 = vecStimTypeLuminance;
	
	% trial timing
	vecTrialStimType1 = vecTrialStimType;
	vecTrialStimRep1 = vecTrialStimRep;
	vecTrialStartSecs1 = vecTrialStartSecs;
	vecStimStartSecs1 = vecStimStartSecs;
	vecStimStopSecs1 = vecStimStopSecs;
	vecTrialEndSecs1 = vecTrialEndSecs;
	
	%spike counts
	matModelResp1 = matModelResp;
	if boolLoadLGN
		matModelRespLGN_ON1 = matModelRespLGN_ON;
		matModelRespLGN_OFF1 = matModelRespLGN_OFF;
	end
	clearvars matModelResp matModelRespLGN_ON matModelRespLGN_OFF
	
	%% check stim types and add types
	vecTrialOris=sData.vecTrialOris; 
	vecTrialOriIdx=sData.vecTrialOriIdx; 
	vecTrialSFs=sData.vecTrialSFs; 
	vecTrialSFIdx=sData.vecTrialSFIdx; 
	vecTrialTFs=sData.vecTrialTFs; 
	vecTrialTFIdx=sData.vecTrialTFIdx; 
	vecTrialContrasts=sData.vecTrialContrasts; 
	vecTrialContrastIdx=sData.vecTrialContrastIdx; 
	vecTrialLuminance=sData.vecTrialLuminance;
	vecTrialLuminanceIdx=sData.vecTrialLuminanceIdx;
	matStimTypeCombos=sData.matStimTypeCombos; 
	vecStimTypeOris=sData.vecStimTypeOris; 
	vecStimTypeSFs=sData.vecStimTypeSFs; 
	vecStimTypeTFs=sData.vecStimTypeTFs;
	vecStimTypeContrasts = sData.vecStimTypeContrasts;
	vecStimTypeLuminance = sData.vecStimTypeLuminance;
	
	%% trial identifiers
	vecTrialStartSecs = sData.vecTrialStartSecs;
	vecStimStartSecs = sData.vecStimStartSecs;
	vecStimStopSecs = sData.vecStimStopSecs;
	vecTrialEndSecs = sData.vecTrialEndSecs;
	vecTrialStimType = sData.vecTrialStimType;
	vecTrialStimRep = sData.vecTrialStimRep;
	
	%% spiking data
	cellSpikeTimesCortex = sData.cellSpikeTimesCortex;
	if boolLoadLGN && isfield(sData,'cellSpikeTimesLGN_ON')
		cellSpikeTimesLGN_ON = sData.cellSpikeTimesLGN_ON;
		cellSpikeTimesLGN_OFF = sData.cellSpikeTimesLGN_OFF;
	elseif isfield(sData,'cellSpikeTimesLGN_ON')
		sData = rmfield(sData,{'cellSpikeTimesLGN_ON','cellSpikeTimesLGN_OFF'});
	end
	
	
elseif isfield(sData,'matStimTypeCombos')
	fprintf('Detected type [matStimTypeCombos] \n')
	%% unpack
	%identifier
	strConnFile = sData.strConnFile;
	 
	%direct from run
	vecOverallT = sData.vecOverallT;
	dblDeltaT = sData.dblDeltaT;
	strStimType = sData.strStimType;
	if isfield(sData,'cellSpikeTimesCortex'),cellSpikeTimesCortex = sData.cellSpikeTimesCortex;end
	if boolLoadLGN && isfield(sData,'cellSpikeTimesLGN_ON')
		cellSpikeTimesLGN_ON = sData.cellSpikeTimesLGN_ON;
		cellSpikeTimesLGN_OFF = sData.cellSpikeTimesLGN_OFF;
	elseif isfield(sData,'cellSpikeTimesLGN_ON')
		sData = rmfield(sData,{'cellSpikeTimesLGN_ON','cellSpikeTimesLGN_OFF'});
	end

	vecTrialStartSecs = sData.vecTrialStartSecs;
	vecStimStartSecs = sData.vecStimStartSecs;
	vecStimStopSecs = sData.vecStimStopSecs;
	vecTrialEndSecs = sData.vecTrialEndSecs;
	vecTrialStimRep = sData.vecTrialStimRep;
	vecTrialStimType = sData.vecTrialStimType;
	
	%sConnParams
	if isfield(sData,'intCellsV1'),intCellsV1 = sData.intCellsV1;else intCellsV1 = 0;end
	if isfield(sData,'intCellsV2'),intCellsV2 = sData.intCellsV2;else intCellsV2 = 0;end
	if isfield(sData,'intCortexCells'),intCortexCells = sData.intCortexCells;else intCortexCells = intCellsV1+intCellsV2;end
	
	
	%sConnectivity
	%synapse-based
	matSynFromTo = sData.matSynFromTo;
	vecSynExcInh = sData.vecSynExcInh;
	vecSynDelay = sData.vecSynDelay;
	vecSynWeight = sData.vecSynWeight;
	vecSynConductance = sData.vecSynConductance;
	vecSynType = sData.vecSynType;
	%LGN
	if isfield(sData,'vecSynConductanceON_to_Cort')
		vecSynConductanceON_to_Cort = sData.vecSynConductanceON_to_Cort;
		vecSynConductanceOFF_to_Cort = sData.vecSynConductanceOFF_to_Cort;
		vecSynWeightON_to_Cort = sData.vecSynWeightON_to_Cort;
		vecSynWeightOFF_to_Cort = sData.vecSynWeightOFF_to_Cort;
		vecSynDelayON_to_Cort = sData.vecSynDelayON_to_Cort;
		vecSynDelayOFF_to_Cort = sData.vecSynDelayOFF_to_Cort;
		matSynConnON_to_Cort = sData.matSynConnON_to_Cort;
		matSynConnOFF_to_Cort = sData.matSynConnOFF_to_Cort;
	end
	
	%cell-based
	vecCellThresh = sData.vecCellThresh; 
	if isfield(sData,'vecTauPeakByType'),vecTauPeakByType = sData.vecTauPeakByType;else vecTauPeakByTypePeakByType=unique(sData.vecCellTauPeak);end
	vecCellTauPeak = sData.vecCellTauPeak; 
	vecCellV_E = sData.vecCellV_E; 
	vecCellV_I = sData.vecCellV_I; 
	vecCellV_AHP = sData.vecCellV_AHP; 
	vecCellV_Leak = sData.vecCellV_Leak; 
	vecCellCm = sData.vecCellCm; 
	vecCellG_Leak = sData.vecCellG_Leak; 
	vecCellG_AHP = sData.vecCellG_AHP; 
	vecCellArea = sData.vecCellArea; 
	vecCellTypes = sData.vecCellTypes; 
	vecPrefPsi = sData.vecPrefPsi; 
	vecPrefOri = sData.vecPrefOri; 
	vecPrefSF = sData.vecPrefSF; 
	vecPrefRF_X = sData.vecPrefRF_X; 
	vecPrefRF_Y = sData.vecPrefRF_Y; 
	dblSigmaX = sData.dblSigmaX; 
	dblSigmaY = sData.dblSigmaY; 
	if isfield(sData,'matPrefGabors')
		matPrefGabors = sData.matPrefGabors; 
		matFieldsV2 = sData.matFieldsV2;
	end
	
	%sStimInputs
	dblTrialDur=sData.dblTrialDur; 
	dblVisSpacing=sData.dblVisSpacing; 
	varDeltaSyn=sData.varDeltaSyn; 
	matBlankLGN_ON=sData.matBlankLGN_ON; 
	matBlankLGN_OFF=sData.matBlankLGN_OFF;
	if isfield(sData,'cellR_ON')
		cellR_ON=sData.cellR_ON; 
		cellR_OFF=sData.cellR_OFF; 
		cellLGN_ON=sData.cellLGN_ON; 
		cellLGN_OFF=sData.cellLGN_OFF; 
	end
	cellContrast=sData.cellContrast; 
	cellLuminance=sData.cellLuminance; 
	vecTrialOris=sData.vecTrialOris; 
	vecTrialOriIdx=sData.vecTrialOriIdx; 
	vecTrialSFs=sData.vecTrialSFs; 
	vecTrialSFIdx=sData.vecTrialSFIdx; 
	vecTrialTFs=sData.vecTrialTFs; 
	vecTrialTFIdx=sData.vecTrialTFIdx; 
	vecTrialContrasts=sData.vecTrialContrasts; 
	vecTrialContrastIdx=sData.vecTrialContrastIdx; 
	if isfield(sData,'vecTrialLuminance'),vecTrialLuminance=sData.vecTrialLuminance; end
	if isfield(sData,'vecTrialLuminanceIdx'),vecTrialLuminanceIdx=sData.vecTrialLuminanceIdx; end
	matStimTypeCombos=sData.matStimTypeCombos; 
	vecStimTypeOris=sData.vecStimTypeOris; 
	vecStimTypeSFs=sData.vecStimTypeSFs; 
	vecStimTypeTFs=sData.vecStimTypeTFs;
	vecStimTypeContrasts = sData.vecStimTypeContrasts;
	if isfield(sData,'vecStimTypeLuminance'),vecStimTypeLuminance = sData.vecStimTypeLuminance;end
elseif strcmp(sData.strConnFile(1:7),'KohnExp')
	fprintf('Detected type [KohnExperiment] \n')
	
	%% put all fields into separate variables
	cellFields = fieldnames(sData);
	for intField=1:numel(cellFields);
		strExpression = strcat(cellFields{intField},' = sData.',cellFields{intField},';');
		eval(strExpression);
	end
else
	fprintf('Error: type not detected!\n');
	error
end
%clear sData;
%% get spiking data
%cellSpikeTimesCortex
if exist('sLoadPrePro','var');
	matModelResp = sLoadPrePro.matModelResp;
	clear sLoadPrePro;
end

intTrials = numel(vecTrialStimType);
vecBinsTime = sort([0 vecStimStartSecs vecStimStopSecs vecStimStopSecs(end)+0.1]);
%vecBinsDur = diff(vecBinsTime);
if ~exist('matModelResp','var')
	intNeurons = numel(cellSpikeTimesCortex);
	hTic = tic;
	matModelResp = zeros(intNeurons,intTrials,'uint16');
	
	for intNeuron=1:intNeurons
		[vecCounts,edges,bin] = histcounts(cellSpikeTimesCortex{intNeuron},vecBinsTime);
		matModelResp(intNeuron,:) = uint16(vecCounts(2:2:end)); %Counts, not Hz!
		
		if toc(hTic) > 5 || intNeuron==1
			hTic = tic;
			fprintf('Now at neuron %d/%d [%s]\n',intNeuron,intNeurons,getTime);
		end
	end
	
	%save
	fprintf('Processing completed; saving prepro data [%s]\n',getTime);
	save([strPath 'Simulation_' strSimulation '_prepro.mat'],'matModelResp');
end
intNeurons = size(matModelResp,1);
%get LGN info
if exist('sLoadPreProLGN','var');
	matModelRespLGN_ON = sLoadPreProLGN.matModelRespLGN_ON;
	matModelRespLGN_OFF = sLoadPreProLGN.matModelRespLGN_OFF;
	clear sLoadPreProLGN;
end
if boolLoadLGN && ~exist('matModelRespLGN_ON','var') && exist('cellSpikeTimesLGN_ON','var')
	intNeuronsLGN = numel(cellSpikeTimesLGN_ON);
	hTic = tic;
	matModelRespLGN_ON = zeros(intNeuronsLGN,intTrials,'uint16');
	matModelRespLGN_OFF = zeros(intNeuronsLGN,intTrials,'uint16');
	
	for intNeuron=1:intNeuronsLGN
		[vecCounts,edges,bin] = histcounts(cellSpikeTimesLGN_ON{intNeuron},vecBinsTime);
		matModelRespLGN_ON(intNeuron,:) = uint16(vecCounts(2:2:end)); %Counts, not Hz!

		[vecCounts,edges,bin] = histcounts(cellSpikeTimesLGN_OFF{intNeuron},vecBinsTime);
		matModelRespLGN_OFF(intNeuron,:) = uint16(vecCounts(2:2:end)); %Counts, not Hz!
		
		if toc(hTic) > 5 || intNeuron==1
			hTic = tic;
			fprintf('Now at LGN neuron %d/%d [%s]\n',intNeuron,intNeuronsLGN,getTime);
		end
	end
	
	%save
	fprintf('Processing completed; saving prepro data [%s]\n',getTime);
	save([strPath 'Simulation_' strSimulation '_preproLGN.mat'],'matModelRespLGN_ON','matModelRespLGN_OFF');
end

%get cortex sync info
if exist('sLoadPreProSyncCort','var')
	matMeanCoincCort = sLoadPreProSyncCort.matMeanCoincCort;
	matDiagCoincCort = sLoadPreProSyncCort.matDiagCoincCort;
	clear sLoadPreProSyncCort;
end
%perform cortex sync calc
if boolLoadSync && ~exist('matMeanCoincCort','var') && exist('cellSpikeTimesCortex','var')
	hTic = tic;
	%
	%cellSpikeTimes=cellSpikeTimesCortex;
	%vecWindow = [vecStimStartSecs(1) vecStimStopSecs(1)];
	
	%transform spike times to integers
	fprintf('Preparing coinc cort calc... [%s]\n',getTime);
	cellSpikeTimes = doSpikeDoubleToInt(cellSpikeTimesCortex,dblStepSize,0);
	vecStimStartInt = doSpikeDoubleToInt(vecStimStartSecs,dblStepSize,0);
	vecStimStopInt = doSpikeDoubleToInt(vecStimStopSecs,dblStepSize,0);

	matMeanCoincCort = nan(intNeurons,intTrials);
	matDiagCoincCort = nan(intNeurons,intTrials);
	for intTrial=1:intTrials
		vecWindow = [vecStimStartInt(intTrial) vecStimStopInt(intTrial)];
		matCoincidence = calcCoincidence(cellSpikeTimes,vecWindow);
		matMeanCoincCort(:,intTrial) = mean(matCoincidence);
		matDiagCoincCort(:,intTrial) = diag(matCoincidence);
		
		if toc(hTic) > 5 || intTrial==1
			hTic = tic;
			fprintf('Coinc cort calc; Now at trial %d/%d [%s]\n',intTrial,intTrials,getTime);
		end
	end
	save([strPath 'Simulation_' strSimulation '_preproSyncCort.mat'],'matMeanCoincCort','matDiagCoincCort');
end

%% concatenate spike count matrices
if exist('matModelResp1','var')
	
	vecTrialOris=cat(2,vecTrialOris1,vecTrialOris);
	vecTrialOriIdx=label2idx(vecTrialOris);
	vecTrialSFs=cat(2,vecTrialSFs1,vecTrialSFs);
	vecTrialSFIdx=label2idx(vecTrialSFs);
	vecTrialTFs=cat(2,vecTrialTFs1, vecTrialTFs);
	vecTrialTFIdx=label2idx(vecTrialTFs);
	vecTrialContrasts=cat(2,vecTrialContrasts1, vecTrialContrasts);
	vecTrialContrastIdx=label2idx(vecTrialContrasts);
	vecTrialLuminance=cat(2,vecTrialLuminance1,vecTrialLuminance);
	vecTrialLuminanceIdx=label2idx(vecTrialLuminance);
	vecStimTypeOris=cat(2,vecStimTypeOris1, vecStimTypeOris);
	vecStimTypeSFs=cat(2,vecStimTypeSFs1,vecStimTypeSFs);
	vecStimTypeTFs=cat(2,vecStimTypeTFs1,vecStimTypeTFs);
	vecStimTypeContrasts = cat(2,vecStimTypeContrasts1,vecStimTypeContrasts);
	vecStimTypeLuminance = cat(2,vecStimTypeLuminance1,vecStimTypeLuminance);
	%stim combos
	cellParamTypes = {vecTrialOriIdx,vecTrialSFIdx,vecTrialTFIdx,vecTrialContrastIdx,vecTrialLuminanceIdx};
	matStimTypeCombos = buildStimCombos(cellParamTypes);
	
	% trial timing & remove excess trial repetitions
	intMaxRep1 = max(vecTrialStimRep1);
	intMaxRep = max(vecTrialStimRep);
	intUseMaxRep = min([intMaxRep intMaxRep1]);
	indUseTrials1 = vecTrialStimRep1<=intUseMaxRep;
	indUseTrials = vecTrialStimRep<=intUseMaxRep;
	intTypes1 = numel(unique(vecTrialStimType1));
	dblEndTime1 = max(vecTrialEndSecs1(indUseTrials1));
	vecTrialStimType = cat(2,vecTrialStimType1(indUseTrials1),vecTrialStimType(indUseTrials)+intTypes1);
	vecTrialStimRep = cat(2,vecTrialStimRep1(indUseTrials1),vecTrialStimRep(indUseTrials));
	vecTrialStartSecs = cat(2,vecTrialStartSecs1(indUseTrials1),vecTrialStartSecs(indUseTrials)+dblEndTime1);
	vecStimStartSecs = cat(2,vecStimStartSecs1(indUseTrials1),vecStimStartSecs(indUseTrials)+dblEndTime1);
	vecStimStopSecs = cat(2,vecStimStopSecs1(indUseTrials1),vecStimStopSecs(indUseTrials)+dblEndTime1);
	vecTrialEndSecs = cat(2,vecTrialEndSecs1(indUseTrials1),vecTrialEndSecs(indUseTrials)+dblEndTime1);
	
	%spike counts
	matModelResp = cat(2,matModelResp1(:,indUseTrials1),matModelResp(:,indUseTrials));
end
if boolLoadLGN && exist('matModelRespLGN_ON1','var')
	%spike counts
	matModelRespLGN_ON = cat(2,matModelRespLGN_ON1(:,indUseTrials1),matModelRespLGN_ON(:,indUseTrials));
	matModelRespLGN_OFF = cat(2,matModelRespLGN_OFF1(:,indUseTrials1),matModelRespLGN_OFF(:,indUseTrials));
end
clearvars matModelResp1 matModelRespLGN_ON1 matModelRespLGN_OFF1;

%% save sData with removal of spike arrays
if ~indLoadedPrePros(end) && ~boolLoadSpikeTimes
	fprintf('Processing completed; saving prepro0 data [%s]\n',getTime);
	sData = reduceData(sData);
	save([strPath 'Simulation_' strSimulation '_prepro0.mat'],'sData');
end

%% build structStim
%{
structStim.Orientation = vecTrialOris;
structStim.TrialNumber = 1:length(vecTrialOris);
structStim.FrameOn = vecStimStartSecs;
structStim.FrameOff = vecStimStopSecs;
vecOrientations = unique(vecTrialOris);
vecOriDegs = rad2deg(vecOrientations);
sTypes = getStimulusTypes(structStim);
cellSelect = getSelectionVectors(structStim,sTypes);
%}
