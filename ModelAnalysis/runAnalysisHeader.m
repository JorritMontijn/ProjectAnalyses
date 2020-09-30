%% different pre-processing!
%{
We counted spikes in 100 ms bins, beginning 160 ms after stimulus onset and spanning a total of
1 s (10 bins per trial). Because our focus is on trial-to-trial variability, we subtracted the peristimulus
time histogram (PSTH) for each neuron and condition, from the binned spike counts,
yielding the residual activity for each neuron on each trial. Neurons which fired less than 0.5
spikes/s were excluded from the analyses. We compared our analysis of V1-V2 interactions to
the results of applying the same analyses to a held-out V1 population (V1-V1). The source
V1 population for the V1-V2 and V1-V1 analyses was identical. The target population in
the V1-V1 analyses was a held-out subset of the originally recorded population, which was
matched in neuron count to the corresponding V2 population. We also matched the firing
rate distribution (mean-matched) to the V2 population separately for each stimulus condition
(as in Churchland et al. [?]). Briefly, we binned the firing rate distribution of the V1 and
V2 populations (for each neuron, the average firing rate was taken across time and trials for
each dataset), and determined the common firing rate distribution (i.e. for each firing rate
bin, we took the minimum neuron count between the two populations). For each firing rate
interval, we then randomly picked this minimum number of neurons from the corresponding bin
in each population, without replacement. Because we had many more V1 than V2 neurons, the
common distribution usually matched the V2 distribution and we selected an equal number of
V1 neurons. The size of the matched populations ranged from 15 to 31 units across all datasets
(mean: 22.3). The V1 neurons that were not selected for the held-out population defined the
source V1 population. V2 neurons that were not selected for the V2 mean-matched population
were not used in the analysis. We repeated the mean- matching procedure 25 times, using
di?erent random, mean-matched subsets of neurons (and consequently producing a di?erent
source population). Results were based on averages across these repeats.
%}

%runModelHeader
%error([mfilename ':Obsolete'],'Please use version 2');

%% load data
indRequiredPrePros = [true; ... %prepro
	false;... %prepro2
	true;...%preproLGN
	false;...%preproSync
	false;...%preprotrialsplit
	true]; %prepro0
indLoadedPrePros = false(size(indRequiredPrePros));
indRequiredPrePros(3) = boolLoadLGN;
indRequiredPrePros(4) = boolLoadSync;
indRequiredPrePros(5) = boolLoadTrialSplit;
if ~exist('boolLoadSpikeTimes','var')
	boolLoadSpikeTimes = false;
end
if strcmpi(strSimulation(end-3:end),'.mat'),error([mfilename ':ExtensionSupplied'],'Filename includes extension; please remove ".mat"');end
if (~exist('sData','var') && ~exist('sSimRun','var')) || ~exist('strLastSim','var') || ~strcmp(strLastSim,strSimulation)
	fprintf('Loading simulation <%s> [%s]\n',strSimulation,getTime);
	%clearvars -except intLoadSim strSimulation strFigDir bool*
	if ~exist('boolLoadLGN','var'),boolLoadLGN=false;end
	strLastSim = strSimulation;
	strPath = 'F:\Data\Results\SimAggregates\';
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
	if boolLoadTrialSplit && exist([strPath 'Simulation_' strSimulation '_preproTrialSplit.mat'],'file')
		indLoadedPrePros(5) = true;
		sLoadPreProTrialSplit = load([strPath 'Simulation_' strSimulation '_preproTrialSplit.mat']);%matModelRespTS
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

% load connectivity separately
if ~isfield(sData,'matSynFromTo')
	load(['F:\Code\Simulations\SimulationsEVS\Connectivity\' sData.strConnFile]);
	sData = catstruct(sData,sConnectivity');
end
	
%% further processing
if isfield(sData,'matSynFromTo') && exist('vecOverallT','var')
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
	cellFields = fieldnames(sData);
	for intField=1:numel(cellFields);
		strField = cellFields{intField};
		if strcmpi(strField(1:length('vecTrial')),'vecTrial') || strcmpi(strField(1:length('vecStimType')),'vecStimType')
			strExpression = strcat(strField,' = sData.',strField,';');
			eval(strExpression);
		end
	end
	matStimTypeCombos=sData.matStimTypeCombos;
	
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
	
	
elseif isfield(sData,'matSynFromTo')
	fprintf('Detected type [matSynFromTo] \n')
	%% unpack
	%identifier
	strConnFile = sData.strConnFile;
	
	%% check stim types and add types
	cellFields = fieldnames(sData);
	for intField=1:numel(cellFields);
		strField = cellFields{intField};
		cellUse = {'vecTrial','vecStimType'};
		for intUse=1:numel(cellUse)
			strUse = cellUse{intUse};
			intLength = length(strUse);
			if length(strField) > intLength && strcmpi(strField(1:intLength),strUse)
				strExpression = strcat(strField,' = sData.',strField,';');
				eval(strExpression);
			end
		end
	end
	matStimTypeCombos=sData.matStimTypeCombos;
	
	%% trial identifiers
	vecTrialStartSecs = sData.vecTrialStartSecs;
	vecStimStartSecs = sData.vecStimStartSecs;
	vecStimStopSecs = sData.vecStimStopSecs;
	vecTrialEndSecs = sData.vecTrialEndSecs;
	vecTrialStimType = sData.vecTrialStimType;
	vecTrialStimRep = sData.vecTrialStimRep;
	
	%% spiking data
	if isfield(sData,'cellSpikeTimesCortex'),cellSpikeTimesCortex = sData.cellSpikeTimesCortex;end
	if boolLoadLGN && isfield(sData,'cellSpikeTimesLGN_ON')
		cellSpikeTimesLGN_ON = sData.cellSpikeTimesLGN_ON;
		cellSpikeTimesLGN_OFF = sData.cellSpikeTimesLGN_OFF;
	elseif isfield(sData,'cellSpikeTimesLGN_ON')
		sData = rmfield(sData,{'cellSpikeTimesLGN_ON','cellSpikeTimesLGN_OFF'});
	end
	
	%sConnParams
	if isfield(sData,'intCellsV1'),intCellsV1 = sData.intCellsV1;else intCellsV1 = 0;end
	if isfield(sData,'intCellsV2'),intCellsV2 = sData.intCellsV2;else intCellsV2 = 0;end
	if isfield(sData,'intCortexCells'),intCortexCells = sData.intCortexCells;else intCortexCells = intCellsV1+intCellsV2;end
	
	
	%% sConnectivity
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
	
	%% sStimInputs
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
	
elseif strcmp(sData.strConnFile(1:7),'KohnExp')
	fprintf('Detected type [KohnExperiment] \n')
	
	%% put all fields into separate variables
	cellFields = fieldnames(sData);
	for intField=1:numel(cellFields);
		strExpression = strcat(cellFields{intField},' = sData.',cellFields{intField},';');
		eval(strExpression);
	end
elseif strcmp(sData.strConnFile(1:8),'RubenExp')
	fprintf('Detected type [RubenExperiment] \n')
	
	%% put all fields into separate variables
	cellFields = fieldnames(sData);
	for intField=1:numel(cellFields);
		strExpression = strcat(cellFields{intField},' = sData.',cellFields{intField},';');
		eval(strExpression);
	end
elseif strcmp(sData.strConnFile(1:7),'AverExp')
	fprintf('Detected type [AverbeckExperiment] \n')
	
	%% put all fields into separate variables
	cellFields = fieldnames(sData);
	for intField=1:numel(cellFields);
		strExpression = strcat(cellFields{intField},' = sData.',cellFields{intField},';');
		eval(strExpression);
	end
elseif strcmp(sData.strConnFile(1:7),'JoGuExp')
	fprintf('Detected type [JoGuExperiment] \n')
	
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
%% get spiking data per trial, split in 100ms blocks
if boolLoadTrialSplit
	%% get spiking data per trial
	if exist('sLoadPreProTrialSplit','var')
		if isfield(sLoadPreProTrialSplit,'matModelRespTS')
			matModelRespTS = sLoadPreProTrialSplit.matModelRespTS;
		end
		if isfield(sLoadPreProTrialSplit,'matModelRespNormTS3')
			matModelRespNormTS3 = sLoadPreProTrialSplit.matModelRespNormTS3;
		end
		if isfield(sLoadPreProTrialSplit,'matModelRespTS3')
			matModelRespTS3 = sLoadPreProTrialSplit.matModelRespTS3;
		end
		clear sLoadPreProTrialSplit;
	end
	%{
	We counted spikes in 100 ms bins, beginning 160 ms after stimulus onset
	and spanning a total of 1 s (10 bins per trial). Because our focus is
	on trial-to-trial variability, we subtracted the peristimulus time
	histogram (PSTH) for each neuron and condition, from the binned spike
	counts, yielding the residual activity for each neuron on each trial
	%}
	
	%get time stamps for bins
	intNeurons = intCortexCells;
	intTrials = numel(vecTrialStimType);
	dblStimDur = unique(roundi(vecStimStopSecs-vecStimStartSecs,3));
	dblBinSizeSecs = 0.1; %100 ms
	if dblStimDur > 0.5
		dblOffsetInitial = 0.16; %160 ms for exp
		intBinsPerTrial = 10;
	else
		dblOffsetInitial = 0; %0 ms for sim
		intBinsPerTrial = (dblStimDur-dblOffsetInitial)/dblBinSizeSecs;
	end
	if ~isint(intBinsPerTrial),error([mfilename ':BinsNotInteger'],'Number of bins per trial is not an integer');end
	vecTrialBinStartTimes = dblOffsetInitial:dblBinSizeSecs:(dblOffsetInitial+intBinsPerTrial*dblBinSizeSecs);
	matT = (bsxfun(@plus,vecStimStartSecs,vecTrialBinStartTimes'));
	vecBinsTime = sort([0 matT(:)' vecStimStopSecs(end)+0.1]);
	vecRemove = 1:(intBinsPerTrial+1):numel(vecBinsTime);
	intTimePoints = intBinsPerTrial*intTrials;
	indKeep = true(1,intTimePoints + intTrials + 1);
	indKeep(vecRemove) = false;
	%trial indicator
	vecTrialIdxTS = reshape(ones(intBinsPerTrial,1)*(1:intTrials),[1 intTimePoints]);
	matTrialIdxTS3 = reshape(repmat(vecTrialIdxTS,[intNeurons 1]),[intNeurons intBinsPerTrial intTrials]);
	%bindicator
	vecTrialBinTS = reshape((1:intBinsPerTrial)'*ones(1,intTrials),[1 intTimePoints]);
	matTrialBinTS3 = reshape(repmat(vecTrialBinTS,[intNeurons 1]),[intNeurons intBinsPerTrial intTrials]);
	
	%run
	if ~exist('matModelRespTS','var')
		fprintf('Running trial-split processing... [%s]\n',getTime);
		hTic = tic;
		matModelRespTS = zeros(intNeurons,intTimePoints,'int8');
		
		for intNeuron=1:intNeurons
			[vecCounts,edges,bin] = histcounts(cellSpikeTimesCortex{intNeuron},vecBinsTime);
			matModelRespTS(intNeuron,:) = int8(vecCounts(indKeep)); %Counts, not Hz!
			
			if toc(hTic) > 5 || intNeuron==1
				hTic = tic;
				fprintf('Now at neuron %d/%d [%s]\n',intNeuron,intNeurons,getTime);
			end
		end
	end
	
	%get 3D
	if ~exist('matModelRespNormTS3','var') || ~exist('matModelRespTS3','var')
		matModelRespTS3 = reshape(matModelRespTS,[intNeurons intBinsPerTrial intTrials]);
		matModelRespNormTS3 = prepTrialData(matModelRespTS3,vecTrialStimType);
		%save
		fprintf('Trial-split processing completed; saving prepro data [%s]\n',getTime);
		save([strPath 'Simulation_' strSimulation '_preproTrialSplit.mat'],'matModelRespTS','matModelRespTS3','matModelRespNormTS3');
	end
end

%% get spiking data per trial
if exist('sLoadPrePro','var');
	matModelResp = sLoadPrePro.matModelResp;
	clear sLoadPrePro;
end

intTrials = numel(vecTrialStimType);
vecBinsTime = sort([0 vecStimStartSecs vecStimStopSecs vecStimStopSecs(end)+0.1]);
%vecBinsDur = diff(vecBinsTime);
if ~exist('matModelResp','var')
	fprintf('Running trial-mean processing... [%s]\n',getTime);
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

%% get LGN spikes
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
	fprintf('Running LGN processing... [%s]\n',getTime);
	
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

%% sync stuff
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

%% prep variables for use in main script
cellIn = strsplit(strSimulation,'_');
strFiller = cellIn{1};
strType = cell2mat(cellIn(2:(end-1)));
strDate = cellIn{end};
strTag = [strType '_' strDate];
if boolDoSplitAnalysis
	strTag = ['TS_' strTag];
	strType = ['TS_' strType];
end

%% calculate which parameter to use
vecRanges = range(matStimTypeCombos(2:end,:),2);
[d,intUseParam]=max(vecRanges); %exclude orientation
intUseParam = intUseParam + 1;
if all(vecRanges==0), intUseParam = 0;vecParamStimType = [0];strParam='None';
	vecParamIdx = ones(1,size(matStimTypeCombos,2));
	vecParamVals = 1;
	intParamNum=1;
	mapC = redbluepurple(intParamNum);
else
	if intUseParam == 1,vecParamStimType = vecStimTypeOris;strParam='Ori';
	elseif intUseParam == 2,vecParamStimType = vecStimTypeSFs;strParam='SF';
	elseif intUseParam == 3,vecParamStimType = vecStimTypeTFs;strParam='TF';
	elseif intUseParam == 4,vecParamStimType = vecStimTypeContrasts;strParam='Contrast';
	elseif intUseParam == 5,vecParamStimType = vecStimTypeLuminance;strParam='Luminance';
	end
	vecParamIdx = matStimTypeCombos(intUseParam,:);
	vecParamVals = unique(vecParamStimType);
	intParamNum=numel(vecParamVals);
	mapC = redbluepurple(intParamNum);
end


%% build comparison matrix for stim types
vecUniqueStimTypes = unique(vecTrialStimType);
intMultiComp = false;
intStimTypes = numel(vecUniqueStimTypes);
if ~exist('matCompareTypes','var')
	if sum(range(matStimTypeCombos,2)>0)>1
		intMultiComp = true;
		matCompareTypes = [ones(1,size(matStimTypeCombos,2)-1)' (2:size(matStimTypeCombos,2))'];
		matCompareTypes = [1 2; 1 3; 1 4];
	else
		
		matCompareTypes = [(1:intStimTypes)' circshift((1:intStimTypes)',-1)];
	end
end

%% get data
if boolLoadTrialSplit
	fprintf('Transforming data [%s]\n',getTime);
	[vecTrialStimTypeUnsplit,matRespUnsplit,vecTrialStimTypeSplit,matRespSplitNorm,matRespSplit,matSplitTrialIdx] = ...
		prepSplitData(matModelRespTS3,vecTrialStimType);
	
	if boolDoSplitAnalysis
		matData = matRespSplitNorm;
		matDataFull = matRespSplit;
		vecTrialStimType = vecTrialStimTypeSplit;
	else
		matData = matRespUnsplit;
		vecTrialStimType = vecTrialStimTypeUnsplit;
	end
else
	matData = double(matModelResp);
	vecTrialStimType = vecTrialStimType;
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
