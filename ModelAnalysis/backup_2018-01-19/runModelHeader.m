%runModelHeader
%% load data
if (~exist('sData','var') && ~exist('sSimRun','var')) || ~exist('strLastSim','var') || ~strcmp(strLastSim,strSimulation)
	fprintf('Loading simulation <%s> [%s]\n',strSimulation,getTime);
	%clearvars -except intLoadSim strSimulation strFigDir bool*
	if ~exist('boolLoadLGN','var'),boolLoadLGN=false;end
	strLastSim = strSimulation;
	strPath = 'A:\SimAggregates\';
	if exist([strPath 'Simulation_' strSimulation '_prepro.mat'],'file')
		load([strPath 'Simulation_' strSimulation '_prepro.mat']);
	end
	if exist([strPath 'Simulation_' strSimulation '_prepro2.mat'],'file')
		load([strPath 'Simulation_' strSimulation '_prepro2.mat']);
	end
	if boolLoadLGN && exist([strPath 'Simulation_' strSimulation '_preproLGN.mat'],'file')
		load([strPath 'Simulation_' strSimulation '_preproLGN.mat']);
	end
	sLoad = load([strPath 'Simulation_' strSimulation '.mat']);
	
end

%% check if vars are saved as fields
if exist('sData','var') %do nothing
elseif isfield(sLoad,'sData')
	sData = sLoad.sData;
	clear sLoad;
elseif isfield(sLoad,'sParallel')
	sParallel = sLoad.sParallel;
	sParams = sLoad.sParams;
	sSimRun = sLoad.sSimRun;
elseif isfield(sLoad,'vecOverallT')
	sData = sLoad;
	clear sLoad;
end

%% unpack
if exist('sParallel','var') %parallel data; reformat
	sData = catstruct(rmfield(sParams,{'sStimParams','sStimInputs','sConnParams','sConnectivity'}),sParams.sStimParams,sParams.sStimInputs,sParams.sConnParams,sParams.sConnectivity);
	sData.vecOverallT = sort(sParallel.matOverallT(:))';
	sData.vecThisV = 'parallel';
	
	%% remove nans
	vecSpikeCounterCortex = zeros(size(sSimRun(1).cellSpikeTimesCortex));
	vecSpikeCounterLGN_ON = zeros(size(sSimRun(1).cellSpikeTimesLGN_ON));
	vecSpikeCounterLGN_OFF = zeros(size(sSimRun(1).cellSpikeTimesLGN_OFF));
	
	hTic=[];
	for intRepRun = 1:numel(sSimRun)
		if isempty(hTic) || toc(hTic) > 5,fprintf('Removing nans for repetition run %d/%d [%s]\n',intRepRun,numel(sSimRun),getTime);drawnow;hTic=tic;end
		
		% remove placeholder spike entries cortex
		cellSTCT = sSimRun(intRepRun).cellSpikeTimesCortex;
		for intN=1:numel(cellSTCT)
			if numel(cellSTCT{intN}) > 1
				cellSTCT{intN} = cellSTCT{intN}(~isnan(cellSTCT{intN}));
			else
				cellSTCT{intN} = [];
			end
		end
		sSimRun(intRepRun).cellSpikeTimesCortex = cellSTCT;
		vecSpikeCounterCortex = vecSpikeCounterCortex + cellfun(@numel,sSimRun(intRepRun).cellSpikeTimesCortex);
		
		% remove placeholder spike entries LGN ON
		cellSTLNT = sSimRun(intRepRun).cellSpikeTimesLGN_ON;
		for intN=1:numel(cellSTLNT)
			if numel(cellSTLNT{intN}) > 1
				cellSTLNT{intN} = cellSTLNT{intN}(~isnan(cellSTLNT{intN}));
			else
				cellSTLNT{intN} = [];
			end
		end
		sSimRun(intRepRun).cellSpikeTimesLGN_ON = cellSTLNT;
		vecSpikeCounterLGN_ON = vecSpikeCounterLGN_ON + cellfun(@numel,sSimRun(intRepRun).cellSpikeTimesLGN_ON);
		
		% remove placeholder spike entries LGN OFF
		cellSTLFT = sSimRun(intRepRun).cellSpikeTimesLGN_OFF;
		for intN=1:numel(cellSTLFT)
			if numel(cellSTLFT{intN}) > 1
				cellSTLFT{intN} = cellSTLFT{intN}(~isnan(cellSTLFT{intN}));
			else
				cellSTLFT{intN} = [];
			end
		end
		sSimRun(intRepRun).cellSpikeTimesLGN_OFF = cellSTLFT;
		vecSpikeCounterLGN_OFF = vecSpikeCounterLGN_OFF + cellfun(@numel,sSimRun(intRepRun).cellSpikeTimesLGN_OFF);
		
	end
	
	%% concatenate spiking arrays
	cellSpikeTimesCortex = cell(size(sSimRun(1).cellSpikeTimesCortex));
	cellSpikeTimesLGN_ON = cell(size(sSimRun(1).cellSpikeTimesLGN_ON));
	cellSpikeTimesLGN_OFF = cell(size(sSimRun(1).cellSpikeTimesLGN_OFF));
	
	hTic=[];
	for intN = 1:numel(sSimRun(1).cellSpikeTimesCortex)
		if isempty(hTic) || toc(hTic) > 5,fprintf('Concatenating cortex neuron %d/%d [%s]\n',intN,numel(sSimRun(1).cellSpikeTimesCortex),getTime);drawnow;hTic=tic;end
		
		%pre-allocate
		cellSpikeTimesCortex{intN} = nan(1,vecSpikeCounterCortex(intN));
		intC = 0;
		
		for intRepRun = 1:numel(sSimRun)
			vecSpikes = sSimRun(intRepRun).cellSpikeTimesCortex{intN};
			intS = numel(vecSpikes);
			cellSpikeTimesCortex{intN}((intC+1):(intC+intS)) = vecSpikes;
			intC = intC + intS;
		end
	end
	for intN = 1:numel(sSimRun(1).cellSpikeTimesLGN_ON)
		if isempty(hTic) || toc(hTic) > 5,fprintf('Concatenating LGN ON neuron %d/%d [%s]\n',intN,numel(sSimRun(1).cellSpikeTimesLGN_ON),getTime);drawnow;hTic=tic;end
		
		%pre-allocate
		cellSpikeTimesLGN_ON{intN} = nan(1,vecSpikeCounterLGN_ON(intN));
		intC = 0;
		
		for intRepRun = 1:numel(sSimRun)
			vecSpikes = sSimRun(intRepRun).cellSpikeTimesLGN_ON{intN};
			intS = numel(vecSpikes);
			cellSpikeTimesLGN_ON{intN}((intC+1):(intC+intS)) = vecSpikes;
			intC = intC + intS;
		end
	end
	for intN = 1:numel(sSimRun(1).cellSpikeTimesLGN_OFF)
		if isempty(hTic) || toc(hTic) > 5,fprintf('Concatenating LGN OFF neuron %d/%d [%s]\n',intN,numel(sSimRun(1).cellSpikeTimesLGN_OFF),getTime);drawnow;hTic=tic;end
		
		%pre-allocate
		cellSpikeTimesLGN_OFF{intN} = nan(1,vecSpikeCounterLGN_OFF(intN));
		intC = 0;
		
		for intRepRun = 1:numel(sSimRun)
			vecSpikes = sSimRun(intRepRun).cellSpikeTimesLGN_OFF{intN};
			intS = numel(vecSpikes);
			cellSpikeTimesLGN_OFF{intN}((intC+1):(intC+intS)) = vecSpikes;
			intC = intC + intS;
		end
	end
	
	%put in sData
	sData.cellSpikeTimesLGN_ON = cellSpikeTimesLGN_ON;
	sData.cellSpikeTimesLGN_OFF = cellSpikeTimesLGN_OFF;
	sData.cellSpikeTimesCortex = cellSpikeTimesCortex;
	sData.vecSpikeCounterLGN_ON = vecSpikeCounterLGN_ON;
	sData.vecSpikeCounterLGN_OFF = vecSpikeCounterLGN_OFF;
	sData.vecSpikeCounterCortex = vecSpikeCounterCortex;
	
	%% transform cellSpikes to uint32 arrays
	dblStepSize = (sData.vecOverallT(end)-sData.vecOverallT(1))/(numel(sData.vecOverallT)-1);
	vecUniqueSteps = unique(diff(sData.vecOverallT));
	vecUniqueOffsets = abs(vecUniqueSteps-dblStepSize);
	if any(vecUniqueOffsets > (dblStepSize/1000))
		warning([mfilename ':StepSizeError'],'Step size offset is larger than acceptable tolerance: %f',max(vecUniqueOffsets));
	end
	fprintf('Transforming Doubles to Int32 timestamps; maximum step-error is %e with step size %e [%s]\n',max(vecUniqueOffsets),dblStepSize,getTime);
	sData.cellSpikeTimesCortex = doSpikeDoubleToInt(sData.cellSpikeTimesCortex,dblStepSize,sData.vecOverallT(1));
	fprintf('   Cortex spikes have been transformed... [%s]\n',getTime);
	sData.cellSpikeTimesLGN_ON = doSpikeDoubleToInt(sData.cellSpikeTimesLGN_ON,dblStepSize,sData.vecOverallT(1));
	fprintf('   LGN ON spikes have been transformed... [%s]\n',getTime);
	sData.cellSpikeTimesLGN_OFF = doSpikeDoubleToInt(sData.cellSpikeTimesLGN_OFF,dblStepSize,sData.vecOverallT(1));
	fprintf('   LGN OFF spikes have been transformed... [%s]\n',getTime);

	%% save
	strTargetDir = 'D:\Data\Processed\V1_LIFmodel\';
	strOutputFile = ['Simulation_' strSimulation '.mat'];
	fprintf('Processing of parallel structure complete, saving data to <%s%s> [%s]\n',strTargetDir,strOutputFile,getTime);
	save([strTargetDir strOutputFile],'-struct','sData','-v7');
	fprintf('Done! [%s]\n',getTime);
end
%check if cellSpike arrays are uint32
if exist('sData','var') && ~isa(sData.cellSpikeTimesCortex{1},'double')
	dblStepSize = (sData.vecOverallT(end)-sData.vecOverallT(1))/(numel(sData.vecOverallT)-1);
	sData.cellSpikeTimesCortex = doSpikeIntToDouble(sData.cellSpikeTimesCortex,dblStepSize,sData.vecOverallT(1));
	fprintf('   Cortex spikes have been transformed... [%s]\n',getTime);
	if isfield(sData,'cellSpikeTimesLGN_ON')
		sData.cellSpikeTimesLGN_ON = doSpikeIntToDouble(sData.cellSpikeTimesLGN_ON,dblStepSize,sData.vecOverallT(1));
		fprintf('   LGN ON spikes have been transformed... [%s]\n',getTime);
		sData.cellSpikeTimesLGN_OFF = doSpikeIntToDouble(sData.cellSpikeTimesLGN_OFF,dblStepSize,sData.vecOverallT(1));
		fprintf('   LGN OFF spikes have been transformed... [%s]\n',getTime);
	end
end

%% further processing
if isfield(sData,'matCortConn') %type 1
	vecOverallT = sData.vecOverallT;
	dblDeltaT = sData.dblDeltaT;
	matCortConn = sData.matCortConn;
	dblSynSpikeMem = sData.dblSynSpikeMem;
	vecCortSynType = sData.vecCortSynType;
	intCortexCells = sData.intCortexCells;
	vecCortDelay = sData.vecCortDelay;
	vecCortConductance = sData.vecCortConductance;
	vecCellThresh = sData.vecCellThresh;
	vecTauPeakByType = sData.vecTauPeakByType;
	vecCellV_E = sData.vecCellV_E;
	vecCellV_I = sData.vecCellV_I;
	vecCellV_AHP = sData.vecCellV_AHP;
	vecCellV_Leak = sData.vecCellV_Leak;
	vecCellCm = sData.vecCellCm;
	vecCellG_Leak = sData.vecCellG_Leak;
	vecCellG_AHP = sData.vecCellG_AHP;
	vecSynConductanceON_to_Cort = sData.vecSynConductanceON_to_Cort;
	vecSynConductanceOFF_to_Cort = sData.vecSynConductanceOFF_to_Cort;
	vecSynWeightON_to_Cort = sData.vecSynWeightON_to_Cort;
	vecSynWeightOFF_to_Cort = sData.vecSynWeightOFF_to_Cort;
	vecSynDelayON_to_Cort = sData.vecSynDelayON_to_Cort;
	vecSynDelayOFF_to_Cort = sData.vecSynDelayOFF_to_Cort;
	matSynConnON_to_Cort = sData.matSynConnON_to_Cort;
	matSynConnOFF_to_Cort = sData.matSynConnOFF_to_Cort;
	matBlankLGN_ON = sData.matBlankLGN_ON;
	matBlankLGN_OFF = sData.matBlankLGN_OFF;
	cellLGN_ON = sData.cellLGN_ON;
	cellLGN_OFF = sData.cellLGN_OFF;
	vecTrialOris = sData.vecTrialOris;
	vecTrialOriIdx = sData.vecTrialOriIdx;
	vecStimStartSecs = sData.vecStimStartSecs+0.25;
	vecTrialEndSecs = sData.vecTrialEndSecs;
	vecTrialStimType = sData.vecTrialStimType;
	vecThisV = sData.vecThisV;
	boolStimPresent = sData.boolStimPresent;
	intPrevTrial = sData.intPrevTrial;
	intTrialT = sData.intTrialT;
	intIter = sData.intIter;
	cellSpikeTimesLGN_ON = sData.cellSpikeTimesLGN_ON;
	cellSpikeTimesLGN_OFF = sData.cellSpikeTimesLGN_OFF;
	cellSpikeTimesCortex = sData.cellSpikeTimesCortex;
	vecSpikeCounterLGN_ON = sData.vecSpikeCounterLGN_ON;
	vecSpikeCounterLGN_OFF = sData.vecSpikeCounterLGN_OFF;
	vecSpikeCounterCortex = sData.vecSpikeCounterCortex;
	intPreAllocationSize = sData.intPreAllocationSize;
	vecCellTypes = sData.vecCellTypes;
	dblStimDur = sData.dblStimDur;
	vecStimStopSecs = vecStimStartSecs+dblStimDur;
	
	vecCellArea = ones(size(vecCellTypes));
	vecPrefPsi = sData.vecPrefPsi;
	vecPrefOri = sData.vecPrefOri;
	vecPrefSF = sData.vecPrefOri;
	vecPrefRF_X = sData.vecPrefRF_X;
	vecPrefRF_Y = sData.vecPrefRF_Y;
	
	vecSynExcInh = vecCortSynType;
	matSynFromTo = matCortConn;
	vecSynConductance = vecCortConductance;
	vecSynWeight = ones(size(vecCortConductance));

elseif isfield(sData,'sStimInputs') && isfield(sData,'sConnParams')
	%% unpack
	%direct from run
	vecOverallT = sData.vecOverallT;
	dblDeltaT = sData.dblDeltaT;
	
	sData = rmfield(sData,{'cellSpikeTimesLGN_ON','cellSpikeTimesLGN_OFF'});
	cellSpikeTimesCortex = sData.cellSpikeTimesCortex;
	
	vecTrialStartSecs = sData.vecTrialStartSecs;
	vecStimStartSecs = sData.vecStimStartSecs;
	vecStimStopSecs = sData.vecStimStopSecs;
	vecTrialEndSecs = sData.vecTrialEndSecs;
	vecTrialStimRep = sData.vecTrialStimRep;
	vecTrialStimType = sData.vecTrialStimType;
	
	%sConnParams
	if isfield(sData.sConnParams,'intCellsV1'),intCellsV1 = sData.sConnParams.intCellsV1;else intCellsV1 = 0;end
	if isfield(sData.sConnParams,'intCellsV2'),intCellsV2 = sData.sConnParams.intCellsV2;else intCellsV2 = 0;end
	if isfield(sData.sConnParams,'intCortexCells'),intCortexCells = sData.sConnParams.intCortexCells;else intCortexCells = intCellsV1+intCellsV2;end
	
	
	%sConnectivity
	vecSynConductanceON_to_Cort = sData.sConnectivity.vecSynConductanceON_to_Cort;
	vecSynConductanceOFF_to_Cort = sData.sConnectivity.vecSynConductanceOFF_to_Cort; 
	vecSynWeightON_to_Cort = sData.sConnectivity.vecSynWeightON_to_Cort; 
	vecSynWeightOFF_to_Cort = sData.sConnectivity.vecSynWeightOFF_to_Cort; 
	vecSynDelayON_to_Cort = sData.sConnectivity.vecSynDelayON_to_Cort; 
	vecSynDelayOFF_to_Cort = sData.sConnectivity.vecSynDelayOFF_to_Cort; 
	matSynConnON_to_Cort = sData.sConnectivity.matSynConnON_to_Cort; 
	matSynConnOFF_to_Cort = sData.sConnectivity.matSynConnOFF_to_Cort; 
	matSynFromTo = sData.sConnectivity.matSynFromTo; 
	vecSynExcInh = sData.sConnectivity.vecSynExcInh; 
	vecSynDelay = sData.sConnectivity.vecSynDelay; 
	vecSynWeight = sData.sConnectivity.vecSynWeight; 
	vecSynConductance = sData.sConnectivity.vecSynConductance; 
	vecSynType = sData.sConnectivity.vecSynType; 
	vecCellThresh = sData.sConnectivity.vecCellThresh; 
	vecTauPeakByType = sData.sConnectivity.vecTauPeakByType; 
	vecCellTauPeak = sData.sConnectivity.vecCellTauPeak; 
	vecCellV_E = sData.sConnectivity.vecCellV_E; 
	vecCellV_I = sData.sConnectivity.vecCellV_I; 
	vecCellV_AHP = sData.sConnectivity.vecCellV_AHP; 
	vecCellV_Leak = sData.sConnectivity.vecCellV_Leak; 
	vecCellCm = sData.sConnectivity.vecCellCm; 
	vecCellG_Leak = sData.sConnectivity.vecCellG_Leak; 
	vecCellG_AHP = sData.sConnectivity.vecCellG_AHP; 
	vecCellArea = sData.sConnectivity.vecCellArea; 
	vecCellTypes = sData.sConnectivity.vecCellTypes; 
	vecPrefPsi = sData.sConnectivity.vecPrefPsi; 
	vecPrefOri = sData.sConnectivity.vecPrefOri; 
	vecPrefSF = sData.sConnectivity.vecPrefSF; 
	vecPrefRF_X = sData.sConnectivity.vecPrefRF_X; 
	vecPrefRF_Y = sData.sConnectivity.vecPrefRF_Y; 
	dblSigmaX = sData.sConnectivity.dblSigmaX; 
	dblSigmaY = sData.sConnectivity.dblSigmaY; 
	matPrefGabors = sData.sConnectivity.matPrefGabors; 
	matFieldsV2 = sData.sConnectivity.matFieldsV2;
	
	%sStimInputs
	vecStimTypeContrasts=sData.sStimInputs.vecStimTypeContrasts; 
	dblTrialDur=sData.sStimInputs.dblTrialDur; 
	dblVisSpacing=sData.sStimInputs.dblVisSpacing; 
	varDeltaSyn=sData.sStimInputs.varDeltaSyn; 
	matBlankLGN_ON=sData.sStimInputs.matBlankLGN_ON; 
	matBlankLGN_OFF=sData.sStimInputs.matBlankLGN_OFF; 
	cellR_ON=sData.sStimInputs.cellR_ON; 
	cellR_OFF=sData.sStimInputs.cellR_OFF; 
	cellLGN_ON=sData.sStimInputs.cellLGN_ON; 
	cellLGN_OFF=sData.sStimInputs.cellLGN_OFF; 
	cellContrast=sData.sStimInputs.cellContrast; 
	cellLuminance=sData.sStimInputs.cellLuminance; 
	vecTrialOris=sData.sStimInputs.vecTrialOris; 
	vecTrialOriIdx=sData.sStimInputs.vecTrialOriIdx; 
	vecTrialSFs=sData.sStimInputs.vecTrialSFs; 
	vecTrialSFIdx=sData.sStimInputs.vecTrialSFIdx; 
	vecTrialTFs=sData.sStimInputs.vecTrialTFs; 
	vecTrialTFIdx=sData.sStimInputs.vecTrialTFIdx; 
	vecTrialContrasts=sData.sStimInputs.vecTrialContrasts; 
	vecTrialContrastIdx=sData.sStimInputs.vecTrialContrastIdx; 
	matStimTypeCombos=sData.sStimInputs.matStimTypeCombos; 
	vecStimTypeOris=sData.sStimInputs.vecStimTypeOris; 
	vecStimTypeSFs=sData.sStimInputs.vecStimTypeSFs; 
	vecStimTypeTFs=sData.sStimInputs.vecStimTypeTFs;
elseif isfield(sData,'matStimTypeCombos')
	%% unpack
	%direct from run
	vecOverallT = sData.vecOverallT;
	dblDeltaT = sData.dblDeltaT;
	strStimType = sData.strStimType;
	%if isfield(sData,'cellSpikeTimesLGN_ON'),sData = rmfield(sData,{'cellSpikeTimesLGN_ON','cellSpikeTimesLGN_OFF'});end
	cellSpikeTimesCortex = sData.cellSpikeTimesCortex;
	if boolLoadLGN && isfield(sData,'cellSpikeTimesLGN_ON')
		cellSpikeTimesLGN_ON = sData.cellSpikeTimesLGN_ON;
		cellSpikeTimesLGN_OFF = sData.cellSpikeTimesLGN_OFF;
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
	vecSynConductanceON_to_Cort = sData.vecSynConductanceON_to_Cort;
	vecSynConductanceOFF_to_Cort = sData.vecSynConductanceOFF_to_Cort; 
	vecSynWeightON_to_Cort = sData.vecSynWeightON_to_Cort; 
	vecSynWeightOFF_to_Cort = sData.vecSynWeightOFF_to_Cort; 
	vecSynDelayON_to_Cort = sData.vecSynDelayON_to_Cort; 
	vecSynDelayOFF_to_Cort = sData.vecSynDelayOFF_to_Cort; 
	matSynConnON_to_Cort = sData.matSynConnON_to_Cort; 
	matSynConnOFF_to_Cort = sData.matSynConnOFF_to_Cort; 
	matSynFromTo = sData.matSynFromTo; 
	vecSynExcInh = sData.vecSynExcInh; 
	vecSynDelay = sData.vecSynDelay; 
	vecSynWeight = sData.vecSynWeight; 
	vecSynConductance = sData.vecSynConductance; 
	vecSynType = sData.vecSynType; 
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
	matPrefGabors = sData.matPrefGabors; 
	matFieldsV2 = sData.matFieldsV2;
	
	%sStimInputs
	dblTrialDur=sData.dblTrialDur; 
	dblVisSpacing=sData.dblVisSpacing; 
	varDeltaSyn=sData.varDeltaSyn; 
	matBlankLGN_ON=sData.matBlankLGN_ON; 
	matBlankLGN_OFF=sData.matBlankLGN_OFF; 
	cellR_ON=sData.cellR_ON; 
	cellR_OFF=sData.cellR_OFF; 
	cellLGN_ON=sData.cellLGN_ON; 
	cellLGN_OFF=sData.cellLGN_OFF; 
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
	
elseif isfield(sData,'matSynFromTo')
	%unpack
	intCortexCells = sData.intCortexCells;
	if isfield(sData,'intCellsV1'),intCellsV1 = sData.intCellsV1;else intCellsV1 = intCortexCells;end
	if isfield(sData,'intCellsV2'),intCellsV2 = sData.intCellsV2;else intCellsV2 = 0;end
	vecOverallT = sData.vecOverallT;
	dblDeltaT = sData.dblDeltaT;
	
	matSynFromTo = sData.matSynFromTo;
	dblSynSpikeMem = sData.dblSynSpikeMem;
	if isfield(sData,'vecSynExcInh'),vecSynExcInh = sData.vecSynExcInh;else vecSynExcInh=sData.vecSynType;end
	vecSynDelay = sData.vecSynDelay;
	vecSynConductance = sData.vecSynConductance;
	if isfield(sData,'vecSynWeight'),vecSynWeight = sData.vecSynWeight;else vecSynWeight = ones(size(vecSynConductance));end
	if isfield(sData,'vecSynType'),vecSynType = sData.vecSynType;else vecSynType=[];end
	
	vecCellThresh = sData.vecCellThresh;
	if isfield(sData,'vecTauPeakByType'),vecTauPeakByType = sData.vecTauPeakByType;else vecTauPeakByTypePeakByType=unique(sData.vecCellTauPeak);end
	vecCellV_E = sData.vecCellV_E;
	vecCellV_I = sData.vecCellV_I;
	vecCellV_AHP = sData.vecCellV_AHP;
	vecCellV_Leak = sData.vecCellV_Leak;
	vecCellCm = sData.vecCellCm;
	vecCellG_Leak = sData.vecCellG_Leak;
	vecCellG_AHP = sData.vecCellG_AHP;
	vecSynConductanceON_to_Cort = sData.vecSynConductanceON_to_Cort;
	vecSynConductanceOFF_to_Cort = sData.vecSynConductanceOFF_to_Cort;
	vecSynWeightON_to_Cort = sData.vecSynWeightON_to_Cort;
	vecSynWeightOFF_to_Cort = sData.vecSynWeightOFF_to_Cort;
	vecSynDelayON_to_Cort = sData.vecSynDelayON_to_Cort;
	vecSynDelayOFF_to_Cort = sData.vecSynDelayOFF_to_Cort;
	matSynConnON_to_Cort = sData.matSynConnON_to_Cort;
	matSynConnOFF_to_Cort = sData.matSynConnOFF_to_Cort;
	matBlankLGN_ON = sData.matBlankLGN_ON;
	matBlankLGN_OFF = sData.matBlankLGN_OFF;
	cellLGN_ON = sData.cellLGN_ON;
	cellLGN_OFF = sData.cellLGN_OFF;
	vecTrialOris = sData.vecTrialOris;
	vecTrialOriIdx = sData.vecTrialOriIdx;
	vecTrialStimType = sData.vecTrialStimType;
	vecStimStartSecs = sData.vecStimStartSecs;
	vecTrialStartSecs = sData.vecTrialStartSecs;
	vecTrialEndSecs = sData.vecTrialEndSecs;
	%vecThisV = sData.vecThisV;
	%boolStimPresent = sData.boolStimPresent;
	%intPrevTrial = sData.intPrevTrial;
	%intTrialT = sData.intTrialT;
	%intIter = sData.intIter;
	cellSpikeTimesLGN_ON = sData.cellSpikeTimesLGN_ON;
	cellSpikeTimesLGN_OFF = sData.cellSpikeTimesLGN_OFF;
	cellSpikeTimesCortex = sData.cellSpikeTimesCortex;
	%vecSpikeCounterLGN_ON = sData.vecSpikeCounterLGN_ON;
	%vecSpikeCounterLGN_OFF = sData.vecSpikeCounterLGN_OFF;
	%vecSpikeCounterCortex = sData.vecSpikeCounterCortex;
	%intPreAllocationSize = sData.intPreAllocationSize;
	
	%cell preferences
	if isfield(sData,'vecCellArea'),vecCellArea = sData.vecCellArea;
	else
		vecCellArea=ones(size(vecCellV_E));
		if isfield(sData,'intCellsV2')
			vecCellArea((end-sData.intCellsV2+1):end) = 2;
		end
	end
	vecCellTypes = sData.vecCellTypes;
	if isfield(sData,'vecPrefPsi'),vecPrefPsi = sData.vecPrefPsi;end %cell's phase offset
	if isfield(sData,'vecPrefOri'),vecPrefOri=sData.vecPrefOri;end %cell's preferred orientation
	if isfield(sData,'vecPrefSF'),vecPrefSF = sData.vecPrefSF;end
	if isfield(sData,'vecPrefRF_X'),vecPrefRF_X = sData.vecPrefRF_X;end
	if isfield(sData,'vecPrefRF_Y'),vecPrefRF_Y = sData.vecPrefRF_Y;end
	if isfield(sData,'dblSigmaX'),dblSigmaX = sData.dblSigmaX;end %length of gabor response
	if isfield(sData,'dblSigmaY'),dblSigmaY = sData.dblSigmaY;end %width of gabor response
	if isfield(sData,'matPrefGabors'),matPrefGabors = sData.matPrefGabors;end
	if isfield(sData,'matFieldsV2'),matFieldsV2 = sData.matFieldsV2;end

	dblStimDur = sData.dblStimDur;
	vecStimStopSecs = vecStimStartSecs+dblStimDur;
end
%clear sData;
%% get spiking data
%cellSpikeTimesCortex
intTrials = numel(vecTrialStimType);
intNeurons = numel(cellSpikeTimesCortex);
vecBinsTime = sort([0 vecStimStartSecs vecStimStopSecs vecStimStopSecs(end)+0.1]);
%vecBinsDur = diff(vecBinsTime);
if ~exist('matModelResp','var')
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
