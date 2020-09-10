
%% define data
strDataPath = 'F:\Data\Processed\PlaidData\';
cellStrSes = cell(2,8);
cellStrSes(:,1) = {'20130612xyt01_ses','20130612xyt02_ses'};
cellStrSes(:,2) = {'20130625xyt01_ses','20130625xyt02_ses'};
cellStrSes(:,3) = {'20131016xyt02_ses','20131016xyt01_ses'};
cellStrSes(:,4) = {'20131022xyt02_ses','20131022xyt01_ses'};
cellStrSes(:,5) = {'20140129xyt02_ses','20140129xyt01_ses'};
cellStrSes(:,6) = {'20140314xyt09_ses','20140314xyt08_ses'};
cellStrSes(:,7) = {'20140423xyt01_ses','20140423xyt02_ses'};
cellStrSes(:,8) = {'20140425xyt10_ses','20140425xyt09_ses'};

cellGratings = cellStrSes(2,:);
cellPlaids = cellStrSes(1,:);
vecPV = [2 3 4 7 8];
cellTypes = {'neuron','PV'};

%% recalc params
dblLowerValThresh = 0.5;
dblSecWindowSize = 30;
dblNeuropilSubtractionFactor = 0.3;

%% load grating data
sAggStim = [];
for intRec=1:numel(cellGratings)
	sLoad = load([strDataPath cellStrSes{2,intRec} '.mat']);
	sGratSes = sLoad.ses;
	sGratSes = doRecalcdFoF(sGratSes,8,[],[],dblLowerValThresh,dblSecWindowSize,dblNeuropilSubtractionFactor);
	matGratResp = getRespMat(sGratSes);
	
	%check for PV
	matGratRespPV = [];
	if isfield(sGratSes,'PV')
		sGratSesPV = sGratSes;
		sGratSesPV.neuron = sGratSesPV.PV;
		sGratSesPV = doRecalcdFoF(sGratSesPV,8,[],[],dblLowerValThresh,dblSecWindowSize,dblNeuropilSubtractionFactor);
		matGratRespPV = getRespMat(sGratSesPV);
		clear sGratSesPV;
	end
	vecOri = sGratSes.structStim.Orientation;
	matPreResp = getRespMat(sGratSes,[],-2);
	vecCellTypeGrat = ismember({sGratSes.neuron.type},cellTypes);
	
	boolPlot = 0;
	[sOut] = getTuningCurves(matGratResp,vecOri,boolPlot);
	vecPrefDeg = rad2deg(sOut.matFittedParams(:,1))/2; %prefDir,kappa,baseline,gain
	dblMaxOriPlot = 135;
	vecPrefDeg(vecPrefDeg>dblMaxOriPlot) = vecPrefDeg(vecPrefDeg>dblMaxOriPlot) - 180;
	
	%% load plaid data
	sLoad = load([strDataPath cellStrSes{1,intRec} '.mat']);
	sPlaidSes = sLoad.ses;
	sPlaidSes = doRecalcdFoF(sPlaidSes,8,[],[],dblLowerValThresh,dblSecWindowSize,dblNeuropilSubtractionFactor);
	%check for PV
	intNumPVs = 0;
	if isfield(sPlaidSes,'PV')
		intNumPVs = numel(sPlaidSes.PV);
		sPlaidSes = doRecalcdFoF(sPlaidSes,8,[],'PV',dblLowerValThresh,dblSecWindowSize,dblNeuropilSubtractionFactor);
	end
	
	vecCellTypePlaid = ismember({sPlaidSes.neuron.type},cellTypes);
	vecCellType = vecCellTypeGrat(1:(min([numel(vecCellTypeGrat) numel(vecCellTypePlaid)])));
	intNeurons = numel(vecCellType);
	vecPrefDeg = vecPrefDeg(1:intNeurons);
	vecP = sOut.vecFitP(1:intNeurons);
	[dummy,dummy,vecP_bh] = fdr_bh(vecP);
	vecInclude = vecP_bh < 0.05;
	
	%% pre-allocate plaid variables
	intStimDur = 77;
	intStimNum = numel(sPlaidSes.structStim.Orientation);
	vecStopFrame = nan(1,intStimNum);
	vecStartFrame = nan(1,intStimNum);
	vecContrast0 = zeros(1,intStimNum);
	vecContrast90 = zeros(1,intStimNum);
	vecPrecedingContrast0 = zeros(1,intStimNum);
	vecPrecedingContrast90 = zeros(1,intStimNum);
	matPlaidResp = nan(intNeurons,intStimNum);
	matPlaidRespPV = nan(intNumPVs,intStimNum);
	
	%% transform stimuli
	for intStim=1:numel(sPlaidSes.structStim.Orientation)
		%% get all stimuli active at start frame
		intStartFrameT = sPlaidSes.structStim.FrameOn(intStim);
		intStopFrameT = sPlaidSes.structStim.FrameOff(intStim);
		vecAllStims = find(sPlaidSes.structStim.FrameOn <= intStartFrameT & sPlaidSes.structStim.FrameOff > intStartFrameT);
		
		%get data
		for intNeuron=1:intNeurons
			matPlaidResp(intNeuron,intStim) = mean(sPlaidSes.neuron(intNeuron).dFoF(intStartFrameT:(intStartFrameT+intStimDur)));
		end
		
		%get data
		for intPV=1:intNumPVs
			matPlaidRespPV(intNeuron,intStim) = mean(sPlaidSes.PV(intPV).dFoF(intStartFrameT:(intStartFrameT+intStimDur)));
		end
		
		%assign on/off
		vecStopFrame(intStim) = intStartFrameT;
		vecStartFrame(intStim) = intStopFrameT;
		
		%assign contrasts
		for intStimIdx=vecAllStims
			dblOri = sPlaidSes.structStim.Orientation(intStimIdx);
			%check current oris
			if dblOri==0
				vecContrast0(intStim) = sPlaidSes.structStim.Contrast(intStimIdx);
			else
				vecContrast90(intStim) = sPlaidSes.structStim.Contrast(intStimIdx);
			end
			%check prior onset
			if sPlaidSes.structStim.FrameOn(intStimIdx) < intStartFrameT
				%check current oris
				if dblOri==0
					vecPrecedingContrast0(intStim) = sPlaidSes.structStim.Contrast(intStimIdx);
				else
					vecPrecedingContrast90(intStim) = sPlaidSes.structStim.Contrast(intStimIdx);
				end
			end
		end
	end
	matStim = cat(1,vecStartFrame,vecStopFrame,vecContrast0,vecContrast90,vecPrecedingContrast0,vecPrecedingContrast90);
	
	%% remove duplicates
	vecDuplicates = find(vecStartFrame(1:(end-1))==vecStartFrame(2:end))+1;
	vecStopFrame(vecDuplicates) = [];
	vecStartFrame(vecDuplicates) = [];
	vecContrast0(vecDuplicates) = [];
	vecContrast90(vecDuplicates) = [];
	vecPrecedingContrast0(vecDuplicates) = [];
	vecPrecedingContrast90(vecDuplicates) = [];
	matPlaidResp(:,vecDuplicates) = [];
	matPlaidRespPV(:,vecDuplicates) = [];
	matStim(:,vecDuplicates) = [];
	
	%% add gratings
	matGratResp = getRespMat(sGratSes);
	vecHorzGrat = find(vecOri == 0);
	vecVertGrat = find(vecOri == 90);
	vecAsgn = (numel(vecStopFrame)+1):(numel(vecStopFrame)+numel(vecHorzGrat)+numel(vecVertGrat));
	vecStopFrame(vecAsgn) = 0;
	vecStartFrame(vecAsgn) = 0;
	vecContrast0(vecAsgn) = 100*[ones(size(vecHorzGrat)) zeros(size(vecHorzGrat))];
	vecContrast90(vecAsgn) = 100*[zeros(size(vecHorzGrat)) ones(size(vecHorzGrat))];
	vecPrecedingContrast0(vecAsgn) = 0;
	vecPrecedingContrast90(vecAsgn) = 0;
	matResp = cat(2,matPlaidResp(1:intNeurons,:),matGratResp(1:intNeurons,[vecHorzGrat vecVertGrat]));
	matRespPV = [];
	if intNumPVs > 0
		matRespPV = cat(2,matPlaidRespPV(1:intNumPVs,:),matGratRespPV(1:intNumPVs,[vecHorzGrat vecVertGrat]));
	end
	sRec(intRec).vecStopFrame = vecStopFrame;
	sRec(intRec).vecStartFrame = vecStartFrame;
	sRec(intRec).vecContrast0 = vecContrast0;
	sRec(intRec).vecContrast90 = vecContrast90;
	sRec(intRec).vecPrecedingContrast0 = vecPrecedingContrast0;
	sRec(intRec).vecPrecedingContrast90 = vecPrecedingContrast90;
	sRec(intRec).matResp = matResp;
	sRec(intRec).matRespPV = matRespPV;
	
	%% add cell vars
	sRec(intRec).vecPrefDeg = vecPrefDeg;
	sRec(intRec).vecCellType = vecCellType;
	sRec(intRec).intNeurons = intNeurons;
	sRec(intRec).vecP = vecP;
	sRec(intRec).vecP_bh = vecP_bh;
	sRec(intRec).vecInclude = vecInclude;
end

%% save
strFileOut = [strDataPath 'PreProAggregatePlaids.mat'];
fprintf('Saving prepro data to %s [%s]...\n',strFileOut,getTime);
save(strFileOut,'sRec','cellStrSes','cellGratings','cellPlaids','dblMaxOriPlot',...
	'vecPV','dblLowerValThresh','dblSecWindowSize','dblNeuropilSubtractionFactor');
fprintf('\b Done! [%s]\n',getTime);