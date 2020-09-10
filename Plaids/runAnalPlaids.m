
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
intRec = 2;
sLoad = load([strDataPath cellStrSes{2,intRec} '.mat']);
sGratSes = sLoad.ses;
sGratSes = doRecalcdFoF(sGratSes,8,[],[],dblLowerValThresh,dblSecWindowSize,dblNeuropilSubtractionFactor);


matGratResp = getRespMat(sGratSes);
vecOri = sGratSes.structStim.Orientation;
matPreResp = getRespMat(sGratSes,[],-2);
vecCellTypeGrat = ismember({sGratSes.neuron.type},cellTypes);

[sOut] = getTuningCurves(matGratResp,vecOri);
vecPrefDeg = rad2deg(sOut.matFittedParams(:,1)/2); %prefDir,kappa,baseline,gain
dblMaxOriPlot = 135;
vecPrefDeg(vecPrefDeg>dblMaxOriPlot) = vecPrefDeg(vecPrefDeg>dblMaxOriPlot) - 180;
vecInclude = sOut.vecFitP < 0.05;

%% load plaid data
sLoad = load([strDataPath cellStrSes{1,intRec} '.mat']);
sPlaidSes = sLoad.ses;
sPlaidSes = doRecalcdFoF(sPlaidSes,8,[],[],dblLowerValThresh,dblSecWindowSize,dblNeuropilSubtractionFactor);

vecCellTypePlaid = ismember({sPlaidSes.neuron.type},cellTypes);
vecCellType = vecCellTypeGrat(1:(min([numel(vecCellTypeGrat) numel(vecCellTypePlaid)])));

%% pre-allocate plaid variables
intStimDur = 77;
intNeurons = numel(vecCellType);
intStimNum = numel(sPlaidSes.structStim.Orientation);
vecStopFrame = nan(1,intStimNum);
vecStartFrame = nan(1,intStimNum);
vecContrast0 = zeros(1,intStimNum);
vecContrast90 = zeros(1,intStimNum);
vecPrecedingContrast0 = zeros(1,intStimNum);
vecPrecedingContrast90 = zeros(1,intStimNum);
matPlaidResp = nan(intNeurons,intStimNum);

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
matStim(:,vecDuplicates) = [];

%% plot stimuli
%{
figure
hold on
vecOriCol = [0 90];
matColor = [1 0 0;0 0 1];
for intStim=1:numel(sPlaidSes.structStim.Orientation)
intC = ismember(vecOriCol,sPlaidSes.structStim.Orientation(intStim));
h = patch([sPlaidSes.structStim.FrameOn(intStim) sPlaidSes.structStim.FrameOff(intStim) sPlaidSes.structStim.FrameOff(intStim) sPlaidSes.structStim.FrameOn(intStim)],[0 0 1 1],matColor(intC,:),'edgecolor','none','FaceAlpha',sPlaidSes.structStim.Contrast(intStim)/100);
%alpha(h,0.5);
end
fixfig
ylim([-1 2])
%}
%% plot population tuning curves grouped by stimulus, preferred ori vs z-resp
boolRemPreceding = false;
indInclPrec = true(size(vecContrast0));
if boolRemPreceding
	indInclPrec = vecPrecedingContrast0 == 0 & vecPrecedingContrast90 == 0;
end

vecHorzC = [0 0 0 0 25 25 25 nan 50 50 nan nan 100 nan nan nan];
vecVertC = [0 25 50 100 0 25 50 nan 0 25 nan nan 0 nan nan nan];

matPlaidZ = zscore(matPlaidResp,[],2);
intPrefOriGroups = 15;
vecPrefOriEdges = linspace(dblMaxOriPlot-180,dblMaxOriPlot,intPrefOriGroups+1);
vecPrefOriCenters = vecPrefOriEdges(2:end)-mean(diff(vecPrefOriEdges))/2;
vecNeuronOriGroup = sum(bsxfun(@gt,vecPrefDeg,vecPrefOriEdges),2);

for intPlot = 2:16
	if isnan(vecHorzC(intPlot)),continue;end
	subplot(4,4,intPlot)
	
	intOri0Contrast = vecHorzC(intPlot);
	intOri90Contrast = vecVertC(intPlot);
	
	vecUseStims = vecContrast0 == intOri0Contrast & vecContrast90 == intOri90Contrast & indInclPrec;
	cellAct = cell(1,intPrefOriGroups);
	for intOriGroup = 1:intPrefOriGroups
		cellAct{intOriGroup} = mean(matPlaidZ(vecNeuronOriGroup==intOriGroup & vecInclude,vecUseStims),1);
	end
	errorbar(vecPrefOriCenters,cellfun(@mean,cellAct),cellfun(@std,cellAct)./sqrt(cellfun(@numel,cellAct)))
end