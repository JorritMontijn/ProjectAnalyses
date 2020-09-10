%% load preprocessed data
clear all;
strDataPath = 'F:\Data\Processed\PlaidData\';
strFileOut = [strDataPath 'PreProAggregatePlaids.mat'];
load(strFileOut);
return
%sRec
%sRec(intRec).vecStopFrame = vecStopFrame;
%sRec(intRec).vecStartFrame = vecStartFrame;
%sRec(intRec).vecContrast0 = vecContrast0;
%sRec(intRec).vecContrast90 = vecContrast90;
%sRec(intRec).vecPrecedingContrast0 = vecPrecedingContrast0;
%sRec(intRec).vecPrecedingContrast90 = vecPrecedingContrast90;
%sRec(intRec).matResp = matResp(neuron,trial);
%sRec(intRec).matRespPV = matRespPV(neuron,trial);

%% plot population tuning curves grouped by stimulus, preferred ori vs z-resp
boolRemPreceding = false;
indInclPrec = true(size(vecContrast0));
if boolRemPreceding
	matAggData = getPlaidData(sRec,'vecPrecedingContrast0 == 0 & vecPrecedingContrast90 == 0');
else
	matAggData = getPlaidData(sRec,'vecPrecedingContrast0 == 0 & vecPrecedingContrast90 == 0');
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