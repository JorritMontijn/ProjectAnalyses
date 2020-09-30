%% load preprocessed data
clear all;
strDataPath = 'F:\Data\Processed\PlaidData\';
strFileOut = [strDataPath 'PreProAggregatePlaids.mat'];
load(strFileOut);
vecPV = [2 3 4 7 8];%example: 7


%% plot population tuning curves grouped by stimulus, preferred ori vs z-resp
vecHorzC = [0 0 0 0 25 25 25 nan 50 50 50 nan 100 nan nan nan];
vecVertC = [0 25 50 100 0 25 50 nan 0 25 50 nan 0 nan nan nan];

intPrefOriGroups = 12;
vecPrefOriEdges = linspace(dblMaxOriPlot-180,dblMaxOriPlot,intPrefOriGroups+1);
vecPrefOriCenters = vecPrefOriEdges(2:end)-mean(diff(vecPrefOriEdges))/2;

cellAggPlotData = cell(numel(vecHorzC),intPrefOriGroups);

%%
cellR_25_25 = {};
cellR_25_50 = {};
cellR_50_25 = {};
cellR_50_50 = {};
vecPV_source = [];
for intRec=vecPV
	%sRec
	%sRec(intRec).vecStopFrame = vecStopFrame;
	%sRec(intRec).vecStartFrame = vecStartFrame;
	%sRec(intRec).vecContrast0 = vecContrast0;
	%sRec(intRec).vecContrast90 = vecContrast90;
	%sRec(intRec).vecPrecedingContrast0 = vecPrecedingContrast0;
	%sRec(intRec).vecPrecedingContrast90 = vecPrecedingContrast90;
	%sRec(intRec).matResp = matResp(neuron,trial);
	%sRec(intRec).matRespPV = matRespPV(neuron,trial);
	%sRec(intRec).vecPrefDeg = vecPrefDeg(neuron);
	
	%retrieve data
	vecInclude = sRec(intRec).vecInclude;
	vecContrast0 = sRec(intRec).vecContrast0;
	vecContrast90 = sRec(intRec).vecContrast90;
	vecPrecedingContrast0 = sRec(intRec).vecPrecedingContrast0;
	vecPrecedingContrast90 = sRec(intRec).vecPrecedingContrast90;
	vecPrefDeg = sRec(intRec).vecPrefDeg;
	matResp = sRec(intRec).matResp;
	matPreResp = sRec(intRec).matPreResp;
	matPreRemResp = matResp;% - matPreResp;
	indGrat = vecContrast0 == 100 | vecContrast90 == 100;
	indPlaid = ~indGrat;
	
	%pv data
	matRespPV = sRec(intRec).matRespPV;
	matPreRespPV = sRec(intRec).matPreRespPV;
	if isempty(matRespPV),continue;end
	matPreRemRespPV = matRespPV;% - matPreRespPV;
	
	%calc vars
	boolRemPreceding = true;
	indInclPrec = true(size(vecContrast0));
	if boolRemPreceding
		indInclPrec = vecPrecedingContrast0 == 0 & vecPrecedingContrast90 == 0;
	end
	
	%calculate COI for all plaids 
	%{
	After subtracting this baseline
	activity, a cross-orientation suppression index is calculated for every condition
	and every trial within conditions where we expect cross-orientation suppression
	to take place (the 25-25, 25-50, 50-25 and 50-50 conditions). Two reference
	conditions are used for this, which are the conditions corresponding to both of
	the gratings in a plaid, without a second grade present (as such, a 25-50 plaid
	trial analysis uses the 25-0 and 0-50 reference conditions). For the analysis we
	use the formula
	COI = (A - R1 - R2)=(R1 + R2)
	, where A is the mean activity of that neuron during the trial,
	%}
	indPlaidTrials = vecContrast0>0 & vecContrast90>0 & indInclPrec;
	indTrials_25_0 = vecContrast0==25 & vecContrast90==0 & indInclPrec;
	indTrials_50_0 = vecContrast0==50 & vecContrast90==0 & indInclPrec;
	indTrials_0_25 = vecContrast0==0 & vecContrast90==25 & indInclPrec;
	indTrials_0_50 = vecContrast0==0 & vecContrast90==50 & indInclPrec;
	
	vecMean_25_0 = mean(matPreRemResp(:,indTrials_25_0),2);
	vecMean_50_0 = mean(matPreRemResp(:,indTrials_50_0),2);
	vecMean_0_25 = mean(matPreRemResp(:,indTrials_0_25),2);
	vecMean_0_50 = mean(matPreRemResp(:,indTrials_0_50),2);
	
	%plaid resp
	matPlaidR_25_25 = matPreRemResp(:,vecContrast0==25 & vecContrast90==25 & indInclPrec);
	matPlaidR_25_50 = matPreRemResp(:,vecContrast0==25 & vecContrast90==50 & indInclPrec);
	matPlaidR_50_25 = matPreRemResp(:,vecContrast0==50 & vecContrast90==25 & indInclPrec);
	matPlaidR_50_50 = matPreRemResp(:,vecContrast0==50 & vecContrast90==50 & indInclPrec);
	
	%coi
	matCOI_25_25 = (matPlaidR_25_25 - vecMean_0_25 - vecMean_25_0) ./ (vecMean_0_25 + vecMean_25_0);
	matCOI_25_50 = (matPlaidR_25_50 - vecMean_0_50 - vecMean_25_0) ./ (vecMean_0_50 + vecMean_25_0);
	matCOI_50_25 = (matPlaidR_50_25 - vecMean_0_25 - vecMean_50_0) ./ (vecMean_0_25 + vecMean_50_0);
	matCOI_50_50 = (matPlaidR_50_50 - vecMean_0_50 - vecMean_50_0) ./ (vecMean_0_50 + vecMean_50_0);
	
	%pv act
	matPV_25_25 = matPreRemRespPV(:,vecContrast0==25 & vecContrast90==25 & indInclPrec);
	matPV_25_50 = matPreRemRespPV(:,vecContrast0==25 & vecContrast90==50 & indInclPrec);
	matPV_50_25 = matPreRemRespPV(:,vecContrast0==50 & vecContrast90==25 & indInclPrec);
	matPV_50_50 = matPreRemRespPV(:,vecContrast0==50 & vecContrast90==50 & indInclPrec);
	
	%mean
	matPV_25_25 = mean(matPV_25_25,1);
	matPV_25_50 = mean(matPV_25_50,1);
	matPV_50_25 = mean(matPV_50_25,1);
	matPV_50_50 = mean(matPV_50_50,1);
	
	%correlation
	intNumPV = size(matPV_25_25,1);
	vecPV_source((end+1):(end+intNumPV)) = intRec;
	for intPV=1:intNumPV
		cellR_25_25(end+1) = {corr(matCOI_25_25',matPV_25_25(intPV,:)')};
		cellR_25_50(end+1) = {corr(matCOI_25_50',matPV_25_50(intPV,:)')};
		cellR_50_25(end+1) = {corr(matCOI_50_25',matPV_50_25(intPV,:)')};
		cellR_50_50(end+1) = {corr(matCOI_50_50',matPV_50_50(intPV,:)')};
	end
end

%% aggregate
vecEdges = -1:0.1:1;
vecCenters = vecEdges(2:end)-(median(diff(vecEdges))/2);
vecCorr_25_25 = cell2vec(cellR_25_25);
vecCorr_25_50 = cell2vec(cellR_25_50);
vecCorr_50_25 = cell2vec(cellR_50_25);
vecCorr_50_50 = cell2vec(cellR_50_50);

[p1]=ranksum(vecCorr_25_25,vecCorr_25_50);
[p2]=ranksum(vecCorr_25_25,vecCorr_50_25);
[p3]=ranksum(vecCorr_25_50,vecCorr_50_50);
[p4]=ranksum(vecCorr_25_50,vecCorr_50_25);
[p5]=ranksum(vecCorr_50_25,vecCorr_50_50);
[p6]=ranksum(vecCorr_25_25,vecCorr_50_50);

figure
subplot(2,3,1);
vecCounts25_25 = histcounts(vecCorr_25_25,vecEdges);
plot(vecCenters,vecCounts25_25);
title(sprintf('25/25, median r=%.3f',median(vecCorr_25_25)));

subplot(2,3,2);
vecCounts25_50 = histcounts(vecCorr_25_50,vecEdges);
plot(vecCenters,vecCounts25_50);
title(sprintf('25/50, median r=%.3f',median(vecCorr_25_50)));

subplot(2,3,4);
vecCounts50_25 = histcounts(vecCorr_50_25,vecEdges);
plot(vecCenters,vecCounts50_25);
title(sprintf('50/25, median r=%.3f',median(vecCorr_50_25)));

subplot(2,3,5);
vecCounts50_50 = histcounts(vecCorr_50_50,vecEdges);
plot(vecCenters,vecCounts50_50);
title(sprintf('50/50, median r=%.3f',median(vecCorr_50_50)));

normaxes;
return
figure
subplot(2,3,1)
plot(vecCenters,vecCounts25_25 + vecCounts50_50);

subplot(2,3,2)
plot(vecCenters,vecCounts25_50 + vecCounts50_25);
normaxes;

