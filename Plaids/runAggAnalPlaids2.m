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
vecR_25_25 = [];
vecR_25_50 = [];
vecR_50_25 = [];
vecR_50_50 = [];
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
	
	%calculate COSI for all plaids
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
	matCOSI_25_25 = (matPlaidR_25_25 - vecMean_0_25 - vecMean_25_0) ./ (vecMean_0_25 + vecMean_25_0);
	matCOSI_25_50 = (matPlaidR_25_50 - vecMean_0_50 - vecMean_25_0) ./ (vecMean_0_50 + vecMean_25_0);
	matCOSI_50_25 = (matPlaidR_50_25 - vecMean_0_25 - vecMean_50_0) ./ (vecMean_0_25 + vecMean_50_0);
	matCOSI_50_50 = (matPlaidR_50_50 - vecMean_0_50 - vecMean_50_0) ./ (vecMean_0_50 + vecMean_50_0);
	
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
		cellR_25_25(end+1) = {corr(matCOSI_25_25',matPV_25_25(intPV,:)')};
		cellR_25_50(end+1) = {corr(matCOSI_25_50',matPV_25_50(intPV,:)')};
		cellR_50_25(end+1) = {corr(matCOSI_50_25',matPV_50_25(intPV,:)')};
		cellR_50_50(end+1) = {corr(matCOSI_50_50',matPV_50_50(intPV,:)')};
	end
end
%calc per rec
intRecs = numel(vecPV);
for intRecIdx=1:intRecs
	intRec = vecPV(intRecIdx);
	vecR_25_25(intRecIdx) = mean(cell2vec(cellR_25_25(vecPV_source==intRec)));
	vecR_25_50(intRecIdx) = mean(cell2vec(cellR_25_50(vecPV_source==intRec)));
	vecR_50_25(intRecIdx) = mean(cell2vec(cellR_50_25(vecPV_source==intRec)));
	vecR_50_50(intRecIdx) = mean(cell2vec(cellR_50_50(vecPV_source==intRec)));
end
vecR_Asym = (vecR_25_50+vecR_50_25)/2;

matPlotDotsY = [vecR_25_25; vecR_Asym; vecR_50_50];
matPlotDotsX = [1*ones(size(vecR_25_25)); 2*ones(size(vecR_Asym)); 3*ones(size(vecR_50_50))];
close all;
figure
subplot(2,3,1)
hold on
scatter(1*ones(size(vecR_25_25)),vecR_25_25,'k')
scatter(2*ones(size(vecR_Asym)),vecR_Asym,'k')
scatter(3*ones(size(vecR_50_50)),vecR_50_50,'k')
bplot(vecR_25_25,1)
bplot(vecR_Asym,2)
bplot(vecR_50_50,3)
[h,pAsymvs25]=ttest(vecR_25_25,vecR_Asym)
[h,pAsymvs50]=ttest(vecR_50_50,vecR_Asym)
xlim([0.5 3.5])
ylabel(sprintf('r(PV,Pyr)'))
set(gca,'xtick',1:3,'xticklabel',{'Iso-25',sprintf('Asym 25/50'),'Iso-50'})
fixfig;

subplot(2,3,2)
plot(matPlotDotsX,matPlotDotsY,'k');
xlim([0.5 3.5])
ylabel(sprintf('r(PV,Pyr)'))
set(gca,'xtick',1:3,'xticklabel',{'Iso-25',sprintf('Asym 25/50'),'Iso-50'})
fixfig;

subplot(2,3,4)
hold on
scatter(1*ones(size(vecR_25_25)),vecR_25_25-vecR_Asym,'k')
scatter(2*ones(size(vecR_50_50)),vecR_50_50-vecR_Asym,'k')
xlim([0.5 2.5]);
ylim([-0.5 0.5])
ylabel(sprintf('%sr(PV,Pyr) rel. to asym',getGreek('Delta')))
set(gca,'xtick',1:3,'xticklabel',{'Iso-25','Iso-50'})
title(sprintf('t-test asym vs 25,p=%.3f; vs 50,p=%.3f,n=%d recs',pAsymvs25,pAsymvs50,intRecs))
fixfig;


subplot(2,3,5)
plot(matPlotDotsX,[vecR_Asym - vecR_Asym;
	vecR_25_25 - vecR_Asym;
	vecR_50_50 - vecR_Asym],'k');
ylabel(sprintf('%sr(PV,Pyr) rel. to asym',getGreek('Delta')))
set(gca,'xtick',1:3,'xticklabel',{sprintf('Asym 25/50'),'Iso-25','Iso-50'})
fixfig;
maxfig;

%% aggregate
vecEdges = -1:0.1:1;
vecCenters = vecEdges(2:end)-(median(diff(vecEdges))/2);
vecCorr_25_25 = cell2vec(cellR_25_25);
vecCorr_25_50 = cell2vec(cellR_25_50);
vecCorr_50_25 = cell2vec(cellR_50_25);
vecCorr_50_50 = cell2vec(cellR_50_50);

figure
subplot(2,3,1);
vecCounts25_25 = histcounts(vecCorr_25_25,vecEdges);
plot(vecCenters,vecCounts25_25);
title(sprintf('25/25, median r=%.3f',median(vecCorr_25_25)));
ylabel('Number of pyr cells')
xlabel('Avg r with PV population')
fixfig;

subplot(2,3,2);
vecCounts25_50 = histcounts(vecCorr_25_50,vecEdges);
plot(vecCenters,vecCounts25_50);
title(sprintf('25/50, median r=%.3f',median(vecCorr_25_50)));
ylabel('Number of pyr cells')
xlabel('Avg r with PV population')
fixfig;

subplot(2,3,4);
vecCounts50_25 = histcounts(vecCorr_50_25,vecEdges);
plot(vecCenters,vecCounts50_25);
title(sprintf('50/25, median r=%.3f',median(vecCorr_50_25)));
ylabel('Number of pyr cells')
xlabel('Avg r with PV population')
fixfig;

subplot(2,3,5);
vecCounts50_50 = histcounts(vecCorr_50_50,vecEdges);
plot(vecCenters,vecCounts50_50);
title(sprintf('50/50, median r=%.3f',median(vecCorr_50_50)));
ylabel('Number of pyr cells')
xlabel('Avg r with PV population')
fixfig;
normaxes;

subplot(2,3,3)
hold on
vecAsymCorr = [vecCorr_25_50 + vecCorr_50_25]/2;
bplot(vecCorr_25_25,1,'outliers')
bplot(vecAsymCorr,2,'outliers')
bplot(vecCorr_50_50,3,'outliers')
[p1]=ranksum(vecCorr_25_25,vecAsymCorr);
[p2]=ranksum(vecCorr_50_50,vecAsymCorr);

set(gca,'xtick',1:3,'xticklabel',{'Iso-25','Asym 25/50','Iso-50'})
ylabel('Correlation PV/Pyr activity')
title(sprintf('MW U-test;Asym vs 25,p=%.2e;Asym vs 50,p=%.2e',p1,p2))
fixfig;

subplot(2,3,6)
hold on
vecAsymCorr = [vecCorr_25_50 + vecCorr_50_25];
[p25_1]=ranksum(vecCorr_25_25,vecCorr_25_50);
[p25_2]=ranksum(vecCorr_25_25,vecCorr_50_25);
[p50_1]=ranksum(vecCorr_50_50,vecCorr_25_50);
[p50_2]=ranksum(vecCorr_50_50,vecCorr_50_25);

bplot(vecCorr_25_25,1,'outliers')
bplot(vecCorr_25_50,2,'outliers')
bplot(vecCorr_50_25,3,'outliers')
bplot(vecCorr_50_50,4,'outliers')
set(gca,'xtick',1:4,'xticklabel',{'Iso-25','25/50','50/25','Iso-50'})
ylabel('Correlation PV/Pyr activity')
title(sprintf('MW;25-A1,p=%.0e;25-A2,p=%.0e;50-A1,p=%.0e;50-A2,p=%.0e;',p25_1,p25_2,p50_1,p50_2))
fixfig;
maxfig;
return

%% per animal