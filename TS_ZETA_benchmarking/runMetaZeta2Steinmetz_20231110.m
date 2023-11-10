clear all;
%close all;
if isfolder('F:\Drive\MontijnHeimel_TimeseriesZeta')
	strPath = 'F:\Drive\MontijnHeimel_TimeseriesZeta';
else
	strPath = 'C:\Drive\MontijnHeimel_TimeseriesZeta';
end
strDataPath = fullfile(strPath,'\Data\');
strFigPath = fullfile(strPath,'\Figs\');
boolGrouping = true;
if boolGrouping
	strGrouping = 'Lumped';
else
	strGrouping = 'Split';
end

cellRunRand = {...
	'',...Rand 1
	'-Rand',...Rand 2
	};


%% prep
strArea = 'Primary visual';%cellUniqueAreas{intArea}; %V1, SC, Retina, Poisson, GCaMP
strStim = 'RunDriftingGratings';


strFile = ['Zeta2Steinmetz' strGrouping '*.mat'];
sFiles = dir(fullpath(strDataPath,strFile));
cellFiles = {sFiles.name};
intFiles=numel(cellFiles);

matTtest2P = [];
matZeta2P = [];
matAnova2P = [];
vecRecording = [];
cellArea = {};
matTrialNum = [];
for intFile=1:intFiles
	strFile = sFiles(intFile).name;
	sLoad=load([strDataPath strFile]);
	intN = size(sLoad.matTtest2P,1);
	vecIdx = size(matTtest2P,1) + (1:intN);
	matTtest2P(vecIdx,1) = nanmean(sLoad.matTtest2P(:,:,1),2);
	matTtest2P(vecIdx,2) = nanmean(sLoad.matTtest2P(:,:,2),2);
	matZeta2P(vecIdx,1) = nanmean(sLoad.matZeta2P(:,:,1),2);
	matZeta2P(vecIdx,2) = nanmean(sLoad.matZeta2P(:,:,2),2);
	matAnova2P(vecIdx,1) = nanmean(sLoad.matAnova2P(:,:,1),2);
	matAnova2P(vecIdx,2) = nanmean(sLoad.matAnova2P(:,:,2),2);
	%cellArea((end+1):(end+intN)) = sLoad.cellArea;
	%matTrialNum((end+1):(end+intN),1) = sLoad.matTrialNum(:,1);
	%matTrialNum((end+1):(end+intN),2) = sLoad.matTrialNum(:,2);
	vecRecording(vecIdx,1) = intFile;
end


%hMegaFig = figure;maxfig;

%% check if recording is above chance
cellRealT = groupby(matTtest2P(:,1),vecRecording);
vecResp = cellfun(@(x) sum(x<0.05),cellRealT);
vecTot = cellfun(@(x) sum(~isnan(x)),cellRealT);
vecRespR = vecResp./vecTot;
pBino=bonf_holm(myBinomTest(vecResp,vecTot,0.05,'two'));
vecRemRecs = find(pBino>0.05);
if ~isempty(vecRemRecs)
	indRemEntries = ismember(vecRecording,vecRemRecs);
	matTtest2P(indRemEntries,:) = [];
	matZeta2P(indRemEntries,:) = [];
	matAnova2P(indRemEntries,:) = [];
	cellArea(indRemEntries) = [];
	matTrialNum(indRemEntries,:) = [];
	vecRecording(indRemEntries) = [];
end

%remove nans
indRem = any(isnan(matZeta2P),2) | any(isnan(matTtest2P),2) | any(isnan(matAnova2P),2);
matTtest2P(indRem,:)=[];
matZeta2P(indRem,:)=[];
matAnova2P(indRem,:)=[];
cellArea(indRem)=[];
matTrialNum(indRem,:)=[];
vecRecording(indRem)=[];

%% plot
matTtest2P = matTtest2P';
matZeta2P = matZeta2P';
matAnova2P = matAnova2P';

matMeanZ = -norminv(matTtest2P'/2);
matZetaZ = -norminv(matZeta2P'/2);
matAnovaZ = -norminv(matAnova2P'/2);

figure
matMeanZ(isinf(matMeanZ(:))) = max(matMeanZ(~isinf(matAnovaZ(:))));
matZetaZ(isinf(matZetaZ(:))) = max(matZetaZ(~isinf(matAnovaZ(:))));
matAnovaZ(isinf(matAnovaZ(:))) = max(matAnovaZ(~isinf(matAnovaZ(:))));

matTtest2P(matTtest2P(:)==0) = 1e-29;
matZeta2P(matZeta2P(:)==0) = 1e-29;
matAnova2P(matAnova2P(:)==0) = 1e-29;
h1 =subplot(2,3,1);
matC = [0.5 0.5 0.5;...
	0 0.8 0;...
	0.8 0 0;...
	0 0 0.8];
vecColor1 = 1 + (matZeta2P(1,:) < 0.05 & matTtest2P(1,:) > 0.05) + 2*(matZeta2P(1,:) > 0.05 & matTtest2P(1,:) < 0.05) + 3*(matZeta2P(1,:) < 0.05 & matTtest2P(1,:) < 0.05);
scatter(matMeanZ(1,:),matZetaZ(1,:),100,vecColor1,'.');
colormap(h1,matC(1:max(vecColor1),:));
%xlim([0 1]);ylim([0 1]);
xlabel('Z-statistic mean-based t-test (\Phi^-^1(1-p/2))')
ylabel('ZETA (\zeta_c)')
title(sprintf('A) Inclusion at %s=0.05: %s=%.3f, %s=%.3f; n=%d',getGreek('alpha'),getGreek('zeta'),sum(matZeta2P(1,:)<0.05)/numel(matZeta2P(1,:)),getGreek('mu'),sum(matTtest2P(1,:)<0.05)/numel(matTtest2P(1,:)),size(matZeta2P,2)))
%set(gca,'xscale','log','yscale','log');

h2=subplot(2,3,2);
vecColor2 = 1 + 1*(matZeta2P(2,:) > 0.05 & matTtest2P(2,:) < 0.05) + 2*(matZeta2P(2,:) < 0.05 & matTtest2P(2,:) > 0.05) + 3*(matZeta2P(2,:) < 0.05 & matTtest2P(2,:) < 0.05);
scatter(matMeanZ(2,:),matZetaZ(2,:),100,vecColor1,'.');
colormap(h2,matC(1:max(vecColor1),:));
%xlim([0 1]);ylim([0 1]);
xlabel('Z-statistic mean-based t-test (\Phi^-^1(1-p/2))')
ylabel('ZETA (\zeta_c)')
title(sprintf('B) False alarms at %s=0.05: %s=%.3f, %s=%.3f',getGreek('alpha'),getGreek('zeta'),sum(matZeta2P(2,:)<0.05)/numel(matZeta2P(2,:)),getGreek('mu'),sum(matTtest2P(2,:)<0.05)/numel(matTtest2P(2,:))))
%set(gca,'xscale','log','yscale','log');
fixfig;maxfig;

dblAlpha = 0.05;
intNumN = size(matZeta2P,2);
vecFP_sortedZ = sort(matZeta2P(2,:));
vecFP_sortedA = sort(matAnova2P(2,:));
dblAlphaAtFpAlphaPercZ = dblAlpha;%vecFP_sortedZ(round(intNumN*dblAlpha));
dblAlphaAtFpAlphaPercA = dblAlpha;%vecFP_sortedA(round(intNumN*dblAlpha));
dblInclusionZ_at_Alpha = sum(matZeta2P(1,:)<dblAlphaAtFpAlphaPercZ)/numel(matZeta2P(1,:));
h4 =subplot(2,3,4);
matC = [0.5 0.5 0.5;...
	0 0.8 0;...
	0.8 0 0;...
	0 0 0.8];
vecColor1 = 1 + (matZeta2P(1,:) < dblAlpha & matAnova2P(1,:) > dblAlpha) + 2*(matZeta2P(1,:) > dblAlpha & matAnova2P(1,:) < dblAlpha) + 3*(matZeta2P(1,:) < dblAlpha & matAnova2P(1,:) < dblAlpha);
scatter(matAnovaZ(1,:),matZetaZ(1,:),100,vecColor1,'.');
colormap(h4,matC(1:max(vecColor1),:));
%xlim([0 1]);ylim([0 1]);
xlabel('Z-statistic ANOVA (\Phi^-^1(1-p/2))')
ylabel('ZETA (\zeta_c)')
title(sprintf('C) Inclusion at FPR=%.3f: %s=%.3f, %s=%.3f; n=%d',dblAlpha,getGreek('zeta'),dblInclusionZ_at_Alpha,'A',sum(matAnova2P(1,:)<dblAlphaAtFpAlphaPercA)/numel(matAnova2P(1,:)),intNumN))
%set(gca,'xscale','log','yscale','log');
fixfig;

h5=subplot(2,3,5);
vecColor2 = 1 + 1*(matZeta2P(2,:) > dblAlpha & matAnova2P(2,:) < dblAlpha) + 2*(matZeta2P(2,:) < dblAlpha & matAnova2P(2,:) > dblAlpha) + 3*(matZeta2P(2,:) < dblAlpha & matAnova2P(2,:) < dblAlpha);
scatter(matAnovaZ(2,:),matZetaZ(2,:),100,vecColor1,'.');
colormap(h5,matC(1:max(vecColor1),:));
%xlim([0 1]);ylim([0 1]);
xlabel('Z-statistic ANOVA (\Phi^-^1(1-p/2))')
ylabel('ZETA (\zeta_c)')
title(sprintf('D) False alarms at %s=%.3f: %s=%.3f, %s=%.3f',getGreek('alpha'),dblAlpha,getGreek('zeta'),sum(matZeta2P(2,:)<dblAlpha)/numel(matZeta2P(2,:)),'A',sum(matAnova2P(2,:)<dblAlpha)/numel(matAnova2P(2,:))))
%set(gca,'xscale','log','yscale','log');
fixfig;maxfig;

%% plot ROC
cellColor = {lines(1),'r','k','b','m'};
%vecH(intResampNpx) = subplot(4,3,intResampNpx);
subplot(2,3,3)
maxfig;
hold on;

for intTest=1:3
	if intTest == 1
		matData = matZeta2P;
	elseif intTest == 2
		matData = matAnova2P;
	elseif intTest == 3
		matData = matTtest2P;
	elseif intTest == 4
		matData = matKsP;
	elseif intTest == 5
		matData = matWilcoxP;
	end
	intCells = size(matData,2);
	vecBothData = cat(2,matData(1,:),matData(2,:));
	vecBothLabels = cat(2,zeros(size(matData(1,:))),ones(size(matData(1,:))));
	vecThresholds = sort(vecBothData);
	vecRealP = matData(1,:);
	vecShuffP = matData(2,:);
	
	vecTP = sum(vecRealP<=vecThresholds',2)/intCells;
	vecFP = sum(vecShuffP<=vecThresholds',2)/intCells;
	
	plot(vecFP,vecTP,'Color',cellColor{intTest});
	
	[dblAUC,Aci,Ase] = getAuc(vecShuffP,vecRealP,0.05,'mann-whitney');
	vecAUC(intTest) = dblAUC;
	vecAUC_se(intTest) = Ase;
end

%% run tests on aucs
AUC_T = vecAUC(3) ;
AUC_A = vecAUC(2);
AUC_Z = vecAUC(1);
Ase_T = vecAUC_se(3) ;
Ase_A = vecAUC_se(2);
Ase_Z = vecAUC_se(1);

%t vs a
m0 = AUC_T - AUC_A;
s0 = (Ase_T + Ase_A)/2;
zTA = abs(m0/s0);
AUC_pTA = normcdf(zTA,'upper')*2;

%t v z
m0 = AUC_T - AUC_Z;
s0 = (Ase_T + Ase_Z)/2;
zTZ = abs(m0/s0);
AUC_pTZ = normcdf(zTZ,'upper')*2;

%a vs z
m0 = AUC_A - AUC_Z;
s0 = (Ase_A + Ase_Z)/2;
zAZ = abs(m0/s0);
AUC_pAZ = normcdf(zAZ,'upper')*2;


%plot
hold off;
xlabel('False positive fraction');
ylabel('Inclusion fraction');
title(sprintf('E) ROC; %s %s %s',strQ,strR,strIndicator));
legend({sprintf('ZETA-test, AUC=%.3f',vecAUC(1)),sprintf('ANOVA, AUC=%.3f',vecAUC(2)),sprintf('t-test, AUC=%.3f',vecAUC(3))},'location','best','interpreter','none');

subplot(2,3,6)
cellLegend = {};
hold on;
for intTest=1:3
	if intTest == 1
		matData = matZeta2P;
		cellLegend(end+1) = {'T-ZETA'};
	elseif intTest == 2
		matData = matAnova2P;
		cellLegend(end+1) = {'ANOVA'};
	elseif intTest == 3
		matData = matTtest2P;
		cellLegend(end+1) = {'T-test'};
	elseif intTest == 4
		matData = matKsP;
		cellLegend(end+1) = {'K-S test'};
	elseif intTest == 5
		matData = matWilcoxP;
		cellLegend(end+1) = {'Wilcoxon'};
	end
	vecRandSorted = sort(matData(2,:));
	vecQuantile = linspace(1/numel(vecRandSorted),1,numel(vecRandSorted));
	plot(vecQuantile,vecRandSorted,'Color',cellColor{intTest});
end
xlabel(sprintf('Significance level %s',getGreek('alpha')));
ylabel(sprintf('P-value threshold required to match empirical FPR'));
set(gca,'xscale','log','yscale','log');
dblMinVal = max(get(gca,'xlim'),get(gca,'ylim'));
plot([dblMinVal 1],[dblMinVal 1],'k--');
cellLegend(end+1) = {'Theoretical norm'};
hold off;
legend(cellLegend,'location','best');
title(sprintf('MW AUC tests; T vs A,p=%.1e; T vs Z,p=%.1e; A vs Z,p=%.1e;',...
	AUC_pTA,AUC_pTZ,AUC_pAZ));

%% save
fixfig;
drawnow;
%export_fig(fullpath(strFigPath,['Zeta2' strQ strR strIndicator '.tif']));
%export_fig(fullpath(strFigPath,['Zeta2' strQ strR strIndicator '.pdf']));

%% add to mega fig
figure(hMegaFig);drawnow;
intPlot = (boolDoOGB) + 1;
subplot(2,3,intPlot);
cellColor = {lines(1),'r','k'};
%vecH(intResampNpx) = subplot(4,3,intResampNpx);
hold on;

for intTest=1:3
	if intTest == 1
		matData = matZeta2P;
	elseif intTest == 2
		matData = matAnova2P;
	elseif intTest == 3
		matData = matTtest2P;
	end
	intCells = size(matData,2);
	vecBothData = cat(2,matData(1,:),matData(2,:));
	vecBothLabels = cat(2,zeros(size(matData(1,:))),ones(size(matData(1,:))));
	vecThresholds = sort(vecBothData);
	vecRealP = matData(1,:);
	vecShuffP = matData(2,:);
	
	vecTP = sum(vecRealP<=vecThresholds',2)/intCells;
	vecFP = sum(vecShuffP<=vecThresholds',2)/intCells;
	
	plot(vecFP,vecTP,'Color',cellColor{intTest});
	
	[dblAUC,Aci] = getAuc(vecShuffP,vecRealP);
	vecAUC(intTest) = dblAUC;
end
hold off;
xlabel('False positive fraction');
ylabel('Inclusion fraction');
title(sprintf('ROC; %s %s %s',strQ,strR,strIndicator));
legend({sprintf('ZETA-test, AUC=%.3f',vecAUC(1)),sprintf('ANOVA, AUC=%.3f',vecAUC(2)),sprintf('t-test, AUC=%.3f',vecAUC(3))},'location','best','interpreter','none');
fixfig;
