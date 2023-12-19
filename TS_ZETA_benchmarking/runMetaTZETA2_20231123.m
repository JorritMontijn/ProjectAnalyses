clear all;
%close all;
if isfolder('F:\Drive\MontijnHeimel_TimeseriesZeta')
	strPath = 'F:\Drive\MontijnHeimel_TimeseriesZeta';
else
	strPath = 'C:\Drive\MontijnHeimel_TimeseriesZeta';
end
strDataPath = fullfile(strPath,'\Data\');
strFigPath = fullfile(strPath,'\Figs\');

intResamps = 1001; %Q1R10000T64 / Q0R250T64
strComp = 'DiffNeurons';%DiffNeurons, DiffStims, PeakHeight, PeakTime
boolDirectQuantile = false;



%% prep
strQ = ['Q' num2str(boolDirectQuantile) ];
strR = ['Resamp' num2str(intResamps)];
strTest = 'TsZeta3';
cellRunRand = {...
	'',...Rand 1
	'-Rand',...Rand 2
	};
vecRunTests = [1:4];
if contains(strComp,'Diff')
	sDirAll1=dir([strDataPath strTest '_' strComp '*' strQ '*ses' strR '.mat']);
	sDirAll2=dir([strDataPath strTest '_' strComp '*' strQ '*ses-Rand' strR '.mat']);
else
	sDirAll1=dir([strDataPath strTest '_' strComp '*' strQ '_' strR '.mat']);
	sDirAll2=dir([strDataPath strTest '_' strComp '*' strQ '_Rand' strR '.mat']);
end
sDirAll = cat(1,sDirAll1,sDirAll2);
vecRand = contains({sDirAll.name},'Rand');
sDirReal = sDirAll(~vecRand);
sDirRand = sDirAll(vecRand);
cellMeanP = [];
cellZetaP = [];
cellAnovaP = [];
cellClustP = [];

matZetaP = [];
matMeanP = [];
matClustP = [];
cellAggClustP = cellfill(cell(1,2),[1,6]);
for intRandType=1:2
	if intRandType == 1
		sDir = sDirReal;
	else
		sDir = sDirRand;
	end
	intFiles=numel(sDir);
	
	
	for intFile=1:intFiles
		strFile = sDir(intFile).name;
		sLoad=load([strDataPath strFile]);
		cellMeanP{intRandType}{intFile} = sLoad.vecTtestP;
		cellZetaP{intRandType}{intFile} = sLoad.vecTsZetaP;
		cellAnovaP{intRandType}{intFile} = sLoad.vecAnovaP;
		if isfield(sLoad,'vecClustP')
			cellClustP{intRandType}{intFile} = sLoad.vecClustP;
		elseif isfield(sLoad,'matClustP')
			for intClust=1:6
				cellAggClustP{intClust}{intRandType} = cat(1,cellAggClustP{intClust}{intRandType},sLoad.matClustP(:,intClust));
			end
		end
	end
end

%% check if recording is above chance
vecRandMeanP = cell2vec(cellMeanP{2});
vecRealMeanP = cell2vec(cellMeanP{1});
vecRandZetaP = cell2vec(cellZetaP{2});
vecRealZetaP = cell2vec(cellZetaP{1});
vecRandAnovaP = cell2vec(cellAnovaP{2});
vecRealAnovaP = cell2vec(cellAnovaP{1});
matRandClustP = [];
matRealClustP = [];
vecCutOffs = 0.05./(2.^(6:-2:-4));%0.0008    0.0031    0.0125    0.0500    0.2000    0.8000
intAlphaNum = numel(vecCutOffs);
if isempty(cellClustP)
	matClustTPRs = nan(3,intAlphaNum);
	matClustFPRs = nan(3,intAlphaNum);
	for intClust=1:intAlphaNum
		matRandClustP(intClust,:) = cellAggClustP{intClust}{2};
		matRealClustP(intClust,:) = cellAggClustP{intClust}{1};
		[dblTPR,vecTPRci] = binofit(sum(matRealClustP(intClust,:)<0.05),numel(matRealClustP(intClust,:)));
		[dblFPR,vecFPRci] = binofit(sum(matRandClustP(intClust,:)<0.05),numel(matRandClustP(intClust,:)));
		matClustTPRs(:,intClust) = cat(1,dblTPR,vecTPRci');
		matClustFPRs(:,intClust) = cat(1,dblFPR,vecFPRci');
	end
	vecRandClustP = matRandClustP(4,:)';
	vecRealClustP = matRealClustP(4,:)';
	matAggClustP = cat(3,matRealClustP,matRandClustP);
else
	vecRandClustP = cell2vec(cellClustP{2});
	vecRealClustP = cell2vec(cellClustP{1});
end
matMeanP = cat(2,vecRealMeanP,vecRandMeanP)';
matZetaP = cat(2,vecRealZetaP,vecRandZetaP)';
matAnovaP = cat(2,vecRealAnovaP,vecRandAnovaP)';
matClustP = cat(2,vecRealClustP,vecRandClustP)';

matMeanZ = cat(2,-norminv(vecRealMeanP/2),-norminv(vecRandMeanP/2))';
matZetaZ = cat(2,-norminv(vecRealZetaP/2),-norminv(vecRandZetaP/2))';
matAnovaZ = cat(2,-norminv(vecRealAnovaP/2),-norminv(vecRandAnovaP/2))';
matClustZ = cat(2,-norminv(vecRealClustP/2),-norminv(vecRandClustP/2))';

%remove nans
indRem = any(isnan(matZetaP),1) | any(isnan(matMeanP),1) | any(isnan(matAnovaZ),1);
matMeanP(:,indRem)=[];
matMeanZ(:,indRem)=[];
matZetaP(:,indRem)=[];
matZetaZ(:,indRem)=[];
matAnovaP(:,indRem)=[];
matAnovaZ(:,indRem)=[];

%rename ks => anova, wilcox=>t
%matAnovaP = matKsP;
%matAnovaZ = matKsZ;
%matMeanP = matWilcoxP;
%matMeanZ = matWilcoxZ;

%% plot
figure
matMeanZ(isinf(matMeanZ(:))) = max(matMeanZ(~isinf(matAnovaZ(:))));
matZetaZ(isinf(matZetaZ(:))) = max(matZetaZ(~isinf(matAnovaZ(:))));
matAnovaZ(isinf(matAnovaZ(:))) = max(matAnovaZ(~isinf(matAnovaZ(:))));
matMeanP(matMeanP(:)==0) = 1e-29;
matZetaP(matZetaP(:)==0) = 1e-29;
matAnovaP(matAnovaP(:)==0) = 1e-29;
h1 =subplot(2,3,1);
matC = [0.5 0.5 0.5;...
	0 0.8 0;...
	0.8 0 0;...
	0 0 0.8];
vecColor1 = 1 + (matZetaP(1,:) < 0.05 & matMeanP(1,:) > 0.05) + 2*(matZetaP(1,:) > 0.05 & matMeanP(1,:) < 0.05) + 3*(matZetaP(1,:) < 0.05 & matMeanP(1,:) < 0.05);
scatter(matMeanZ(1,:),matZetaZ(1,:),100,vecColor1,'.');
colormap(h1,matC(1:max(vecColor1),:));
%xlim([0 1]);ylim([0 1]);
xlabel('Z-statistic mean-based t-test (\Phi^-^1(1-p/2))')
ylabel('ZETA (\zeta_c)')
title(sprintf('A) Inclusion at %s=0.05: %s=%.3f, %s=%.3f; n=%d',getGreek('alpha'),getGreek('zeta'),sum(matZetaP(1,:)<0.05)/numel(matZetaP(1,:)),getGreek('mu'),sum(matMeanP(1,:)<0.05)/numel(matMeanP(1,:)),size(matZetaP,2)))
%set(gca,'xscale','log','yscale','log');

h2=subplot(2,3,2);
vecColor2 = 1 + 1*(matZetaP(2,:) > 0.05 & matMeanP(2,:) < 0.05) + 2*(matZetaP(2,:) < 0.05 & matMeanP(2,:) > 0.05) + 3*(matZetaP(2,:) < 0.05 & matMeanP(2,:) < 0.05);
scatter(matMeanZ(2,:),matZetaZ(2,:),100,vecColor1,'.');
colormap(h2,matC(1:max(vecColor1),:));
%xlim([0 1]);ylim([0 1]);
xlabel('Z-statistic mean-based t-test (\Phi^-^1(1-p/2))')
ylabel('ZETA (\zeta_c)')
title(sprintf('B) False alarms at %s=0.05: %s=%.3f, %s=%.3f',getGreek('alpha'),getGreek('zeta'),sum(matZetaP(2,:)<0.05)/numel(matZetaP(2,:)),getGreek('mu'),sum(matMeanP(2,:)<0.05)/numel(matMeanP(2,:))))
%set(gca,'xscale','log','yscale','log');
fixfig;maxfig;

dblAlpha = 0.05;
intNumN = size(matZetaP,2);
vecFP_sortedZ = sort(matZetaP(2,:));
vecFP_sortedA = sort(matAnovaP(2,:));
dblAlphaAtFpAlphaPercZ = dblAlpha;%vecFP_sortedZ(round(intNumN*dblAlpha));
dblAlphaAtFpAlphaPercA = dblAlpha;%vecFP_sortedA(round(intNumN*dblAlpha));
dblInclusionZ_at_Alpha = sum(matZetaP(1,:)<dblAlphaAtFpAlphaPercZ)/numel(matZetaP(1,:));
h4 =subplot(2,3,4);
matC = [0.5 0.5 0.5;...
	0 0.8 0;...
	0.8 0 0;...
	0 0 0.8];
vecColor1 = 1 + (matZetaP(1,:) < dblAlpha & matAnovaP(1,:) > dblAlpha) + 2*(matZetaP(1,:) > dblAlpha & matAnovaP(1,:) < dblAlpha) + 3*(matZetaP(1,:) < dblAlpha & matAnovaP(1,:) < dblAlpha);
scatter(matAnovaZ(1,:),matZetaZ(1,:),100,vecColor1,'.');
colormap(h4,matC(1:max(vecColor1),:));
%xlim([0 1]);ylim([0 1]);
xlabel('Z-statistic ANOVA (\Phi^-^1(1-p/2))')
ylabel('ZETA (\zeta_c)')
title(sprintf('C) Inclusion at FPR=%.3f: %s=%.3f, %s=%.3f; n=%d',dblAlpha,getGreek('zeta'),dblInclusionZ_at_Alpha,'A',sum(matAnovaP(1,:)<dblAlphaAtFpAlphaPercA)/numel(matAnovaP(1,:)),intNumN))
%set(gca,'xscale','log','yscale','log');
fixfig;

h5=subplot(2,3,5);
vecColor2 = 1 + 1*(matZetaP(2,:) > dblAlpha & matAnovaP(2,:) < dblAlpha) + 2*(matZetaP(2,:) < dblAlpha & matAnovaP(2,:) > dblAlpha) + 3*(matZetaP(2,:) < dblAlpha & matAnovaP(2,:) < dblAlpha);
scatter(matAnovaZ(2,:),matZetaZ(2,:),100,vecColor1,'.');
colormap(h5,matC(1:max(vecColor1),:));
%xlim([0 1]);ylim([0 1]);
xlabel('Z-statistic ANOVA (\Phi^-^1(1-p/2))')
ylabel('ZETA (\zeta_c)')
title(sprintf('D) FPR at %s=%.3f: %s=%.3f, %s=%.3f, %s=%.3f',getGreek('alpha'),dblAlpha,getGreek('zeta'),sum(matZetaP(2,:)<dblAlpha)/numel(matZetaP(2,:)),'A',sum(matAnovaP(2,:)<dblAlpha)/numel(matAnovaP(2,:)),'T',sum(matMeanP(2,:)<dblAlpha)/numel(matMeanP(2,:))))
%set(gca,'xscale','log','yscale','log');
fixfig;maxfig;

%% plot ROC
cellColor = {lines(1),'r','k','b','m'};
%vecH(intResampNpx) = subplot(4,3,intResampNpx);
subplot(2,3,3)
maxfig;
hold on;
cellNames = {'TZETA2','ANOVA','T-test','Clust','ANOVA-b'};
cellLegendAuc = {};

for intTest=vecRunTests
	if intTest == 1
		matData = matZetaP;
	elseif intTest == 2
		matData = matAnovaP;
	elseif intTest == 3
		matData = matMeanP;
	elseif intTest == 4
		matData = matClustP;
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
	
	cellLegendAuc(end+1) = {sprintf('%s, AUC=%.3f',cellNames{intTest},dblAUC)};
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
title(sprintf('E) ROC; %s %s %s',strQ,strR,strComp));
legend(cellLegendAuc,'location','best','interpreter','none');

subplot(2,3,6)
cellLegend = {};
hold on;
for intTest=vecRunTests
	if intTest == 1
		matData = matZetaP;
		cellLegend(end+1) = {'TZETA2'};
	elseif intTest == 2
		matData = matAnovaP;
		cellLegend(end+1) = {'ANOVA'};
	elseif intTest == 3
		matData = matMeanP;
		cellLegend(end+1) = {'T-test'};
	elseif intTest == 4
		matData = matClustP;
		cellLegend(end+1) = {'Clust'};
	end
	vecRandSorted = sort(matData(2,:));
	vecQuantile = linspace(1/numel(vecRandSorted),1,numel(vecRandSorted));
	
	vecFPR = sum(vecRandSorted'<vecQuantile,1)/numel(vecRandSorted);
	plot(vecQuantile,vecFPR,'Color',cellColor{intTest});
end

xlabel(sprintf('Significance level %s',getGreek('alpha')));
ylabel(sprintf('False positive rate'));
set(gca,'xscale','log','yscale','log');
dblMinVal = max(get(gca,'xlim'),get(gca,'ylim'));
plot([dblMinVal 1],[dblMinVal 1],'k--');
cellLegend(end+1) = {'Theoretical norm'};
hold off;
legend(cellLegend,'location','best');
title(sprintf('MW AUC tests; T vs A,p=%.1e; T vs Z,p=%.1e; A vs Z,p=%.1e;',...
	AUC_pTA,AUC_pTZ,AUC_pAZ));
xlim([1e-3 1]);
ylim([1e-3 1]);
fixfig;
drawnow;

%% save
export_fig(fullpath(strFigPath,[strTest strComp strQ strR '.tif']));
export_fig(fullpath(strFigPath,[strTest strComp strQ strR '.pdf']));

%% cluster plots
intN = size(matZetaP,2);

%remove nans
vecClusterAlphas = vecCutOffs;
intAlphaNum = numel(vecClusterAlphas);
matColor = redblack(intAlphaNum);
figure;maxfig;
hAx1=subplot(2,3,1);hold on;
hAx4=subplot(2,3,4);hold on;
cellLegend2 = {};
for intAlphaIdx=1:intAlphaNum
	dblAlpha = vecClusterAlphas(intAlphaIdx);
	matP = squeeze(matAggClustP(intAlphaIdx,:,:));
	cellLegend2 = cat(1,cellLegend2,{num2str(dblAlpha)});

	% AUCs
	axes(hAx1);
	vecCol=matColor(intAlphaIdx,:);
	
	intCells = size(matP,1);
	vecBothData = cat(1,matP(:,1),matP(:,2));
	vecBothLabels = cat(1,zeros(size(matP(:,1))),ones(size(matP(:,1))));
	vecThresholds = sort(vecBothData);
	vecThresholds(isnan(vecThresholds))=1;
	vecRealP = matP(:,1);
	vecShuffP = matP(:,2);
	%remove if both==1
	indRem = vecRealP==1 & vecShuffP==1;
	vecRealP(indRem)=[];
	vecShuffP(indRem)=[];
	
	
	vecRealSorted = sort(matP(:,1));
	vecQuantile = linspace(1/numel(vecRealSorted),1,numel(vecRealSorted));
	%plot(vecQuantile,vecRandSorted,'Color',cellColor{intTest});
	vecAlphas = (1:numel(vecRealSorted))/numel(vecRealSorted);
	vecTPR = sum(vecRealSorted<vecAlphas)/numel(vecRealSorted);
	plot(vecAlphas,vecTPR,'Color',vecCol);
	
	%fprs
	axes(hAx4);
	
	vecRandSorted = sort(matP(:,2));
	vecQuantile = linspace(1/numel(vecRandSorted),1,numel(vecRandSorted));
	%plot(vecQuantile,vecRandSorted,'Color',cellColor{intTest});
	vecAlphas = (1:numel(vecRandSorted))/numel(vecRandSorted);
	vecFPR = sum(vecRandSorted<vecAlphas)/numel(vecRandSorted);
	plot(vecAlphas,vecFPR,'Color',vecCol);
	
	
end

%finish
axes(hAx1);
matZ_P = matZetaP';
matC_P = matClustP';

%zeta
vecRealSorted = sort(matZ_P(:,1));
vecQuantile = linspace(1/numel(vecRealSorted),1,numel(vecRealSorted));
vecAlphas = (1:numel(vecRealSorted))/numel(vecRealSorted);
vecFPR = sum(vecRealSorted<vecAlphas)/numel(vecRealSorted);
plot(vecAlphas,vecFPR,'Color',lines(1));
%clust
vecRealSorted = sort(matC_P(:,1));
vecQuantile = linspace(1/numel(vecRealSorted),1,numel(vecRealSorted));
vecAlphas = (1:numel(vecRealSorted))/numel(vecRealSorted);
vecFPR = sum(vecRealSorted<vecAlphas)/numel(vecRealSorted);
%plot(vecAlphas,vecFPR,'Color','r');


%finish
xlabel(sprintf('Significance level %s',getGreek('alpha')));
ylabel(sprintf('Fraction of included cells'));
set(gca,'xscale','log');%,'yscale','log');
dblMinVal = max(get(gca,'xlim'),get(gca,'ylim'));
hold off;
xlim([1e-2 1]);
%ylim([1e-2 1]);
legend(cellLegend2,'location','best');

%fprs
axes(hAx4);
%zeta
vecRandSorted = sort(matZ_P(:,2));
vecQuantile = linspace(1/numel(vecRandSorted),1,numel(vecRandSorted));
vecAlphas = (1:numel(vecRandSorted))/numel(vecRandSorted);
vecFPR = sum(vecRandSorted<vecAlphas)/numel(vecRandSorted);
plot(vecAlphas,vecFPR,'Color',lines(1));
%clust
vecRandSorted = sort(matC_P(:,2));
vecQuantile = linspace(1/numel(vecRandSorted),1,numel(vecRandSorted));
vecAlphas = (1:numel(vecRandSorted))/numel(vecRandSorted);
vecFPR = sum(vecRandSorted<vecAlphas)/numel(vecRandSorted);
%plot(vecAlphas,vecFPR,'Color','r');


%finish
xlabel(sprintf('Significance level %s',getGreek('alpha')));
ylabel(sprintf('Fraction of false positives'));
set(gca,'xscale','log','yscale','log');
dblMinVal = max(get(gca,'xlim'),get(gca,'ylim'));
plot([dblMinVal 1],[dblMinVal 1],'k--');
hold off;
xlim([1e-2 1]);
ylim([1e-2 1]);

boolPatch = false;
%plot inclusions/fprs
intN = size(matZetaP,2);
subplot(2,3,2);hold on;
colormap(redblack);
[dblTPR_Z,vecTPR_Z] = binofit( sum(matZetaP(1,:)<0.05),intN);
[dblFPR_Z,vecFPR_Z] = binofit( sum(matZetaP(2,:)<0.05),intN);
cline(vecCutOffs,matClustTPRs(1,:),1:6,boolPatch);
cline(vecCutOffs,matClustTPRs(2,:),1:6,boolPatch);
cline(vecCutOffs,matClustTPRs(3,:),1:6,boolPatch);
plot(vecCutOffs([1 end]),dblTPR_Z*[1 1],'Color',lines(1));
plot(vecCutOffs([1 end]),vecTPR_Z(1)*[1 1],'--','Color',lines(1));
plot(vecCutOffs([1 end]),vecTPR_Z(2)*[1 1],'--','Color',lines(1));
set(gca,'xscale','log');
xlim([5e-4 1]);
ylim([0 0.3]);
set(gca,'xtick',[1e-3 1e-2 1e-1 1e0]);
title('Alpha=0.05');
xlabel(sprintf('Cluster threshold'));
ylabel(sprintf('Fraction of included cells'));


subplot(2,3,5);hold on;
cline(vecCutOffs,matClustFPRs(1,:),1:6,boolPatch);
cline(vecCutOffs,matClustFPRs(2,:),1:6,boolPatch);
cline(vecCutOffs,matClustFPRs(3,:),1:6,boolPatch);
plot(vecCutOffs([1 end]),dblFPR_Z*[1 1],'Color',lines(1));
plot(vecCutOffs([1 end]),vecFPR_Z(1)*[1 1],'--','Color',lines(1));
plot(vecCutOffs([1 end]),vecFPR_Z(2)*[1 1],'--','Color',lines(1));
set(gca,'xscale','log')
xlim([5e-4 1]);
ylim([0 0.3]);
set(gca,'xtick',[1e-3 1e-2 1e-1 1e0]);
title('Alpha=0.05');
xlabel(sprintf('Cluster threshold'));
ylabel(sprintf('False positive fraction'));

fixfig([],[],2,16)

%% save
export_fig(fullpath(strFigPath,[strTest strComp strQ strR 'ClustTest.tif']));
export_fig(fullpath(strFigPath,[strTest strComp strQ strR 'ClustTest.pdf']));
