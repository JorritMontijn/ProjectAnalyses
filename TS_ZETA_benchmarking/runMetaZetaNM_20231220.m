clear all;
%close all;
if isfolder('F:\Drive\MontijnHeimel_TimeseriesZeta')
	strPath = 'F:\Drive\MontijnHeimel_TimeseriesZeta';
else
	strPath = 'C:\Drive\MontijnHeimel_TimeseriesZeta';
end
strDataPath = fullfile(strPath,'\Data\');
strFigPath = fullfile(strPath,'\Figs\');

%% prep
intResamps = 500;
boolDirectQuantile = false;
strR = ['Resamp' num2str(intResamps)];
%zeta
strFile=fullpath(strDataPath, 'ZetaDataV1RunNaturalMovieResamp500.mat');
sLoad = load(strFile);
matAnovaP = sLoad.matAnovaP;
matAnovaP_optimal = sLoad.matAnovaP_optimal;
matNumSpikes = sLoad.matNumSpikes;
matTtestP = sLoad.matTtestP;
matZetaP_NoStitch = sLoad.matZetaP_NoStitch;
matZetaP_Old = sLoad.matZetaP_Old;
matZetaP_Stitch = sLoad.matZetaP_Stitch;


%z
matTtestZ = -norminv(matTtestP/2);
matZetaZ = -norminv(matZetaP_Stitch/2);
matZeta2Z = -norminv(matZetaP_NoStitch/2);
matAnovaZ = -norminv(matAnovaP_optimal/2);


%plot ROC
matAUCp = [];
matAUC_dprime = [];

%remove neurons
indRem = any(matTtestZ==0,2) | any(matZetaZ==0,2) | any(matZeta2Z==0,2) | any(matAnovaZ==0,2);
%indRem = all(matTtestZ==0,2) & all(matZetaZ==0,2) & all(matZeta2Z==0,2) & all(matAnovaZ==0,2);
matTtestZ(indRem,:) = [];
matZetaZ(indRem,:) = [];
matZeta2Z(indRem,:) = [];
matAnovaZ(indRem,:) = [];

%get p
matTtestP = 2-2*normcdf(matTtestZ);
matZetaP = 2-2*normcdf(matZetaZ);
matZeta2P = 2-2*normcdf(matZeta2Z);
matAnovaP = 2-2*normcdf(matAnovaZ);

figure
dblAlpha = 0.05;
matTtestZ(isinf(matTtestZ(:))) = max(matTtestZ(~isinf(matAnovaZ(:))));
matZetaZ(isinf(matZetaZ(:))) = max(matZetaZ(~isinf(matAnovaZ(:))));
matAnovaZ(isinf(matAnovaZ(:))) = max(matAnovaZ(~isinf(matAnovaZ(:))));
matTtestP(matTtestP(:)==0) = 1e-29;
matZetaP(matZetaP(:)==0) = 1e-29;
matAnovaP(matAnovaP(:)==0) = 1e-29;
h1 =subplot(2,3,1);
matC = [0.5 0.5 0.5;...
	0 0.8 0;...
	0.8 0 0;...
	0 0 0.8];
vecColor1 = 1 + (matZetaP(:,1) < dblAlpha & matTtestP(:,1) > dblAlpha) + 2*(matZetaP(:,1) > dblAlpha & matTtestP(:,1) < dblAlpha) + 3*(matZetaP(:,1) < dblAlpha & matTtestP(:,1) < dblAlpha);
scatter(matTtestZ(:,1),matZetaZ(:,1),100,vecColor1,'.');
colormap(h1,matC(1:max(vecColor1),:));
%xlim([0 1]);ylim([0 1]);
xlabel('Z-statistic mean-based t-test (\Phi^-^1(1-p/2))')
ylabel('ZETA (\zeta_c)')
strTit = sprintf('A) Inclusion at %s=%.3f: %s=%.3f, %s=%.3f; n=%d',getGreek('alpha'),dblAlpha,getGreek('zeta'),sum(matZetaP(:,1)<dblAlpha)/numel(matZetaP(:,1)),getGreek('mu'),sum(matTtestP(:,1)<dblAlpha)/numel(matTtestP(:,1)),size(matZetaP,1));
title(strTit)	%set(gca,'xscale','log','yscale','log');

h2=subplot(2,3,2);
vecColor2 = 1 + 1*(matZetaP(:,2) > dblAlpha & matTtestP(:,2) < dblAlpha) + 2*(matZetaP(:,2) < dblAlpha & matTtestP(:,2) > dblAlpha) + 3*(matZetaP(:,2) < dblAlpha & matTtestP(:,2) < dblAlpha);
scatter(matTtestZ(:,2),matZetaZ(:,2),100,vecColor1,'.');
colormap(h2,matC(1:max(vecColor1),:));
%xlim([0 1]);ylim([0 1]);
xlabel('Z-statistic mean-based t-test (\Phi^-^1(1-p/2))')
ylabel('ZETA (\zeta_c)')
title(sprintf('B) False alarms at %s=%.3f: %s=%.3f, %s=%.3f',getGreek('alpha'),dblAlpha,getGreek('zeta'),sum(matZetaP(:,2)<dblAlpha)/numel(matZetaP(:,2)),getGreek('mu'),sum(matTtestP(:,2)<dblAlpha)/numel(matTtestP(:,2))))
%set(gca,'xscale','log','yscale','log');

intNumN = size(matZetaP,1);
vecFP_sortedZ = sort(matZetaP(:,2));
vecFP_sortedA = sort(matAnovaP(:,2));
dblAlphaAtFpAlphaPercZ = vecFP_sortedZ(round(intNumN*dblAlpha));
dblAlphaAtFpAlphaPercA = vecFP_sortedA(round(intNumN*dblAlpha));
dblInclusionZ_at_Alpha = sum(matZetaP(:,1)<dblAlphaAtFpAlphaPercZ)/numel(matZetaP(:,1));
h4 =subplot(2,3,4);
matC = [0.5 0.5 0.5;...
	0 0.8 0;...
	0.8 0 0;...
	0 0 0.8];
vecColor1 = 1 + (matZetaP(:,1) < dblAlpha & matAnovaP(:,1) > dblAlpha) + 2*(matZetaP(:,1) > dblAlpha & matAnovaP(:,1) < dblAlpha) + 3*(matZetaP(:,1) < dblAlpha & matAnovaP(:,1) < dblAlpha);
scatter(matAnovaZ(:,1),matZetaZ(:,1),100,vecColor1,'.');
colormap(h4,matC(1:max(vecColor1),:));
%xlim([0 1]);ylim([0 1]);
xlabel('Z-statistic ANOVA (\Phi^-^1(1-p/2))')
ylabel('ZETA (\zeta_c)')
title(sprintf('C) Inclusion at FPR=%.3f: %s=%.3f, %s=%.3f; n=%d',dblAlpha,getGreek('zeta'),dblInclusionZ_at_Alpha,'A',sum(matAnovaP(:,1)<dblAlphaAtFpAlphaPercA)/numel(matAnovaP(:,1)),intNumN))
%set(gca,'xscale','log','yscale','log');

h5=subplot(2,3,5);
vecColor2 = 1 + 1*(matZetaP(:,2) > dblAlpha & matAnovaP(:,2) < dblAlpha) + 2*(matZetaP(:,2) < dblAlpha & matAnovaP(:,2) > dblAlpha) + 3*(matZetaP(:,2) < dblAlpha & matAnovaP(:,2) < dblAlpha);
scatter(matAnovaZ(:,2),matZetaZ(:,2),100,vecColor1,'.');
colormap(h5,matC(1:max(vecColor1),:));
%xlim([0 1]);ylim([0 1]);
xlabel('Z-statistic ANOVA (\Phi^-^1(1-p/2))')
ylabel('ZETA (\zeta_c)')
title(sprintf('D) False alarms at %s=%.3f: %s=%.3f, %s=%.3f',getGreek('alpha'),dblAlpha,getGreek('zeta'),sum(matZetaP(:,2)<dblAlpha)/numel(matZetaP(:,2)),'A',sum(matAnovaP(:,2)<dblAlpha)/numel(matAnovaP(:,2))))
%set(gca,'xscale','log','yscale','log');

%% plot ROC
cellColor = {lines(1),'r','k','b','m'};
%vecH(intResampNpx) = subplot(4,3,intResampNpx);
subplot(2,3,3)
maxfig;
hold on;

cellLegend = {};
hold on;
for intTest=1:4
	if intTest == 1
		matData = matZetaP;
		cellLegend(end+1) = {'ZETA'};
	elseif intTest == 2
		matData = matAnovaP;
		cellLegend(end+1) = {'ANOVA'};
	elseif intTest == 3
		matData = matTtestP;
		cellLegend(end+1) = {'T-test'};
	elseif intTest == 4
		matData = matZeta2P;
		cellLegend(end+1) = {'ZETA-NS'};
	end
	intCells = size(matData,1);
	vecBothData = cat(1,matData(:,1),matData(:,2));
	vecBothLabels = cat(1,zeros(size(matData(:,1))),ones(size(matData(:,1))));
	vecThresholds = [0; sort(vecBothData); 1];
	vecRealP = matData(:,1);
	vecShuffP = matData(:,2);
	
	vecTP = sum(vecRealP<=vecThresholds',1)/intCells;
	vecFP = sum(vecShuffP<=vecThresholds',1)/intCells;
	
	plot(vecFP,vecTP,'Color',cellColor{intTest});
	
	[dblAUC,Aci,Ase] = getAuc(vecShuffP,vecRealP,0.05,'mann-whitney');
	vecAUC(intTest) = dblAUC;
	vecAUC_se(intTest) = Ase;
	
	cellLegend{end} = [cellLegend{end} sprintf(', AUC=%.3f',dblAUC)];
end

%% run tests on aucs
AUC_T = vecAUC(3) ;
AUC_A = vecAUC(2);
AUC_Z = vecAUC(1);
AUC_Zns = vecAUC(4);
Ase_T = vecAUC_se(3) ;
Ase_A = vecAUC_se(2);
Ase_Z = vecAUC_se(1);
Ase_Zns = vecAUC_se(4);

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

%z vs zna
m0 = AUC_Zns - AUC_Z;
s0 = (Ase_Zns + Ase_Z)/2;
zZZns = abs(m0/s0);
AUC_pZZns = normcdf(zZZns,'upper')*2;


%plot
hold off;
xlabel('False positive fraction');
ylabel('Inclusion fraction');
title(sprintf('E) NM ROC; %s',strR));
legend(cellLegend,'location','best');
%legend({sprintf('ZETA-test, AUC=%.3f',vecAUC(1)),sprintf('ANOVA, AUC=%.3f',vecAUC(2)),sprintf('t-test, AUC=%.3f',vecAUC(3))},'location','best','interpreter','none');

subplot(2,3,6)
cellLegend = {};
hold on;
for intTest=1:4
	if intTest == 1
		matData = matZetaP;
		cellLegend(end+1) = {'ZETA'};
	elseif intTest == 2
		matData = matAnovaP;
		cellLegend(end+1) = {'ANOVA'};
	elseif intTest == 3
		matData = matTtestP;
		cellLegend(end+1) = {'T-test'};
	elseif intTest == 4
		matData = matZeta2P;
		cellLegend(end+1) = {'ZETA-NS'};
	end
	vecRandSorted = sort(matData(:,2));
	%vecRandSorted(vecRandSorted<1e-5)=1e-5;
	%plot(vecQuantile,vecRandSorted,'Color',cellColor{intTest});
	vecAlphas = (1:numel(vecRandSorted))/numel(vecRandSorted);
	vecFPR = sum(vecRandSorted<vecAlphas)/numel(vecRandSorted);
	plot(vecAlphas,vecFPR,'Color',cellColor{intTest});
end
xlabel(sprintf('Significance level %s',getGreek('alpha')));
ylabel(sprintf('False positive fraction'));
set(gca,'xscale','log','yscale','log');
plot([1e-3 1],[1e-3 1],'k--');
xlim([1e-3 1]);
ylim([1e-3 1]);
cellLegend(end+1) = {'Theoretical norm'};
hold off;
legend(cellLegend,'location','best');
title(sprintf('MW AUC tests; Z-Zns,p=%.1e; Z-T,p=%.1e; Z-A,p=%.1e;',...
	AUC_pZZns,AUC_pTZ,AUC_pAZ));
fixfig;maxfig;

%% save
drawnow;
export_fig(fullpath(strFigPath,['ZetaNM' strR '.tif']));
export_fig(fullpath(strFigPath,['ZetaNM' strR '.pdf']));

