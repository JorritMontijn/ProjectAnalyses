
%% set recording
close all;
clear all;
if isfolder('F:\Drive\MontijnHeimel_TimeseriesZeta')
	strPath = 'F:\Drive\MontijnHeimel_TimeseriesZeta';
	strDataSourcePath = 'F:\Data\Processed\Neuropixels\';
else
	strPath = 'C:\Drive\MontijnHeimel_TimeseriesZeta';
	strDataSourcePath = 'E:\DataPreProcessed\';
end
strDataTargetPath = fullfile(strPath,'\Data\');
strFigPath = fullfile(strPath,'\Figs\');
vecRunTypes = [1 2];
vecResamps = 100;%[100 200 500 1000 2000];
intResamps= numel(vecResamps);
%vecResamps = [5000 10000];
boolSave = true;

%% prep
strArea = 'PoissonPeak';
strQ = 'Q0';
strR = 'R100';
strFileSearch = ['Zeta2DataAnova' strArea 'Resamp250.mat'];
sDir = dir(fullpath(strDataTargetPath,strFileSearch));
strFile = sDir(1).name;
sLoad = load(fullpath(sDir(1).folder,strFile));

cellNeuron = sLoad.cellNeuron;
matAnovaP = sLoad.matAnova2;
matAnovaP_ub = sLoad.matAnova2_unbalanced;
matMeanP = sLoad.matTtest2;
matZetaP = sLoad.matZeta2;
matZetaP_old = sLoad.matZeta2_old;

%% plot
%flatten
matMeanZ = -norminv(matMeanP/2);
matZetaZ = -norminv(matZetaP/2);
matZetaZ_old = -norminv(matZetaP_old/2);
matAnovaZ = -norminv(matAnovaP/2);
matAnovaZ_ub = -norminv(matAnovaP_ub/2);

%plot ROC
matAUCp = [];
matAUC_dprime = [];
if size(matMeanZ,1) >= 1
	figure
	dblAlpha = 0.05;
	matMeanZ(isinf(matMeanZ(:))) = max(matMeanZ(~isinf(matAnovaZ(:))));
	matZetaZ(isinf(matZetaZ(:))) = max(matZetaZ(~isinf(matAnovaZ(:))));
	matZetaZ_old(isinf(matZetaZ_old(:))) = max(matZetaZ(~isinf(matAnovaZ(:))));
	matAnovaZ(isinf(matAnovaZ(:))) = max(matAnovaZ(~isinf(matAnovaZ(:))));
	matAnovaZ_ub(isinf(matAnovaZ_ub(:))) = max(matAnovaZ(~isinf(matAnovaZ(:))));
	
	matMeanP(matMeanP(:)==0) = 1e-29;
	matZetaP(matZetaP(:)==0) = 1e-29;
	matAnovaP(matAnovaP(:)==0) = 1e-29;
	matZetaP_old(matZetaP_old(:)==0) = 1e-29;
	matAnovaP_ub(matAnovaP_ub(:)==0) = 1e-29;
	
	h1 =subplot(2,3,1);
	matC = [0.5 0.5 0.5;...
		0 0.8 0;...
		0.8 0 0;...
		0 0 0.8];
	vecColor1 = 1 + (matZetaP(:,1) < dblAlpha & matMeanP(:,1) > dblAlpha) + 2*(matZetaP(:,1) > dblAlpha & matMeanP(:,1) < dblAlpha) + 3*(matZetaP(:,1) < dblAlpha & matMeanP(:,1) < dblAlpha);
	scatter(matMeanZ(:,1),matZetaZ(:,1),100,vecColor1,'.');
	colormap(h1,matC(1:max(vecColor1),:));
	%xlim([0 1]);ylim([0 1]);
	xlabel('Z-statistic mean-based t-test (\Phi^-^1(1-p/2))')
	ylabel('ZETA (\zeta_c)')
	strTit = sprintf('A) Inclusion at %s=%.3f: %s=%.3f, %s=%.3f; n=%d',getGreek('alpha'),dblAlpha,getGreek('zeta'),sum(matZetaP(:,1)<dblAlpha)/numel(matZetaP(:,1)),getGreek('mu'),sum(matMeanP(:,1)<dblAlpha)/numel(matMeanP(:,1)),size(matZetaP,1));
	title(strTit)	%set(gca,'xscale','log','yscale','log');
	
	h2=subplot(2,3,2);
	vecColor2 = 1 + 1*(matZetaP(:,2) > dblAlpha & matMeanP(:,2) < dblAlpha) + 2*(matZetaP(:,2) < dblAlpha & matMeanP(:,2) > dblAlpha) + 3*(matZetaP(:,2) < dblAlpha & matMeanP(:,2) < dblAlpha);
	scatter(matMeanZ(:,2),matZetaZ(:,2),100,vecColor1,'.');
	colormap(h2,matC(1:max(vecColor1),:));
	%xlim([0 1]);ylim([0 1]);
	xlabel('Z-statistic mean-based t-test (\Phi^-^1(1-p/2))')
	ylabel('ZETA (\zeta_c)')
	title(sprintf('B) False alarms at %s=%.3f: %s=%.3f, %s=%.3f',getGreek('alpha'),dblAlpha,getGreek('zeta'),sum(matZetaP(:,2)<dblAlpha)/numel(matZetaP(:,2)),getGreek('mu'),sum(matMeanP(:,2)<dblAlpha)/numel(matMeanP(:,2))))
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
	
	for intTest=1:5
		if intTest == 1
			matData = matZetaP;
		elseif intTest == 2
			matData = matAnovaP;
		elseif intTest == 3
			matData = matMeanP;
		elseif intTest == 4
			matData = matZetaP_old;
		elseif intTest == 5
			matData = matAnovaP_ub;
		end
		intCells = size(matData,1);
		vecBothData = cat(1,matData(:,1),matData(:,2));
		vecBothLabels = cat(1,zeros(size(matData(:,1))),ones(size(matData(:,1))));
		vecThresholds = sort(vecBothData);
		vecRealP = matData(:,1);
		vecShuffP = matData(:,2);
		
		vecTP = sum(vecRealP<=vecThresholds',1)/intCells;
		vecFP = sum(vecShuffP<=vecThresholds',1)/intCells;
		
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
	title(sprintf('E) ROC; %s %s %s',strQ,strR,strArea));
	legend({sprintf('ZETA-test, AUC=%.3f',vecAUC(1)),sprintf('ANOVA, AUC=%.3f',vecAUC(2)),sprintf('t-test, AUC=%.3f',vecAUC(3))},'location','best','interpreter','none');
	
	subplot(2,3,6)
	cellLegend = {};
	hold on;
	for intTest=1:5
		if intTest == 1
			matData = matZetaP;
			cellLegend(end+1) = {'ZETA'};
		elseif intTest == 2
			matData = matAnovaP;
			cellLegend(end+1) = {'ANOVA'};
		elseif intTest == 3
			matData = matMeanP;
			cellLegend(end+1) = {'T-test'};
		elseif intTest == 4
			matData = matZetaP_old;
			cellLegend(end+1) = {'ZETA-old'};
		elseif intTest == 5
			matData = matAnovaP_ub;
			cellLegend(end+1) = {'ANOVA-ub'};
		end
		vecRandSorted = sort(matData(:,2));
		vecQuantile = linspace(1/numel(vecRandSorted),1,numel(vecRandSorted));
		plot(vecQuantile,vecRandSorted,'Color',cellColor{intTest});
	end
	xlabel(sprintf('Significance level %s',getGreek('alpha')));
	ylabel(sprintf('P-value cut-off needed for FPR=%s',getGreek('alpha')));
	set(gca,'xscale','log','yscale','log');
	dblMinVal = max(get(gca,'xlim'),get(gca,'ylim'));
	plot([dblMinVal 1],[dblMinVal 1],'k--');
	cellLegend(end+1) = {'Theoretical norm'};
	hold off;
	legend(cellLegend,'location','best');
	title(sprintf('MW AUC tests; T-A,p=%.1e; T-Z,p=%.1e; A-Z,p=%.1e;',...
		AUC_pTA,AUC_pTZ,AUC_pAZ));
	fixfig;maxfig;
	
	%% save
	drawnow;
	export_fig(fullpath(strFigPath,['Zeta2' strQ strR strArea '.tif']));
	export_fig(fullpath(strFigPath,['Zeta2' strQ strR strArea '.pdf']));
end
