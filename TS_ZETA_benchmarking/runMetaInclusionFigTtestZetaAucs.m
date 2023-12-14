%% file locs
clear all;
%close all;
if isfolder('F:\Drive\MontijnHeimel_TimeseriesZeta')
	strPath = 'F:\Drive\MontijnHeimel_TimeseriesZeta';
else
	strPath = 'C:\Drive\MontijnHeimel_TimeseriesZeta';
end
strDataPath = fullfile(strPath,'\Data\');
strFigPath = fullfile(strPath,'\Figs\');

strIndicator = 'GCaMP';

%% load
cellTests = {'Spike times','Time-series','Two-sample spike t','Two-sample t-s'};
cellTests = {'Npx','GCaMP6f','2-s Npx','2-s GCaMP6f'};
figure;maxfig;
for intCase=1:4
	if intCase==1
		%zeta
		strFile=fullpath(strDataPath, 'ZetaDataAnovaPrimary visualRunDriftingGratingsResamp250.mat');
		sLoad = load(strFile);
		strFile=fullpath(strDataPath, 'ZetaDataAnovaPrimary visual-RandRunDriftingGratingsResamp250.mat');
		sLoadRand = load(strFile);
		matMeanP = cat(1,sLoad.vecTtestP,sLoadRand.vecTtestP);
		matZetaP = cat(1,sLoad.vecZetaP_UniNoStitch,sLoadRand.vecZetaP_UniNoStitch); %with replacement
		matAnovaP = cat(1,sLoad.vecAnovaP_optimal,sLoadRand.vecAnovaP_optimal);
	elseif intCase==2
		%ts-zeta
		[matMeanP,matZetaP,matAnovaP] = loadTsZeta(strIndicator);
	elseif intCase==3
		%2-sample zeta
		%strFile = fullpath(strDataPath, 'Zeta2DataAnovaV1RunDriftingGratingsResamp500.mat');
		strFile = fullpath(strDataPath, 'Zeta2DataStimDiffV1RunDriftingGratingsResamp500Q0.mat');
		%strFile = fullpath(strDataPath, 'Zeta2DataShiftRespV1RunDriftingGratingsResamp500Q0.mat');
		
		sLoad = load(strFile);
		matMeanP = sLoad.matTtest2';
		matZetaP = sLoad.matZeta2';
		matAnovaP = sLoad.matAnova2_optimal';
	elseif intCase==4
		%2-sample ts-zeta
		[matMeanP,matZetaP,matAnovaP] = loadTsZeta2('');
	end
	%remove nans
	indRem = any(isnan(matMeanP) | isnan(matZetaP) | isnan(matAnovaP),1);
	matMeanP(:,indRem)=[];
	matZetaP(:,indRem)=[];
	matAnovaP(:,indRem)=[];
	
	intN = size(matZetaP,2);
	[dblTPRZ,vecTPRZci] = binofit(sum(matZetaP(1,:)<0.05),intN);
	[dblTPRT,vecTPRTci] = binofit(sum(matMeanP(1,:)<0.05),intN);
	[dblTPRA,vecTPRAci] = binofit(sum(matAnovaP(1,:)<0.05),intN);
	
	[dblFPRZ,vecFPRZci] = binofit(sum(matZetaP(2,:)<0.05),intN);
	[dblFPRT,vecFPRTci] = binofit(sum(matMeanP(2,:)<0.05),intN);
	[dblFPRA,vecFPRAci] = binofit(sum(matAnovaP(2,:)<0.05),intN);
	
	dblSensitivityZ = dblTPRZ/dblFPRZ;
	dblSensitivityT = dblTPRT/dblFPRT;
	dblSensitivityA = dblTPRA/dblFPRA;
	
	%% plot
	% AUCs
	subplot(2,4,intCase);
	vecColZ = lines(1);
	vecColT = [0 0 0];
	vecColA = [0.8 0 0];
	cellColor = {vecColZ,vecColA,vecColT};
	hold on
	vecAUC = [];
	vecAUC_se = [];
	for intTest=1:3
		if intTest == 1
			matData = matZetaP';
			vecCol=vecColZ;
		elseif intTest == 2
			matData = matAnovaP';
			vecCol=vecColA;
		else
			matData = matMeanP';
			vecCol=vecColT;
		end
		intCells = size(matData,1);
		vecBothData = cat(1,matData(:,1),matData(:,2));
		vecBothLabels = cat(1,zeros(size(matData(:,1))),ones(size(matData(:,1))));
		vecThresholds = sort(vecBothData);
		vecThresholds(isnan(vecThresholds))=1;
		vecRealP = matData(:,1);
		vecShuffP = matData(:,2);
		%remove if both==1
		indRem = vecRealP==1 & vecShuffP==1;
		vecRealP(indRem)=[];
		vecShuffP(indRem)=[];
		
		vecTP = sum(vecRealP<=vecThresholds',1)/sum(~isnan(vecRealP));
		vecFP = sum(vecShuffP<=vecThresholds',1)/sum(~isnan(vecShuffP));
		dblAlpha = 0.05;
		[dblTPR, dblTPR_se] = binofit(sum(vecRealP<dblAlpha),sum(~isnan(vecRealP)),0.05);
		[dblFPR, dblFPR_se] = binofit(sum(vecShuffP<dblAlpha),sum(~isnan(vecShuffP)),0.05);
		
		plot(vecFP,vecTP,'Color',vecCol);
		[dblAUC,Aci,Ase] = getAuc(vecShuffP,vecRealP,0.05,'mann-whitney');
		vecAUC(intTest) = dblAUC;
		vecAUC_se(intTest) = Ase;
	end
	xlabel('Fraction of false positives');
	ylabel('Fraction of included cells');
	title(sprintf('%s:t=%.3f;a=%.3f;z=%.3f',cellTests{intCase},vecAUC(3),vecAUC(2),vecAUC(1)))
		
		
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
	
	
	%fprs
	subplot(2,4,4+intCase);
	cellLegend = {};
	hold on;
	
	for intTest=1:3
		if intTest == 1
			matData = matZetaP';
			cellLegend(end+1) = {'ZETA'};
		elseif intTest == 2
			matData = matAnovaP';
			cellLegend(end+1) = {'ANOVA'};
		elseif intTest == 3
			matData = matMeanP';
			cellLegend(end+1) = {'T-test'};
		end
		vecRandSorted = sort(matData(:,2));
		vecQuantile = linspace(1/numel(vecRandSorted),1,numel(vecRandSorted));
		%plot(vecQuantile,vecRandSorted,'Color',cellColor{intTest});
		vecAlphas = (1:numel(vecRandSorted))/numel(vecRandSorted);
		vecFPR = sum(vecRandSorted<vecAlphas)/numel(vecRandSorted);
		plot(vecAlphas,vecFPR,'Color',cellColor{intTest});
	end
	xlabel(sprintf('Significance level %s',getGreek('alpha')));
	ylabel(sprintf('Fraction of false positives'));
	%ylabel(sprintf('P-value cut-off needed for FPR=%s',getGreek('alpha')));
	set(gca,'xscale','log','yscale','log');
	dblMinVal = max(get(gca,'xlim'),get(gca,'ylim'));
	plot([dblMinVal 1],[dblMinVal 1],'k--');
	cellLegend(end+1) = {'Theoretical norm'};
	hold off;
	title(sprintf('MW; Z-T,p=%.1e; Z-A,p=%.1e;',...
		AUC_pTZ,AUC_pAZ));
	xlim([1e-2 1]);
	ylim([1e-2 1]);
end
legend(cellLegend,'location','best');
fixfig([],[],2,16)

%% save
strFile = 'MetaInclusionAucs';
export_fig(fullpath(strFigPath,[strFile '.pdf']));
export_fig(fullpath(strFigPath,[strFile '.png']));