clear all;
%close all;
strPath = 'D:\Data\Processed\TraceZeta\data\';
strFigPath = 'D:\Data\Results\TraceZeta\';
intResamps = 500;
boolDirectQuantile = false;
strQ = ['Q' num2str(boolDirectQuantile) ];
strR = ['Resamp' num2str(intResamps)];
cellRunRand = {...
	'',...Rand 1
	'-Rand',...Rand 2
	};


%% prep
strArea = 'Primary visual';%cellUniqueAreas{intArea}; %V1, SC, Retina, Poisson, GCaMP
strStim = 'RunDriftingGratings';

hMegaFig = figure;maxfig;

for boolDoOGB = [false true]
	%% load data
	if boolDoOGB
		strIndicator = 'OGB';
	else
		strIndicator = 'GCaMP';
	end
	sDirAll1=dir([strPath 'TsZetaTwoSample' strIndicator '*' strQ '*ses' strR '.mat']);
	sDirAll2=dir([strPath 'TsZetaTwoSample' strIndicator '*' strQ '*ses-Rand' strR '.mat']);
	sDirAll = cat(1,sDirAll1,sDirAll2);
	vecRand = contains({sDirAll.name},'Rand');
	sDirReal = sDirAll(~vecRand);
	sDirRand = sDirAll(vecRand);
	cellMeanP = [];
	cellZetaP = [];
	cellAnovaP = [];
	cellZetaDur = [];
	cellAnovaDur = [];
	
	matZetaP = [];
	matMeanP = [];

	for intRandType=1:2
		%swap data because of mistake in other script
		if intRandType == 2
			sDir = sDirReal;
		else
			sDir = sDirRand;
		end
		intFiles=numel(sDir);
		
		
		for intFile=1:intFiles
			strFile = sDir(intFile).name;
			sLoad=load([strPath strFile]);
			cellMeanP{intRandType}{intFile} = sLoad.vecTtestP;
			cellZetaP{intRandType}{intFile} = sLoad.vecTsZetaP;
			cellAnovaP{intRandType}{intFile} = sLoad.vecAnovaP;
		end
	end
	
	%% check if recording is above chance
	vecRandMeanP = cell2vec(cellMeanP{2});
	vecRealMeanP = cell2vec(cellMeanP{1});
	vecRandZetaP = cell2vec(cellZetaP{2});
	vecRealZetaP = cell2vec(cellZetaP{1});
	vecRandAnovaP = cell2vec(cellAnovaP{2});
	vecRealAnovaP = cell2vec(cellAnovaP{1});
	
	matMeanP = cat(2,vecRealMeanP,vecRandMeanP)';
	matZetaP = cat(2,vecRealZetaP,vecRandZetaP)';
	matAnovaP = cat(2,vecRealAnovaP,vecRandAnovaP)';
	
	matMeanZ = cat(2,norminv(1-vecRealMeanP/2),norminv(1-vecRandMeanP/2))';
	matZetaZ = cat(2,norminv(1-vecRealZetaP/2),norminv(1-vecRandZetaP/2))';
	matAnovaZ = cat(2,norminv(1-vecRealAnovaP/2),norminv(1-vecRandAnovaP/2))';
	
	%remove nans
	indRem = any(isnan(matZetaP),1) | any(isnan(matMeanP),1) | any(isnan(matAnovaZ),1);
	matMeanP(:,indRem)=[];
	matZetaP(:,indRem)=[];
	matMeanZ(:,indRem)=[];
	matZetaZ(:,indRem)=[];
	matAnovaP(:,indRem)=[];
	matAnovaZ(:,indRem)=[];
	
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
	fixfig;
	
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
	
	h4 =subplot(2,3,4);
	matC = [0.5 0.5 0.5;...
		0 0.8 0;...
		0.8 0 0;...
		0 0 0.8];
	vecColor1 = 1 + (matZetaP(1,:) < 0.05 & matAnovaP(1,:) > 0.05) + 2*(matZetaP(1,:) > 0.05 & matAnovaP(1,:) < 0.05) + 3*(matZetaP(1,:) < 0.05 & matAnovaP(1,:) < 0.05);
	scatter(matAnovaZ(1,:),matZetaZ(1,:),100,vecColor1,'.');
	colormap(h4,matC(1:max(vecColor1),:));
	%xlim([0 1]);ylim([0 1]);
	xlabel('Z-statistic ANOVA (\Phi^-^1(1-p/2))')
	ylabel('ZETA (\zeta_c)')
	title(sprintf('C) Inclusion at %s=0.05: %s=%.3f, %s=%.3f; n=%d',getGreek('alpha'),getGreek('zeta'),sum(matZetaP(1,:)<0.05)/numel(matZetaP(1,:)),'A',sum(matAnovaP(1,:)<0.05)/numel(matAnovaP(1,:)),size(matZetaP,2)))
	%set(gca,'xscale','log','yscale','log');
	fixfig;
	
	h5=subplot(2,3,5);
	vecColor2 = 1 + 1*(matZetaP(2,:) > 0.05 & matAnovaP(2,:) < 0.05) + 2*(matZetaP(2,:) < 0.05 & matAnovaP(2,:) > 0.05) + 3*(matZetaP(2,:) < 0.05 & matAnovaP(2,:) < 0.05);
	scatter(matAnovaZ(2,:),matZetaZ(2,:),100,vecColor1,'.');
	colormap(h5,matC(1:max(vecColor1),:));
	%xlim([0 1]);ylim([0 1]);
	xlabel('Z-statistic ANOVA (\Phi^-^1(1-p/2))')
	ylabel('ZETA (\zeta_c)')
	title(sprintf('D) False alarms at %s=0.05: %s=%.3f, %s=%.3f',getGreek('alpha'),getGreek('zeta'),sum(matZetaP(2,:)<0.05)/numel(matZetaP(2,:)),'A',sum(matAnovaP(2,:)<0.05)/numel(matAnovaP(2,:))))
	%set(gca,'xscale','log','yscale','log');
	fixfig;maxfig;
	
	%% plot ROC
	cellColor = {lines(1),'r','k'};
	%vecH(intResampNpx) = subplot(4,3,intResampNpx);
	subplot(2,3,3)
	maxfig;
	hold on;
	
	for intTest=1:3
		if intTest == 1
			matData = matZetaP;
		elseif intTest == 2
			matData = matAnovaP;
		elseif intTest == 3
			matData = matMeanP;
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
	fixfig;
	
	subplot(2,3,6)
	title(sprintf('Mann-Whitney AUC tests; T vs A,p=%.1e; T vs Z,p=%.1e; A vs Z,p=%.1e;',...
		AUC_pTA,AUC_pTZ,AUC_pAZ));
	
	%% save
	drawnow;
	export_fig(fullpath(strFigPath,['TsZeta' strQ strR strIndicator '.tif']));
	export_fig(fullpath(strFigPath,['TsZeta' strQ strR strIndicator '.pdf']));
	
	%% add to mega fig
	figure(hMegaFig);drawnow;
	intPlot = (boolDoOGB) + 1;
	subplot(2,3,intPlot);
	cellColor = {lines(1),'r','k'};
	%vecH(intResampNpx) = subplot(4,3,intResampNpx);
	hold on;
	
	for intTest=1:3
		if intTest == 1
			matData = matZetaP;
		elseif intTest == 2
			matData = matAnovaP;
		elseif intTest == 3
			matData = matMeanP;
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
	
end

%% save
figure(hMegaFig)
drawnow;
export_fig(fullpath(strFigPath,['TsZetaSummary2' strQ strR '.tif']));
export_fig(fullpath(strFigPath,['TsZetaSummary2' strQ strR '.pdf']));
