clear all;
%close all;
if isfolder('F:\Drive\MontijnHeimel_TimeseriesZeta')
	strPath = 'F:\Drive\MontijnHeimel_TimeseriesZeta';
else
	strPath = 'C:\Drive\MontijnHeimel_TimeseriesZeta';
end
strDataPath = fullfile(strPath,'\Data\');
strFigPath = fullfile(strPath,'\Figs\');


dblSuperResFactor = 100;
intResamps = 1001; %Q1R10000T64 / Q0R250T64
vecRunTests = [1:3];
boolDirectQuantile = false;


strQ = ['Q' num2str(boolDirectQuantile) ];
strR = ['Resamp' num2str(intResamps)];
strTest = 'TsZeta2';


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
	sDirAll1=dir([strDataPath strTest strIndicator '*' strQ '*ses' strR '.mat']);
	sDirAll2=dir([strDataPath strTest strIndicator '*' strQ '*ses-Rand' strR '.mat']);
	sDirAll = cat(1,sDirAll1,sDirAll2);
	vecRand = contains({sDirAll.name},'Rand');
	sDirReal = sDirAll(~vecRand);
	sDirRand = sDirAll(vecRand);
	cellMeanP = [];
	cellZetaP = [];
	cellAnovaP = [];
	cellZetaP_wr = [];
	
	matZetaP = [];
	matMeanP = [];
	matZetawrP = [];
	
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
			cellZetaP_wr{intRandType}{intFile} = sLoad.vecTsZetaP_wr;
		end
	end
	if isempty(cellMeanP),continue;end
	%% check if recording is above chance
	vecResp = cellfun(@(x) sum(x<0.05),cellMeanP{1});
	vecTot = cellfun(@numel,cellMeanP{1});
	vecRespR = vecResp./vecTot;
	pBino=bonf_holm(myBinomTest(vecResp,vecTot,0.05,'two'));
	indRemRecs = pBino>0.05;
	cellMeanP{1}(indRemRecs) = [];
	cellMeanP{2}(indRemRecs) = [];
	cellZetaP{1}(indRemRecs) = [];
	cellZetaP{2}(indRemRecs) = [];
	cellAnovaP{1}(indRemRecs) = [];
	cellAnovaP{2}(indRemRecs) = [];
	cellZetaP_wr{1}(indRemRecs) = [];
	cellZetaP_wr{2}(indRemRecs) = [];
	
	vecRandMeanP = cell2vec(cellMeanP{2});
	vecRealMeanP = cell2vec(cellMeanP{1});
	vecRandZetaP = cell2vec(cellZetaP{2});
	vecRealZetaP = cell2vec(cellZetaP{1});
	vecRandAnovaP = cell2vec(cellAnovaP{2});
	vecRealAnovaP = cell2vec(cellAnovaP{1});
	vecRandZetawrP = cell2vec(cellZetaP_wr{2});
	vecRealZetawrP = cell2vec(cellZetaP_wr{1});
	
	matMeanP = cat(2,vecRealMeanP,vecRandMeanP)';
	matZetaP = cat(2,vecRealZetaP,vecRandZetaP)';
	matAnovaP = cat(2,vecRealAnovaP,vecRandAnovaP)';
	matZetawrP = cat(2,vecRealZetawrP,vecRandZetawrP)';
	
	matMeanZ = cat(2,-norminv(vecRealMeanP/2),-norminv(vecRandMeanP/2))';
	matZetaZ = cat(2,-norminv(vecRealZetaP/2),-norminv(vecRandZetaP/2))';
	matAnovaZ = cat(2,-norminv(vecRealAnovaP/2),-norminv(vecRandAnovaP/2))';
	matZetawrZ = cat(2,-norminv(vecRealZetawrP/2),-norminv(vecRandZetawrP/2))';
	
	%remove nans
	indRem = any(isnan(matZetaP),1) | any(isnan(matMeanP),1) | any(isnan(matAnovaZ),1);
	matMeanP(:,indRem)=[];
	matMeanZ(:,indRem)=[];
	matZetaP(:,indRem)=[];
	matZetaZ(:,indRem)=[];
	matAnovaP(:,indRem)=[];
	matAnovaZ(:,indRem)=[];
	matZetawrP(:,indRem)=[];
	matZetawrZ(:,indRem)=[];
	
	%rename ks => anova, wilcox=>t
	%matAnovaP = matKsP;
	%matAnovaZ = matKsZ;
	%matMeanP = matWilcoxP;
	%matMeanZ = matWilcoxZ;
	
	%% plot
	figure
	matMeanZ(isinf(matMeanZ(:))) = max(matMeanZ(~isinf(matAnovaZ(:))));
	matZetaZ(isinf(matZetaZ(:))) = max(matZetaZ(~isinf(matAnovaZ(:))));
	matZetawrZ(isinf(matZetawrZ(:))) = max(matZetawrZ(~isinf(matAnovaZ(:))));
	matAnovaZ(isinf(matAnovaZ(:))) = max(matAnovaZ(~isinf(matAnovaZ(:))));
	matMeanP(matMeanP(:)==0) = 1e-29;
	matZetaP(matZetaP(:)==0) = 1e-29;
	matZetawrP(matZetawrP(:)==0) = 1e-29;
	matAnovaP(matAnovaP(:)==0) = 1e-29;
	h1 =subplot(2,3,1);
	matC = [0.5 0.5 0.5;...
		0 0.8 0;...
		0.8 0 0;...
		0 0 0.8];
	vecColor1 = 1 + (matZetawrP(1,:) < 0.05 & matMeanP(1,:) > 0.05) + 2*(matZetawrP(1,:) > 0.05 & matMeanP(1,:) < 0.05) + 3*(matZetawrP(1,:) < 0.05 & matMeanP(1,:) < 0.05);
	scatter(matMeanZ(1,:),matZetaZ(1,:),100,vecColor1,'.');
	colormap(h1,matC(1:max(vecColor1),:));
	%xlim([0 1]);ylim([0 1]);
	xlabel('Z-statistic mean-based t-test (\Phi^-^1(1-p/2))')
	ylabel('ZETA (\zeta_c)')
	title(sprintf('A) Inclusion at %s=0.05: %s=%.3f, %s=%.3f; n=%d',getGreek('alpha'),getGreek('zeta'),sum(matZetawrP(1,:)<0.05)/numel(matZetawrP(1,:)),getGreek('mu'),sum(matMeanP(1,:)<0.05)/numel(matMeanP(1,:)),size(matZetawrP,2)))
	%set(gca,'xscale','log','yscale','log');
	
	h2=subplot(2,3,2);
	vecColor2 = 1 + 1*(matZetawrP(2,:) > 0.05 & matMeanP(2,:) < 0.05) + 2*(matZetawrP(2,:) < 0.05 & matMeanP(2,:) > 0.05) + 3*(matZetawrP(2,:) < 0.05 & matMeanP(2,:) < 0.05);
	scatter(matMeanZ(2,:),matZetaZ(2,:),100,vecColor1,'.');
	colormap(h2,matC(1:max(vecColor1),:));
	%xlim([0 1]);ylim([0 1]);
	xlabel('Z-statistic mean-based t-test (\Phi^-^1(1-p/2))')
	ylabel('ZETA (\zeta_c)')
	title(sprintf('B) False alarms at %s=0.05: %s=%.3f, %s=%.3f',getGreek('alpha'),getGreek('zeta'),sum(matZetawrP(2,:)<0.05)/numel(matZetawrP(2,:)),getGreek('mu'),sum(matMeanP(2,:)<0.05)/numel(matMeanP(2,:))))
	%set(gca,'xscale','log','yscale','log');
	fixfig;maxfig;
	
    dblAlpha = 0.05;
    intNumN = size(matZetawrP,2);
	vecFP_sortedZ = sort(matZetawrP(2,:));
	vecFP_sortedA = sort(matAnovaP(2,:));
	dblAlphaAtFpAlphaPercZ = dblAlpha;%vecFP_sortedZ(round(intNumN*dblAlpha));
	dblAlphaAtFpAlphaPercA = dblAlpha;%vecFP_sortedA(round(intNumN*dblAlpha));
	dblInclusionZ_at_Alpha = sum(matZetawrP(1,:)<dblAlphaAtFpAlphaPercZ)/numel(matZetawrP(1,:));
	h4 =subplot(2,3,4);
	matC = [0.5 0.5 0.5;...
		0 0.8 0;...
		0.8 0 0;...
		0 0 0.8];
	vecColor1 = 1 + (matZetawrP(1,:) < dblAlpha & matAnovaP(1,:) > dblAlpha) + 2*(matZetawrP(1,:) > dblAlpha & matAnovaP(1,:) < dblAlpha) + 3*(matZetawrP(1,:) < dblAlpha & matAnovaP(1,:) < dblAlpha);
	scatter(matAnovaZ(1,:),matZetaZ(1,:),100,vecColor1,'.');
	colormap(h4,matC(1:max(vecColor1),:));
	%xlim([0 1]);ylim([0 1]);
	xlabel('Z-statistic ANOVA (\Phi^-^1(1-p/2))')
	ylabel('ZETA (\zeta_c)')
	title(sprintf('C) Inclusion at FPR=%.3f: %s=%.3f, %s=%.3f; n=%d',dblAlpha,getGreek('zeta'),dblInclusionZ_at_Alpha,'A',sum(matAnovaP(1,:)<dblAlphaAtFpAlphaPercA)/numel(matAnovaP(1,:)),intNumN))
	%set(gca,'xscale','log','yscale','log');
	fixfig;
	
	h5=subplot(2,3,5);
	vecColor2 = 1 + 1*(matZetawrP(2,:) > dblAlpha & matAnovaP(2,:) < dblAlpha) + 2*(matZetawrP(2,:) < dblAlpha & matAnovaP(2,:) > dblAlpha) + 3*(matZetawrP(2,:) < dblAlpha & matAnovaP(2,:) < dblAlpha);
	scatter(matAnovaZ(2,:),matZetaZ(2,:),100,vecColor1,'.');
	colormap(h5,matC(1:max(vecColor1),:));
	%xlim([0 1]);ylim([0 1]);
	xlabel('Z-statistic ANOVA (\Phi^-^1(1-p/2))')
	ylabel('ZETA (\zeta_c)')
	title(sprintf('D) FPR at %s=%.3f: %s=%.3f, %s=%.3f, %s=%.3f',getGreek('alpha'),dblAlpha,getGreek('zeta'),sum(matZetawrP(2,:)<dblAlpha)/numel(matZetawrP(2,:)),'A',sum(matAnovaP(2,:)<dblAlpha)/numel(matAnovaP(2,:)),'T',sum(matMeanP(2,:)<dblAlpha)/numel(matMeanP(2,:))))
	%set(gca,'xscale','log','yscale','log');
	fixfig;maxfig;
	
	%% plot ROC
	cellColor = {lines(1),'r','k','b','m'};
	%vecH(intResampNpx) = subplot(4,3,intResampNpx);
	subplot(2,3,3)
    maxfig;
    hold on;
	cellNames = {'TS-ZETA-wr','ANOVA','T-test','ZETA','ANOVA-b'};
	cellLegendAuc = {};
	
    for intTest=vecRunTests
        if intTest == 1
            matData = matZetawrP;
        elseif intTest == 2
            matData = matAnovaP;
        elseif intTest == 3
            matData = matMeanP;
		 elseif intTest == 4
            matData = matZetaP;
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
	title(sprintf('E) ROC; %s %s %s',strQ,strR,strIndicator));
	legend(cellLegendAuc,'location','best','interpreter','none');

    subplot(2,3,6)
    cellLegend = {};
    hold on;
   for intTest=vecRunTests
        if intTest == 1
            matData = matZetawrP;
            cellLegend(end+1) = {'TS-ZETA-wr'};
        elseif intTest == 2
            matData = matAnovaP;
            cellLegend(end+1) = {'ANOVA'};
        elseif intTest == 3
            matData = matMeanP;
            cellLegend(end+1) = {'T-test'};
		elseif intTest == 4
			matData = matZetaP;
			cellLegend(end+1) = {'TS-ZETA'};
		end
        vecRandSorted = sort(matData(2,:));
        vecQuantile = linspace(1/numel(vecRandSorted),1,numel(vecRandSorted));
        plot(vecQuantile,vecRandSorted,'Color',cellColor{intTest});
    end
    xlabel(sprintf('Significance level %s',getGreek('alpha')));
    ylabel(sprintf('P-value threshold required to match empirical FPR'));
    set(gca,'xscale','log','yscale','log');
    xlim([0.001 1]);
    dblMinVal = max(get(gca,'xlim'),get(gca,'ylim'));
	plot([dblMinVal 1],[dblMinVal 1],'k--');
    cellLegend(end+1) = {'Theoretical norm'};
    hold off;
    legend(cellLegend,'location','best');
    title(sprintf('M-W; T vs A,p=%.1e; T vs Z,p=%.1e; A vs Z,p=%.1e;',...
        AUC_pTA,AUC_pTZ,AUC_pAZ));

	%% save
    fixfig;
	drawnow;
	export_fig(fullpath(strFigPath,[strTest strQ strR strIndicator '.tif']));
	export_fig(fullpath(strFigPath,[strTest strQ strR strIndicator '.pdf']));
	
	%% add to mega fig
	figure(hMegaFig);drawnow;
	intPlot = (boolDoOGB) + 1;
	subplot(2,3,intPlot);
	cellColor = {lines(1),'r','k','b'};
	%vecH(intResampNpx) = subplot(4,3,intResampNpx);
	hold on;
	
	for intTest=vecRunTests
		if intTest == 1
			matData = matZetawrP;
		elseif intTest == 2
			matData = matAnovaP;
		elseif intTest == 3
			matData = matMeanP;
		elseif intTest == 4
			matData = matZetaP;
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
	legend(cellLegendAuc,'location','best','interpreter','none');
	fixfig;
	
end

%% save
figure(hMegaFig)
drawnow;
export_fig(fullpath(strFigPath,[strTest 'Summary2' strQ strR '.tif']));
export_fig(fullpath(strFigPath,[strTest 'Summary2' strQ strR '.pdf']));
