clear all;
%close all;
if isfolder('F:\Drive\MontijnHeimel_TimeseriesZeta')
	strPath = 'F:\Drive\MontijnHeimel_TimeseriesZeta';
else
	strPath = 'C:\Drive\MontijnHeimel_TimeseriesZeta';
end
strDataPath = fullfile(strPath,'\Data\');
strFigPath = fullfile(strPath,'\Figs\');

%% load data
sDir=dir([strDataPath 'ZetaDataKsResamp1000.mat']);
strFile = sDir(1).name;
sLoad=load([strDataPath strFile]);

%% ROC
figure;maxfig;
%prep figure
hAx1=subplot(2,3,1);hold on; %inclusion rate at alpha=0.05
hAx2=subplot(2,3,2);hold on; %fpr rate at alpha=0.05
hAx3=subplot(2,3,3);hold on; %auc

% analysis
cellColor = {lines(1),'m'};
%vecH(intResampNpx) = subplot(4,3,intResampNpx);
cellNames = {'ZETA2','KS2','Ttest2'};
cellLegend = {};
intTrialNum = numel(sLoad.vecRunTrialNums);
vecRunTests = 1:3;
matTPR = nan(numel(vecRunTests),intTrialNum);
matTPR_se = nan(numel(vecRunTests),intTrialNum,2);
matFPR = nan(numel(vecRunTests),intTrialNum);
matFPR_se = nan(numel(vecRunTests),intTrialNum,2);
matAUC = nan(numel(vecRunTests),intTrialNum);
matAUC_se = nan(numel(vecRunTests),intTrialNum);
hold on;
dblAlpha = 0.05;
for intTest=vecRunTests
	if intTest == 1
		matFullData = sLoad.matZetaP;
	elseif intTest == 2
		matFullData = sLoad.matKsP;
	elseif intTest == 3
		matFullData = sLoad.matKsP;
	end
	for intTrialIdx = 1:intTrialNum
		intTrials = sLoad.vecRunTrialNums(intTrialIdx);
		matData = squeeze(matFullData(intTrialIdx,:,:));
		intCells = size(matData,1);
		vecBothData = cat(1,matData(:,1),matData(:,2));
		vecBothLabels = cat(1,zeros(size(matData(:,1))),ones(size(matData(:,1))));
		vecThresholds = sort(vecBothData);
		vecThresholds(isnan(vecThresholds))=1;
		vecRealP = matData(:,1);
		vecShuffP = matData(:,2);
		vecTP = sum(vecRealP<=vecThresholds',1)/sum(~isnan(vecRealP));
		vecFP = sum(vecShuffP<=vecThresholds',1)/sum(~isnan(vecShuffP));
		
		[dblTPR, dblTPR_se] = binofit(sum(vecRealP<dblAlpha),sum(~isnan(vecRealP)),0.05);
		[dblFPR, dblFPR_se] = binofit(sum(vecShuffP<dblAlpha),sum(~isnan(vecShuffP)),0.05);
		matTPR(intTest,intTrialIdx) = dblTPR;
		matTPR_se(intTest,intTrialIdx,:) = dblTPR_se;
		matFPR(intTest,intTrialIdx) = dblFPR;
		matFPR_se(intTest,intTrialIdx,:) = dblFPR_se;
		
		%plot(vecFP,vecTP,'Color',cellColor{intTest});
		
		[dblAUC,Aci,Ase] = getAuc(vecShuffP,vecRealP,0.05,'mann-whitney');
		matAUC(intTest,intTrialIdx) = dblAUC;
		matAUC_se(intTest,intTrialIdx) = Ase;
		
		%cellLegend(end+1) = {sprintf('%s, AUC=%.3f',cellNames{intTest},dblAUC)};
	end
	% plot
	plot(hAx1,sLoad.vecRunTrialNums,matTPR(intTest,:),'color',cellColor{intTest});
	plot(hAx1,sLoad.vecRunTrialNums,matTPR_se(intTest,:,1),'--','color',cellColor{intTest},'HandleVisibility','off');
	plot(hAx1,sLoad.vecRunTrialNums,matTPR_se(intTest,:,2),'--','color',cellColor{intTest},'HandleVisibility','off');
	
	plot(hAx2,sLoad.vecRunTrialNums,matFPR(intTest,:),'color',cellColor{intTest});
	plot(hAx2,sLoad.vecRunTrialNums,matFPR_se(intTest,:,1),'--','color',cellColor{intTest},'HandleVisibility','off');
	plot(hAx2,sLoad.vecRunTrialNums,matFPR_se(intTest,:,2),'--','color',cellColor{intTest},'HandleVisibility','off');
	
	plot(hAx3,sLoad.vecRunTrialNums,matAUC(intTest,:),'color',cellColor{intTest});
	plot(hAx3,sLoad.vecRunTrialNums,matAUC(intTest,:)+matAUC_se(intTest,:),'--','color',cellColor{intTest},'HandleVisibility','off');
	plot(hAx3,sLoad.vecRunTrialNums,matAUC(intTest,:)-matAUC_se(intTest,:),'--','color',cellColor{intTest},'HandleVisibility','off');
end
legend(hAx1,{cellNames{1},cellNames{2}},'location','best');
ylabel(hAx1,sprintf('TPR at %s = %.3f',getGreek('alpha'),dblAlpha));
xlabel(hAx1,'Number of trials');

plot(hAx2,sLoad.vecRunTrialNums([1 end]),[dblAlpha dblAlpha],'--','color',[0.5 0.5 0.5]);
legend(hAx2,{cellNames{1},cellNames{2},'Norm'},'location','best');
ylabel(hAx2,sprintf('FPR at %s = %.3f',getGreek('alpha'),dblAlpha));
xlabel(hAx2,'Number of trials');

plot(hAx3,sLoad.vecRunTrialNums([1 end]),[0.5 0.5],'--','color',[0.5 0.5 0.5]);
legend(hAx3,{cellNames{1},cellNames{2},'Chance'},'location','best');
ylabel(hAx3,'AUC');
xlabel(hAx3,'Number of trials');
fixfig;
