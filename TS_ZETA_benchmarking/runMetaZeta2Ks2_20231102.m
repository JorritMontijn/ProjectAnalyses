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
sDir=dir([strDataPath 'Zeta2DataKsResamp1000.mat']);
strFile = sDir(1).name;
sLoad=load([strDataPath strFile]);

%% ROC
figure;maxfig;
%prep figure
hAx1=subplot(2,3,1);hold on; %inclusion rate at alpha=0.05
hAx2=subplot(2,3,2);hold on; %fpr rate at alpha=0.05
hAx3=subplot(2,3,3);hold on; %auc

% analysis
cellColor = {lines(1),'m','k'};
%vecH(intResampNpx) = subplot(4,3,intResampNpx);
cellNames = {'ZETA2','KS2','Ttest2'};
cellLegend = {};
vecRunTests = 1:3;
indKeep = ~any(all(isnan(sLoad.matZetaP),2),3);
vecRunTrialNums = sLoad.vecRunTrialNums(indKeep);
intTrialNum = numel(vecRunTrialNums);

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
		matFullData = sLoad.matZetaP(indKeep,:,:);
	elseif intTest == 2
		matFullData = sLoad.matKsP(indKeep,:,:);
	elseif intTest == 3
		matFullData = sLoad.matTtestP(indKeep,:,:);
	end
	for intTrialIdx = 1:intTrialNum
		intTrials = vecRunTrialNums(intTrialIdx);
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
	plot(hAx1,vecRunTrialNums,matTPR(intTest,:),'color',cellColor{intTest});
	plot(hAx1,vecRunTrialNums,matTPR_se(intTest,:,1),'--','color',cellColor{intTest},'HandleVisibility','off');
	plot(hAx1,vecRunTrialNums,matTPR_se(intTest,:,2),'--','color',cellColor{intTest},'HandleVisibility','off');
	
	plot(hAx2,vecRunTrialNums,matFPR(intTest,:),'color',cellColor{intTest});
	plot(hAx2,vecRunTrialNums,matFPR_se(intTest,:,1),'--','color',cellColor{intTest},'HandleVisibility','off');
	plot(hAx2,vecRunTrialNums,matFPR_se(intTest,:,2),'--','color',cellColor{intTest},'HandleVisibility','off');
	
	plot(hAx3,vecRunTrialNums,matAUC(intTest,:),'color',cellColor{intTest});
	plot(hAx3,vecRunTrialNums,matAUC(intTest,:)+matAUC_se(intTest,:),'--','color',cellColor{intTest},'HandleVisibility','off');
	plot(hAx3,vecRunTrialNums,matAUC(intTest,:)-matAUC_se(intTest,:),'--','color',cellColor{intTest},'HandleVisibility','off');
end
legend(hAx1,cellNames,'location','best');
ylabel(hAx1,sprintf('TPR at %s = %.3f',getGreek('alpha'),dblAlpha));
xlabel(hAx1,'Number of trials');

plot(hAx2,vecRunTrialNums([1 end]),[dblAlpha dblAlpha],'--','color',[0.5 0.5 0.5]);
legend(hAx2,[cellNames(:);{'Norm'}],'location','best');
ylabel(hAx2,sprintf('FPR at %s = %.3f',getGreek('alpha'),dblAlpha));
xlabel(hAx2,'Number of trials');

plot(hAx3,vecRunTrialNums([1 end]),[0.5 0.5],'--','color',[0.5 0.5 0.5]);
legend(hAx3,[cellNames(:);{'Chance'}],'location','best');
ylabel(hAx3,'AUC');
xlabel(hAx3,'Number of trials');
fixfig;

%save
drawnow;
export_fig(fullpath(strFigPath,['TZETA2_KS2.tif']));
export_fig(fullpath(strFigPath,['TZETA2_KS2.pdf']));

%% plot fprs
fZ = @(x) -norminv(x/2);
matPH0Z_Z = abs(fZ(sLoad.matZetaP(indKeep,:,2))-0);
matPH0Z_K = fZ(sLoad.matKsP(indKeep,:,2));
matPH0Z_T = fZ(sLoad.matTtestP(indKeep,:,2));

fP = @(x) 2-2*normcdf(x);
matPH0_Z = fP(matPH0Z_Z);
matPH0_K = fP(matPH0Z_K);
matPH0_T = fP(matPH0Z_T);

figure
subplot(2,3,1)
histogram(matPH0Z_Z(2:end,:));

subplot(2,3,2)
histogram(matPH0Z_K(2:end,:));

subplot(2,3,3)
histogram(matPH0Z_T(2:end,:));


subplot(2,3,4)
histogram(matPH0_Z(2:end,:));

subplot(2,3,5)
histogram(matPH0_K(2:end,:));

subplot(2,3,6)
histogram(matPH0_T(2:end,:));