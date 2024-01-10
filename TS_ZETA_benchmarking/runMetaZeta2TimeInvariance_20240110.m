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
sDir=dir([strDataPath 'Zeta2TimeInvariance.mat']);
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
cellNames = {'ZETA2','Ttest2'};
cellLegend = {};
vecRunTests = 1:2;
indKeep = ~any(all(isnan(sLoad.matZeta2P),2),3);
vecTimeShifts = sLoad.vecTimeShifts(indKeep);
intTimeShiftNum = numel(vecTimeShifts);

matTPR = nan(numel(vecRunTests),intTimeShiftNum);
matTPR_se = nan(numel(vecRunTests),intTimeShiftNum,2);
matFPR = nan(numel(vecRunTests),intTimeShiftNum);
matFPR_se = nan(numel(vecRunTests),intTimeShiftNum,2);
matAUC = nan(numel(vecRunTests),intTimeShiftNum);
matAUC_se = nan(numel(vecRunTests),intTimeShiftNum);
hold on;
dblAlpha = 0.05;
for intTest=vecRunTests
	if intTest == 1
		matFullData = sLoad.matZeta2P(indKeep,:,:);
	elseif intTest == 2
		matFullData = sLoad.matTtest2P(indKeep,:,:);
	end
	for intShiftIdx = 1:intTimeShiftNum
		dblShift = vecTimeShifts(intShiftIdx);
		matData = squeeze(matFullData(intShiftIdx,:,:));
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
		matTPR(intTest,intShiftIdx) = dblTPR;
		matTPR_se(intTest,intShiftIdx,:) = dblTPR_se;
		matFPR(intTest,intShiftIdx) = dblFPR;
		matFPR_se(intTest,intShiftIdx,:) = dblFPR_se;
		
		%plot(vecFP,vecTP,'Color',cellColor{intTest});
		
		[dblAUC,Aci,Ase] = getAuc(vecShuffP,vecRealP,0.05,'mann-whitney');
		matAUC(intTest,intShiftIdx) = dblAUC;
		matAUC_se(intTest,intShiftIdx) = Ase;
		
		%cellLegend(end+1) = {sprintf('%s, AUC=%.3f',cellNames{intTest},dblAUC)};
	end
	% plot
	plot(hAx1,vecTimeShifts,matTPR(intTest,:),'color',cellColor{intTest});
	plot(hAx1,vecTimeShifts,matTPR_se(intTest,:,1),'--','color',cellColor{intTest},'HandleVisibility','off');
	plot(hAx1,vecTimeShifts,matTPR_se(intTest,:,2),'--','color',cellColor{intTest},'HandleVisibility','off');
	
	plot(hAx2,vecTimeShifts,matFPR(intTest,:),'color',cellColor{intTest});
	plot(hAx2,vecTimeShifts,matFPR_se(intTest,:,1),'--','color',cellColor{intTest},'HandleVisibility','off');
	plot(hAx2,vecTimeShifts,matFPR_se(intTest,:,2),'--','color',cellColor{intTest},'HandleVisibility','off');
	
	plot(hAx3,vecTimeShifts,matAUC(intTest,:),'color',cellColor{intTest});
	plot(hAx3,vecTimeShifts,matAUC(intTest,:)+matAUC_se(intTest,:),'--','color',cellColor{intTest},'HandleVisibility','off');
	plot(hAx3,vecTimeShifts,matAUC(intTest,:)-matAUC_se(intTest,:),'--','color',cellColor{intTest},'HandleVisibility','off');
end
legend(hAx1,cellNames,'location','best');
ylabel(hAx1,sprintf('TPR at %s = %.3f',getGreek('alpha'),dblAlpha));
xlabel(hAx1,'Time shift (s)');

plot(hAx2,vecTimeShifts([1 end]),[dblAlpha dblAlpha],'--','color',[0.5 0.5 0.5]);
legend(hAx2,[cellNames(:);{'Norm'}],'location','best');
ylabel(hAx2,sprintf('FPR at %s = %.3f',getGreek('alpha'),dblAlpha));
xlabel(hAx2,'Time shift (s)');

plot(hAx3,vecTimeShifts([1 end]),[0.5 0.5],'--','color',[0.5 0.5 0.5]);
legend(hAx3,[cellNames(:);{'Chance'}],'location','best');
ylabel(hAx3,'AUC');
xlabel(hAx3,'Time shift (s)');
fixfig;

%save
drawnow;
export_fig(fullpath(strFigPath,['ZETA2_TimeInv.jpg']));
export_fig(fullpath(strFigPath,['ZETA2_TimeInv.pdf']));

%% plot fprs
fZ = @(x) -norminv(x/2);
matPH0Z_Z = abs(fZ(sLoad.matZeta2P(indKeep,:,2))-0);
matPH0Z_T = fZ(sLoad.matTtest2P(indKeep,:,2));

fP = @(x) 2-2*normcdf(x);
matPH0_Z = fP(matPH0Z_Z);
matPH0_T = fP(matPH0Z_T);

figure
subplot(2,3,1)
histogram(matPH0Z_Z(2:end,:));

subplot(2,3,3)
histogram(matPH0Z_T(2:end,:));


subplot(2,3,4)
histogram(matPH0_Z(2:end,:));

subplot(2,3,6)
histogram(matPH0_T(2:end,:));