
%start
clear all;
close all;

strDataDir = 'D:\Data\Results\Data8Top50\';
sFiles = dir([strDataDir '*Subsp*.mat']);
intFiles = numel(sFiles);

%pre-allocate data
intMaxDim = inf;
cellPhiAgg = cell(1,intFiles);
cellRaw = cell(1,intFiles);
cellShuff = cell(1,intFiles);
cellExpName = cell(1,intFiles);
vecRemDataSets = false(1,intFiles);
vecReps = zeros(1,intFiles);
%loop through data sets
for intFileIdx=1:numel(sFiles)
	if strfind(sFiles(intFileIdx).name,'Noise063')
		intBaselineData = intFileIdx;
		vecRemDataSets(intBaselineData) = true;
	end
	
	%load data
	strFile = sFiles(intFileIdx).name;
	sLoad = load([strDataDir strFile]);
	
	%get relative information for raw and shuffled data
	matRelRaw = sLoad.matRawCos';
	matRelShuff = sLoad.matShuffCos';
	if size(matRelRaw,1) == size(matRelShuff,1)
		%normalize to shuffled
		matPhi = (matRelRaw - matRelShuff);% ./ (1 - matRelShuff + 1/1000);
		
	else
		disp a
		%normalize to shuffled
		vecS = (1 - mean(matRelShuff,1));
		matR = (1 - matRelRaw);
		matPhi = (bsxfun(@minus,vecS,matR));% ./ matR;
	end
	
	%check if min dim
	intSize = size(matPhi,2);
	if intSize < intMaxDim && ~isinf(intMaxDim),
		vecRemDataSets(intFileIdx) = true;
		continue;
	end
	
	%process
	%matPhi = (sLoad.matPredI  - sLoad.matRandI) ./ (sLoad.matFullI - sLoad.matRandI);
	cellPhiAgg{intFileIdx} = sum(matPhi(:,1:min(intMaxDim,intSize)),2)/min(intMaxDim,intSize);
	cellRaw{intFileIdx} = matRelRaw;
	cellShuff{intFileIdx} = matRelShuff;
	cellExpName{intFileIdx} = sLoad.strName;
	
end

%get baseline
vecBaselineData = cellPhiAgg{intBaselineData};
cellPhiAgg(vecRemDataSets) = [];
cellExpName(vecRemDataSets) = [];

%plot
vecMean = cellfun(@mean,cellPhiAgg);
vecSD = cellfun(@std,cellPhiAgg);
vecIters = cellfun(@numel,cellPhiAgg);
vecSEM = vecSD./sqrt(vecIters);
[h,dblP_zero]=ttest(vecMean);

dblBaseMean = mean(vecBaselineData);
[h,dblP_base]=ttest(vecMean,dblBaseMean);

%% plot
intN = numel(cellPhiAgg);
figure
%subplot(2,3,[1 2 4 5])
hold on
plot([0.5 intN+0.5],[0 0],'k--');
errorbar(1:intN,repmat(mean(vecBaselineData),[1 intN]),repmat(std(vecBaselineData),[1 intN]));
title(['Top 50% neurons, \phi ' sprintf('(AUC) up to subspace rank %d; mean=%.3f; t-test vs 0,p=%.6f; t-test vs baseline (%.3f),p=%.6f',intMaxDim,mean(vecMean),dblP_zero,dblBaseMean,dblP_base)]);
legend({'DiCo=0','DiCo=behav'});
errorbar(1:intN,vecMean,vecSD,'xk');
fixfig;
set(gca,'xtick',1:intN,'xticklabel',cellExpName,'XTickLabelRotation',45)
ylabel('Norm. diff. corr. strength (\phi)')

%save fig
drawnow;
jFig = get(handle(gcf), 'JavaFrame');
jFig.setMaximized(true);
figure(gcf);
drawnow;
strFig = ['AggregateGraphDim' num2str(intMaxDim) '_' getDate];
export_fig([strDataDir strFig '.tif']);
export_fig([strDataDir strFig '.pdf']);

%% second plot
%{
vecAllPoints = cell2mat(cellPhiAgg');
[h,dblP]=ttest(vecAllPoints);
figure
dblBinSize = 0.1;
vecBins = [-1:dblBinSize:1];
vecPlot = vecBins(2:end) - dblBinSize/2;
[vecCounts,edges] = histcounts(vecAllPoints,vecBins);
plot(vecPlot,vecCounts);
hold on
plot(mean(vecAllPoints)*[1 1],[0 max(get(gca,'ylim'))],'b--');
hold off
ylabel('Number of data points (count)');
xlabel('Norm. diff. corr. strength (\phi)')
title(['\phi ' sprintf('(AUC) up to rank %d; mean=%.3f; t-test,p=%.3f',intMaxDim,mean(vecAllPoints),dblP)]);
fixfig
%}

