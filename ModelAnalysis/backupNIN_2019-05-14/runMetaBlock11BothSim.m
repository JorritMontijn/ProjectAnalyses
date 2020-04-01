
%start
clear all;
%close all;

strFigType = 'Block11';
strDataDir = 'D:\Data\Results\Data11\';
strDataDir2 = 'D:\Data\Results\Data11Merge\';
sFiles = dir([strDataDir '*Subsp*.mat']);
intFiles = numel(sFiles);

%pre-allocate data
intMaxDim = inf;
cellPhiAgg = cell(1,intFiles);
cellPhiMerge = cell(1,intFiles);
cellShuff = cell(1,intFiles);
cellExpName = cell(1,intFiles);
vecRemDataSets = false(1,intFiles);
vecReps = zeros(1,intFiles);
vecdTheta = zeros(1,intFiles);

%loop through data sets
for intFileIdx=1:numel(sFiles)
	if strfind(sFiles(intFileIdx).name,'Noise088')
		intBaselineData = intFileIdx;
		vecRemDataSets(intBaselineData) = true;
	elseif isempty(strfind(sFiles(intFileIdx).name,'Noise'))
		vecRemDataSets(intFileIdx) = true; %remove simulations
		continue;
	%elseif strfind(sFiles(intFileIdx).name,'JoGu')
	%	vecRemDataSets(intFileIdx) = true; %remove simulations
	%	continue;
	end
	
	%load data
	strFile = sFiles(intFileIdx).name;
	sLoad = load([strDataDir strFile]);
	vecdTheta(intFileIdx) = sLoad.dTheta;
	
	%get relative information for raw and shuffled data
	matRawFisher = sLoad.matRawFisher';
	matShuffFisher = sLoad.matShuffFisher';
	if size(matRawFisher,1) == size(matShuffFisher,1)
		%normalize to shuffled
		matPhi = (matRawFisher - matShuffFisher);% ./ (1 - matRelShuff + 1/1000);
		
	else
		disp a
		%normalize to shuffled
		vecS = (1 - mean(matShuffFisher,1));
		matR = (1 - matRawFisher);
		matPhi = (bsxfun(@minus,vecS,matR));% ./ matR;
	end
	
	%check if min dim
	intSize = size(matPhi,2);
	if intSize < intMaxDim && ~isinf(intMaxDim),
		vecRemDataSets(intFileIdx) = true;
		continue;
	end
	
	%get relative information for raw and shuffled data
	matRawFisherMerged = sLoad.matRawFisherMerged';
	matShuffFisherMerged = sLoad.matShuffFisherMerged';
	if size(matRawFisherMerged,1) == size(matShuffFisherMerged,1)
		%normalize to shuffled
		matPhiMerge = (matRawFisherMerged - matShuffFisherMerged);% ./ (1 - matRelShuff + 1/1000);
		
	else
		disp a
		%normalize to shuffled
		vecS = (1 - mean(matShuffFisherMerged,1));
		matR = (1 - matRawFisherMerged);
		matPhiMerge = (bsxfun(@minus,vecS,matR));% ./ matR;
	end
	
	%check if min dim
	intSizeMerge = size(matPhiMerge,2);
	if intSizeMerge < intMaxDim && ~isinf(intMaxDim),
		vecRemDataSets(intFileIdx) = true;
		continue;
	end
	%}
	%process
	%matPhi = (sLoad.matPredI  - sLoad.matRandI) ./ (sLoad.matFullI - sLoad.matRandI);
	cellPhiAgg{intFileIdx} = nansum(matPhi(:,1:min(intMaxDim,intSize)),2)/min(intMaxDim,intSize);
	cellPhiMerge{intFileIdx} = nansum(matPhiMerge(:,1:min(intMaxDim,intSizeMerge)),2)/min(intMaxDim,intSizeMerge);
	cellExpName{intFileIdx} = sLoad.strName;
	
end

%get baseline
vecBaselineData = cellPhiAgg{intBaselineData};
vecBaselineMerge = cellPhiMerge{intBaselineData};
dThetaBase = vecdTheta(intBaselineData);
dblBaseEquivSD = roundi(dThetaBase*(sqrt(2)/2),1);
cellPhiAgg(vecRemDataSets) = [];
cellPhiMerge(vecRemDataSets) = [];
cellExpName(vecRemDataSets) = [];
vecdTheta(vecRemDataSets) = [];

%% reorder
[cellExpName,vecReorder] = sort(cellExpName);
cellPhiAgg = cellPhiAgg(vecReorder);
cellPhiMerge = cellPhiMerge(vecReorder);
vecdTheta = vecdTheta(vecReorder);
celldTheta = cellfun(@num2str,num2cell(vecdTheta)','UniformOutput',false);

%transform dTheta to equivalent wiggling; x*sqrt(2)/2
cellEquivalentSD = cellfun(@num2str,num2cell(roundi(vecdTheta'*(sqrt(2)/2),1))','UniformOutput',false);

%normal
vecMean = cellfun(@mean,cellPhiAgg);
vecSD = cellfun(@std,cellPhiAgg);
vecIters = cellfun(@numel,cellPhiAgg);
vecSEM = vecSD./sqrt(vecIters);
[h,dblP_zero]=ttest(vecMean);

dblBaseMean = mean(vecBaselineData);
[h,dblP_base]=ttest(vecMean,dblBaseMean);

%merge
vecMeanMerge = cellfun(@mean,cellPhiMerge);
vecSDMerge = cellfun(@std,cellPhiMerge);
vecItersMerge = cellfun(@numel,cellPhiMerge);
vecSEMMerge = vecSDMerge./sqrt(vecItersMerge);
[h,dblP_zeroMerge]=ttest(vecMeanMerge);

dblBaseMeanMerge = mean(vecBaselineMerge);
[h,dblP_baseMerge]=ttest(vecMeanMerge,dblBaseMeanMerge);


%% plot
intN = numel(cellPhiAgg);
figure
%subplot(2,3,[1 2 4 5])
hold on
plot([0.5 intN+0.5],[0 0],'k--');
errorbar((1:intN)-0.1,repmat(mean(vecBaselineData),[1 intN]),repmat(std(vecBaselineData),[1 intN]),'c--');
errorbar((1:intN)+0.1,repmat(mean(vecBaselineMerge),[1 intN]),repmat(std(vecBaselineMerge),[1 intN]),'m--');
text((intN)+0.2,mean(vecBaselineMerge),sprintf(' %.1f',dblBaseEquivSD),'Color','magenta','FontSize',14);
title(['Fisher I, data&merged, \phi ' sprintf('(AUC) up to subspace rank %d; mean=%.3f; t-test vs 0,p=%.6f; t-test vs baseline (%.3f),p=%.6f',intMaxDim,mean(vecMean),dblP_zero,dblBaseMean,dblP_base)]);
errorbar((1:intN)-0.1,vecMean,vecSD,'xb');
errorbar((1:intN)+0.1,vecMeanMerge,vecSDMerge,'xr');
legend({'DiCo=0','DiCo=behav','DiCo=behav,merge','DiCo=data','DiCo=merge'},'Location','best');
text((1:intN)+0.2,vecMeanMerge,cellEquivalentSD,'Color','red','FontSize',14);
fixfig;
set(gca,'xtick',1:intN,'xticklabel',cellExpName,'XTickLabelRotation',45)
ylabel('Norm. diff. corr. strength (\phi)')

%save fig
drawnow;
jFig = get(handle(gcf), 'JavaFrame');
jFig.setMaximized(true);
figure(gcf);
drawnow;
strFig = [strFigType 'AggregateGraphDim'  num2str(intMaxDim) '_' getDate];
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

