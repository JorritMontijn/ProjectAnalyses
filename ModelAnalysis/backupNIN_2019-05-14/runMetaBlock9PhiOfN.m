
%start
clearvars -except strDataDir intSubSample* vecNeuronN vecDoShuff vecRunSims intLoadSim bool* vecRunAreas
	%close all;

strFigType = 'Block9';
if ~exist('strDataDir','var'),strDataDir = 'D:\Data\Results\Data9\';end
strMetaFigDir = 'D:\Data\Results\Data9\';
sFiles = dir([strDataDir '*Block9*.mat']);
intFiles = numel(sFiles);

%pre-allocate data
intMaxDim = inf;
matPhiAgg = [];
matPhiAggMerged = [];
vecRemDataSets = false(1,intFiles);
vecNeurons = zeros(1,intFiles);
vecdTheta = zeros(1,intFiles);

%loop through data sets
intBaselineData = [];
for intFileIdx=1:numel(sFiles)
	
	%load data
	strFile = sFiles(intFileIdx).name;
	sLoad = load([strDataDir strFile]);
	strName = sLoad.strName;
	vecdTheta(intFileIdx) = sLoad.dTheta;
	vecNeurons(intFileIdx) = sLoad.intNeurons;
	
	%get relative information for raw and shuffled data
	matRelRaw = sLoad.matRawCos';
	matRelShuff = sLoad.matShuffCos';
	matRelRawMerged = sLoad.matRawCosMerged';
	matRelShuffMerged = sLoad.matShuffCosMerged';
	
	%normalize to shuffled
	matPhi = (matRelRaw - matRelShuff);% ./ (1 - matRelShuff + 1/1000);
	matPhiMerged = (matRelRawMerged - matRelShuffMerged);% ./ (1 - matRelShuff + 1/1000);
		
	%check if min dim
	intSize = size(matPhi,2);
	if intSize < intMaxDim && ~isinf(intMaxDim),
		vecRemDataSets(intFileIdx) = true;
		continue;
	end

	%}
	%process
	matPhiAgg(intFileIdx,:) = nansum(matPhi(:,1:min(intMaxDim,intSize)),2)/min(intMaxDim,intSize);
	matPhiAggMerged(intFileIdx,:) = nansum(matPhiMerged(:,1:min(intMaxDim,intSize)),2)/min(intMaxDim,intSize);
	
end

%% remove data sets
matPhiAgg(vecRemDataSets,:) = [];
matPhiAggMerged(vecRemDataSets,:) = [];
vecdTheta(vecRemDataSets) = [];
vecNeurons(vecRemDataSets) = [];

%resort
[vecNeurons,vecReorder] = sort(vecNeurons,'ascend');
matPhiAgg = matPhiAgg(vecReorder,:);
matPhiAggMerged = matPhiAggMerged(vecReorder,:);

%transform dTheta to equivalent wiggling; x*sqrt(2)/2
cellEquivalentSD = cellfun(@num2str,num2cell(roundi(vecdTheta'*(sqrt(2)/2),1))','UniformOutput',false);

%normal
vecMean = xmean(matPhiAgg,2);
vecSD = xstd(matPhiAgg,2);
vecIters = size(matPhiAgg,2);
vecSEM = vecSD./sqrt(vecIters);
[h,dblP_zero]=ttest(vecMean);

%merge
vecMeanMerged = xmean(matPhiAggMerged,2);
vecSDMerged = xstd(matPhiAggMerged,2);
vecItersMerged = size(matPhiAggMerged,2);
vecSEMMerged = vecSDMerged./sqrt(vecIters);
[h,dblP_zeroMerged]=ttest(vecMeanMerged);

%% plot
intN = numel(vecMeanMerged);
figure
%subplot(2,3,[1 2 4 5])
hold on
title(['f'' in \Sigma, Noise ' getFlankedBy(strName,'Noise','')]);
errorbar(vecNeurons-1,vecMean,vecSD,'xb-');
errorbar(vecNeurons+1,vecMeanMerged,vecSDMerged,'xr-');
xlabel('Subsampled population size (N)');
ylabel('Norm. diff. corr. strength (\phi)')
legend({'Unmerged','Merged'},'Location','Best');
fixfig;

%% save fig
drawnow;
strFig = [strFigType 'VaryingN_' strName '_' getDate];
export_fig([strMetaFigDir strFig '.tif']);
export_fig([strMetaFigDir strFig '.pdf']);

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

