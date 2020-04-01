
%start
clear all;
%close all;
boolPlotSims = true; %true for sims, false for monkeys
strFigType = 'Block8';
strDataDir = 'D:\Data\Results\Data8\';
sFiles = dir([strDataDir '*Block8*.mat']);
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
intBaselineData = [];
vecNumNeurons = [];
vecNN2 = [];
vecNoise = [];
for intFileIdx=1:numel(sFiles)
	if strfind(sFiles(intFileIdx).name,'Noise063')
		intBaselineData = intFileIdx;
		vecRemDataSets(intBaselineData) = true;
	elseif boolPlotSims && strfind(sFiles(intFileIdx).name,'Noise')
		vecRemDataSets(intFileIdx) = false; %remove simulations
		%continue;
	elseif ~boolPlotSims && ~isempty(strfind(sFiles(intFileIdx).name,'monyet')) || ~isempty(strfind(sFiles(intFileIdx).name,'cadet'))
		vecRemDataSets(intFileIdx) = false; %remove data
		%continue;
	else
		vecRemDataSets(intFileIdx) = true; %remove simulations
		continue;
	end
	
	%load data
	strFile = sFiles(intFileIdx).name;
	sLoad = load([strDataDir strFile]);
	vecdTheta(intFileIdx) = sLoad.dTheta;
	vecNN2(intFileIdx) = sLoad.intNeurons;
	if sLoad.intNeurons < 60
		vecRemDataSets(intFileIdx) = true;
		continue
	end
	vecNoise(intFileIdx) = sLoad.dblNoise;
	
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
	
	%% merged
	%get relative information for raw and shuffled data
	matRawCosMerged = sLoad.matRawCosMerged';
	matShuffCosMerged = sLoad.matShuffCosMerged';
	if size(matRawCosMerged,1) == size(matShuffCosMerged,1)
		%normalize to shuffled
		matPhiMerge = (matRawCosMerged - matShuffCosMerged);% ./ (1 - matRelShuff + 1/1000);
		
	else
		disp a
		%normalize to shuffled
		vecS = (1 - mean(matShuffCosMerged,1));
		matR = (1 - matRawCosMerged);
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
	vecNumNeurons(intFileIdx) = intSize;
end

%get baseline
if ~isempty(intBaselineData)
	vecBaselineData = cellPhiAgg{intBaselineData};
	vecBaselineMerge = cellPhiMerge{intBaselineData};
	dThetaBase = vecdTheta(intBaselineData);
	dblBaseEquivSD = roundi(dThetaBase*(sqrt(2)/2),1);
end
cellPhiAgg(vecRemDataSets) = [];
cellPhiMerge(vecRemDataSets) = [];
cellExpName(vecRemDataSets) = [];
vecdTheta(vecRemDataSets) = [];
vecNumNeurons(vecRemDataSets) = [];
vecNN2(vecRemDataSets) = [];
vecNoise(vecRemDataSets) = [];

%% reorder
[vecNoise,vecReorder] = sort(vecNoise);
cellExpName = cellExpName(vecReorder);
cellPhiAgg = cellPhiAgg(vecReorder);
cellPhiMerge = cellPhiMerge(vecReorder);
vecdTheta = vecdTheta(vecReorder);
celldTheta = cellfun(@num2str,num2cell(vecdTheta)','UniformOutput',false);
vecNumNeurons = vecNumNeurons(vecReorder);
vecNN2 = vecNN2(vecReorder);


%transform dTheta to equivalent wiggling; x*sqrt(2)/2
cellEquivalentSD = cellfun(@num2str,num2cell(roundi(vecdTheta'*(sqrt(2)/2),1))','UniformOutput',false);

%normal
vecMean = cellfun(@mean,cellPhiAgg);
vecSD = cellfun(@std,cellPhiAgg);
vecIters = cellfun(@numel,cellPhiAgg);
vecSEM = vecSD./sqrt(vecIters);
[h,dblP_zero]=ttest(vecMean);

%merge
vecMeanMerge = cellfun(@mean,cellPhiMerge);
vecSDMerge = cellfun(@std,cellPhiMerge);
vecItersMerge = cellfun(@numel,cellPhiMerge);
vecSEMMerge = vecSDMerge./sqrt(vecItersMerge);
[h,dblP_zeroMerge]=ttest(vecMeanMerge);

%baseline
if ~isempty(intBaselineData)
	dblBaseMean = mean(vecBaselineData);
	[h,dblP_base]=ttest(vecMean,dblBaseMean);
	
	dblBaseMeanMerge = mean(vecBaselineMerge);
	[h,dblP_baseMerge]=ttest(vecMeanMerge,dblBaseMeanMerge);
end

%% plot
vecX = vecNoise;
intN = numel(cellPhiAgg);
figure
%subplot(2,3,[1 2 4 5])
hold on
plot([min(vecX) max(vecX)],[0 0],'k--');
if ~isempty(intBaselineData)
errorbar(vecX-0.01,repmat(mean(vecBaselineData),[1 intN]),repmat(std(vecBaselineData),[1 intN]),'c--');
errorbar(vecX+0.01,repmat(mean(vecBaselineMerge),[1 intN]),repmat(std(vecBaselineMerge),[1 intN]),'m--');
%text((intN)+0.02,mean(vecBaselineMerge),sprintf(' %.1f',dblBaseEquivSD),'Color','magenta','FontSize',14);
title(['f'' in \Sigma, data&merged, \phi ' sprintf('(AUC) up to subspace rank %d; mean=%.3f; t-test vs 0,p=%.6f; t-test vs baseline (%.3f),p=%.6f',intMaxDim,mean(vecMean),dblP_zero,dblBaseMean,dblP_base)]);
legend({'DiCo=0','DiCo=behav','DiCo=behav,merge','DiCo=data','DiCo=merge'},'Location','Best');
else
	title(['f'' in \Sigma, data&merged, \phi ' sprintf('(AUC) up to subspace rank %d; mean=%.3f; t-test vs 0,p=%.6f',intMaxDim,mean(vecMean),dblP_zero)]);
end
errorbar(vecX-0.01,vecMean,vecSD,'xb');
%errorbar(vecX+0.01,vecMeanMerge,vecSDMerge,'xr');
%text(vecX+0.02,vecMeanMerge,cellEquivalentSD,'Color','red','FontSize',14);
fixfig;
%set(gca,'xtick',vecX,'xticklabel',cellExpName,'XTickLabelRotation',45)
ylabel('Norm. diff. corr. strength (\phi)')
xlabel('Noise injection (degs)');
xlim([min(vecX)-0.1 max(vecX)+0.1]);

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

vecAllPoints = cell2mat(cellPhiAgg');
[h,dblP]=ttest(vecMean);
figure
dblBinSize = 0.05;
vecBins = [(-0.2-dblBinSize/2):dblBinSize:(0.2+dblBinSize/2)];
vecPlot = vecBins(2:end) - dblBinSize/2;
[vecCounts,edges] = histcounts(vecMean,vecBins);
hold on
%plot([0 0],[0 4],'k--');
plot(vecPlot,vecCounts,'x-');
hold off
ylabel('Number of data points (count)');
xlabel('Norm. diff. corr. strength (\phi)')
title(['\phi ' sprintf('(AUC) up to rank %d; mean=%.3f; t-test,p=%.3f',intMaxDim,mean(vecMean),dblP)]);
fixfig

drawnow;
strFig = [strFigType 'AggregateHistoGraphDim'  num2str(intMaxDim) '_' getDate];
export_fig([strDataDir strFig '.tif']);
export_fig([strDataDir strFig '.pdf']);
