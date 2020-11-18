
%start
clear all;
%close all;

strFigType = 'Block8';
strDataDir = ['F:\Data\Results\SimFigs\Data8\'];
strDataDir2 = 'D:\Data\Results\Data8Merge\';
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
vecNoise = zeros(1,intFiles);

%loop through data sets
intBaselineData = [];
for intFileIdx=1:numel(sFiles)
	%if ~isempty(strfind(sFiles(intFileIdx).name,'Noise0_')) || ~isempty(strfind(sFiles(intFileIdx).name,'Noise063_')) || ~isempty(strfind(sFiles(intFileIdx).name,'Noise5_'))
	if ~isempty(strfind(sFiles(intFileIdx).name,'Noise'))
		if ~isempty(strfind(sFiles(intFileIdx).name,'Noise063_')) || ...
				...%~isempty(strfind(sFiles(intFileIdx).name,'Noise0_')) || ...
				~isempty(strfind(sFiles(intFileIdx).name,'Noise088_')) || ...
				~isempty(strfind(sFiles(intFileIdx).name,'Noise5_'))
			
			sFiles(intFileIdx).name
			vecRemDataSets(intFileIdx) = true; %remove non-simulations
			continue;
		end
	else
		vecRemDataSets(intFileIdx) = true; %remove non-simulations
		continue;
	end
	
	%load data
	strFile = sFiles(intFileIdx).name;
	sLoad = load([strDataDir strFile]);
	vecdTheta(intFileIdx) = sLoad.dTheta;
	
	%get noise
	strNoise = getFlankedBy(sLoad.strName,'Noise','');
	vecNum = arrayfun(@str2double,strNoise);
	intPos= find(isnan(vecNum),1,'first');
	if isempty(intPos),intPos = numel(vecNum)+1;end
	dblNoise = str2double(strNoise(1:(intPos-1)));
	if dblNoise > 9,dblNoise=dblNoise/10;end
	if vecNum(1) == 0 && ~all(vecNum==0)
		dblNoise = dblNoise/(10^ceil(log10(dblNoise)));
	end
	vecNoise(intFileIdx) = dblNoise;
	
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
	
end
vecRemDataSets(cellfun(@isempty,cellPhiAgg)) = true;

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
vecNoise(vecRemDataSets) = [];

%% reorder
[cellExpName,vecReorder] = sort(cellExpName);
cellPhiAgg = cellPhiAgg(vecReorder);
cellPhiMerge = cellPhiMerge(vecReorder);
vecdTheta = vecdTheta(vecReorder);
vecNoise = vecNoise(vecReorder);
celldTheta = cellfun(@num2str,num2cell(vecdTheta)','UniformOutput',false);

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
%fit logistic function
vecP0 = [0.5 1 0.5 0];
vecX = vecNoise*2.56;
vecY = vecMean;
sOptions = optimoptions('lsqcurvefit','MaxFunctionEvaluations',5000,'MaxIterations',5000);

[vecFitP] = lsqcurvefit(@logisticfitx0,vecP0,vecX,vecY,[],[],sOptions);
vecPlotFit = 0:0.01:max(vecX);
vecFitY = logisticfitx0(vecFitP,vecPlotFit);

%calc equivalence
dblData =0.1359;
L = vecFitP(1);
k = vecFitP(2);
x0 = vecFitP(3);
y0 = vecFitP(4);
dblEqNoise = x0 + log(1/(0.5+ ((dblData-y0)/L)/2) -1)/-k;

%% plot
intN = numel(cellPhiAgg);
figure
%subplot(2,3,[1 2 4 5])
hold on
if exist('vecPlotFit','var') && ~isempty(vecPlotFit),plot(vecPlotFit,vecFitY);end
if ~isempty(intBaselineData)
errorbar(vecX,repmat(mean(vecBaselineData),[1 intN]),repmat(std(vecBaselineData),[1 intN]),'c--');
title(['f'' in \Sigma, equiv=' sprintf('%.3f',dblEqNoise) ',\phi ' sprintf('(AUC) up to subspace rank %d; mean=%.3f; t-test vs 0,p=%.6f; t-test vs baseline (%.3f),p=%.6f',intMaxDim,mean(vecMean),dblP_zero,dblBaseMean,dblP_base)]);
legend({'Logistic fit','DiCo=behav','DiCo=data'},'Location','Best');
else
	title(['f'' in \Sigma, data&merged, \phi ' sprintf('(AUC) up to subspace rank %d; mean=%.3f; t-test vs 0,p=%.6f',intMaxDim,mean(vecMean),dblP_zero)]);
end
%errorbar(vecX,vecY,vecE,'xb');
scatter(vecX,vecY,'xb');
%errorbar(vecX,vecMeanMerge,vecSDMerge,'xr');
scatter(dblEqNoise,dblData,'xk');
%text(vecX,vecMeanMerge,cellEquivalentSD,'Color','red','FontSize',14);
fixfig;
[vecSortX,vecReorder]=sort(vecX);
%set(gca,'xtick',vecSortX,'xticklabel',cellExpName(vecReorder),'XTickLabelRotation',45)
ylabel('Norm. diff. corr. strength (\phi)')
xlabel('Orientation Discrimination Threshold (degs)')

%save fig
drawnow;
jFig = get(handle(gcf), 'JavaFrame');
jFig.setMaximized(true);
figure(gcf);
drawnow;
strFig = [strFigType 'AggregateGraphSimDim'  num2str(intMaxDim) '_' getDate];
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

