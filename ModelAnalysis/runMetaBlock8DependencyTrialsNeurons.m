
%start
clear all;
close all;

strFigType = 'Dep8';
strDataDir = 'D:\Data\Results\Data8\';
sFiles = dir([strDataDir '*Block*.mat']);
intFiles = numel(sFiles);

%pre-allocate data
intMaxDim = inf;
cellPhiAgg = cell(1,intFiles);
cellPhiMerge = cell(1,intFiles);
cellRaw = cell(1,intFiles);
cellShuff = cell(1,intFiles);
cellExpName = cell(1,intFiles);
vecRemDataSets = false(1,intFiles);
vecRepetitions = zeros(1,intFiles);
vecNeurons = zeros(1,intFiles);
vecdTheta = zeros(1,intFiles);

%loop through data sets
for intFileIdx=1:numel(sFiles)
	if strfind(sFiles(intFileIdx).name,'Noise088')
		intBaselineData = intFileIdx;
		vecRemDataSets(intBaselineData) = true;
	elseif strfind(sFiles(intFileIdx).name,'Noise')
		vecRemDataSets(intFileIdx) = true; %remove simulations
		continue;
	end
	
	%load data
	strFile = sFiles(intFileIdx).name;
	sLoad = load([strDataDir strFile]);
	vecdTheta(intFileIdx) = sLoad.dTheta;
	
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
	vecRepetitions(intFileIdx) = sLoad.intRepetitions;
	vecNeurons(intFileIdx) = sLoad.intNeurons;
	vecdTheta(intFileIdx) = sLoad.dTheta;
end

%get baseline
vecBaselineData = cellPhiAgg{intBaselineData};
dThetaBase = vecdTheta(intBaselineData);
cellPhiAgg(vecRemDataSets) = [];
cellPhiMerge(vecRemDataSets) = [];
cellExpName(vecRemDataSets) = [];
vecdTheta(vecRemDataSets) = [];
vecNeurons(vecRemDataSets) = [];
vecRepetitions(vecRemDataSets) = [];

%% reorder
[cellExpName,vecReorder] = sort(cellExpName);
cellPhiAgg = cellPhiAgg(vecReorder);
cellPhiMerge = cellPhiMerge(vecReorder);
vecdTheta = vecdTheta(vecReorder);
vecNeurons = vecNeurons(vecReorder);
vecRepetitions = vecRepetitions(vecReorder);
celldTheta = cellfun(@num2str,num2cell(vecdTheta)','UniformOutput',false);

%normal
vecMean = cellfun(@mean,cellPhiAgg);
vecSD = cellfun(@std,cellPhiAgg);
vecIters = cellfun(@numel,cellPhiAgg);
vecSEM = vecSD./sqrt(vecIters);
[h,dblP_zero]=ttest(vecMean);

dblBaseMean = mean(vecBaselineData);
[h,dblP_base]=ttest(vecMean,dblBaseMean);

%merged
vecMeanMerged = cellfun(@mean,cellPhiMerge);
vecSDMerged = cellfun(@std,cellPhiMerge);
vecItersMerged = cellfun(@numel,cellPhiMerge);
vecSEMMerged = vecSDMerged./sqrt(vecItersMerged);
[h,dblP_zero]=ttest(vecMeanMerged);

%% plot
intN = numel(cellPhiAgg);
figure
subplot(2,2,1)
[r,dblP]=corr(vecdTheta(:),vecMean(:));
[vecP,s] = polyfit(vecdTheta,vecMean,1);
vecX1 = min(vecdTheta):(range(vecdTheta)/100):max(vecdTheta);
hold on
plot(vecX1,vecX1*vecP(1)+vecP(2));
[y,delta] = polyval(vecP,vecX1,s);
plot(vecX1,y,'r');
plot(vecX1,y-delta,'r--');
plot(vecX1,y+delta,'r--');
scatter(vecdTheta,vecMean,'bo')
hold off
xlim([0 max(get(gca,'xlim'))]);
title(['f'' in \Sigma, \phi (AUC) up to rank' sprintf('%d; Pearson corr, r=%.3f, p=%.3f',intMaxDim,r,dblP)])
xlabel('Difference in orientation (\delta\theta)')
ylabel('Diff. corr. strength, \phi (AUC)');
fixfig

subplot(2,2,2)
[r,dblP]=corr(vecNeurons(:),vecMean(:));
[vecP,s] = polyfit(vecNeurons,vecMean,1);
vecX2 = min(vecNeurons):(range(vecNeurons)/100):max(vecNeurons);
hold on
plot(vecX2,vecX2*vecP(1)+vecP(2));
[y,delta] = polyval(vecP,vecX2,s);
plot(vecX2,y,'r');
plot(vecX2,y-delta,'r--');
plot(vecX2,y+delta,'r--');
scatter(vecNeurons,vecMean,'bo')
hold off
xlim([0 max(get(gca,'xlim'))]);
title(['f'' in \Sigma, \phi (AUC) up to rank' sprintf('%d; Pearson corr, r=%.3f, p=%.3f',intMaxDim,r,dblP)])
xlabel('Number of neurons')
ylabel('Diff. corr. strength, \phi (AUC)');
fixfig


subplot(2,2,3)
[r,dblP]=corr(vecRepetitions(:),vecMean(:));
[vecP,s] = polyfit(vecRepetitions,vecMean,1);
vecX3 = min(vecRepetitions):(range(vecRepetitions)/100):max(vecRepetitions);
hold on
plot(vecX3,vecX3*vecP(1)+vecP(2));
[y,delta] = polyval(vecP,vecX3,s);
plot(vecX3,y,'r');
plot(vecX3,y-delta,'r--');
plot(vecX3,y+delta,'r--');
scatter(vecRepetitions,vecMean,'bo')
hold off
xlim([0 max(get(gca,'xlim'))]);
title(['f'' in \Sigma, \phi (AUC) up to rank' sprintf('%d; Pearson corr, r=%.3f, p=%.3f',intMaxDim,r,dblP)])
xlabel('Number of repetitions')
ylabel('Diff. corr. strength, \phi (AUC)');
fixfig

subplot(2,2,4)
[rUM,dblPUM]=corr(vecMean(:),vecMeanMerged(:));
[vecP,s] = polyfit(vecMean,vecMeanMerged,1);
vecX3 = min(vecMean):(range(vecMean)/100):max(vecMean);
hold on
plot(vecX3,vecX3*vecP(1)+vecP(2));
[y,delta] = polyval(vecP,vecX3,s);
plot(vecX3,y,'r');
plot(vecX3,y-delta,'r--');
plot(vecX3,y+delta,'r--');
scatter(vecMean,vecMeanMerged,'bo')
hold off
title(['Merged vs unmerged, \phi (AUC); ' sprintf('Pearson corr, r=%.3f, p=%.3f',rUM,dblPUM)])
xlabel('Unmerged \phi (AUC)')
ylabel('Merged \phi (AUC)');
fixfig



%save fig
drawnow;
jFig = get(handle(gcf), 'JavaFrame');
jFig.setMaximized(true);
figure(gcf);
drawnow;

strFig = ['AggregateGraphDim' num2str(intMaxDim) '_' strFigType '_' getDate];
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

