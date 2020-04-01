
%start
clear all;
close all;

strDataDir = 'D:\Data\Results\Data4\';
sFiles = dir([strDataDir '*Shuffnone*.mat']);
intFiles = numel(sFiles);

%pre-allocate data
intMaxDim = 10;
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
	
	%load shuffled data
	cellSplit = strsplit(strFile,'_');
	sShuffFile = dir([strDataDir '*' cellSplit{2} '*' cellSplit{3} '*Shuffreps*.mat']);
	if isempty(sShuffFile)
		vecRemDataSets(intFileIdx) = true;
		continue;
	end
	sLoadS = load([strDataDir sShuffFile(1).name]);
	
	%get relative information for raw and shuffled data
	matRelRaw = sLoad.matPredI ./ sLoad.matFullI;
	matRelShuff = sLoadS.matPredI ./ sLoadS.matFullI;
	if size(matRelRaw,1) == size(matRelShuff,1)
		%normalize to shuffled
		matPhi = (matRelRaw - matRelShuff);% ./ (1 - matRelShuff);
		
	else
		disp a
		%normalize to shuffled
		vecS = (1 - mean(matRelShuff,1));
		matR = (1 - matRelRaw);
		matPhi = (bsxfun(@minus,vecS,matR));% ./ matR;
	end
	
	%check if min dim
	if size(matPhi,2) < intMaxDim,
		vecRemDataSets(intFileIdx) = true;
		continue;
	end
	
	%process
	%matPhi = (sLoad.matPredI  - sLoad.matRandI) ./ (sLoad.matFullI - sLoad.matRandI);
	cellPhiAgg{intFileIdx} = sum(matPhi(:,1:intMaxDim),2)/intMaxDim;
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

%plot
intN = numel(cellPhiAgg);
figure
subplot(2,3,[1 2 4 5])
hold on
plot([0.5 intN+0.5],[0 0],'k--');
errorbar(1:intN,repmat(mean(vecBaselineData),[1 intN]),repmat(std(vecBaselineData),[1 intN]));
title(['\phi ' sprintf('(AUC) up to subspace rank %d',intMaxDim)]);
legend({'DiCo=0','DiCo=behav'});
fixfig;
matAggPhi = cell2mat(cellPhiAgg);
boxplot(matAggPhi);
set(gca,'xticklabel',cellExpName,'XTickLabelRotation',45)
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
