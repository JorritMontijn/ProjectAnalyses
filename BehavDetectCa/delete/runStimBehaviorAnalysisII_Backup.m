%function vecPredictorAssemblyZscores = doPlotActivityDistributions(cellContrastDataAD,sMetaDataAD,sPC,sParams)

%if ~exist('sParams','var'),sParams = struct;end
%if isfield(sParams,'intStartSes'),intStartSes=sParams.intStartSes;else intStartSes=1;end
%if isfield(sParams,'intStopSes'),intStopSes=sParams.intStopSes;else intStopSes=length(cellContrastDataAD);end
close all
clear all
%get block data
intMouse = 5;
if ~exist('cellContrastDataAD','var')
	if intMouse == 1
		strSes = '20140207';
	elseif intMouse == 2
		strSes = '20140314';
	elseif intMouse == 3
		strSes = '20140425';
	elseif intMouse == 4
		%return; %exclude?; bad behavior, weird signals
		strSes = '20140430';
	elseif intMouse == 5
		strSes = '20140507';
	elseif intMouse == 6
		strSes = '20140530';
	elseif intMouse == 7
		strSes = '20140604';
	end
	load(['D:\Data\Results\stimdetection\dataRawPre_aggregate' strSes '.mat']);
end

%vecBlock = [1 1 1 2 2 2 2 2]; %define whether neurons are in same population or not
vecBlockTypes = unique(vecBlock);
intNumBlocks = length(vecBlockTypes);
%vecNeuronNum = zeros(1,intNumBlocks);
%cellKeepList = cell(1,intNumBlocks);
%#ok<*ASGLU>
%#ok<*AGROW>
clear cellSes sLoad sSesAggregate ses
sParams.strFigDir = ['D:\Data\Results\stimdetection' filesep strSes filesep];
sParams.boolSavePlots = true;
sParams.boolSaveData = true;
strOldDir = cd(sParams.strFigDir);

%define bins
%vecBins = -0.055:0.005:0.105;
intBinStepSize = 0.1;
intBinStart = -1; %-1
intBinStop = 2; %2
vecBins = (intBinStart-intBinStepSize):intBinStepSize:(intBinStop+intBinStepSize);
matPR = zeros(length(vecBins),length(cellContrastDataAD{intStartSes}.cellContrastPR)-1);
matNPR = zeros(length(vecBins),length(cellContrastDataAD{intStartSes}.cellContrastNPR)-1);
matPN = zeros(length(vecBins),length(cellContrastDataAD{intStartSes}.cellContrastPN)-1);
matNPN = zeros(length(vecBins),length(cellContrastDataAD{intStartSes}.cellContrastNPN)-1);
intNumNeurons = 0;


sAgg = struct;
%transform data structure
for intContrast=2:length(cellContrastDataAD{intStartSes}.cellContrastPR)
	%pre-allocate
	vecPR = [];
	vecNPR = [];
	vecPN = [];
	vecNPN = [];
	vecPR_index = [];
	vecNPR_index = [];
	vecPN_index = [];
	vecNPN_index = [];
	
	for intPopulation=vecBlockTypes
		%pre-allocate
		vecPopPR = [];
		vecPopNPR = [];
		vecPopPN = [];
		vecPopNPN = [];
		vecPopPR_index = [];
		vecPopNPR_index = [];
		vecPopPN_index = [];
		vecPopNPN_index = [];
		
		intSubStartSes = vecFirstBlock(intPopulation);
		intSubStopSes = vecLastBlock(intPopulation);
		
		%aggregate
		for intSes=intSubStartSes:intSubStopSes
			vecPopPR = [vecPopPR zscore(cellContrastDataAD{intSes}.cellContrastPR{intContrast})];
			vecPopNPR = [vecPopNPR zscore(cellContrastDataAD{intSes}.cellContrastNPR{intContrast})];
			vecPopPN = [vecPopPN zscore(cellContrastDataAD{intSes}.cellContrastPN{intContrast})];
			vecPopNPN = [vecPopNPN zscore(cellContrastDataAD{intSes}.cellContrastNPN{intContrast})];
			vecPopPR_index = [vecPopPR_index cellContrastDataAD{intSes}.cellContrastPR_index{intContrast}];
			vecPopNPR_index = [vecPopNPR_index cellContrastDataAD{intSes}.cellContrastNPR_index{intContrast}];
			vecPopPN_index = [vecPopPN_index cellContrastDataAD{intSes}.cellContrastPN_index{intContrast}];
			vecPopNPN_index = [vecPopNPN_index cellContrastDataAD{intSes}.cellContrastNPN_index{intContrast}];
		end
		
		%put in aggregate structure
		sAgg(intPopulation).cellContrastPR{intContrast} = vecPopPR;
		sAgg(intPopulation).cellContrastNPR{intContrast} = vecPopNPR;
		sAgg(intPopulation).cellContrastPN{intContrast} = vecPopPN;
		sAgg(intPopulation).cellContrastNPN{intContrast} = vecPopNPN;
		sAgg(intPopulation).cellContrastPR_index{intContrast} = vecPopPR_index;
		sAgg(intPopulation).cellContrastNPR_index{intContrast} = vecPopNPR_index;
		sAgg(intPopulation).cellContrastPN_index{intContrast} = vecPopPN_index;
		sAgg(intPopulation).cellContrastNPN_index{intContrast} = vecPopNPN_index;
		
		%put in multi vec
		vecPR = [vecPR vecPopPR];
		vecNPR = [vecNPR vecPopNPR];
		vecPN = [vecPN vecPopPN];
		vecNPN = [vecNPN vecPopNPN];
		vecPR_index = [vecPR_index vecPopPR_index];
		vecNPR_index = [vecNPR_index vecPopNPR_index];
		vecPN_index = [vecPN_index vecPopPN_index];
		vecNPN_index = [vecNPN_index vecPopNPN_index];
		intNumNeurons = max([intNumNeurons vecPR_index vecNPR_index vecPN_index vecNPN_index]);
	end
	
	%make into bins and add to matrix
	vecCounts = hist(vecPR,vecBins)';
	matPR(:,intContrast-1) = matPR(:,intContrast-1) + vecCounts(1:(end));
	
	vecCounts = hist(vecNPR,vecBins)';
	matNPR(:,intContrast-1) = matNPR(:,intContrast-1) + vecCounts(1:(end));
	
	vecCounts = hist(vecPN,vecBins)';
	matPN(:,intContrast-1) = matPN(:,intContrast-1) + vecCounts(1:(end));
	
	vecCounts = hist(vecNPN,vecBins)';
	matNPN(:,intContrast-1) = matNPN(:,intContrast-1) + vecCounts(1:(end));
end

% =====>
%normalize raw df/f values to range [0 1] per contrast
%or z-score normalize raw df/f values
matPR = conv2(matPR,(1/9)*ones(3,3), 'same');
matNPR = conv2(matNPR,(1/9)*ones(3,3), 'same');
matPN = conv2(matPN,(1/9)*ones(3,3), 'same');
matNPN = conv2(matNPN,(1/9)*ones(3,3), 'same');

boolNorm = true;
if boolNorm
	%normalize counts per contrast
	matNormPR = matPR./repmat(sum(matPR,1),[size(matPR,1) 1]);
	matNormNPR = matNPR./repmat(sum(matNPR,1),[size(matNPR,1) 1]);
	matNormPN = matPN./repmat(sum(matPN,1),[size(matPN,1) 1]);
	matNormNPN = matNPN./repmat(sum(matNPN,1),[size(matNPN,1) 1]);
	matDiffP = (matNormPR) ./ (matNormPR + matNormPN);
	matDiffNP = (matNormNPR) ./ (matNormNPR + matNormNPN);
	%matDiffP(matNormPR==0 & matNormPN == 0) = nan;
	%matDiffNP(matNormNPR==0 & matNormNPN == 0) = nan;
else
	%don't normalize counts per contrast
	matNormPR = matPR;
	matNormNPR = matNPR;
	matNormPN = matPN;
	matNormNPN = matNPN;
	matDiffP = (matNormPR) ./ (matNormPR + matNormPN);
	matDiffNP = (matNormNPR) ./ (matNormNPR + matNormNPN);
	matDiffP(matNormPR==0 & matNormPN == 0) = nan;
	matDiffNP(matNormNPR==0 & matNormNPN == 0) = nan;
end

%transform difference matrices
matDiffP = conv2(matDiffP,(1/3)*ones(3,1), 'same');
matDiffNP = conv2(matDiffNP,(1/3)*ones(3,1), 'same');

%remove upper/lower rows
matNormPR = matNormPR(2:(end-1),:);
matNormNPR = matNormNPR(2:(end-1),:);
matNormPN = matNormPN(2:(end-1),:);
matNormNPN = matNormPR(2:(end-1),:);
matDiffP = matDiffP(2:(end-1),:);
matDiffNP = matDiffNP(2:(end-1),:);

%matNormPR = matNormNPR;
%matNormPN = matNormNPN;
%matDiffP = matDiffNP;

%general vars
vecTickY = [1 11 21 31];%[6 16 26 36];%[1 6 11 16 20];
vecTickLabelY = vecBins(vecTickY+1);
vecTickLabelX = sMetaDataAD{intStartSes}.vecContrasts(2:end)*100;
strLabelY = sprintf('Z-scored dF/F (sigma)');

%plot over contrasts
handleFig1 = figure;
set(handleFig1,'Color',[1 1 1]);
figure(handleFig1);
colormap(jet(256))

subplot(2,2,1)
imagesc(matNormPR)
cellSaveMatrices{1} = matNormPR;
axis xy
title('Pref / Resp; z=normalized count')
set(gca,'YTick',vecTickY,'YTickLabel',vecTickLabelY)
set(gca,'XTick',1:length(vecTickLabelX),'XTickLabel',vecTickLabelX)
ylabel(strLabelY)
xlabel('Contrast (%)')
colorbar

subplot(2,2,2)
imagesc(matNormPN)
cellSaveMatrices{2} = matNormPN;
axis xy
title('Pref / NonResp; z=normalized count')
set(gca,'YTick',vecTickY,'YTickLabel',vecTickLabelY)
set(gca,'XTick',1:length(vecTickLabelX),'XTickLabel',vecTickLabelX)
ylabel(strLabelY)
xlabel('Contrast (%)')
colorbar

subplot(2,2,3)
imagesc(matDiffP,[0.35-eps 0.65+eps])
cellSaveMatrices{3} = matDiffP;
axis xy
title('Diff Pref')
set(gca,'YTick',vecTickY,'YTickLabel',vecTickLabelY)
set(gca,'XTick',1:length(vecTickLabelX),'XTickLabel',vecTickLabelX)
ylabel(strLabelY)
xlabel('Contrast (%)')
colorbar

%select which bins to plot
intTestType = 3;
intBins = length(matDiffP(:,1));
matH = false(size(matDiffP));
matP = nan(size(matDiffP));
if intTestType < 3
	for intContrastIndex=1:size(matDiffP,2)
		if intTestType == 1
			%bonferroni corrected ttest
			alpha=0.05/intBins;
			for intBin=1:length(matDiffP(:,intContrastIndex))
				[h,p]=ttest(matDiffP(:,intContrastIndex),matDiffP(intBin,intContrastIndex),alpha);
				matP(intBin,intContrastIndex) = p;
				matH(intBin,intContrastIndex) = h;
			end
		elseif intTestType == 2
			%z-score bin count percentage, then take  >2sd's
			matP(:,intContrastIndex) = zscore(matDiffP(:,intContrastIndex));
			matH(:,intContrastIndex) = matP(:,intContrastIndex)>1.5;
		end
	end
elseif intTestType == 3
	matP = reshape(zscore(matDiffP(:)),size(matDiffP));
	matH = matP>1.5;
end

%select for calculation all bins that are significant on any contrast
matH(vecBins<1,:) = 0;
indFitBins = false(intBins,1);
for intContrastIndex=1:size(matDiffP,2)
	indFitBins = matH(:,intContrastIndex) == 1 | indFitBins;
end
vecFitBins = find(indFitBins');
intFitBins = numel(vecFitBins);

matPlotP = ones(size(matDiffP));
matPlotP(matH) = matP(matH);
subplot(2,2,4)
imagesc(matH)
cellSaveMatrices{4} = matH;
axis xy
title('Significant outlier bins')
set(gca,'YTick',vecTickY,'YTickLabel',vecTickLabelY)
set(gca,'XTick',1:length(vecTickLabelX),'XTickLabel',vecTickLabelX)
ylabel(strLabelY)
xlabel('Contrast (%)')
colorbar
if sParams.boolSavePlots
	drawnow;
	jFig = get(handle(gcf), 'JavaFrame');
	jFig.setMaximized(true);
	figure(gcf);
	drawnow;
	strFig = sprintf('%sagg_detectcorrelated_heatmaps_raw',strSes);
	export_fig([strFig '.tif']);
	export_fig([strFig '.pdf']);
end

%% plot R^2 for fits of all selected z-score bins
handleFig2 = figure;
set(handleFig2,'Color',[1 1 1]);
figure(handleFig2);

%run analysis on those bins
vecR2 = nan(1,intFitBins);
vecP = nan(1,intFitBins);
vecSlopes = nan(1,intFitBins);
intCounter = 0;
vecPercentageSignal = (sPC.vecCorrect(end) - sPC.vecCorrect(2:end));%./(1 - sPC.vecCorrect(2:end));
for intFitBin=vecFitBins
	intCounter = intCounter + 1;
	%get data
	vecSpecificAssembly = (mean(matDiffP(intFitBin,:),1)-0.5)*100;
	
	%do regression
	sStats=regstats(vecPercentageSignal,vecSpecificAssembly,'linear');
	vecSlopes(intCounter) = sStats.beta(2);
	vecR2(intCounter) = sStats.rsquare;
	vecP(intCounter) = sStats.tstat.pval(2);
end

%plot
subplot(2,2,2)
vecZscore = vecBins(vecFitBins);
vecSlopes(vecP > 0.05) = 1;
scatter(vecZscore,vecR2,[],vecP)
ylabel('Behavioral prediction (R^2) of assembly activation')
xlabel('Z-score activity bin')
xlim([vecBins(2) vecBins(end-1)])
ylim([0 1])

%plot regression at highest R^2
subplot(2,2,1)
[dummy,intMaxIndex] = max(vecR2);
dblCorrBinZ = vecZscore(intMaxIndex);
vecCorrBin = find(vecBins == dblCorrBinZ);
vecPercentageSignal = (sPC.vecCorrect(end) - sPC.vecCorrect(2:end));%./(1 - sPC.vecCorrect(2:end));
vecSpecificAssembly = (mean(matDiffP(vecCorrBin,:),1)-0.5)*100;
scatter(vecPercentageSignal,vecSpecificAssembly)
xlabel('Expected perceptual signal strength from behavior')
ylabel('Stimulus detection specific assembly activation (%)')
%ylim([0 1])

%regression
sStats=regstats(vecSpecificAssembly,vecPercentageSignal,'linear');
vecX = get(gca,'XLim');
vecY = polyval(sStats.beta([2 1]),vecX);
hold on
plot(vecX,vecY,'r')
hold off
title(sprintf(['Correlation of behavioral stimulus detection signal [x] \n'...
	'with presence of highly active neuronal subassembly [y]\nat z-score=%.1f; Linear regression: R^2=%.2f'],...
	dblCorrBinZ,sStats.rsquare))
%% <=================================================================================================================================================================================
%% switch x and y for above plot and perform shuffle of expected behavioral stimulus detection signal over neuronal subassembly activity; save data +  mean of shuffle, then perform ttest as meta analysis to check if this  measure predicts more than shuffled

%% check whether neuron identities are more consistent than expected
%select significant bins
matHeroded = conv2(double(matH),ones(3,3),'same') > 0;
[vecRow,vecCol] = find(matHeroded);
vecZscoreIndices = vecBins(vecRow);
vecContrastIndices = vecCol+1;

%pre-allocate
cellPR = cell(intStopSes-intStartSes+1,length(vecZscoreIndices));
cellPN = cell(intStopSes-intStartSes+1,length(vecZscoreIndices));
cellDiff = cell(intStopSes-intStartSes+1,length(vecZscoreIndices));
matPR_index = zeros(intNumNeurons,intStopSes-intStartSes+1);
matPN_index = zeros(intNumNeurons,intStopSes-intStartSes+1);
matDiff_index = zeros(intNumNeurons,intStopSes-intStartSes+1);

%select neuron identities belonging to this bin per session
intSesCounter = 0;
for intSes=intStartSes:intStopSes
	intSesCounter = intSesCounter + 1;
	%loop through all significant bins
	for intBin=1:length(vecZscoreIndices)
		%pre-allocate
		intContrastIndex = vecContrastIndices(intBin);
		vecZscoreRange = [vecZscoreIndices(intBin)-intBinStepSize/2 vecZscoreIndices(intBin)+intBinStepSize/2];
		
		%retrieve data
		vecPR = zscore(cellContrastDataAD{intSes}.cellContrastPR{intContrastIndex});
		vecPN = zscore(cellContrastDataAD{intSes}.cellContrastPN{intContrastIndex});
		vecPR_index = cellContrastDataAD{intSes}.cellContrastPR_index{intContrastIndex};
		vecPN_index = cellContrastDataAD{intSes}.cellContrastPN_index{intContrastIndex};
		
		%remove from response vector all neurons that are also in no-response
		indSelectPR = vecPR > vecZscoreRange(1) & vecPR < vecZscoreRange(2);
		indSelectPN = vecPN > vecZscoreRange(1) & vecPN < vecZscoreRange(2);
		vecIdentitiesPR = vecPR_index(indSelectPR);
		vecIdentitiesPN = vecPN_index(indSelectPN);
		vecIdentitiesDiff = vecIdentitiesPR(~ismember(vecIdentitiesPR,vecIdentitiesPN));
		
		%put in output
		cellPR{intSesCounter,intBin} = vecIdentitiesPR;
		cellPN{intSesCounter,intBin} = vecIdentitiesPN;
		cellDiff{intSesCounter,intBin} = vecIdentitiesDiff;
		matPR_index(vecIdentitiesPR,intSesCounter) = matPR_index(vecIdentitiesPR,intSesCounter) + 1;
		matPN_index(vecIdentitiesPN,intSesCounter) = matPN_index(vecIdentitiesPN,intSesCounter) + 1;
		matDiff_index(vecIdentitiesDiff,intSesCounter) = matDiff_index(vecIdentitiesDiff,intSesCounter) + 1;
	end
end
matDiff_index = logical(matDiff_index);

%plot
subplot(2,2,3)
imagesc(matDiff_index)
ylabel('Neuron ID')
xlabel('Recording block')
colorbar

%correlation
subplot(2,2,4)
matAssemblyCorr = corr(matDiff_index);
matPlotCorr = nan(size(matAssemblyCorr));
for intPopulation=vecBlockTypes
	
	%prepare this block
	vecRecordings = find(vecBlock==intPopulation);
	intRecsInBlock = length(vecRecordings);
	intValues = (intRecsInBlock^2-intRecsInBlock)/2;
	
	%loop through sessions and put correlation values in vector
	for intRecCounter1=1:intRecsInBlock
		intRec1 = vecRecordings(intRecCounter1);
		for intRecCounter2=intRecCounter1:intRecsInBlock
			intRec2 = vecRecordings(intRecCounter2);
			matPlotCorr(intRec2,intRec1) = matAssemblyCorr(intRec2,intRec1);
			matPlotCorr(intRec1,intRec2) = matAssemblyCorr(intRec1,intRec2);
		end
	end
end
imagesc(matPlotCorr)
ylabel('Recording block')
xlabel('Recording block')
nancolorbar([-1 1])
if sParams.boolSavePlots
	drawnow;
	jFig = get(handle(gcf), 'JavaFrame');
	jFig.setMaximized(true);
	figure(gcf);
	drawnow;
	strFig = sprintf('%sagg_detectcorrelated_consistency_raw',strSes);
	export_fig([strFig '.tif']);
	export_fig([strFig '.pdf']);
end

%% plot consistency of response specific subassembly
%loop through blocks and calculate correlation per block
vecCorrBlocks = zeros(1,length(vecBlockTypes));
matCorrConfInterval = zeros(2,length(vecBlockTypes));
cellCorrelationVector = cell(1,length(vecBlockTypes));
cellShuffledCorrMatrix = cell(1,length(vecBlockTypes));
cellShuffledCorrVector = cell(1,length(vecBlockTypes));
intIterations = 1000;
dblAlpha = 0.05;
for intPopulation=vecBlockTypes
	
	%prepare this block
	vecRecordings = find(vecBlock==intPopulation);
	intRecsInBlock = length(vecRecordings);
	intValues = (intRecsInBlock^2-intRecsInBlock)/2;
	cellCorrelationVector{intPopulation} = nan(1,intValues);
	
	%loop through sessions and put correlation values in vector
	intValueCounter = 0;
	for intRecCounter1=1:intRecsInBlock
		intRec1 = vecRecordings(intRecCounter1);
		for intRecCounter2=(intRecCounter1+1):intRecsInBlock
			intValueCounter = intValueCounter + 1;
			intRec2 = vecRecordings(intRecCounter2);
			cellCorrelationVector{intPopulation}(intValueCounter) = matAssemblyCorr(intRec1,intRec2);
		end
	end
	vecCorrBlocks(intPopulation) = mean(cellCorrelationVector{intPopulation});
	
	%shuffle correlations over recordings per block
	cellShuffledCorrMatrix{intPopulation} = nan(intIterations,intValues);
	vecNeuronsPerRecording = sum(matDiff_index(:,vecRecordings),1);
	intTotalAssemblyNeurons = sum(vecNeuronsPerRecording);
	intAllNeurons = vecNeuronNum(intPopulation);
	for intIter=1:intIterations
		%create random vectors
		matShuffled = zeros(intAllNeurons,intRecsInBlock);
		vecNeuronsPerRecShuffled = hist(ceil(rand(1,intTotalAssemblyNeurons)*intRecsInBlock),intRecsInBlock);
		for intRec=1:intRecsInBlock
			vecRandNeurons = randperm(2*intAllNeurons,vecNeuronsPerRecShuffled(intRec));
			matShuffled(vecRandNeurons<intAllNeurons,intRec) = 1;
			if any(vecRandNeurons>intAllNeurons)
				vecDouble = vecRandNeurons(vecRandNeurons>intAllNeurons) - intAllNeurons;
				matShuffled(vecDouble,intRec) = matShuffled(vecDouble,intRec) + 1;
			end
		end
		
		%calculate correlations
		matShuffledAssemblyCorr = corr(matShuffled);
		
		%loop through sessions and put correlation values in vector
		intValueCounter = 0;
		for intRec1=1:intRecsInBlock
			for intRec2=(intRec1+1):intRecsInBlock
				intValueCounter = intValueCounter + 1;
				
				%put correlation into matrix
				cellShuffledCorrMatrix{intPopulation}(intIter,intValueCounter) = matShuffledAssemblyCorr(intRec1,intRec2);
			end
		end
	end
	
	%calculate if shuffled correlation is different from actual
	%correlation
	cellShuffledCorrVector{intPopulation} = mean(cellShuffledCorrMatrix{intPopulation},2);
	vecShuffSorted = sort(cellShuffledCorrVector{intPopulation},'ascend');
	
	%get confidence intervals
	intLow = round(intIterations*dblAlpha);
	intHigh = round(intIterations*(1-dblAlpha));
	
	%put correlation into matrix
	dblLowVal = vecShuffSorted(intLow);
	dblHighVal =  vecShuffSorted(intHigh);
	matCorrConfInterval(:,intPopulation) = [dblLowVal; dblHighVal];
	
	
end

%plot
handleFig3 = figure;
set(handleFig3,'Color',[1 1 1]);
figure(handleFig3);

%plot confidence intervals of correlation values
%subplot(2,2,1)
intCounter = 0;
vecX = nan(1,length(vecBlockTypes));
vecY = nan(1,length(vecBlockTypes));
vecL = nan(1,length(vecBlockTypes));
vecH = nan(1,length(vecBlockTypes));
cellLabelsX = cell(1,length(vecBlockTypes));
for intPopulation = vecBlockTypes
	intCounter = intCounter + 1;
	dblMeanAct = mean(cellShuffledCorrVector{intPopulation});
	vecShuffSorted = sort(cellShuffledCorrVector{intPopulation},'ascend');
	dblLowVal = vecShuffSorted(intLow);
	dblHighVal =  vecShuffSorted(intHigh);
	
	vecX(intCounter) = intPopulation;
	vecY(intCounter) = dblMeanAct;
	vecL(intCounter) = dblMeanAct-dblLowVal;
	vecH(intCounter) = dblHighVal-dblMeanAct;
	cellLabelsX{intPopulation} = [num2str(vecFirstBlock(intPopulation)) ' - ' num2str(vecLastBlock(intPopulation))];
	
	%put in output
	cellSaveHCARNeuronConsistency{intPopulation} = [vecCorrBlocks(intPopulation); vecShuffSorted];
end
errorbar(vecX,vecY,vecL,vecH,'LineStyle','none','Marker','o')
hold on
scatter(vecBlockTypes,vecCorrBlocks, 'rx', 'SizeData',60)
plot(get(gca,'XLim'),[0 0],'k--')
hold off
legend('Shuffled','Data','Location','Best')
set(gca,'XTick',vecBlockTypes)
set(gca,'XTickLabel',cellLabelsX)
xlabel('Recording blocks')
ylabel('Pearson correlation')
title('Visual response specific assembly consistency')
ylim([-1 1])
if sParams.boolSavePlots
	drawnow;
	strFig = sprintf('%sagg_detectcorrelated_consistency2_raw',strSes);
	export_fig([strFig '.tif']);
	export_fig([strFig '.pdf']);
end
%end

%% plot properties of response specific subassembly
%create neuron assembly selection indices per block
cellNeuronIndicesPerBlock = cell(1,intNumBlocks);
cellSortingIndexPerBlock = cell(1,intNumBlocks);
cellNeuronIndicesReorderedPerBlock = cell(1,intNumBlocks);
%go through blocks
for intPopulation = vecBlockTypes
	handleFig = figure;
	set(handleFig,'Color',[1 1 1]);
	figure(handleFig);
	
	%get selection vectors&data
	vecSelectBlocks = vecFirstBlock(intPopulation):vecLastBlock(intPopulation);
	cellNeuronIndicesPerBlock{intPopulation} = logical(sum(matDiff_index(1:vecNeuronNum(intPopulation),vecSelectBlocks),2));
	[cellNeuronIndicesReorderedPerBlock{intPopulation},cellSortingIndexPerBlock{intPopulation}]=sort(cellNeuronIndicesPerBlock{intPopulation},'descend');
	sCorr = calcStimCorrs(cellMultiSes{intPopulation});
	intBorderSize = ceil(log10(vecNeuronNum(intPopulation))+1);
	
	%signal correlations
	matReorderedSC = sCorr.matSignalCorrs(cellSortingIndexPerBlock{intPopulation},cellSortingIndexPerBlock{intPopulation});
	matPlotSC = [matReorderedSC; nan(intBorderSize,size(matReorderedSC,intPopulation))];
	matPlotSC = [matPlotSC nan(size(matPlotSC,1),intBorderSize)];
	subplot(2,2,1)
	imagesc(matPlotSC);
	nancolorbar(matPlotSC);
	
	cmap = colormap(jet(2));
	colormap('default');
	
	hold on;
	intCounter = 0;
	for intAssemblyMember = [1 0] %plot pref stim legend bar
		intCounter = intCounter + 1;
		intStartNeuron = find(cellNeuronIndicesReorderedPerBlock{intPopulation}==intAssemblyMember,1,'first');
		intStopNeuron = find(cellNeuronIndicesReorderedPerBlock{intPopulation}==intAssemblyMember,1,'last');
		if ~isempty(intStartNeuron)
			fill([intStartNeuron intStartNeuron intStopNeuron intStopNeuron],[vecNeuronNum(intPopulation)+1 vecNeuronNum(intPopulation)+intBorderSize vecNeuronNum(intPopulation)+intBorderSize vecNeuronNum(intPopulation)+1],cmap(intCounter,:),'EdgeColor',cmap(intCounter,:));
			fill([vecNeuronNum(intPopulation)+1 vecNeuronNum(intPopulation)+intBorderSize vecNeuronNum(intPopulation)+intBorderSize vecNeuronNum(intPopulation)+1],[intStartNeuron intStartNeuron intStopNeuron intStopNeuron],cmap(intCounter,:),'EdgeColor',cmap(intCounter,:));
		end
	end
	hold off;
	title(sprintf('Signal correlation for block %d; blue=response correlated',intPopulation));
	
	xlabel('Neuron # [grouped by assembly membership]');
	ylabel('Neuron # [grouped by assembly membership]');
	
	%plot mean signal correlation within assembly, between assembly/non-assembly and within non-assembly
	subplot(2,2,2)
	[matMeanSC,matStdSC,cellValsSC] = getMatrixBlockMeans(matReorderedSC,cellNeuronIndicesReorderedPerBlock{intPopulation});
	vecY = [matMeanSC(1,1) matMeanSC(2,1) matMeanSC(2,2)];%within assembly/between/within non-assembly
	vecE = [matStdSC(1,1)/sqrt(length(cellValsSC{1,1})) matStdSC(2,1)/sqrt(length(cellValsSC{2,1})) matStdSC(2,2)/sqrt(length(cellValsSC{2,2}))];
	[h,dblP_AN] = ttest2(cellValsSC{1,1},cellValsSC{2,2}); %assembly/non-assembly
	[h,dblP_AB] = ttest2(cellValsSC{1,1},cellValsSC{2,1});%assembly/between
	[h,dblP_BN] = ttest2(cellValsSC{2,1},cellValsSC{2,2});%between/non-assembly
	errorbar(vecY,vecE,'ob','LineStyle','none');
	set(gca,'XTick',[1 2 3],'XTickLabel',{'Assembly','Between','Non-assembly'})
	ylabel('Mean Signal Correlation')
	title(sprintf('Mean +/- st err of signal correlation; A-N, p=%.3f; A-B, p=%.3f; B-N, p=%.3f',dblP_AN,dblP_AB,dblP_BN))
	
	%put in output
	cellSaveSignalCorrs{1,intPopulation} = cellValsSC{1,1};%signal correlations [3 (within-HCAR, between, within-non-HCAR) x n (animals)] with every cell = vector of pairwise correlation values
	cellSaveSignalCorrs{2,intPopulation} = cellValsSC{2,1};
	cellSaveSignalCorrs{3,intPopulation} = cellValsSC{2,2};
	
	%noise correlations
	subplot(2,2,3)
	matReorderedNC = sCorr.matNoiseCorrs(cellSortingIndexPerBlock{intPopulation},cellSortingIndexPerBlock{intPopulation});
	matPlotNC = [matReorderedNC; nan(intBorderSize,size(matReorderedNC,intPopulation))];
	matPlotNC = [matPlotNC nan(size(matPlotNC,1),intBorderSize)];
	imagesc(matPlotNC);
	nancolorbar(matPlotNC);
	
	cmap = colormap(jet(2));
	colormap('default');
	
	hold on;
	intCounter = 0;
	for intAssemblyMember = [1 0] %plot pref stim legend bar
		intCounter = intCounter + 1;
		intStartNeuron = find(cellNeuronIndicesReorderedPerBlock{intPopulation}==intAssemblyMember,1,'first');
		intStopNeuron = find(cellNeuronIndicesReorderedPerBlock{intPopulation}==intAssemblyMember,1,'last');
		if ~isempty(intStartNeuron)
			fill([intStartNeuron intStartNeuron intStopNeuron intStopNeuron],[vecNeuronNum(intPopulation)+1 vecNeuronNum(intPopulation)+intBorderSize vecNeuronNum(intPopulation)+intBorderSize vecNeuronNum(intPopulation)+1],cmap(intCounter,:),'EdgeColor',cmap(intCounter,:));
			fill([vecNeuronNum(intPopulation)+1 vecNeuronNum(intPopulation)+intBorderSize vecNeuronNum(intPopulation)+intBorderSize vecNeuronNum(intPopulation)+1],[intStartNeuron intStartNeuron intStopNeuron intStopNeuron],cmap(intCounter,:),'EdgeColor',cmap(intCounter,:));
		end
	end
	hold off;
	title(sprintf('Noise correlation for block %d; blue=response correlated',intPopulation));
	
	xlabel('Neuron # [grouped by assembly membership]');
	ylabel('Neuron # [grouped by assembly membership]');
	
	%plot mean noise correlation within assembly, between assembly/non-assembly and within non-assembly
	subplot(2,2,4)
	[matMeanNC,matStdNC,cellValsNC] = getMatrixBlockMeans(matReorderedNC,cellNeuronIndicesReorderedPerBlock{intPopulation});
	vecY = [matMeanNC(1,1) matMeanNC(2,1) matMeanNC(2,2)];%within assembly/between/within non-assembly
	vecE = [matStdNC(1,1)/sqrt(length(cellValsNC{1,1})) matStdNC(2,1)/sqrt(length(cellValsNC{2,1})) matStdNC(2,2)/sqrt(length(cellValsNC{2,2}))];
	[h,dblP_AN] = ttest2(cellValsNC{1,1},cellValsNC{2,2}); %assembly/non-assembly
	[h,dblP_AB] = ttest2(cellValsNC{1,1},cellValsNC{2,1});%assembly/between
	[h,dblP_BN] = ttest2(cellValsNC{2,1},cellValsNC{2,2});%between/non-assembly
	errorbar(vecY,vecE,'ob','LineStyle','none');
	set(gca,'XTick',[1 2 3],'XTickLabel',{'Assembly','Between','Non-assembly'})
	ylabel('Mean Noise Correlation')
	title(sprintf('Mean +/- st err of noise correlation; A-N, p=%.3f; A-B, p=%.3f; B-N, p=%.3f',dblP_AN,dblP_AB,dblP_BN))
	
	%put in output
	cellSaveNoiseCorrs{1,intPopulation} = cellValsNC{1,1};%signal correlations [3 (within-HCAR, between, within-non-HCAR) x n (animals)] with every cell = vector of pairwise correlation values
	cellSaveNoiseCorrs{2,intPopulation} = cellValsNC{2,1};
	cellSaveNoiseCorrs{3,intPopulation} = cellValsNC{2,2};
	
	if sParams.boolSavePlots
		drawnow;
		jFig = get(handle(gcf), 'JavaFrame');
		jFig.setMaximized(true);
		figure(gcf);
		drawnow;
		strFig = sprintf('%sagg_HCARneurons_correlations_pop%d_raw',strSes,intPopulation);
		export_fig([strFig '.tif']);
		export_fig([strFig '.pdf']);
	end
end

%perform dF/F difference comparison
clear cellAggPN;
clear cellAggPR;
clear indSelectHitCorrelatedNeurons;
clear sMeta;
clear sRespDiff;
strFile = sprintf('dataRaw_aggregate%s.mat',strSes);
if exist(['D:\Data\Results\stimdetection\' strFile],'file')
	fprintf('\nLoading aggregate data, please wait...\n')
	load(['D:\Data\Results\stimdetection\' strFile]);
	intNeurons = numel(cellMultiSes{1}.neuron);
else
	fprintf('\nRetrieving aggregate data, please wait...\n')
	for intPopulation = vecBlockTypes
		[sAggSub,sMetaDataAgg] = getStimDetectDataAD(cellMultiSes{intPopulation},cellMultiSes{intPopulation}.structStim);
		sAgg(intPopulation) = sAggSub;
		sMeta(intPopulation) = sMetaDataAgg;
	end
	
	%pre-allocate
	for intContrastIndex = 2:length(sAgg(1).cellContrastPR_index)
		cellAggPR{intContrastIndex} = [];
		cellAggPN{intContrastIndex} = [];
	end
	
	for intPopulation = vecBlockTypes
		intNeurons = numel(cellMultiSes{intPopulation}.neuron);
		matSubAllPR = nan(length(sAgg(intPopulation).cellContrastPR_index),intNeurons);
		matSubAllPN = nan(length(sAgg(intPopulation).cellContrastPN_index),intNeurons);
		for intContrastIndex = 2:length(sAgg(1).cellContrastPR_index)
			vecNeuronsPR = unique(sAgg(intPopulation).cellContrastPR_index{intContrastIndex});
			vecNeuronsPN = unique(sAgg(intPopulation).cellContrastPN_index{intContrastIndex});
			vecNeuronsBoth = vecNeuronsPR(ismember(vecNeuronsPR,vecNeuronsPN));
			
			intNeuronsBoth = length(vecNeuronsBoth);
			
			vecSubPR = nan(1,intNeuronsBoth);
			vecSubPN = nan(1,intNeuronsBoth);
			
			intNeuronCounter = 0;
			for intNeuron=vecNeuronsBoth
				intNeuronCounter = intNeuronCounter + 1;
				indSelectPR = sAgg(intPopulation).cellContrastPR_index{intContrastIndex} == intNeuron;
				indSelectPN = sAgg(intPopulation).cellContrastPN_index{intContrastIndex} == intNeuron;
				
				vecTempPR = sAgg(intPopulation).cellContrastPR{intContrastIndex}(indSelectPR);
				vecTempPN = sAgg(intPopulation).cellContrastPN{intContrastIndex}(indSelectPN);
				
				vecSubPR(intNeuronCounter) = mean(vecTempPR);
				vecSubPN(intNeuronCounter) = mean(vecTempPN);
				matSubAllPR(intContrastIndex,intNeuron) = mean(vecTempPR);
				matSubAllPN(intContrastIndex,intNeuron) = mean(vecTempPN);
			end
			cellAggPR{intContrastIndex} = [cellAggPR{intContrastIndex} vecSubPR];
			cellAggPN{intContrastIndex} = [cellAggPN{intContrastIndex} vecSubPN];
			
			
			sRespDiff(intPopulation).matSubAllPR = matSubAllPR;
			sRespDiff(intPopulation).matSubAllPN = matSubAllPN;
		end
	end
end

%plot
figure
intCounter = 0;
intMin = -5;
intMax = 5;
intStep = 2/3;
vecX = (intMin+intStep/2):intStep:(intMax-intStep/2);
matNormN = zeros(5,length(vecX));
matResiduals = zeros(5,length(vecX));
vecAggDiffRespNormP = [];
matAggDiffRespNormP = zeros(5,intNeurons);
for intContrastIndex = 2:length(cellAggPR)
	intCounter = intCounter + 1;
	
	subplot(2,5,intCounter)
	vecDiffRespP = cellAggPR{intContrastIndex} - cellAggPN{intContrastIndex};
	dblMu = mean(vecDiffRespP);
	dblSigma = std(vecDiffRespP);
	[h,p] = ttest(vecDiffRespP);
	
	%normalize by dividing by sigma
	vecDiffRespNormP = vecDiffRespP/dblSigma;
	vecAggDiffRespNormP = [vecAggDiffRespNormP vecDiffRespNormP];
	[vecN]= hist(vecDiffRespNormP,vecX);
	vecNormN = vecN/max(vecN);
	matNormN(intCounter,:) = vecNormN;
	
	%plot difference from normally distributed
	vecGaussY = normpdf(vecX,0,std(vecDiffRespNormP));
	vecGaussY = (vecGaussY/sum(vecGaussY))*sum(vecNormN);
	vecResiduals = vecNormN-vecGaussY;
	matResiduals(intCounter,:) = vecResiduals;
	bar(vecX,vecNormN)
	hold on
	plot(vecX,vecGaussY,'r')
	hold off
	title(sprintf('C=%.3f, Mu=%.3f; mu sd=%.3f',sMeta(1).vecContrasts(intContrastIndex),mean(vecDiffRespP),mean(vecDiffRespP)/std(vecDiffRespP)))
	xlabel('Normalized increase in d(dF/F) for response trials (sd''s)')
	ylabel('Normalized number of neurons (count)')
	
	%put in output
	cellSaveDCAE{intCounter} = vecNormN; %normalized increase in dF/F for detection trials [5 (c) x n (animals)] cell array with each cell = vector of neurons with normalized dF/F increase for detection trials
	
	subplot(2,5,intCounter+5)
	bar(vecX,vecResiduals)
	ylabel('Residuals after Gaussian fit at mu=0 (normalized count)')
	xlabel('Normalized increase in d(dF/F) for response trials (sd''s)')
	ylim([-0.3 0.5])
end
if sParams.boolSavePlots
	drawnow;
	jFig = get(handle(gcf), 'JavaFrame');
	jFig.setMaximized(true);
	figure(gcf);
	drawnow;
	strFig = sprintf('%sagg_detectcorrelated_actEnhancementAll_raw',strSes);
	export_fig([strFig '.tif']);
	export_fig([strFig '.pdf']);
end

figure
subplot(1,2,1)

%compute measures on aggregate over contrasts
[vecN]= hist(vecAggDiffRespNormP,vecX);
vecNormN = vecN/max(vecN);
vecGaussY = normpdf(vecX,0,std(vecAggDiffRespNormP));
vecGaussY = (vecGaussY/sum(vecGaussY))*sum(vecNormN);
vecResiduals = vecNormN-vecGaussY;
[h,p] = ttest(vecAggDiffRespNormP);

%plot
bar(vecX,vecNormN)
hold on
plot(vecX,vecGaussY,'r')
hold off


title(sprintf('Mean over all contrasts; mu=%.3f; p=%.4f',mean(vecAggDiffRespNormP),p))
xlabel('Normalized increase in d(dF/F) for response trials (sd''s)')
ylabel('Normalized number of neurons (count)')


subplot(1,2,2)
bar(vecX,vecResiduals)
ylabel('Residuals after Gaussian fit at mu=0 (normalized count)')
xlabel('Normalized increase in d(dF/F) for response trials (sd''s)')
ylim([-0.3 0.5])

if sParams.boolSavePlots
	drawnow;
	strFig = sprintf('%sagg_detectcorrelated_actEnhancementMean_raw',strSes);
	export_fig([strFig '.tif']);
	export_fig([strFig '.pdf']);
end
cellSaveSignalCorrsBiDir = cell(3,length(vecBlockTypes));
cellSaveNoiseCorrsBiDir = cell(3,length(vecBlockTypes));

for intPopulation = vecBlockTypes	
	%% calculation based on HCAR neurons [hit-correlated activity region]
	%calculate neuronal stimulus-detection correlated firing enhancement
	vecPrefStim = mod(sMeta(intPopulation).vecPrefStim,4);
	vecPrefStim(vecPrefStim == 0) = 4;
	intNeurons = length(vecPrefStim);
	intTrials = length(cellMultiSes{intPopulation}.structStim.FrameOn);
	matRespDiff = sRespDiff(intPopulation).matSubAllPR-sRespDiff(intPopulation).matSubAllPN;
	matRespDiffNorm = matRespDiff./repmat(nanstd(matRespDiff,[],2),[1 intNeurons]);
	vecStimDetectActInc = nanmean(matRespDiffNorm(2:3,:),1); %contrast 2&3 (0.005 & 0.02 because significant)
	
	%get neuronal responses per trial
	[matTrialResponse,cellSelectContrasts] = getTrialResponseData(cellMultiSes{intPopulation},cellMultiSes{intPopulation}.structStim);
	
	%get orientation-based trial selection vectors
	sTypesOri = getStimulusTypes(cellMultiSes{intPopulation},{'Orientation'});
	cellSelectOri = getSelectionVectors(cellMultiSes{intPopulation}.structStim,sTypesOri);
	
	%group per pref stim
	vecReorder = [];
	vecReorderedPrefStim = [];
	for intStim=unique(vecPrefStim)
		%neurons
		vecSubPop = find(vecPrefStim==intStim);
		[vecVals,vecSubSort] = sort(vecStimDetectActInc(vecSubPop),'descend');
		vecSortOrig = vecSubPop(vecSubSort);
		vecReorder = [vecReorder vecSortOrig];
		vecReorderedPrefStim = [vecReorderedPrefStim ones(1,length(vecSubPop))*intStim];
	end
	
	%loop through contrasts & stim types for trial reordering
	vecTrialOriVector = [];
	vecTrialContrastVector = [];
	
	for intStim=unique(vecPrefStim)
		for intContrastIndex=1:length(cellSelectContrasts)
			%trials
			indSelectOri = cellSelectOri{intStim} | cellSelectOri{intStim+4};
			vecTrialOriVector = [vecTrialOriVector find(indSelectOri)];
			
			vecSelectContrastTrials = cellSelectContrasts{intContrastIndex} & indSelectOri;
			vecTrialContrastVector = [vecTrialContrastVector find(vecSelectContrastTrials)];
		end
	end
	
	%normalize per contrast
	matRespNormPerContrast = nan(size(matTrialResponse));
	for intContrastIndex=1:length(cellSelectContrasts)
		vecSelectContrastTrials = cellSelectContrasts{intContrastIndex};
		matRespNormPerContrast(:,vecSelectContrastTrials) = zscore(matTrialResponse(:,vecSelectContrastTrials),[],2);
	end
	
	%sort neurons by detection correlated activity enhancement
	[dummy,vecSortNeurons] = sort(vecStimDetectActInc,'descend');
	
	%plot
	figure
	colormap(hot(256))
	%imagesc(matRespNormPerContrast)
	imagesc(matRespNormPerContrast(vecSortNeurons,vecTrialContrastVector))
	set(gca,'XTick',1:sum(vecSelectContrastTrials):length(vecTrialContrastVector))
	colorbar
	title('Z-scored activation level per neuron per trial')
	ylabel('Neuron number sorted by HCAR')
	xlabel('Trial number, grouped by contrast, sorted by stimulus ori')
	%save plot
	if sParams.boolSavePlots
		drawnow;
		strFig = sprintf('%sagg_detectcorrelated_actdissimilarity_pop%d_raw',strSes,intPopulation);
		export_fig([strFig '.tif']);
		export_fig([strFig '.pdf']);
	end
	
	%pre-allocate
	%HCAR/low HCAR
	indSelectHitCorrelatedNeurons = cellNeuronIndicesPerBlock{intPopulation}';
	indSelectHitUncorrelatedNeurons = ~cellNeuronIndicesPerBlock{intPopulation}';
	
	%high/low DCAE
	%indSelectHitCorrelatedNeurons = vecStimDetectActInc>(mean(vecStimDetectActInc)+std(vecStimDetectActInc));
	%indSelectHitUncorrelatedNeurons = vecStimDetectActInc<(mean(vecStimDetectActInc)-std(vecStimDetectActInc));
	%indSelectHitCorrelatedNeurons = vecStimDetectActInc>(mean(vecStimDetectActInc));
	%indSelectHitUncorrelatedNeurons = vecStimDetectActInc<(mean(vecStimDetectActInc));
	intNeuronsHigh = sum(indSelectHitCorrelatedNeurons);
	intNeuronsLow = sum(indSelectHitUncorrelatedNeurons);
	intValuesHigh = (intNeuronsHigh*intNeuronsHigh-intNeuronsHigh)/2;
	intValuesLow = (intNeuronsLow*intNeuronsLow-intNeuronsLow)/2;
	intValuesTotal = (intNeurons*intNeurons-intNeurons)/2;
	matHeterogeneity = nan(intTrials,intValuesTotal);
	matRespDistNoRemHitCorrelated = nan(intTrials,intValuesHigh);
	matRespDistNoRemHitUncorrelated = nan(intTrials,intValuesLow);
	
	matRespDistLoZRemHitCorrelated = nan(intTrials,intValuesHigh);
	matRespDistLoZRemHitUncorrelated = nan(intTrials,intValuesLow);
	
	matRespDistHiZRemHitCorrelated = nan(intTrials,intValuesHigh);
	matRespDistHiZRemHitUncorrelated = nan(intTrials,intValuesLow);
	
	matRespDistLoARemHitCorrelated = nan(intTrials,intValuesHigh);
	matRespDistLoARemHitUncorrelated = nan(intTrials,intValuesLow);
	
	matRespDistHiARemHitCorrelated = nan(intTrials,intValuesHigh);
	matRespDistHiARemHitUncorrelated = nan(intTrials,intValuesLow);
	
	matTrialIndexHigh = nan(intTrials,intValuesHigh);
	matTrialIndexLow = nan(intTrials,intValuesLow);
	%calculate inverse correlation matrix per trial
	for intTrial=1:intTrials
		%general
		matTrialIndexHigh(intTrial,:) = intTrial;
		matTrialIndexLow(intTrial,:) = intTrial;
		
		%% STANDARD
		%perform calculation for all neurons
		vecActivity = matRespNormPerContrast(:,intTrial);
		matZ1 = repmat(vecActivity,[1 intNeurons]);
		matZ2 = repmat(vecActivity',[intNeurons 1]);
		matDistAll = abs(matZ1 - matZ2);
		
		%do calculation for high response correlated neurons
		vecActivity = matRespNormPerContrast(indSelectHitCorrelatedNeurons,intTrial);
		matZ1 = repmat(vecActivity,[1 intNeuronsHigh]);
		matZ2 = repmat(vecActivity',[intNeuronsHigh 1]);
		matDistHigh = abs(matZ1 - matZ2);
		
		%do calculation for low response correlated neurons
		vecActivity = matRespNormPerContrast(indSelectHitUncorrelatedNeurons,intTrial);
		matZ1 = repmat(vecActivity,[1 intNeuronsLow]);
		matZ2 = repmat(vecActivity',[intNeuronsLow 1]);
		matDistLow = abs(matZ1 - matZ2);
		
		%save data as vectors in matrix
		matSelectAll = tril(true(size(matDistAll)),-1);
		matSelectHigh = tril(true(size(matDistHigh)),-1);
		matSelectLow = tril(true(size(matDistLow)),-1);
		vecRespDistHighHCAR = matDistHigh(matSelectHigh);
		vecRespDistLowHCAR = matDistLow(matSelectLow);
		vecRespDistAll = matDistAll(matSelectAll);
		
		matRespDistNoRemHitCorrelated(intTrial,:) = vecRespDistHighHCAR;
		matRespDistNoRemHitUncorrelated(intTrial,:) = vecRespDistLowHCAR;
		matHeterogeneity(intTrial,:) = vecRespDistAll;
		
		%% REMOVE MOST ACTIVE NEURONS (z-score)
		%do calculation for high response correlated neurons
		vecActivity = matRespNormPerContrast(indSelectHitCorrelatedNeurons,intTrial);
		vecActivity = findmin(vecActivity,round(length(vecActivity)*0.8))'; %remove 20% most active cells
		matZ1 = repmat(vecActivity,[1 length(vecActivity)]);
		matZ2 = repmat(vecActivity',[length(vecActivity) 1]);
		matDistHigh = abs(matZ1 - matZ2);
		
		%do calculation for low response correlated neurons
		vecActivity = matRespNormPerContrast(indSelectHitUncorrelatedNeurons,intTrial);
		vecActivity = findmin(vecActivity,round(length(vecActivity)*0.8))'; %remove 20% most active cells
		matZ1 = repmat(vecActivity,[1 length(vecActivity)]);
		matZ2 = repmat(vecActivity',[length(vecActivity) 1]);
		matDistLow = abs(matZ1 - matZ2);
		
		%save data as vectors in matrix
		matSelectHigh = tril(true(size(matDistHigh)),-1);
		matSelectLow = tril(true(size(matDistLow)),-1);
		vecRespDistHighHCAR = matDistHigh(matSelectHigh);
		vecRespDistLowHCAR = matDistLow(matSelectLow);
		
		matRespDistHiZRemHitCorrelated(intTrial,:) = [vecRespDistHighHCAR;nan(intValuesHigh-length(vecRespDistHighHCAR),1)];
		matRespDistHiZRemHitUncorrelated(intTrial,:) = [vecRespDistLowHCAR;nan(intValuesLow-length(vecRespDistLowHCAR),1)];
		
		%% REMOVE LEAST ACTIVE NEURONS (z-score)
		%do calculation for high response correlated neurons
		vecActivity = matRespNormPerContrast(indSelectHitCorrelatedNeurons,intTrial);
		vecActivity = findmax(vecActivity,round(length(vecActivity)*0.8))'; %remove 20% least active cells
		matZ1 = repmat(vecActivity,[1 length(vecActivity)]);
		matZ2 = repmat(vecActivity',[length(vecActivity) 1]);
		matDistHigh = abs(matZ1 - matZ2);
		
		%do calculation for low response correlated neurons
		vecActivity = matRespNormPerContrast(indSelectHitUncorrelatedNeurons,intTrial);
		vecActivity = findmax(vecActivity,round(length(vecActivity)*0.8))'; %remove 20% least active cells
		matZ1 = repmat(vecActivity,[1 length(vecActivity)]);
		matZ2 = repmat(vecActivity',[length(vecActivity) 1]);
		matDistLow = abs(matZ1 - matZ2);
		
		%save data as vectors in matrix
		matSelectHigh = tril(true(size(matDistHigh)),-1);
		matSelectLow = tril(true(size(matDistLow)),-1);
		vecRespDistHighHCAR = matDistHigh(matSelectHigh);
		vecRespDistLowHCAR = matDistLow(matSelectLow);
		
		matRespDistLoZRemHitCorrelated(intTrial,:) = [vecRespDistHighHCAR;nan(intValuesHigh-length(vecRespDistHighHCAR),1)];
		matRespDistLoZRemHitUncorrelated(intTrial,:) = [vecRespDistLowHCAR;nan(intValuesLow-length(vecRespDistLowHCAR),1)];
		
		%% REMOVE MOST ACTIVE NEURONS (df/f)
		%do calculation for high response correlated neurons
		vecRawActivity = matTrialResponse(indSelectHitCorrelatedNeurons,intTrial);
		vecActivity = matRespNormPerContrast(indSelectHitCorrelatedNeurons,intTrial);
		[dummy,vecIndex] = findmin(vecRawActivity,round(length(vecActivity)*0.8)); %remove 20% most active cells
		vecActivity = vecActivity(vecIndex);
		matZ1 = repmat(vecActivity,[1 length(vecActivity)]);
		matZ2 = repmat(vecActivity',[length(vecActivity) 1]);
		matDistHigh = abs(matZ1 - matZ2);
		
		%do calculation for low response correlated neurons
		vecRawActivity = matTrialResponse(indSelectHitUncorrelatedNeurons,intTrial);
		vecActivity = matRespNormPerContrast(indSelectHitUncorrelatedNeurons,intTrial);
		[dummy,vecIndex] = findmin(vecRawActivity,round(length(vecActivity)*0.8)); %remove 20% most active cells
		vecActivity = vecActivity(vecIndex);
		matZ1 = repmat(vecActivity,[1 length(vecActivity)]);
		matZ2 = repmat(vecActivity',[length(vecActivity) 1]);
		matDistLow = abs(matZ1 - matZ2);
		
		%save data as vectors in matrix
		matSelectHigh = tril(true(size(matDistHigh)),-1);
		matSelectLow = tril(true(size(matDistLow)),-1);
		vecRespDistHighHCAR = matDistHigh(matSelectHigh);
		vecRespDistLowHCAR = matDistLow(matSelectLow);
		
		matRespDistHiARemHitCorrelated(intTrial,:) = [vecRespDistHighHCAR;nan(intValuesHigh-length(vecRespDistHighHCAR),1)];
		matRespDistHiARemHitUncorrelated(intTrial,:) = [vecRespDistLowHCAR;nan(intValuesLow-length(vecRespDistLowHCAR),1)];
		
		%% REMOVE LEAST ACTIVE NEURONS (df/f)
		%do calculation for high response correlated neurons
		vecRawActivity = matTrialResponse(indSelectHitCorrelatedNeurons,intTrial);
		vecActivity = matRespNormPerContrast(indSelectHitCorrelatedNeurons,intTrial);
		[dummy,vecIndex] = findmax(vecRawActivity,round(length(vecActivity)*0.8)); %remove 20% most active cells
		vecActivity = vecActivity(vecIndex);
		matZ1 = repmat(vecActivity,[1 length(vecActivity)]);
		matZ2 = repmat(vecActivity',[length(vecActivity) 1]);
		matDistHigh = abs(matZ1 - matZ2);
		
		%do calculation for low response correlated neurons
		vecRawActivity = matTrialResponse(indSelectHitUncorrelatedNeurons,intTrial);
		vecActivity = matRespNormPerContrast(indSelectHitUncorrelatedNeurons,intTrial);
		[dummy,vecIndex] = findmax(vecRawActivity,round(length(vecActivity)*0.8)); %remove 20% most active cells
		vecActivity = vecActivity(vecIndex);
		matZ1 = repmat(vecActivity,[1 length(vecActivity)]);
		matZ2 = repmat(vecActivity',[length(vecActivity) 1]);
		matDistLow = abs(matZ1 - matZ2);
		
		%save data as vectors in matrix
		matSelectHigh = tril(true(size(matDistHigh)),-1);
		matSelectLow = tril(true(size(matDistLow)),-1);
		vecRespDistHighHCAR = matDistHigh(matSelectHigh);
		vecRespDistLowHCAR = matDistLow(matSelectLow);
		
		matRespDistLoARemHitCorrelated(intTrial,:) = [vecRespDistHighHCAR;nan(intValuesHigh-length(vecRespDistHighHCAR),1)];
		matRespDistLoARemHitUncorrelated(intTrial,:) = [vecRespDistLowHCAR;nan(intValuesLow-length(vecRespDistLowHCAR),1)];
		
	end
	
	for intRemType=1:5
		if intRemType == 1
			%no rem
			matRespDistHitCorrelated = matRespDistNoRemHitCorrelated;
			matRespDistHitUncorrelated = matRespDistNoRemHitUncorrelated;
			strRemType = 'No';
		elseif intRemType == 2
			%rem hi
			matRespDistHitCorrelated = matRespDistHiZRemHitCorrelated;
			matRespDistHitUncorrelated = matRespDistHiZRemHitUncorrelated;
			strRemType = 'High Z';
		elseif intRemType == 3
			%rem lo
			matRespDistHitCorrelated = matRespDistLoZRemHitCorrelated;
			matRespDistHitUncorrelated = matRespDistLoZRemHitUncorrelated;
			strRemType = 'Low Z';
		elseif intRemType == 4
			%rem hi
			matRespDistHitCorrelated = matRespDistHiARemHitCorrelated;
			matRespDistHitUncorrelated = matRespDistHiARemHitUncorrelated;
			strRemType = 'High A';
		elseif intRemType == 5
			%rem lo
			matRespDistHitCorrelated = matRespDistLoARemHitCorrelated;
			matRespDistHitUncorrelated = matRespDistLoARemHitUncorrelated;
			strRemType = 'Low A';
		end
		
		%{
	for intTrial=1:intTrials
		%do calculation for high response correlated neurons
		vecActivity = matRespNormPerContrast(:,intTrial);
		matZ1 = repmat(vecActivity,[1 intNeurons]);
		matZ2 = repmat(vecActivity',[intNeurons 1]);
		matDist = abs(matZ1 - matZ2);
		
		figure
		imagesc(matDist)
		colormap(flipud(jet))
		colorbar
		xlabel('Neuron ID')
		ylabel('Neuron ID')
		title(sprintf('Neuronal pairwise distance in z-scored activation level\nExample trial (# %d)',intTrial))
		pause
		close
	end
		%}
		
		%plot activity dissimilarity over trials
		figure
		subplot(2,2,1)
		%hit correlated neurons
		vecX = 1:intTrials;
		vecInvX = intTrials:-1:1;
		vecMeanH = nanmean(matRespDistHitCorrelated,2)';
		vecErrH = nanstd(matRespDistHitCorrelated,[],2)'/sqrt(intValuesHigh);
		vecMinTrace = vecMeanH-vecErrH;
		vecMaxTrace = vecMeanH+vecErrH;
		
		fill([vecX vecInvX],[vecMaxTrace vecMinTrace(vecInvX)],[1.0 0.7 0.7],'EdgeColor','none');
		hold on
		plot(vecX,vecMeanH,'-','LineWidth',2,'Color',[1 0 0]);
		
		%hit uncorrelated neurons
		vecMeanL = nanmean(matRespDistHitUncorrelated,2)';
		vecErrL = nanstd(matRespDistHitUncorrelated,[],2)'/sqrt(intValuesLow);
		vecMinTrace = vecMeanL-vecErrL;
		vecMaxTrace = vecMeanL+vecErrL;
		
		fill([vecX vecInvX],[vecMinTrace vecMaxTrace(vecInvX)],[0.7 0.7 1],'EdgeColor','none');
		plot(vecX,vecMeanL,'-','LineWidth',2,'Color',[0 0 1]);
		hold off
		[h,dblP] = ttest2(matRespDistHitCorrelated(:),matRespDistHitUncorrelated(:));
		title(sprintf('%s rem; High HCAR (red): mean=%.3f; Low HCAR (blue): mean=%.3f; ttest p=%.3f',strRemType,mean(matRespDistHitCorrelated(:)),mean(matRespDistHitUncorrelated(:)),dblP))
		ylabel(sprintf('Mean population activity dissimilarity'))
		xlabel('Trial number')
		ylim([0 2])
		
		%split data set
		indSelectContrasts = cellMultiSes{intPopulation}.structStim.Contrast == 0.005 | cellMultiSes{intPopulation}.structStim.Contrast == 0.02;
		indSelectResp = cellMultiSes{intPopulation}.structStim.vecTrialResponse;
		
		%hit trials only
		subplot(2,2,3)
		indSelectRespTrials = indSelectContrasts & indSelectResp;
		%hit correlated neurons
		intRespTrials = sum(indSelectRespTrials);
		matRespDistHitcorrelatedSub = matRespDistHitCorrelated(indSelectRespTrials,:);
		vecX = 1:intRespTrials;
		vecInvX = intRespTrials:-1:1;
		vecMeanRespH = nanmean(matRespDistHitcorrelatedSub,2)';
		vecErrRespH = nanstd(matRespDistHitcorrelatedSub,[],2)'/sqrt(intValuesHigh);
		vecMinTraceRespH = vecMeanRespH-vecErrRespH;
		vecMaxTraceRespH = vecMeanRespH+vecErrRespH;
		
		fill([vecX vecInvX],[vecMaxTraceRespH vecMinTraceRespH(vecInvX)],[1.0 0.7 0.7],'EdgeColor','none');
		hold on
		plot(vecX,vecMeanRespH,'-','LineWidth',2,'Color',[1 0 0]);
		
		%hit uncorrelated neurons
		matRespDistHituncorrelatedSub = matRespDistHitUncorrelated(indSelectRespTrials,:);
		vecX = 1:intRespTrials;
		vecInvX = intRespTrials:-1:1;
		vecMeanRespL = nanmean(matRespDistHituncorrelatedSub,2)';
		vecErrRespL = nanstd(matRespDistHituncorrelatedSub,[],2)'/sqrt(intValuesLow);
		vecMinTraceRespL = vecMeanRespL-vecErrRespL;
		vecMaxTraceRespL = vecMeanRespL+vecErrRespL;
		
		fill([vecX vecInvX],[vecMaxTraceRespL vecMinTraceRespL(vecInvX)],[0.7 0.7 1],'EdgeColor','none');
		plot(vecX,vecMeanRespL,'-','LineWidth',2,'Color',[0 0 1]);
		hold off
		[h,dblP] = ttest2(matRespDistHituncorrelatedSub(:),matRespDistHitcorrelatedSub(:));
		dblMeanRH = nanmean(matRespDistHitcorrelatedSub(:));
		dblSdRH = nanstd(matRespDistHitcorrelatedSub(:));
		intValsRH = length(matRespDistHitcorrelatedSub(:));
		dblMeanRL = nanmean(matRespDistHituncorrelatedSub(:));
		dblSdRL = nanstd(matRespDistHituncorrelatedSub(:));
		intValsRL = length(matRespDistHituncorrelatedSub(:));
		title(sprintf('%s rem; Resp; High HCAR (red): mean=%.3f; Low HCAR (blue): mean=%.3f; ttest p=%.3f',strRemType,dblMeanRH,dblMeanRL,dblP))
		ylabel(sprintf('Mean population activity dissimilarity'))
		xlabel('Trial number')
		ylim([0 2])
		
		%miss trials only
		indSelectNoRespTrials = indSelectContrasts & ~indSelectResp;
		subplot(2,2,4)
		%hit correlated neurons
		intNoRespTrials = sum(indSelectNoRespTrials);
		matNoRespHigh = matRespDistHitCorrelated(indSelectNoRespTrials,:);
		vecX = 1:intNoRespTrials;
		vecInvX = intNoRespTrials:-1:1;
		vecMeanNoRespH = nanmean(matNoRespHigh,2)';
		vecErrNoRespH = nanstd(matNoRespHigh,[],2)'/sqrt(intValuesHigh);
		vecMinTraceNoRespH = vecMeanNoRespH-vecErrNoRespH;
		vecMaxTraceNoRespH = vecMeanNoRespH+vecErrNoRespH;
		
		fill([vecX vecInvX],[vecMaxTraceNoRespH vecMinTraceNoRespH(vecInvX)],[1.0 0.7 0.7],'EdgeColor','none');
		hold on
		plot(vecX,vecMeanNoRespH,'-','LineWidth',2,'Color',[1 0 0]);
		
		%hit uncorrelated neurons
		matNoRespLow = matRespDistHitUncorrelated(indSelectNoRespTrials,:);
		vecX = 1:intNoRespTrials;
		vecInvX = intNoRespTrials:-1:1;
		vecMeanNoRespL = nanmean(matNoRespLow,2)';
		vecErrNoRespL = nanstd(matNoRespLow,[],2)'/sqrt(intValuesLow);
		vecMinTraceNoRespL = vecMeanNoRespL-vecErrNoRespL;
		vecMaxTraceNoRespL = vecMeanNoRespL+vecErrNoRespL;
		
		fill([vecX vecInvX],[vecMaxTraceNoRespL vecMinTraceNoRespL(vecInvX)],[0.7 0.7 1],'EdgeColor','none');
		plot(vecX,vecMeanNoRespL,'-','LineWidth',2,'Color',[0 0 1]);
		hold off
		[h,dblP] = ttest2(matNoRespLow(:),matNoRespHigh(:));
		dblMeanNRH = nanmean(matNoRespHigh(:));
		dblSdNRH = nanstd(matNoRespHigh(:));
		intValsNRH = length(matNoRespHigh(:));
		dblMeanNRL = nanmean(matNoRespLow(:));
		dblSdNRL = nanstd(matNoRespLow(:));
		intValsNRL = length(matNoRespLow(:));
		
		title(sprintf('%s rem; No Resp; High HCAR (red): mean=%.3f; Low HCAR (blue): mean=%.3f; ttest p=%.3f',strRemType,dblMeanNRH,dblMeanNRL,dblP))
		ylabel(sprintf('Mean population activity dissimilarity'))
		xlabel('Trial number')
		ylim([0 2])
		
		%comparison of hit/miss and hit corr/uncorr
		subplot(2,2,2)
		dblErrNRH = dblSdNRH/sqrt(intValsNRH);
		dblErrNRL = dblSdNRL/sqrt(intValsNRL);
		dblErrRH = dblSdRH/sqrt(intValsRH);
		dblErrRL = dblSdRL/sqrt(intValsRL);
		errorbar(1:4,[dblMeanNRH dblMeanNRL dblMeanRH dblMeanRL],[dblErrNRH dblErrNRL dblErrRH dblErrRL],'Linestyle','none','Marker','x');
		set(gca,'XTick',1:4,'XTickLabel',{'No Resp High','No Resp Low','Resp High','Resp Low'})
		xlim([0.5 4.5])
		ylabel('Mean within-group z-scored activation dissimilarity')
		title(sprintf('%s rem; Block %d',strRemType,intPopulation));
		
		%put in output
		cellSaveNormActDissim{1,1,intPopulation,intRemType} = matRespDistHitcorrelatedSub(:); %within-group z-scored activation dissimilarity [2 (high/low HCAR) x 2 (detect/no-detect)] x n (animals) with every cell = vector of dissimilarity values
		cellSaveNormActDissim{1,2,intPopulation,intRemType} = matNoRespHigh(:);
		cellSaveNormActDissim{2,1,intPopulation,intRemType} = matRespDistHituncorrelatedSub(:);
		cellSaveNormActDissim{2,2,intPopulation,intRemType} = matNoRespLow(:);
		
		if sParams.boolSavePlots
			drawnow;
			jFig = get(handle(gcf), 'JavaFrame');
			jFig.setMaximized(true);
			figure(gcf);
			drawnow;
			strFig = sprintf('%sagg_detectcorrelated_actdissimilarity_pop%d_%srem_raw',strSes,intPopulation,strRemType);
			export_fig([strFig '.tif']);
			export_fig([strFig '.pdf']);
		end
	end
	%no rem
	matRespDistHitCorrelated = matRespDistNoRemHitCorrelated;
	matRespDistHitUncorrelated = matRespDistNoRemHitUncorrelated;
	strRemType = 'No';
	
	%% calculate the face of god
	%group per pref stim
	[vecPrefStimC,vecReorderC] = sort(vecPrefStim(indSelectHitCorrelatedNeurons),'ascend');
	[vecPrefStimU,vecReorderU] = sort(vecPrefStim(indSelectHitUncorrelatedNeurons),'ascend');
	
	%loop through contrasts & stim types for trial reordering
	vecTrialOriVector = [];
	vecTrialContrastVector = [];
	
	intStimType = 0;
	vecStimType = zeros(1,length(cellSelectOri{1}));
	for intStim=unique(vecPrefStim)
		for intContrastIndex=2:length(cellSelectContrasts)
			intStimType = intStimType + 1;
			%trials
			indSelectTrials = (cellSelectOri{intStim} | cellSelectOri{intStim+4}) & cellSelectContrasts{intContrastIndex};
			vecStimType(indSelectTrials) = intStimType;
		end
	end
	
	intCorrType = 2;
	if intCorrType == 1 %all neurons
		%calculate trial-by-trial correlations of population activity
		matCorrC=corr(matRespNormPerContrast(indSelectHitCorrelatedNeurons,:));
		matCorrU=corr(matRespNormPerContrast(indSelectHitUncorrelatedNeurons,:));
	else %only pref-pop
		for intPrefType=1:4
			indPrefNeurons = vecPrefStim == intPrefType;
			indSelectC = indSelectHitCorrelatedNeurons & indPrefNeurons;
			indSelectU = indSelectHitUncorrelatedNeurons & indPrefNeurons;
			
			%calculate trial-by-trial correlations of population activity
			if sum(indSelectC) > 0,cellCorrC{intPrefType}=corr(matRespNormPerContrast(indSelectC,:));else cellCorrC{intPrefType} = nan(intTrials,intTrials);end
			if sum(indSelectU) > 0,cellCorrU{intPrefType}=corr(matRespNormPerContrast(indSelectU,:));else cellCorrU{intPrefType} = nan(intTrials,intTrials);end
		end
		matCorrU = nan(size(cellCorrU{1}));
		matCorrC = nan(size(cellCorrC{1}));
		
		for intStimType=vecStimType
			intPrefType = mod(intStimType,4);
			if intPrefType == 0,intPrefType = 4;end

			%hit
			indSelectHitTrials = vecStimType == intStimType & indSelectResp;

			%get trials
			matCorrU(indSelectHitTrials,indSelectHitTrials) = cellCorrU{intPrefType}(indSelectHitTrials,indSelectHitTrials);
			matCorrC(indSelectHitTrials,indSelectHitTrials) = cellCorrC{intPrefType}(indSelectHitTrials,indSelectHitTrials);
			
			%miss
			indSelectMissTrials = vecStimType == intStimType & ~indSelectResp;

			%get trials
			matCorrU(indSelectMissTrials,indSelectMissTrials) = cellCorrU{intPrefType}(indSelectMissTrials,indSelectMissTrials);
			matCorrC(indSelectMissTrials,indSelectMissTrials) = cellCorrC{intPrefType}(indSelectMissTrials,indSelectMissTrials);
		end
	end
	
	
	%compute absolute trial distance
	matTrialDist = abs(repmat(1:intTrials,[intTrials 1])-repmat((1:intTrials)',[1 intTrials]));
	%matSelect = tril(true(size(matTrialDist)),-1);
	
	%separate across-trial correlation correlations for different trial types (1-4)
	%also split for hit/miss trials
	%select only 0.5 & 2%
	intHitTrials = sum(indSelectResp);
	intMissTrials = sum(~indSelectResp);
	vecDistHit = [];
	vecDistMiss = [];
	vecCorrHitC = [];
	vecCorrHitU = [];
	vecCorrMissC = [];
	vecCorrMissU = [];
	for intStimType=vecStimType
		%% hit; select trials
		indSelectHitTrials = vecStimType == intStimType & indSelectResp;
		
		%get trials
		matDistSubHit = matTrialDist(indSelectHitTrials,indSelectHitTrials);
		matCorrSubHitLow = matCorrU(indSelectHitTrials,indSelectHitTrials);
		matCorrSubHitHigh = matCorrC(indSelectHitTrials,indSelectHitTrials);
		
		%subselection matrix hit
		matSubSelectHit = tril(true(size(matDistSubHit)),-1);
		
		%subselect
		vecDistHit = [vecDistHit matDistSubHit(matSubSelectHit)'];
		vecCorrHitU = [vecCorrHitU matCorrSubHitLow(matSubSelectHit)'];
		vecCorrHitC = [vecCorrHitC matCorrSubHitHigh(matSubSelectHit)'];
		
		%% miss; select trials
		indSelectMissTrials = vecStimType == intStimType & ~indSelectResp;
		
		%get trials
		matDistSubMiss = matTrialDist(indSelectMissTrials,indSelectMissTrials);
		matCorrSubMissLow = matCorrU(indSelectMissTrials,indSelectMissTrials);
		matCorrSubMissHigh = matCorrC(indSelectMissTrials,indSelectMissTrials);
		
		%subselection matrix hit
		matSubSelectMiss = tril(true(size(matDistSubMiss)),-1);
		
		%subselect
		vecDistMiss = [vecDistMiss matDistSubMiss(matSubSelectMiss)'];
		vecCorrMissU = [vecCorrMissU matCorrSubMissLow(matSubSelectMiss)'];
		vecCorrMissC = [vecCorrMissC matCorrSubMissHigh(matSubSelectMiss)'];
	end
	
	%plot
	figure
	subplot(2,2,1)
	imagesc(matCorrU);
	colorbar;
	title(sprintf('Correlation between trials of z-scored activation distance (Low HCAR; block %d)',intPopulation))
	xlabel('Trial')
	ylabel('Trial')
	
	subplot(2,2,2)
	imagesc(matCorrC);
	colorbar;
	title(sprintf('Correlation between trials of z-scored activation distance (High HCAR; block %d)',intPopulation))
	xlabel('Trial')
	ylabel('Trial')
	
	subplot(2,2,3)
	%high HCAR
	intStep = 25;
	intMax = ceil(intTrials/intStep)*intStep;
	vecBinX = 0:intStep:intMax;
	vecPlotX = (intStep/2):intStep:(intMax-intStep/2);
	[nVec,meanVec,stdVec,cellVals,cellIDs] = makeBins(vecDistHit,vecCorrHitC,vecBinX);
	errorbar(vecPlotX,meanVec,stdVec./nVec, 'Color', 'r')
	sStatsH=regstats(vecCorrHitC,vecDistHit,'linear',{'beta','rsquare','tstat'});
	vecX = get(gca,'XLim');
	vecY = polyval(sStatsH.beta([2 1]),vecX);
	hold on
	plot(vecX,vecY,'r')

	%low HCAR
	vecBinX = 0:intStep:intMax;
	vecPlotX = (intStep/2):intStep:(intMax-intStep/2);
	[nVec,meanVec,stdVec,cellVals,cellIDs] = makeBins(vecDistHit,vecCorrHitU,vecBinX);
	errorbar(vecPlotX,meanVec,stdVec./nVec, 'Color', 'b')
	sStatsL=regstats(vecCorrHitU,vecDistHit,'linear',{'beta','rsquare','tstat'});
	vecX = get(gca,'XLim');
	vecY = polyval(sStatsL.beta([2 1]),vecX);
	plot(vecX,vecY,'b')
	hold off
	title(sprintf('Hits; Lin reg high HCAR (red): slope=%.3f; p=%.3f; intercept=%.3f; p=%.3f; R^2=%.3f\nLin reg low HCAR (blue): slope=%.3f; p=%.3f; intercept=%.3f; p=%.3f; R^2=%.3f',...
		sStatsH.beta(2),sStatsH.tstat.pval(2),sStatsH.beta(1),sStatsH.tstat.pval(1),sStatsH.rsquare,sStatsL.beta(2),sStatsL.tstat.pval(2),sStatsL.beta(1),sStatsL.tstat.pval(1),sStatsL.rsquare))
	ylabel(sprintf('Assembly consistency over time\n(Inter-trial correlation of mean population activity)'))
	xlabel('Inter-trial distance')

	subplot(2,2,4)
	%high HCAR
	vecBinX = 0:intStep:intMax;
	vecPlotX = (intStep/2):intStep:(intMax-intStep/2);
	[nVec,meanVec,stdVec,cellVals,cellIDs] = makeBins(vecDistMiss,vecCorrMissC,vecBinX);
	errorbar(vecPlotX,meanVec,stdVec./nVec, 'Color', 'm')
	sStatsH=regstats(vecCorrMissC,vecDistMiss,'linear',{'beta','rsquare','tstat'});
	vecX = get(gca,'XLim');
	vecY = polyval(sStatsH.beta([2 1]),vecX);
	hold on
	plot(vecX,vecY,'m')

	%low HCAR
	vecBinX = 0:intStep:intMax;
	vecPlotX = (intStep/2):intStep:(intMax-intStep/2);
	[nVec,meanVec,stdVec,cellVals,cellIDs] = makeBins(vecDistMiss,vecCorrMissU,vecBinX);
	errorbar(vecPlotX,meanVec,stdVec./nVec, 'Color', 'c')
	sStatsL=regstats(vecCorrMissU,vecDistMiss,'linear',{'beta','rsquare','tstat'});
	vecX = get(gca,'XLim');
	vecY = polyval(sStatsL.beta([2 1]),vecX);
	plot(vecX,vecY,'c')
	hold off
	title(sprintf('Misses; Lin reg high HCAR (magenta): slope=%.3f; p=%.3f; intercept=%.3f; p=%.3f; R^2=%.3f\nLin reg low HCAR (cyan): slope=%.3f; p=%.3f; intercept=%.3f; p=%.3f; R^2=%.3f',...
		sStatsH.beta(2),sStatsH.tstat.pval(2),sStatsH.beta(1),sStatsH.tstat.pval(1),sStatsH.rsquare,sStatsL.beta(2),sStatsL.tstat.pval(2),sStatsL.beta(1),sStatsL.tstat.pval(1),sStatsL.rsquare))
	ylabel(sprintf('Assembly consistency over time\n(Inter-trial correlation of mean population activity)'))
	xlabel('Inter-trial distance')
	
	%%%% OUTPUT
	cellSaveCorrITD{1,1,intPopulation} = vecDistHit;
	cellSaveCorrITD{1,2,intPopulation} = vecCorrHitC;
	cellSaveCorrITD{1,3,intPopulation} = vecCorrHitU;
	cellSaveCorrITD{2,1,intPopulation} = vecDistMiss;
	cellSaveCorrITD{2,2,intPopulation} = vecCorrMissC;
	cellSaveCorrITD{2,3,intPopulation} = vecCorrMissU;
	
	%plot
	if sParams.boolSavePlots
		drawnow;
		jFig = get(handle(gcf), 'JavaFrame');
		jFig.setMaximized(true);
		figure(gcf);
		drawnow;
		strFig = sprintf('%sagg_HCARneurons_dissimilaritycorrelations_pop%d_%srem_raw',strSes,intPopulation,strRemType);
		export_fig([strFig '.tif']);
		export_fig([strFig '.pdf']);
	end
	
	
	%comparison of means
	figure
	dblErrCCHH = std(vecCorrHitC)/sqrt(length(vecCorrHitC));
	dblErrCCHL = std(vecCorrHitU)/sqrt(length(vecCorrHitU));
	dblErrCCMH = std(vecCorrMissC)/sqrt(length(vecCorrMissC));
	dblErrCCML = std(vecCorrMissU)/sqrt(length(vecCorrMissU));
	errorbar(1:4,[mean(vecCorrMissC) mean(vecCorrMissU) mean(vecCorrHitC) mean(vecCorrHitU)],[dblErrCCMH dblErrCCML dblErrCCHH dblErrCCHL],'Linestyle','none','Marker','x');
	set(gca,'XTick',1:4,'XTickLabel',{'Miss HCARN','Miss Other','Hit HCARN','Hit Other'})
	xlim([0.5 4.5])
	ylabel('Mean inter-trial correlation of activation dissimilarity')
	title(sprintf('Block %d',intPopulation));
	
	if sParams.boolSavePlots
		drawnow;
		strFig = sprintf('%sagg_HCARneurons_intertrialcorr_pop%d_%srem_raw',strSes,intPopulation,strRemType);
		export_fig([strFig '.tif']);
		export_fig([strFig '.pdf']);
	end
	
	%plot debiased activity distance increases
	figure
	subplot(1,2,1)
	dblMeanNRHd = dblMeanNRH / dblMeanNRH; %debiased
	dblErrNRHd = dblErrNRH / dblMeanNRH;
	dblMeanRHd = dblMeanRH / dblMeanNRH;
	dblErrRHd = dblErrRH/ dblMeanNRH;
	dblMeanNRLd = dblMeanNRL / dblMeanNRL;
	dblErrNRLd = dblErrNRL / dblMeanNRL;
	dblMeanRLd = dblMeanRL / dblMeanNRL;
	dblErrRLd = dblErrRL/ dblMeanNRL;
	
	errorbar(1:4,[dblMeanNRHd dblMeanNRLd dblMeanRHd dblMeanRLd],[dblErrNRHd dblErrNRLd dblErrRHd dblErrRLd],'Linestyle','none','Marker','x');
	set(gca,'XTick',1:4,'XTickLabel',{'No Resp High','No Resp Low','Resp High','Resp Low'})
	xlim([0.5 4.5])
	ylabel('Mean debiased within-group z-scored activation dissimilarity')
	title(sprintf('%s rem; Block %d',strRemType,intPopulation))
	%relative changes
	
	subplot(1,2,2)
	vecHighResp = matRespDistHitcorrelatedSub(:);
	vecLowResp = matRespDistHituncorrelatedSub(:);
	vecHighRespInv = ((vecHighResp - dblMeanNRH) / dblMeanNRH)*100;
	vecLowRespInv = ((vecLowResp - dblMeanNRL) / dblMeanNRL)*100;
	dblMeanHRI = mean(vecHighRespInv);
	dblErrHRI = std(vecHighRespInv) / sqrt(length(vecHighRespInv));
	dblMeanLRI = mean(vecLowRespInv);
	dblErrLRI = std(vecLowRespInv) / sqrt(length(vecLowRespInv));
	
	errorbar(1:2,[dblMeanHRI dblMeanLRI],[dblErrHRI dblErrLRI],'Linestyle','none','Marker','x');
	set(gca,'XTick',1:2,'XTickLabel',{'High HCAR neurons','Low HCAR neurons'})
	xlim([0.5 2.5])
	ylabel('% increase of within-group z-scored distances during detection trials')
	if sParams.boolSavePlots
		drawnow;
		jFig = get(handle(gcf), 'JavaFrame');
		jFig.setMaximized(true);
		figure(gcf);
		drawnow;
		strFig = sprintf('%sagg_debiased_dissimilarity_pop%d_raw',strSes,intPopulation);
		export_fig([strFig '.tif']);
		export_fig([strFig '.pdf']);
	end
	
	%plot dependency on reaction times
	vecRTs = cellMultiSes{intPopulation}.structStim.vecTrialRespSecs(indSelectRespTrials) - cellMultiSes{intPopulation}.structStim.SecsOn(indSelectRespTrials); %get RTs for hit trials and selected contrasts
	[vecRTsSorted,vecRTSortIndex] = sort(vecRTs,'ascend');
	matHitcorrRTSorted = matRespDistHitcorrelatedSub(vecRTSortIndex,:);
	matHituncorrRTSorted = matRespDistHituncorrelatedSub(vecRTSortIndex,:);
	vecHitcorr = mean(matHitcorrRTSorted,2);
	vecHituncorr = mean(matHituncorrRTSorted,2);
	
	%Z activity
	matActZHCAR = matRespNormPerContrast(indSelectHitCorrelatedNeurons,indSelectRespTrials)';
	matActZNonHCAR = matRespNormPerContrast(indSelectHitUncorrelatedNeurons,indSelectRespTrials)';
	
	matZHCARRTSorted = matActZHCAR(vecRTSortIndex,:);
	matZNonHCARRTSorted = matActZNonHCAR(vecRTSortIndex,:);
	vecZHCAR = mean(matZHCARRTSorted,2);
	vecZNonHCAR = mean(matZNonHCARRTSorted,2);
	
	%dF/F
	matActHCAR = matTrialResponse(indSelectHitCorrelatedNeurons,indSelectRespTrials)';
	matActNonHCAR = matTrialResponse(indSelectHitUncorrelatedNeurons,indSelectRespTrials)';
	
	matAHCARRTSorted = matActHCAR(vecRTSortIndex,:);
	matANonHCARRTSorted = matActNonHCAR(vecRTSortIndex,:);
	vecAHCAR = mean(matAHCARRTSorted,2);
	vecANonHCAR = mean(matANonHCARRTSorted,2);
	
	%save
	cellSaveRTDependency{1,intPopulation} = vecRTsSorted;
	cellSaveRTDependency{2,intPopulation} = vecHitcorr;
	cellSaveRTDependency{3,intPopulation} = vecHituncorr;
	cellSaveRTDependency{4,intPopulation} = vecZHCAR;
	cellSaveRTDependency{5,intPopulation} = vecZNonHCAR;
	cellSaveRTDependency{6,intPopulation} = vecAHCAR;
	cellSaveRTDependency{7,intPopulation} = vecANonHCAR;
	
	%plot
	subplot(2,2,1);
	scatter(vecRTsSorted,vecHitcorr,'r')
	hold on
	scatter(vecRTsSorted,vecHituncorr,'b')
	drawnow
	
	%perform regressions
	sStatsC=regstats(vecHitcorr,vecRTsSorted,'linear',{'beta','rsquare','tstat'});
	vecX = get(gca,'XLim');
	vecY = polyval(sStatsC.beta([2 1]),vecX);
	plot(vecX,vecY,'r')
	
	sStatsU=regstats(vecHituncorr,vecRTsSorted,'linear',{'beta','rsquare','tstat'});
	vecY = polyval(sStatsU.beta([2 1]),vecX);
	
	%plot reaction time vs activity dissimilarity correlation
	plot(vecX,vecY,'b')
	hold off
	
	title(sprintf('Pop heterogeneity; Lin reg high HCAR (blue): slope=%.3f; p=%.3f; intercept=%.3f; p=%.3f; R^2=%.3f\nLin reg low HCAR (red): slope=%.3f; p=%.3f; intercept=%.3f; p=%.3f; R^2=%.3f',...
		sStatsC.beta(2),sStatsC.tstat.pval(2),sStatsC.beta(1),sStatsC.tstat.pval(1),sStatsC.rsquare,sStatsU.beta(2),sStatsU.tstat.pval(2),sStatsU.beta(1),sStatsU.tstat.pval(1),sStatsU.rsquare))
	
	xlabel('Reaction Time (s)')
	ylabel('Mean population activity dissimilarity')
	
	%z-scored activity
	subplot(2,2,2);
	scatter(vecRTsSorted,vecZHCAR,'r')
	hold on
	scatter(vecRTsSorted,vecZNonHCAR,'b')
	drawnow
	
	%perform regressions
	sStatsC=regstats(vecZHCAR,vecRTsSorted,'linear',{'beta','rsquare','tstat'});
	vecX = get(gca,'XLim');
	vecY = polyval(sStatsC.beta([2 1]),vecX);
	plot(vecX,vecY,'r')
	
	sStatsU=regstats(vecZNonHCAR,vecRTsSorted,'linear',{'beta','rsquare','tstat'});
	vecY = polyval(sStatsU.beta([2 1]),vecX);
	
	%plot reaction time vs activity dissimilarity correlation
	plot(vecX,vecY,'b')
	hold off
	
	title(sprintf('Z-scored act; Lin reg high HCAR (blue): slope=%.3f; p=%.3f; intercept=%.3f; p=%.3f; R^2=%.3f\nLin reg low HCAR (red): slope=%.3f; p=%.3f; intercept=%.3f; p=%.3f; R^2=%.3f',...
		sStatsC.beta(2),sStatsC.tstat.pval(2),sStatsC.beta(1),sStatsC.tstat.pval(1),sStatsC.rsquare,sStatsU.beta(2),sStatsU.tstat.pval(2),sStatsU.beta(1),sStatsU.tstat.pval(1),sStatsU.rsquare))
	
	xlabel('Reaction Time (s)')
	ylabel('Mean z-scored population activity')
	
	%dF/F activity
	subplot(2,2,3);
	scatter(vecRTsSorted,vecAHCAR,'r')
	hold on
	scatter(vecRTsSorted,vecANonHCAR,'b')
	drawnow
	
	%perform regressions
	sStatsC=regstats(vecAHCAR,vecRTsSorted,'linear',{'beta','rsquare','tstat'});
	vecX = get(gca,'XLim');
	vecY = polyval(sStatsC.beta([2 1]),vecX);
	plot(vecX,vecY,'r')
	
	sStatsU=regstats(vecANonHCAR,vecRTsSorted,'linear',{'beta','rsquare','tstat'});
	vecY = polyval(sStatsU.beta([2 1]),vecX);
	
	%plot reaction time vs activity dissimilarity correlation
	plot(vecX,vecY,'b')
	hold off
	
	title(sprintf('Z-scored act; Lin reg HCAR (red): slope=%.3f; p=%.3f; intercept=%.3f; p=%.3f; R^2=%.3f\nLin reg non-HCAR (blue): slope=%.3f; p=%.3f; intercept=%.3f; p=%.3f; R^2=%.3f',...
		sStatsC.beta(2),sStatsC.tstat.pval(2),sStatsC.beta(1),sStatsC.tstat.pval(1),sStatsC.rsquare,sStatsU.beta(2),sStatsU.tstat.pval(2),sStatsU.beta(1),sStatsU.tstat.pval(1),sStatsU.rsquare))
	
	xlabel('Reaction Time (s)')
	ylabel('Mean dF/F population activity')
	
	if sParams.boolSavePlots
		drawnow;
		jFig = get(handle(gcf), 'JavaFrame');
		jFig.setMaximized(true);
		figure(gcf);
		drawnow;
		strFig = sprintf('%sagg_behaviorallycorrelated_act_z_dissimilarity_pop%d_raw',strSes,intPopulation);
		export_fig([strFig '.tif']);
		export_fig([strFig '.pdf']);
	end
	
	%% heterogeneity for different contrasts, split for hit/miss (like dF/F)
	%figure
	hHetCon = figure;
	set(hHetCon,'Color',[1 1 1]);
	figure(hHetCon);
	
	%pre-allocate
	vecContrasts = sMetaDataAD{1}.vecContrasts(2:end)*100;
	intContrasts = length(vecContrasts);
	vecP = nan(1,intContrasts);
	vecMetaHitY = nan(1,intContrasts);
	vecMetaHitE = nan(1,intContrasts);
	vecMetaMissY = nan(1,intContrasts);
	vecMetaMissE = nan(1,intContrasts);
	intContrastCounter = 0;
	for intContrastIndex=1:intContrasts
		%get contrast
		dblContrast = vecContrasts(intContrastIndex);
		if dblContrast == 0,continue;end
		intContrastCounter = intContrastCounter + 1;
		
		%get resp
		vecHeteroHit = matHeterogeneity(cellSelectContrasts{intContrastCounter+1} & indSelectResp,:);
		vecHeteroHit = vecHeteroHit(:);
		vecHeteroMiss = matHeterogeneity(cellSelectContrasts{intContrastCounter+1} & ~indSelectResp,:);
		vecHeteroMiss = vecHeteroMiss(:);
		
		%perform analyses per contrast
		[h,p,ci] = ttest2(vecHeteroHit,vecHeteroMiss);
		vecY = [mean(vecHeteroHit) mean(vecHeteroMiss)];
		
		%put in output
		cellSaveMatContHetero{intPopulation}(intContrastCounter,1) = vecY(1);
		cellSaveMatContHetero{intPopulation}(intContrastCounter,2) = vecY(2);
		
		%put in meta vector
		vecP(intContrastIndex) = p;
		vecMetaHitY(intContrastIndex) = vecY(1);
		vecMetaHitE(intContrastIndex) = std(vecHeteroHit)/sqrt(length(vecHeteroHit));
		vecMetaMissY(intContrastIndex) = vecY(2);
		vecMetaMissE(intContrastIndex) = std(vecHeteroMiss)/sqrt(length(vecHeteroMiss));
	end
	
	%pre-compute variables
	vecWindow = [1 5];
	vecWindowSelect = vecWindow(1):vecWindow(end);
	intWL = length(vecWindowSelect);
	vecWindowInv = intWL:-1:1;
	vecWindowPlotInv = vecWindow(end):-1:vecWindow(1);
	vecX = [vecContrasts(vecWindowSelect) vecContrasts(vecWindowPlotInv)];
	vecLineX = vecContrasts(vecWindowSelect);
	
	%get data
	for intResp=[0 1]
		if intResp == 1
			vecMeanTrace = vecMetaHitY(vecWindowSelect);
			vecSE = vecMetaHitE(vecWindowSelect);
			vecColorFill = [0.7 1 0.7];
			vecColorLine = [0 1 0];
		else
			vecColorLine = [1 0 0];
			vecColorFill = [1 0.7 0.7];
			vecMeanTrace = vecMetaMissY(vecWindowSelect);
			vecSE = vecMetaMissE(vecWindowSelect);
		end
		
		vecMinTrace = vecMeanTrace-vecSE;
		vecMaxTrace = vecMeanTrace+vecSE;
		vecY = [vecMinTrace vecMaxTrace(vecWindowInv)];
		
		%plot
		hold on
		fill(vecX,vecY,vecColorFill,'EdgeColor','none');
		plot(vecLineX,vecMeanTrace,'-','LineWidth',2,'Color',vecColorLine);
		hold off
	end
	set(gca,'XScale','log','YScale','linear')
	set(gca,'XTick',vecLineX,'XTickLabel',vecLineX)
	title(sprintf('Pop het [%d]; p-vals: 0.5: %.3f; 2: %.3f; 8: %.3f; 32: %.3f; 100: %.3f',intPopulation,vecP))
	grid on
	xlabel('Contrast')
	ylabel('Population response heterogeneity')
	xlim(vecContrasts(vecWindow))
	%ylim([-0.01 0.06])
	legend({'SEM','Miss','SEM','Hit'},'Location','Best')
	
	drawnow;
	
	if sParams.boolSavePlots
		drawnow;
		strFig = sprintf('%s_heterogeneity_over_contrasts_pop%d_raw',strSes,intPopulation);
		export_fig([strFig '.tif']);
		export_fig([strFig '.pdf']);
	end
	
	
	
	%% dF/F over time + heterogeneity over time; also split for fast/slow trials
	sParamsHet.intWindowLength = 1;
	[vecHeterogeneity,vecActivity_dFoF] = calcSlidingHeteroGen(cellMultiSes{intPopulation},sParamsHet);
	
	%get stim data
	cellFieldsC = {'Contrast'};
	sTypesC = getStimulusTypes(cellMultiSes{intPopulation},cellFieldsC);
	cellSelectC = getSelectionVectors(cellMultiSes{intPopulation}.structStim,sTypesC);
	vecC = sTypesC.matTypes;
	
	%get resp data
	indMiss = isnan(cellMultiSes{intPopulation}.structStim.vecTrialRespSecs);
	indFast = cellMultiSes{intPopulation}.structStim.vecTrialRespSecs < nanmean(cellMultiSes{intPopulation}.structStim.vecTrialRespSecs);
	indSlow = cellMultiSes{intPopulation}.structStim.vecTrialRespSecs >= nanmean(cellMultiSes{intPopulation}.structStim.vecTrialRespSecs);
	
	%general
	vecWindowSecs = [-3 5];
	vecWindow = round(vecWindowSecs*cellMultiSes{intPopulation}.samplingFreq);
	vecLineX = (vecWindow(1):vecWindow(end))/cellMultiSes{intPopulation}.samplingFreq;
	vecWindowInv = length(vecLineX):-1:1;
	vecX = [vecLineX vecLineX(vecWindowInv)];
	
	%plot
	cellSaveActTime = {};
	cellSaveHetTime = {};
	hHetTime = figure;
	hActTime = figure;
	for intC=1:length(vecC)
		indTrials = cellSelectC{intC};
		for intAct=[0 1]
			if intAct == 0,figure(hHetTime);else figure(hActTime);end
			subplot(2,3,intC);
			for intType=1:3
				if intType == 1 %miss
					vecTrials = find(indTrials&indMiss);
					strType = 'Miss';
					vecColorLine = [1 0 0];
					vecColorFill = [1 0.7 0.7];
				elseif intType == 2 %slow
					vecTrials = find(indTrials&indSlow);
					strType = 'Slow';
					vecColorLine = [1 1 0];
					vecColorFill = [1 1 0.7];
				else %fast
					vecTrials = find(indTrials&indFast);
					strType = 'Fast';
					vecColorLine = [0 1 0];
					vecColorFill = [0.7 1 0.7];
				end
				matRawData = zeros(length(vecTrials),length(vecWindow(1):vecWindow(end)));
				
				%get data
				vecStimOn = cellMultiSes{intPopulation}.structStim.FrameOn(vecTrials);
				for intTrial=1:length(vecTrials)
					intStart = vecStimOn(intTrial)+vecWindow(1);
					intStop = vecStimOn(intTrial)+vecWindow(end);
					if intAct == 0,
						matRawData(intTrial,:) = vecHeterogeneity(intStart:intStop);
					else
						matRawData(intTrial,:) = vecActivity_dFoF(intStart:intStop);
					end
					
				end
				cellType{intType} = matRawData;
				
				vecMeanTrace = mean(matRawData,1);
				vecSE = std(matRawData,[],1);
				vecMinTrace = vecMeanTrace-vecSE;
				vecMaxTrace = vecMeanTrace+vecSE;
				vecY = [vecMinTrace vecMaxTrace(vecWindowInv)];
				
				%plot
				hold on
				%fill(vecX,vecY,vecColorFill,'EdgeColor','none');
				plot(vecLineX,vecMeanTrace,'-','LineWidth',2,'Color',vecColorLine);
				hold off
			end
			
			%save data&labels
			if intAct == 0,
				cellSaveHetTime{intC,intPopulation} = cellType;
				title(sprintf('Pop Het [%d]; Contrast %.1f',intPopulation,vecC(intC)*100))
				ylabel('Population response heterogeneity')
			else
				cellSaveActTime{intC,intPopulation} = cellType;
				title(sprintf('Pop Act [%d]; Contrast %.1f',intPopulation,vecC(intC)*100))
				ylabel('Population dF/F')
			end
			grid on
			xlabel('Time after stim onset (s)')
			xlim(vecWindowSecs)
			%ylim([-0.01 0.06])
			drawnow;
		end
	end
	
	if sParams.boolSavePlots
		figure(hHetTime)
		drawnow;
		strFig = sprintf('%s_heterogeneity_over_time_pop%d_raw',strSes,intPopulation);
		export_fig([strFig '.tif']);
		export_fig([strFig '.pdf']);
		
		figure(hActTime)
		drawnow;
		strFig = sprintf('%s_activity_over_time_pop%d_raw',strSes,intPopulation);
		export_fig([strFig '.tif']);
		export_fig([strFig '.pdf']);
	end
	
	%% decode if trial will be fast, slow or miss based on dF/F and heterogeneity
	%{
	per animal, sum miss, slow and fast trials; trake min(m,s,f)-1 trials for
	likelihood construction and take one other random trial per condition to
	decode [make procedure work over all time points]
	=> then bootstrap procedure
	%}
	
	%get data
	if ~exist('vecHeterogeneity','var') || isempty(vecHeterogeneity)
		sParamsHet.intWindowLength = 1;%round(cellMultiSes{1}.samplingFreq);
		[vecHeterogeneity,vecActivity_dFoF] = calcSlidingHeteroGen(cellMultiSes{intPopulation},sParamsHet);
	end
	
	%get resp data
	indMiss = isnan(cellMultiSes{intPopulation}.structStim.vecTrialRespSecs);
	indFast = cellMultiSes{intPopulation}.structStim.vecTrialRespSecs < nanmean(cellMultiSes{intPopulation}.structStim.vecTrialRespSecs);
	indSlow = cellMultiSes{intPopulation}.structStim.vecTrialRespSecs >= nanmean(cellMultiSes{intPopulation}.structStim.vecTrialRespSecs);
	
	%general
	vecWindowSecs = [-3 5];
	vecWindow = round(vecWindowSecs*cellMultiSes{intPopulation}.samplingFreq);
	vecLineX = (vecWindow(1):vecWindow(end))/cellMultiSes{intPopulation}.samplingFreq;
	vecSelectSecs = [-1 0];
	vecSelect = round(vecSelectSecs*cellMultiSes{intPopulation}.samplingFreq);
	intStart = find(vecLineX>=vecSelectSecs(1),1);
	intStop = find(vecLineX>=vecSelectSecs(2),1);
	intFramesPerStim = length(vecWindow(1):vecWindow(end));
	vecSelectFrames = intStart:intStop;
	
	%save vars
	cellSaveActDecodeTime = {};
	cellSaveHetDecodeTime = {};
	hDecRespTypeTime = figure;
	
	%get nr of trials
	intMiss = sum(indMiss);
	vecMiss = find(indMiss);
	intFast = sum(indFast);
	vecFast = find(indFast);
	intSlow = sum(indSlow);
	vecSlow = find(indSlow);
	intNrTrials = min([intMiss intFast intSlow]);
	
	%transform data to easy-access format
	%pre-allocate
	matRawDataAct = zeros(length(indMiss),intFramesPerStim);
	matRawDataHet = zeros(length(indMiss),intFramesPerStim);
	
	%get data
	vecStimOn = cellMultiSes{intPopulation}.structStim.FrameOn;
	for intTrial=1:length(indMiss)
		intStart = vecStimOn(intTrial)+vecWindow(1);
		intStop = vecStimOn(intTrial)+vecWindow(end);
		
		matRawDataHet(intTrial,:) = vecHeterogeneity(intStart:intStop);
		matRawDataAct(intTrial,:) = vecActivity_dFoF(intStart:intStop);
	end
	
	
	%bootstrap decoding
	intIters = 1000;
	
	%pre-allocate output
	matPostAct = nan(3,3,intIters); %[actual stim] x [decoded stim] x [time point] x [trials]; value = probability; order= miss / slow / fast
	matPostHet = nan(3,3,intIters); %[actual stim] x [decoded stim] x [time point] x [trials]; value = probability; order= miss / slow / fast
	
	%loop
	for intIter=1:intIters
		%get random trials
		vecMissRand = randperm(intMiss,intNrTrials);
		vecMissLikelihood = vecMiss(vecMissRand(1:(end-1)));
		intDecodeTrialMiss = vecMiss(vecMissRand(end));
		
		vecFastRand = randperm(intFast,intNrTrials);
		vecFastLikelihood = vecFast(vecFastRand(1:(end-1)));
		intDecodeTrialFast = vecFast(vecFastRand(end));
		
		vecSlowRand = randperm(intSlow,intNrTrials);
		vecSlowLikelihood = vecSlow(vecSlowRand(1:(end-1)));
		intDecodeTrialSlow = vecSlow(vecSlowRand(end));
		
		%get likelihoods df/f
		dblLikeActMissMean = mean(mean(matRawDataAct(vecMissLikelihood,vecSelectFrames),2));
		dblLikeActMissSD = std(mean(matRawDataAct(vecMissLikelihood,vecSelectFrames),2));
		
		dblLikeActFastMean = mean(mean(matRawDataAct(vecFastLikelihood,vecSelectFrames),2));
		dblLikeActFastSD = std(mean(matRawDataAct(vecFastLikelihood,vecSelectFrames),2));
		
		dblLikeActSlowMean = mean(mean(matRawDataAct(vecSlowLikelihood,vecSelectFrames),2));
		dblLikeActSlowSD = std(mean(matRawDataAct(vecSlowLikelihood,vecSelectFrames),2));
		
		%get likelihoods heterogeneity
		dblLikeHetMissMean = mean(mean(matRawDataHet(vecMissLikelihood,vecSelectFrames),2));
		dblLikeHetMissSD = std(mean(matRawDataHet(vecMissLikelihood,vecSelectFrames),2));
		
		dblLikeHetFastMean = mean(mean(matRawDataHet(vecFastLikelihood,vecSelectFrames),2));
		dblLikeHetFastSD = std(mean(matRawDataHet(vecFastLikelihood,vecSelectFrames),2));
		
		dblLikeHetSlowMean = mean(mean(matRawDataHet(vecSlowLikelihood,vecSelectFrames),2));
		dblLikeHetSlowSD = std(mean(matRawDataHet(vecSlowLikelihood,vecSelectFrames),2));
		

		%get act probabilities for miss trial
		dblActMissTrial = mean(matRawDataAct(intDecodeTrialMiss,vecSelectFrames));
		matPostAct(1,1,intIter) = normpdf(dblActMissTrial,dblLikeActMissMean,dblLikeActMissSD);
		matPostAct(1,2,intIter) = normpdf(dblActMissTrial,dblLikeActSlowMean,dblLikeActSlowSD);
		matPostAct(1,3,intIter) = normpdf(dblActMissTrial,dblLikeActFastMean,dblLikeActFastSD);
		
		%get act probabilities for slow trial
		dblActSlowTrial = mean(matRawDataAct(intDecodeTrialSlow,vecSelectFrames));
		matPostAct(2,1,intIter) = normpdf(dblActSlowTrial,dblLikeActMissMean,dblLikeActMissSD);
		matPostAct(2,2,intIter) = normpdf(dblActSlowTrial,dblLikeActSlowMean,dblLikeActSlowSD);
		matPostAct(2,3,intIter) = normpdf(dblActSlowTrial,dblLikeActFastMean,dblLikeActFastSD);
		
		%get act probabilities for fast trial
		dblActFastTrial = mean(matRawDataAct(intDecodeTrialFast,vecSelectFrames));
		matPostAct(3,1,intIter) = normpdf(dblActFastTrial,dblLikeActMissMean,dblLikeActMissSD);
		matPostAct(3,2,intIter) = normpdf(dblActFastTrial,dblLikeActSlowMean,dblLikeActSlowSD);
		matPostAct(3,3,intIter) = normpdf(dblActFastTrial,dblLikeActFastMean,dblLikeActFastSD);
		
		
		%get het probabilities for miss trial
		vecHetMissTrial = mean(matRawDataHet(intDecodeTrialMiss,vecSelectFrames));
		matPostHet(1,1,intIter) = normpdf(vecHetMissTrial,dblLikeHetMissMean,dblLikeHetMissSD);
		matPostHet(1,2,intIter) = normpdf(vecHetMissTrial,dblLikeHetSlowMean,dblLikeHetSlowSD);
		matPostHet(1,3,intIter) = normpdf(vecHetMissTrial,dblLikeHetFastMean,dblLikeHetFastSD);
		
		%get het probabilities for slow trial
		vecHetSlowTrial = mean(matRawDataHet(intDecodeTrialSlow,vecSelectFrames));
		matPostHet(2,1,intIter) = normpdf(vecHetSlowTrial,dblLikeHetMissMean,dblLikeHetMissSD);
		matPostHet(2,2,intIter) = normpdf(vecHetSlowTrial,dblLikeHetSlowMean,dblLikeHetSlowSD);
		matPostHet(2,3,intIter) = normpdf(vecHetSlowTrial,dblLikeHetFastMean,dblLikeHetFastSD);
		
		%get het probabilities for fast trial
		vecHetFastTrial = mean(matRawDataHet(intDecodeTrialFast,vecSelectFrames));
		matPostHet(3,1,intIter) = normpdf(vecHetFastTrial,dblLikeHetMissMean,dblLikeHetMissSD);
		matPostHet(3,2,intIter) = normpdf(vecHetFastTrial,dblLikeHetSlowMean,dblLikeHetSlowSD);
		matPostHet(3,3,intIter) = normpdf(vecHetFastTrial,dblLikeHetFastMean,dblLikeHetFastSD);
	end
	
	%normalize probabilities
	matPostHet(:,1,:) = matPostHet(:,1,:)/mean(reshape(matPostHet(:,1,:),[1 numel(matPostHet(:,1,:))]));
	matPostHet(:,2,:) = matPostHet(:,2,:)/mean(reshape(matPostHet(:,2,:),[1 numel(matPostHet(:,2,:))]));
	matPostHet(:,3,:) = matPostHet(:,3,:)/mean(reshape(matPostHet(:,3,:),[1 numel(matPostHet(:,3,:))]));
	matPostAct(:,1,:) = matPostAct(:,1,:)/mean(reshape(matPostAct(:,1,:),[1 numel(matPostAct(:,1,:))]));
	matPostAct(:,2,:) = matPostAct(:,2,:)/mean(reshape(matPostAct(:,2,:),[1 numel(matPostAct(:,2,:))]));
	matPostAct(:,3,:) = matPostAct(:,3,:)/mean(reshape(matPostAct(:,3,:),[1 numel(matPostAct(:,3,:))]));
	
	%get fraction correct
	%act
	[dummy,matDecodedActMiss]=max(matPostAct(1,:,:),[],2);
	vecMissActCorrect = sum(squeeze(matDecodedActMiss == 1))/intIters;
	[dummy,matDecodedActSlow]=max(matPostAct(2,:,:),[],2);
	vecSlowActCorrect = sum(squeeze(matDecodedActSlow == 2))/intIters;
	[dummy,matDecodedActFast]=max(matPostAct(3,:,:),[],2);
	vecFastActCorrect = sum(squeeze(matDecodedActFast == 3))/intIters;
	
	%het
	[dummy,matDecodedHetMiss]=max(matPostHet(1,:,:),[],2);
	vecMissHetCorrect = sum(squeeze(matDecodedHetMiss == 1))/intIters;
	[dummy,matDecodedHetSlow]=max(matPostHet(2,:,:),[],2);
	vecSlowHetCorrect = sum(squeeze(matDecodedHetSlow == 2))/intIters;
	[dummy,matDecodedHetFast]=max(matPostHet(3,:,:),[],2);
	vecFastHetCorrect = sum(squeeze(matDecodedHetFast == 3))/intIters;
	
	
	%% plot decoding performance
	%{
	make angle for m/s/f separated by 2/3pi. Calc resultant vector angle of
	triplet decoding, graph aspolar plot
	%}
	cellSaveRespDecode{intPopulation,1} = zeros(2,2,3);%vectors; [act/het] x [theta/rho] x [miss/slow/fast]
	cellSaveRespDecode{intPopulation,2} = zeros(2,3);%distance; [act/het] x [miss/slow/fast]
	for intRespType=1:3
		if intRespType == 1
			strTitle = 'Miss';
			matDecodedAct = matDecodedActMiss;
			matDecodedHet = matDecodedHetMiss;
			thetaCorr = (4/3)*pi;
			cellColor = {'r','k','k'};
		elseif intRespType == 2
			strTitle = 'Slow';
			matDecodedAct = matDecodedActSlow;
			matDecodedHet = matDecodedHetSlow;
			thetaCorr = (0/3)*pi;
			cellColor = {'k','y','k'};
		elseif intRespType == 3
			strTitle = 'Fast';
			matDecodedAct = matDecodedActFast;
			matDecodedHet = matDecodedHetFast;
			thetaCorr = (2/3)*pi;
			cellColor = {'k','k','g'};
		end
		subplot(2,2,intRespType)
		%calc corr vector
		[xCorr,yCorr] = pol2cart(thetaCorr,1);
		
		%calc resultant vector act
		[x1,y1] = pol2cart((4/3)*pi,sum(squeeze(matDecodedAct == 1))/intIters);%miss
		[x2,y2] = pol2cart((0/3)*pi,sum(squeeze(matDecodedAct == 2))/intIters);%slow
		[x3,y3] = pol2cart((2/3)*pi,sum(squeeze(matDecodedAct == 3))/intIters);%fast
		x=x1+x2+x3;
		y=y1+y2+y3;
		[thetaAct,rhoAct]=cart2pol(x,y);
		[thetaDiffAct,rhoDiffAct] = cart2pol(xCorr-x,yCorr-y);
		
		%resultant vector het
		[x1,y1] = pol2cart((4/3)*pi,sum(squeeze(matDecodedHet == 1))/intIters);%miss
		[x2,y2] = pol2cart((0/3)*pi,sum(squeeze(matDecodedHet == 2))/intIters);%slow
		[x3,y3] = pol2cart((2/3)*pi,sum(squeeze(matDecodedHet == 3))/intIters);%fast
		x=x1+x2+x3;
		y=y1+y2+y3;
		[thetaHet,rhoHet]=cart2pol(x,y);
		[thetaDiffHet,rhoDiffHet] = cart2pol(xCorr-x,yCorr-y);
		
		%save data
		cellSaveRespDecode{intPopulation,1}(1,1,intRespType) = thetaAct;%vectors;[act/het] x [theta/rho] x [miss/slow/fast]
		cellSaveRespDecode{intPopulation,1}(1,2,intRespType) = rhoAct;
		cellSaveRespDecode{intPopulation,1}(2,1,intRespType) = thetaHet;
		cellSaveRespDecode{intPopulation,1}(2,2,intRespType) = rhoHet;
		cellSaveRespDecode{intPopulation,2}(1,intRespType) = rhoDiffAct;%distance; [act/het] x [miss/slow/fast]
		cellSaveRespDecode{intPopulation,2}(2,intRespType) = rhoDiffHet;%distance; [act/het] x [miss/slow/fast]
		
		%make polar background
		[x,y] = pol2cart((4/3)*pi,0.9);%miss
		text(x,y,'Miss','Color',cellColor{1});
		hold on
		[x,y] = pol2cart((0/3)*pi,0.7);%slow
		text(x,y,'Slow','Color',cellColor{2});
		[x,y] = pol2cart((2/3)*pi,0.9);%fast
		text(x,y,'Fast','Color',cellColor{3});
		[x,y] = pol2cart((1/3)*pi,1);%slow/fast
		plot([0 x],[0 y],'k');
		[x,y] = pol2cart((3/3)*pi,1);%miss/fast
		plot([0 x],[0 y],'k');
		[x,y] = pol2cart((5/3)*pi,1);%slow/miss
		plot([0 x],[0 y],'k');
		ang=0:0.001:2*pi; 
		x=cos(ang);
		y=sin(ang);
		plot(x,y, 'k');

		%plot data
		[x,y]=pol2cart(thetaAct,rhoAct);
		scatter(x,y,'bx');
		[x,y]=pol2cart(thetaHet,rhoHet);
		scatter(x,y,'kx');
		hold off;
		xlim([-1 1]);
		ylim([-1 1]);
		title(sprintf('%s trials',strTitle))
	end
	
	%plot summary
	subplot(2,2,4);
	matPerformance = imnorm(1-abs(cellSaveRespDecode{intPopulation,2}));
	vecActPerfomance = matPerformance(1,:);
	vecHetPerfomance = matPerformance(2,:);
	dblOffset=0.1;
	scatter((1-dblOffset)*ones(size(vecActPerfomance)),vecActPerfomance,'bx')
	hold on
	dblMeanAct=mean(vecActPerfomance);
	dblSDAct = std(vecActPerfomance);
	errorbar((1+dblOffset),dblMeanAct,dblSDAct/sqrt(length(vecActPerfomance)),'xb')
	dblMeanHet=mean(vecHetPerfomance);
	dblSDHet = std(vecHetPerfomance);
	scatter((2-dblOffset)*ones(size(vecHetPerfomance)),vecHetPerfomance,'kx')
	errorbar((2+dblOffset),dblMeanHet,dblSDHet/sqrt(length(vecHetPerfomance)),'xk')
	hold off
	title('Normalized resp type decoding')
	ylabel('Norm decod perf')
	set(gca,'XTick',[1 2],'XTickLabel',{'dF/F0', 'Heterogen'})
	
	%save
	if sParams.boolSavePlots
		drawnow;
		strFig = sprintf('%s_resptypedecoding_over_time_pop%d_raw',strSes,intPopulation);
		export_fig([strFig '.tif']);
		export_fig([strFig '.pdf']);
	end
end

%% save data structures
if sParams.boolSaveData
	save(['D:\Data\Results\stimdetection\' strFile],'sMeta','cellAggPR','cellAggPN','sRespDiff','indSelectHitCorrelatedNeurons','-v7.3');
	vecClock = fix(clock);
	strDate = num2str(vecClock(1:3));
	strRecs = num2str(vecRecordings);
	strFile = ['data_aggregateAD' strSes '_' strrep(strDate,'    ','_') '_' strRecs(~isspace(strRecs))];
	%save(['D:\Data\Results\stimdetection\' strFile],'-v7.3');
	%if exist(['D:\Data\Results\stimdetection\' strFile '.mat'],'file')
	%	load(['D:\Data\Results\stimdetection\' strFile]);
	%end
	save(['D:\Data\Results\stimdetection\' strFile],'cellSave*','-v7.3');
	cd(strOldDir);
end

%{
figures to save and take mean over animals
- cellSavePreStim: dF/F detect/no-detect preceding stimulus frames [5 (c) x 2 (d/nd) cell array] x n (animals) with every cell = vector of values
- cellSaveBehavDetect: behavioral detection [5 (c) x n (animals)] cell array with every cell = vector of detect/no detect values (1/0)
- cellSaveBehavRT: behavioral RT [5 (c) x n (animals)] cell array with every cell = vector of RTs
- cellSaveBaseSubtrStimAct: dF/F detect/no-detect baseline-subtracted activity [5 (c) x 2 (d/nd) cell array] x n (animals) with every cell = vector of values
- cellSaveDuringStim: dF/F detect/no-detect during stimulus frames [5 (c) x 2 (d/nd) cell array] x n (animals) with every cell = vector of values
- cellSaveTraceAct: dF/F trace detect/no-detect pref/non-pref [2 (pref/no-pref) x 2 (detect/no-detect) x 6 (c)] x n (animals) with every cell = [Nreps x Nframes] matrix of trace dF/F values
-  AD: cellSaveMatrices: blob matrices in [1 x 4 (detect/no-detect/diff/significance)] cell array with [5 (c) x 31 (bins)] matrix
-  AD: cellSaveHCARNeuronConsistency: [1 x i (neuronal population)] cell array with vector of shuffled correlation values, where first value is actual data
-  AD: cellSaveSignalCorrs: signal correlations [3 (within-HCAR, between, within-non-HCAR) x n (animals/blocks)] with every cell = vector of pairwise correlation values
-  AD: cellSaveNoiseCorrs: noise correlations [3 (within-HCAR, between, within-non-HCAR) x n (animals/blocks)] with every cell = vector of pairwise correlation values
-  AD: cellSaveHCAR: normalized increase in dF/F for detection trials [5 (c) x n (animals/blocks)] cell array with each cell = vector of neurons with normalized dF/F increase for detection trials
-  AD: cellSaveSignalCorrsBiDir: signal correlations [3 (within-HCAR, between, within-non-HCAR) x n (animals/blocks)] with every cell = vector of pairwise correlation values
-  AD: cellSaveNoiseCorrsBiDir: noise correlations [3 (within-HCAR, between, within-non-HCAR) x n (animals/blocks)] with every cell = vector of pairwise correlation values
-  AD: cellSaveNormActDissim: within-group z-scored activation dissimilarity [2 (high/low HCAR) x 2 (detect/no-detect)] x n (animals/blocks) with every cell = vector of dissimilarity values
-  AD: cellSaveDissimCorrITD: inter-trial-distance dependence of assembly consistency [2 (high/low HCAR) x 2 (hits/misses) x n (animals/blocks)] with every cell = [p x 2] matrix with intertrial distance and correlation values for each trial pair
	- extra graph of above with means comparison of hits/misses high/low HCAR

%}