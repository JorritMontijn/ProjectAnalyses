%% aim
%{
https://www.nature.com/articles/s41598-021-95037-z
https://www.nature.com/articles/s42003-021-02437-y

q1: is coding better in trials with higher or lower firing rates?

q2: can we define peaks in the IFR as population events, and find which cells spike in the beginning
or end? does this ordering differ between orientations?

%}

%% define qualifying areas
%close all;
clear all;
boolSaveFigs = true;
boolHome = true;
if boolHome
	strDataPath = 'F:\Data\Processed\Neuropixels\';
	strFigurePath = 'F:\Data\Results\PopTimeCoding\figures\';
	strTargetDataPath = 'F:\Data\Results\PopTimeCoding\data\';
else
	strDataPath = 'E:\DataPreProcessed\';
	strFigurePath = 'D:\Data\Results\PopTimeCoding\figures\';
	strTargetDataPath = 'D:\Data\Results\PopTimeCoding\data\';
end
dblBimoThreshold = inf;%0.5;%0.4
dblDevThreshold = 0.7;%0.7;%0.7
intPlotMu = 2;
if intPlotMu == 1
	strMu = '\DeltaMu';
else
	strMu = 'Pop mean';
end

%% find data
strStim = 'NM';%DG/NM
%cellTypes = {'Real','Shuff','Poiss','UniStretch','VarFixed','Saturating','TuneFixed'};
cellTypes = {'Real','Shuff','Poiss','UniStretch','VarFixed','Saturating'};
sDir = dir([strTargetDataPath 'TimeCodingAggQC1ABI*.mat']); %or ABA if old
indUseRecs = contains({sDir.name},['ABI_' strStim]);
sDir(~indUseRecs) = [];

cellAggLRActPerQ = cell([0 0 0 0 0]);
cellAggPopMuPerQ = cell([0 0 0 0 0]);
matBC = nan([0 0]);
matMDV = nan([0 0]);
vecCounter = zeros(1,numel(cellTypes));
for intFile=1:numel(sDir)
	%% load data
	strFolder = sDir(intFile).folder;
	strFile = sDir(intFile).name;
	strType = strrep(strrep(getFlankedBy(strFile,'AggQC1ABI','_','first'),'ABI',''),['_' strStim],'');
	intType = find(ismember(cellTypes,strType));
	sData = load(fullpath(strFolder,strFile));
	if size(sData.cellLRActPerQ,2) < 8 || any(flat(cellfun(@(x) any(isnan(x(:))),sData.cellLRActPerQ)))
		continue;
	end
	
	if sData.dblBC > dblBimoThreshold || sData.dblMaxDevFrac > dblDevThreshold
		continue;
	end
	
	vecCounter(intType) = vecCounter(intType) + 1;
	matMeanRate = sData.matMeanRate;
	[intQuantiles,intStimNum,intStimCompN] = size(sData.cellLRActPerQ);
	[intNeuronNum,intTrialNum] = size(matMeanRate);
	
	%% aggregate data
	matBC(intType,vecCounter(intType)) = sData.dblBC;
	matMDV(intType,vecCounter(intType)) = sData.dblMaxDevFrac;
	cellAggLRActPerQ(:,:,:,intType,vecCounter(intType)) = sData.cellLRActPerQ;
	cellAggPopMuPerQ(:,:,:,intType,vecCounter(intType)) = sData.cellPopMuPerQ;
end

%% plot
%pre-allocate
%intUseRec = 1:31;
cellUseLRActPerQ = cellAggLRActPerQ;%(:,:,:,:,intUseRec);
cellUsePopMuPerQ = cellAggPopMuPerQ;
intRecs = size(cellUseLRActPerQ,5);
matR_Discr=nan(3,intRecs);
matR_MuVar=nan(3,intRecs);
hAggFig = figure;
hAggAx = axes();
hold on;
for intType=1:numel(cellTypes)
	strType = cellTypes{intType};
	
	%% average over all orthogonal (or adjacent?) stimuli
	if intType == 3
		dblStep = 0.1;
		vecBinE = -4:dblStep:4;
	else
		dblStep = 0.1;
		vecBinE = -4:dblStep:4;
	end
	vecBinC = vecBinE(2:end)-dblStep/2;
	figure;maxfig;
	subplot(2,3,1);
	hold on
	intMidQ = ceil(intQuantiles/2);
	dblAggMu1 = mean(flat(cellfun(@mean,cellUseLRActPerQ(:,:,1,:,:))));
	dblAggMu2 = mean(flat(cellfun(@mean,cellUseLRActPerQ(:,:,2,:,:))));
	cellMidMu1 = num2cell(repmat(cellfun(@mean,cellUseLRActPerQ(intMidQ,:,1,:,:)),[5 1 1 1 1]));
	cellMidMu2 = num2cell(repmat(cellfun(@mean,cellUseLRActPerQ(intMidQ,:,2,:,:)),[5 1 1 1 1]));
	cellMidSd1 = num2cell(repmat(cellfun(@std,cellUseLRActPerQ(intMidQ,:,1,:,:)),[5 1 1 1 1]));
	cellMidSd2 = num2cell(repmat(cellfun(@std,cellUseLRActPerQ(intMidQ,:,2,:,:)),[5 1 1 1 1]));
	
	%cellUseLRActPerQ(:,:,1,:,:)  = cellfun(@(x,m,s) dblAggMu1+((x-m)/s),cellUseLRActPerQ(:,:,1,:,:),cellMidMu1,cellMidSd1,'UniformOutput',false);
	%cellUseLRActPerQ(:,:,2,:,:)  = cellfun(@(x,m,s) dblAggMu2+((x-m)/s),cellUseLRActPerQ(:,:,2,:,:),cellMidMu2,cellMidSd2,'UniformOutput',false);
	
	for intQ=1:intQuantiles
		%plot distros
		vecAct1 = cell2vec(cellUseLRActPerQ(intQ,:,1,intType,:));
		vecAct2 = cell2vec(cellUseLRActPerQ(intQ,:,2,intType,:));
		vecCounts1 = histcounts(vecAct1,vecBinE);
		vecCounts2 = histcounts(vecAct2,vecBinE);
		plot(vecBinC,0.8*(vecCounts1/max(vecCounts1))+intQ,'Color',[1 0 0]);
		plot(vecBinC,0.8*(vecCounts2/max(vecCounts2))+intQ,'Color',[0 0 1]);
	end
	%finish plot
	hold off;
	set(gca,'ytick',0.5 + (1:intQuantiles),'yticklabel',1:5);
	ylabel('Quantile; y=trials per quantile');
	xlabel('LR activation');
	title(sprintf('%s; Mean over adj stim pairs',strType));
	fixfig;grid off
	
	%plot d', variance and distance in mean
	matPooledSd = nan(intQuantiles,intStimNum,intRecs);
	matMeanD = nan(intQuantiles,intStimNum,intRecs);
	matPopMu = nan(intQuantiles,intStimNum,intRecs);
	
	matQ = nan(intQuantiles,intStimNum,intRecs);
	for intRec=1:intRecs
		for intQ=1:intQuantiles
			for intOriIdx = 1:intStimNum
				vecR1 = cellAggLRActPerQ{intQ,intOriIdx,1,intType,intRec};
				vecR2 = cellAggLRActPerQ{intQ,intOriIdx,2,intType,intRec};
				%matPooledVar(intQ,intOriIdx,intRec) = (var(vecR1) + var(vecR2))/2;
				matPooledSd(intQ,intOriIdx,intRec) = ((std(vecR1) + std(vecR2))/2);
				matMeanD(intQ,intOriIdx,intRec)  = abs(mean(vecR1) - mean(vecR2));
				matQ(intQ,intOriIdx,intRec) = intQ;
				
				dblPopMu1 = mean(cellUsePopMuPerQ{intQ,intOriIdx,1,intType,intRec});
				dblPopMu2 = mean(cellUsePopMuPerQ{intQ,intOriIdx,2,intType,intRec});
				matPopMu(intQ,intOriIdx,intRec) = (dblPopMu1 + dblPopMu2)/2;
			end
		end
	end
	matDprime = matMeanD ./ matPooledSd;
	matColMap = redbluepurple(intQuantiles);
	matColor2 = matColMap(matQ(:),:);
	
	h=subplot(2,3,2);
	colormap(h,matColMap);
	%scatter(mean(matPooledSd,2),mean(matMeanD,2),[],matColMap)
	%calc mean+sem per q
	matSdV = mean(matPooledSd,2);
	matDp = mean(matDprime,2);
	if intPlotMu == 1
		matMuV = mean(matMeanD,2);
	else
		matMuV = mean(matPopMu,2);
	end
	indRem = any(matDp > 100,1) | any(isnan(matSdV),1) | any(isnan(matDp),1) | any(isnan(matMuV),1);
	matSdV(:,:,indRem) = [];
	matDp(:,:,indRem) = [];
	matMuV(:,:,indRem) = [];
	intRecs = sum(~indRem);
	
	vecMeanDprime = mean(matDp,3);
	vecSemDprime = std(matDp,[],3)./sqrt(intRecs);
	vecMeanSd = mean(matSdV,3);
	vecSemSdL = (std(matSdV,[],3)./sqrt(intRecs))/2;
	vecSemSdR = (std(matSdV,[],3)./sqrt(intRecs))/2;
	vecMeanMu = mean(matMuV,3);
	vecSemMu = std(matMuV,[],3)./sqrt(intRecs);
	vecSemMuL = vecSemMu/2;
	vecSemMuR = vecSemMu/2;
	hold on
	cline(h,vecMeanMu',vecMeanSd',[],1:5);
	for intQ=1:intQuantiles
		errorbar(vecMeanMu(intQ),vecMeanSd(intQ),vecSemSdL(intQ),vecSemSdR(intQ),vecSemMuL(intQ),vecSemMuR(intQ),'x','color',matColMap(intQ,:));
	end
	hold off
	ylabel('Pooled sd over trials');
	xlabel([strMu ' over trials']);
	title('Point = stim+quantile mu+/-sem');
	xlim([min(0,min(get(gca,'xlim'))) max(get(gca,'xlim'))]);
	ylim([min(0,min(get(gca,'ylim'))) max(get(gca,'ylim'))]);
	fixfig;
	
	h=subplot(2,3,3);
	colormap(h,matColMap);
	hold on
	cline(h,vecMeanMu,vecMeanDprime,[],1:5);
	for intQ=1:intQuantiles
		errorbar(vecMeanMu(intQ),vecMeanDprime(intQ),vecSemDprime(intQ)/2,vecSemDprime(intQ)/2,vecSemMuL(intQ),vecSemMuR(intQ),'x','color',matColMap(intQ,:));
	end
	hold off
	xlabel([strMu ' over trials']);
	ylabel('Discriminability (d'')');
	title('Point = stim+quantile mu+/-sem');
	xlim([min(0,min(get(gca,'xlim'))) max(get(gca,'xlim'))]);
	ylim([min(0,min(get(gca,'ylim'))) max(get(gca,'ylim'))]);
	fixfig;
	
	
	%corrs per rec
	matFanoV = squeeze(matDp);
	matSdV = squeeze(matSdV);
	matMuV = squeeze(matMuV);
	for intRec=1:intRecs
		matR_Discr(intType,intRec) = corr(flat(matMeanD(:,:,intRec)),flat(matDprime(:,:,intRec)));
		matR_MuVar(intType,intRec) = corr(flat(matMeanD(:,:,intRec)),flat(matPooledSd(:,:,intRec)));
	end
	
	subplot(2,3,4)
	dblStep = 0.1;
	vecBinsE = -1:dblStep:1;
	vecBinsC = vecBinsE(2:end)-dblStep/2;
	vecBinsMV = histcounts(matR_MuVar(intType,:),vecBinsE);
	plot(vecBinsC,vecBinsMV)
	ylabel('# of recordings');
	xlabel('Pearson correlation mean/sd');
	title('muvar');
	fixfig;
	
	subplot(2,3,5)
	dblStep = 0.25;
	vecBinsE = -1:dblStep:1;
	vecBinsC = vecBinsE(2:end)-dblStep/2;
	vecBinsMD = histcounts(matR_Discr(intType,:),vecBinsE);
	[h,pD]=ttest(matR_Discr(intType,:));
	plot(vecBinsC,vecBinsMD)
	ylabel('# of recordings');
	xlabel('Pearson correlation mean vs d''');
	title(sprintf('Pearson r(mu,d''), mu=%.3f, p=%.2e',nanmean(matR_Discr(intType,:)),pD));
	fixfig;
	
	if intType == 2
		subplot(2,3,6);
		[h,p_RealShuff]=ttest(matR_Discr(1,:),matR_Discr(2,:));
		title(sprintf('T-test Pearson r(real) vs r(shuff),p=%.4f',p_RealShuff));
		fixfig;
	end
	if boolSaveFigs
		%% save fig
		export_fig(fullpath(strFigurePath,sprintf('2C_ABI_%s_Agg3_GainInvariance_%s.tif',strStim,strType)));
		export_fig(fullpath(strFigurePath,sprintf('2C_ABI_%s_Agg3_GainInvariance_%s.pdf',strStim,strType)));
	end
	
	%% add plot to agg fig
	axes(hAggAx);
	colormap(hAggAx,matColMap);
	cline(hAggAx,vecMeanMu,vecMeanDprime,[],1:5);
	for intQ=1:intQuantiles
		errorbar(hAggAx,vecMeanMu(intQ),vecMeanDprime(intQ),vecSemDprime(intQ)/2,vecSemDprime(intQ)/2,vecSemMuL(intQ),vecSemMuR(intQ),'x','color',matColMap(intQ,:));
	end
	text(hAggAx,vecMeanMu(end),vecMeanDprime(end),strType);
end
% finish agg fig
axes(hAggAx);
hold off
xlabel([strMu ' over trials']);
ylabel('Discriminability (d'')');
title('Point = stim+quantile mu+/-sem');
xlim([min(0,min(get(gca,'xlim'))) max(get(gca,'xlim'))]);
ylim([min(0,min(get(gca,'ylim'))) max(get(gca,'ylim'))]);
fixfig;

if boolSaveFigs
	%% save fig
	export_fig(fullpath(strFigurePath,sprintf('2C_ABI_%s_Agg4_GainInvariance.tif',strStim)));
	export_fig(fullpath(strFigurePath,sprintf('2C_ABI_%s_Agg4_GainInvariance.pdf',strStim)));
end
