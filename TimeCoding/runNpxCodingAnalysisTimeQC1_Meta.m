%% aim
%{
https://www.nature.com/articles/s41598-021-95037-z
https://www.nature.com/articles/s42003-021-02437-y

q1: is coding better in trials with higher or lower firing rates?

q2: can we define peaks in the IFR as population events, and find which cells spike in the beginning
or end? does this ordering differ between orientations?

%}

%% define qualifying data
%close all;
clear all;
strRunType = 'Npx'; %ABI or Npx or all?
strRunStim = 'DG';%DG or NM
intPlotMu = 3; %1=delta mu of LR axes + sd of LR axes, 2= pop mean + pop sd, 3= pop mean + sd of LR axes
boolSaveFigs = true;

%folders
if isfolder('F:\Drive\PopTimeCoding') && isfolder('F:\Data\Processed\Neuropixels\')
	strDataPath = 'F:\Data\Processed\Neuropixels\';
	strFigurePathSR = 'F:\Drive\PopTimeCoding\single_recs';
	strFigurePath = 'F:\Drive\PopTimeCoding\figures\';
	strTargetDataPath = 'F:\Drive\PopTimeCoding\data\';
else
	strDataPath = 'E:\DataPreProcessed\';
	strFigurePathSR = 'C:\Drive\PopTimeCoding\single_recs';
	strFigurePath = 'C:\Drive\PopTimeCoding\figures\';
	strTargetDataPath = 'C:\Drive\PopTimeCoding\data\';
end


%% find data
%cellTypes = {'Real','Shuff','Poiss','UniStretch','VarFixed','Saturating','TuneFixed'};
cellTypes = {'Real','TShuff','TPoiss','TSdScaling'};
if strcmp(strRunType,'ABI')
	strFindStrType = '_ABI';
elseif strcmp(strRunType,'Npx')
	strFindStrType = '_Rec';
else
	strFindStrType = '_';
end
if intPlotMu == 1
	strMu = '\DeltaMu on LR axis';
else
	strMu = 'Pop mean';
end
sDir = dir([strTargetDataPath 'QC1Data*' strFindStrType '*.mat']); %or ABA if old
%indUseRecs = contains({sDir.name},['ABI_' strStim]);
%sDir(~indUseRecs) = [];

matDecPerf_TrainAll = nan(5,0,0);
matDecPerf_TrainOnQ = nan(5,0,0);
cellAggLRActPerQ = cell([0 0 0 0 0]); %[quantile x stim x (this stim/adja stim) x type x rec]
cellAggPopMuPerQ = cell([0 0 0 0 0]);
vecCounter = zeros(1,numel(cellTypes));
for intFile=1:numel(sDir)
	%% load data
	strFolder = sDir(intFile).folder;
	strFile = sDir(intFile).name;
	cellSplit = strsplit(strFile,'_');
	sData = load(fullpath(strFolder,strFile));
	intType = find(ismember(cellTypes,sData.strType));
	if ~strcmp(sData.strRunStim,strRunStim) || size(sData.cellLRActPerQ,2) < 8 || any(flat(cellfun(@(x) any(isnan(x(:))),sData.cellLRActPerQ)))
		continue;
	end

	vecCounter(intType) = vecCounter(intType) + 1;
	matMeanRate = sData.matMeanRate;
	[intQuantiles,intStimNum,intStimCompN] = size(sData.cellLRActPerQ);
	[intNeuronNum,intTrialNum] = size(matMeanRate);

	vecSplitPerfMu = mean(sData.matSplitPerf,2)';%Train once on all, test per quantile
	vecQuantilePerf = sData.vecQuantilePerf;%Train+test per quantile

	matDecPerf_TrainAll(:,intType,vecCounter(intType)) = vecSplitPerfMu;
	matDecPerf_TrainOnQ(:,intType,vecCounter(intType)) = vecQuantilePerf;

	%% aggregate data
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
		dblStep = 1/4;
		vecBinE = -3:dblStep:3;
	else
		dblStep = 1/4;
		vecBinE = -3:dblStep:3;
	end
	dblIntegralFactor = numel(vecBinE)/range(vecBinE);
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
		plot(vecBinC,1*dblIntegralFactor*(vecCounts1/sum(vecCounts1))+intQ,'Color',[1 0 0]);
		plot(vecBinC,1*dblIntegralFactor*(vecCounts2/sum(vecCounts2))+intQ,'Color',[0 0 1]);
	end
	%finish plot
	hold off;
	set(gca,'ytick',0.5 + (1:intQuantiles),'yticklabel',1:5);
	ylabel('Quantile; y=trials per quantile');
	xlabel('LR activation');
	title(sprintf('%s; Mean over orth stim pairs',strType));
	fixfig;grid off

	%plot d', variance and distance in mean
	matPooledSd = nan(intQuantiles,intStimNum,intRecs);
	matMeanD = nan(intQuantiles,intStimNum,intRecs);
	matPopMu = nan(intQuantiles,intStimNum,intRecs);
	matPopSd = nan(intQuantiles,intStimNum,intRecs);

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

				dblPopSd1 = std(cellUsePopMuPerQ{intQ,intOriIdx,1,intType,intRec});
				dblPopSd2 = std(cellUsePopMuPerQ{intQ,intOriIdx,2,intType,intRec});
				matPopSd(intQ,intOriIdx,intRec) = (dblPopSd1 + dblPopSd2)/2;
			end
		end
	end
	matColMap = redbluepurple(intQuantiles);
	matColor2 = matColMap(matQ(:),:);

	h=subplot(2,3,2);
	colormap(h,matColMap);
	%scatter(mean(matPooledSd,2),mean(matMeanD,2),[],matColMap)
	%calc mean+sem per q
	if intPlotMu == 1
		matDprime = matMeanD ./ matPooledSd;
		matMuV = squeeze(mean(matMeanD,2));
		matSdV = squeeze(mean(matPooledSd,2));
		matDp = squeeze(mean(matDprime,2));
		
		indRem = any(matDp > 100,1) | any(isnan(matSdV),1) | any(isnan(matDp),1) | any(isnan(matMuV),1);
		matSdV(:,indRem) = [];
		matDp(:,indRem) = [];
		matMuV(:,indRem) = [];
	elseif intPlotMu == 2
		matDprime = matPopMu ./ matPopSd;
		matMuV = squeeze(mean(matPopMu,2));
		matSdV = squeeze(mean(matPopSd,2));
		matDp = squeeze(mean(matDprime,2));

		indRem = any(matDp > 100,1) | any(isnan(matSdV),1) | any(isnan(matDp),1) | any(isnan(matMuV),1);
		matSdV(:,indRem) = [];
		matDp(:,indRem) = [];
		matMuV(:,indRem) = [];
	elseif intPlotMu == 3
		matDprime = matMeanD ./ matPooledSd;
		matMuV = squeeze(mean(matPopMu,2));
		matSdV = squeeze(mean(matPooledSd,2));
		matDp = squeeze(mean(matDprime,2));

		indRem = any(matDp > 100,1) | any(isnan(matSdV),1) | any(isnan(matDp),1) | any(isnan(matMuV),1);
		matSdV(:,indRem) = [];
		matDp(:,indRem) = [];
		matMuV(:,indRem) = [];
	end
	intRecs = sum(~indRem);
%error is d prime correct? sd/mean in plot (2,3,2) seems different for for example SdFixed

	vecMeanDprime = mean(matDp,2);
	vecSemDprime = std(matDp,[],2)./sqrt(intRecs);
	vecMeanSd = mean(matSdV,2);
	vecSemSdL = (std(matSdV,[],2)./sqrt(intRecs))/2;
	vecSemSdR = (std(matSdV,[],2)./sqrt(intRecs))/2;
	vecMeanMu = mean(matMuV,2);
	vecSemMu = std(matMuV,[],2)./sqrt(intRecs);
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
		[h,p_RealShuff]=ttest(matR_MuVar(1,:),matR_MuVar(2,:));
		title(sprintf('T-test Pearson r(real) vs r(shuff),p=%.4f',p_RealShuff));
		fixfig;
	end
	if boolSaveFigs
		%% save fig
		export_fig(fullpath(strFigurePath,sprintf('QC1a_GainInvariance_%s_%s_%s.tif',strRunType,strRunStim,strType)));
		export_fig(fullpath(strFigurePath,sprintf('QC1a_GainInvariance_%s_%s_%s.pdf',strRunType,strRunStim,strType)));
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
	export_fig(fullpath(strFigurePath,sprintf('QC1b_GainInvariance_%s_%s.tif',strRunType,strRunStim)));
	export_fig(fullpath(strFigurePath,sprintf('QC1b_GainInvariance_%s_%s.pdf',strRunType,strRunStim)));
end

%% plot decoding
figure;maxfig;
dblChance = 1/numel(unique(sData.vecStimIdx));
h1=subplot(2,3,1);hold on;plot([1 5],[dblChance dblChance],'--','color',[0.5 0.5 0.5]);xlabel('Activity quantile');ylabel('Decoding accuracy');text(h1,5,dblChance,'Chance');
h2=subplot(2,3,2);hold on;plot([1 5],[dblChance dblChance],'--','color',[0.5 0.5 0.5]);xlabel('Activity quantile');ylabel('Decoding accuracy');text(h2,5,dblChance,'Chance');
h3=subplot(2,3,3);hold on;xlabel('Activity quantile');ylabel('\DeltaDecoding accuracy');
for intType=1%:numel(cellTypes)
	strType = cellTypes{intType};

	vecTrainAllMu = mean(matDecPerf_TrainAll(:,intType,:),3);
	vecTrainAllSem = std(matDecPerf_TrainAll(:,intType,:),[],3)./sqrt(intRecs);
	vecTrainOnQMu = mean(matDecPerf_TrainOnQ(:,intType,:),3);
	vecTrainOnQSem = std(matDecPerf_TrainOnQ(:,intType,:),[],3)./sqrt(intRecs);
	
	vecDiffMu = mean(matDecPerf_TrainAll(:,intType,:)-matDecPerf_TrainOnQ(:,intType,:),3);
	vecDiffSem = std(matDecPerf_TrainAll(:,intType,:)-matDecPerf_TrainOnQ(:,intType,:),[],3)./sqrt(intRecs);
	
	%% add plot to agg fig
	axes(h1);
	colormap(h1,matColMap);
	cline(h1,1:5,vecTrainAllMu,[],1:5);
	for intQ=1:intQuantiles
		errorbar(h1,intQ,vecTrainAllMu(intQ),vecTrainAllSem(intQ),'x','color',matColMap(intQ,:));
	end
	text(h1,5,vecTrainAllMu(end),strType);

	axes(h2);
	colormap(h2,matColMap);
	cline(h2,1:5,vecTrainOnQMu,[],1:5);
	for intQ=1:intQuantiles
		errorbar(h2,intQ,vecTrainOnQMu(intQ),vecTrainOnQSem(intQ),'x','color',matColMap(intQ,:));
	end
	text(h2,5,vecTrainOnQMu(end),strType);

	axes(h3);
	colormap(h3,matColMap);
	cline(h3,1:5,vecDiffMu,[],1:5);
	for intQ=1:intQuantiles
		errorbar(h3,intQ,vecDiffMu(intQ),vecDiffSem(intQ),'x','color',matColMap(intQ,:));
	end
	text(h3,5,vecDiffMu(end),strType);
end

title(h1,'Train on all, test per Q');
title(h2,'Train + test per quantile');
title(h3,'Generalization benefit');
fixfig;

if boolSaveFigs
	%% save fig
	export_fig(fullpath(strFigurePath,sprintf('QC1c_DecodingPerQuantile_%s_%s.tif',strRunType,strRunStim)));
	export_fig(fullpath(strFigurePath,sprintf('QC1c_DecodingPerQuantile_%s_%s.pdf',strRunType,strRunStim)));
end