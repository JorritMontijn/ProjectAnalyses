%% aim
%{
https://www.nature.com/articles/s41598-021-95037-z
https://www.nature.com/articles/s42003-021-02437-y

q1: is coding better in trials with higher or lower firing rates?

q2: can we define peaks in the IFR as population events, and find which cells spike in the beginning
or end? does this ordering differ between orientations?

%}
%% define qualifying areas
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

%% find data
strStim = 'DG';%DG/NM
cellTypes = {'Real','Shuff','Poiss'};
sDir = dir([strTargetDataPath 'TimeCodingAggQC1ABI*.mat']);
indOri = contains({sDir.name},'ABI_Ori');
if strcmp(strStim,'DG')
	sDir(~indOri) = [];
elseif strcmp(strStim,'NM')
	sDir(indOri) = [];
else
	error;
end
cellAggLRActPerQ = cell([0 0 0 0 0]);
vecCounter = [0 0 0];
for intFile=1:numel(sDir)
	%% load data
	strFolder = sDir(intFile).folder;
	strFile = sDir(intFile).name;
	strType = strrep(strrep(getFlankedBy(strFile,'QC1ABI','_','first'),'ABI',''),'_Ori','');
	intType = find(ismember(cellTypes,strType));
	sData = load(fullpath(strFolder,strFile));
	if size(sData.cellLRActPerQ,2) < 8 || any(flat(cellfun(@(x) any(isnan(x(:))),sData.cellLRActPerQ)))
		continue;
	end
	vecCounter(intType) = vecCounter(intType) + 1;
	matMeanRate = sData.matMeanRate;
	[intQuantiles,intStimNum,intStimCompN] = size(sData.cellLRActPerQ);
	[intNeuronNum,intTrialNum] = size(matMeanRate);

	%% aggregate data
	cellAggLRActPerQ(:,:,:,intType,vecCounter(intType)) = sData.cellLRActPerQ;
end

%pre-allocate
intRecs = size(cellAggLRActPerQ,5);
matR_Discr=nan(3,intRecs);
matR_MuVar=nan(3,intRecs);
for intType=1:3
	strType = cellTypes{intType};
	
	%% average over all orthogonal (or adjacent?) stimuli
	dblStep = 0.5;
	vecBinE = -10:dblStep:10;
	vecBinC = vecBinE(2:end)-dblStep/2;
	figure;maxfig;
	subplot(2,3,1);
	hold on
	for intQ=1:intQuantiles
		%plot distros
		vecAct1 = cell2vec(cellAggLRActPerQ(intQ,:,1,intType,:));
		vecAct2 = cell2vec(cellAggLRActPerQ(intQ,:,2,intType,:));
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
	title(sprintf('%s; Mean over ~orth stim pairs',strType));
	fixfig;grid off
	
	%plot d', variance and distance in mean
	matDprime = nan(intQuantiles,intStimNum,intRecs);
	matPooledVar = nan(intQuantiles,intStimNum,intRecs);
	matMeanD = nan(intQuantiles,intStimNum,intRecs);
	matQ = nan(intQuantiles,intStimNum,intRecs);
	for intRec=1:intRecs
		for intQ=1:intQuantiles
			for intOriIdx = 1:intStimNum
				vecR1 = cellAggLRActPerQ{intQ,intOriIdx,1,intType,intRec};
				vecR2 = cellAggLRActPerQ{intQ,intOriIdx,2,intType,intRec};
				matDprime(intQ,intOriIdx,intRec) = abs(mean(vecR1) - mean(vecR2))/((var(vecR1) + var(vecR2))*0.5);
				%matDprime(intQ,intOriIdx,intRec) = abs(getdprime2(vecR1,vecR2));
				matPooledVar(intQ,intOriIdx,intRec) = (var(vecR1) + var(vecR2))/2;
				matMeanD(intQ,intOriIdx,intRec)  = abs(mean(vecR1) - mean(vecR2));
				matQ(intQ,intOriIdx,intRec) = intQ;
			end
		end
	end
	matColMap = redbluepurple(intQuantiles);
	matColor2 = matColMap(matQ(:),:);
	
	h=subplot(2,3,2);
	colormap(h,matColMap);
	%scatter(mean(matPooledSd,2),mean(matMeanD,2),[],matColMap)
	%calc mean+sem per q
	matVarV = mean(matPooledVar,2);
	matDp = mean(matDprime,2);
	matMuV = mean(matMeanD,2);
	indRem = any(isnan(matVarV),1) | any(isnan(matDp),1) | any(isnan(matMuV),1);
	matVarV(:,:,indRem) = [];
	matDp(:,:,indRem) = [];
	matMuV(:,:,indRem) = [];
	intRecs = sum(~indRem);
	
	vecMeanDprime = mean(matDp,3);
	vecSemDprime = std(matDp,[],3)./sqrt(intRecs);
	vecMeanSd = mean(matVarV,3);
	vecSemSdL = (std(matVarV,[],3)./sqrt(intRecs))/2;
	vecSemSdR = (std(matVarV,[],3)./sqrt(intRecs))/2;
	vecMeanMu = mean(matMuV,3);
	vecSemMu = std(matMuV,[],3)./sqrt(intRecs);
	hold on
	cline(h,vecMeanSd',vecMeanMu',[],1:5);
	for intQ=1:intQuantiles
		errorbar(vecMeanSd(intQ),vecMeanMu(intQ),vecSemMu(intQ)/2,vecSemMu(intQ)/2,vecSemSdL(intQ),vecSemSdR(intQ),'x','color',matColMap(intQ,:));
	end
	hold off
	xlabel('Pooled var over trials');
	ylabel('\DeltaMean over trials');
	title('Point = stim+quantile mu+/-sem');
	xlim([min(0,min(get(gca,'xlim'))) max(get(gca,'xlim'))]);
	ylim([min(0,min(get(gca,'ylim'))) max(get(gca,'ylim'))]);
	fixfig;
	
	h=subplot(2,3,3);
	colormap(h,matColMap);
	hold on
	cline(h,vecMeanSd,vecMeanDprime,[],1:5);
	for intQ=1:intQuantiles
		errorbar(vecMeanSd(intQ),vecMeanDprime(intQ),vecSemDprime(intQ)/2,vecSemDprime(intQ)/2,vecSemSdL(intQ),vecSemSdR(intQ),'x','color',matColMap(intQ,:));
	end
	hold off
	xlabel('Pooled var over trials');
	ylabel('Discriminability (mean/var)');
	title('Point = stim+quantile mu+/-sem');
	xlim([min(0,min(get(gca,'xlim'))) max(get(gca,'xlim'))]);
	ylim([min(0,min(get(gca,'ylim'))) max(get(gca,'ylim'))]);
	fixfig;
	
	%corrs per rec
	matFanoV = squeeze(matDp);
	matVarV = squeeze(matVarV);
	matMuV = squeeze(matMuV);
	for intRec=1:intRecs
		matR_Discr(intType,intRec) = corr(matMuV(:,intRec),matFanoV(:,intRec));
		matR_MuVar(intType,intRec) = corr(matMuV(:,intRec),matVarV(:,intRec));
	end
	
	subplot(2,3,4)
	dblStep = 0.1;
	vecBinsE = -1:dblStep:1;
	vecBinsC = vecBinsE(2:end)-dblStep/2;
	vecBinsMV = histcounts(matR_MuVar(intType,:),vecBinsE);
	plot(vecBinsC,vecBinsMV)
	ylabel('# of recordings');
	xlabel('Pearson correlation mean/var');
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
	xlabel('Pearson correlation mean vs (mean/var)');
	title(sprintf('Pearson r(mu,fano-1), mu=%.3f, p=%.2e',nanmean(matR_Discr(intType,:)),pD));
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
end
