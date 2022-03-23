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
sDir = dir([strTargetDataPath 'TimeCodingAggQC1*.mat']);
for intFile=1:numel(sDir)
	%% load data
	strFolder = sDir(intFile).folder;
	strFile = sDir(intFile).name;
	strType = getFlankedBy(strFile,'QC1','_','first');
	
	sData = load(fullpath(strFolder,strFile));
	cellAggRate = sData.cellAggRate;
	[intQuantiles,intStimNum,intIterNum] = size(sData.cellAggLRActPerQ,[1 2 4]);
	intCounter = 1;
	for intIter=1:intIterNum
		[intNeuronNum,intTrialNum] = size(cellAggRate{intIter});
		matMeanRate(intCounter:(intCounter+intNeuronNum-1),:) = cellAggRate{intIter};
		intCounter = intCounter + intNeuronNum;
	end
	
	%% make plot 1
	%real pop mean+sd
	vecPopMean = mean(matMeanRate,1);
	vecPopSd = std(matMeanRate,[],1);
	vecPopCv = vecPopSd./vecPopMean;
	vecPopFano = (vecPopSd.^2)./vecPopMean;
	
	vecFitX = linspace(min(vecPopMean),max(vecPopMean),100);
	[fitobject,gof] = fit(vecPopMean',vecPopSd','poly1');
	matCoefInt = confint(fitobject);
	dblSlope = fitobject.p1;
	vecSlopeCI = matCoefInt(:,1);
	dblIntercept = fitobject.p2;
	vecIntCI = matCoefInt(:,2);
	dblAdjR2 = gof.adjrsquare;
	[vecPredCI,vecPredY] = predint(fitobject,vecFitX);
	
	%real neuron mean+sd
	vecNeuronMean = mean(matMeanRate,2);
	vecNeuronSd = std(matMeanRate,[],2);
	vecNeuronCv = vecNeuronSd./vecNeuronMean;
	vecNeuronFano = (vecNeuronSd.^2)./vecNeuronMean;
	
	vecNeuronFitX = linspace(min(vecNeuronMean),max(vecNeuronMean),100);
	[fitobjectNeuron,gofNeuron] = fit(vecNeuronMean,vecNeuronSd,'poly1');
	matNeuronCoefInt = confint(fitobjectNeuron);
	dblNeuronSlope = fitobjectNeuron.p1;
	vecNeuronSlopeCI = matNeuronCoefInt(:,1);
	dblNeuronIntercept = fitobjectNeuron.p2;
	vecNeuronIntCI = matNeuronCoefInt(:,2);
	dblNeuronAdjR2 = gofNeuron.adjrsquare;
	[vecNeuronPredCI,vecNeuronPredY] = predint(fitobjectNeuron,vecNeuronFitX);
	
	figure;maxfig
	subplot(3,4,1)
	hold on
	scatter(vecPopMean,vecPopSd,[],lines(1),'.');
	plot(vecFitX,vecPredY,'k--');
	%plot(vecFitX,vecPredCI(:,1),'--','color',[0.5 0.5 0.5]);
	%plot(vecFitX,vecPredCI(:,2),'--','color',[0.5 0.5 0.5]);;
	hold off
	title(sprintf('%s,R^2-adj=%.3f,slope=%.3f +/- %.3f',strType,dblAdjR2,dblSlope,mean(abs(vecSlopeCI-dblSlope))));
	xlabel('Population mean rate (Hz)');
	ylabel('Pop sd (Hz)');fixfig;
	
	subplot(3,4,5)
	scatter(vecPopMean,vecPopCv,'.')
	xlabel('Population mean rate (Hz)');
	ylabel('Pop cv');fixfig;
	
	subplot(3,4,9)
	scatter(vecPopMean,vecPopFano,'.')
	xlabel('Population mean rate (Hz)');
	ylabel('Pop Fano factor');fixfig;
	
	subplot(3,4,2)
	hold on
	scatter(vecNeuronMean,vecNeuronSd,[],lines(1),'.');
	plot(vecNeuronFitX,vecNeuronPredY,'k--');
	hold off
	title(sprintf('R^2-adj=%.3f,slope=%.3f +/- %.3f',dblNeuronAdjR2,dblNeuronSlope,mean(abs(vecNeuronSlopeCI-dblNeuronSlope))));
	xlabel('Neuron mean rate (Hz)');
	ylabel('Neuron sd (Hz)');fixfig;
	
	subplot(3,4,6)
	scatter(vecNeuronMean,vecNeuronCv,'.')
	xlabel('Neuron mean rate (Hz)');
	ylabel('Neuron cv');fixfig;
	
	subplot(3,4,10)
	scatter(vecNeuronMean,vecNeuronFano,'.')
	xlabel('Neuron mean rate (Hz)');
	ylabel('Neuron Fano factor');fixfig;
	
	if boolSaveFigs
		%% save fig
		export_fig(fullpath(strFigurePath,sprintf('2CAgg1_PopSpikeStatistics_%s.tif',strType)));
		export_fig(fullpath(strFigurePath,sprintf('2CAgg1_PopSpikeStatistics_%s.pdf',strType)));
	end
	
	%% make plot 2
	vecQuantilePerfMu = mean(sData.matQuantilePerf,2);
	vecQuantilePerfSd = std(sData.matQuantilePerf,[],2);
	
	vecSplitQuantilePerfMu = mean(sData.matSplitQuantilePerf,2);
	vecSplitQuantilePerfSd = std(sData.matSplitQuantilePerf,[],2);
	
	figure;maxfig;
	subplot(2,3,4)
	hold on
	plot([1 numel(vecQuantilePerfMu)],(1/intStimNum)*[1 1],'--','color',[0.5 0.5 0.5]);
	errorbar(1:intQuantiles,vecSplitQuantilePerfMu,vecSplitQuantilePerfSd,'color',lines(1));
	xlabel('Pop. activity quantile');
	ylabel('CV decoding accuracy');
	title(sprintf('%s,Train once on all, test per quantile',strType));
	fixfig;
	
	%split train, output split
	subplot(2,3,5)
	hold on
	plot([1 numel(vecQuantilePerfMu)],(1/intStimNum)*[1 1],'--','color',[0.5 0.5 0.5]);
	errorbar(1:intQuantiles,vecQuantilePerfMu,vecQuantilePerfSd,'color',lines(1));
	xlabel('Pop. activity quantile');
	ylabel('CV decoding accuracy');
	title('Train+test per quantile');
	fixfig;
	normaxes;
	
	%difference
	%[phat3,pci3] = binofit((phat-vecQuantilePerf)*intTrialsPerQ,intTrialsPerQ,dblAlphaEquivOfSd);
	
	subplot(2,3,6)
	hold on
	%plot([1 numel(vecQuantilePerf)],(1/intStimNr)*[1 1],'--','color',[0.5 0.5 0.5]);
	plot(1:intQuantiles,(vecSplitQuantilePerfMu-vecQuantilePerfMu));
	xlabel('Pop. activity quantile');
	ylabel('Generalization penalty (\Deltaaccuracy)');
	title('Accuracy difference');
	fixfig;
	
	if boolSaveFigs
		%% save fig
		export_fig(fullpath(strFigurePath,sprintf('2CAgg2_QuantileDecoding_%s.tif',strType)));
		export_fig(fullpath(strFigurePath,sprintf('2CAgg2_QuantileDecoding_%s.pdf',strType)));
	end
	
	%% plot 3
	
	%% average over all orthogonal (or adjacent?) stimuli
	cellAggLRActPerQ = sData.cellAggLRActPerQ;
	cellLRActPerQ = cell(intQuantiles,intStimNum,2);
	dblStep = 0.25;
	vecBinE = -2:dblStep:2;
	vecBinC = vecBinE(2:end)-dblStep/2;
	figure;maxfig;
	subplot(2,3,4);
	hold on
	for intQ=1:intQuantiles
		%plot distros
		vecAct1 = cell2vec(cellAggLRActPerQ(intQ,:,1,:));
		vecAct2 = cell2vec(cellAggLRActPerQ(intQ,:,2,:));
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
	matDprime = nan(intQuantiles,intStimNum,intIterNum);
	matPooledSd = nan(intQuantiles,intStimNum,intIterNum);
	matMeanD = nan(intQuantiles,intStimNum,intIterNum);
	matQ = nan(intQuantiles,intStimNum,intIterNum);
	for intIter=1:intIterNum
		for intQ=1:intQuantiles
			for intOriIdx = 1:intStimNum
				matDprime(intQ,intOriIdx,intIter) = abs(getdprime2(cellAggLRActPerQ{intQ,intOriIdx,1,intIter},cellAggLRActPerQ{intQ,intOriIdx,2,intIter}));
				matPooledSd(intQ,intOriIdx,intIter) = (std(cellAggLRActPerQ{intQ,intOriIdx,1,intIter}) + std(cellAggLRActPerQ{intQ,intOriIdx,2,intIter}))/2;
				matMeanD(intQ,intOriIdx,intIter)  = abs(mean(cellAggLRActPerQ{intQ,intOriIdx,1,intIter}) - mean(cellAggLRActPerQ{intQ,intOriIdx,2,intIter}));
				matQ(intQ,intOriIdx,intIter) = intQ;
			end
		end
	end
	
	matColMap = redbluepurple(intQuantiles);
	matColor2 = matColMap(matQ(:),:);
	
	h=subplot(2,3,5);
	colormap(h,matColMap);
	%scatter(mean(matPooledSd,2),mean(matMeanD,2),[],matColMap)
	%calc mean+sem per q
	vecMeanDprime = mean(mean(matDprime,3),2);
	vecSemDprime = std(mean(matDprime,3),[],2)./sqrt(intStimNum);
	vecMeanSd = mean(mean(matPooledSd,3),2);
	vecSemSd = std(mean(matPooledSd,3),[],2)./sqrt(intStimNum);
	vecMeanMu = mean(mean(matMeanD,3),2);
	vecSemMu = std(mean(matMeanD,3),[],2)./sqrt(intStimNum);
	hold on
	cline(h,vecMeanSd,vecMeanMu,[],1:5);
	for intQ=1:intQuantiles
		errorbar(vecMeanSd(intQ),vecMeanMu(intQ),vecSemSd(intQ)/2,vecSemSd(intQ)/2,vecSemMu(intQ)/2,vecSemMu(intQ)/2,'x','color',matColMap(intQ,:));
	end
	hold off
	xlabel('Sd over trials');
	ylabel('Mean over trials');
	title('Point = stim+quantile mu+/-sem');
	fixfig;
	
	h=subplot(2,3,6);
	colormap(h,matColMap);
	hold on
	cline(h,vecMeanSd,vecMeanDprime,[],1:5);
	for intQ=1:intQuantiles
		errorbar(vecMeanSd(intQ),vecMeanDprime(intQ),vecSemSd(intQ)/2,vecSemSd(intQ)/2,vecSemDprime(intQ)/2,vecSemDprime(intQ)/2,'x','color',matColMap(intQ,:));
	end
	hold off
	xlabel('Sd over trials');
	ylabel('d''');
	title('Point = stim+quantile mu+/-sem');
	fixfig;
	
	if boolSaveFigs
		%% save fig
		export_fig(fullpath(strFigurePath,sprintf('2CAgg3_DynamicCoding_%s.tif',strType)));
		export_fig(fullpath(strFigurePath,sprintf('2CAgg3_DynamicCoding_%s.pdf',strType)));
	end
end
