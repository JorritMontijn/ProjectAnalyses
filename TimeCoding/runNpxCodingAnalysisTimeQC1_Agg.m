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
cellUseAreas = {...
	'Primary visual area',...
	...'posteromedial visual area',...
	};

boolSaveData = true;
boolSaveFigs = true;
boolHome = true;
if boolHome
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


%% select all neurons in LP and drifting grating stimuli
if ~exist('sAggStim','var') || isempty(sAggStim)
	[sAggStim,sAggNeuron,sAggSources]=loadDataNpx('','driftinggrating',strDataPath);
end
sAggNeuron(strcmpi({sAggNeuron.SubjectType},'DBA')) = [];

%% pre-allocate matrices
intAreas = numel(cellUseAreas);
matDecPerf = [];
dblStartT = 0;

%% go through recordings
tic
for intRandomize=1:3
	for intRec=1:numel(sAggStim)
		% get matching recording data
		strRec = sAggStim(intRec).Exp;
		sThisRec = sAggStim(strcmpi(strRec,{sAggStim(:).Exp}));
		sThisSource = sAggSources(strcmpi(strRec,{sAggSources(:).Exp}));
		
		%prep grating data
		[sUseNeuron,vecStimOnTime,vecStimOffTime,vecOrientation] = NpxPrepGrating(sAggNeuron,sThisRec,cellUseAreas);
		[vecOriIdx,vecUniqueOris,vecRepNum,cellSelect,vecTrialRepetition] = val2idx(vecOrientation);
		intTrialNum = numel(vecStimOnTime);
		intOriNum = numel(unique(vecOrientation));
		intRepNum = intTrialNum/intOriNum;
		if numel(sUseNeuron) == 0, continue;end
		
		%change name
		if intRandomize == 1
			strRec = ['Real_' strRec];
		elseif intRandomize == 2
			strRec = ['Shuff_' strRec];
		elseif intRandomize == 3
			strRec = ['Poiss_' strRec];
		end
		
		%% select area 1
		for intArea=1:intAreas
			strArea = cellUseAreas{intArea};
			indArea1Neurons = contains({sUseNeuron.Area},strArea,'IgnoreCase',true);
			if sum(indArea1Neurons) == 0, continue;end
			%% get orientation responses & single-trial population noise
			sArea1Neurons = sUseNeuron(indArea1Neurons);
			cellSpikeTimes = {sArea1Neurons.SpikeTimes};
			[matData,indTuned,indResp,cellSpikeTimes,sOut] = NpxPrepData(cellSpikeTimes,vecStimOnTime,vecStimOffTime,vecOrientation);
			vecOri180 = mod(vecOrientation,180)*2;
			intTunedN = sum(indTuned);
			intRespN = sum(indResp);
			
			%% stim timings
			dblStimDur = roundi(median(vecStimOffTime - vecStimOnTime),1,'ceil');
			dblPreTime = -dblStartT;%0.3;
			dblPostTime = 0;%0.3;
			dblMaxDur = dblStimDur+dblPreTime+dblPostTime;
			dblBinWidth = 0.05;
			vecBinEdges = 0:dblBinWidth:dblMaxDur;
			vecStimTime = vecBinEdges(2:end)-dblBinWidth/2 - dblPreTime;
			indStimBins = vecStimTime > 0 & vecStimTime < dblStimDur;
			intBinNum = numel(vecBinEdges)-1;
			matBNSR = nan(intBinNum,intRespN,intOriNum,intRepNum);
			matBNT = nan(intBinNum,intRespN,intTrialNum);
			%matBNT_shifted = nan(intBinNum,intRespN,intTrialNum);
			
			vecRepCounter = zeros(1,intOriNum);
			%get spikes per trial per neuron
			cellSpikeTimesPerCellPerTrial = cell(intRespN,intTrialNum);
			for intN=1:intRespN
				% build pseudo data, stitching stimulus periods
				[vecPseudoSpikeTimes,vecPseudoEventT] = getPseudoSpikeVectors(cellSpikeTimes{intN},vecStimOnTime-dblPreTime,dblMaxDur);
				
				%real
				[vecTrialPerSpike,vecTimePerSpike] = getSpikesInTrial(vecPseudoSpikeTimes,vecPseudoEventT,dblMaxDur);
				for intTrial=1:intTrialNum
					vecSpikeT = vecTimePerSpike(vecTrialPerSpike==intTrial);
					cellSpikeTimesPerCellPerTrial{intN,intTrial} = vecSpikeT;
				end
			end
			vecStimOnStitched = vecPseudoEventT;
			
			%% define quantiles and remove zero-variance neurons
			%constants
			[vecTrialTypeIdx,vecUnique,vecPriorDistribution,cellSelect,vecRepetition] = val2idx(vecOri180);
			intStimNr = numel(vecUnique);
			dblLambda = 1;%1
			intTypeCV = 2;
			dblUseStartT = 0;
			dblUseMaxDur = dblMaxDur-dblUseStartT;
			intUseMax = inf;
			intReps = mean(vecPriorDistribution);
			
			%remove zero-variance neurons
			matSpikeCounts_pre = cellfun(@(x) sum(x>dblUseStartT),cellSpikeTimesPerCellPerTrial);
			matMeanRate_pre = matSpikeCounts_pre./dblUseMaxDur;
			vecPopRate_pre = sum(matMeanRate_pre,1);
			
			intQuantiles = 5;
			vecStartTrials = round(intReps*linspace(1/intReps,(1+1/intReps),intQuantiles+1));
			vecStartTrials(end)=[];
			intSplitTrialsPerOri = min(floor(cellfun(@sum,cellSelect)/intQuantiles));
			indZeroVarNeurons = false(intRespN,1);
			for intQ=1:intQuantiles
				vecUseTrialsTemp = vecStartTrials(intQ):(vecStartTrials(intQ)+intSplitTrialsPerOri-1);
				for intOri=1:intStimNr
					vecThisOri = find(cellSelect{intOri});
					[vecSorted,vecReorder]=sort(vecPopRate_pre(vecThisOri));
					vecQualifyingTrials = vecThisOri(vecReorder(vecUseTrialsTemp));
					indZeroVarNeurons = indZeroVarNeurons | (var(matMeanRate_pre(:,vecQualifyingTrials),[],2) == 0);
				end
			end
			indZeroVarNeurons = false(size(indZeroVarNeurons));
			vecUseNeurons = find(~indZeroVarNeurons);
			vecRemNeurons = find(indZeroVarNeurons);
			cellSpikeTimesPerCellPerTrial(vecRemNeurons,:) = [];
			intRespN = numel(vecUseNeurons);
			
			%simple "rate code"
			matSpikeCounts = cellfun(@(x) sum(x>dblUseStartT),cellSpikeTimesPerCellPerTrial);
			matMeanRate = matSpikeCounts./dblUseMaxDur;
			
			%randomize per orientation
			if intRandomize > 1
				for intStim=1:intStimNr
					vecUseT = find(vecTrialTypeIdx==intStim);
					for intN=1:intRespN
						if intRandomize == 2
							%shuffle spikes
							matMeanRate(intN,vecUseT) = matMeanRate(intN,vecUseT(randperm(numel(vecUseT))));
						elseif intRandomize == 3
							%generate spikes
							dblMean = mean(cellfun(@(x) sum(x>dblUseStartT),cellSpikeTimesPerCellPerTrial(intN,vecUseT)));
							vecRates = poissrnd(dblMean,size(vecUseT));
							matMeanRate(intN,vecUseT) = vecRates;
						else error
							
						end
					end
				end
			end
			%matMeanRate = log(1+matMeanRate);
			
			%% plot population mean + sd over neurons, compare with neuron mean+sd over trials
			intTrials = size(matMeanRate,2);
			
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
			
			%pop mean+sd, generated according to poisson cells
			matGenRate = poissrnd(repmat(vecNeuronMean,[1 intTrials]));
			vecGenPopMean = mean(matGenRate,1);
			vecGenPopSd = std(matGenRate,[],1);
			vecGenPopCv = vecGenPopSd./vecGenPopMean;
			vecGenPopFano = (vecGenPopSd.^2)./vecGenPopMean;
			
			vecGenFitX = linspace(min(vecGenPopMean),max(vecGenPopMean),100);
			[fitobjectGen,gofGen] = fit(vecGenPopMean',vecGenPopSd','poly1');
			matGenCoefInt = confint(fitobjectGen);
			dblGenSlope = fitobjectGen.p1;
			vecGenSlopeCI = matGenCoefInt(:,1);
			dblGenIntercept = fitobjectGen.p2;
			vecGenIntCI = matGenCoefInt(:,2);
			dblGenAdjR2 = gofGen.adjrsquare;
			[vecGenPredCI,vecGenPredY] = predint(fitobjectGen,vecGenFitX);
			
			%neuron mean+sd, generated according to poisson cells
			vecGenNeuronMean = mean(matGenRate,2);
			vecGenNeuronSd = std(matGenRate,[],2);
			vecGenNeuronCv = vecGenNeuronSd./vecGenNeuronMean;
			vecGenNeuronFano = (vecGenNeuronSd.^2)./vecGenNeuronMean;
			
			vecGenNeuronFitX = linspace(min(vecGenNeuronMean),max(vecGenNeuronMean),100);
			[fitobjectGenNeuron,gofGenNeuron] = fit(vecGenNeuronMean,vecGenNeuronSd,'poly1');
			matGenNeuronCoefInt = confint(fitobjectGenNeuron);
			dblGenNeuronSlope = fitobjectGenNeuron.p1;
			vecGenNeuronSlopeCI = matGenNeuronCoefInt(:,1);
			dblGenNeuronIntercept = fitobjectGenNeuron.p2;
			vecGenNeuronIntCI = matGenNeuronCoefInt(:,2);
			dblGenNeuronAdjR2 = gofGenNeuron.adjrsquare;
			[vecGenNeuronPredCI,vecGenNeuronPredY] = predint(fitobjectGenNeuron,vecGenNeuronFitX);
			
			figure;maxfig
			subplot(3,4,1)
			hold on
			scatter(vecPopMean,vecPopSd,[],lines(1),'.');
			plot(vecFitX,vecPredY,'k--');
			%plot(vecFitX,vecPredCI(:,1),'--','color',[0.5 0.5 0.5]);
			%plot(vecFitX,vecPredCI(:,2),'--','color',[0.5 0.5 0.5]);;
			hold off
			title(sprintf('R^2-adj=%.3f,slope=%.3f +/- %.3f',dblAdjR2,dblSlope,mean(abs(vecSlopeCI-dblSlope))));
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
			
			subplot(3,4,3)
			hold on
			scatter(vecGenPopMean,vecGenPopSd,'.')
			plot(vecGenFitX,vecGenPredY,'k--');
			hold off
			title(sprintf('R^2-adj=%.3f,slope=%.3f +/- %.3f',dblGenAdjR2,dblGenSlope,mean(abs(vecGenSlopeCI-dblGenSlope))));
			xlabel('Gen Population mean rate (Hz)');
			ylabel('Gen Pop sd (Hz)');fixfig;
			
			subplot(3,4,7)
			scatter(vecGenPopMean,vecGenPopCv,'.')
			xlabel('Gen Population mean rate (Hz)');
			ylabel('Gen Pop cv');fixfig;
			
			subplot(3,4,11)
			scatter(vecGenPopMean,vecGenPopFano,'.')
			xlabel('Gen Population mean rate (Hz)');
			ylabel('Gen Pop Fano factor');fixfig;
			
			subplot(3,4,4)
			hold on
			scatter(vecGenNeuronMean,vecGenNeuronSd,[],lines(1),'.');
			plot(vecGenNeuronFitX,vecGenNeuronPredY,'k--');
			hold off
			title(sprintf('R^2-adj=%.3f,slope=%.3f +/- %.3f',dblGenNeuronAdjR2,dblGenNeuronSlope,mean(abs(vecGenNeuronSlopeCI-dblGenNeuronSlope))));
			xlabel('Gen Neuron mean rate (Hz)');
			ylabel('Gen Neuron sd (Hz)');fixfig;
			
			subplot(3,4,8)
			scatter(vecGenNeuronMean,vecGenNeuronCv,'.')
			xlabel('Gen Neuron mean rate (Hz)');
			ylabel('Gen Neuron cv');fixfig;
			
			subplot(3,4,12)
			scatter(vecGenNeuronMean,vecGenNeuronFano,'.')
			xlabel('Gen Neuron mean rate (Hz)');
			ylabel('Gen Neuron Fano factor');fixfig;
			
			if boolSaveFigs
				%% save fig
				export_fig(fullpath(strFigurePath,sprintf('2C1_PopSpikeStatistics_%s.tif',strRec)));
				export_fig(fullpath(strFigurePath,sprintf('2C1_PopSpikeStatistics_%s.pdf',strRec)));
			end
			
			%% run decoding
			%run samples with balanced trials
			intTrialsPerQ = intSplitTrialsPerOri*intStimNr;
			intSplitTrialsPerOri = min(floor(cellfun(@sum,cellSelect)/intQuantiles));
			intIters = 100;
			matSplitPerf = nan(intQuantiles,intIters);
			for intIter=1:intIters
				intIter
				%select balanced trials
				vecPopRate = sum(matMeanRate,1);
				vecUseBalancedTrials = nan(1,intSplitTrialsPerOri*intStimNr*intQuantiles);
				vecBalancedQ = nan(1,intSplitTrialsPerOri*intStimNr*intQuantiles);
				intCounterBT = 1;
				vecPriorDistributionSplit = ones(size(vecPriorDistribution))*intSplitTrialsPerOri;
				vecStartTrials = round(intReps*linspace(1/intReps,(1+1/intReps),intQuantiles+1));
				vecStartTrials(end)=[];
				matUseTrials = nan(intStimNr,intSplitTrialsPerOri,intQuantiles);
				for intQ=1:intQuantiles
					vecUseTrialsTemp = vecStartTrials(intQ):(vecStartTrials(intQ)+intSplitTrialsPerOri-1);
					for intOri=1:intStimNr
						vecThisOri = find(cellSelect{intOri});
						[vecSorted,vecReorder]=sort(vecPopRate(vecThisOri));
						vecQualifyingTrials = vecThisOri(vecReorder(vecUseTrialsTemp));
						matUseTrials(intOri,:,intQ) = vecQualifyingTrials;
						
						vecUseBalancedTrials(intCounterBT:(intCounterBT+intSplitTrialsPerOri-1)) = vecQualifyingTrials(randperm(intSplitTrialsPerOri,intSplitTrialsPerOri));
						vecBalancedQ(intCounterBT:(intCounterBT+intSplitTrialsPerOri-1)) = intQ;
						intCounterBT = intCounterBT + intSplitTrialsPerOri;
					end
				end
				
				%perform overall decoding
				[vecBalancedTrialTypeIdx,vecBalancedUnique,vecBalancedPriorDistribution,cellBalancedSelect,vecBalancedRepetition] = val2idx(vecOri180(vecUseBalancedTrials));
				[dblPerf,vecDecodedIndexCV,matPosteriorProbability] = ...
					doCrossValidatedDecodingLR(matMeanRate(:,vecUseBalancedTrials),vecOri180(vecUseBalancedTrials),intTypeCV,[],dblLambda);
				vecConfidence = nan(size(vecDecodedIndexCV));
				for intTrial=1:numel(vecBalancedTrialTypeIdx)
					vecConfidence(intTrial) = matPosteriorProbability(vecBalancedTrialTypeIdx(intTrial),intTrial);
				end
				
				%assign per quantile
				matCorrPerQ = nan(intTrialsPerQ,intQuantiles);
				vecCorr = vecDecodedIndexCV(:) == vecBalancedTrialTypeIdx(:);
				for intQ=1:intQuantiles
					matCorrPerQ(:,intQ) = vecCorr(vecBalancedQ==intQ);
				end
				
				%save data
				matSplitPerf(:,intIter) = sum(matCorrPerQ,1)/intTrialsPerQ;
			end
			
			%calc mean and sem
			vecSplitPerfMu = mean(matSplitPerf,2)';
			vecSplitPerfSem = std(matSplitPerf,[],2)';
			
			%perform decoding per quantile
			vecQuantilePerf = nan(1,intQuantiles);
			matQuantileConfusion = nan(intStimNr,intStimNr,intQuantiles);
			%split trials into quantiles
			figure;maxfig;
			vecPlotQ = [];%[1 ceil(intQuantiles/2) intQuantiles];
			vecTrialQuantile = zeros(1,intTrials);
			for intQ=1:intQuantiles
				%divide into quantiles
				vecUseTrials = sort(flat(matUseTrials(:,:,intQ)));
				vecTrialQuantile(vecUseTrials) = intQ;
				
				matUseRate = matMeanRate(:,vecUseTrials);
				vecOriUse = vecOri180(vecUseTrials);
				
				%[dblPerfQ,vecDecodedIndexRateCV,matPosteriorProbabilityRate,dblMeanErrorDegsRate,matConfusion,matWeightsRate] = ...
				%	doCrossValidatedDecodingLR(matUseRate,vecOriUse,intTypeCV,vecPriorDistributionSplit,dblLambda);
				[dblPerfQ,vecDecodedIndexCV,matPosteriorProbability,dblMeanErrorDegs,matConfusion,matWeights,matAggActivation] = ...
					doCrossValidatedDecodingLR(matUseRate,vecOriUse,intTypeCV,[],dblLambda);
				vecQuantilePerf(intQ) = dblPerfQ;
				matQuantileConfusion(:,:,intQ) = matConfusion;
				intPlot = find(vecPlotQ==intQ);
				if ~isempty(intPlot)
					subplot(2,3,intPlot)
					imagesc(matConfusion)
					xlabel('Real Orientation')
					ylabel('Decoded Orientation')
					title(sprintf('Quantile %d',intQ))
					fixfig;grid off;
				end
			end
			
			% overall train, output split
			clf;
			subplot(2,3,4)
			hold on
			plot([1 numel(vecQuantilePerf)],(1/intStimNr)*[1 1],'--','color',[0.5 0.5 0.5]);
			errorbar(1:intQuantiles,vecSplitPerfMu,vecSplitPerfSem,'color',lines(1));
			xlabel('Pop. activity quantile');
			ylabel('CV decoding accuracy');
			title('Train once on all, test per quantile');
			fixfig;
			
			%split train, output split
			dblAlphaEquivOfSd = normcdf(1)-normcdf(-1);
			[phat2,pci2] = binofit(vecQuantilePerf*intTrialsPerQ,intTrialsPerQ,dblAlphaEquivOfSd);
			subplot(2,3,5)
			hold on
			plot([1 numel(vecQuantilePerf)],(1/intStimNr)*[1 1],'--','color',[0.5 0.5 0.5]);
			errorbar(1:intQuantiles,phat2,phat2'-pci2(:,1),phat2'-pci2(:,2),'color',lines(1));
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
			plot(1:intQuantiles,(vecSplitPerfMu-vecQuantilePerf));
			xlabel('Pop. activity quantile');
			ylabel('Generalization penalty (\Deltaaccuracy)');
			title('Accuracy difference');
			fixfig;
			
			if boolSaveFigs
				%% save fig
				export_fig(fullpath(strFigurePath,sprintf('2C2_QuantileDecoding_%s.tif',strRec)));
				export_fig(fullpath(strFigurePath,sprintf('2C2_QuantileDecoding_%s.pdf',strRec)));
			end
			
			%% plot PCA
			vecStim1 = find(vecOri180==0);
			vecStim2 = find(vecOri180==90);
			vecUseTrials = cat(2,vecStim1,vecStim2);
			vecTQR = vecTrialQuantile(vecUseTrials);
			vecUseTrials(vecTQR==0)=[];
			vecTQR = vecTrialQuantile(vecUseTrials);
			vecOriR = vecOri180(vecUseTrials);
			vecPMR = vecPopMean(vecUseTrials);
			matForPCA = matMeanRate(:,vecUseTrials)';
			[coeff,score,latent,tsquared,explained,mu] = pca(matForPCA);
			vecFracExplained = cumsum(explained)/sum(explained);
			matReduced = matForPCA(:,vecFracExplained<0.9);
			
			%% plot decoding as a function of PCs
			%{
		vecPerfPerDim = size(score,2);
		for intDim=1:size(score,2)
			[dblPerfP,vecDecodedIndexRateCV,matPosteriorProbabilityBin,dblMeanErrorDegsRate,matConfusion,matWeightsBin,matActivation] = ...
				doCrossValidatedDecodingLR(score(:,1:intDim),vecOriR,intTypeCV,[],dblLambda);%dblLambda);
			vecPerfPerDim(intDim) = dblPerfP;
		end
		[matW, dblMSE, intRank, sSuppOut] = doRidgeRRR(matX,matY,intRank,dblLambda,intUseVersion);
			%}
			
			%% plot PCA
			matColMap = redbluepurple(intQuantiles);
			matColor = matColMap(vecTQR,:);
			%matColor(vecOriR==0,2) = 0.5;
			indOri1 = vecOriR==0;
			indOri2 = vecOriR==90;
			
			figure;maxfig;
			subplot(2,3,1);
			colormap(redbluepurple(intQuantiles));
			hold on
			scatter(score(indOri1,1),score(indOri1,2),[],matColor(indOri1,:),'d');
			scatter(score(indOri2,1),score(indOri2,2),[],matColor(indOri2,:),'*');
			hold off
			xlabel('PC1');
			ylabel('PC2');
			legend({'0 degs','90 degs'},'location','best');
			title('point=trial, color=pop activity quantile');
			fixfig;
			
			%make plots of distributions of activation per quantile
			subplot(2,3,3)
			dblStep = 0.25;
			vecBinE = -1.5:dblStep:1.5;
			vecBinC = vecBinE(2:end)-dblStep/2;
			hold on;
			vecAllAct = nan(1,numel(vecTQR));
			vecAbsW = nan(1,intQuantiles);
			vecBinaryPerf = nan(1,intQuantiles);
			beta0 = [0;0];
			cellLRActPerQ = cell(2,intQuantiles);
			for intQ=1:intQuantiles
				vecUseTrialsQ = vecTQR==intQ;
				matUseResp = matForPCA(vecUseTrialsQ,:)';
				vecUseOri = val2idx(vecOriR(vecUseTrialsQ));
				[dblPerfP,vecDecodedIndexRateCV,matPosteriorProbabilityBin,dblMeanErrorDegsRate,matConfusion,matWeightsBin,matActivation] = ...
					doCrossValidatedDecodingLR(matUseResp,vecUseOri,intTypeCV,[],1);%dblLambda);
				
				vecBinaryPerf(intQ) = dblPerfP;
				vecAbsW(intQ) = sum(abs(matWeightsBin(:,1)));
				vecAllAct(vecUseTrialsQ) = matActivation(1,1:sum(vecUseTrialsQ))/vecAbsW(intQ);
				
				%split by group & plot
				vecAct = matActivation(1,1:sum(vecUseTrialsQ))/vecAbsW(intQ);
				vecCounts1 = histcounts(vecAct(vecUseOri==1),vecBinE);
				vecCounts2 = histcounts(vecAct(vecUseOri==2),vecBinE);
				plot(vecBinC,0.8*(vecCounts1/max(vecCounts1))+intQ,'Color',[1 0 0]);
				plot(vecBinC,0.8*(vecCounts2/max(vecCounts2))+intQ,'Color',[0 0 1]);
				cellLRActPerQ{1,intQ} = vecAct(vecUseOri==1);
				cellLRActPerQ{2,intQ} = vecAct(vecUseOri==2);
			end
			hold off;
			set(gca,'ytick',0.5 + (1:intQuantiles),'yticklabel',1:5);
			ylabel('Quantile; y=trials per quantile');
			xlabel('LR activation');
			title('0 vs 90 degrees');
			fixfig;grid off
			
			%plot mean & decision axis
			subplot(2,3,2)
			hold on;
			scatter(vecPMR(indOri1),vecAllAct(indOri1),[],matColor(indOri1,:),'d');
			scatter(vecPMR(indOri2),vecAllAct(indOri2),[],matColor(indOri2,:),'*');
			hold off
			xlabel('Mean pop activity');
			ylabel('LR activation');
			title('point=trial, color=pop activity quantile');
			fixfig;
			
			%% average over all orthogonal (or adjacent?) stimuli
			intOris = numel(vecUnique);
			matBinaryPerf = nan(intQuantiles,intOris);
			cellLRActPerQ = cell(intQuantiles,intOris,2);
			
			[vecTrialTypeIdx,vecUnique,vecPriorDistribution,cellSelect,vecRepetition] = val2idx(vecOri180);
			dblStep = 0.25;
			vecBinE = -2:dblStep:2;
			vecBinC = vecBinE(2:end)-dblStep/2;
			subplot(2,3,4);
			hold on
			for intQ=1:intQuantiles
				for intOriIdx1 = 1:intOris
					intOriIdx2 = intOriIdx1+floor(intOris/2)-1;
					intOriIdx2 = modx(intOriIdx2,intOris);
					
					vecUseTrials = vecTrialQuantile==intQ & (vecTrialTypeIdx == intOriIdx1 | vecTrialTypeIdx == intOriIdx2);
					
					matUseResp = matMeanRate(:,vecUseTrials)';
					vecUseOri = vecTrialTypeIdx(vecUseTrials);
					[dblPerfP,vecDecodedIndexRateCV,matPosteriorProbabilityBin,dblMeanErrorDegsRate,matConfusion,matWeightsBin,matActivation] = ...
						doCrossValidatedDecodingLR(matUseResp,vecUseOri,intTypeCV,[],dblLambda);
					
					
					%split by group & plot
					matBinaryPerf(intQ,intOriIdx1) = dblPerfP;
					vecAct = matActivation(1,:)/sum(abs(matWeightsBin(:,1)));
					vecAct1 = vecAct(vecUseOri==intOriIdx1);
					vecAct2 = vecAct(vecUseOri==intOriIdx2);
					if mean(vecAct1) > mean(vecAct2)
						[vecAct2,vecAct1]= swap(vecAct1,vecAct2);
					end
					cellLRActPerQ{intQ,intOriIdx1,1} = vecAct1;
					cellLRActPerQ{intQ,intOriIdx1,2} = vecAct2;
				end
				
				%plot distros
				vecAct1 = cell2vec(cellLRActPerQ(intQ,:,1));
				vecAct2 = cell2vec(cellLRActPerQ(intQ,:,2));
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
			title('Mean over ~orth stim pairs');
			fixfig;grid off
			
			%plot d', variance and distance in mean
			matDprime = nan(intQuantiles,intOris);
			matPooledSd = nan(intQuantiles,intOris);
			matMeanD = nan(intQuantiles,intOris);
			matQ = nan(intQuantiles,intOris);
			for intQ=1:intQuantiles
				for intOriIdx = 1:intOris
					matDprime(intQ,intOriIdx) = abs(getdprime2(cellLRActPerQ{intQ,intOriIdx,1},cellLRActPerQ{intQ,intOriIdx,2}));
					matPooledSd(intQ,intOriIdx) = (std(cellLRActPerQ{intQ,intOriIdx,1}) + std(cellLRActPerQ{intQ,intOriIdx,2}))/2;
					matMeanD(intQ,intOriIdx)  = abs(mean(cellLRActPerQ{intQ,intOriIdx,1}) - mean(cellLRActPerQ{intQ,intOriIdx,2}));
					matQ(intQ,intOriIdx) = intQ;
				end
			end
			
			matColMap = redbluepurple(intQuantiles);
			matColor2 = matColMap(matQ(:),:);
			%{
		subplot(2,3,4)
		scatter(matDprime(:),matBinaryPerf(:),[],matColor2)
		xlabel('d''');
		ylabel('Ori-pair decoding accuracy');
		title('Point = stim+quantile combination');
		fixfig;
			%}
			%{
		subplot(2,3,5)
		scatter(matPooledSd(:),matBinaryPerf(:),[],matColor2)
		xlabel('Sd over trials');
		ylabel('Ori-pair decoding accuracy');
		title('Point = stim+quantile combination');
		fixfig;
			%}
			
			h=subplot(2,3,5);
			colormap(h,matColMap);
			%scatter(mean(matPooledSd,2),mean(matMeanD,2),[],matColMap)
			%calc mean+sem per q
			vecMeanDprime = mean(matDprime,2);
			vecSemDprime = std(matDprime,[],2)./sqrt(intOris);
			vecMeanSd = mean(matPooledSd,2);
			vecSemSd = std(matPooledSd,[],2)./sqrt(intOris);
			vecMeanMu = mean(matMeanD,2);
			vecSemMu = std(matMeanD,[],2)./sqrt(intOris);
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
				export_fig(fullpath(strFigurePath,sprintf('2C3_DynamicCoding_%s.tif',strRec)));
				export_fig(fullpath(strFigurePath,sprintf('2C3_DynamicCoding_%s.pdf',strRec)));
			end
			
			%% save data
			if boolSaveData
				save([strTargetDataPath 'TimeCodingAggQC1' strRec '.mat'],...
					'strRec','strArea','matMeanRate','matBinaryPerf','vecQuantilePerf','matSplitPerf','cellLRActPerQ');
			end
		end
	end
end
toc
