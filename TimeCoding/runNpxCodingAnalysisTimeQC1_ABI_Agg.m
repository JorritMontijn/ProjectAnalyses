%% aim
%{
https://www.nature.com/articles/s41598-021-95037-z
https://www.nature.com/articles/s42003-021-02437-y

q1: can single-neuron variability be explained as population-level gain multiplication?
> estimate tuning curve from real data, then apply trial-by-trial gain to all neurons
>model predicting firing rate as combination of mean-rate multiplied by tuning curve => what do
residuals look like?
A: gain axis and mean-rate axis do not align, but ori tuning is distributed as conic manifold around
gain axis. Using stim-specific gain coupling calculated cross-validated for each neuron, allows pop
response to be predicted with R^2 of 0.72

%see articles:
https://elifesciences.org/articles/8998
https://www.nature.com/articles/nn.3711
https://www.jneurosci.org/content/39/37/7344.abstract
etc

%}
%% define qualifying areas
clearvars -except sAggABI;
%{'VISal'}    {'VISam'}    {'VISl'}    {'VISp'}    {'VISpm'}    {'VISrl'}
strRunArea = 'VISp';%'posteromedial visual area' 'Primary visual area'
cellUseAreas = {strRunArea};
strRunStim = 'NM'; %DG,NM

vecRandomize = 1:3; %1=real data, 2=shuffled, 3=generated
boolMakeFigs = true;
boolSaveFigs = true;
boolSaveData = true;
boolHome = true;
if boolHome
	strDataPath = 'F:\Data\Processed\AllenBrainVisualEphys\Aggregates\';
	strFigurePath = 'F:\Data\Results\PopTimeCoding\figures\';
	strTargetDataPath = 'F:\Data\Results\PopTimeCoding\data\';
else
	strDataPath = 'D:\Data\Processed\AllenBrainVisualEphys\Aggregates\';
	strFigurePath = 'D:\Data\Results\PopTimeCoding\figures\';
	strTargetDataPath = 'D:\Data\Results\PopTimeCoding\data\';
end

%% select all neurons in LP and drifting grating stimuli
%strUseRec = '20191211_MP3_RunDriftingGratingsR01_g0_t0';
%strUseRec = '20191216_MP3_RunNaturalMovieR01_g0_t0';
if ~exist('sAggABI','var') || isempty(sAggABI)
	fprintf('Loading Allen Brain EcEphys data... [%s]\n',getTime);
	sLoad = load([strDataPath 'AggSes2022-03-28.mat']);
	sAggABI = sLoad.sSes;
	clear sLoad;
end
vecUseRec = 1:numel(sAggABI);%find(contains({sAggABA.Exp},'MP'));

%% pre-allocate matrices
intAreas = numel(cellUseAreas);
dblStartT = 0;
dblBimoThreshold = 0.4;
dblDevThreshold = 0.7;
intProjType = 2; %1=train+project per quantile; 2=train overall + project per quantile

%% go through recordings
tic
for intRec=1:numel(vecUseRec)
	% get matching recording data
	sRec = sAggABI(intRec);
	strRecOrig = sRec.Exp;
	
	strField = ['s' strRunStim];
	if strcmp(strRunStim,'DG') && isfield(sRec.structStimAgg,strField)
		%% DG
		% concatenate stimulus structures
		structStim = sRec.structStimAgg.sDG; %sDG, sNM, sNS
		
		%generate stimulus vectors
		vecStimLabels = mod(structStim.orientation(:)',360);
		vecStimOnTime = structStim.startT(:)';
		vecStimOffTime = structStim.stopT(:)';
		%remove nans
		indRem = isnan(vecStimLabels);
		vecStimLabels(indRem) = [];
		vecStimOnTime(indRem) = [];
		vecStimOffTime(indRem) = [];
		
		
	elseif strcmp(strRunStim,'NM') && isfield(sRec.structStimAgg,strField)
		%% NM
		% concatenate stimulus structures
		structStim = sRec.structStimAgg.sNM; %sDG, sNM, sNS
		vecOrigStartIdx = [1; 1+find(diff(structStim.frame)<0)];
		vecOrigStimOnTime = flat(structStim.startT(vecOrigStartIdx))';
		dblDur = median(diff(vecOrigStimOnTime));
		vecOrigStimOffTime = vecOrigStimOnTime+dblDur;
		intFramesInMovie = max(structStim.frame)+1;
		dblBinAvg = 10;
		intUseBinsInMovie = intFramesInMovie/dblBinAvg;
		dblBinRate = round(intUseBinsInMovie/dblDur);
		dblBinDur = 1/dblBinRate;
		vecBinEdges = linspace(0,dblBinDur*intUseBinsInMovie,intUseBinsInMovie+1);
		
		%generate fake stimulus vectors
		vecUniqueStims = 1:intUseBinsInMovie;
		vecFrameIdx = repmat(vecUniqueStims,[1 numel(vecOrigStimOnTime)]);
		vecStimLabels = (vecFrameIdx/max(vecFrameIdx))*180;
		vecStimOnTime = flat(vecBinEdges(1:(end-1))' + vecOrigStimOnTime)';
		vecStimOffTime = vecStimOnTime + dblBinDur;
	else
		continue;
	end
	%general prep
	[vecStimIdx,vecUniqueStims,vecRepNum,cellSelect,vecTrialRepetition] = val2idx(vecStimLabels);
	intRepNum = min(vecRepNum);
	indRem2 = vecTrialRepetition>intRepNum;
	vecStimLabels(indRem2) = [];
	vecStimOnTime(indRem2) = [];
	vecStimOffTime(indRem2) = [];
	[vecStimIdx,vecUniqueStims,vecRepNum,cellSelect,vecTrialRepetition] = val2idx(vecStimLabels);
	intTrialNum = numel(vecStimOnTime);
	intStimNum = numel(vecUniqueStims);
	
	%% select area 1
	for intArea=1:intAreas
		strArea = cellUseAreas{intArea};
		indArea1Neurons = contains(sRec.cellAreas,strArea,'IgnoreCase',true);
		if sum(indArea1Neurons) == 0, continue;end
		fprintf('Preparing %d/%d: %s, area %s [%s]\n',intRec,numel(vecUseRec),strRecOrig,strArea,getTime);
		
		%% prep recording
		runPopCodingPrep;
		if intNeuronNum < 25
			fprintf('Number of neurons is %d for %s: skipping... [%s]\n',intNeuronNum,strRecOrig,getTime);
			continue;
		end
		
		%% run data types
		%change name
		for intRandomize = vecRandomize
			if intRandomize == 1
				strType = 'Real';
			elseif intRandomize == 2
				strType = 'Shuff';
			elseif intRandomize == 3
				strType = 'Poiss';
			end
			strRec = [strType '_' strRecOrig];
			close all;
			fprintf('Running data type %s for %s [%s]\n',strType,strRecOrig,getTime);
			
			%simple "rate code"
			matSpikeCounts = cellfun(@(x) sum(x>dblUseStartT),cellSpikeTimesPerCellPerTrial);
			matMeanRate = matSpikeCounts./dblUseMaxDur;
			
			%randomize per orientation
			if intRandomize > 1
				for intStim=1:intStimNum
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
						else
							error
						end
					end
				end
			end
			
			%run approximation by mean-rate multiplied with tuning curve
			[matRespNSR,vecStimTypes,vecUnique] = getStimulusResponses(matMeanRate,vecStimIdx);
			matTuningCurves = mean(matRespNSR,3);
			vecMeanPopRate = mean(matMeanRate,1);
			vecNormPopRate = vecMeanPopRate ./ mean(vecMeanPopRate);
			
			%% plot population mean + sd over neurons, compare with neuron mean+sd over trials
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
			matGenRate = poissrnd(repmat(vecNeuronMean,[1 intTrialNum]));
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
			
			if boolMakeFigs
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
					export_fig(fullpath(strFigurePath,sprintf('2C1%s_PopSpikeStatistics_%s.tif',strRunStim,strRec)));
					export_fig(fullpath(strFigurePath,sprintf('2C1%s_PopSpikeStatistics_%s.pdf',strRunStim,strRec)));
				end
			end
			
			%% prep analysis
			%select balanced trials
			vecPopRate = sum(matMeanRate,1);
			vecStartTrials = round(intReps*linspace(1/intReps,(1+1/intReps),intQuantiles+1));
			vecStartTrials(end)=[];
			matUseTrials = nan(intStimNum,intSplitTrialsPerOri,intQuantiles);
			vecTrialQuantile = zeros(1,intTrialNum);
			for intQ=1:intQuantiles
				vecUseTrialsTemp = vecStartTrials(intQ):(vecStartTrials(intQ)+intSplitTrialsPerOri-1);
				for intStim=1:intStimNum
					vecThisStim = find(cellSelect{intStim});
					[vecSorted,vecReorder]=sort(vecPopRate(vecThisStim));
					vecQualifyingTrials = vecThisStim(vecReorder(vecUseTrialsTemp));
					matUseTrials(intStim,:,intQ) = vecQualifyingTrials;
				end
				%divide into quantiles
				vecUseTrials = sort(flat(matUseTrials(:,:,intQ)));
				vecTrialQuantile(vecUseTrials) = intQ;
			end
			
			%% example
			%for intUseStim=1:intStimNum
			intUseStim = 1;
			vecStim1 = find(vecStimIdx==intUseStim);
			vecStim2 = find(vecStimIdx==modx(round(intUseStim+intStimNum/2),intStimNum));
			vecUseTrials = cat(2,vecStim1,vecStim2);
			vecTQR = vecTrialQuantile(vecUseTrials);
			vecUseTrials(vecTQR==0)=[];
			vecTQR = vecTrialQuantile(vecUseTrials);
			vecStimR = vecStimIdx(vecUseTrials);
			vecPMR = vecPopMean(vecUseTrials);
			matStimPair = matMeanRate(:,vecUseTrials)';
			indStim1=vecStimR==min(vecStimR);
			indStim2=vecStimR==max(vecStimR);
			
			%% plot
			if boolMakeFigs && intUseStim == 1
				
				matColMap = redbluepurple(intQuantiles);
				matColor = matColMap(vecTQR,:);
				
				figure;maxfig;
				%make plots of distributions of activation per quantile
				subplot(2,3,1)
				title(sprintf('%s, %s', strRec,strRunArea),'interpreter','none');
				
				subplot(2,3,3)
				hold on;
			end
			dblStep = 1;
			vecBinE = (-10:dblStep:10)/10;
			vecBinC = vecBinE(2:end)-dblStep/2;
			vecAllAct = nan(1,numel(vecTQR));
			vecAbsW = nan(1,intQuantiles);
			vecBinaryPerf = nan(1,intQuantiles);
			beta0 = [0;0];
			cellLRActPerQ = cell(2,intQuantiles);
			
			if intProjType == 2
				matUseResp = matStimPair';
				vecUseOri = val2idx(vecStimR);
				[dblPerfP,vecDecodedIndexRateCV,matPosteriorProbabilityBin,dblMeanErrorDegsRate,matConfusion,matWeightsBin,matActivation,matAggWeights,vecAggRep] = ...
					doCrossValidatedDecodingLR(matUseResp,vecUseOri,intTypeCV,[],1);%dblLambda);
				%calculate normalization factors
				vecBinIdx = val2idx(vecUseOri);
				vecRepWeights = squeeze(sum(abs(matAggWeights(:,1,:)),1));
				vecNormFactors = vecRepWeights(vecAggRep);
			end
			
			for intQ=1:intQuantiles
				vecUseTrialsQ = vecTQR==intQ;
				vecUseOriQ = val2idx(vecStimR(vecUseTrialsQ));
				if intProjType == 1
					%% proj 1
					matUseResp = matStimPair(vecUseTrialsQ,:)';
					[dblPerfP,vecDecodedIndexRateCV,matPosteriorProbabilityBin,dblMeanErrorDegsRate,matConfusion,matWeightsBin,matActivation] = ...
						doCrossValidatedDecodingLR(matUseResp,vecUseOriQ,intTypeCV,[],1);%dblLambda);
					
					vecBinaryPerf(intQ) = dblPerfP;
					vecAbsW(intQ) = sum(abs(matWeightsBin(:,1)));
					vecAllAct(vecUseTrialsQ) = matActivation(1,1:sum(vecUseTrialsQ))/vecAbsW(intQ);
					
					%split by group & plot
					vecAct = matActivation(1,1:sum(vecUseTrialsQ))/vecAbsW(intQ);
				else
					%% proj 2
					%split by group & plot
					
					vecBinaryPerf(intQ) = sum(vecDecodedIndexRateCV(vecUseTrialsQ) == vecBinIdx(vecUseTrialsQ)) / sum(vecUseTrialsQ);
					vecAbsW(intQ) = mean(vecNormFactors(vecUseTrialsQ));
					vecAllAct(vecUseTrialsQ) = matActivation(1,1:sum(vecUseTrialsQ))/vecAbsW(intQ);
					
					%split by group & plot
					vecAct = matActivation(1,vecUseTrialsQ)'./vecNormFactors(vecUseTrialsQ);
				end
				%% plot
				vecCounts1 = histcounts(vecAct(vecUseOriQ==1),vecBinE);
				vecCounts2 = histcounts(vecAct(vecUseOriQ==2),vecBinE);
				if boolMakeFigs && intUseStim == 1
					plot(vecBinC,0.8*(vecCounts1/max(vecCounts1))+intQ,'Color',[1 0 0]);
					plot(vecBinC,0.8*(vecCounts2/max(vecCounts2))+intQ,'Color',[0 0 1]);
				end
				cellLRActPerQ{1,intQ} = vecAct(vecUseOriQ==1);
				cellLRActPerQ{2,intQ} = vecAct(vecUseOriQ==2);
			end
			
			if boolMakeFigs && intUseStim == 1
				hold off;
				set(gca,'ytick',0.5 + (1:intQuantiles),'yticklabel',1:5);
				ylabel('Quantile; y=trials per quantile');
				xlabel('LR activation');
				title('0 vs 90 degrees');
				fixfig;grid off
				
				%plot mean & decision axis
				subplot(2,3,2)
				hold on;
				scatter(vecPMR(indStim1),vecAllAct(indStim1),[],matColor(indStim1,:),'d');
				scatter(vecPMR(indStim2),vecAllAct(indStim2),[],matColor(indStim2,:),'*');
				hold off
				xlabel('Mean pop activity');
				ylabel('LR activation');
				title('point=trial, color=pop activity quantile');
				fixfig;
			end
			
			%% average over all orthogonal (or adjacent?) stimuli
			matBinaryPerf = nan(intQuantiles,intStimNum);
			cellLRActPerQ = cell(intQuantiles,intStimNum,2);
			
			[vecTrialTypeIdx,vecUnique,vecPriorDistribution,cellSelect,vecRepetition] = val2idx(vecStimIdx);
			
			%% get projection onto activation axis, type 1 or 2
			if intProjType == 1
				%% get axes and projection per quantile
				for intQ=1:intQuantiles
					for intOriIdx1 = 1:intStimNum
						intOriIdx2 = intOriIdx1+1;
						%intOriIdx2 = intOriIdx1+floor(intStimNum/2)-1;
						intOriIdx2 = modx(intOriIdx2,intStimNum);
						
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
				end
			else
				%% get axes from all data, then project per quantile
				for intOriIdx1 = 1:intStimNum
					intOriIdx2 = intOriIdx1+1;
					%intOriIdx2 = intOriIdx1+floor(intStimNum/2)-1;
					intOriIdx2 = modx(intOriIdx2,intStimNum);
					
					indUseTrainTrials = (vecTrialTypeIdx == intOriIdx1 | vecTrialTypeIdx == intOriIdx2);
					matUseResp = matMeanRate(:,indUseTrainTrials)';
					vecUseOri = vecTrialTypeIdx(indUseTrainTrials);
					[dblPerformanceCV,vecDecodedIndexCV,matPosteriorProbability,dblMeanErrorDegs,matConfusion,matWeights,matAggActivation,matAggWeights,vecAggRep] = ...
						doCrossValidatedDecodingLR(matUseResp,vecUseOri,intTypeCV,[],dblLambda);
					%calculate normalization factors
					vecBinIdx = val2idx(vecTrialTypeIdx(indUseTrainTrials));
					vecQuantInOri = vecTrialQuantile(indUseTrainTrials);
					vecRepWeights = squeeze(sum(abs(matAggWeights(:,1,:)),1));
					vecNormFactors = vecRepWeights(vecAggRep);
					
					for intQ=1:intQuantiles
						indUseTrials = vecQuantInOri==intQ;
						%split by group & plot
						matBinaryPerf(intQ,intOriIdx1) = sum(vecDecodedIndexCV(indUseTrials) == vecBinIdx(indUseTrials)) / sum(indUseTrials);
						vecAct = matAggActivation(1,indUseTrials)'./vecNormFactors(indUseTrials);
						vecAct1 = vecAct(vecUseOri(indUseTrials)==intOriIdx1);
						vecAct2 = vecAct(vecUseOri(indUseTrials)==intOriIdx2);
						if mean(vecAct1) > mean(vecAct2)
							[vecAct2,vecAct1]= swap(vecAct1,vecAct2);
						end
						cellLRActPerQ{intQ,intOriIdx1,1} = vecAct1;
						cellLRActPerQ{intQ,intOriIdx1,2} = vecAct2;
					end
				end
			end
			
			if boolMakeFigs
				dblStep = 1;
				vecBinE = -10:dblStep:10;
				vecBinC = vecBinE(2:end)-dblStep/2;
				
				subplot(2,3,4);
				hold on
				for intQ=1:intQuantiles
					%plot distros
					vecAct1 = cell2vec(cellLRActPerQ(intQ,:,1));
					vecAct2 = cell2vec(cellLRActPerQ(intQ,:,2));
					vecCounts1 = histcounts(vecAct1,vecBinE);
					vecCounts2 = histcounts(vecAct2,vecBinE);
					if boolMakeFigs
						plot(vecBinC,0.8*(vecCounts1/max(vecCounts1))+intQ,'Color',[1 0 0]);
						plot(vecBinC,0.8*(vecCounts2/max(vecCounts2))+intQ,'Color',[0 0 1]);
					end
				end
				%finish plot
				hold off;
				set(gca,'ytick',0.5 + (1:intQuantiles),'yticklabel',1:5);
				ylabel('Quantile; y=trials per quantile');
				xlabel('LR activation');
				title('Mean over ~orth stim pairs');
				fixfig;grid off
			
				%plot d', variance and distance in mean
				matDprime = nan(intQuantiles,intStimNum);
				matPooledSd = nan(intQuantiles,intStimNum);
				matMeanD = nan(intQuantiles,intStimNum);
				matQ = nan(intQuantiles,intStimNum);
				for intQ=1:intQuantiles
					for intOriIdx = 1:intStimNum
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
				vecSemDprime = std(matDprime,[],2)./sqrt(intStimNum);
				vecMeanSd = mean(matPooledSd,2);
				vecSemSd = std(matPooledSd,[],2)./sqrt(intStimNum);
				vecMeanMu = mean(matMeanD,2);
				vecSemMu = std(matMeanD,[],2)./sqrt(intStimNum);
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
					export_fig(fullpath(strFigurePath,sprintf('2C3%s_DynamicCoding_%s.tif',strRunStim,strRec)));
					export_fig(fullpath(strFigurePath,sprintf('2C3%s_DynamicCoding_%s.pdf',strRunStim,strRec)));
				end
			end
			
			%% save data
			if boolSaveData
				save([strTargetDataPath 'TimeCodingAggQC1ABI_' strRunStim strRec '_' strRunArea '.mat'],'strRec','strRunArea','strRunStim','matMeanRate','cellLRActPerQ','dblBC','dblMaxDevFrac');
			end
		end
	end
end
toc
