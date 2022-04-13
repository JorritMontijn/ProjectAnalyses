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

% new to do:
- non-stationary data exclusion
- sanity check: minder trials in eigen data
- change spiking rate to spike count inclusion threshold
- gain axes over alle data berekenen, dan projecteren per quantile
- probeer adjacent stimuli voor gain invariance analyse

%}
%% define qualifying areas
clearvars -except sAggABI;
%{'VISal'}    {'VISam'}    {'VISl'}    {'VISp'}    {'VISpm'}    {'VISrl'}
strRunArea = 'VISp';%'posteromedial visual area' 'Primary visual area'
strRunStim = 'NM'; %DG, NM
cellUseAreas = {strRunArea};

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
			
			%% get CV model prediction
			fprintf('   Running model predictions [%s]\n',getTime);
			%leave one repetition of one neuron out
			matTunePredCV = nan(size(matMeanRate));
			matMeanPredCV = nan(size(matMeanRate));
			matGainPredCV = nan(size(matMeanRate));
			matGain1PredCV = nan(size(matMeanRate));
			for intTestRep = 1:intRepNum
				indTestTrials = vecTrialRepetition==intTestRep;
				vecTestTrials = find(indTestTrials);
				vecTrainTrials = find(~indTestTrials);
				matTrainR = matMeanRate(:,vecTrainTrials);
				vecTrainS = vecStimIdx(vecTrainTrials);
				vecTestS = vecStimIdx(vecTestTrials);
				matTestR = matMeanRate(:,vecTestTrials);
				dblMeanPopRate = mean(matTrainR(:));
				%get tuning curves
				[matRespNSR,vecStimTypes,vecUniqueDegs] = getStimulusResponses(matTrainR,vecTrainS);
				matTrainTuning = mean(matRespNSR,3);
				%get gain for training trials
				for intTestNeuron=1:intNeuronNum
					%% get vars
					indOtherNeurons = ~ismember(1:intNeuronNum,intTestNeuron);
					matOtherR = matMeanRate(indOtherNeurons,vecTestTrials);
					vecTrainOtherRate = mean(matTrainR(indOtherNeurons,:),1)';
					vecTestOtherRate = mean(matOtherR,1)';
					matTrainOtherR = matTrainR(indOtherNeurons,:);
					matTestOtherR = matTestR(indOtherNeurons,:);
					
					%% simple prediction: only tuning curve
					matTunePredCV(intTestNeuron,vecTestTrials) = matTrainTuning(intTestNeuron,vecTestS);
					
					%% classic prediction: tuning curve times population mean
					%get pop mean per trial, and average pop mean over trials
					vecPopMeanTrain = mean(matTrainOtherR,1);
					dblAvgPopMeanTrain = mean(vecPopMeanTrain);
					
					%get test pop mean
					vecPopMeanTest =  mean(matTestOtherR,1);
					
					%get mean resp during train
					vecTrainR = matTrainTuning(intTestNeuron,vecTestS);
					
					%predict: R(phi,i,t)= (pop_mean(not i,t) / avg(pop_mean(not i))) * mu(phi,i)
					matMeanPredCV(intTestNeuron,vecTestTrials) = (vecPopMeanTest ./ dblAvgPopMeanTrain) .* vecTrainR;
					
					%{
					%% old
					%calculate optimal coupling
					vecTrainSelfAct = (matTrainR(~indOtherNeurons,:)-matTrainTuning(intTestNeuron,vecTrainS)) ./ (matTrainTuning(intTestNeuron,vecTrainS)+eps);
					vecMeanOffset = cat(2,vecTrainOtherRate,ones(size(vecTrainOtherRate)))\vecTrainSelfAct';
					
					%optimal pred
					vecMeanPred = matTrainTuning(intTestNeuron,vecTestS)' + matTrainTuning(intTestNeuron,vecTestS)'.*(vecMeanOffset(1).*vecTestOtherRate + vecMeanOffset(2));
					matMeanPredCV(intTestNeuron,indTestTrials) = vecMeanPred;
					
					
					%direct pred
					%vecOtherPopNormRate = mean(matOtherR)./dblMeanPopRate;
					%matMeanPredCV(intTestNeuron,indTestTrials) = matTrainTuning(intTestNeuron,vecTestOri) .* vecOtherPopNormRate;
					%}
					
					%% same gain for all stims
					%calculate gain axis for this stimulus from training trials
					vecTrainGain1Ax = mean(matTrainOtherR,2);
					
					%get population test responses and project onto training gain axis
					vecOtherGain1=getProjOnLine(matTestOtherR,vecTrainGain1Ax);
					
					%get mean resp during train
					vecTrainR = matTrainTuning(intTestNeuron,vecTestS);
					
					%predict: R(phi,i,t)= (gain1(not i,t) / norm(mu(not i))) * mu(phi,i)
					matGain1PredCV(intTestNeuron,vecTestTrials) = (vecOtherGain1 ./ norm(vecTrainGain1Ax)) .* vecTrainR';
					
					%{
					vecGain1Ax = mean(matTrainTuning(indOtherNeurons,:),2);
					vecTrainOtherGain1=getProjOnLine(matTrainR(indOtherNeurons,:),vecGain1Ax);
					vecTestOtherGain1=getProjOnLine(matOtherR,vecGain1Ax);
					
					%calculate optimal coupling
					vecTrainSelfAct = (matTrainR(~indOtherNeurons,:)-matTrainTuning(intTestNeuron,vecTrainS)) ./ (matTrainTuning(intTestNeuron,vecTrainS)+eps);
					vecGainOffset = cat(2,vecTrainOtherGain1,ones(size(vecTrainOtherGain1)))\vecTrainSelfAct';
					
					vecGain1Pred = matTrainTuning(intTestNeuron,vecTestOri)' + matTrainTuning(intTestNeuron,vecTestOri)'.*(vecGainOffset(1).*vecTestOtherGain1 + vecGainOffset(2));
					matGain1PredCV(intTestNeuron,indTestTrials) = vecGain1Pred;
					
					%dblCoupling=corr(vecTrainOtherGain1Norm',vecTrainSelfActNorm');
					%matGain1PredCV(intTestNeuron,indTestTrials) = matTrainTuning(intTestNeuron,vecTestOri) .* vecTestOtherGain1Norm*dblCoupling;
					%}
					
					%% gain prediction: stim-specific coupling between population gain and predicted neuron
					%R(phi,i,t)= (gain(phi,not i,t) / norm(mu(phi,not i))) * mu(phi,i)
					for intStimIdx=1:intStimNum
						%calculate gain axis for this stimulus from training trials
						indTrainT = vecTrainS==intStimIdx;
						matTrainOtherR = matTrainR(indOtherNeurons,indTrainT);
						vecTrainGainAx = mean(matTrainOtherR,2);
						
						%get population test responses and project onto training gain axis
						indTestT = vecTestS==intStimIdx;
						matTestOtherR = matTestR(indOtherNeurons,indTestT);
						vecOtherGain=getProjOnLine(matTestOtherR,vecTrainGainAx);
						
						%get mean response for this stimulus for this neuron over training trials
						dblTrainSelf = mean(matTrainR(~indOtherNeurons,indTrainT));
						
						%predict: R(phi,i,t)= (gain(phi,not i,t) / norm(mu(phi,not i))) * mu(phi,i)
						matGainPredCV(intTestNeuron,vecTestTrials(indTestT)) = (vecOtherGain ./ norm(vecTrainGainAx)) * dblTrainSelf;
					end
					
					%{
					%% old
					%get stimulus-specific gain
					vecOthersGain = nan(1,numel(vecTrainS));
					vecGain = nan(1,intStimNum);
					vecOffset = nan(1,intStimNum);
					for intOriIdx=1:intStimNum
						indUseT = vecTrainS==intOriIdx;
						matTrainOriR = matTrainR(indOtherNeurons,indUseT);
						vecTrainGainAx = mean(matTrainOriR,2);
						vecOtherGain=getProjOnLine(matTrainOriR,vecTrainGainAx);
						
						%linear fit
						vecGainOffset = cat(2,vecOtherGain,ones(size(vecOtherGain)))\matTrainR(intTestNeuron,indUseT)';
						vecGain(intOriIdx) = vecGainOffset(1);
						vecOffset(intOriIdx) = vecGainOffset(2);
					end
					
					%calculate location along gain axis
					vecGainAx = mean(matTrainTuning(indOtherNeurons,:),2);
					vecGainV=getProjOnLine(matOtherR,vecGainAx)';
					vecGainPred = vecGain(vecTestOri).*vecGainV + vecOffset(vecTestOri);
					matGainPredCV(intTestNeuron,vecTestTrials) = vecGainPred;
					%}
				end
			end
			
			[dblR2_Tune,dblSS_tot,dblSS_res,dblT_Tune,dblP_Tune,dblR2_adjusted_Tune,dblR2_SE_Tune] = getR2(matMeanRate(:),matTunePredCV(:),numel(matTrainTuning));
			[dblR2_Mean,dblSS_tot,dblSS_res,dblT_Mean,dblP_Mean,dblR2_adjusted_Mean,dblR2_SE_Mean] = getR2(matMeanRate(:),matMeanPredCV(:),numel(matTrainTuning));
			[dblR2_Gain,dblSS_tot,dblSS_res,dblT_Gain,dblP_Gain,dblR2_adjusted_Gain,dblR2_SE_Gain] = getR2(matMeanRate(:),matGainPredCV(:),numel(matTrainTuning));
			[dblR2_Gain1,dblSS_tot,dblSS_res,dblT_Gain1,dblP_Gain1,dblR2_adjusted_Gain1,dblR2_SE_Gain1] = getR2(matMeanRate(:),matGain1PredCV(:),numel(matTrainTuning));
			
			%prediction per neuron
			vecR2_Tune = nan(1,intNeuronNum);
			vecR2_Mean = nan(1,intNeuronNum);
			vecR2_Gain = nan(1,intNeuronNum);
			vecR2_Gain1 = nan(1,intNeuronNum);
			for intN=1:intNeuronNum
				%save prediction
				[dblR2_TuneN,dblSS_tot,dblSS_res,dblT,dblP,dblR2_adjusted_TuneN] = getR2(matMeanRate(intN,:),matTunePredCV(intN,:),numel(matTrainTuning(intN,:)));
				[dblR2_MeanN,dblSS_tot,dblSS_res,dblT,dblP,dblR2_adjusted_MeanN] = getR2(matMeanRate(intN,:),matMeanPredCV(intN,:),numel(matTrainTuning(intN,:)));
				[dblR2_GainN,dblSS_tot,dblSS_res,dblT,dblP,dblR2_adjusted_GainN] = getR2(matMeanRate(intN,:),matGainPredCV(intN,:),numel(matTrainTuning(intN,:)));
				[dblR2_Gain1N,dblSS_tot,dblSS_res,dblT,dblP,dblR2_adjusted_Gain1N] = getR2(matMeanRate(intN,:),matGain1PredCV(intN,:),numel(matTrainTuning(intN,:)));
				vecR2_Tune(intN) = dblR2_TuneN;
				vecR2_Mean(intN) = dblR2_MeanN;
				vecR2_Gain(intN) = dblR2_GainN;
				vecR2_Gain1(intN) = dblR2_Gain1N;
			end
			dblPredPerNeuronTune = mean(vecR2_Tune);
			dblPredPerNeuronMean = mean(vecR2_Mean);
			dblPredPerNeuronGain = mean(vecR2_Gain);
			dblPredPerNeuronGain1 = mean(vecR2_Gain1);
			
			%stim decoding; warning, takes very long (~half an hour)
			%[dblPerformanceCV,vecDecodedIndexCV,matPosteriorProbability,dblMeanErrorDegs,matConfusion,matWeights,matAggActivation] = ...
			%	doCrossValidatedDecodingLR(matMeanRate,vecStimIdx,intTypeCV,vecPriorDistribution,dblLambda);
			
			%%
			dblStep = 4;
			vecX = matMeanRate(:);
			vecY = matMeanPredCV(:);
			vecBinsX = 0:dblStep:81 - 0.5*dblStep;
			vecBinsY = 0:dblStep:80 - 0.5*dblStep;
			matCounts = histcounts2(vecX,vecY,vecBinsX,vecBinsY);
			vecLimC = [0 80];
			vecLimC_d = [-40 40];
			if boolMakeFigs
				figure;maxfig
				%real
				colormap(redwhite);
				subplot(2,4,1)
				imagesc(matMeanRate,vecLimC)
				colorbar;
				ylabel('Neuron');
				xlabel('Trial');
				title(sprintf('%s; Recorded spiking',strType));
				h=colorbar;
				ylabel(h,'Spiking rate (Hz)');
				fixfig;grid off;
				
				%summary
				subplot(2,4,5)
				errorbar(1:4,[dblR2_Tune dblR2_Mean dblR2_Gain1 dblR2_Gain],[dblR2_SE_Tune dblR2_SE_Mean dblR2_SE_Gain1 dblR2_SE_Gain],'x','CapSize',20)
				set(gca,'xtick',1:4,'xticklabel',{'Tune','Tune*Mean','Tune*Gain1','Tune*Gain'});
				xtickangle(20);
				ylabel('CV pop resp predictability (R^2)');
				xlim([0 5]);
				title(sprintf('%s',strRec),'interpreter','none');
				fixfig;
				
				%tuning only
				subplot(2,4,2)
				imagesc(matTunePredCV,vecLimC)
				colorbar;
				ylabel('Neuron');
				xlabel('Trial');
				title(sprintf('CV model: R = ori tuning'));
				h=colorbar;
				ylabel(h,'CV-predicted rate (Hz)');
				fixfig;grid off;
				
				h=subplot(2,4,6);
				imagesc(matMeanRate-matTunePredCV,vecLimC_d)
				h.Colormap=redblue;
				colorbar;
				ylabel('Neuron');
				xlabel('Trial');
				title(sprintf('Resids; adj-R^2=%.3f; T=%.1f',dblR2_adjusted_Tune,dblT_Tune));
				h=colorbar;
				ylabel(h,'Prediction error (Hz)');
				fixfig;grid off;
				
				%tune*mean
				subplot(2,4,3)
				imagesc(matMeanPredCV,vecLimC)
				colorbar;
				ylabel('Neuron');
				xlabel('Trial');
				title(sprintf('CV model: R = ori tuning * pop-mean'));
				h=colorbar;
				ylabel(h,'CV-predicted rate (Hz)');
				fixfig;grid off;
				
				h=subplot(2,4,7);
				imagesc(matMeanRate-matMeanPredCV,vecLimC_d)
				h.Colormap=redblue;
				colorbar;
				ylabel('Neuron');
				xlabel('Trial');
				title(sprintf('Resids; adj-R^2=%.3f; T=%.1f',dblR2_adjusted_Mean,dblT_Mean));
				h=colorbar;
				ylabel(h,'Prediction error (Hz)');
				fixfig;grid off;
				
				%tune*gain
				subplot(2,4,4)
				imagesc(matGainPredCV,vecLimC)
				colorbar;
				ylabel('Neuron');
				xlabel('Trial');
				title(sprintf('CV model: R = ori tuning * pop-gain'));
				h=colorbar;
				ylabel(h,'CV-predicted rate (Hz)');
				fixfig;grid off;
				
				h=subplot(2,4,8);
				imagesc(matMeanRate-matGainPredCV,vecLimC_d)
				h.Colormap=redblue;
				colorbar;
				ylabel('Neuron');
				xlabel('Trial');
				title(sprintf('Resids; adj-R^2=%.3f; T=%.1f',dblR2_adjusted_Gain,dblT_Gain));
				h=colorbar;
				ylabel(h,'Prediction error (Hz)');
				fixfig;grid off;
				
				if boolSaveFigs
					%% save fig
					export_fig(fullpath(strFigurePath,sprintf('2B1O_PopRespPrediction_%s.tif',strRec)));
					export_fig(fullpath(strFigurePath,sprintf('2B1O_PopRespPrediction_%s.pdf',strRec)));
				end
			end
			
			%% use correlations to predict residuals
			%{
		matResids = matMeanRate-matGainPredCV;
		matCorr = corr(matResids');
		matCorr(diag(diag(true(size(matCorr))))) = nan;
		[matResids_hat,dblR2_CV,matB] = doCrossValidatedDecodingRR(matCorr,matResids,1,1);
		
		figure
		subplot(2,3,1)
		imagesc(matCorr)
		
		subplot(2,3,2)
		imagesc(matB)
			%}
			
			%% project
			%for a given cloud of pop responses (n=trials), what fraction of noise aligns with the
			%mean-rate axis? what proportion with f'? how does this compare to a random direction (=spectral decomposition of covariance matrix)?
			
			%run projections
			fprintf('   Running multi-dim analysis [%s]\n',getTime);
			vecGainAx = mean(matTuningCurves,2);
			matGainResp = nan(size(matMeanRate));
			matCov = nan(intNeuronNum,intNeuronNum,intStimNum);
			vecSdProjMean = nan(1,intStimNum);
			vecSdProjGain1 = nan(1,intStimNum);
			vecSdProjGain = nan(1,intStimNum);
			vecSdProjRand = nan(1,intStimNum);
			vecSdProjOrth = nan(1,intStimNum);
			vecSdProjAdja = nan(1,intStimNum);
			vecDistFromOrigin = nan(1,intStimNum);
			vecDistFromOrthOri = nan(1,intStimNum);
			vecDistFromAdjaOri = nan(1,intStimNum);
			%reflections
			vecReflDistMuOrth = nan(1,intStimNum);
			vecReflDistMuAdja = nan(1,intStimNum);
			vecReflDistGainOrth = nan(1,intStimNum);
			vecReflDistGainAdja = nan(1,intStimNum);
			%cossim
			vecAngleMeanAndGain = nan(1,intStimNum);
			vecAngleMeanAndGainRand = nan(1,intStimNum);
			
			for intStimIdx=1:intStimNum
				vecOriTrials = find(vecStimIdx==intStimIdx);
				indAdjaOriTrials = vecStimIdx==modx(intStimIdx+1,intStimNum);
				indAdja2OriTrials = vecStimIdx==modx(intStimIdx-1,intStimNum);
				indOrthOriTrials = vecStimIdx==modx(intStimIdx+ceil(intStimNum/2)-1,intStimNum);
				matPoints =matMeanRate(:,vecOriTrials);
				vecOriR = mean(matPoints,2);
				vecAdjaOriR = mean(matMeanRate(:,indAdjaOriTrials),2);
				vecAdja2OriR = mean(matMeanRate(:,indAdja2OriTrials),2);
				vecOrthOriR = mean(matMeanRate(:,indOrthOriTrials),2);
				matGainResp(:,vecOriTrials) = vecOriR * vecNormPopRate(vecOriTrials);
				matCov(:,:,intStimIdx) = cov(matMeanRate(:,vecOriTrials)');
				
				vecProjMean=getProjOnLine(matPoints,ones(intNeuronNum,1));
				vecSdProjMean(intStimIdx) = std(vecProjMean);
				vecProjRand=getProjOnLine(matPoints,(rand(intNeuronNum,1)-0.5)*2);
				vecSdProjRand(intStimIdx) = std(vecProjRand);
				
				%Project activity unto single gain axis Vs one per stim and compare SDs
				%for each stim separately
				vecProjGain=getProjOnLine(matPoints,vecOriR);
				vecSdProjGain(intStimIdx) = std(vecProjGain);
				
				%for avg-gain
				vecProjGain1=getProjOnLine(matPoints,vecGainAx);
				vecSdProjGain1(intStimIdx) = std(vecProjGain1);
				
				vecProjOrth=getProjOnLine(matPoints-vecOrthOriR,vecOriR-vecOrthOriR);
				vecSdProjOrth(intStimIdx) = std(vecProjOrth);
				vecProjAdja=getProjOnLine(matPoints-vecAdjaOriR,vecOriR-vecAdjaOriR);
				vecSdProjAdja(intStimIdx) = std(vecProjAdja);
				
				
				vecDistFromOrigin(intStimIdx) = norm(vecOriR);
				vecDistFromOrthOri(intStimIdx) = norm(vecOriR-vecOrthOriR);
				vecDistFromAdjaOri(intStimIdx) = (norm(vecOriR-vecAdjaOriR) + norm(vecOriR-vecAdja2OriR))/2;
				
				vecGainDir = vecOriR;
				vecMeanDir = ones(size(vecOriR));
				[dblCosSim,dblAngSim] = cossim(vecGainDir,vecMeanDir);
				vecRandSim = nan(1,1);
				for i=1:1
					%[dblCosSimRand,dblAngSimRand] = cossim(vecOrthDir(randperm(intNeuronNum)),vecRefDir);
					%[dblCosSimRand,dblAngSimRand] = cossim((rand(intNeuronNum,1)-0.5)*2,vecMeanDir);
					[dblCosSimRand,dblAngSimRand] = cossim(rand(intNeuronNum,1),vecMeanDir);
					vecRandSim(i) = dblCosSimRand;
				end
				vecAngleMeanAndGain(intStimIdx) = dblCosSim;
				vecAngleMeanAndGainRand(intStimIdx) = mean(vecRandSim);
				
				%reflect around mean-axis
				intIters=100;
				vecTempReflDistMu = nan(1,intIters);
				vecTempReflDistMuAdja = nan(1,intIters);
				vecTempReflDistG = nan(1,intIters);
				vecTempReflDistGAdja = nan(1,intIters);
				for i=1:100
					%define variables
					vecReorder = randperm(numel(vecOriTrials));
					%vecOriTrials1 = vecOriTrials;%vecReorder(1:(floor(numel(vecReorder)/2)));
					%vecOriTrials2 = indAdjaOriTrials;%vecReorder(((ceil(numel(vecReorder)/2))+1):end);
					vecOriTrials1 = vecOriTrials(vecReorder(1:(floor(numel(vecReorder)/2))));
					vecOriTrials2 = vecOriTrials(vecReorder(((ceil(numel(vecReorder)/2))+1):end));
					indOrthOriTrials = vecStimIdx==modx(intStimIdx+ceil(intStimNum/2),intStimNum);
					vecX = mean(matMeanRate(:,vecOriTrials1),2);
					vecY = mean(matMeanRate(:,indOrthOriTrials),2);
					vecZ = mean(matMeanRate(:,vecOriTrials2),2);
					
					%reflect around mean
					[dummy,vecXprime]=getProjOnLine(vecX,ones(intNeuronNum,1));
					vecXrefl = 2*vecXprime - vecX;
					%compare with orthogonal
					dblReflDistMu = norm(vecXrefl-vecY);
					%compare with adjacent
					dblReflDistMuAdja = norm(vecXrefl-vecZ);
					
					
					%reflect around gain axis
					[dummy,vecXprimeG]=getProjOnLine(vecX,vecGainAx);
					vecXreflG = 2*vecXprimeG - vecX;
					%compare with orthogonal
					dblReflDistG = norm(vecXreflG-vecY);
					%compare with adjacent
					dblReflDistGAdja = norm(vecXreflG-vecZ);
					
					%save
					vecTempReflDistMu(i) = dblReflDistMu;
					vecTempReflDistMuAdja(i) = dblReflDistMuAdja;
					vecTempReflDistG(i) =dblReflDistG;
					vecTempReflDistGAdja(i) = dblReflDistGAdja;
				end
				
				%mean over resamples
				vecReflDistMuOrth(intStimIdx) = mean(vecTempReflDistMu(i));
				vecReflDistMuAdja(intStimIdx) = mean(vecTempReflDistMuAdja(i));
				vecReflDistGainOrth(intStimIdx) = mean(vecTempReflDistG(i));
				vecReflDistGainAdja(intStimIdx) = mean(vecTempReflDistGAdja(i));
			end
			
			%symmetry
			vecGainSymmetry = 100 - (vecReflDistGainOrth./vecReflDistGainAdja)*100;
			vecMeanSymmetry = 100 - (vecReflDistMuOrth./vecReflDistMuAdja)*100;
			
			% plot
			[h,pMeanGain]=ttest(vecSdProjMean,vecSdProjGain);
			[h,pMeanOrth]=ttest(vecSdProjMean,vecSdProjOrth);
			[h,pGainOrth]=ttest(vecSdProjGain,vecSdProjOrth);
			if boolMakeFigs
				figure;maxfig;
				subplot(2,3,1)
				hold on
				vecLoc = [0.2 0.5 0.8 1.1];
				dblJit = 0.05;
				scatter(vecLoc(1)+(rand(1,intStimNum)-0.5)*2*dblJit,vecSdProjMean,[200],[0.5 0.5 0.5],'.')
				scatter(vecLoc(2)+(rand(1,intStimNum)-0.5)*2*dblJit,vecSdProjGain1,[200],[0.5 0.5 0.5],'.')
				scatter(vecLoc(3)+(rand(1,intStimNum)-0.5)*2*dblJit,vecSdProjGain,[200],[0.5 0.5 0.5],'.')
				%scatter(vecLoc(3)+(rand(1,intStimNum)-0.5)*2*dblJit,vecSdProjRand,[200],[0.5 0.5 0.5],'.')
				%scatter(vecLoc(3)+(rand(1,intStimNum)-0.5)*2*dblJit,vecSdProjAdja,[200],[0.5 0.5 0.5],'.')
				scatter(vecLoc(4)+(rand(1,intStimNum)-0.5)*2*dblJit,vecSdProjOrth,[200],[0.5 0.5 0.5],'.')
				errorbar(vecLoc(1),mean(vecSdProjMean),std(vecSdProjMean)/sqrt(intStimNum),'bx','CapSize',20)
				errorbar(vecLoc(2),mean(vecSdProjGain1),std(vecSdProjGain1)/sqrt(intStimNum),'kx','CapSize',20)
				errorbar(vecLoc(3),mean(vecSdProjGain),std(vecSdProjGain)/sqrt(intStimNum),'kx','CapSize',20)
				%errorbar(vecLoc(3),mean(vecSdProjRand),std(vecSdProjRand)/sqrt(intStimNum),'rx','CapSize',20)
				%errorbar(vecLoc(3),mean(vecSdProjAdja),std(vecSdProjAdja)/sqrt(intStimNum),'kx','CapSize',20)
				errorbar(vecLoc(4),mean(vecSdProjOrth),std(vecSdProjOrth)/sqrt(intStimNum),'rx','CapSize',20)
				hold off;
				set(gca,'xtick',vecLoc,'xticklabel',{'Pop mean','Pop gain1','Pop gain','Stim-Orth'});
				xlabel('Projection axis');
				ylabel('Range of spiking rates (\sigmaHz)');
				title(sprintf('%s; MG,p=%.2e; MO,p=%.3f; GO,p=%.2e',strType,pMeanGain,pMeanOrth,pGainOrth));
				fixfig;
				ylim(gca,[0 max(get(gca,'ylim'))]);
				xlim([vecLoc(1)-0.2 vecLoc(end)+0.2]);
				
				%what is distance from origin to mean R compared with distance between ori Rs?
				subplot(2,3,2);
				hold on
				vecLoc = [0.2 0.5 0.8];
				dblJit = 0.05;
				scatter(vecLoc(1)+(rand(1,intStimNum)-0.5)*2*dblJit,vecDistFromOrigin,[200],[0.5 0.5 0.5],'.')
				errorbar(vecLoc(1),mean(vecDistFromOrigin),std(vecDistFromOrigin)/sqrt(intStimNum),'kx','CapSize',20)
				scatter(vecLoc(2)+(rand(1,intStimNum)-0.5)*2*dblJit,vecDistFromOrthOri,[200],[0.5 0.5 0.5],'.')
				errorbar(vecLoc(2),mean(vecDistFromOrthOri),std(vecDistFromOrthOri)/sqrt(intStimNum),'bx','CapSize',20)
				scatter(vecLoc(3)+(rand(1,intStimNum)-0.5)*2*dblJit,vecDistFromAdjaOri,[200],[0.5 0.5 0.5],'.')
				errorbar(vecLoc(3),mean(vecDistFromAdjaOri),std(vecDistFromAdjaOri)/sqrt(intStimNum),'rx','CapSize',20)
				hold off
				set(gca,'xtick',vecLoc,'xticklabel',{'Origin','Orth. ori','Adja. ori'});
				%ylim([0 60]);
				xlim([vecLoc(1)-0.2 vecLoc(end)+0.2]);
				xlabel('Distance of pop. stim resp to');
				ylabel('Pop. spiking rate distance (\DeltaHz)');
				title('Mu +/- sem, n=12 orientations');
				fixfig;
				
				%cos sim of gain with mean vs rand with mean
				%note that random vectors do not lead to cossim of 0 because the space is limited to
				%positive values only (Hz > 0)
				[h,pCosSim] = ttest2(vecAngleMeanAndGain,vecAngleMeanAndGainRand);
				subplot(2,3,3);
				hold on
				vecLoc = [0.2 0.8];
				dblJit = 0.05;
				scatter(vecLoc(1)+(rand(1,intStimNum)-0.5)*2*dblJit,vecAngleMeanAndGain,[200],[0.5 0.5 0.5],'.')
				errorbar(vecLoc(1),mean(vecAngleMeanAndGain),std(vecAngleMeanAndGain)/sqrt(intStimNum),'kx','CapSize',20)
				scatter(vecLoc(2)+(rand(1,intStimNum)-0.5)*2*dblJit,vecAngleMeanAndGainRand,[200],[0.5 0.5 0.5],'.')
				errorbar(vecLoc(2),mean(vecAngleMeanAndGainRand),std(vecAngleMeanAndGainRand)/sqrt(intStimNum),'bx','CapSize',20)
				hold off
				set(gca,'xtick',vecLoc,'xticklabel',{'Gain','Rand'});
				xlim([vecLoc(1)-0.2 vecLoc(end)+0.2]);
				xlabel('Neural axis');
				ylabel('Alignment with mean-axis (cos sim)');
				title(sprintf('T-test, cos-sim mean with gain vs rand: p=%.2e',pCosSim));
				fixfig;
				
				%symmetry around mean or gain axis
				vecLoc = [0.2 0.5 0.8 1.1];
				dblJit = 0.05;
				subplot(2,3,4);
				hold on
				errorbar(vecLoc(1),mean(vecReflDistMuOrth),std(vecReflDistMuOrth)/sqrt(intStimNum),'kx','CapSize',20)
				errorbar(vecLoc(2),mean(vecReflDistMuAdja),std(vecReflDistMuAdja)/sqrt(intStimNum),'bx','CapSize',20)
				errorbar(vecLoc(3),mean(vecReflDistGainOrth),std(vecReflDistGainOrth)/sqrt(intStimNum),'kx','CapSize',20)
				errorbar(vecLoc(4),mean(vecReflDistGainAdja),std(vecReflDistGainAdja)/sqrt(intStimNum),'bx','CapSize',20)
				hold off
				set(gca,'xtick',vecLoc,'xticklabel',{'Mu/Orth','Mu/Adj','Gain/Orth','Gain/Adj'});
				xlim([vecLoc(1)-0.2 vecLoc(end)+0.2]);
				xlabel('Reflection over/stim ori');
				ylabel('Prediction error after reflection (\DeltaHz)');
				fixfig;
				
				% symmetry around gain axis
				vecLoc = [0.2 0.8];
				subplot(2,3,5);
				[h,pSym]=ttest2(vecMeanSymmetry,vecGainSymmetry);
				hold on
				errorbar(0.2,mean(vecMeanSymmetry),std(vecMeanSymmetry)/sqrt(intStimNum),'bd','CapSize',20)
				errorbar(0.8,mean(vecGainSymmetry),std(vecGainSymmetry)/sqrt(intStimNum),'ko','CapSize',20)
				hold off
				ylabel('Manifold symmetry (%)');
				set(gca,'xtick',vecLoc,'xticklabel',{'Pop. mean','Pop. gain'});
				xlabel('Axis of symmetry');
				xlim([vecLoc(1)-0.2 vecLoc(end)+0.2]);
				title(sprintf('T-test, p=%.2e',pSym));
				fixfig;
				
				%make example figures of the above
				intS1 = 1;
				intS1Adja = 2;
				intS2 = 1+round(intStimNum/2);
				
				dblDistToGainCenter = mean(vecDistFromOrigin);
				dblAlignment=cossim(vecGainAx,ones(size(vecGainAx)));
				dblAng = acos(dblAlignment)/2+deg2rad(45);
				vecGainCenter = [cos(dblAng) sin(dblAng)]*dblDistToGainCenter;
				
				dblAng1 = dblAng + atan((mean(vecDistFromOrthOri(intS1)))/dblDistToGainCenter);
				dblAng1A = dblAng + atan((mean(vecDistFromOrthOri(intS1)) - mean(vecDistFromAdjaOri(intS1))/4)/dblDistToGainCenter);
				dblAng2 = dblAng - atan((mean(vecDistFromOrthOri(intS2)))/dblDistToGainCenter);
				
				vecCenter1 = [cos(dblAng1) sin(dblAng1)]*vecDistFromOrigin(intS1);
				vecCenter1A = [cos(dblAng1A) sin(dblAng1A)]*vecDistFromOrigin(intS1Adja);
				vecCenter2 = [cos(dblAng2) sin(dblAng2)]*vecDistFromOrigin(intS2);
				
				%reflect
				vecX = vecCenter1';
				[dummy,vecXprimeG]=getProjOnLine(vecX,vecGainCenter');
				vecXreflG = 2*vecXprimeG - vecX;
				[dummy,vecXprime]=getProjOnLine(vecX,ones(size(vecX)));
				vecXrefl = 2*vecXprime - vecX;
				
				subplot(2,3,6)%[3 4 7 8]);
				cla;
				hold on
				%mean
				dblMaxLim = ceil((max([vecCenter1 vecCenter2]) + max(vecSdProjGain))/5)*5;
				intOffset = 1;
				plot([0 dblMaxLim],[0 dblMaxLim],'-','color',[0.5 0.5 0.5])
				text(dblMaxLim*0.2+intOffset,dblMaxLim*0.2-intOffset,'Pop. mean axis','FontSize',16,'Color',[0.5 0.5 0.5],'Rotation',45)
				%gain
				plot([0 vecGainCenter(1)*1.4],[0 vecGainCenter(2)*1.4],'-','color','k')
				text(vecGainCenter(1)*0.3+intOffset,vecGainCenter(2)*0.3-intOffset,'Pop. gain axis','FontSize',16,'Color','k','Rotation',rad2deg(dblAng))
				
				%centers
				scatter(vecCenter1(1),vecCenter1(2),'rx');
				scatter(vecCenter2(1),vecCenter2(2),'bx');
				scatter(vecCenter1A(1),vecCenter1A(2),[],[0.5 0 0.5],'x');
				
				
				%covars
				ellipse(vecCenter1(1),vecCenter1(2),mean(vecSdProjGain(intS1)),mean(vecSdProjOrth(intS1)),dblAng1,'Color','r','LineStyle','-');
				ellipse(vecCenter2(1),vecCenter2(2),mean(vecSdProjGain(intS2)),mean(vecSdProjOrth(intS2)),dblAng2,'Color','b','LineStyle','-');
				
				%draw reflections
				plot([vecCenter1(1) vecXreflG(1)],[vecCenter1(2) vecXreflG(2)],'k--');
				plot([vecCenter1(1) vecXrefl(1)],[vecCenter1(2) vecXrefl(2)],'--','color',[0.5 0.5 0.5]);
				scatter(vecXreflG(1),vecXreflG(2),'ro');
				scatter(vecXrefl(1),vecXrefl(2),'rd');
				
				%text
				text(vecCenter1(1)+intOffset,vecCenter1(2)+intOffset,'Ref','FontSize',16,'Color','r')
				text(vecCenter2(1)+intOffset,vecCenter2(2)+intOffset,'Orth','FontSize',16,'Color','b')
				text(vecCenter1A(1)+intOffset,vecCenter1A(2)+intOffset,'Adja','FontSize',16,'Color',[0.5 0 0.5])
				text(vecXreflG(1)+intOffset,vecXreflG(2)-intOffset,'Gain-reflect','FontSize',16,'Color','r')
				text(vecXrefl(1)+intOffset,vecXrefl(2)+intOffset,'Mean-reflect','FontSize',16,'Color','r')
				
				%finish fig
				hold off
				axis equal;
				%xlim([0 dblMaxLim*(2/3)]);
				%ylim([0 dblMaxLim]);
				xlabel('Spiking rate axis X (Hz)')
				ylabel('Spiking rate axis Y (Hz)')
				fixfig;grid off;
				
				if boolSaveFigs
					%% save fig
					export_fig(fullpath(strFigurePath,sprintf('2B2O_Projections_%s.tif',strRec)));
					export_fig(fullpath(strFigurePath,sprintf('2B2O_Projections_%s.pdf',strRec)));
				end
			end
			
			%% predict gain from pupil
			%get pupil size per trial
			if 0%(isfield(sRec,'matPupilXY') && ~isempty(sRec.matPupilXY)) || (isfield(sRec,'vecRunningTime') && ~isempty(sRec.vecRunningTime))
				figure;maxfig;
			end
			if 0%isfield(sRec,'matPupilXY') && ~isempty(sRec.matPupilXY)
				sPupil = sRec.Pupil;
				matPupilXY = sRec.matPupilXY;
				vecPupilT = sRec.vecPupilT;
				vecTrialSpeed = sRec.vecPupilSize;
				
				
				vecPupilStimOn = structStim.vecPupilStimOn;
				vecPupilStimOff = structStim.vecPupilStimOff;
				vecPupilStimOn(indRem) = [];
				vecPupilStimOff(indRem) = [];
				
				vecTime = sPupil.vecTime;
				vecVals = sPupil.vecRadius;
				vecEvents = vecPupilStimOn;
				vecWindow  =[0 dblDur];
				vecPupMov = [0 sqrt(diff(sPupil.vecCenterX).^2 + diff(sPupil.vecCenterY).^2)];
				vecTrialSpeed = mean(getRespMat(sPupil.vecTime,sPupil.vecRadius,vecPupilStimOn-0.5,[0 dblDur]),2);
				vecPupilMov = mean(getRespMat(sPupil.vecTime,vecPupMov,vecPupilStimOn-0.5,[0 dblDur]),2);
				vecTrialGain = getProjOnLine(matMeanRate,vecGainAx);
				vecTrialMean = getProjOnLine(matMeanRate,ones(size(vecGainAx)));
				vecTrialMean2 = mean(matMeanRate,1);
				rGainPupilS = nancorr(vecTrialGain(:),vecTrialSpeed(:));
				rMeanPupilS = nancorr(vecTrialMean(:),vecTrialSpeed(:));
				indUseTrials = ~isnan(vecTrialSpeed) & ~isnan(vecTrialMean) & ~isnan(vecTrialGain);
				
				%run predictions
				vecTrialGain_hat = doCrossValidatedDecodingRR(zscore(vecTrialSpeed(indUseTrials)),zscore(vecTrialGain(indUseTrials)),1,1);
				vecTrialMean_hat = doCrossValidatedDecodingRR(zscore(vecTrialSpeed(indUseTrials)),zscore(vecTrialMean(indUseTrials)),1,1);
				
				[dblR2G,dblSS_totG,dblSS_resG,dblTG,dblPG,dblR2_adjustedG,dblR2_SEG] = getR2(zscore(vecTrialGain(indUseTrials)),vecTrialGain_hat,1);
				[dblR2M,dblSS_totM,dblSS_resM,dblTM,dblPM,dblR2_adjustedM,dblR2_SEM] = getR2(zscore(vecTrialMean(indUseTrials)),vecTrialMean_hat,1);
				
				subplot(2,3,1);
				scatter(zscore(vecTrialGain(indUseTrials)),vecTrialGain_hat,'b.');
				xlabel('Z-scored trial gain');
				ylabel('CV predicted norm trial gain');
				title('Gain prediction from pupil size');
				fixfig;
				
				subplot(2,3,2);
				scatter(zscore(vecTrialMean(indUseTrials)),vecTrialMean_hat,'k.');
				xlabel('Z-scored trial mean');
				ylabel('CV predicted norm trial mean');
				title('Mean prediction from pupil size');
				fixfig;
				
				subplot(2,3,3);
				hold on
				errorbar(0.2,dblR2G,dblR2_SEG,'bx');
				errorbar(0.8,dblR2M,dblR2_SEM,'kx');
				set(gca,'xtick',[0.2 0.8],'xticklabel',{'Gain','Mean'});
				ylabel('CV predictability (R^2)');
				xlim([0 1]);
				title(sprintf('T-tests vs 0: gain-p=%.2e,mean-p=%.2e',dblPG,dblPM));
				fixfig;
				
				if boolSaveFigs
					%% save fig
					export_fig(fullpath(strFigurePath,sprintf('2B3O_PupilPredGain_%s.tif',strRec)));
					export_fig(fullpath(strFigurePath,sprintf('2B3O_PupilPredGain_%s.pdf',strRec)));
				end
			end
			
			%% use locomotion
			if 0%isfield(sRec,'vecRunningTime') && ~isempty(sRec.vecRunningTime)
				vecRunningTime = sRec.vecRunningTime;
				vecRunningSpeed = sRec.vecRunningSpeed;
				
				
				vecTime = vecRunningTime;
				vecVals = vecRunningSpeed;
				vecEvents = vecStimOnTime;
				vecWindow  =[0 dblDur];
				
				vecTrialSpeed = mean(getRespMat(vecTime,vecVals,vecEvents,vecWindow),2);
				
				vecTrialGain = getProjOnLine(matMeanRate,vecGainAx);
				vecTrialMean = getProjOnLine(matMeanRate,ones(size(vecGainAx)));
				vecTrialMean2 = mean(matMeanRate,1);
				rGainPupilS = nancorr(vecTrialGain(:),vecTrialSpeed(:));
				rMeanPupilS = nancorr(vecTrialMean(:),vecTrialSpeed(:));
				indUseTrials = ~isnan(vecTrialSpeed) & ~isnan(vecTrialMean) & ~isnan(vecTrialGain);
				
				%run predictions
				matX = vecTrialSpeed(indUseTrials);
				dblMaxMean = max(abs(mean(matX) ./ range(matX)));
				
				vecTrialGain_hat = doCrossValidatedDecodingRR(zscore(vecTrialSpeed(indUseTrials)),zscore(vecTrialGain(indUseTrials)),1,1);
				vecTrialMean_hat = doCrossValidatedDecodingRR(zscore(vecTrialSpeed(indUseTrials)),zscore(vecTrialMean(indUseTrials)),1,1);
				
				[dblR2G,dblSS_totG,dblSS_resG,dblTG,dblPG,dblR2_adjustedG,dblR2_SEG] = getR2(zscore(vecTrialGain(indUseTrials)),vecTrialGain_hat,1);
				[dblR2M,dblSS_totM,dblSS_resM,dblTM,dblPM,dblR2_adjustedM,dblR2_SEM] = getR2(zscore(vecTrialMean(indUseTrials)),vecTrialMean_hat,1);
				
				subplot(2,3,1);
				scatter(zscore(vecTrialGain(indUseTrials)),vecTrialGain_hat,'b.');
				xlabel('Z-scored trial gain');
				ylabel('CV predicted norm trial gain');
				title('Gain prediction from pupil size');
				fixfig;
				
				subplot(2,3,2);
				scatter(zscore(vecTrialMean(indUseTrials)),vecTrialMean_hat,'k.');
				xlabel('Z-scored trial mean');
				ylabel('CV predicted norm trial mean');
				title('Mean prediction from pupil size');
				fixfig;
				
				subplot(2,3,3);
				hold on
				errorbar(0.2,dblR2G,dblR2_SEG,'bx');
				errorbar(0.8,dblR2M,dblR2_SEM,'kx');
				set(gca,'xtick',[0.2 0.8],'xticklabel',{'Gain','Mean'});
				ylabel('CV predictability (R^2)');
				xlim([0 1]);
				title(sprintf('T-tests vs 0: gain-p=%.2e,mean-p=%.2e',dblPG,dblPM));
				fixfig;
				
			end
			
			%% save data
			if boolSaveData
				%save data
				%% predictions
				sPrediction = struct;
				sPrediction.matMeanRate = matMeanRate;
				sPrediction.matTunePredCV = matTunePredCV;
				sPrediction.matMeanPredCV = matMeanPredCV;
				sPrediction.matGainPredCV = matGainPredCV;
				sPrediction.matGain1PredCV = matGain1PredCV;
				
				sPrediction.dblR2_Tune = dblR2_Tune;
				sPrediction.dblR2_Mean = dblR2_Mean;
				sPrediction.dblR2_Gain = dblR2_Gain;
				sPrediction.dblR2_Gain1 = dblR2_Gain1;
				
				sPrediction.vecR2_Tune = vecR2_Tune;
				sPrediction.vecR2_Mean = vecR2_Mean;
				sPrediction.vecR2_Gain = vecR2_Gain;
				sPrediction.vecR2_Gain1 = vecR2_Gain1;
				
				sPrediction.dblPredPerNeuronTune = dblPredPerNeuronTune;
				sPrediction.dblPredPerNeuronMean = dblPredPerNeuronMean;
				sPrediction.dblPredPerNeuronGain = dblPredPerNeuronGain;
				sPrediction.dblPredPerNeuronGain1 = dblPredPerNeuronGain1;
				
				
				%% projections
				sProjection = struct;
				sProjection.vecReflDistMuOrth = vecReflDistMuOrth;
				sProjection.vecReflDistMuAdja = vecReflDistMuAdja;
				sProjection.vecReflDistGainOrth = vecReflDistGainOrth;
				sProjection.vecReflDistGainAdja = vecReflDistGainAdja;
				
				sProjection.vecGainSymmetry = vecGainSymmetry;
				sProjection.vecMeanSymmetry = vecMeanSymmetry;
				
				sProjection.vecAngleMeanAndGain = vecAngleMeanAndGain;
				sProjection.vecAngleMeanAndGainRand = vecAngleMeanAndGainRand;
				
				sProjection.vecSdProjMean = vecSdProjMean;
				sProjection.vecSdProjGain1 = vecSdProjGain1;
				sProjection.vecSdProjGain = vecSdProjGain;
				sProjection.vecSdProjRand = vecSdProjRand;
				sProjection.vecSdProjOrth = vecSdProjOrth;
				sProjection.vecSdProjAdja = vecSdProjAdja;
				
				sProjection.vecDistFromOrigin = vecDistFromOrigin;
				sProjection.vecDistFromOrthOri = vecDistFromOrthOri;
				sProjection.vecDistFromAdjaOri = vecDistFromAdjaOri;
				
				save([strTargetDataPath 'TimeCodingAggQC3ABI_' strRunStim strRec '_' strRunArea '.mat'],'strRec','strRunArea','sPrediction','sProjection','dblBC','dblMaxDevFrac');
			end
		end
	end
end
toc
