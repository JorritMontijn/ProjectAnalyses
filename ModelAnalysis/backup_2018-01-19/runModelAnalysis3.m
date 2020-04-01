%% load data
%close all
if ~exist('sData','var')
	strSimulation = 'ScaledUp2016-12-11';
	load(['D:\Data\Results\V1_LIFmodel\Simulation_' strSimulation '.mat']);
end
%% unpack
vecOverallT = sData.vecOverallT;
dblDeltaT = sData.dblDeltaT;
matCortConn = sData.matCortConn;
dblSynSpikeMem = sData.dblSynSpikeMem;
vecCortSynType = sData.vecCortSynType;
intCortexCells = sData.intCortexCells;
vecCortDelay = sData.vecCortDelay;
vecCortConductance = sData.vecCortConductance;
vecCellThresh = sData.vecCellThresh;
vecTauPeakByType = sData.vecTauPeakByType;
vecCellV_E = sData.vecCellV_E;
vecCellV_I = sData.vecCellV_I;
vecCellV_AHP = sData.vecCellV_AHP;
vecCellV_Leak = sData.vecCellV_Leak;
vecCellCm = sData.vecCellCm;
vecCellG_Leak = sData.vecCellG_Leak;
vecCellG_AHP = sData.vecCellG_AHP;
vecSynConductanceON_to_Cort = sData.vecSynConductanceON_to_Cort;
vecSynConductanceOFF_to_Cort = sData.vecSynConductanceOFF_to_Cort;
vecSynWeightON_to_Cort = sData.vecSynWeightON_to_Cort;
vecSynWeightOFF_to_Cort = sData.vecSynWeightOFF_to_Cort;
vecSynDelayON_to_Cort = sData.vecSynDelayON_to_Cort;
vecSynDelayOFF_to_Cort = sData.vecSynDelayOFF_to_Cort;
matSynConnON_to_Cort = sData.matSynConnON_to_Cort;
matSynConnOFF_to_Cort = sData.matSynConnOFF_to_Cort;
matBlankLGN_ON = sData.matBlankLGN_ON;
matBlankLGN_OFF = sData.matBlankLGN_OFF;
cellLGN_ON = sData.cellLGN_ON;
cellLGN_OFF = sData.cellLGN_OFF;
vecTrialOris = sData.vecTrialOris;
vecTrialOriIdx = sData.vecTrialOriIdx;
vecStimStartSecs = sData.vecStimStartSecs;
vecTrialEndSecs = sData.vecTrialEndSecs;
vecThisV = sData.vecThisV;
boolStimPresent = sData.boolStimPresent;
intPrevTrial = sData.intPrevTrial;
intTrialT = sData.intTrialT;
intIter = sData.intIter;
cellSpikeTimesLGN_ON = sData.cellSpikeTimesLGN_ON;
cellSpikeTimesLGN_OFF = sData.cellSpikeTimesLGN_OFF;
cellSpikeTimesCortex = sData.cellSpikeTimesCortex;
vecSpikeCounterLGN_ON = sData.vecSpikeCounterLGN_ON;
vecSpikeCounterLGN_OFF = sData.vecSpikeCounterLGN_OFF;
vecSpikeCounterCortex = sData.vecSpikeCounterCortex;
intPreAllocationSize = sData.intPreAllocationSize;
%end
%% get data
%load('Simulation2016-10-04.mat')
vecStimStartSecs = vecStimStartSecs+0.25;
dblStimDur = 0.5;
vecStimStopSecs = vecStimStartSecs+dblStimDur;
%vecTrialEndSecs = vecStimStartSecs+dblTrialDur;
intColumns = 252; %252
vecCellTypesPerColumn = [1 1 1 1 2]; %[1=pyramid 2=interneuron]
vecCellTypes = repmat(vecCellTypesPerColumn',[intColumns 1]);

%% get spiking data
cellSpikeTimesCortex = cat(1,cellSpikeTimesLGN_ON,cellSpikeTimesLGN_OFF);
vecCellSize = size(cellSpikeTimesCortex);
intTrials = numel(vecTrialOris);
intNeurons = numel(cellSpikeTimesCortex);
if ~exist('matModelResp','var')
	matModelResp = nan(intNeurons,intTrials);
	parfor intTrial=1:intTrials
		dblStartT = vecStimStartSecs(intTrial);
		dblStopT = vecStimStopSecs(intTrial);
		
		matModelResp(:,intTrial) = cellfun(@sum,...
			cellfun(@and,cellfun(@gt,cellSpikeTimesCortex,cellfill(dblStartT,vecCellSize),'UniformOutput',false),...
			cellfun(@lt,cellSpikeTimesCortex,cellfill(dblStopT,vecCellSize),'UniformOutput',false),...
			'UniformOutput',false));
		if mod(intTrial,100) == 0
			fprintf('Now at trial %d/%d [%s]\n',intTrial,intTrials,getTime);
		end
	end
end

%% build structStim
structStim.Orientation = vecTrialOris;
structStim.TrialNumber = 1:length(vecTrialOris);
%structStim.FrameOn = vecStimStartSecs;
%structStim.FrameOff = vecStimStopSecs;
vecOrientations = unique(vecTrialOris);
vecOriDegs = rad2deg(vecOrientations);
sTypes = getStimulusTypes(structStim);
cellSelect = getSelectionVectors(structStim,sTypes);
vecTrialOriIdx = label2idx(vecTrialOris);
vecTrialOriIndex = vecTrialOriIdx;


%%
%{
matData = matModelResp(1:20,:);
[dblPerformanceCV,vecDecodedIndexCV,matPosteriorProbabilityCV,matWeights] = doCrossValidatedDecodingLR(matData,vecTrialTypes,0.1);
dblPerformanceCV

vecW1 = matWeights(:,1);
vecdMu12 = [xmean(matData(:,vecTrialOriIdx==1),2); 1] - [xmean(matData(:,vecTrialOriIdx==2),2); 1];

vecdMu12'*matWeights(:,2)

matDataPlusLin = [matData; ones(1,size(matData,2))];
matActivation = matWeights'*matDataPlusLin;
	
vecDeltaMu

matCov1 = cov(matData(:,vecTrialOriIdx==1)');
%}
%% reproduce fig. 4 from Series et al.
%matModelRespP = matModelResp(vecCellTypes==1,:);
matModelRespP = matModelResp;
intNeurons = size(matModelRespP,1);
vecGroupSizes = [2.^[0:9] intNeurons];
%vecGroupSizes = 2.^[0:8];
dblLambda = 1/100;

intGroups = numel(vecGroupSizes);
vecI_Cross = nan(1,intGroups);
vecI_Same = nan(1,intGroups);
vecI_shuff_Cross = nan(1,intGroups);
vecI_shuff_Same = nan(1,intGroups);
vecI_CV = nan(1,intGroups);
vecI_CV_shuff = nan(1,intGroups);
vecI_shuff_bc = nan(1,intGroups);

vecPerf = nan(1,intGroups);
vecPerfPred = nan(1,intGroups);

vecPerf_shuff = nan(1,intGroups);
vecPerfPred_shuff = nan(1,intGroups);

vecPerf_shuff_bc = nan(1,intGroups);
vecPerfPred_shuff_bc = nan(1,intGroups);
for intGroupSizeIdx=1:intGroups;
	intGroupSize = vecGroupSizes(intGroupSizeIdx);
	vecNeurons = linspace(1,intNeurons,intGroupSize+1);
	vecNeurons = ceil(vecNeurons(1:(end-1)));
	%vecNeurons = 1:intGroupSize;
	
	%select only classes 1 and 2
	matData = matModelRespP(vecNeurons,:);
	indClasses12 = (vecTrialOriIdx==1 | vecTrialOriIdx==2);
	matData12 = matData(:,indClasses12);
	intTrials = size(matData12,2);
	indFirstHalf = true(1,intTrials);
	indFirstHalf((intTrials/2)+1:end) = false;
	indSecondHalf = ~indFirstHalf;
	vecClasses01 = vecTrialOriIdx(indClasses12);
	vecClasses01(vecClasses01==2)=0;
	
	%first half
	matData01FH = matData12(:,indFirstHalf);
	vecClasses01FH = vecClasses01(indFirstHalf);
	
	%second half
	matData01SH = matData12(:,indSecondHalf);
	vecClasses01SH = vecClasses01(indSecondHalf);

	%check if some need to be removed
	vecVar01_FH = ([xvar(matData01FH(:,vecClasses01FH==1),2); 1] + [xvar(matData01FH(:,vecClasses01FH==0),2); 1])/2;
	vecVar01_SH = ([xvar(matData01SH(:,vecClasses01SH==1),2); 1] + [xvar(matData01SH(:,vecClasses01SH==0),2); 1])/2;
	indRem = (vecVar01_SH==0 | vecVar01_FH==0);
	matData01FH(indRem,:) = [];
	matData01SH(indRem,:) = [];
	matData12(indRem,:) = [];
	
	%get correction factors
	dblDiffTheta = range(vecClasses01FH);
	dblSubFac =(2*intNeurons)/(intTrials*(dblDiffTheta.^2));
	dblProdFacRaw = ((2*intTrials-intNeurons-3)/(2*intTrials-2));
	
	% get weights for binary logistic regression
	[vecWeights, vecLLH] = logitBin(matData01FH, vecClasses01FH, dblLambda);
	matDPL = [matData01FH;ones(1,size(matData01FH,2))];
	y = ~round(exp(-log1pexp(vecWeights'*matDPL)));
	dblPerformanceSame = sum(y==vecClasses01FH)/numel(vecClasses01FH);
	
	%cross
	matDPL = [matData01SH;ones(1,size(matData01SH,2))];
	vecOutputA = vecWeights'*matDPL;
	vecClass0 = vecOutputA(vecClasses01SH==0);
	vecClass1 = vecOutputA(vecClasses01SH==1);
	getdprime2(vecClass0,vecClass1)
	
	
	y = ~round(exp(-log1pexp(vecOutputA)));
	dblPerformanceCross = sum(y==vecClasses01SH)/numel(vecClasses01SH);
	
	%get information for same-half
	vecdMu01_FH = [xmean(matData01FH(:,vecClasses01FH==1),2); 1] - [xmean(matData01FH(:,vecClasses01FH==0),2); 1];
	vecVar01_FH = ([xvar(matData01FH(:,vecClasses01FH==1),2); 1] + [xvar(matData01FH(:,vecClasses01FH==0),2); 1])/2;
	dblI_Same = vecdMu01_FH'*vecWeights;
	
	%get information for unseen half
	vecdMu01_SH = [xmean(matData01SH(:,vecClasses01SH==1),2); 1] - [xmean(matData01SH(:,vecClasses01SH==0),2); 1];
	dblI_Cross = vecdMu01_SH'*vecWeights;
	dblPredA_Cross = nansum((1/2)*erf((sqrt([dblI_Cross dblI_Cross])/2)/sqrt(2)));
	dblPredA_Cross = dblPredA_Cross + (1-dblPredA_Cross)*(1/2);
	
	%calculate directly
	[dblPredA_same,matPredA,dblD2_same,dblD2mat,dblD2_diag] = getSeparation(matData01FH,vecClasses01FH,0);
	dblI_CV_unbiased = dblD2_same*dblProdFacRaw-dblSubFac;
	
	%save
	vecI_CV(intGroupSizeIdx) = dblI_CV_unbiased;
	vecI_Same(intGroupSizeIdx) = dblI_Same;
	vecI_Cross(intGroupSizeIdx) = dblI_Cross;
	vecPerf(intGroupSizeIdx) = dblPerformanceCross;
	vecPerfPred(intGroupSizeIdx) = dblPredA_Cross;
	
	%% shuffled
	%get information for shuffled
	dblI_shuff = sum(((vecdMu01_FH.^2 / dblDiffTheta) ./ vecVar01_FH));
	dblI_shuff_unbiased = ...
		sum(((vecdMu01_FH.^2 / dblDiffTheta) ./ vecVar01_FH)) ...
		* ((intTrials - 2) / (intTrials - 1))...
		- dblSubFac;
	
	if dblI_shuff_unbiased > eps
		dblPredA_shuff_unbiased = nansum((1/2)*erf((sqrt([dblI_shuff_unbiased dblI_shuff_unbiased])/2)/sqrt(2)));
		dblPredA_shuff_unbiased = dblPredA_shuff_unbiased + (1-dblPredA_shuff_unbiased)*(1/2);
	else
		dblPredA_shuff_unbiased = 0.5;
	end
	
	%get shuffled information directly v2
	matShuffledData12 = nan(size(matData12));
	for intClass = [0 1];
		vecThisClass = find(vecClasses01==intClass);
		for intThisNeuron=1:size(matShuffledData12,1)
			matShuffledData12(intThisNeuron,vecThisClass) = matData12(intThisNeuron,circshift(vecThisClass,intThisNeuron-1,2));
		end
	end
	
	%first half
	matData01FH_shuff = matShuffledData12(:,indFirstHalf);
	
	%second half
	matData01SH_shuff = matShuffledData12(:,indSecondHalf);
	
	%calculate directly
	[dblPredA_shuff,matPredA_shuff,dblD2_shuff] = getSeparation(matData01FH_shuff,vecClasses01FH,0);
	dblI_CV_shuff_unbiased = dblD2_shuff*dblProdFacRaw-dblSubFac;
	
	%get decoding for shuffled
	[vecWeights_shuff, vecLLH_shuff] = logitBin(matData01FH_shuff, vecClasses01FH, dblLambda);
	matDPL_shuff = [matData01FH_shuff;ones(1,size(matData01FH_shuff,2))];
	y = ~round(exp(-log1pexp(vecWeights_shuff'*matDPL_shuff)));
	dblPerformanceSame_shuff = sum(y==vecClasses01FH)/numel(vecClasses01FH);
	
	%get decoding for shuffled
	matDPL_shuff = [matData01SH_shuff;ones(1,size(matData01SH_shuff,2))];
	y = ~round(exp(-log1pexp(vecWeights_shuff'*matDPL_shuff)));
	dblPerformanceCross_shuff = sum(y==vecClasses01SH)/numel(vecClasses01SH);
	
	%get information for same-half
	dblI_Same_shuff = vecdMu01_FH'*vecWeights_shuff;
	
	%get information for unseen half
	dblI_Cross_shuff = vecdMu01_SH'*vecWeights_shuff;
	
	
	%save
	vecI_CV_shuff(intGroupSizeIdx) = dblI_CV_shuff_unbiased;
	vecI_shuff_Cross(intGroupSizeIdx) = dblI_Cross_shuff;
	vecI_shuff_bc(intGroupSizeIdx) = dblI_shuff_unbiased;
	vecPerf_shuff(intGroupSizeIdx) = dblPerformanceCross_shuff;
	vecPerfPred_shuff(intGroupSizeIdx) = dblPredA_shuff;
	vecPerfPred_shuff_bc(intGroupSizeIdx) = dblPredA_shuff_unbiased;
	
	%% diag: train on shuffled, test on original
	%get decoding for shuffled
	[vecWeights_shuff, vecLLH_shuff] = logitBin(matData01FH_shuff, vecClasses01FH, dblLambda);
	matDPL_shuff = [matData01FH;ones(1,size(matData01FH,2))];
	y = ~round(exp(-log1pexp(vecWeights_shuff'*matDPL_shuff)));
	dblPerformanceSame_diag = sum(y==vecClasses01FH)/numel(vecClasses01FH);
	
	
	%get decoding for shuffled
	[vecWeights_shuff, vecLLH_shuff] = logitBin(matData01FH_shuff, vecClasses01FH, dblLambda);
	matDPL_shuff = [matData01SH;ones(1,size(matData01SH,2))];
	y = ~round(exp(-log1pexp(vecWeights_shuff'*matDPL_shuff)));
	dblPerformanceCross_diag = sum(y==vecClasses01SH)/numel(vecClasses01SH);
	
end

% plot
% information
%logistic regression non-shuffled, CV 50/50 (cross)
%vecI

%logistic regression shuffled, CV 50/50 (cross)
%vecI_shuff

%bias corrected
%vecI_shuff_bc %direct
figure
subplot(2,2,1)
plot(vecGroupSizes,vecI_Cross,'xb-');
hold on
plot(vecGroupSizes,vecI_shuff_Cross,'xr-');
hold off
xlabel('Sample size')
ylabel('Fisher information (d''^2)')
title('Log reg weight based, blue=raw, red=shuffled')
fixfig

subplot(2,2,2)
plot(vecGroupSizes,vecI_CV,'xb-');
hold on
plot(vecGroupSizes,vecI_shuff_bc,'xm-');
plot(vecGroupSizes,vecI_CV_shuff,'xr-');
hold off
xlabel('Sample size')
ylabel('Fisher information (d''^2)')
title('Direct; blue=I(CV), magenta=I_s_h_u_f_f(mu/var), red=I_s_h_u_f_f(CV)')
fixfig

% decoding accuracy
%vecPerf
%vecPerfPred


%vecPerf_shuff
%vecPerfPred_shuff
%vecPerfPred_shuff_bc

subplot(2,2,3)
plot(vecGroupSizes,vecPerf,'xb-');
hold on
plot(vecGroupSizes,vecPerfPred,'xb--');
hold off
xlabel('Sample size')
ylabel('Decoding accuracy')
set(gca,'xscale','log')
title('Solid=real, dashed=predicted')
fixfig

subplot(2,2,4)
plot(vecGroupSizes,vecPerf_shuff,'xr-');
hold on
plot(vecGroupSizes,vecPerfPred_shuff,'xr--');
hold off
xlabel('Sample size')
ylabel('Decoding accuracy')
set(gca,'xscale','log')
title('Solid=real, dashed=predicted')
fixfig
%{
subplot(2,3,6)
plot(vecGroupSizes,vecPerfPred_shuff_bc,'xm--');
hold off
xlabel('Sample size')
ylabel('Decoding accuracy')
set(gca,'xscale','log')
title('Solid=real, dashed=predicted')
fixfig
%}
%%
return
%save figure
export_fig(['D:\Data\Results\V1_LIFmodel\figures\' strSimulation '_InformationAndDecoding_LGN_Lambda1_100.tif']);
export_fig(['D:\Data\Results\V1_LIFmodel\figures\' strSimulation '_InformationAndDecoding_LGN_Lambda1_100.pdf']);


%%
figure
subplot(2,2,1)
plot(vecGroupSizes,vecI_Cross,'xb-');
hold on
plot(vecGroupSizes,vecI_shuff_Cross,'xr-');
hold off
xlabel('Sample size')
ylabel('Fisher information (d''^2)')
title('Blue=raw, red=shuffled')
fixfig

subplot(2,2,2)
plot(vecGroupSizes,vecPerf,'xb-');
hold on
plot(vecGroupSizes,vecPerf_shuff,'xr-');
hold off
xlabel('Sample size')
ylabel('Decoding accuracy')
set(gca,'xscale','log')
title('Blue=raw, red=shuffled')
fixfig


%% how to calculate fisher information; PLOS Comp Bio (2015)
%I = train and test on raw data (CV, 50/50)

%I_shuffle = CV (50/50) of shuffled data
%I_shuffle = sum(F_prime(intNeuron).^2 / var(intNeuron)) across neurons
%I_shuffle,unbiased = SUM[{(dMu(intNeuron).^2 / d(Theta)) / var(intNeuron)} * {(T - 2) /
%(T - 1)}] - {(2 * N) / (T * d(Theta).^2)}
% Here, T=intTrials, N=intNeurons

%I_diag = train on shuffled, test on original (CV 100/100 => use 50/50)




return
%% get relationship predicted/actual
intIters=10;
vecSizes = [2.^(0:8)];
intSizes = numel(vecSizes);
matDprimeSquared = nan(intIters,intSizes);
matPredictedAcc = nan(intIters,intSizes);
matActualAcc = nan(intIters,intSizes);
intNeuronCounter = 0;
indOriPair = vecTrialOriIdx==1 | vecTrialOriIdx==2;
matModelRespPair = matModelResp(:,indOriPair);
vecOriPair = vecTrialOriIdx(indOriPair);
for intSize=1:intSizes
	intGroupSize = vecSizes(intSize);
	fprintf('Running group size %d/%d (%d) [%s]\n',intSize,intSizes,intGroupSize,getTime);
	parfor intIter=1:intIters
		vecNeurons = randperm(intNeurons,intGroupSize);
		
		[dblPredA,matPredA,dblDprimeSquaredOffDiagonal] = getSeparation(matModelResp(vecNeurons,:),vecTrialOriIdx);
		matDprimeSquared(intIter,intSize) = dblDprimeSquaredOffDiagonal;
		matPredictedAcc(intIter,intSize) = dblPredA(end,1);
		matActualAcc(intIter,intSize) = doCrossValidatedDecodingLR(matModelResp(vecNeurons,:),vecTrialOriIdx,0.01);
	end
end

%

plot([0 1],[0 1],'k--')
hold on
scatter(mean(matActualAcc,1),mean(matPredictedAcc,1))
hold off
%xlim([0.5 1])
%ylim([0.5 1])

%% plot
clf
subplot(2,2,1);
plot(sort(vecDprimeSquared));
xlabel('Number of neurons');
ylabel('D''^2');
title('Greedy classifier');
xlim([0 intNeurons]);

subplot(2,2,2);
plot(sort(vecPredictedAcc));
xlabel('Number of neurons');
ylabel('Predicted decoding accuracy');
xlim([0 intNeurons]);
ylim([0 1]);
title('Decoding accuracy predicted from d''^2' )

subplot(2,2,3);
plot(vecPerformance);
xlabel('Number of neurons');
ylabel('Real decoding accuracy');
xlim([0 intNeurons]);
ylim([0 1]);
title('Real decoding accuracy, neurons sorted by d''^2 contribution' )

subplot(2,2,4)
title('Model')
axes off
return

export_fig(['D:\Data\Results\V1_LIFmodel\' strSimulation '_greedyDecoding2.tif']);
export_fig(['D:\Data\Results\V1_LIFmodel\' strSimulation '_greedyDecoding2.pdf']);



%% get neuron group
vecUseNeurons = [1 2 3];
matDataPoints = matModelResp(vecUseNeurons,:)';

%calc orth/para
[vecSD_Orth,vecSD_Para] = calcOrthPara(matDataPoints,vecTrialOriIndex);
[vecShuffSD_Orth,vecShuffSD_Para] = calcOrthPara(matDataPoints,vecTrialOriIndex,true);
clc
dblOrthRel = mean(vecSD_Orth)/mean(vecShuffSD_Orth)
dblParaRel = mean(vecSD_Para)/mean(vecShuffSD_Para)

%transform matrix
intStimTypes = numel(vecOrientations);
intReps = intTrials/intStimTypes;
matStimResp = nan(intNeurons,intStimTypes,intReps);
for intStimType=1:intStimTypes
	matStimResp(:,intStimType,:) = matModelResp(:,vecTrialOriIdx==intStimType);
end

mat3Resp = matStimResp(vecUseNeurons,:,:);
plotTube(mat3Resp)

%% investigate orientation coding per spike

