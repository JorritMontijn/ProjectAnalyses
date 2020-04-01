

%% get simulation name [strSimulation] from [intLoadSim]
clear all;close all;
strSimulation = 'Simulation_xAreaDistributed_C17_2019-07-12';
strSimDataPath = 'D:\Data\SimAggregates\';
load([strSimDataPath strSimulation '.mat']);
strFigPath = 'D:\Data\ResultsModelXarea\';

%% RUN: #header
dblSimStep = median(diff(vecOverallT));
if strcmp(class(cellSpikeTimesCortex{1}),'int32')
	cellSpikeTimesCortex = doSpikeIntToDouble(cellSpikeTimesCortex,dblSimStep,vecOverallT(1)-dblSimStep);
end

cellIn = strsplit(strSimulation,'_');
strFiller = cellIn{1};
strType = cell2mat(cellIn(2:(end-1)));
strDate = cellIn{end};
strTag = [strType '_' strDate];

%% split population
intPopSizeV1 = 1200;
vecAllNeuronsV1 = find(vecCellArea==1);
vecSelectNeuronsV1 = sort(vecAllNeuronsV1(randperm(numel(vecAllNeuronsV1),intPopSizeV1)));
intPopSizeV2 = 1200;
vecAllNeuronsV2 = find(vecCellArea==2);
vecSelectNeuronsV2 = sort(vecAllNeuronsV2(randperm(numel(vecAllNeuronsV2),intPopSizeV2)));
vecUseNeurons = cat(2,vecSelectNeuronsV1,vecSelectNeuronsV2);
vecCellAreaV1V2 = cat(2,ones(size(vecSelectNeuronsV1)),2*ones(size(vecSelectNeuronsV2)));
cellSpikesV1V2 = cellSpikeTimesCortex(vecUseNeurons);

%% get stimulus types
vecUniqueStimTypes = unique(vecTrialStimType);
intStimTypes = numel(vecUniqueStimTypes);
intTrials = numel(vecStimStartSecs);

 %% run trial-based analysis: is spiking rate and classification performance in V2 associated with V1 synchronization?
%build trial by neuron matrix
vecTrialEdges = zeros(1,size(vecTrialStartSecs,2)*2);
vecTrialEdges(1:2:end) = vecStimStartSecs;
vecTrialEdges(2:2:end) = vecStimStopSecs;
vecTrialEdges(end+1) = vecTrialEdges(end) + median(vecStimStartSecs(2:end) - vecStimStopSecs(1:(end-1)));
dblITIDur = median(vecStimStartSecs(2:end) - vecStimStopSecs(1:(end-1)));
dblTrialDur = median(vecStimStopSecs - vecStimStartSecs);
%get V1
vecCellsV1 = find(vecCellAreaV1V2==1);
intNeuronsV1 = numel(vecCellsV1);
matRespV1 = zeros(intTrials,intNeuronsV1);
for intNeuronIdx=1:intNeuronsV1
	intCellV1 = vecCellsV1(intNeuronIdx);
	vecSpikes = cellSpikesV1V2{intCellV1};
	vecTrialResp = histcounts(vecSpikes,vecTrialEdges);
	vecTrialResp = vecTrialResp(1:2:end);
	matRespV1(:,intNeuronIdx) = vecTrialResp;
end
%get V2
vecCellsV2 = find(vecCellAreaV1V2==2);
intNeuronsV2 = numel(vecCellsV2);
matRespV2 = zeros(intTrials,intNeuronsV2);
for intNeuronIdx=1:intNeuronsV2
	intCellV2 = vecCellsV2(intNeuronIdx);
	vecSpikes = cellSpikesV1V2{intCellV2};
	vecTrialResp = histcounts(vecSpikes,vecTrialEdges);
	vecTrialResp = vecTrialResp(1:2:end);
	matRespV2(:,intNeuronIdx) = vecTrialResp;
end

%% go through contrasts and compare with 0
[vecTrialTypes,vecContrasts,vecCounts,cellSelect,vecRepetition] = label2idx(vecTrialContrasts);
vecH = zeros(1,numel(vecContrasts)-1);
vecM = zeros(1,numel(vecContrasts)-1);
vecFA = zeros(1,numel(vecContrasts)-1);
vecCR = zeros(1,numel(vecContrasts)-1);
vecHV1 = zeros(1,numel(vecContrasts)-1);
for intC=find(vecContrasts>0)
	intC
	vecUseContrasts = [0 vecContrasts(intC)];
	vecUseIdxC = find(ismember(vecContrasts,vecUseContrasts));
	indUseTrialC = ismember(vecTrialTypes,vecUseIdxC);
	vecUseTrialTypes = vecTrialTypes(indUseTrialC);
	dblLambda = 1;
	intTypeCV = 2;
	
	%V1
	matUseResp = matRespV1(indUseTrialC,:);
	[dblPerformanceV1,vecDecodedIndexV1,matPosteriorProbabilityV1,matWeightsV1,dblMeanErrorDegs,matConfusionV1] = doCrossValidatedDecodingLR(matUseResp,vecUseTrialTypes,intTypeCV,dblLambda);
	vecHV1(intC-1) = matConfusionV1(2,2);
	
	%V2
	matUseResp = matRespV2(indUseTrialC,:);
	[dblPerformanceV2,vecDecodedIndexV2,matPosteriorProbabilityV2,matWeightsV2,dblMeanErrorDegs,matConfusionV2] = doCrossValidatedDecodingLR(matUseResp,vecUseTrialTypes,intTypeCV,dblLambda);
	
	%assign
	vecH(intC-1) = matConfusionV2(2,2);
	vecCR(intC-1) = matConfusionV2(1,1);
	vecM(intC-1) = matConfusionV2(1,2);
	vecFA(intC-1) = matConfusionV2(2,1);
end

%plot
[vecPV1,vecCIV1] = binofit(vecHV1,vecCounts(2:end));
[vecPV2,vecCIV2] = binofit(vecH,vecCounts(2:end));
figure;
hold on
plot([vecContrasts(2) vecContrasts(end)],[0.5 0.5],'k--')
errorbar(vecContrasts(2:end),vecPV1,vecCIV1(:,1)-vecPV1',vecCIV1(:,2)-vecPV1','b');
errorbar(vecContrasts(2:end),vecPV2,vecCIV2(:,1)-vecPV2',vecCIV2(:,2)-vecPV2','r');
hold off
set(gca,'xscale','log')
xlabel('Visual contrast (%)');
ylabel('Stimulus detection V2 (hit rate)')
title('LR decoder V1 (blue) and V2 (red) population')
fixfig

%% decode one contrast in V1 & V2
[dummy,intUseC] = min(abs(vecPV2-0.75));
intUseC = intUseC + 1;
dblUseContrast = vecContrasts(intUseC);
vecUseContrasts = [0 dblUseContrast];
vecUseIdxC = find(ismember(vecContrasts,vecUseContrasts));
indUseTrialC = ismember(vecTrialTypes,vecUseIdxC);
vecUseTrialTypes = vecTrialTypes(indUseTrialC);
dblLambda = 1;
intTypeCV = 2;
%decode V1
matUseRespV1 = matRespV1(indUseTrialC,:);
[dblPerformance,vecDecodedIndex,matPosteriorProbability,matWeightsV1,dblMeanErrorDegs,matConfusion] = doCrossValidatedDecodingLR(matUseRespV1,vecUseTrialTypes,intTypeCV,dblLambda);
%decode V2
matUseRespV2 = matRespV2(indUseTrialC,:);
[dblPerformance,vecDecodedIndex,matPosteriorProbability,matWeightsV2,dblMeanErrorDegs,matConfusion] = doCrossValidatedDecodingLR(matUseRespV2,vecUseTrialTypes,intTypeCV,dblLambda);
%select V1 neurons
intSelectN = 100;
vecAbsW = matWeightsV1(1:(end-1),1)-matWeightsV1(1:(end-1),2);
[vecMaxW,vecReO] = sort(vecAbsW,'descend');
vecSelectNeuronsV1 = sort(vecReO(1:intSelectN));

%select V2 neurons
vecAbsW = matWeightsV2(1:(end-1),1)-matWeightsV2(1:(end-1),2);
[vecMaxW,vecReO] = sort(vecAbsW,'descend');
vecSelectNeuronsV2 = sort(vecReO(1:intSelectN))+intNeuronsV1;

%% rerun decoding with subset
[vecTrialTypes,vecContrasts,vecCounts,cellSelect,vecRepetition] = label2idx(vecTrialContrasts);
vecSubH = zeros(1,numel(vecContrasts)-1);
vecSubM = zeros(1,numel(vecContrasts)-1);
vecSubFA = zeros(1,numel(vecContrasts)-1);
vecSubCR = zeros(1,numel(vecContrasts)-1);
for intC=find(vecContrasts>0)
	intC
	vecUseContrasts = [0 vecContrasts(intC)];
	vecUseIdxC = find(ismember(vecContrasts,vecUseContrasts));
	indUseTrialC = ismember(vecTrialTypes,vecUseIdxC);
	vecUseTrialTypes = vecTrialTypes(indUseTrialC);
	matUseResp = matRespV2(indUseTrialC,vecSelectNeuronsV2-intNeuronsV1);
	dblLambda = 1;
	intTypeCV = 2;
	[dblPerformance,vecDecodedIndex,matPosteriorProbability,matWeightsV2,dblMeanErrorDegs,matConfusion] = doCrossValidatedDecodingLR(matUseResp,vecUseTrialTypes,intTypeCV,dblLambda);
	%assign
	vecSubH(intC-1) = matConfusion(2,2);
	vecSubCR(intC-1) = matConfusion(1,1);
	vecSubM(intC-1) = matConfusion(1,2);
	vecSubFA(intC-1) = matConfusion(2,1);
end

% decode one contrast in V1 & V2
[vecSubP,vecSubCI] = binofit(vecSubH,vecCounts(2:end));
[dummy,intUseC] = min(abs(vecP-0.75));
intUseC = intUseC + 1;
dblUseContrast = vecContrasts(intUseC);
vecUseContrasts = [0 dblUseContrast];
vecUseIdxC = find(ismember(vecContrasts,vecUseContrasts));
indUseTrialC = ismember(vecTrialTypes,vecUseIdxC);
vecUseTrialTypes = vecTrialTypes(indUseTrialC);

%decode V2
matUseRespV2 = matRespV2(indUseTrialC,vecSelectNeuronsV2-intNeuronsV1);
[dblPerformance,vecDecodedIndex,matPosteriorProbability,matWeightsV2,dblMeanErrorDegs,matConfusion] = doCrossValidatedDecodingLR(matUseRespV2,vecUseTrialTypes,intTypeCV,dblLambda);
%get performance V2
dblSubH = matConfusion(2,2);
dblSubCR = matConfusion(1,1);
dblSubM = matConfusion(1,2);
dblSubFA = matConfusion(2,1);

%select only stimulation trials
indSubHit = vecUseTrialTypes==intUseC & vecDecodedIndex'==2;
indSubCRs = vecUseTrialTypes==1 & vecDecodedIndex'==1;
indSubMiss = vecUseTrialTypes==intUseC & vecDecodedIndex'==1;
indSubFA = vecUseTrialTypes==1 & vecDecodedIndex'==2;

