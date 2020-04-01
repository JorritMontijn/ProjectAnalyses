

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

%% subsample population
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

%% build synchronization vectors
dblBinSize_ms = 1;
dblJitter_ms = 5;
dblJitterSecs = dblJitter_ms/1000;
dblBinSizeSecs = dblBinSize_ms/1000;
%build timestamps
dblFirstSpike = min(cellfun(@min,cellSpikesV1V2));
dblLastSpike = max(cellfun(@max,cellSpikesV1V2));
dblStart = dblFirstSpike - dblBinSizeSecs - dblJitterSecs;
dblStop = dblLastSpike + dblBinSizeSecs + dblJitterSecs*2;
vecBinEdges = dblStart:dblBinSizeSecs:dblStop;
vecTimestamps = vecBinEdges(2:end) - dblBinSizeSecs/2;

%get neuron numbers
intNeuronsV1 = sum(vecCellAreaV1V2==1);
intNeuronsV2 = sum(vecCellAreaV1V2==2);
%V1
[vecPopSyncV1,vecTimestampsV1,intJitterBins] = getPopulationSynchronization(cellSpikesV1V2(vecCellAreaV1V2==1),dblJitter_ms,vecTimestamps);
%V2
%[vecPopSyncV2,vecTimestampsV2] = getPopulationSynchronization(cellSpikeTimesCortex(vecCellArea==2),dblJitter_ms,vecTimestamps);

%% analysis
ptrTic = tic;
%get V1 events >2 sd
vecPopSyncZ_V1 = zscore(vecPopSyncV1);
indCutOff = vecPopSyncZ_V1>2;
vecSelectedSyncV1 = vecTimestampsV1(indCutOff);
vecOnOff = diff([0 indCutOff 0]);
vecOn = find(vecOnOff==1);
vecOnSecs = vecTimestamps(vecOn)-dblBinSizeSecs/2;
vecOff = find(vecOnOff==-1)-1;
vecOffSecs = vecTimestamps(vecOff)+dblBinSizeSecs/2;
vecSyncDurSecs = (vecOffSecs - vecOnSecs);

%get spikes
vecSpikesV1 = sort(cell2mat(cellSpikesV1V2(vecCellAreaV1V2==1)));
vecSpikesV2 = sort(cell2mat(cellSpikesV1V2(vecCellAreaV1V2==2)));

%% loop through windows and count V1 sync epoch dur, # of V1 spikes, and # of V2 spikes
%build windows
dblWindow_ms = 25;
dblWindowSecs = dblWindow_ms/1000;
dblWindowStepSecs = 1/1000;
vecWindowEdges = vecTimestampsV1(1):dblWindowSecs:(vecTimestampsV1(end)-dblWindowSecs);
vecWindowStartOffsets = 0:dblWindowStepSecs:(dblWindowSecs-dblWindowStepSecs);
vecWindowStarts = sort(flat(bsxfun(@plus,vecWindowEdges,vecWindowStartOffsets')))';
intWindowOffsets = numel(vecWindowStartOffsets);
intWindowNum = (numel(vecWindowEdges)-1)*intWindowOffsets;
%pre-allocate
vecWinSyncEpochDurV1 = zeros(1,intWindowNum);
vecWinSpikeCountV1 = zeros(1,intWindowNum);
vecWinSpikeCountV2 = zeros(1,intWindowNum);
%loop
for intWinOffset=1:intWindowOffsets
	%% get output indices
	vecAssign = intWinOffset:intWindowOffsets:intWindowNum;
	dblWindowOffset = vecWindowStartOffsets(intWinOffset);
	
	%% get sync
	vecWinSyncEpochDurV1(vecAssign) = histcounts(vecSelectedSyncV1,vecWindowEdges+dblWindowOffset);
	
	%% get V1 & V2 spikes
	vecWinSpikeCountV1(vecAssign) = histcounts(vecSpikesV1,vecWindowEdges+dblWindowOffset);
	vecWinSpikeCountV2(vecAssign) = histcounts(vecSpikesV2,vecWindowEdges+dblWindowOffset);
	
	%% msg
	if toc(ptrTic) > 10
		fprintf('Now at window offset %d/%d [%s]\n',intWinOffset,intWindowOffsets,getTime);
		ptrTic = tic;
	end
end

%% transform to Hz
vecWinSpikeHzV1 = (vecWinSpikeCountV1 / dblWindowSecs) / intNeuronsV1;
vecWinSpikeHzV2 = (vecWinSpikeCountV2 / dblWindowSecs) / intNeuronsV2;
vecWinSyncHzV1 = ((((vecWinSyncEpochDurV1)*(dblBinSizeSecs/dblWindowStepSecs))/dblWindowSecs) / intNeuronsV1) / intJitterBins;

%% prep counts
%bin spike V1
dblMaxCountV1 = max(vecWinSpikeHzV1);
dblMinCountV1 = 0;
if dblMaxCountV1 < 10
	dblStepCountV1 = dblMaxCountV1/10;
else
	dblStepCountV1 = 1;
end
vecEdgesCountV1 = dblMinCountV1:dblStepCountV1:dblMaxCountV1;

%bin sync V1
[vecUniques,dummy,vecValIdx]=unique(vecWinSyncHzV1);
intUniques = numel(vecUniques);
vecValCounts = zeros(1,intUniques);
for intU=1:intUniques
	vecValCounts(intU) = sum(vecValIdx==intU);
end
vecCutOffCount = 10;
intCUV = find(vecValCounts < vecCutOffCount);
if isempty(intCUV),intCUV = intUniques;end
dblSyncStep = vecUniques(2) - vecUniques(1);
vecEdgesSyncV1 = (vecUniques(1) - dblSyncStep/2):dblSyncStep:(vecUniques(intCUV) + dblSyncStep/2);

%% get counts
%build matrix
[matCountsOverall,matValMeansOverall,matValSDsOverall,cellValsOverall,cellIDsOverall] = makeBins2(vecWinSpikeHzV1,vecWinSyncHzV1,vecWinSpikeHzV2,vecEdgesCountV1,vecEdgesSyncV1);
matValMeansOverall(1,1) = nan;
matCountsOverall(1,1) = 0;

% plot heat map
matPlotMeans = matValMeansOverall;
matPlotMeans(matCountsOverall<100) = nan;
figure,plotSyncMap(vecEdgesCountV1,vecEdgesSyncV1,matPlotMeans,dblWindow_ms);
return
%% save
%strFile = sprintf('NeuralMarkerStimDetection3_V1Synchrony_Window%dms',round(dblWindow_ms));
%export_fig([strFile '.tif']);
%export_fig([strFile '.pdf']);
%if 1,return;end

%% prepare contrast analysis
%build trial id vector
intTrials = numel(vecStimStartSecs);
vecWinTrialID = zeros(size(vecWinSpikeCountV1),'single');
for intTrial=1:intTrials
	dblStart = vecStimStartSecs(intTrial);
	dblStop = vecStimStopSecs(intTrial);
	vecAssign = vecTimestampsV1 > dblStart & vecTimestampsV1 < dblStop;
	vecWinTrialID(vecAssign) = intTrial;
end
%% run contrast analysis v2
intNumC = numel(vecContrasts);
for intIdxC = []%1:intNumC
	%% get data
	dblContrast = vecStimTypeContrasts(intIdxC);
	vecUseTrials = find(vecTrialContrasts==dblContrast);
	indUseBins = ismember(vecWinTrialID,vecUseTrials);
	
	%retrieve entries
	vecTheseWinSpikeHzV1 = vecWinSpikeHzV1(indUseBins);
	vecTheseWinSpikeHzV2 = vecWinSpikeHzV2(indUseBins);
	vecTheseWinSyncHzV1 = vecWinSyncHzV1(indUseBins);
	
	%% prep counts
	%bin spike V1
	dblMaxCountV1 = max(vecTheseWinSpikeHzV1);
	dblMinCountV1 = 0;
	if dblMaxCountV1 < 10
		dblStepCountV1 = dblMaxCountV1/10;
	else
		dblStepCountV1 = 1;
	end
	vecEdgesCountV1 = dblMinCountV1:dblStepCountV1:dblMaxCountV1;
	
	%bin sync V1
	[vecUniques,dummy,vecValIdx]=unique(vecTheseWinSyncHzV1);
	intUniques = numel(vecUniques);
	vecValCounts = zeros(1,intUniques);
	for intU=1:intUniques
		vecValCounts(intU) = sum(vecValIdx==intU);
	end
	vecCutOffCount = 10;
	intCUV = find(vecValCounts < vecCutOffCount);
	if isempty(intCUV),intCUV = intUniques;end
	dblSyncStep = vecUniques(2) - vecUniques(1);
	vecEdgesSyncV1 = (vecUniques(1) - dblSyncStep/2):dblSyncStep:(vecUniques(intCUV) + dblSyncStep/2);
	
	%% get counts
	%build matrix
	[matCounts,matValMeans] = makeBins2(vecTheseWinSpikeHzV1,vecTheseWinSyncHzV1,vecTheseWinSpikeHzV2,vecEdgesCountV1,vecEdgesSyncV1);
	matValMeans(1,1) = nan;
	matCounts(1,1) = 0;

	% plot heat map
	matPlotMeans = matValMeans;
	matPlotMeans(matCounts<25) = nan;
	figure,plotSyncMap(vecEdgesCountV1,vecEdgesSyncV1,matPlotMeans,dblWindow_ms);
	subplot(2,2,4)
	title(sprintf('Contrast %.2f',vecContrasts(intIdxC)));
	fixfig
	axis off
end

%% run trial-based analysis: is spiking rate and classification performance in V2 associated with V1 synchronization?
%build trial by neuron matrix
matResp = zeros(intTrials,intNeuronsV2);
vecTrialEdges = zeros(1,size(vecTrialStartSecs,2)*2);
vecTrialEdges(1:2:end) = vecStimStartSecs;
vecTrialEdges(2:2:end) = vecStimStopSecs;
vecTrialEdges(end+1) = vecTrialEdges(end) + median(vecStimStartSecs(2:end) - vecStimStopSecs(1:(end-1)));
dblITIDur = median(vecStimStartSecs(2:end) - vecStimStopSecs(1:(end-1)));
dblTrialDur = median(vecStimStopSecs - vecStimStartSecs);
vecCellsV2 = find(vecCellAreaV1V2==2);
for intNeuronIdx=1:intNeuronsV2
	intCellV2 = vecCellsV2(intNeuronIdx);
	vecSpikes = cellSpikesV1V2{intCellV2};
	vecTrialResp = histcounts(vecSpikes,vecTrialEdges);
	vecTrialResp = vecTrialResp(1:2:end);
	matResp(:,intNeuronIdx) = vecTrialResp;
end

%% go through contrasts and compare with 0
cellRealIndex = cell(1,numel(vecContrasts)-1);
cellDecodedIndex = cell(1,numel(vecContrasts)-1);
cellRelPostProb = cell(1,numel(vecContrasts)-1);
[vecTrialTypes,vecContrasts,vecCounts,cellSelect,vecRepetition] = label2idx(vecTrialContrasts);
vecH = zeros(1,numel(vecContrasts)-1);
vecM = zeros(1,numel(vecContrasts)-1);
vecFA = zeros(1,numel(vecContrasts)-1);
vecCR = zeros(1,numel(vecContrasts)-1);
for intC=find(vecContrasts>0)
	intC
	vecUseContrasts = [0 vecContrasts(intC)];
	vecUseIdxC = find(ismember(vecContrasts,vecUseContrasts));
	indUseTrialC = ismember(vecTrialTypes,vecUseIdxC);
	vecUseTrialTypes = vecTrialTypes(indUseTrialC);
	matUseResp = matResp(indUseTrialC,:);
	dblLambda = 1;
	intTypeCV = 2;
	[dblPerformance,vecDecodedIndex,matPosteriorProbability,matWeights,dblMeanErrorDegs,matConfusion] = doCrossValidatedDecodingLR(matUseResp,vecUseTrialTypes,intTypeCV,dblLambda);
	%assign
	cellRealIndex{intC-1} = label2idx(vecUseTrialTypes);
	cellDecodedIndex{intC-1} = vecDecodedIndex;
	cellRelPostProb{intC-1} = matPosteriorProbability(:,2) - matPosteriorProbability(:,1);
	vecH(intC-1) = matConfusion(2,2);
	vecCR(intC-1) = matConfusion(1,1);
	vecM(intC-1) = matConfusion(1,2);
	vecFA(intC-1) = matConfusion(2,1);
end

%plot
[vecP,vecCI] = binofit(vecH,vecCounts(2:end));
figure;
hold on
plot([vecContrasts(2) vecContrasts(end)],[0.5 0.5],'k--')
errorbar(vecContrasts(2:end),vecP,vecCI(:,1)-vecP',vecCI(:,2)-vecP');
hold off
set(gca,'xscale','log')
xlabel('Visual contrast (%)');
ylabel('Stimulus detection V2 (hit rate)')
title('LR decoder V2 population')
fixfig

%% analyze one contrast
[dummy,intUseC] = min(abs(vecP-0.75));
intUseC = 16;%intUseC + 1;
dblUseContrast = vecContrasts(intUseC);
vecUseContrasts = [0 dblUseContrast];
vecUseIdxC = find(ismember(vecContrasts,vecUseContrasts));
indUseTrialC = ismember(vecTrialTypes,vecUseIdxC);
vecUseTrialTypes = vecTrialTypes(indUseTrialC);
matUseResp = matResp(indUseTrialC,:);
dblLambda = 1;
intTypeCV = 2;
[dblPerformance,vecDecodedIndex,matPosteriorProbability,matWeights,dblMeanErrorDegs,matConfusion] = doCrossValidatedDecodingLR(matUseResp,vecUseTrialTypes,intTypeCV,dblLambda);
dblH = matConfusion(2,2);
dblCR = matConfusion(1,1);
dblM = matConfusion(1,2);
dblFA = matConfusion(2,1);


indHit = vecUseTrialTypes==intUseC & vecDecodedIndex'==2;
indCRs = vecUseTrialTypes==1 & vecDecodedIndex'==1;
indMiss = vecUseTrialTypes==intUseC & vecDecodedIndex'==1;
indFA = vecUseTrialTypes==1 & vecDecodedIndex'==2;

%% select only stimulation trials
vecRealIndex = label2idx(vecUseTrialTypes);
vecDecodedIndex = vecDecodedIndex;
vecRelPostProb = matPosteriorProbability(:,2) - matPosteriorProbability(:,1);

vecUseStimuli = 16;
vecRealIndex = cell2mat(cellRealIndex(vecUseStimuli));
vecDecodedIndex = cell2mat(cellDecodedIndex(vecUseStimuli)')';
vecRelPostProb = cell2mat(cellRelPostProb(vecUseStimuli)')';
vecUseTrials = [];
for intS=vecUseStimuli
	vecUseTrials = cat(2,vecUseTrials,find(ismember(vecTrialTypes,[1 intS+1])));
end

%% get trial-based sync data
dblUseTrialOffset = 0;%0.01;
dblUseTrialDur = 0.5;%0.05
%select trial windows
vecTrialEdges = zeros(1,size(vecTrialStartSecs,2)*2);
vecTrialEdges(1:2:end) = vecStimStartSecs+dblUseTrialOffset;
vecTrialEdges(2:2:end) = vecStimStartSecs+dblUseTrialOffset+dblUseTrialDur;
vecTrialEdges(end+1) = vecTrialEdges(end) + median(vecStimStartSecs(2:end) - vecStimStopSecs(1:(end-1)));

%subselect data
vecOrigTrialIdx = find(indUseTrialC);
%vecUseTrials = vecOrigTrialIdx;
vecUseDecodedIndex = vecDecodedIndex;
indUseBins = ismember(vecWinTrialID,vecUseTrials);

%retrieve entries
vecTheseTimestamps = vecTimestampsV1(indUseBins);
vecTheseWinSpikeHzV1 = vecWinSpikeHzV1(indUseBins);
vecTheseWinSpikeHzV2 = vecWinSpikeHzV2(indUseBins);
vecTheseWinSyncHzV1 = vecWinSyncHzV1(indUseBins);
	
%get data
[vecCountsSpikeHzV2,vecMeansSpikeHzV2] = ...
	makeBins(vecTheseTimestamps,vecTheseWinSpikeHzV2,vecTrialEdges);
[vecCountsSyncHzV1,vecMeansSyncHzV1] = ...
	makeBins(vecTheseTimestamps,vecTheseWinSyncHzV1,vecTrialEdges);
[vecCountsSpikeHzV1,vecMeansSpikeHzV1] = ...
	makeBins(vecTheseTimestamps,vecTheseWinSpikeHzV1,vecTrialEdges);

% get trial data
%spike V2
	vecTheseMeansSpikeHzV2 = vecMeansSpikeHzV2(1:2:end);
	vecTheseMeansSpikeHzV2 = vecTheseMeansSpikeHzV2(vecUseTrials);
%sync V1
	vecTheseMeansSyncHzV1 = vecMeansSyncHzV1(1:2:end);
	vecTheseMeansSyncHzV1 = vecTheseMeansSyncHzV1(vecUseTrials);
%spike V1
	vecTheseMeansSpikeHzV1 = vecMeansSpikeHzV1(1:2:end);
	vecTheseMeansSpikeHzV1 = vecTheseMeansSpikeHzV1(vecUseTrials);
%reliability detection
	vecTheseDetectProbs = vecRelPostProb;%matPosteriorProbability(:,2);%matPosteriorProbability(:,2);%vecRelPostProb

	%% regress
	[vecUseX,vecReorder] = sort(vecTheseMeansSyncHzV1(vecRealIndex==2));
	vecUseY = vecTheseDetectProbs(vecRealIndex==2)';
	vecUseY = vecUseY(vecReorder);
	matX = [ones(size(vecUseX)) vecUseX];
	
	mdl = fitlm(vecUseX, vecUseY, 'linear');
	xnew = linspace(min(vecUseX), max(vecUseX), 1000)';
	[ypred, yci] = predict(mdl, xnew);
	figure
	scatter(vecUseX, vecUseY, 'b.')
	hold on
	plot(xnew, ypred, '-r')
	plot(xnew, yci(:,1), '-g')
	plot(xnew, yci(:,2), '-g')
	hold off

	%%
%plot
figure
subplot(2,2,1);
[RSp2B,PSp2B,RLSp2B,RUSp2B] = corrcoef(vecTheseMeansSpikeHzV2(vecRealIndex==1),vecTheseDetectProbs(vecRealIndex==1));
[RSp2S,PSp2S,RLSp2S,RUSp2S] = corrcoef(vecTheseMeansSpikeHzV2(vecRealIndex==2),vecTheseDetectProbs(vecRealIndex==2));
hold on
scatter(vecTheseMeansSpikeHzV2(vecRealIndex==1),vecTheseDetectProbs(vecRealIndex==1),'rx');
scatter(vecTheseMeansSpikeHzV2(vecRealIndex==2),vecTheseDetectProbs(vecRealIndex==2),'gx');
hold off
title(sprintf('Stim (g), r=%.3f (p=%.3f); Base (r) r=%.3f (p=%.3f)',RSp2S(1,2),PSp2S(1,2),RSp2B(1,2),PSp2B(1,2)));
xlabel('V2 spiking rate per trial per neuron (Hz)');
ylabel('Decoder certainty on stim. presence');
%set(gca,'yscale','log')
fixfig;

subplot(2,2,2);
[RSy1B,PSy1B,RLSy1B,RUSy1B] = corrcoef(vecTheseMeansSyncHzV1(vecRealIndex==1),vecTheseDetectProbs(vecRealIndex==1));
[RSy1S,PSy1S,RLSy1S,RUSy1S] = corrcoef(vecTheseMeansSyncHzV1(vecRealIndex==2),vecTheseDetectProbs(vecRealIndex==2));
hold on
scatter(vecTheseMeansSyncHzV1(vecRealIndex==1),vecTheseDetectProbs(vecRealIndex==1),'rx');
scatter(vecTheseMeansSyncHzV1(vecRealIndex==2),vecTheseDetectProbs(vecRealIndex==2),'gx');
hold off
title(sprintf('Stim (g), r=%.3f (p=%.3f); Base (r) r=%.3f (p=%.3f)',RSy1S(1,2),PSy1S(1,2),RSy1B(1,2),PSy1B(1,2)));
xlabel('V1 sync. rate per trial per neuron (Hz)');
ylabel('Decoder certainty on stim. presence');
%set(gca,'yscale','log')
fixfig;

subplot(2,2,3);
[RSp1B,PSp1B,RLSp1B,RUSp1B] = corrcoef(vecTheseMeansSpikeHzV1(vecRealIndex==1),vecTheseDetectProbs(vecRealIndex==1));
[RSp1S,PSp1S,RLSp1S,RUSp1S] = corrcoef(vecTheseMeansSpikeHzV1(vecRealIndex==2),vecTheseDetectProbs(vecRealIndex==2));
hold on
scatter(vecTheseMeansSpikeHzV1(vecRealIndex==1),vecTheseDetectProbs(vecRealIndex==1),'rx');
scatter(vecTheseMeansSpikeHzV1(vecRealIndex==2),vecTheseDetectProbs(vecRealIndex==2),'gx');
hold off
title(sprintf('Stim (g), r=%.3f (p=%.3f); Base (r) r=%.3f (p=%.3f)',RSp1S(1,2),PSp1S(1,2),RSp1B(1,2),PSp1B(1,2)));
xlabel('V1 spiking rate per trial per neuron (Hz)');
ylabel('Decoder certainty on stim. presence');
%set(gca,'yscale','log')
fixfig;

%% make colormap plot
strContrastsUsed = sprintf('%.0f; ',vecContrasts(vecUseStimuli+1));
close all
figure
vecSubDetP = vecUseY;
vecSubSync = vecUseX;

dblSubSyncStep = 0.0005/10;
vecSubSyncEdges = 0:dblSubSyncStep:0.01;
dblSubDecStep = 0.05/10;
vecSubDecEdges = 0:dblSubDecStep:1;
matC=histcounts2(vecSubDetP(:)',vecSubSync(:)',vecSubDecEdges,vecSubSyncEdges);

%smooth
matFilt = normpdf(-20:20,0,10)'*normpdf(-20:20,0,10);
matFilt = matFilt / sum(matFilt(:));
matSmoothC = conv2(matC,matFilt,'same');
vecPlotX = vecSubSyncEdges(2:end) - dblSubSyncStep/2;
vecPlotY = vecSubDecEdges(2:end) - dblSubDecStep/2;
imagesc(vecPlotX,vecPlotY,matSmoothC);colormap(grey);freezeColors;
hold on
contour(vecPlotX,vecPlotY,matSmoothC);colormap(redwhite);freezeColors;colorbar
%scatter(vecUseX, vecUseY, 'b.')
%plot(xnew, ypred, '-r')
%plot(xnew, yci(:,1), '-g')
%plot(xnew, yci(:,2), '-g')
hold off
xlabel('V1 sync. rate per trial per neuron (Hz)');
ylabel('Decoder certainty on stim. presence');
title(sprintf('Contrasts %s',strContrastsUsed));
axis xy
xlim([min(vecSubSyncEdges) max(vecSubSyncEdges)]);ylim([min(vecSubDecEdges) max(vecSubDecEdges)]);fixfig;grid off;



figure
%smooth
matFilt = normpdf(-20:20,0,10)'*normpdf(-20:20,0,10);
matFilt = matFilt / sum(matFilt(:));
matSmoothC = conv2(matC,matFilt,'same');
vecPlotX = vecSubSyncEdges(2:end) - dblSubSyncStep/2;
vecPlotY = vecSubDecEdges(2:end) - dblSubDecStep/2;
contour(vecPlotX,vecPlotY,matSmoothC);colormap(redwhite);freezeColors;colorbar
hold on
	plot(xnew, ypred, '-r')
	plot(xnew, yci(:,1), '-g')
	plot(xnew, yci(:,2), '-g')
	hold off
xlabel('V1 sync. rate per trial per neuron (Hz)');
ylabel('Decoder certainty on stim. presence');
title(sprintf('Contrasts %s',strContrastsUsed));
axis xy
xlim([min(vecSubSyncEdges) max(vecSubSyncEdges)]);ylim([min(vecSubDecEdges) max(vecSubDecEdges)]);fixfig;grid off;

figure
%smooth
matFilt = normpdf(-20:20,0,10)'*normpdf(-20:20,0,10);
matFilt = matFilt / sum(matFilt(:));
matSmoothC = conv2(matC,matFilt,'same');
vecPlotX = vecSubSyncEdges(2:end) - dblSubSyncStep/2;
vecPlotY = vecSubDecEdges(2:end) - dblSubDecStep/2;
imagesc(vecPlotX,vecPlotY,matSmoothC);colormap(grey);freezeColors;colorbar
xlabel('V1 sync. rate per trial per neuron (Hz)');
ylabel('Decoder certainty on stim. presence');
title(sprintf('Contrasts %s',strContrastsUsed));
axis xy
xlim([min(vecSubSyncEdges) max(vecSubSyncEdges)]);ylim([min(vecSubDecEdges) max(vecSubDecEdges)]);fixfig;grid off;
