

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
matPlotMeans(matCountsOverall<250) = nan;
figure,plotSyncMap(vecEdgesCountV1,vecEdgesSyncV1,matPlotMeans,dblWindow_ms);

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
	cellRelPostProb{intC-1} = matPosteriorProbability(:,1) - matPosteriorProbability(:,2);
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
[dummy,intUseC] = min(abs(vecP-1));
intUseC = intUseC + 1;
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

%select only stimulation trials
vecRelPostProb = matPosteriorProbability(:,2) - matPosteriorProbability(:,1);
indHit = vecUseTrialTypes==intUseC & vecDecodedIndex'==2;
indCRs = vecUseTrialTypes==1 & vecDecodedIndex'==1;
indMiss = vecUseTrialTypes==intUseC & vecDecodedIndex'==1;
indFA = vecUseTrialTypes==1 & vecDecodedIndex'==2;


%% get trial-based sync data
%select neurons
vecAbsW = matWeights(1:(end-1),1)-matWeights(1:(end-1),2);
[vecMaxW,vecReO] = sort(vecAbsW,'descend');
intSelectN = 100;
vecSelectNeurons = vecReO(1:intSelectN);
indCorr = indHit;
indIncorr = indCRs;

dblUseTrialOffset = 0;%0.01;
dblUseTrialDur = 0.5;%0.05
%select trial windows
vecTrialEdges = zeros(1,size(vecTrialStartSecs,2)*2);
vecTrialEdges(1:2:end) = vecStimStartSecs+dblUseTrialOffset;
vecTrialEdges(2:2:end) = vecStimStartSecs+dblUseTrialOffset+dblUseTrialDur;
vecTrialEdges(end+1) = vecTrialEdges(end) + median(vecStimStartSecs(2:end) - vecStimStopSecs(1:(end-1)));

%subselect data
vecOrigTrialIdx = find(indUseTrialC);
vecIncorrect = vecOrigTrialIdx(indIncorr);
vecCorrect = vecOrigTrialIdx(indCorr);
vecUseTrials = sort(cat(2,vecCorrect,vecIncorrect));
vecUseDecodedIndex = vecDecodedIndex(indCorr | indIncorr);
indUseBins = ismember(vecWinTrialID,vecUseTrials);

%retrieve entries
vecTheseTimestamps = vecTimestampsV1(indUseBins);
vecTheseWinSpikeHzV1 = vecWinSpikeHzV1(indUseBins);
vecTheseWinSpikeHzV2 = vecWinSpikeHzV2(indUseBins);
vecTheseWinSyncHzV1 = vecWinSyncHzV1(indUseBins);
	
%get data
[vecCountsSpikeHzV2,vecMeansSpikeHzV2,vecSDsSpikeHzV2,cellValsSpikeHzV2,cellIDsSpikeHzV2] = ...
	makeBins(vecTheseTimestamps,vecTheseWinSpikeHzV2,vecTrialEdges);
[vecCountsSyncHzV1,vecMeansSyncHzV1,vecSDsSyncHzV1,cellValsSyncHzV1,cellIDsSyncHzV1] = ...
	makeBins(vecTheseTimestamps,vecTheseWinSyncHzV1,vecTrialEdges);
[vecCountsSpikeHzV1,vecMeansSpikeHzV1,vecSDsSpikeHzV1,cellValsSpikeHzV1,cellIDsSpikeHzV1] = ...
	makeBins(vecTheseTimestamps,vecTheseWinSpikeHzV1,vecTrialEdges);

% get trial data
intSwitch = 1;
if intSwitch == 1
	hFunc = @mean;
elseif intSwitch == 2
	hFunc = @max;
elseif intSwitch == 3
	hFunc = @median;
elseif intSwitch == 4
	hFunc = @min;
elseif intSwitch == 5
	hFunc = @std;
end
%spike V2
	vecTheseMeansSpikeHzV2 = cellValsSpikeHzV2(1:2:end);
	vecTheseMeansSpikeHzV2 = vecTheseMeansSpikeHzV2(vecUseTrials);
	vecTheseMeansSpikeHzV2 = cellfun(hFunc,vecTheseMeansSpikeHzV2);
%sync V1
	vecTheseMeansSyncHzV1 = cellValsSyncHzV1(1:2:end);
	vecTheseMeansSyncHzV1 = vecTheseMeansSyncHzV1(vecUseTrials);
	vecTheseMeansSyncHzV1 = cellfun(hFunc,vecTheseMeansSyncHzV1);
%spike V1

	vecTheseMeansSpikeHzV1 = cellValsSpikeHzV1(1:2:end);
	vecTheseMeansSpikeHzV1 = vecTheseMeansSpikeHzV1(vecUseTrials);
	vecTheseMeansSpikeHzV1 = cellfun(hFunc,vecTheseMeansSpikeHzV1);

% run correct/incorrect trials
close all
%set edges
%bin spike V1
vecEdgesCountV1 = 0:dblStepCountV1:25;

%bin sync V1
dblSyncStep = vecUniques(2) - vecUniques(1);
vecEdgesSyncV1 = (0 - dblSyncStep/2):dblSyncStep:(0.045 + dblSyncStep/2);
	
%edges for distros
dblStepSp1 = 0.5;
vecBinSpikeV1 = 0:dblStepSp1:30;
vecPlotSpikeV1 = (dblStepSp1/2):dblStepSp1:(30-dblStepSp1/2);
vecCountsSpikeV1 = cell(1,2);

dblStepSp2 = 1;
vecBinSpikeV2 = 0:dblStepSp2:60;
vecPlotSpikeV2 = (dblStepSp2/2):dblStepSp2:(60-dblStepSp2/2);
vecCountsSpikeV2 = cell(1,2);

dblStepSy1 = dblSyncStep;
vecBinSyncV1 = (vecUniques(1) - dblSyncStep/2):dblSyncStep:(vecUniques(end) + dblSyncStep/2);
vecPlotSyncV1 = vecUniques(1):dblStepSy1:vecUniques(end);
vecCountsSyncV1 = cell(1,2);

%run
figure;
hAxScat1 = subplot(2,2,1);hold on;
hAxScat2 = subplot(2,2,2);hold on;
hAxScat3 = subplot(2,2,3);hold on;
figure;
hAxTrialScat1 = subplot(2,2,1);hold on;
hAxTrialScat2 = subplot(2,2,2);hold on;
hAxTrialScat3 = subplot(2,2,3);hold on;
for intCorr=[1 2]
	if intCorr == 1
		strCorrIncorr = 'Miss trials';
		vecUseTheseTrials = vecIncorrect;
		vecColor = [0.7 0 0];
	else
		strCorrIncorr = 'Hit trials';
		vecUseTheseTrials = vecCorrect;
		vecColor = [0 0.7 0];
	end
	indUseTheseBins = ismember(vecWinTrialID,vecUseTheseTrials);
	
	%retrieve entries
	vecTheseWinSpikeHzV1 = vecWinSpikeHzV1(indUseTheseBins);
	vecTheseWinSpikeHzV2 = vecWinSpikeHzV2(indUseTheseBins);
	vecTheseWinSyncHzV1 = vecWinSyncHzV1(indUseTheseBins);
	
	%% get counts
	%build matrix
	[matCounts,matValMeans] = makeBins2(vecTheseWinSpikeHzV1,vecTheseWinSyncHzV1,vecTheseWinSpikeHzV2,vecEdgesCountV1,vecEdgesSyncV1);
	matValMeans(1,1) = nan;
	matCounts(1,1) = 0;

	% plot heat map
	matPlotMeans = matValMeans;
	matPlotMeans(matCounts<25) = nan;
	figure,plotSyncMap(vecEdgesCountV1,vecEdgesSyncV1,matPlotMeans,dblWindow_ms,false);
	subplot(2,2,4)
	title(sprintf('Contrast %.2f; %s',dblUseContrast, strCorrIncorr));
	fixfig
	axis off
	
	%% plot trial-mean distros
	indUseEntries = vecUseDecodedIndex==intCorr;
	intTotBins = sum(indUseEntries);
	%vecCountsSpikeV1{intCorr} = histcounts(vecTheseWinSpikeHzV1,vecBinSpikeV1)/intTotBins;
	%vecCountsSpikeV2{intCorr} = histcounts(vecTheseWinSpikeHzV2,vecBinSpikeV2)/intTotBins;
	%vecCountsSyncV1{intCorr} = histcounts(vecTheseWinSyncHzV1,vecBinSyncV1)/intTotBins;
	vecCountsSpikeV1{intCorr} = histcounts(vecTheseMeansSpikeHzV1(indUseEntries),vecBinSpikeV1)/intTotBins;
	vecCountsSpikeV2{intCorr} = histcounts(vecTheseMeansSpikeHzV2(indUseEntries),vecBinSpikeV2)/intTotBins;
	vecCountsSyncV1{intCorr} = histcounts(vecTheseMeansSyncHzV1(indUseEntries),vecBinSyncV1)/intTotBins;

	%get trial data
	axes(hAxScat1);
	plot(vecPlotSpikeV1,vecCountsSpikeV1{intCorr},'Color',vecColor);
	xlabel('Average V1 spiking rate (Hz)');
	ylabel('Normalized count');
	fixfig;
	
	axes(hAxScat2);
	plot(vecPlotSpikeV2,vecCountsSpikeV2{intCorr},'Color',vecColor);
	xlabel('Average V2 spiking rate (Hz)');
	ylabel('Normalized count');
	fixfig;
	
	axes(hAxScat3);
	plot(vecPlotSyncV1,vecCountsSyncV1{intCorr},'Color',vecColor);
	xlabel('Av. sync. rate per V1 neuron (Hz)');
	ylabel('Normalized count');
	fixfig;
	
	%% plot trial-mean distros
	%get trial data
	axes(hAxTrialScat1);
	scatter(vecTheseMeansSpikeHzV1(vecUseDecodedIndex==intCorr),vecTheseMeansSyncHzV1(vecUseDecodedIndex==intCorr),[],vecColor);
	xlabel('Average V1 spiking rate (Hz)');
	ylabel('Av. sync. rate per V1 neuron (Hz)');
	title('Trial average')
	fixfig;
	
	axes(hAxTrialScat2);
	scatter(vecTheseMeansSpikeHzV1(vecUseDecodedIndex==intCorr),vecTheseMeansSpikeHzV2(vecUseDecodedIndex==intCorr),[],vecColor);
	xlabel('Average V1 spiking rate (Hz)');
	ylabel('Average V2 spiking rate (Hz)');
	title('Trial average')
	fixfig;
	
	axes(hAxTrialScat3);
	scatter(vecTheseMeansSyncHzV1(vecUseDecodedIndex==intCorr),vecTheseMeansSpikeHzV2(vecUseDecodedIndex==intCorr),[],vecColor);
	xlabel('Av. sync. rate per V1 neuron (Hz)');
	ylabel('Average V2 spiking rate (Hz)');
	title('Trial average')
	fixfig;
end
% plotdistro diffs
	%get trial data
	figure
	subplot(2,2,1)
	plot(vecPlotSpikeV1,vecCountsSpikeV1{2}-vecCountsSpikeV1{1},'Color','b');
	xlabel('Average V1 spiking rate (Hz)');
	ylabel('Normalized count');
	
	[h,pSp1]=ttest2(vecTheseMeansSpikeHzV1(vecUseDecodedIndex==2),vecTheseMeansSpikeHzV1(vecUseDecodedIndex==1));
	dblMissMeanSp1 = mean(vecTheseMeansSpikeHzV1(vecUseDecodedIndex==1));
	dblHitMeanSp1 = mean(vecTheseMeansSpikeHzV1(vecUseDecodedIndex==2));
	title(sprintf('%s spikes V1, Miss=%.3f; Hit=%.3f, p=%.3f',char(hFunc),dblMissMeanSp1,dblHitMeanSp1,pSp1));
	fixfig;
	
	subplot(2,2,2)
	plot(vecPlotSpikeV2,vecCountsSpikeV2{2}-vecCountsSpikeV2{1},'Color','b');
	xlabel('Average V2 spiking rate (Hz)');
	ylabel('Normalized count');
	[h,pSp2]=ttest2(vecTheseMeansSpikeHzV2(vecUseDecodedIndex==2),vecTheseMeansSpikeHzV2(vecUseDecodedIndex==1));
	dblMissMeanSp2 = mean(vecTheseMeansSpikeHzV2(vecUseDecodedIndex==1));
	dblHitMeanSp2 = mean(vecTheseMeansSpikeHzV2(vecUseDecodedIndex==2));
	title(sprintf('%s spikes V2, Miss=%.3f; Hit=%.3f, p=%.3f',char(hFunc),dblMissMeanSp2,dblHitMeanSp2,pSp2));
	fixfig;
	
	subplot(2,2,3)
	plot(vecPlotSyncV1,vecCountsSyncV1{2}-vecCountsSyncV1{1},'Color','b');
	xlabel('Av. sync. rate per V1 neuron (Hz)');
	ylabel('Normalized count');
	[h,pSy1]=ttest2(vecTheseMeansSyncHzV1(vecUseDecodedIndex==2),vecTheseMeansSyncHzV1(vecUseDecodedIndex==1));
	dblMissMeanSy1 = mean(vecTheseMeansSyncHzV1(vecUseDecodedIndex==1));
	dblHitMeanSy1 = mean(vecTheseMeansSyncHzV1(vecUseDecodedIndex==2));
	title(sprintf('%s sync V1 x100, Miss=%.3f; Hit=%.3f, p=%.3f',char(hFunc),dblMissMeanSy1*100,dblHitMeanSy1*100,pSy1));
	fixfig;
	
%% run analysis to check per V2 spike per cell what the preceding and following spikes in V1 and V2 are
% what is the ratio in predictability of V1 and V2?
% are population events in V1 more common during V2 spikes?
% are there specific temporal patterns of activation?
if 1,return,end

%% loop through stimulus types
boolMakePlots = true;
vecFitOriDegs = 0:359;
cellArea = {'V1','V2'};
intNeurons = numel(cellSpikeTimesCortex);
dblStepT = mean(diff(vecOverallT));
for intStim=1:intStimTypes
	%% get stim types
	intStimIdx = vecUniqueStimTypes(intStim);
	indUseTrials = vecTrialStimType==intStimIdx;
	intRepetitions = sum(indUseTrials);
	vecUseTrials = find(indUseTrials);
	
	[vecTrialPerSpike,vecTimePerSpike] = getSpikesInTrial(vecSpikes,vecTrialStarts);
	%% loop through trials
	for intRep=1:intRepetitions
		%% get trial and timing
		intTrial = vecUseTrials(intRep);
		dblStimOnSecs = vecStimStartSecs(intTrial);
		dblStimOffSecs = vecStimStopSecs(intTrial);
		return
		%% analyze markers in V1 and quantify effect on V2
		
	end
end