%% load data
strSimulation = 'ScaledUp2016-11-04';
load(['D:\Simulations\Results\Simulation_' strSimulation '.mat']);

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
%cellSpikeTimesCortex
vecCellSize = size(cellSpikeTimesCortex);
intTrials = numel(vecTrialOris);
intNeurons = numel(cellSpikeTimesCortex);
matResp = nan(intNeurons,intTrials);
for intTrial=1:intTrials
	dblStartT = vecStimStartSecs(intTrial);
	dblStopT = vecStimStopSecs(intTrial);
	
	matResp(:,intTrial) = cellfun(@sum,...
		cellfun(@and,cellfun(@gt,cellSpikeTimesCortex,cellfill(dblStartT,vecCellSize),'UniformOutput',false),...
		cellfun(@lt,cellSpikeTimesCortex,cellfill(dblStopT,vecCellSize),'UniformOutput',false),...
		'UniformOutput',false));
	if mod(intTrial,100) == 0
		fprintf('Now at trial %d/%d [%s]\n',intTrial,intTrials,getTime);
	end
end

%%
structStim.Orientation = vecTrialOris;
structStim.TrialNumber = 1:length(vecTrialOris);
structStim.FrameOn = vecStimStartSecs;
structStim.FrameOff = vecStimStopSecs;
vecOrientations = unique(vecTrialOris);
vecOriDegs = rad2deg(vecOrientations);

%%
sTypes = getStimulusTypes(structStim);
cellSelect = getSelectionVectors(structStim,sTypes);

%% transform to binned data
dblBinSize = 0.05; %50 ms bins
vecStart= 0:dblBinSize:max(vecOverallT);
matSpikeCounts = getSpikeCounts(cellSpikeTimesCortex,vecStart,dblBinSize);

%plot
sEvents.vecOn = vecStimStartSecs/dblBinSize;
sEvents.cellSelect = cellSelect;
sEvents.sTypes = sTypes;
sEvents.vecWindow = [-0.5/dblBinSize 1/dblBinSize];
sEvents.dblFrameRate = 1/dblBinSize;

for intNeuron=1:size(matSpikeCounts,1)
[cellHandles,sOut] = doPEP(sEvents,matSpikeCounts(intNeuron,:));
pause
close all
end
	%doPEP Performs Peri-Event Plot of supplied trace
	%syntax: cellHandles = doPEP(sEvents,vecTrace)
	%	input:
	%	- sEvents, structure containing the following fields:
	%		- handleFig; handle to figure
	%		- vecOn; vector containing starts of events
	%		- vecOff; vector containing stops of events [default; same as vecOn]
	%		- cellSelect; combo from getStimulusTypes /
	%			getSelectionVectors supplying event types [optional input]
	%			- if cellSelect is defined; sTypes is an optional field for
	%				supplying variable names
	%		- vecTypes; vector supplying event types [optional input]
	%			- if vecTypes is defined; varNames is an optional field for
	%				supplying variable names
	%		- vecWindow; 2-element vector specifying which window to plot
	%			around events [default; [-75 125]]
	%		- dblFrameRate; frame rate [optional; used for setting x-axis]
	%	- vecTrace; trace containing data to be plotted
	%
	

%% 
vecPlotN = randi(intCortexCells,[1 9]);
for intN=1:9%[5 432 887]%1:size(matResp,1)
	subplot(3,3,intN)
doPlotTuningCurve(matResp(vecPlotN(intN),:),structStim,false)
title(sprintf('neuron %d',vecPlotN(intN)))
xlabel('stimulus orientation (radians)')
ylabel('mean firing rate')
end
%export figure
%export_fig(['D:\Results\V1SubPopulations\simulations\' strSimulation '_ExampleTuningCurves.tif']);
%export_fig(['D:\Results\V1SubPopulations\simulations\' strSimulation '_ExampleTuningCurves.pdf']);

%%
clf;
[dblPerformanceML,vecDecodedIndexML,matPosteriorProbabilityML,dblMeanErrorDegsML,matConfusionML] = doCrossValidatedDecodingML(matResp',vecTrialOris,1);

%%
vecTrialOriIdx = label2idx(vecTrialOris);
vecCorrect = vecDecodedIndexML' == vecTrialOriIdx;
vecIncorrect = ~vecCorrect;

vecCorMeanR = mean(matResp(:,vecCorrect),2);
vecIncMeanR = mean(matResp(:,vecIncorrect),2);

vecRespPyrC=vecCorMeanR(vecCellTypes==1);
vecRespPyrI=vecIncMeanR(vecCellTypes==1);

vecRespIntC=vecCorMeanR(vecCellTypes==2);
vecRespIntI=vecIncMeanR(vecCellTypes==2);

vecCorMeanR-vecIncMeanR
%{
%plot
clf
subplot(2,3,1)
plot([0 20],[0 20],'k--')
hold on
scatter(vecRespPyrC,vecRespPyrI,'g^');
scatter(vecRespIntC,vecRespIntI,'ro');
xlim([0 20])
ylim([0 20])
xlabel('Mean response correctly decoded trials')
ylabel('Mean response incorrectly decoded trials')
title('Green=pyramid; red=interneuron')


% >> why are some trials incorrect?
mean(vecRespIntC-vecRespIntI)
mean(vecRespPyrC-vecRespPyrI)

subplot(2,3,4)
dblBinSize = 0.5;
vecBins = 0:dblBinSize:20;
intBins = length(vecBins);
vecPC = min(max(round(vecRespPyrC/dblBinSize),1),intBins);
vecPI = min(max(round(vecRespPyrI/dblBinSize),1),intBins);

vecIC = min(max(round(vecRespIntC/dblBinSize),1),intBins);
vecII = min(max(round(vecRespIntI/dblBinSize),1),intBins);
vecII(vecII>intBins)=intBins;

matProbDensPyrCI= getFillGrid(zeros(intBins),vecPI,vecPC);
matProbDensIntCI = getFillGrid(zeros(intBins),vecII,vecIC);
matProbDensCI(:,:,1) = -log(matProbDensPyrCI+1);
matProbDensCI(:,:,2) = -log(matProbDensIntCI+1);
matIm = 1-imnorm(matProbDensCI);
matIm(:,:,3) = matIm(:,:,2)+matIm(:,:,1);

image([0 20],[0 20],1-matIm);
axis xy
hold on
plot([0 intBins],[0 intBins],'k--')
xlabel('Mean response correctly decoded trials')
ylabel('Mean response incorrectly decoded trials')
title('log prob density')

%export figure
export_fig(['D:\Results\V1SubPopulations\simulations\' strSimulation '_CorrIncorrDecoding.tif']);
export_fig(['D:\Results\V1SubPopulations\simulations\' strSimulation '_CorrIncorrDecoding.pdf']);
%}
%% plot rest
figure
subplot(2,3,1)
imagesc(vecOriDegs,vecOriDegs,matConfusionML,[0 1]);colormap(hot);colorbar;freezeColors;
set(gca,'xtick',vecOriDegs(1:2:end),'ytick',vecOriDegs(1:2:end));
xlabel('Presented orientation');
ylabel('Decoded orientation');
title(sprintf('naive Bayes decoder, correct=%.3f',dblPerformanceML))
fixfig;grid off;box on;
return
%decode shuffled
intIters=10;
intOris = numel(vecOriDegs);
vecShuffledPerformanceML = nan(1,intIters);
matShuffledConfusionML = zeros(intOris,intOris);
for i=1:intIters
	fprintf('Decoding ori, shuffle iter %d/%d [%s]\n',i,intIters,getTime);
	[dblShuffledPerformanceML,x,y,z,matThisShuffledConfusionML] = doCrossValidatedDecodingML(matResp,vecTrialOris(randperm(numel(vecTrialOris))),1);
	vecShuffledPerformanceML(i) = dblShuffledPerformanceML;
	matShuffledConfusionML = matShuffledConfusionML + matThisShuffledConfusionML;
end
matShuffledConfusionML = matShuffledConfusionML / intIters;

subplot(2,3,4)
imagesc(vecOriDegs,vecOriDegs,matShuffledConfusionML,[0 1]);colormap(hot);colorbar;freezeColors;
set(gca,'xtick',vecOriDegs(1:2:end),'ytick',vecOriDegs(1:2:end));
xlabel('Presented orientation');
ylabel('Decoded orientation');
title(sprintf('Shuffled; naive Bayes decoder, correct=%.3f',mean(vecShuffledPerformanceML)))
fixfig;grid off;box on;

% tuning curves
structOut = calcStimCorrsRespMat(matResp,cellSelect);
intNeurons = size(structOut.matSignalCorrs,1);
subplot(2,3,2)
imagesc(structOut.matSignalCorrs,[-1 1]);colormap(redblue);colorbar;freezeColors;
title('signal corrs')
set(gca,'xtick',[1 intNeurons/2 intNeurons],'ytick',[1 intNeurons/2 intNeurons])
xlabel('Neuron #, sorted by pref. ori.');
ylabel('Neuron #, sorted by pref. ori.');
fixfig;grid off;

subplot(2,3,5)
matSelectLT = logical(tril(ones(intNeurons),-1));
histx(structOut.matSignalCorrs(matSelectLT));
xlim([-1 1]);
title('signal corrs')
ylabel('Number of neuronal pairs');
xlabel('Signal correlation (Pearsons r)');
fixfig;
h = findobj(gca,'Type','patch');
h.EdgeColor = 'none';
h.FaceColor = 'k';

subplot(2,3,3)
imagesc(structOut.matNoiseCorrs,[-1 1]);colormap(redblue);colorbar;freezeColors;
title('noise corrs')
set(gca,'xtick',[1 intNeurons/2 intNeurons],'ytick',[1 intNeurons/2 intNeurons])
xlabel('Neuron #, sorted by pref. ori.');
ylabel('Neuron #, sorted by pref. ori.');
fixfig;grid off;


subplot(2,3,6)
histx(structOut.matNoiseCorrs(matSelectLT));
xlim([-1 1]);
title('noise corrs')
ylabel('Number of neuronal pairs');
xlabel('Noise correlation (Pearsons r)');
fixfig;
h = findobj(gca,'Type','patch');
h.EdgeColor = 'none';
h.FaceColor = 'k';


return

%export figure
export_fig(['D:\Results\V1SubPopulations\simulations\' strSimulation '_PopStatistics.tif']);
export_fig(['D:\Results\V1SubPopulations\simulations\' strSimulation '_PopStatistics.pdf']);
