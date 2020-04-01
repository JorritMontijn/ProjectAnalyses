

%% get simulation name [strSimulation] from [intLoadSim]
clear all;close all;
strSimulation = 'Simulation_xAreaDistributed_SG18_2019-07-04';
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

%% get stimulus types
vecUniqueStimTypes = unique(vecTrialStimType);
intStimTypes = numel(vecUniqueStimTypes);
 
%% build synchronization vectors
dblBinSize_ms = 1;
dblJitter_ms = 2;
dblJitterSecs = dblJitter_ms/1000;
dblBinSizeSecs = dblBinSize_ms/1000;
%build timestamps
dblFirstSpike = min(cellfun(@min,cellSpikeTimesCortex));
dblLastSpike = max(cellfun(@max,cellSpikeTimesCortex));
dblStart = dblFirstSpike - dblBinSizeSecs - dblJitterSecs;
dblStop = dblLastSpike + dblBinSizeSecs + dblJitterSecs*2;
vecBinEdges = dblStart:dblBinSizeSecs:dblStop;
vecTimestamps = vecBinEdges(2:end) - dblBinSizeSecs/2;

%get neuron numbers
intNeuronsV1 = sum(vecCellArea==1);
intNeuronsV2 = sum(vecCellArea==2);
%V1
[vecPopSyncV1,vecTimestampsV1] = getPopulationSynchronization(cellSpikeTimesCortex(vecCellArea==1),dblJitter_ms,vecTimestamps);
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
vecSpikesV1 = cell2mat(cellSpikeTimesCortex(vecCellArea==1));
vecSpikesV2 =  cell2mat(cellSpikeTimesCortex(vecCellArea==2));

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

%%
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
%transform sync epochs 
%% transform to Hz
vecWinSpikeHzV1 = (vecWinSpikeCountV1 / dblWindowSecs) / intNeuronsV1;
vecWinSpikeHzV2 = (vecWinSpikeCountV2 / dblWindowSecs) / intNeuronsV2;
vecWinSyncHzV1 = (((vecWinSyncEpochDurV1)*(dblBinSizeSecs/dblWindowStepSecs))/dblWindowSecs) / intNeuronsV1;

%% prep counts
%bin spike V1
dblMaxCountV1 = mean(vecWinSpikeHzV1) + 3*std(vecWinSpikeHzV1);
dblMinCountV1 = 0;
dblStepCountV1 = 1;
vecEdgesCountV1 = dblMinCountV1:dblStepCountV1:dblMaxCountV1;

%bin sync V1
[vecUniques,dummy,vecValIdx]=unique(vecWinSyncHzV1);
intUniques = numel(vecUniques);
vecValCounts = zeros(1,intUniques);
for intU=1:intUniques
	vecValCounts(intU) = sum(c==intU);
end
vecCutOffCount = 100;
intCUV = find(vecValCounts < vecCutOffCount);
if isempty(intCUV),intCUV = intUniques;end
dblSyncStep = vecUniques(2) - vecUniques(1);
vecEdgesSyncV1 = (vecUniques(1) - dblSyncStep/2):dblSyncStep:(vecUniques(intCUV) + dblSyncStep/2);

%define x,y,z
vecValsX = vecWinSpikeHzV1;
vecValsY = vecWinSyncHzV1;
vecValsZ = vecWinSpikeHzV2;

%define edges x,y,z
vecEdgesX = vecEdgesCountV1;
vecEdgesY = vecEdgesSyncV1;

%% get counts
boolCellVals = false;
intValNumX = (numel(vecEdgesX)-1);
intValNumY = (numel(vecEdgesY)-1);
matCounts = zeros(intValNumY,intValNumX);
matValMeans = zeros(intValNumY,intValNumX);
matValSDs = zeros(intValNumY,intValNumX);
cellVals = cell(intValNumY,intValNumX);
for intEdgeX=1:intValNumX
	dblStartEdgeX = vecEdgesX(intEdgeX);
	dblStopEdgeX = vecEdgesX(intEdgeX+1);
	%pre-select X bin
	indTheseValsX = vecValsX > dblStartEdgeX & vecValsX < dblStopEdgeX;
	vecTheseX = vecValsX(indTheseValsX);
	vecTheseY = vecValsY(indTheseValsX);
	vecTheseZ = vecValsZ(indTheseValsX);
	
	for intEdgeY=1:intValNumY
		dblStartEdgeY = vecEdgesY(intEdgeY);
		dblStopEdgeY = vecEdgesY(intEdgeY+1);
		
		%get values for bin
		indTheseValsXY = vecTheseY > dblStartEdgeY & vecTheseY < dblStopEdgeY;
		vecTheseZValsXY = vecTheseZ(indTheseValsXY);
		
		%assign
		matCounts(intEdgeY,intEdgeX) = sum(indTheseValsXY);
		matValMeans(intEdgeY,intEdgeX) = mean(vecTheseZValsXY);
		matValSDs(intEdgeY,intEdgeX) = std(vecTheseZValsXY);
		if boolCellVals,cellVals{intEdgeY,intEdgeX} = vecTheseZValsXY;end
	end
end
%get counts
[N,Xedges,Yedges,binX,binY] = histcounts2(vecValsX,vecValsY,vecEdgesX,vecEdgesY);

%% plot heat map
close all
figure
vecPlotX = vecEdgesX(2:end)-median(diff(vecEdgesX))/2;
vecPlotY = vecEdgesY(2:end)-median(diff(vecEdgesY))/2;
matPlotZ = matValMeans;
matPlotZ(matCounts<250) = nan;

subplot(2,2,2)
matC = colormap(parula(numel(vecPlotX)));
hold all;
for intL=1:size(matC,1)
	plot(vecPlotY,matPlotZ(:,intL),'Color',matC(intL,:));
end
hold off;
xlabel('Synchronization rate (Hz)');
ylabel('V2 spiking rate (Hz)');
title('Lines: various V1 spiking rate bins');
fixfig;

subplot(2,2,3)
matC = colormap(parula(numel(vecPlotY)));
hold all;
for intL=1:size(matC,1)
	plot(vecPlotX,matPlotZ(intL,:),'Color',matC(intL,:));
end
hold off
ylabel('V2 spiking rate (Hz)');
xlabel('Synchronization rate (Hz)');
title('Lines: various synchronization rate bins');
fixfig;

subplot(2,2,1)
imagesc(vecPlotX,vecPlotY,matPlotZ);
axis xy;
nancolorbar(matPlotZ,[],parula(256));
title(sprintf('Average V2 spiking rate (Hz); window: %dms',dblWindow_ms));
xlabel('Average V1 spiking rate (Hz)');
ylabel('Average sync. rate per neuron (Hz)');
fixfig;
grid off;

%% save
strFile = sprintf('NeuralMarker3_V1Synchrony_Window%dms',round(dblWindow_ms));
export_fig([strFile '.tif']);
export_fig([strFile '.pdf']);
return
%% run analysis to check per V2 spike per cell what the preceding and following spikes in V1 and V2 are
% what is the ratio in predictability of V1 and V2?
% are population events in V1 more common during V2 spikes?
% are there specific temporal patterns of activation?

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