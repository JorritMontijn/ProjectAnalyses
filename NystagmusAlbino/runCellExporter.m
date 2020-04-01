%% load spiking data & plot tuning curves
%{
%20190315:
2 (loc1, cort): some possibly visual responses, not very convincing
4 (loc1, subcort), awesome cells, very strong visual responses
7 (loc1, subcort), several strongly visual cells
10 (loc2, cort), multiple strongly visual cells
12 (loc2, subcort, end error), possibly some visual cells, not very strong
%}

%% set recording
clear all;close all;
strSelectMouseType = 'WT';
strSelectArea = 'V1'; 
strSelectStim = 'DG'; %DG=drifting grating, RF=receptive field mapping, RF-DG=merged block
intSelectPopulation	= 1; %V1:1-3, NOT:1, SC:1-5
for intSelectPopulation=1:6
	clearvars -except strSelectMouseType strSelectArea strSelectStim intSelectPopulation
%% load data
loadNystagmusData;

%% run header
runNystagmusHeader;
%{
%% get ISI
for intCell=1:intNumSU
	sOut = getClusterQuality(SU_st{intCell},vecTrialStartTime);
	title(sprintf('Cell %d (orig %d), %sB%s',intCell,vecSingleUnitClusters(intCell),strDate,strBlock));
	pause;
	close;
end

%}
%non-stationarity limit: -0.3 / +0.3
%ISI violation index: both <0.05

%{
V1/WT/DG pop 1, 20190314
vecExclude = [1 7 8 14];

V1/WT/DG pop 2, 20190515
vecExclude = [1 6 10];

V1/WT/DG pop 3, 20190516
vecExclude = [];


SC/WT/DG pop 1, 20190315B4
vecExclude = [171 176 203];

SC/WT/DG pop 2, 20190315B10
vecExclude = [];

SC/WT/DG pop 3, 20190315B7
vecExclude = [85];

SC/WT/DG pop 4, 20190508B11-12
vecExclude = [126? ];

SC/WT/DG pop 5, 20190510B15-16-17
vecExclude = [49 ];

%}

%% calculate Isolation Distance and L_ratio
%https://www.sciencedirect.com/science/article/abs/pii/S0306452204008425?via%3Dihub
%Redish, neuroscience (2005)

%% process eye-tracking
%plot(vecEyeTimestamps,matEyeData(5,:))
%vecEyeTimestamps
%matEyeData

%% get orientations of DG trials
vecOriObjects = find(cellfun(@isfield,structEP.cellStimObject,cellfill('Orientation',size(structEP.cellStimObject))));
indOriTrials = ismember(structEP.ActStimType,vecOriObjects);
indOriTrials(numel(vecStimOnTime):end) = false;
cellOriObjects = structEP.cellStimObject(structEP.ActStimType(indOriTrials));
vecStimOriDegrees = cellfun(@getfield,cellOriObjects,cellfill('Orientation',size(cellOriObjects)));
[vecOriIdx,vecOriTypes,vecReps] = label2idx(vecStimOriDegrees);
vecAngles = deg2rad(vecStimOriDegrees);

%% transform to spikes to response matrix
vecStimOnTime = vecStimOnTime(indOriTrials);
vecStimOffTime = vecStimOffTime(indOriTrials);
vecStimOffTime = vecStimOffTime(1:numel(vecStimOnTime));
dblStimDur = mean(vecStimOffTime - vecStimOnTime);
vecTrialStartTime = vecTrialStartTime(1:numel(vecStimOnTime));
matStimCounts = getSpikeCounts(SU_st,vecStimOnTime,vecStimOffTime);
matStimResp = bsxfun(@rdivide,matStimCounts,(vecStimOffTime-vecStimOnTime)'); %transform to Hz
matBaseCounts = getSpikeCounts(SU_st,vecTrialStartTime,vecStimOnTime);
matBaseResp = bsxfun(@rdivide,matBaseCounts,(vecStimOnTime-vecTrialStartTime)'); %transform to Hz
matResp = matStimResp - matBaseResp;

%% export data
boolSave = true;
if boolSave
	vecSaveSU = [];
	if strcmpi(strSelectArea,'NOT')
		vecSaveSU = [4 7 12 15];
	else
		vecSaveSU = 1:intNumSU;
	end
	for intSaveSU=vecSaveSU
		if ~isempty(intSaveSU)
			runExportNeuronDG;
		end
		clear intSaveSU;
	end
	vecPlotNeurons = vecSaveSU;
else
	vecPlotNeurons = 1:intNumSU;
end
end
if 1,return;end

%% make PSTH
%gather data
close all;
if numel(vecOriTypes) == 8
	intVertPlot = 3;
	intHorPlot = 3;
elseif numel(vecOriTypes) == 20
	intVertPlot = 6;
	intHorPlot = 6;
end
% set plot order; start at right, then counter-clockwise
vecPlotRight = intHorPlot*((intVertPlot-1):-1:2);
vecPlotRightTop = vecPlotRight((floor(end/2)+1):1:end);
vecPlotRightBottom = vecPlotRight(~ismember(vecPlotRight,vecPlotRightTop));
vecPlotTop = intHorPlot:-1:1;
vecPlotLeft = (intHorPlot+1):intHorPlot:(intHorPlot*(intVertPlot-2)+1);
vecPlotBottom = ((intVertPlot-1)*intHorPlot+1):1:(intVertPlot*intHorPlot);
vecUseSubPlot = cat(2,vecPlotRightTop,vecPlotTop,vecPlotLeft,vecPlotBottom,vecPlotRightBottom);
vecMiddle = find(~ismember(1:(intVertPlot*intHorPlot),vecUseSubPlot));

sOptions.handleFig = -1;
dblStep = 0.1;
vecWindow = -0.5:dblStep:2.5;
vecWindowBinCenters = vecWindow(2:end)-dblStep/2;
indUseStimBins = vecWindowBinCenters > 0 & vecWindowBinCenters < dblStimDur;
matPEP = nan(numel(vecWindow)-1,numel(vecOriTypes),2,intNumSU); %[time x stim x mean/sd x neuron]
for intSU = vecPlotNeurons
	%% plot neuron
	vecSpikeTimes = SU_st{intSU};
	figure;
	if numel(vecSpikeTimes) < 50,continue;end
	%get data
	for intStimType=1:numel(vecOriTypes)
		vecEvents = vecStimOnTime(vecOriIdx==intStimType);
		[vecMean,vecSEM] = doPEP(vecSpikeTimes,vecWindow,vecEvents,sOptions);
		matPEP(:,intStimType,1,intSU) = vecMean;
		matPEP(:,intStimType,2,intSU) = vecSEM;
	end
	
	%plot
	dblMax = max(flat(matPEP(:,:,1,intSU) + matPEP(:,:,2,intSU)));
	dblMin = min(flat(matPEP(:,:,1,intSU) - matPEP(:,:,2,intSU)));
	vecLimY = [min([0 dblMin]) dblMax];
	for intStimType=1:numel(vecOriTypes)
		subplot(intVertPlot,intHorPlot,vecUseSubPlot(intStimType))
		errorfill(vecWindowBinCenters,matPEP(:,intStimType,1,intSU),matPEP(:,intStimType,2,intSU));
		ylim(vecLimY);
		xlabel('Time from stim on (s)');
		ylabel('Spiking rate (Hz)');
		title(sprintf('Ori %d (%.1f degs)',intStimType,vecOriTypes(intStimType)));
	end
	
	%plot radial tuning curve
	vecMeanResp = mean(matPEP(indUseStimBins,:,1,intSU),1);
	subplot(intVertPlot,intHorPlot,vecMiddle);
	polarplot(deg2rad(vecOriTypes([1:end 1])),vecMeanResp([1:end 1]));
	title(sprintf('SU %d (orig %d); depth=%.1f%sm (ch %.1f)',intSU,vecSingleUnitClusters(intSU),...
		SU_depth_micron(intSU),getGreek(12,'lower'),SU_depth_ch(intSU)));
	%pause
end
 
%plot across neurons
matMeanAcrossN = mean(matPEP(:,:,1,:),4); %[time x stim]
matSEMAcrossN = std(matPEP(:,:,1,:),[],4)./sqrt(intNumSU); %[time x stim]
%plot
dblMax = max(flat(matMeanAcrossN + matSEMAcrossN));
dblMin = min(flat(matMeanAcrossN - matSEMAcrossN));
vecLimY = [min([0 dblMin]) dblMax];
figure
for intStimType=1:numel(vecOriTypes)
	subplot(3,3,intStimType)
	errorfill(vecWindowBinCenters,matMeanAcrossN(:,intStimType),matSEMAcrossN(:,intStimType));
	ylim(vecLimY);
	xlabel('Time from stim on (s)');
	ylabel('Spiking rate (Hz)');
	title(sprintf('Mean across %d units; ori %d (%.1f degs)',intNumSU,intStimType,vecOriTypes(intStimType)));
end

%define NOT comparison
vecLeftRight = [0 180];
return
%% save figure
strFigPath = 'D:\Data\ResultsNystagmus\';
strFig = strcat('ExampleCell_',strSelectArea,'_',strMouse,'_',strDate,'B',strBlock);
export_fig([strFigPath strFig '.tif'])
export_fig([strFigPath strFig '.pdf'])