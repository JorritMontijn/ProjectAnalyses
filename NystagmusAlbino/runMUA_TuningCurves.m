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
clear all;%close all;
strSelectMouseType = 'WT';
strSelectArea = 'V1';
strSelectStim = 'DG'; %DG=drifting grating, RF=receptive field mapping, RF-DG=merged block
intSelectPopulation	= 2;
strFigPath = 'D:\Data\ResultsNystagmus\';

%% load data
loadNystagmusData;

%% run header
runNystagmusHeader;

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
intTrials = sum(indOriTrials);

%% transform to spikes to response matrix
vecStimOnTime = vecStimOnTime(indOriTrials);
vecStimOffTime = vecStimOffTime(indOriTrials);
vecStimOffTime = vecStimOffTime(1:numel(vecStimOnTime));
dblStimDur = mean(vecStimOffTime - vecStimOnTime);
vecTrialStartTime = vecTrialStartTime(1:numel(vecStimOnTime));

%% calc OT estimate
%update variables
[vecOriIdx,vecOriTypes,vecReps,cellSelect,vecRepetition] = label2idx(vecStimOriDegrees);
vecAngles = deg2rad(vecStimOriDegrees);

%prep window
intChNr = size(matEnvData,1);
dblStep = 0.2;
vecWindow = -0.5:dblStep:2.5;
vecWindowBinCenters = vecWindow(2:end)-dblStep/2;
indUseStimBins = vecWindowBinCenters > 0 & vecWindowBinCenters < dblStimDur;
matPEP = nan(numel(vecWindow)-1,numel(vecOriTypes),2,intChNr); %[time x stim x mean/sd x neuron]
sOptions.handleFig = -1;

%% transform data to new time steps
vecBinnedTimestamps = 0:dblStep:ceil(vecEnvTimestamps(end));
vecX = vecEnvTimestamps';
vecY = single(matEnvData)';
vecBins = vecBinnedTimestamps';
[vecCounts,matDataResampled,vecSDs,cellVals,cellIDs] = makeBins(vecX,vecY,vecBins);

%base, stim
matRespBase = nan(intChNr,intTrials);
matRespStim = nan(intChNr,intTrials);
vecStimTypes = nan(1,intTrials);
vecStimOriDeg = nan(1,intTrials);
%go through objects and assign to matrices

%get data
for intStimType=1:numel(vecOriTypes)
	vecEvents = vecStimOnTime(vecOriIdx==intStimType);
	intReps = numel(vecEvents);
	matThisData = nan(size(matPEP,1),intReps,intChNr);
	for intRep=1:intReps
		dblEventSecs = vecEvents(intRep);
		vecSelectSecs = vecWindow + dblEventSecs;
		
		intStart = find(vecBinnedTimestamps>vecSelectSecs(1),1);
		intStop = find(vecBinnedTimestamps>vecSelectSecs(end),1)-1;
		if isempty(intStop),intStop=numel(vecBinnedTimestamps);end
		
		
		vecSelectEntries = intStart:intStop;
		indKeep = ~(vecSelectEntries<1 | vecSelectEntries>size(matDataResampled,1));
		vecSelectBins = vecSelectEntries(indKeep);
		
		for intCh=1:intChNr
			matThisData(:,intRep,intCh) = matDataResampled(vecSelectBins,intCh);
		end
	end
	
	
	for intCh=1:intChNr
		
		matPEP(:,intStimType,1,intCh) = nanmean(matThisData(:,:,intCh),2);
		matPEP(:,intStimType,2,intCh) = nanstd(matThisData(:,:,intCh),[],2)/sqrt(intReps);
	end
end
vecMeanPerChannel = nanmean(matDataResampled,1);
vecSDPerChannel = nanstd(matDataResampled,[],1);
%% get tuning
%[sOut] = getTuningCurves(matStimResp,vecStimOriDegrees);
%vecDeltaPrime = getDeltaPrime(matResp,vecAngles);

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

%get all-channel average
vecShiftDegs = circshift(1:numel(vecOriTypes),[0 round(size(vecOriTypes,2)/2)]);
vecAverage = mean(mean(matPEP(indUseStimBins,:,1,:),1),4);
matMeanOriCh = nan(numel(vecOriTypes),intChNr);
boolManyPlots = false;
for intCh = 1:intChNr
	%calculate mean response
	vecMeanResp = mean(matPEP(indUseStimBins,:,1,intCh),1);
	vecMeanResp = (vecMeanResp - vecMeanPerChannel(intCh)) ./ vecSDPerChannel(intCh);
	matMeanOriCh(:,intCh) = vecMeanResp;
	if boolManyPlots
	%make figure
	figure;
	
	%plot
	dblMax = max(flat(matPEP(:,:,1,intCh) + matPEP(:,:,2,intCh)));
	dblMin = min(flat(matPEP(:,:,1,intCh) - matPEP(:,:,2,intCh)));
	vecLimY = [dblMin dblMax];
	for intStimType=1:numel(vecOriTypes)
		subplot(intVertPlot,intHorPlot,vecUseSubPlot(intStimType))
		errorfill(vecWindowBinCenters,matPEP(:,intStimType,1,intCh),matPEP(:,intStimType,2,intCh));
		ylim(vecLimY);
		xlabel('Time from stim on (s)');
		ylabel('MUA (a.u.)');
		title(sprintf('Ch%02d; ori %d (%.1f degs)',intCh,intStimType,vecOriTypes(intStimType)));
	end
	
	%plot radial tuning curve
	subplot(intVertPlot,intHorPlot,vecMiddle);
	
	%polar plot?
	%polarplot(deg2rad(vecOriTypes([1:end 1])),vecMeanResp([1:end 1]) - min(vecMeanResp(:)));
	
	%or regular plot?
	plot(vecOriTypes,vecMeanResp(vecShiftDegs));
	ylim([-0.5 1]);
	xlim([0 360]);
	vecUseTicks = get(gca,'xtick');
	set(gca,'xtick',(0:90:360),'xticklabel',mod((0:90:360)+180,360));
	%pause
	end
end
 
%plot across channels
matMeanAcrossN = mean(matPEP(:,:,1,:),4); %[time x stim]
matSEMAcrossN = std(matPEP(:,:,1,:),[],4)./sqrt(intNumSU); %[time x stim]
%plot
dblMax = max(flat(matMeanAcrossN + matSEMAcrossN));
dblMin = min(flat(matMeanAcrossN - matSEMAcrossN));
vecLimY = [min([0 dblMin]) dblMax];
figure
for intStimType=1:numel(vecOriTypes)
	subplot(intVertPlot,intHorPlot,vecUseSubPlot(intStimType))
	errorfill(vecWindowBinCenters,matMeanAcrossN(:,intStimType),matSEMAcrossN(:,intStimType));
	ylim(vecLimY);
	xlabel('Time from stim on (s)');
	ylabel('MUA (a.u.)');
	title(sprintf('Mean across %d units; ori %d (%.1f degs)',intNumSU,intStimType,vecOriTypes(intStimType)));
end

%plot radial tuning curve
vecMeanMeanResp = mean(matMeanAcrossN,1);
subplot(intVertPlot,intHorPlot,vecMiddle);

%polar plot?
%polarplot(deg2rad(vecOriTypes([1:end 1])),vecMeanResp([1:end 1]) - min(vecMeanResp(:)));

%or regular plot?
plot(vecOriTypes,vecMeanMeanResp(vecShiftDegs));
%ylim([-0.5 1]);
xlim([0 360]);
vecUseTicks = get(gca,'xtick');
set(gca,'xtick',(0:90:360),'xticklabel',mod((0:90:360)+180,360));
%pause

%define NOT comparison
vecLeftRight = [0 180];

%% make heat map of MUA as function of depth and orientation
figure
vecMedian = median(matMeanOriCh,2);
matMeanNormOriCh = bsxfun(@minus,matMeanOriCh,vecMedian);
imagesc(1:numel(vecOriTypes),vecChannelDepth,matMeanNormOriCh(vecShiftDegs,:)',[-0.1 0.8]);colormap(hot)
set(gca,'xtick',linspace(0.5,numel(vecOriTypes)+0.5,5),'xticklabel',mod((0:90:360)+180,360));
xlabel('Stimulus orientation (degs)');
ylabel('Channel depth from dura (micron)');
title(sprintf('Norm. MUA, %s_%s, B%s (%s)',strMouse,strDate,strBlock,strSelectArea),'interpreter','none');
fixfig;
colorbar;
grid off;
drawnow;

%% save figure
strFig = strcat(strSelectArea,'_',strMouse,'_',strDate,'B',strBlock);
export_fig([strFigPath strFig '.tif'])
export_fig([strFigPath strFig '.pdf'])