%function getPupilPlotsNoEphys(sPupil,cellStim,cellRecNames,sFiles)
%clearvars -except sPupil cellStim cellRecNames sFiles strDir

%% prepro
%set output location
boolMakeFigure = true;
strTargetDir = 'F:\Data\Results\EyeTracking\';
strFile = strrep(sPupil.name,'Processed','Results');

%raw and control vals
vecPupilStimOn = sPupil.vecPupilStimOn;
vecPupilTime = sPupil.vecPupilTime;
vecPupilAbsVidLum = sPupil.vecPupilAbsVidLum;
vecPupilEdgeHardness = sPupil.vecPupilEdgeHardness;
vecPupilSyncLum = sPupil.vecPupilSyncLum;
vecPupilCenterX = sPupil.vecPupilCenterX;
vecPupilCenterY = sPupil.vecPupilCenterY;

%fixed vals
%vecPupilFixedRadius = sPupil.vecPupilFixedRadius;
%vecPupilFixedCenterX = sPupil.vecPupilFixedCenterX;
%vecPupilFixedCenterY = sPupil.vecPupilFixedCenterY;
vecPupilFixedRadius = sPupil.vecPupilRadius;
vecPupilFixedCenterX = sPupil.vecPupilCenterX;
vecPupilFixedCenterY = sPupil.vecPupilCenterY;

dblLatencyCorrection = 0;
dblDetectionRate = 1/median(diff(vecPupilTime));

%calc removals
indRem = vecPupilEdgeHardness<0.1 | abs(zscore(vecPupilCenterX)) > 10 | abs(zscore(vecPupilCenterY)) > 10;
vecRemTimes = vecPupilTime(indRem);

%location
vecPupilMovement = [0 sqrt(diff(vecPupilFixedCenterX).^2 + diff(vecPupilFixedCenterY).^2)];

%% high-pass filter to remove drift
dblLowPass = 0.1/dblDetectionRate;
[fb,fa] = butter(2,dblLowPass,'high');
vecFiltX = filtfilt(fb,fa, vecPupilFixedCenterX);
vecFiltY = filtfilt(fb,fa, vecPupilFixedCenterY);
vecFiltR = filtfilt(fb,fa, vecPupilFixedRadius);
vecFiltAbsVidLum = filtfilt(fb,fa, vecPupilAbsVidLum);
[boolBlinking,dblCritValBlink]=DP_GetUpDown(vecFiltAbsVidLum,0.9,0.995);

%spread the blink
boolBlinking = conv(boolBlinking,ones(1,5),'same')>0;

%find longest period without stimuli
vecFiltSyncLum = filtfilt(fb,fa, vecPupilSyncLum);
[boolStimStart,dblCritValStim]=DP_GetUpDown(vecFiltSyncLum,0.5,1);
vecStimT = [1 find(boolStimStart) numel(boolStimStart)];
vecInterStimT = diff(vecStimT);
[dblDur,intStartIdx]=max(vecInterStimT);
dblStartBlank = vecPupilTime(vecStimT(intStartIdx))+5;
dblStopBlank = vecPupilTime(vecStimT(intStartIdx+1))-5;
dblStartBlankIdx = find(vecPupilTime>dblStartBlank,1);
dblStopBlankIdx = find(vecPupilTime>dblStopBlank,1);

for intRec=1:numel(cellStim)
	if isfield(cellStim{intRec}.structEP,'ActOnPupil')
		cellAllPupilOn{intRec} = cellStim{intRec}.structEP.ActOnPupil;
	else
		cellAllPupilOn{intRec} = [];
	end
end
vecAllPupilOn = cell2vec(cellAllPupilOn);

%plot blink
%{
plot(vecPupilTimeNI,vecFiltAbsVidLum)
hold on
scatter(vecPupilTimeNI(find(boolBlinking==1)),dblCritValBlink*ones(1,sum(boolBlinking)));
scatter(vecPupilTimeNI(find(indRem==1)),dblCritValBlink*ones(1,sum(indRem)));
hold off
%}

%plot sync
%{
plot(vecPupilTimeNI,vecFiltSyncLum)
hold on
plot([dblStartBlankNI dblStopBlankNI],dblCritValStim*[1 1]);
hold off
%}
%% smooth movement
dblFiltWidth = 0.1;
intSize = round(2*dblFiltWidth*dblDetectionRate);
vecFilt = normpdf(-intSize:intSize,0,intSize/2);
vecFilt = vecFilt ./ sum(vecFilt);
vecMoveX = abs(imfilt([0 diff(vecPupilFixedCenterX)],vecFilt));
vecMoveY = abs(imfilt([0 diff(vecPupilFixedCenterY)],vecFilt));
vecMoveR = abs(imfilt([0 diff(vecPupilFixedRadius)],vecFilt));

%% set unlikely values to nan
vecFiltX(boolBlinking | indRem) = nan;
vecFiltY(boolBlinking | indRem) = nan;
vecFiltR(boolBlinking | indRem) = nan;
vecMoveX(boolBlinking | indRem) = nan;
vecMoveY(boolBlinking | indRem) = nan;
vecMoveR(boolBlinking | indRem) = nan;

%% retrieve pupil parameters outside stimuli
vecFiltX_Base = vecFiltX(dblStartBlankIdx:dblStopBlankIdx);
vecFiltY_Base = vecFiltY(dblStartBlankIdx:dblStopBlankIdx);
vecFiltR_Base = vecFiltR(dblStartBlankIdx:dblStopBlankIdx);
vecMoveX_Base = vecMoveX(dblStartBlankIdx:dblStopBlankIdx);
vecMoveY_Base = vecMoveY(dblStartBlankIdx:dblStopBlankIdx);
vecMoveR_Base = vecMoveR(dblStartBlankIdx:dblStopBlankIdx);

vecLocX_Norm = (vecFiltX - nanmean(vecFiltX_Base)) ./ nanstd(vecFiltX_Base);
vecLocY_Norm = (vecFiltY - nanmean(vecFiltY_Base)) ./ nanstd(vecFiltY_Base);
vecLocR_Norm = (vecFiltR - nanmean(vecFiltR_Base)) ./ nanstd(vecFiltR_Base);
vecMoveX_Norm = (vecMoveX - nanmean(vecMoveX_Base)) ./ nanstd(vecMoveX_Base);
vecMoveY_Norm = (vecMoveY - nanmean(vecMoveY_Base)) ./ nanstd(vecMoveY_Base);
vecMoveR_Norm = (vecMoveR - nanmean(vecMoveR_Base)) ./ nanstd(vecMoveR_Base);

%%
if 0
	subplot(2,3,1)
	scatter(sPupil.vecPupilMeanPupilLum,sPupil.vecPupilAbsVidLum,[],sPupil.vecPupilEdgeHardness);
	colormap('redbluepurple');
	
	subplot(2,3,2)
	scatter(zscore(vecFiltX),zscore(vecFiltY),[],sPupil.vecPupilEdgeHardness);
	colormap('redbluepurple');
end

%% concatenate drifting grating trials
vecOrientation = [];
vecPupilStimOnTime = [];
vecGlitchesPerTrials = [];
vecStimOnTime = [];
intCounter = 0;
dblSampFreq = median(diff(vecPupilTime));
dblT0=-1;
dblTotDur = 3;
vecRefT = 0:dblSampFreq:(dblTotDur+dblSampFreq/2);
intBins = numel(vecRefT);
dblRealDur = range(vecRefT)+dblSampFreq/2;
vecInterpT = -0.5:0.01:1.5;
vecLimX = [vecInterpT(1) vecInterpT(end)];

%for each stim set
for intStim=1:numel(cellStim)
	if isfield(cellStim{intStim}.structEP,'strExpType')
		strType = cellStim{intStim}.structEP.strExpType;
	else
		strType = cellStim{intStim}.structEP.strFile;
	end
	sOpt = struct;
	if strcmp(strType,'RunNaturalMovie') || strcmp(strType,'RunReceptiveFieldMapping')
		sOpt.vecWindow = [-0.5 10];
		continue;
	elseif strcmp(strType,'RunDriftingGratings')
		sOpt.vecWindow = [-1 1.5];
		vecSubOrientation = cellStim{intStim}.structEP.Orientation;
		if numel(unique(vecSubOrientation)) ~= 24;continue;end
	elseif strcmp(strType,'')
		strType
		sOpt.vecWindow = [-0.5 1];
		
	else
		strType
		error
	end
	
	%build output name
	cellFile = strsplit(strFile,'.');
	strDataFile = [strjoin(cellFile(1:(end-1)),'.') '_Stim' num2str(intStim)];
	
	%assign time vector
	vecOnTime = cellStim{intStim}.structEP.ActOnPupil;
	vecOffTime = vecOnTime - cellStim{intStim}.structEP.ActOnSecs + cellStim{intStim}.structEP.ActOffSecs;
	if vecOnTime(end) < vecPupilTime(1) || vecOnTime(1) > vecPupilTime(end)
		continue;
	end
	%location
	[dummy,matHorzLocPerTrial] = getTraceInTrial(vecPupilTime,vecLocX_Norm,vecOnTime+dblT0,dblSampFreq,dblRealDur);
	[dummy,matVertLocPerTrial] = getTraceInTrial(vecPupilTime,vecLocY_Norm,vecOnTime+dblT0,dblSampFreq,dblRealDur);
	[dummy,matRadiusLocPerTrial] = getTraceInTrial(vecPupilTime,vecLocR_Norm,vecOnTime+dblT0,dblSampFreq,dblRealDur);
	%movement
	[dummy,matHorzMovePerTrial] = getTraceInTrial(vecPupilTime,vecMoveX_Norm,vecOnTime+dblT0,dblSampFreq,dblRealDur);
	[dummy,matVertMovePerTrial] = getTraceInTrial(vecPupilTime,vecMoveY_Norm,vecOnTime+dblT0,dblSampFreq,dblRealDur);
	[dummy,matRadiusMovePerTrial] = getTraceInTrial(vecPupilTime,vecMoveR_Norm,vecOnTime+dblT0,dblSampFreq,dblRealDur);
	[vecTraceT,matSyncLumPerTrial] =  getTraceInTrial(vecPupilTime,vecFiltSyncLum,vecOnTime+dblT0,dblSampFreq,dblRealDur);
	vecRefT = vecTraceT + dblT0;
	
	%delay between matlab and LED
	vecMeanSyncLum = mean(matSyncLumPerTrial,1);
	[boolSyncLum,dblCritValSL] = DP_GetUpDown(vecMeanSyncLum);
	dblDelayCorrection = vecRefT(find(diff(boolSyncLum)==1,1));
	vecCorrT = vecRefT - dblDelayCorrection;
	
	%vecFiltX
	%vecFiltY
	%vecFiltR
	
	%interpolate to standard time
	intTrials = size(matHorzLocPerTrial,1);
	intNumT = numel(vecInterpT);
	matX = nan(intTrials,intNumT);
	matY = nan(intTrials,intNumT);
	matR = nan(intTrials,intNumT);
	matXderiv = nan(intTrials,intNumT);
	matYderiv = nan(intTrials,intNumT);
	matRderiv = nan(intTrials,intNumT);
	matL = nan(intTrials,intNumT);
	for intTrial=1:intTrials
		matX(intTrial,:) = interp1(vecCorrT,matHorzLocPerTrial(intTrial,:),vecInterpT);
		matY(intTrial,:) = interp1(vecCorrT,matVertLocPerTrial(intTrial,:),vecInterpT);
		matR(intTrial,:) = interp1(vecCorrT,matRadiusLocPerTrial(intTrial,:),vecInterpT);
		matXderiv(intTrial,:) = interp1(vecCorrT,matHorzMovePerTrial(intTrial,:),vecInterpT);
		matYderiv(intTrial,:) = interp1(vecCorrT,matVertMovePerTrial(intTrial,:),vecInterpT);
		matRderiv(intTrial,:) = interp1(vecCorrT,matRadiusMovePerTrial(intTrial,:),vecInterpT);
		matL(intTrial,:) = interp1(vecCorrT,matSyncLumPerTrial(intTrial,:),vecInterpT);
	end
	
	%plot
	intCounter = intCounter + 1;
	indUseTrials = ~all(isnan(matR),2);
	intTrials = sum(indUseTrials);
	vecMeanRad = nanmean(matR(indUseTrials,:),1);
	if numel(vecMeanRad) < intBins,vecMeanRad((end+1):intBins)=vecMeanRad(end);end
	matMeanRadius(1:intBins,intCounter) = vecMeanRad(1:intBins);
	
	if boolMakeFigure
	h=figure;
	maxfig;
	hAx=gca;
	plot(hAx,vecRefT,matSyncLumPerTrial');
	hAx.ColorOrder = redbluepurple;
	subplot(3,3,5);
	hold on;
	plot(vecInterpT,zscore(mean(matL,1))/50,'k');
	errorbar(vecInterpT, vecMeanRad, nanstd(matR,[],1)/sqrt(intTrials),'color',lines(1))
	hold off
	title(sprintf('Blue=radius,Black=sync'));
	fixfig;
	xlim(vecLimX);
	
	intPlotOris = 8;
	vecPlotOrder = [6 3 2 1 4 7 8 9];
	vecHandles = nan(size(vecPlotOrder));
	dblOriStep = (pi/(intPlotOris/2));
	for intPlotOri=1:8
		%get which oris to include
		dblCenterOri = (intPlotOri-1)*dblOriStep;
		vecOriRad = deg2rad(vecSubOrientation);
		indUseOriTrials = abs(circ_dist(dblCenterOri,vecOriRad)) < (dblOriStep*0.49);
		%indIncludeTrials = vecGlitchesPerTrials<dblCutOff;
		indUseT = indUseOriTrials(:) & indUseTrials(:);
		%plot
		
		subplot(3,3,vecPlotOrder(intPlotOri));
		hold on;
		matXZ = (matXderiv - nanmean(matXderiv(:)))./nanstd(matXderiv(:));
		matYZ = (matYderiv - nanmean(matYderiv(:)))./nanstd(matYderiv(:));
		errorbar(vecInterpT, nanmean(matXZ(indUseT,:),1), nanstd(matXZ(indUseT,:),[],1)/sqrt(intTrials),'Color','r');
		errorbar(vecInterpT, nanmean(matYZ(indUseT,:),1), nanstd(matYZ(indUseT,:),[],1)/sqrt(intTrials),'Color','b');
		hold off
		fixfig;
		title(sprintf('X=r,Y=b,Ori=%.1f deg',rad2deg(dblCenterOri)));
		xlim(vecLimX);
		
		continue;
		if numel(vecMean) < intBins,vecMean((end+1):intBins)=vecMean(end);end
		matMeanMoveX(1:intBins,intPlotOri,intCounter) = vecMean(1:intBins);
		
		if numel(vecMean) < intBins,vecMean((end+1):intBins)=vecMean(end);end
		matMeanMoveY(1:intBins,intPlotOri,intCounter) = vecMean(1:intBins);
	end
	
	%save figure
	drawnow;
	export_fig([strTargetDir strDataFile 'Loc.tif']);
	export_fig([strTargetDir strDataFile 'Loc.pdf']);
	end
	%% save data
	% assign output
	sPupilResults = struct;
	if isfield(cellStim{intStim}.structEP.sStimParams,'strRecording')
		sPupilResults.strRec = cellStim{intStim}.structEP.sStimParams.strRecording;
	else
		sPupilResults.strRec = '';
	end
	sPupilResults.strBlock = cellRecNames{intStim};
	sPupilResults.strVidFile = sPupil.strVidFile;
	sPupilResults.strVidPath = sPupil.strVidPath;
	sPupilResults.strProcFile = sPupil.name;
	sPupilResults.Orientation = cellStim{intStim}.structEP.Orientation;
	sPupilResults.vecT = vecInterpT;
	sPupilResults.matX = matX;
	sPupilResults.matY = matY;
	sPupilResults.matR = matR;
	sPupilResults.matXderiv = matXderiv;
	sPupilResults.matYderiv = matYderiv;
	sPupilResults.matRderiv = matRderiv;
	sPupilResults.matL = matL;
	sPupilResults.structEP = cellStim{intStim}.structEP;
	sPupilResults.sPupil = sPupil;
	
	save([strTargetDir strDataFile '.mat'],'sPupilResults');
end
