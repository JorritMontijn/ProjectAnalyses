%prepro
if ~exist('boolSkipSpikes','var') || isempty(boolSkipSpikes)
	boolSkipSpikes = false;
end
sPupil = sAP.sPupil;
vecPupilRadius = sPupil.vecPupilRadius;
vecPupilTime = sPupil.vecPupilTime;
vecPupilFixedRadius = sPupil.vecPupilFixedRadius;
vecPupilFixedCenterX = sPupil.vecPupilFixedCenterX;
vecPupilFixedCenterY = sPupil.vecPupilFixedCenterY;

dblLatencyCorrection = -0.25;
vecPupilTime = vecPupilTime + dblLatencyCorrection;
dblSamplingRate = 1/median(diff(vecPupilTime));

%calc removals
indRem = sPupil.vecPupilEdgeHardness<0.1 | abs(zscore(sPupil.vecPupilCenterX)) > 10 | abs(zscore(sPupil.vecPupilCenterY)) > 10;
vecRemTimes = vecPupilTime(indRem);

%location
vecPupilMovement = [0 sqrt(diff(vecPupilFixedCenterX).^2 + diff(vecPupilFixedCenterY).^2)];

%% high-pass filter to remove drift
dblLowPass = 0.1/dblSamplingRate;
[fb,fa] = butter(2,dblLowPass,'high');
vecFiltX = filtfilt(fb,fa, vecPupilFixedCenterX);
vecFiltY = filtfilt(fb,fa, vecPupilFixedCenterY);
vecFiltR = filtfilt(fb,fa, vecPupilFixedRadius);

%% smooth movement
dblFiltWidth = 0.1;
intSize = round(2*dblFiltWidth*dblSamplingRate);
vecFilt = normpdf(-intSize:intSize,0,intSize/2);
vecFilt = vecFilt ./ sum(vecFilt);
vecMoveX = imfilt([0 diff(vecPupilFixedCenterX)],vecFilt);
vecMoveY = imfilt([0 diff(vecPupilFixedCenterY)],vecFilt);
vecMoveR = imfilt([0 diff(vecPupilFixedRadius)],vecFilt);

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
%for each stim set
for intStim=1:numel(sAP.cellStim)
	strType = sAP.cellStim{intStim}.structEP.strFile;
	sOpt = struct;
	if strcmp(strType,'RunNaturalMovie')
		sOpt.vecWindow = [-0.5 10];
		continue;
	elseif strcmp(strType,'RunDriftingGratings')
		sOpt.vecWindow = [-1 1.5];
		vecSubOrientation = sAP.cellStim{intStim}.structEP.Orientation;
		if numel(unique(vecSubOrientation)) ~= 24;continue;end
	elseif strcmp(strType,'')
		strType
		sOpt.vecWindow = [-0.5 1];
		
	else
		strType
		error
	end
	if ~isfield(sAP.cellStim{intStim}.structEP,'vecPupilStimOnTime'),continue;end
	vecSubPupilStimOnTime = sAP.cellStim{intStim}.structEP.vecPupilStimOnTime;
	vecSubStimOnTime = sAP.cellStim{intStim}.structEP.vecStimOnTime;
	vecSubGlitchesPerTrials = zeros(size(vecSubPupilStimOnTime));
	intTrials = numel(vecSubPupilStimOnTime);
	%remove trials
	for intTrial=1:intTrials
		dblStartT = vecSubPupilStimOnTime(intTrial) + sOpt.vecWindow(1);
		dblStopT = vecSubPupilStimOnTime(intTrial) + sOpt.vecWindow(2);
		vecSubGlitchesPerTrials(intTrial) = sum(vecRemTimes > dblStartT & vecRemTimes < dblStopT);
	end
	
	%concatenate
	vecOrientation = cat(2,vecOrientation,vecSubOrientation);
	vecPupilStimOnTime = cat(2,vecPupilStimOnTime,vecSubPupilStimOnTime);
	vecGlitchesPerTrials = cat(2,vecGlitchesPerTrials,vecSubGlitchesPerTrials);
	vecStimOnTime = cat(2,vecStimOnTime,vecSubStimOnTime);
end

%% get data
cellClustAreas = {sAP.sCluster(:).Area};
cellSpikes = {sAP.sCluster(:).SpikeTimes};
intNeurons = numel(cellSpikes);
vecDecodeWindow = [0 1];
dblCutOff = 2;
vecUseTrialOnsets = vecPupilStimOnTime(vecGlitchesPerTrials<dblCutOff);
intUseTrials = numel(vecUseTrialOnsets);
vecPupilSize = nan(1,intUseTrials);
vecPupilLocX = nan(1,intUseTrials);
vecPupilLocY = nan(1,intUseTrials);
matMeanAct = nan(intNeurons,intUseTrials);
for intTrial=1:intUseTrials
	%pupil
	dblStartPupilT = vecPupilStimOnTime(intTrial) + vecDecodeWindow(1);
	dblStopPupilT = vecPupilStimOnTime(intTrial) + vecDecodeWindow(2);
	vecPupilSize(intTrial) = nanmean(vecFiltR(vecPupilTime > dblStartPupilT & vecPupilTime < dblStopPupilT));
	vecPupilLocX(intTrial) = nanmean(vecFiltX(vecPupilTime > dblStartPupilT & vecPupilTime < dblStopPupilT));
	vecPupilLocY(intTrial) = nanmean(vecFiltY(vecPupilTime > dblStartPupilT & vecPupilTime < dblStopPupilT));
	
	%spikes
	dblStartSpikeT = vecStimOnTime(intTrial) + vecDecodeWindow(1);
	dblStopSpikeT = vecStimOnTime(intTrial) + vecDecodeWindow(2);
	if ~boolSkipSpikes
		parfor intNeuron=1:intNeurons
			matMeanAct(intNeuron,intTrial) = sum(cellSpikes{intNeuron} > dblStartSpikeT & cellSpikes{intNeuron} < dblStopSpikeT);
		end
	end
end
clear boolSkipSpikes;
matMeanAct = matMeanAct./range(vecDecodeWindow);

%get corrected pupil time
mdlr = fitlm(vecPupilStimOnTime,vecStimOnTime,'RobustOpts','on');
beta = mdlr.Coefficients.Estimate;
vecFixedPupilTime = vecPupilTime*beta(2)+beta(1);
%scatter(vecStimOnTime,vecPupilStimOnTime*beta(2)+beta(1)-vecStimOnTime)

%output
vecOrientation;
vecPupilStimOnTime;
vecGlitchesPerTrials;
vecStimOnTime;
vecPupilSize;
vecPupilLocX;
vecPupilLocY;
