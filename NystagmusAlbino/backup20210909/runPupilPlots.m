
%get data
strDataSourcePath = 'F:\Data\Processed\Neuropixels\';
sFiles = dir([strDataSourcePath '*.mat']);
cellFiles = {sFiles(:).name}';
close all;

%% pre-allocate
dblFrameDur = 0.0099;
vecPlotWindow = [-1 1.5];
vecBinsT = vecPlotWindow(1):dblFrameDur:vecPlotWindow(2);
matMeanMoveX = nan(numel(vecBinsT),1);
matMeanMoveY = nan(numel(vecBinsT),1);
matMeanRadius = nan(numel(vecBinsT),1);

%% go through files
for intFile=1:numel(cellFiles)
	strTarget = cellFiles{intFile};
	strRec = strTarget((end-10):(end-7));
	sLoad = load(fullfile(strDataSourcePath,strTarget));
	sAP = sLoad.sAP;
	%prepro
	sPupil = sAP.sPupil;
	vecPupilTime = sPupil.vecPupilTime;
	vecPupilFixedRadius = sPupil.vecPupilRadius;
	vecPupilFixedCenterX = sPupil.vecPupilCenterX;
	vecPupilFixedCenterY = sPupil.vecPupilCenterY;
	
	dblLatencyCorrection = -0.25;
	vecPupilTime = vecPupilTime + dblLatencyCorrection;
	dblSamplingRate = 1/median(diff(vecPupilTime));
	sPupil.vecPupilEdgeHardness;
	sPupil.vecPupilMeanPupilLum;
	sPupil.vecPupilAbsVidLum;
	sPupil.vecPupilSdPupilLum;
	sPupil.vecPupilApproxConfidence;
	sPupil.vecPupilApproxRadius;
	sPupil.vecPupilApproxRoundness;
	
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
	
	%%
	if 0
	subplot(2,3,1)
	scatter(sPupil.vecPupilMeanPupilLum,sPupil.vecPupilAbsVidLum,[],sPupil.vecPupilEdgeHardness);
	colormap('redbluepurple');
	
	subplot(2,3,2)
	scatter(zscore(vecFiltX),zscore(vecFiltY),[],sPupil.vecPupilEdgeHardness);
	colormap('redbluepurple');
	end
	
	%%
	%for each stim set
	for intStim=1:numel(sAP.cellStim)
		strType = sAP.cellStim{intStim}.structEP.strFile;
		sOpt = struct;
		if strcmp(strType,'RunNaturalMovie')
			sOpt.vecWindow = [-0.5 10];
			continue;
		elseif strcmp(strType,'RunDriftingGratings')
			sOpt.vecWindow = [-1 1.5];
			vecOrientation = sAP.cellStim{intStim}.structEP.Orientation;
			if numel(unique(vecOrientation)) > 24;continue;end
		elseif strcmp(strType,'')
			strType
			sOpt.vecWindow = [-0.5 1];
			
		else
			strType
			error
		end
		if ~isfield(sAP.cellStim{intStim}.structEP,'vecPupilStimOnTime'),continue;end
		vecPupilStimOnTime = sAP.cellStim{intStim}.structEP.vecPupilStimOnTime;
		vecGlitchesPerTrials = zeros(size(vecPupilStimOnTime));
		intTrials = numel(vecPupilStimOnTime);
		%remove trials
		for intTrial=1:intTrials
			dblStartT = vecPupilStimOnTime(intTrial) + sOpt.vecWindow(1);
			dblStopT = vecPupilStimOnTime(intTrial) + sOpt.vecWindow(2);
			vecGlitchesPerTrials(intTrial) = sum(vecRemTimes > dblStartT & vecRemTimes < dblStopT);
		end
		dblCutOff = 2;
		
		%%
		error
		%plot
		figure;
		hRadius = subplot(3,3,5);
		sOpt.handleFig = hRadius;
		[vecMean,vecSEM,vecWindowBinCenters,matPET] = doPEP(vecPupilTime,vecFiltR,vecPupilStimOnTime(vecGlitchesPerTrials<dblCutOff),sOpt);
		title(sprintf('%sB%d:%s; used %d/%d',strRec,intStim,strType,numel(vecPupilStimOnTime(vecGlitchesPerTrials<dblCutOff)),intTrials));
		ylabel('Pupil radius');
		xlabel('Time after stim onset (s)');
		fixfig;
		
		intPlotOris = 8;
		vecPlotOrder = [6 3 2 1 4 7 8 9];
		vecHandles = nan(size(vecPlotOrder));
		dblOriStep = (pi/(intPlotOris/2));
		for intPlotOri=1:8
			%get which oris to include
			dblCenterOri = (intPlotOri-1)*dblOriStep;
			vecOriRad = deg2rad(vecOrientation);
			indUseOriTrials = abs(circ_dist(dblCenterOri,vecOriRad)) < (dblOriStep*0.49);
			indIncludeTrials = vecGlitchesPerTrials<dblCutOff;
			indUseT = indUseOriTrials & indIncludeTrials;
			%plot
			vecHandles(intPlotOri) = subplot(3,3,vecPlotOrder(intPlotOri));
			sOpt.handleFig = -1;
			[vecMean,vecSEM,vecWindowBinCenters,matPET] = doPEP(vecPupilTime,zscore(vecMoveX),vecPupilStimOnTime(indUseT),sOpt);
			matMeanMoveX(:,end+1) = nan(numel(vecBinsT),1);
matMeanMoveY = nan(numel(vecBinsT),1);
matMeanRadius = nan(numel(vecBinsT),1);

			[vecMean,vecSEM,vecWindowBinCenters,matPET] = doPEP(vecPupilTime,zscore(vecMoveY),vecPupilStimOnTime(indUseT),sOpt);
			title(sprintf('used %d/%d',sum(indUseT),sum(indUseOriTrials)));
			ylabel('Pupil location');
			xlabel('Time after stim onset (s)');
			fixfig;
		end
		maxfig;
	end
end

% plot grand mean
%plot
		figure;
		hRadius = subplot(3,3,5);
		sOpt.handleFig = hRadius;
		[vecMean,vecSEM,vecWindowBinCenters,matPET] = doPEP(vecPupilTime,vecFiltR,vecPupilStimOnTime(vecGlitchesPerTrials<dblCutOff),sOpt);
		title(sprintf('%sB%d:%s; used %d/%d',strRec,intStim,strType,numel(vecPupilStimOnTime(vecGlitchesPerTrials<dblCutOff)),intTrials));
		ylabel('Pupil radius');
		xlabel('Time after stim onset (s)');
		fixfig;
		
		intPlotOris = 8;
		vecPlotOrder = [6 3 2 1 4 7 8 9];
		vecHandles = nan(size(vecPlotOrder));
		dblOriStep = (pi/(intPlotOris/2));
		for intPlotOri=1:8
			%get which oris to include
			dblCenterOri = (intPlotOri-1)*dblOriStep;
			vecOriRad = deg2rad(vecOrientation);
			indUseOriTrials = abs(circ_dist(dblCenterOri,vecOriRad)) < (dblOriStep*0.49);
			indIncludeTrials = vecGlitchesPerTrials<dblCutOff;
			indUseT = indUseOriTrials & indIncludeTrials;
			%plot
			vecHandles(intPlotOri) = subplot(3,3,vecPlotOrder(intPlotOri));
			sOpt.handleFig = vecHandles(intPlotOri);
			sOpt.vecColor = [0 0 1];
			[vecMean,vecSEM,vecWindowBinCenters,matPET] = doPEP(vecPupilTime,zscore(vecMoveX),vecPupilStimOnTime(indUseT),sOpt);
			sOpt.vecColor = [1 0 0];
			[vecMean,vecSEM,vecWindowBinCenters,matPET] = doPEP(vecPupilTime,zscore(vecMoveY),vecPupilStimOnTime(indUseT),sOpt);
			title(sprintf('used %d/%d',sum(indUseT),sum(indUseOriTrials)));
			ylabel('Pupil location');
			xlabel('Time after stim onset (s)');
			fixfig;
		end
		maxfig;
