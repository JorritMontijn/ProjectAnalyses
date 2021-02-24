
%get data
strDataSourcePath = 'F:\Data\Processed\Neuropixels\';
sFiles = dir([strDataSourcePath '*.mat']);
cellFiles = {sFiles(:).name}';
close all;

%% pre-allocate
dblFrameDur = 0.0099;
intPlotOris = 8;
vecPlotOrder = [6 3 2 1 4 7 8 9];
intCounter = 0;
vecPlotWindow = [-1 1.5];
vecBinsT = vecPlotWindow(1):dblFrameDur:vecPlotWindow(2);
intBins = numel(vecBinsT);
matMeanMoveX = nan(intBins,intPlotOris,1);
matMeanMoveY = nan(intBins,intPlotOris,1);
matMeanRadius = nan(intBins,1);

%% go through files
for intFile=1:numel(cellFiles)
	strTarget = cellFiles{intFile};
	strRec = strTarget((end-10):(end-7));
	sLoad = load(fullfile(strDataSourcePath,strTarget));
	sAP = sLoad.sAP;
	%prepro
	sPupil = sAP.sPupil;
	vecPupilRadius = sPupil.vecPupilRadius;
	vecPupilTime = sPupil.vecPupilTime;
	vecPupilFixedRadius = sPupil.vecPupilFixedRadius;
	vecPupilFixedCenterX = sPupil.vecPupilFixedCenterX;
	vecPupilFixedCenterY = sPupil.vecPupilFixedCenterY;
	vecPupilFixedCenterY = max(vecPupilFixedCenterY) - vecPupilFixedCenterY; %invert y
	
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
	%vecMoveX = imfilt([0 diff(vecPupilFixedCenterX)],vecFilt);
	%vecMoveY = imfilt([0 diff(vecPupilFixedCenterY)],vecFilt);
	vecMoveX = vecFiltX;
	vecMoveY = vecFiltY;
	
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
		
		%plot
		intCounter = intCounter + 1;
		sOpt.handleFig = -1;
		[vecMean,vecSEM,vecWindowBinCenters,matPET] = doPEP(vecPupilTime,vecFiltR,vecPupilStimOnTime(vecGlitchesPerTrials<dblCutOff),sOpt);
		if numel(vecMean) < intBins,vecMean((end+1):intBins)=vecMean(end);end
		matMeanRadius(1:intBins,intCounter) = vecMean(1:intBins);
		
		dblOriStep = (pi/(intPlotOris/2));
		for intPlotOri=1:8
			%get which oris to include
			dblCenterOri = (intPlotOri-1)*dblOriStep;
			vecOriRad = deg2rad(vecOrientation);
			indUseOriTrials = abs(circ_dist(dblCenterOri,vecOriRad)) < (dblOriStep*0.49);
			indIncludeTrials = vecGlitchesPerTrials<dblCutOff;
			indUseT = indUseOriTrials & indIncludeTrials;
			%plot
			sOpt.handleFig = -1;
			[vecMean,vecSEM,vecWindowBinCenters,matPET] = doPEP(vecPupilTime,zscore(vecMoveX),vecPupilStimOnTime(indUseT),sOpt);
			if numel(vecMean) < intBins,vecMean((end+1):intBins)=vecMean(end);end
			matMeanMoveX(1:intBins,intPlotOri,intCounter) = vecMean(1:intBins);
			[vecMean,vecSEM,vecWindowBinCenters,matPET] = doPEP(vecPupilTime,zscore(vecMoveY),vecPupilStimOnTime(indUseT),sOpt);
			if numel(vecMean) < intBins,vecMean((end+1):intBins)=vecMean(end);end
			matMeanMoveY(1:intBins,intPlotOri,intCounter) = vecMean(1:intBins);
		end
	end
end

%% plot grand mean
matGrandMeanX = matMeanMoveX - mean(mean(matMeanMoveX,3),2);
matGrandMeanY = matMeanMoveY - mean(mean(matMeanMoveY,3),2);

%plot
figure;
hRadius = subplot(3,3,5);
errorfill(vecBinsT,mean(matMeanRadius,2),std(matMeanRadius,[],2)/sqrt(size(matMeanRadius,2)));
ylabel('Pupil radius');
xlabel('Time after stim onset (s)');
fixfig;

vecLimY = [-0.2 0.2];
strLabelY = 'Pupil location';
intPlotOris = 8;
vecPlotOrder = [6 3 2 1 4 7 8 9];
vecHandles = nan(size(vecPlotOrder));
dblOriStep = (pi/(intPlotOris/2));
for intPlotOri=1:8
	%get which oris to include
	
	%plot
	vecHandles(intPlotOri) = subplot(3,3,vecPlotOrder(intPlotOri));
	vecColor = [0 0 1];
	errorfill(vecBinsT,mean(matGrandMeanX(:,intPlotOri,:),3),std(matGrandMeanX(:,intPlotOri,:),[],3)/sqrt(size(matMeanRadius,2)),vecColor);
	vecColor = [1 0 0];
	errorfill(vecBinsT,mean(matGrandMeanY(:,intPlotOri,:),3),std(matGrandMeanY(:,intPlotOri,:),[],3)/sqrt(size(matMeanRadius,2)),vecColor);
	ylabel(strLabelY);
	xlabel('Time after stim onset (s)');
	title('b=horz,r=vert');
	ylim(vecLimY);
	fixfig;
end
maxfig;
