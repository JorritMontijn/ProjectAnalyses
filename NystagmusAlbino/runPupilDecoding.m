
%get data
strDataSourcePath = 'F:\Data\Processed\Neuropixels\';
sFiles = dir([strDataSourcePath '*.mat']);
cellFiles = {sFiles(:).name}';
close all;

%% pre-allocate
cellRunAreas = {...
	'lateral geniculate',...Area 1
	'Primary visual',...Area 
	'Lateral posterior nucleus',...Area 
	'Anterior pretectal nucleus',...Area 
	'Nucleus of the optic tract',...Area 
	'Superior colliculus',...Area 
	'Anteromedial visual',...Area 
	'posteromedial visual',...Area 
	'Anterolateral visual',...Area 
	'Lateral visual',...Area 
	'Rostrolateral area',...Area 
	'Anterior area',...Area 
	'Subiculum',...Area 
	'Field CA1',...Area 
	'Field CA2',...Area 
	'Field CA3',...Area 
	'Dentate gyrus',...Area 
	'Retrosplenial'...Area 
	};

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
	vecPupilRadius = sPupil.vecPupilRadius;
	vecPupilTime = sPupil.vecPupilTime;
	vecPupilFixedRadius = sPupil.vecPupilFixedRadius;
	vecPupilFixedCenterX = sPupil.vecPupilFixedCenterX;
	vecPupilFixedCenterY = sPupil.vecPupilFixedCenterY;
	
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
		%error
		%plot
		
		%go through areas
		cellClustAreas = {sAP.sCluster(:).Area};
		intArea = 8;%pm
		%for intArea = 8%1:numel(cellUniqueAreas)
		%% get data
		indUseNeurons = contains(cellClustAreas,cellRunAreas{intArea});
		indUseNeurons = true(size(cellClustAreas));
		cellSpikes = {sAP.sCluster(indUseNeurons).SpikeTimes};
		intNeurons = numel(cellSpikes);
		vecDecodeWindow = [0 1];
		vecUseTrialOnsets = vecPupilStimOnTime(vecGlitchesPerTrials<dblCutOff);
		intUseTrials = numel(vecUseTrialOnsets);
		vecPupilSize = nan(1,intUseTrials);
		vecPupilLocX = nan(1,intUseTrials);
		vecPupilLocY = nan(1,intUseTrials);
		matMeanAct = nan(intNeurons,intUseTrials);
		vecSpikeOnTime = sAP.cellStim{intStim}.structEP.vecStimOnTime;
		for intTrial=1:intUseTrials
			%pupil
			dblStartPupilT = vecPupilStimOnTime(intTrial) + vecDecodeWindow(1);
			dblStopPupilT = vecPupilStimOnTime(intTrial) + vecDecodeWindow(2);
			vecPupilSize(intTrial) = nanmean(vecFiltR(vecPupilTime > dblStartPupilT & vecPupilTime < dblStopPupilT));
			vecPupilLocX(intTrial) = nanmean(vecFiltX(vecPupilTime > dblStartPupilT & vecPupilTime < dblStopPupilT));
			vecPupilLocY(intTrial) = nanmean(vecFiltY(vecPupilTime > dblStartPupilT & vecPupilTime < dblStopPupilT));

			%spikes
			dblStartSpikeT = vecSpikeOnTime(intTrial) + vecDecodeWindow(1);
			dblStopSpikeT = vecSpikeOnTime(intTrial) + vecDecodeWindow(2);
			parfor intNeuron=1:intNeurons
				matMeanAct(intNeuron,intTrial) = sum(cellSpikes{intNeuron} > dblStartSpikeT & cellSpikes{intNeuron} < dblStopSpikeT);
			end
		end
		matMeanAct = matMeanAct./range(vecDecodeWindow);
		
		%% decode pupil size
		matMeanZ = zscore(matMeanAct,[],1)';
		vecPupilZ = zscore(vecPupilSize)';
		matY = vecPupilZ;%(randperm(numel(vecPupilZ)));
		varTypeCV = 1;
		dblLambda = 200;
		[matY_hat,dblR2_CV,matB] = doCrossValidatedDecodingRR(matMeanZ,matY,varTypeCV,dblLambda);
		
		% plot
		figure
		maxfig;
		dblLimAx = max(abs(cat(1,matY_hat,matY)));
		subplot(2,3,1);
		h2=plot(dblLimAx.*[-1 1],dblLimAx.*[-1 1],'--','Color',[0.5 0.5 0.5]);
		hold on
		h1=scatter(matY,matY_hat,[],lines(1),'.');
		hold off
		xlabel('Real pupil size (z-score)');
		ylabel('CV-decoded pupil size (z-score)');
		title(sprintf('%s; Trial-average pupil size; CV R^2=%.3f',strRec,dblR2_CV));
		fixfig;
		h1.SizeData=72;
		h2.LineWidth=1;
		
		subplot(2,3,4)
		[cellUniqueAreas,a,vecAreaIdx] = unique(cellClustAreas);
		cellGroupedAreas = cellfun(@(x,y) find(contains(x,y,'IgnoreCase',true)),cellfill(cellUniqueAreas,size(cellRunAreas)),cellRunAreas,'UniformOutput',false);
		intUniqueAreas = numel(cellUniqueAreas);
		vecGroupInto = nan(1,intUniqueAreas);
		vecRunAreaIdx = vecAreaIdx;
		for intArea=1:intUniqueAreas
			intNewIdx = find(cellfun(@ismember,cellfill(intArea,size(cellGroupedAreas)),cellGroupedAreas));
			if isempty(intNewIdx),intNewIdx=0;end
			vecGroupInto(intArea) = intNewIdx;
			vecRunAreaIdx(vecAreaIdx==intArea)=intNewIdx;
		end
		[u,a,b]=unique(vecRunAreaIdx)
		[varDataOut,vecUnique,vecCounts,cellSelect,vecRepetition] = label2idx(vecRunAreaIdx)
		%% decode horizontal pupil location
		vecHorzLocZ = zscore(vecPupilLocX)';
		varTypeCV = 1;
		%dblLambda = 1;
		[matY_hat,dblR2_CV,matB] = doCrossValidatedDecodingRR(matMeanZ,vecHorzLocZ,varTypeCV,dblLambda);
		
		%plot
		dblLimAx = max(abs(cat(1,matY_hat,vecHorzLocZ)));
		subplot(2,3,2);
		h2=plot(dblLimAx.*[-1 1],dblLimAx.*[-1 1],'--','Color',[0.5 0.5 0.5]);
		hold on
		h1=scatter(vecHorzLocZ,matY_hat,[],lines(1),'.');
		hold off
		xlabel('Real pupil x-location (z-score)');
		ylabel('CV-decoded pupil x-location (z-score)');
		title(sprintf('%s; Trial-average pupil x-loc; CV R^2=%.3f',strRec,dblR2_CV));
		fixfig;
		h1.SizeData=72;
		h2.LineWidth=1;
		
		%% decode vertical pupil location
		vecVertLocZ = zscore(vecPupilLocY)';
		varTypeCV = 1;
		%dblLambda = 1;
		matX = matMeanZ;
		matY = vecVertLocZ;%(randperm(numel(vecVertLocZ)));
		
		%%% test
		%dblLambda = 10;
		%matY = normrnd(0,2,[100 1]);
		%matX = bsxfun(@plus,matY,normrnd(0,1,[100 10]));
		[matY_hat,dblR2_CV,matB] = doCrossValidatedDecodingRR(matX,matY,varTypeCV,dblLambda);
		
		%plot
		dblLimAx = max(abs(cat(1,matY_hat,matY)));
		subplot(2,3,3);
		h2=plot(dblLimAx.*[-1 1],dblLimAx.*[-1 1],'--','Color',[0.5 0.5 0.5]);
		hold on
		h1=scatter(matY,matY_hat,[],lines(1),'.');
		hold off
		xlabel('Real pupil y-location (z-score)');
		ylabel('CV-decoded pupil y-location (z-score)');
		title(sprintf('%s; Trial-average pupil y-loc; CV R^2=%.3f',strRec,dblR2_CV));
		fixfig;
		h1.SizeData=72;
		h2.LineWidth=1;
		
		%get R^2
		vecMu = mean(matY,1);
	dblSSRes_ridge = sum(sum((matY - matY_hat).^2));
	dblSSTot = sum(sum(bsxfun(@minus,matY,vecMu).^2));
	dblR2_CV = 1 - dblSSRes_ridge / dblSSTot;


	end
end
