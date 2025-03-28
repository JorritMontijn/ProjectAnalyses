%% runPlaceDecoding Runs place decoding + GLM of calcium imaging data
%Uses a 6-factor model, including EyePosX, EyePosY, EyeSize, MotLoc, MotSpeed, Licking
%
%	Version History:
%	2019-04-08	Created by Jorrit Montijn

%% load data
clear all;
strPath = 'D:\Data\Processed\imagingKoen\';
strFile = 'Pogba_20190226_001_normcorr_SPSIG';
sLoadData = load([strPath strFile '.mat']);
sLoadMeta = load([strPath strFile '_Res.mat']);
intFoldK = 10;

%% get meta data
%get data time-stamps
vecFrameTimes =  sLoadMeta.info.Frametimes;
dblFrameDur = mean(diff(vecFrameTimes));
dblSamplingFreq = 1/dblFrameDur;

%get sync info
%info.frame [trials x 1] is sync frames for position=0
vecTrialStartFrames = sLoadMeta.info.frame;
vecTrialStartTime = vecFrameTimes(vecTrialStartFrames);
intTrialReps = numel(vecTrialStartFrames);

%get motion log [location-entry x 3]
%col1: running speed
%col2: time-stamp
%col3: position (0=start, 0.5=end)
matMotionLog = sLoadMeta.info.Stim.log.motionlog;
vecMotSpeed = matMotionLog(:,1);
vecMotTime = matMotionLog(:,2);
vecMotLoc = matMotionLog(:,3);

%re-align to logged trial starts
intPointNum = numel(vecMotLoc);
vecFramePos0 = [1;find(diff(vecMotLoc)<0)];
vecTimePos0 = vecMotTime(vecFramePos0);
vecOffset = vecTimePos0-vecTrialStartTime;
vecResOffset = get(resample(timeseries([vecOffset(1);vecOffset;vecOffset(end)],[1;vecFramePos0;intPointNum]),1:intPointNum),'Data');
vecResMotTime = vecMotTime-vecResOffset;
vecResTimePos0 = vecResMotTime(vecFramePos0);

%stimulus log
%stimlog = sLoadMeta.info.Stim.log.stimlog;
vecLickTimes = sLoadMeta.info.Stim.log.licklog(:,2);
vecLickFrames = nan(size(vecLickTimes));
for intLick=1:numel(vecLickFrames)
	vecLickFrames(intLick) = sum(vecFrameTimes<vecLickTimes(intLick));
end

%eye data
vecEyePosX = sLoadMeta.info.eye.pos(:,1);
vecEyePosY = sLoadMeta.info.eye.pos(:,2);
vecEyeArea = sLoadMeta.info.eye.area;
vecEyeTime = sLoadMeta.info.eye.time;

%% get neuronal data
%sigraw = F
%sig = F/F0
%decon=deconvolved F/F0
%den=fitted F/F0
%no pre-processed dF/F0
matAct = sLoadData.sig; %[frame x neuron]

%% transform predictors to same time bins
dblLickFilterSecs = 0.2;
dblFiltWidth = dblLickFilterSecs/dblFrameDur;
vecFiltX = ceil(2*dblFiltWidth)*[-1 1];
vecFilterLicking = normpdf(vecFiltX(1):vecFiltX(end),0,dblFiltWidth);
vecFilterLicking = vecFilterLicking/sum(vecFilterLicking(:));
cellPredictors = {'EyePosX','EyePosY','EyeSize','MotLoc','Licking','MotSpeed'};
intPredictors = numel(cellPredictors);
sResampled = struct;
vecResExtMotTime = [min([vecFrameTimes(1) vecResMotTime(1)]);vecResMotTime;max([vecFrameTimes(end) vecResMotTime(end)])];
sResampled.EyePosX = get(resample(timeseries(vecEyePosX,vecEyeTime),vecFrameTimes),'Data');
sResampled.EyePosY = get(resample(timeseries(vecEyePosY,vecEyeTime),vecFrameTimes),'Data');
sResampled.EyeSize = get(resample(timeseries(vecEyeArea,vecEyeTime),vecFrameTimes),'Data');
sResampled.MotLoc = get(resample(timeseries([0;vecMotLoc;vecMotLoc(end)],vecResExtMotTime),vecFrameTimes),'Data');
sResampled.Licking = false(size(vecFrameTimes));
sResampled.Licking(vecLickFrames) = true;
sResampled.Licking = conv(sResampled.Licking,vecFilterLicking,'same');
sResampled.MotSpeed = get(resample(timeseries([0;vecMotSpeed;vecMotSpeed(end)],vecResExtMotTime),vecFrameTimes),'Data');
%assign
matX = nan(size(matAct,1),intPredictors);
for intPredictor=1:intPredictors
	matX(:,intPredictor) = sResampled.(cellPredictors{intPredictor});
end
%remove nans and remove means
indGood = ~any(isnan(matX),2);
matX = matX(indGood,:);
matAct = matAct(indGood,:);
vecFrameTimes = vecFrameTimes(indGood);

%build trial repetition vector
vecTrialRepetition = zeros(size(vecFrameTimes));
for intRep=1:intTrialReps
	vecTrialRepetition(vecTrialStartFrames(intRep):end) = intRep;
end

%% build generative encoding model
intTimeBins = size(matAct,1);
intNeurons = size(matAct,2);
matR_squared = nan(intNeurons,intPredictors+1); %first predictor is all
matPredY = nan(intTimeBins,intNeurons,intPredictors+1); %first predictor is all
cellCoefficients = cell(intNeurons,intPredictors+1); %first predictor is bias
for intNeuron = 1:intNeurons
	%get data
	fprintf('Processing neuron %d/%d [%s]\n',intNeuron,intNeurons,getTime);
	vecY = matAct(:,intNeuron);
	
	%set starting parameters
	cellCoeffs0 = {0,0,0,[0,0,1],0,0};
	%CV?
	if intFoldK < 2
		%% non-cv
		%full model
		[cellTheseCoeffs,matOutX,cellFunctions] = gnmfit(matX,vecY,cellCoeffs0,'gnmgauss',4,'times_mult',6);
		cellCoefficients(intNeuron,:) = cellTheseCoeffs;
		%predict
		vecPredY = gnmval(cellTheseCoeffs,matOutX,cellFunctions);
		
	else
		%% cv
		%full model
		[vecPredY_CV,cellTheseCoeffs,matOutX,cellFunctions,cellFoldCoeffs,cellFoldTestObs] = gnmkfold(intFoldK,vecTrialRepetition,matX,vecY,cellCoeffs0,'gnmgauss',4,'times_mult',6);
		%predict
		cellCoefficients(intNeuron,:) = cellTheseCoeffs;
		vecPredY = vecPredY_CV;
	end
	%save data
	matPredY(:,intNeuron,1) = vecPredY;
	
	%get R^2
	warning('off','getR2:NaNsDetected');
	[dblR2,dblSS_tot,dblSS_res] = getR2(vecY,vecPredY);
	matR_squared(intNeuron,1) = dblR2;
	warning('on','getR2:NaNsDetected');
	
	%% remove each component to calculate its contribution
	for intPredictor=1:intPredictors
		%remove feature
		matTempX = matOutX;
		matTempX(:,intPredictor) = [];
		%prepare arguments
		cellTheseCoeffs0 = cellTheseCoeffs;
		cellTheseCoeffs0(intPredictor) = [];
		if intPredictor == 4
			cellInputs = {'times_mult',5,'constant',intPredictors};
		elseif intPredictor < 4
			cellInputs = {'gnmgauss',3,'times_mult',5,'constant',intPredictors};
		elseif intPredictor > 4
			if intPredictor == 5
				cellInputs = {'gnmgauss',4,'times_mult',5,'constant',intPredictors};
			elseif intPredictor == 6
				cellInputs = {'gnmgauss',4,'constant',intPredictors};
			end
		end
		
		%recalculate
		if intFoldK < 2
			%get model
			[cellTheseCoeffsRem,matTempOutX,cellTempFunctions] = gnmfit(matTempX,vecY,cellTheseCoeffs0,cellInputs{:});
			%predict
			vecTempPredY = gnmval(cellTheseCoeffsRem,matTempOutX,cellTempFunctions);
			
		else
			%get k-fold CV
			[vecTempPredY,cellMeanCoeffs,matX,cellFunctions,cellFoldCoeffs,cellFoldTestObs] = gnmkfold(intFoldK,vecTrialRepetition,matX,vecY,cellTheseCoeffs0,cellInputs{:});
		end
		
		%save prediction
		matPredY(:,intNeuron,intPredictor+1) = vecTempPredY;
		%get R^2
		warning('off','getR2:NaNsDetected');
		[dblR2,dblSS_tot,dblSS_res] = getR2(vecY,vecTempPredY);
		warning('on','getR2:NaNsDetected');
		matR_squared(intNeuron,intPredictor+1) = dblR2;
	end
end
%% calculate unique contributions
matR2_uniqueContrib = bsxfun(@minus,matR_squared(:,1),matR_squared(:,2:end));
cellPredictorsRemoved = [{''} cellPredictors];

%% plot gnm
if numel(strFile) > 20
	strTitFile = [strFile(1:17) '...'];
else
	strTitFile = strFile;
end
%%
figure
for intPredictor=1:intPredictors
	subplot(2,3,intPredictor);
	hold all;
	histx(matR_squared(:,intPredictor+1));
	ylabel('Number of neurons (count)');
	xlabel('R^2');
	xlim([0 0.6])
	title([cellPredictorsRemoved{intPredictor+1} ' removed'])
	fixfig;
	set(get(gca,'Children'),'Linewidth',1);
end
%%
figure
for intPredictor=1:intPredictors
	subplot(2,3,intPredictor);
	hold all;
	histx(matR2_uniqueContrib(:,intPredictor));
	ylabel('Number of neurons (count)');
	xlabel('R^2');
	xlim([0 0.6])
	title(['Unique contrib of ' cellPredictorsRemoved{intPredictor+1}])
	fixfig;
	set(get(gca,'Children'),'Linewidth',1);
end
%%
vecPlot = 1:5000;
figure
subplot(2,2,1)
dblBinStep = 0.1;
vecBins = 0:dblBinStep:0.6;
vecBinPlotX = vecBins(2:end)-dblBinStep/2;
vecCounts = histcounts(matR_squared(:,1),vecBins);
bar(vecBinPlotX,vecCounts);
ylabel('Number of neurons (count)');
xlabel('R^2');
xlim([0 0.6])
title(sprintf('K-fold CV=%d, Full model, all cells, %s',intFoldK,strTitFile),'Interpreter','none');
fixfig;
set(get(gca,'Children'),'Linewidth',1);

[dblBestR2,intBestNeuron]=max(matR_squared(:,1));
vecPredictedBestNeuron = matPredY(:,intBestNeuron,1);
vecRealBestNeuron = matAct(:,intBestNeuron);

subplot(2,2,2)
scatter(vecRealBestNeuron(vecPlot),vecPredictedBestNeuron(vecPlot),'xk');
dblR = nancorr(vecRealBestNeuron(vecPlot),vecPredictedBestNeuron(vecPlot));
ylabel('Predicted activity (F/F0)');
xlabel('Real activity (F/F0)');
title(sprintf('First %dk frames, best neuron %d; r=%.3f',numel(vecPlot)/1000,intBestNeuron,dblR));
fixfig;


subplot(2,2,[3 4])
plot(vecFrameTimes(vecPlot),vecRealBestNeuron(vecPlot),'k');
hold on
plot(vecFrameTimes(vecPlot),vecPredictedBestNeuron(vecPlot),'r');
hold off
ylabel('Activity (F/F0)');
xlabel('Time (s)');
title(sprintf('K-fold CV=%d, prediction best neuron %d (overall R^2=%.3f), first %dk frames; black=real, red=predicted',intFoldK,intBestNeuron,dblBestR2,numel(vecPlot)/1000));
fixfig;
%set(get(gca,'Children'),'Linewidth',1);

%% subtract activity from other factors
matActClean = zeros(size(matAct));
for intNeuron=1:intNeurons
	cellTheseCoeffs = cellCoefficients(intNeuron,:);
	vecY = matAct(:,intNeuron);
	cellTheseCoeffs{4}(1)=0;
	matActClean(:,intNeuron) = vecY - gnmval(cellTheseCoeffs,matOutX,cellFunctions);
end

%% decode
%bin location data
dblBinStep = 0.01;
vecBins = (-dblBinStep/2):dblBinStep:1;
vecX = matOutX(:,4);
vecResps = vecY;
[vecCounts,vecMeans,vecSDs,cellVals,cellIDs] = makeBins(vecX,vecX,vecBins);
intFirstRemBin = find(vecBins>0.5);%find(vecCounts < 100,1);
cellRems = cellIDs(intFirstRemBin:end);
cellBins = cellIDs(1:(intFirstRemBin-1));
intLocBins = numel(cellBins)+1;
vecLocBin = nan(size(vecX));
for intLocBin=1:(intLocBins-1)
	vecLocBin(cellBins{intLocBin}) = intLocBin;
end
vecLocBin(cell2vec(cellRems)) = intLocBins;

%make figure
figure;
for intDecodeType=1:2
	if intDecodeType==1
		matData=matAct;
		strTitle = strTitFile;
	else
		matData=matActClean;
		strTitle = 'gnm-subtracted';
	end
	%decode
	[dblPerformance,vecDecodedIndexCV,matMahalDistsCV,dblMeanErrorDegs,matConfusion] = ...
		doCrossValidatedDecodingMD(matData,vecLocBin,vecTrialRepetition);
	
	% plot
	subplot(2,2,1+(intDecodeType-1)*2);
	imagesc(matConfusion);colormap(hot)
	ylabel('Decoded location');
	xlabel('Real location');
	title(sprintf('Raw pdf, %s',strTitle),'interpreter','none');
	colorbar;
	fixfig;
	grid off;
	
	%normalize
	vecObservedCounts = accumarray(vecLocBin,ones(size(vecLocBin)));
	matConfusionNorm = bsxfun(@rdivide,matConfusion,vecObservedCounts');
	vecConfCount = sum(matConfusionNorm,1); %should all be one
	matCN2 = matConfusionNorm;
	matCN2(end,end) = 0;
	dblMax = max(matCN2(:));
	
	%calculate MSE
	intBins = size(matConfusion,1);
	matErrDist = (abs((1:intBins) - 1 * (1:intBins)'));
	matError = matErrDist .* matConfusionNorm;
	dblMSE = mean(matError(:).^2);
	
	subplot(2,2,2+(intDecodeType-1)*2);
	imagesc(matConfusionNorm,[0 dblMax]);colormap(hot)
	ylabel('Decoded location');
	xlabel('Real location');
	title(sprintf('MSE=%.3f, Norm pdf, %s',dblMSE,strTitle),'interpreter','none');
	colorbar;
	fixfig;
	grid off;
end