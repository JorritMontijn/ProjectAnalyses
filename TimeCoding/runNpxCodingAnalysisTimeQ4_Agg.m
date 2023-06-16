%% aim
%{
is population gain correlated with pupil size?
%}

%% define qualifying areas
clear all;
cellUseAreas = {...
	'Primary visual area',...
	...'posteromedial visual area',...
	};
if isfolder('F:\Drive\PopTimeCoding') && isfolder('F:\Data\Processed\Neuropixels\')
	strDataPath = 'F:\Data\Processed\Neuropixels\';
	strFigurePathSR = 'F:\Drive\PopTimeCoding\single_recs';
	strFigurePath = 'F:\Drive\PopTimeCoding\figures\';
	strTargetDataPath = 'F:\Drive\PopTimeCoding\data\';
else
	strDataPath = 'E:\DataPreProcessed\';
	strFigurePathSR = 'C:\Drive\PopTimeCoding\single_recs';
	strFigurePath = 'C:\Drive\PopTimeCoding\figures\';
	strTargetDataPath = 'C:\Drive\PopTimeCoding\data\';
end


%% select all neurons in LP and drifting grating stimuli
if ~exist('sAggStim','var') || isempty(sAggStim) || isempty(sAggNeuron)
	[sAggStim,sAggNeuron]=loadDataNpx('','natural',strDataPath);
end
indRemDBA = strcmpi({sAggNeuron.SubjectType},'DBA');
fprintf('Removing %d cells of DBA animals; %d remaining [%s]\n',sum(indRemDBA),sum(~indRemDBA),getTime);
sAggNeuron(indRemDBA) = [];

%% pre-allocate matrices
intAreas = numel(cellUseAreas);

%% go through recordings
tic
for intRec=1:numel(sAggStim) %19 || weird: 11
	% get matching recording data
	strRec = sAggStim(intRec).Exp;
	sThisRec = sAggStim(strcmpi(strRec,{sAggStim(:).Exp}));
	
	%prep nm data
	[sUseNeuron,vecStimOnTime,vecStimOffTime,vecStimIdx] = NpxPrepMovies(sAggNeuron,sThisRec,cellUseAreas);
	[vecFrameIdx,vecUniqueFrames,vecRepNum,cellSelect,vecTrialRepetition] = val2idx(vecStimIdx);
	intTrialNum = numel(vecStimOnTime);
	intStimNum = numel(unique(vecStimIdx));
	intRepNum = intTrialNum/intStimNum;
	
	%original movie time
	vecOrigStimOnTime = sThisRec.cellBlock{1}.vecStimOnTime;
	vecOrigStimOffTime = sThisRec.cellBlock{1}.vecStimOffTime;
	dblStimDur = min(diff(vecOrigStimOnTime));
	intOrigTrialNum = numel(vecOrigStimOnTime);
	
	%% select area 1
	for intArea=1:intAreas
		strArea = cellUseAreas{intArea};
		indArea1Neurons = contains({sUseNeuron.Area},strArea,'IgnoreCase',true);
		if sum(indArea1Neurons) == 0, continue;end
		
		%% get orientation responses & single-trial population noise
		sArea1Neurons = sUseNeuron(indArea1Neurons);
		
		%% prep data
		%get data matrix
		cellSpikeTimesRaw = {sArea1Neurons.SpikeTimes};
		[matMeanRate,cellSpikeTimes] = ...
			NpxPrepMovieData(cellSpikeTimesRaw,vecStimOnTime,vecStimOffTime,vecFrameIdx);
		
		%events
		dblStartEpoch = vecOrigStimOnTime(1)-10;
		dblStopEpoch = vecOrigStimOnTime(end)+dblStimDur+10;
		dblEpochDur = dblStopEpoch - dblStartEpoch;
		intNumN = numel(cellSpikeTimes);
		
		%% calc gain
		%calc gain
		dblWindow = 1; %secs
		vecGainE = dblStartEpoch:dblWindow:dblStopEpoch;
		vecGainT = vecGainE(2:end) - dblWindow/2;
		vecStimE = 0:dblWindow:8;
		vecStimC = vecStimE(2:end) - dblWindow/2;
		[~,vecSEM,vecWindowBinCenters,matPET] = doPEP(cellSpikeTimes,vecStimE,vecOrigStimOnTime);
		matPET = shiftdim(matPET,1);
		matStimR = mean(matPET,3);
		matPET_Noise = bsxfun(@minus,matPET,matStimR);
		intBinNum = numel(vecStimC)*size(matPET,3);
		matBinnedRate = nan(intNumN,intBinNum);
		matBinnedNoise = nan(intNumN,intBinNum);
		for intN=1:intNumN
			matBinnedRate(intN,:)= flat(matPET(:,intN,:));
			matBinnedNoise(intN,:)= flat(matPET_Noise(:,intN,:));
		end
		matBinnedStim = matBinnedRate-matBinnedNoise;
		vecGainAx = mean(matBinnedRate,2);
		
		vecGain=getProjOnLine(matBinnedNoise,vecGainAx)';
		vecNormGain = vecGain./norm(vecGainAx);
		vecMean=mean(matBinnedRate,1);
		vecNormMean = vecMean./mean(vecMean(:));
		
		%% are fluctuations in the direction of the mean-axis or gain-axis?
		%calculate residuals after subtracting gain fluctuations or mean fluctuations
		%matBinnedRate = matBinnedRate./std(matBinnedRate,[],2);
		
		%not cv
		%gain
		matPredGain = vecGainAx*vecNormGain;
		matResidGain = matBinnedNoise-matPredGain;
		intK_G = intNumN;
		intK_M = 1;
		
		%mean
		matPredMean = vecGainAx*vecNormMean;
		matResidMean = matBinnedNoise-matPredMean;
		
		%prep
		matPredGainCV = nan(intNumN,intBinNum);
		matPredMeanCV = nan(intNumN,intBinNum);
		matPredRegCV = nan(intNumN,intBinNum);
		%cross-validate by leaving one neuron out at a given time-point
		for intBin=1:intBinNum
			%get gain axis from all other bins
			indOtherBins = true(1,intBinNum);
			indOtherBins(intBin) = false;
			vecGainAxCv = mean(matBinnedRate(:,indOtherBins),2);
			
			for intN=1:intNumN
				%select other neurons
				indOtherNeurons = true(1,intNumN);
				indOtherNeurons(intN) = false;
				%get cv gain axis for other neurons
				vecOtherGainAxCv = vecGainAxCv(indOtherNeurons);
				%get response of other neurons in this trial
				vecTestOtherR = matBinnedNoise(indOtherNeurons,intBin);
				
				%get gain of rest of population in this bin
				dblOtherGain=getProjOnLine(vecTestOtherR,vecOtherGainAxCv)/norm(vecOtherGainAxCv);
				dblOtherMean=mean(vecTestOtherR);
				
				%predict this repsonse
				matPredGainCV(intN,intBin) = matBinnedStim(intN,intBin)+vecGainAxCv(intN)*dblOtherGain;
				matPredMeanCV(intN,intBin) = matBinnedStim(intN,intBin)+dblOtherMean;
				
				%% regression
				%select other neurons
				matX_Train = matBinnedNoise(indOtherNeurons,indOtherBins)';
				matY_Train = matBinnedNoise(intN,indOtherBins)';
				matX_Test = matBinnedNoise(indOtherNeurons,intBin);
				matY_Test = matBinnedNoise(intN,intBin);
			
				[matW, dblMSE, intRankT, sSuppOut] = doRdRankReg(matX_Train, matY_Train, 1);
				dblOtherR=getProjOnLine(matX_Test,matW)/norm(matX_Test);
				matPredRegCV(intN,intBin) = matBinnedStim(intN,intBin)+vecGainAxCv(intN)*dblOtherR;
				
			end
		end
		matResidMeanCV = matBinnedRate-matPredMeanCV;
		matResidGainCV = matBinnedRate-matPredGainCV;
		matResidRegCV = matBinnedRate-matPredRegCV;
		
			
		%fComp = @(x) 1-getR2(matBinnedRate,((x*matPredMeanCV+(1-x)*matPredGainCV)));
		%dblMeanWeightFrac = fminsearch(fComp,0.5);
		%matPredCompCV = dblMeanWeightFrac*matPredMeanCV+(1-dblMeanWeightFrac)*matPredGainCV;
		%matResidCompCV = matBinnedRate-matPredCompCV;
		%[dblR2_C,dblSS_tot_C,dblSS_res_C,dblT_C,dblP_C,dblR2_adjusted_C,dblR2_SE_C] = getR2(matBinnedRate,matPredCompCV,intK_M);
		%dblR2_C2 = 1-fComp(dblMeanWeightFrac);
		
		%get predictability
		[dblR2_G,dblSS_tot_G,dblSS_res_G,dblT_G,dblP_G,dblR2_adjusted_G,dblR2_SE_G] = getR2(matBinnedRate-matBinnedStim,matPredGainCV-matBinnedStim,intK_G);
		[dblR2_M,dblSS_tot_M,dblSS_res_M,dblT_M,dblP_M,dblR2_adjusted_M,dblR2_SE_M] = getR2(matBinnedRate-matBinnedStim,matPredMeanCV-matBinnedStim,intK_M);
		[dblR2_R,dblSS_tot_R,dblSS_res_R,dblT_R,dblP_R,dblR2_adjusted_R,dblR2_SE_R] = getR2(matBinnedRate-matBinnedStim,matPredRegCV-matBinnedStim,intK_G*intK_G);
		
		%plot
		close all;
		figure;maxfig;
		h1=subplot(2,3,1);
		imagesc(h1,matBinnedNoise)
		colorbar(h1);
		colormap(h1,redblue);
		set(h1,'clim',max(abs(get(gca,'clim')))*[-1 1]);
		ylabel('Neuron #');
		xlabel('Time (s)');
		title('Real fluctuations (Hz)');
		
		h3=subplot(2,3,2);
		imagesc(h3,matResidMeanCV);
		colorbar(h3);
		title(h3,sprintf('Residuals isotropic model (Mean), R2=%.3f',dblR2_adjusted_M));
		colormap(h3,redblue);
		set(h3,'clim',max(abs(get(gca,'clim')))*[-1 1]);
		ylabel('Neuron #');
		xlabel('Time (s)');
		
		h6=subplot(2,3,3);
		imagesc(h6,matResidGainCV);
		colorbar(h6);
		title(h6,sprintf('Residuals anisotropic model (Gain), R2=%.3f',dblR2_adjusted_G));
		set(h6,'clim',max(abs(get(gca,'clim')))*[-1 1]);
		colormap(h6,redblue);
		ylabel('Neuron #');
		xlabel('Time (s)');
		
		h6=subplot(2,3,4);
		imagesc(h6,matResidRegCV);
		colorbar(h6);
		title(h6,sprintf('Residuals optimal linear prediction (Reg), R2=%.3f',dblR2_adjusted_R));
		set(h6,'clim',max(abs(get(gca,'clim')))*[-1 1]);
		colormap(h6,redblue);
		ylabel('Neuron #');
		xlabel('Time (s)');
		
		
		
		h4=subplot(2,3,6);
		errorbar([1 2 3],[dblR2_adjusted_M dblR2_adjusted_G dblR2_adjusted_R],[dblR2_SE_M dblR2_SE_G dblR2_SE_R],'x');
		xlim([0.5 3.5]);
		set(gca,'xtick',[1 2 3],'xticklabel',{'Mean','Gain','Reg'});
		ylabel('Fluctuation prediction, adjusted R^2');
		title(sprintf('Mean=%.3f; Gain=%.3f',dblR2_adjusted_M,dblR2_adjusted_G))
		
		fixfig;
		
		%save figure
		export_fig(fullpath(strFigurePathSR,sprintf('Q4_NoisePredictionGainMean_%s.tif',strRec)));
		export_fig(fullpath(strFigurePathSR,sprintf('Q4_NoisePredictionGainMean_%s.pdf',strRec)));
		
		%% is population gain correlated with pupil size?
		if isfield(sThisRec,'Pupil')
			
			%get pupil data
			sPupil = sThisRec.Pupil;
			x0 = sThisRec.cellBlock{1}.vecPupilStimOn(1)-vecOrigStimOnTime(1);
			dblCorrection = getPupilOffset(cellSpikeTimes,x0,sPupil);
			vecPupilStimOn = vecOrigStimOnTime-dblCorrection;
			
			%bin
			vecTime = sPupil.vecTime;
			vecVals = sPupil.vecRadius;
			[matTE,vecWindowBinCenters] = getRespMat(vecTime,vecVals,vecPupilStimOn,vecStimE);
			vecPupilSize = flat(matTE');
			indUseTrials = ~isnan(vecPupilSize(:)) & ~isnan(vecNormMean(:)) & ~isnan(vecNormGain(:));
			vecG = vecNormGain(:);
			vecM = vecNormMean(:);
			vecP = vecPupilSize(:);
			[R_G,P_G,RL_G,RU_G] = corrcoef(vecG(indUseTrials),vecP(indUseTrials));
			[R_M,P_M,RL_M,RU_M] = corrcoef(vecM(indUseTrials),vecP(indUseTrials));
			dblCorrGainPupil = R_G(1,2);
			dblCorrMeanPupil = R_M(1,2);
			
			%run predictions
			[vecTrialGain_PupilHat,dblR2_CV_G,matB_G] = doCrossValidatedDecodingRR(zscore(vecPupilSize(indUseTrials)),zscore(vecNormGain(indUseTrials))',1,1);
			[vecTrialMean_PupilHat,dblR2_CV_M,matB_M] = doCrossValidatedDecodingRR(zscore(vecPupilSize(indUseTrials)),zscore(vecNormMean(indUseTrials))',1,1);
			
			[dblR2G,dblSS_totG,dblSS_resG,dblTG,dblPG,dblR2_adjustedG,dblR2_SEG] = getR2(zscore(vecNormGain(indUseTrials))',vecTrialGain_PupilHat,1);
			[dblR2M,dblSS_totM,dblSS_resM,dblTM,dblPM,dblR2_adjustedM,dblR2_SEM] = getR2(zscore(vecNormMean(indUseTrials))',vecTrialMean_PupilHat,1);
			
			
			figure;maxfig;
			subplot(2,3,1)
			scatter(zscore(vecPupilSize),zscore(vecNormGain))
			xlabel('Z-scored pupil size');
			ylabel('Z-scored pop gain');
			
			subplot(2,3,2)
			errorbar([1 2],[R_M(1,2) R_G(1,2)],[RL_M(1,2) RL_G(1,2)],[RU_M(1,2) RU_G(1,2)],'x')
			set(gca,'xtick',[1 2],'xticklabel',{'Mean-Pupil','Gain-Pupil'});
			ylabel('Pearson correlation (r)');
			title(sprintf('Corr, mean=%.3f, gain=%.3f',dblCorrMeanPupil,dblCorrGainPupil))
			xlim([0.5 2.5]);
			
			fixfig;
		end
		%save figure
		export_fig(fullpath(strFigurePathSR,sprintf('Q4_PupilPredGainMean_%s.tif',strRec)));
		export_fig(fullpath(strFigurePathSR,sprintf('Q4_PupilPredGainMean_%s.pdf',strRec)));
		
		%% save
		%save data
		save(fullpath(strTargetDataPath,sprintf('Q4Data_%s',strRec)),...
			'strRec',...
			'dblP_M',...
			'dblR2_adjusted_M',...
			'dblR2_adjusted_G',...
			'dblR2_adjusted_R',...
			'dblR2_M',...
			'dblR2_G',...
			'dblR2_R',...
			'dblCorrMeanPupil',...
			'dblCorrGainPupil'...
			);
	end
end
toc
