%% aim

%% set parameters
cellDataTypes = {'Npx','Sim','ABI','SWN'};%topo, model, allen, nora
intRunDataType = 1;
strRunStim = 'DG';%DG or NM? => superseded to WS by SWN
dblRemOnset = 0; %remove onset period in seconds; 0.125 for sim, 0.25 for npx
runHeaderPopTimeCoding;

%% specific parameters
cellRandomize = {'Real','TShuff','TPoiss','xPerm'};% {'Real','TShuff','TPoiss','TUniStretch','TSdFixed','TSaturating','TSdScaling','TSdLinear','TSdQuad','TSdCude'};
boolSaveData = true;
boolMakeFigs = true;
boolSaveFigs = true;

%% go through recordings
tic
for intRec=1:intRecNum
	close all;
	%% prep ABI or Npx data
	if strcmp(strRunType,'ABI')
		error to be updated
		runRecPrepABI;
		strThisRec = strRec;
	elseif strcmp(strRunType,'Sim')
		%load
		runRecPrepSim;
		
		%edit vars
		strThisRec = strRec;
		strDataPathT0=strDataPathSimT0;
		vecOri180 = mod(vecOrientation,180)*2;
		vecStimIdx = vecOri180;
		
		%% move onset
		%remove first x ms
		vecOrigStimOnTime = vecStimOnTime + dblRemOnset;
		
		%get data matrix
		[matData,indTuned,cellSpikeTimes,sOut,cellSpikeTimesPerCellPerTrial,indResp] = ...
			SimPrepData(cellSpikeTimesRaw,vecOrigStimOnTime,vecStimOffTime,vecOrientation);
		
		%get cell props
		%vecNeuronPrefOri = pi-[sLoad.sNeuron.PrefOri];
		vecNeuronType = [sLoad.sNeuron.Types]; %1=pyr,2=interneuron
		
	elseif strcmp(strRunType,'Npx') || strcmp(strRunType,'SWN')
		%prep
		runRecPrepNpx;
		strThisRec = strRec;
		strDataPathT0 = strTargetDataPath;
		
		%% move onset
		%remove first x ms
		vecOrigStimOnTime = vecStimOnTime + dblRemOnset;
		
		%get data matrix
		[matData,indTuned,cellSpikeTimes,sTuning24,cellSpikeTimesPerCellPerTrial] = ...
			NpxPrepData(cellSpikeTimesRaw,vecOrigStimOnTime,vecStimOffTime,vecOrientation);
		
		%get cell props
		vecNeuronType = ones(size(indTuned)); %1=pyr,2=interneuron
		%narrow vs broad not done
	else
		error impossible
	end
	
	%get ori vars
	cellOrigSpikeTimes = cellSpikeTimes;
	intTunedN = sum(indTuned);
	intRespN = size(matData,1);
	intNumN = numel(cellSpikeTimes);
	intTrialNum = numel(vecOrigStimOnTime);
	vecOri180 = mod(vecOrientation,180);
	[vecStimIdx,vecUnique,vecPriorDistribution,cellSelect,vecRepetition] = val2idx(vecOri180);
	sTuning = getTuningCurves(matData,vecOri180,0);
	vecNeuronPrefOri = sTuning.matFittedParams(:,1);
	vecNeuronBandwidth = real(sTuning.matBandwidth);
	
	intOriNum = numel(vecUnique);
	intRepNum = min(vecPriorDistribution);
	dblStimDur = median(vecStimOffTime - vecOrigStimOnTime);
	if mean(sum(matData)) < 90%90 / 50
		fprintf('Avg # of spikes per trial was %.1f for %s; skipping...\n',mean(sum(matData)),strThisRec);
		continue;
	else
		fprintf('Running %s... [%s]\n',strThisRec,getTime);
	end
	
	%% take only period during stimuli
	cellSpikeTimes = cell(size(cellOrigSpikeTimes));
	cellSpikeTimesPerCellPerTrial = cell(numel(cellOrigSpikeTimes),intTrialNum);
	for intN=1:numel(cellOrigSpikeTimes)
		[vecPseudoSpikeTimes,vecPseudoEventT] = getPseudoSpikeVectors(cellOrigSpikeTimes{intN},vecOrigStimOnTime,dblStimDur,true);
		
		dblStart = vecPseudoEventT(1);
		dblEnd = vecPseudoEventT(end)+dblStimDur;
		dblTotDur = dblEnd-dblStart;
		vecPseudoSpikeTimes(vecPseudoSpikeTimes<dblStart | vecPseudoSpikeTimes>dblEnd) = [];
		vecPseudoSpikeTimes = vecPseudoSpikeTimes-dblStart;
		dblLastSpike = max(vecPseudoSpikeTimes);
		
		cellSpikeTimes{intN} = vecPseudoSpikeTimes;
		
		%real
		[vecTrialPerSpike,vecTimePerSpike] = getSpikesInTrial(vecPseudoSpikeTimes,vecPseudoEventT,dblStimDur);
		for intTrial=1:intTrialNum
			vecSpikeT = vecTimePerSpike(vecTrialPerSpike==intTrial);
			cellSpikeTimesPerCellPerTrial{intN,intTrial} = vecSpikeT;
		end
	end
	
	%pupil
	if isfield(structStim,'vecPupilStimOn')
		vecPupilStimOn = structStim.vecPupilStimOn;
		vecPupilStimOff = structStim.vecPupilStimOff;
	else
		vecPupilStimOn =  [];
		vecPupilStimOff = [];
	end
	
	%running
	if isfield(structStim,'vecRunningTime')
		vecRunningTime = structStim.vecRunningTime;
		vecRunningSpeed = structStim.vecRunningSpeed;
	else
		vecRunningTime =  [];
		vecRunningSpeed = [];
	end
	
	
	%% run analysis
	for intRandomize=1:numel(cellRandomize)
		strType = cellRandomize{intRandomize};
		strTypeRec = [strType '_' strRec];
	
		%% define quantiles and remove zero-variance neurons
		%remove zero-variance neurons
		matSpikeCounts_pre = cellfun(@(x) sum(x>dblUseStartT),cellSpikeTimesPerCellPerTrial);
		matMeanRate_pre = matSpikeCounts_pre./dblUseMaxDur;
		vecPopRate_pre = sum(matMeanRate_pre,1);
		
		intQuantiles = 5;
		vecStartTrials = round(intRepNum*linspace(1/intRepNum,(1+1/intRepNum),intQuantiles+1));
		vecStartTrials(end)=[];
		intSplitTrialsPerOri = min(floor(cellfun(@sum,cellSelect)/intQuantiles));
		indZeroVarNeurons = false(intNumN,1);
		for intQ=1:intQuantiles
			vecUseTrialsTemp = vecStartTrials(intQ):(vecStartTrials(intQ)+intSplitTrialsPerOri-1);
			for intOri=1:intStimNum
				vecThisOri = find(cellSelect{intOri});
				[vecSorted,vecReorder]=sort(vecPopRate_pre(vecThisOri));
				vecQualifyingTrials = vecThisOri(vecReorder(vecUseTrialsTemp));
				indZeroVarNeurons = indZeroVarNeurons | (var(matMeanRate_pre(:,vecQualifyingTrials),[],2) == 0);
			end
		end
		indZeroVarNeurons = false(size(indZeroVarNeurons));
		vecUseNeurons = find(~indZeroVarNeurons);
		vecRemNeurons = find(indZeroVarNeurons);
		cellSpikeTimesPerCellPerTrial(vecRemNeurons,:) = [];
		intNeuronNum = numel(vecUseNeurons);
		if intNeuronNum < 25,continue;end
		
		%simple "rate code"
		matSpikeCounts = cellfun(@(x) sum(x>dblUseStartT),cellSpikeTimesPerCellPerTrial);
		matMeanRate = matSpikeCounts./dblUseMaxDur;
		
		%randomize per orientation
		if strcmp(strType,'TShuff') || strcmp(strType,'TPoiss')
			for intStim=1:intStimNum
				vecUseT = find(vecStimIdx==intStim);
				for intN=1:intRespN
					if strcmp(strType,'TShuff')
						%shuffle spikes
						matMeanRate(intN,vecUseT) = matMeanRate(intN,vecUseT(randperm(numel(vecUseT))));
					elseif strcmp(strType,'TPoiss')
						%generate spikes
						dblMean = mean(cellfun(@(x) sum(x>dblUseStartT),cellSpikeTimesPerCellPerTrial(intN,vecUseT)));
						vecRates = poissrnd(dblMean,size(vecUseT));
						matMeanRate(intN,vecUseT) = vecRates;
					else
						error
					end
				end
			end
		elseif strcmp(strType,'xPerm')
			for intN=1:intNumN
				matMeanRate(intN,:) = matMeanRate(intN,randperm(intTrialNum));
			end
		end
		
		%run approximation by mean-rate multiplied with tuning curve
		sOut = getTuningCurves(matMeanRate,vecOri180,0);
		matTuningCurves = sOut.matMeanResp;
		vecMeanPopRate = mean(matMeanRate,1);
		vecNormPopRate = vecMeanPopRate ./ mean(vecMeanPopRate);
		
		%% get CV model prediction
		fprintf('   Running model predictions [%s]\n',getTime);
		%leave one repetition of one neuron out
		matTunePredCV = nan(size(matMeanRate));
		matMeanPredCV = nan(size(matMeanRate));
		matGainPredCV = nan(size(matMeanRate));
		matGain1PredCV = nan(size(matMeanRate));
		for intTestRep = 1:intRepNum
			indTestTrials = vecTrialRepetition==intTestRep;
			vecTestTrials = find(indTestTrials);
			vecTrainTrials = find(~indTestTrials);
			matTrainR = matMeanRate(:,vecTrainTrials);
			vecTrainS = vecStimIdx(vecTrainTrials);
			vecTestS = vecStimIdx(vecTestTrials);
			matTestR = matMeanRate(:,vecTestTrials);
			dblMeanPopRate = mean(matTrainR(:));
			%get tuning curves
			[matRespNSR,vecStimTypes,vecUniqueDegs] = getStimulusResponses(matTrainR,vecTrainS);
			matTrainTuning = mean(matRespNSR,3);
			%get gain for training trials
			for intTestNeuron=1:intNeuronNum
				%% get vars
				indOtherNeurons = ~ismember(1:intNeuronNum,intTestNeuron);
				matOtherR = matMeanRate(indOtherNeurons,vecTestTrials);
				vecTrainOtherRate = mean(matTrainR(indOtherNeurons,:),1)';
				vecTestOtherRate = mean(matOtherR,1)';
				matTrainOtherR = matTrainR(indOtherNeurons,:);
				matTestOtherR = matTestR(indOtherNeurons,:);
				
				%% simple prediction: only tuning curve
				matTunePredCV(intTestNeuron,vecTestTrials) = matTrainTuning(intTestNeuron,vecTestS);
				
				%% classic prediction: tuning curve times population mean
				%get pop mean per trial, and average pop mean over trials
				vecPopMeanTrain = mean(matTrainOtherR,1);
				dblAvgPopMeanTrain = mean(vecPopMeanTrain);
				
				%get test pop mean
				vecPopMeanTest =  mean(matTestOtherR,1);
				
				%get mean resp during train
				vecTrainR = matTrainTuning(intTestNeuron,vecTestS);
				
				%predict: R(phi,i,t)= (pop_mean(not i,t) / avg(pop_mean(not i))) * mu(phi,i)
				matMeanPredCV(intTestNeuron,vecTestTrials) = (vecPopMeanTest ./ dblAvgPopMeanTrain) .* vecTrainR;
				
				%% same gain for all stims
				%calculate gain axis for this stimulus from training trials
				vecTrainGain1Ax = mean(matTrainOtherR,2);
				
				%get population test responses and project onto training gain axis
				vecOtherGain1=getProjOnLine(matTestOtherR,vecTrainGain1Ax);
				
				%get mean resp during train
				vecTrainR = matTrainTuning(intTestNeuron,vecTestS);
				
				%predict: R(phi,i,t)= (gain1(not i,t) / norm(mu(not i))) * mu(phi,i)
				matGain1PredCV(intTestNeuron,vecTestTrials) = (vecOtherGain1 ./ norm(vecTrainGain1Ax)) .* vecTrainR';
				
				%% gain prediction: stim-specific coupling between population gain and predicted neuron
				%R(phi,i,t)= (gain(phi,not i,t) / norm(mu(phi,not i))) * mu(phi,i)
				for intStimIdx=1:intStimNum
					%calculate gain axis for this stimulus from training trials
					indTrainT = vecTrainS==intStimIdx;
					matTrainOtherR = matTrainR(indOtherNeurons,indTrainT);
					vecTrainGainAx = mean(matTrainOtherR,2);
					
					%get population test responses and project onto training gain axis
					indTestT = vecTestS==intStimIdx;
					matTestOtherR = matTestR(indOtherNeurons,indTestT);
					vecOtherGain=getProjOnLine(matTestOtherR,vecTrainGainAx);
					
					%get mean response for this stimulus for this neuron over training trials
					dblTrainSelf = mean(matTrainR(~indOtherNeurons,indTrainT));
					
					%predict: R(phi,i,t)= (gain(phi,not i,t) / norm(mu(phi,not i))) * mu(phi,i)
					matGainPredCV(intTestNeuron,vecTestTrials(indTestT)) = (vecOtherGain ./ norm(vecTrainGainAx)) * dblTrainSelf;
				end
				
			end
		end
		
		
		[dblR2_Tune,dblSS_tot,dblSS_res,dblT_Tune,dblP_Tune,dblR2_adjusted_Tune,dblR2_SE_Tune] = getR2(matMeanRate(:),matTunePredCV(:),numel(matTrainTuning));
		[dblR2_Mean,dblSS_tot,dblSS_res,dblT_Mean,dblP_Mean,dblR2_adjusted_Mean,dblR2_SE_Mean] = getR2(matMeanRate(:),matMeanPredCV(:),numel(matTrainTuning));
		[dblR2_Gain,dblSS_tot,dblSS_res,dblT_Gain,dblP_Gain,dblR2_adjusted_Gain,dblR2_SE_Gain] = getR2(matMeanRate(:),matGainPredCV(:),numel(matTrainTuning));
		[dblR2_Gain1,dblSS_tot,dblSS_res,dblT_Gain1,dblP_Gain1,dblR2_adjusted_Gain1,dblR2_SE_Gain1] = getR2(matMeanRate(:),matGain1PredCV(:),numel(matTrainTuning));
		
		%prediction per neuron
		vecR2_Tune = nan(1,intNeuronNum);
		vecR2_Mean = nan(1,intNeuronNum);
		vecR2_Gain = nan(1,intNeuronNum);
		vecR2_Gain1 = nan(1,intNeuronNum);
		for intN=1:intNeuronNum
			%save prediction
			[dblR2_TuneN,dblSS_tot,dblSS_res,dblT,dblP,dblR2_adjusted_TuneN] = getR2(matMeanRate(intN,:),matTunePredCV(intN,:),numel(matTrainTuning(intN,:)));
			[dblR2_MeanN,dblSS_tot,dblSS_res,dblT,dblP,dblR2_adjusted_MeanN] = getR2(matMeanRate(intN,:),matMeanPredCV(intN,:),numel(matTrainTuning(intN,:)));
			[dblR2_GainN,dblSS_tot,dblSS_res,dblT,dblP,dblR2_adjusted_GainN] = getR2(matMeanRate(intN,:),matGainPredCV(intN,:),numel(matTrainTuning(intN,:)));
			[dblR2_Gain1N,dblSS_tot,dblSS_res,dblT,dblP,dblR2_adjusted_Gain1N] = getR2(matMeanRate(intN,:),matGain1PredCV(intN,:),numel(matTrainTuning(intN,:)));
			vecR2_Tune(intN) = dblR2_TuneN;
			vecR2_Mean(intN) = dblR2_MeanN;
			vecR2_Gain(intN) = dblR2_GainN;
			vecR2_Gain1(intN) = dblR2_Gain1N;
		end
		dblPredPerNeuronTune = mean(vecR2_Tune);
		dblPredPerNeuronMean = mean(vecR2_Mean);
		dblPredPerNeuronGain = mean(vecR2_Gain);
		dblPredPerNeuronGain1 = mean(vecR2_Gain1);
		
		%%
		dblStep = 4;
		vecSelfR = matMeanRate(:);
		vecOrthR = matMeanPredCV(:);
		vecBinsX = 0:dblStep:81 - 0.5*dblStep;
		vecBinsY = 0:dblStep:80 - 0.5*dblStep;
		matCounts = histcounts2(vecSelfR,vecOrthR,vecBinsX,vecBinsY);
		vecLimC = [0 80];
		vecLimC_d = [-40 40];
		if boolMakeFigs
			figure;maxfig
			%real
			colormap(redwhite);
			subplot(2,4,1)
			imagesc(matMeanRate,vecLimC)
			colorbar;
			ylabel('Neuron');
			xlabel('Trial');
			title(sprintf('%s; Recorded spiking',strType));
			h=colorbar;
			ylabel(h,'Spiking rate (Hz)');
			fixfig;grid off;
			
			%summary
			subplot(2,4,5)
			errorbar(1:4,[dblR2_Tune dblR2_Mean dblR2_Gain1 dblR2_Gain],[dblR2_SE_Tune dblR2_SE_Mean dblR2_SE_Gain1 dblR2_SE_Gain],'x','CapSize',20)
			set(gca,'xtick',1:4,'xticklabel',{'Tune','Tune*Mean','Tune*Gain1','Tune*Gain'});
			xtickangle(20);
			ylabel('CV pop resp predictability (R^2)');
			xlim([0 5]);
			title(sprintf('%s',strTypeRec),'interpreter','none');
			fixfig;
			
			%tuning only
			subplot(2,4,2)
			imagesc(matTunePredCV,vecLimC)
			colorbar;
			ylabel('Neuron');
			xlabel('Trial');
			title(sprintf('CV model: R = ori tuning'));
			h=colorbar;
			ylabel(h,'CV-predicted rate (Hz)');
			fixfig;grid off;
			
			h=subplot(2,4,6);
			imagesc(matMeanRate-matTunePredCV,vecLimC_d)
			h.Colormap=redblue;
			colorbar;
			ylabel('Neuron');
			xlabel('Trial');
			title(sprintf('Resids; adj-R^2=%.3f; T=%.1f',dblR2_adjusted_Tune,dblT_Tune));
			h=colorbar;
			ylabel(h,'Prediction error (Hz)');
			fixfig;grid off;
			
			%tune*mean
			subplot(2,4,3)
			imagesc(matMeanPredCV,vecLimC)
			colorbar;
			ylabel('Neuron');
			xlabel('Trial');
			title(sprintf('CV model: R = ori tuning * pop-mean'));
			h=colorbar;
			ylabel(h,'CV-predicted rate (Hz)');
			fixfig;grid off;
			
			h=subplot(2,4,7);
			imagesc(matMeanRate-matMeanPredCV,vecLimC_d)
			h.Colormap=redblue;
			colorbar;
			ylabel('Neuron');
			xlabel('Trial');
			title(sprintf('Resids; adj-R^2=%.3f; T=%.1f',dblR2_adjusted_Mean,dblT_Mean));
			h=colorbar;
			ylabel(h,'Prediction error (Hz)');
			fixfig;grid off;
			
			%tune*gain
			subplot(2,4,4)
			imagesc(matGainPredCV,vecLimC)
			colorbar;
			ylabel('Neuron');
			xlabel('Trial');
			title(sprintf('CV model: R = ori tuning * pop-gain'));
			h=colorbar;
			ylabel(h,'CV-predicted rate (Hz)');
			fixfig;grid off;
			
			h=subplot(2,4,8);
			imagesc(matMeanRate-matGainPredCV,vecLimC_d)
			h.Colormap=redblue;
			colorbar;
			ylabel('Neuron');
			xlabel('Trial');
			title(sprintf('Resids; adj-R^2=%.3f; T=%.1f',dblR2_adjusted_Gain,dblT_Gain));
			h=colorbar;
			ylabel(h,'Prediction error (Hz)');
			fixfig;grid off;
			
			if boolSaveFigs
				%% save fig
				export_fig(fullpath(strFigurePathSR,sprintf('QC2_PopRespPrediction_%s.tif',strTypeRec)));
				export_fig(fullpath(strFigurePathSR,sprintf('QC2_PopRespPrediction_%s.pdf',strTypeRec)));
			end
		end
		
		
		%% project
		%for a given cloud of pop responses (n=trials), what fraction of noise aligns with the
		%mean-rate axis? what proportion with f'? how does this compare to a random direction (=spectral decomposition of covariance matrix)?
		
		%run projections
		fprintf('   Running multi-dim analysis [%s]\n',getTime);
		vecGainAx = mean(matMeanRate,2);
		vecGainAx = vecGainAx./norm(vecGainAx);
		
		matGainResp = nan(size(matMeanRate));
		matCov = nan(intNeuronNum,intNeuronNum,intStimNum);
		vecSdProjMean = nan(1,intStimNum);
		vecSdProjGain1 = nan(1,intStimNum);
		vecSdProjGain = nan(1,intStimNum);
		vecSdProjRand = nan(1,intStimNum);
		vecSdProjOrth = nan(1,intStimNum);
		vecSdProjAdja = nan(1,intStimNum);
		vecDistFromOrigin = nan(1,intStimNum);
		vecDistFromOrthOri = nan(1,intStimNum);
		vecDistFromAdjaOri = nan(1,intStimNum);
		%reflections
		vecReflMuDistToOrth = nan(1,intStimNum);
		vecReflMuDistToSelf = nan(1,intStimNum);
		vecReflGainDistToOrth = nan(1,intStimNum);
		vecReflGainDistToSelf = nan(1,intStimNum);
		vecReflRandDistToOrth = nan(1,intStimNum);
		vecReflRandDistToSelf = nan(1,intStimNum);
		
		%cossim
		vecAngleMeanAndGain = nan(1,intStimNum);
		vecAngleMeanAndGainRand = nan(1,intStimNum);
		
		for intStimIdx=1:intStimNum
			%take stim i
			vecOriTrials = find(vecStimIdx==intStimIdx);
			%i+1 = counterclockwise adjacent
			indAdjaOriTrials = vecStimIdx==modx(intStimIdx+1,intStimNum);
			%i-1 = clockwise adjacent
			indAdja2OriTrials = vecStimIdx==modx(intStimIdx-1,intStimNum);
			%orthogonal (n/2 is opposite, but stim vector is on 180 degrees half unit circle)
			indOrthOriTrials = vecStimIdx==modx(intStimIdx+ceil(intStimNum/2)-1,intStimNum);
			
			%get responses for stim i
			matPoints =matMeanRate(:,vecOriTrials);
			%get mean responses over repetitions for stims i, adjacent and orthogonal
			vecOriR = mean(matPoints,2);
			vecAdjaOriR = mean(matMeanRate(:,indAdjaOriTrials),2);
			vecAdja2OriR = mean(matMeanRate(:,indAdja2OriTrials),2);
			vecOrthOriR = mean(matMeanRate(:,indOrthOriTrials),2);
			
			%project responses onto mean-axis (diagonal)
			vecProjMean=getProjOnLine(matPoints,ones(intNeuronNum,1));
			%get variability (sd) over trials along this mean-axis
			vecSdProjMean(intStimIdx) = std(vecProjMean);
			%project responses onto random axis taken from the nonnegative orthant (hyperquadrant)
			vecProjRand=getProjOnLine(matPoints,(rand(intNeuronNum,1)-0.5)*2);
			%get variability (sd) over trials along this random-axis
			vecSdProjRand(intStimIdx) = std(vecProjRand);
			
			%Project activity unto single gain axis (Gain1) versus one per stim (Gain)
			vecProjGain=getProjOnLine(matPoints,vecOriR);
			vecSdProjGain(intStimIdx) = std(vecProjGain);
			
			%for avg-gain; is the same for all stimuli (average over all stims)
			vecProjGain1=getProjOnLine(matPoints,vecGainAx);
			vecSdProjGain1(intStimIdx) = std(vecProjGain1);
			
			%project onto axis connecting i with adja-i or orth-i: how variable is resp to stim
			%i along f' (vecOriR-vecAdjaOriR) compared with axis connecting (vecOriR-vecOrthOriR)
			vecProjAdja=getProjOnLine(matPoints-vecAdjaOriR,vecOriR-vecAdjaOriR);
			vecSdProjAdja(intStimIdx) = std(vecProjAdja);
			vecProjOrth=getProjOnLine(matPoints-vecOrthOriR,vecOriR-vecOrthOriR);
			vecSdProjOrth(intStimIdx) = std(vecProjOrth);
			
			%what is the distance of the average response of stim i to the origin?
			vecDistFromOrigin(intStimIdx) = norm(vecOriR);
			%what is the distance of the avg resp of stim i to the avg resp of orth-i ?
			vecDistFromOrthOri(intStimIdx) = norm(vecOriR-vecOrthOriR);
			%what is the distance of the avg resp of stim i to the avg resp of adja-i ?
			vecDistFromAdjaOri(intStimIdx) = norm(vecOriR-vecAdjaOriR);%/2 + norm(vecOriR-vecAdja2OriR)/2;
			
			%what is the alignment of the gain axis with the mean axis?
			vecGainDir = vecOriR;
			vecMeanDir = ones(size(vecOriR));
			[dblCosSim,dblAngSim] = cossim(vecGainDir,vecMeanDir);
			vecAngleMeanAndGain(intStimIdx) = dblCosSim;
			
			%what is the alignment of a random axis with the mean axis?
			vecRandSim = nan(1,100);
			for intIter=1:numel(vecRandSim)
				[dblCosSimRand,dblAngSimRand] = cossim(rand(intNeuronNum,1),vecMeanDir);
				vecRandSim(intIter) = dblCosSimRand;
			end
			vecAngleMeanAndGainRand(intStimIdx) = mean(vecRandSim);
			
			%% reflect repetition-averaged population responses around gain-axis and mean-axis
			%is the manifold symmetric around the mean or around the gain?
			intIters=100;
			vecTempReflMuDistToOrth = nan(1,intIters);
			vecTempReflMuDistToSelf = nan(1,intIters);
			vecTempReflGainDistToOrth = nan(1,intIters);
			vecTempReflGainDistToSelf = nan(1,intIters);
			vecTempReflRandDistToOrth = nan(1,intIters);
			vecTempReflRandDistToSelf = nan(1,intIters);
			for intIter=1:100
				%% select n/2 trials from stim i
				vecReorder = randperm(numel(vecOriTrials));
				%vecOriTrials1 = vecOriTrials;%vecReorder(1:(floor(numel(vecReorder)/2)));
				%vecOriTrials2 = indAdjaOriTrials;%vecReorder(((ceil(numel(vecReorder)/2))+1):end);
				%random 50% of stim i
				vecOriTrials1 = vecOriTrials(vecReorder(1:(floor(numel(vecReorder)/2))));
				vecSelfR = mean(matMeanRate(:,vecOriTrials1),2);
				%other 50% of stim i
				vecOriTrials2 = vecOriTrials(vecReorder(((ceil(numel(vecReorder)/2))+1):end));
				vecSelfR_CV = mean(matMeanRate(:,vecOriTrials2),2);
				%all trials of orth-i
				indOrthOriTrials = vecStimIdx==modx(intStimIdx+ceil(intStimNum/2),intStimNum);
				vecOrthR = mean(matMeanRate(:,indOrthOriTrials),2);
				
				%reflect around mean
				[dummy,vecXprime]=getProjOnLine(vecSelfR,ones(intNeuronNum,1));
				vecXrefl = 2*vecXprime - vecSelfR;
				%compare with orthogonal
				dblReflMuDistToOrth = norm(vecXrefl-vecOrthR);
				%compare with self
				dblReflMuDistToSelfCV = norm(vecXrefl-vecSelfR_CV);
				
				
				%reflect around gain axis
				[dummy,vecXprimeG]=getProjOnLine(vecSelfR,vecGainAx);
				vecXreflG = 2*vecXprimeG - vecSelfR;
				%compare with orthogonal
				dblReflGainDistToOrth = norm(vecXreflG-vecOrthR);
				%compare with self
				dblReflGainDistToSelf = norm(vecXreflG-vecSelfR_CV);
				
				%reflect around random axis
				vecRandAx = rand(intNeuronNum,1);
				[dummy,vecXprimeR]=getProjOnLine(vecSelfR,vecRandAx);
				vecXreflR = 2*vecXprimeR - vecSelfR;
				%compare with orthogonal
				dblReflRandDistToOrth = norm(vecXreflR-vecOrthR);
				%compare with self
				dblReflRandDistToSelf = norm(vecXreflR-vecSelfR_CV);
				
				%save
				vecTempReflMuDistToOrth(intIter) = dblReflMuDistToOrth;
				vecTempReflMuDistToSelf(intIter) = dblReflMuDistToSelfCV;
				vecTempReflGainDistToOrth(intIter) = dblReflGainDistToOrth;
				vecTempReflGainDistToSelf(intIter) = dblReflGainDistToSelf;
				vecTempReflRandDistToOrth(intIter) = dblReflRandDistToOrth;
				vecTempReflRandDistToSelf(intIter) = dblReflRandDistToSelf;
			end
			
			%mean over resamples
			vecReflMuDistToOrth(intStimIdx) = mean(vecTempReflMuDistToOrth);
			vecReflMuDistToSelf(intStimIdx) = mean(vecTempReflMuDistToSelf);
			vecReflGainDistToOrth(intStimIdx) = mean(vecTempReflGainDistToOrth);
			vecReflGainDistToSelf(intStimIdx) = mean(vecTempReflGainDistToSelf);
			vecReflRandDistToOrth(intStimIdx) = mean(vecTempReflRandDistToOrth);
			vecReflRandDistToSelf(intStimIdx) = mean(vecTempReflRandDistToSelf);
			
			
			%% save data as scaled ori response & covariance
			matGainResp(:,vecOriTrials) = vecOriR * vecNormPopRate(vecOriTrials);
			matCov(:,:,intStimIdx) = cov(matMeanRate(:,vecOriTrials)');
		end
		
		%symmetry
		vecGainSymmetry = 100 - (vecReflGainDistToOrth./vecReflGainDistToSelf)*100;
		vecMeanSymmetry = 100 - (vecReflMuDistToOrth./vecReflMuDistToSelf)*100;
		vecRandSymmetry = 100 - (vecReflRandDistToOrth./vecReflRandDistToSelf)*100;
		
		%% plot
		[h,pMeanGain]=ttest(vecSdProjMean,vecSdProjGain);
		[h,pMeanOrth]=ttest(vecSdProjMean,vecSdProjOrth);
		[h,pGainOrth]=ttest(vecSdProjGain,vecSdProjOrth);
		if boolMakeFigs
			figure;maxfig;
			subplot(2,3,1)
			hold on
			vecLoc = [0.2 0.5 0.8 1.1];
			dblJit = 0.05;
			scatter(vecLoc(1)+(rand(1,intStimNum)-0.5)*2*dblJit,vecSdProjMean,[200],[0.5 0.5 0.5],'.')
			scatter(vecLoc(2)+(rand(1,intStimNum)-0.5)*2*dblJit,vecSdProjGain1,[200],[0.5 0.5 0.5],'.')
			scatter(vecLoc(3)+(rand(1,intStimNum)-0.5)*2*dblJit,vecSdProjGain,[200],[0.5 0.5 0.5],'.')
			%scatter(vecLoc(3)+(rand(1,intStimNum)-0.5)*2*dblJit,vecSdProjRand,[200],[0.5 0.5 0.5],'.')
			%scatter(vecLoc(3)+(rand(1,intStimNum)-0.5)*2*dblJit,vecSdProjAdja,[200],[0.5 0.5 0.5],'.')
			scatter(vecLoc(4)+(rand(1,intStimNum)-0.5)*2*dblJit,vecSdProjOrth,[200],[0.5 0.5 0.5],'.')
			errorbar(vecLoc(1),mean(vecSdProjMean),std(vecSdProjMean)/sqrt(intStimNum),'bx','CapSize',20)
			errorbar(vecLoc(2),mean(vecSdProjGain1),std(vecSdProjGain1)/sqrt(intStimNum),'kx','CapSize',20)
			errorbar(vecLoc(3),mean(vecSdProjGain),std(vecSdProjGain)/sqrt(intStimNum),'kx','CapSize',20)
			%errorbar(vecLoc(3),mean(vecSdProjRand),std(vecSdProjRand)/sqrt(intStimNum),'rx','CapSize',20)
			%errorbar(vecLoc(3),mean(vecSdProjAdja),std(vecSdProjAdja)/sqrt(intStimNum),'kx','CapSize',20)
			errorbar(vecLoc(4),mean(vecSdProjOrth),std(vecSdProjOrth)/sqrt(intStimNum),'rx','CapSize',20)
			hold off;
			set(gca,'xtick',vecLoc,'xticklabel',{'Pop mean','Pop gain1','Pop gain','Stim-Orth'});
			xlabel('Projection axis');
			ylabel('Range of spiking rates (\sigmaHz)');
			title(sprintf('%s; MG,p=%.2e; MO,p=%.3f; GO,p=%.2e',strType,pMeanGain,pMeanOrth,pGainOrth));
			fixfig;
			ylim(gca,[0 max(get(gca,'ylim'))]);
			xlim([vecLoc(1)-0.2 vecLoc(end)+0.2]);
			
			%what is distance from origin to mean R compared with distance between ori Rs?
			subplot(2,3,2);
			hold on
			vecLoc = [0.2 0.5 0.8];
			dblJit = 0.05;
			scatter(vecLoc(1)+(rand(1,intStimNum)-0.5)*2*dblJit,vecDistFromOrigin,[200],[0.5 0.5 0.5],'.')
			errorbar(vecLoc(1),mean(vecDistFromOrigin),std(vecDistFromOrigin)/sqrt(intStimNum),'kx','CapSize',20)
			scatter(vecLoc(2)+(rand(1,intStimNum)-0.5)*2*dblJit,vecDistFromOrthOri,[200],[0.5 0.5 0.5],'.')
			errorbar(vecLoc(2),mean(vecDistFromOrthOri),std(vecDistFromOrthOri)/sqrt(intStimNum),'bx','CapSize',20)
			scatter(vecLoc(3)+(rand(1,intStimNum)-0.5)*2*dblJit,vecDistFromAdjaOri,[200],[0.5 0.5 0.5],'.')
			errorbar(vecLoc(3),mean(vecDistFromAdjaOri),std(vecDistFromAdjaOri)/sqrt(intStimNum),'rx','CapSize',20)
			hold off
			set(gca,'xtick',vecLoc,'xticklabel',{'Origin','Orth. ori','Adja. ori'});
			%ylim([0 60]);
			xlim([vecLoc(1)-0.2 vecLoc(end)+0.2]);
			xlabel('Distance of pop. stim resp to');
			ylabel('Pop. spiking rate distance (\DeltaHz)');
			title('Mu +/- sem, n=12 orientations');
			fixfig;
			
			%cos sim of gain with mean vs rand with mean
			%note that random vectors do not lead to cossim of 0 because the space is limited to
			%positive values only (Hz > 0)
			[h,pCosSim] = ttest2(vecAngleMeanAndGain,vecAngleMeanAndGainRand);
			subplot(2,3,3);
			hold on
			vecLoc = [0.2 0.8];
			dblJit = 0.05;
			scatter(vecLoc(1)+(rand(1,intStimNum)-0.5)*2*dblJit,vecAngleMeanAndGain,[200],[0.5 0.5 0.5],'.')
			errorbar(vecLoc(1),mean(vecAngleMeanAndGain),std(vecAngleMeanAndGain)/sqrt(intStimNum),'kx','CapSize',20)
			scatter(vecLoc(2)+(rand(1,intStimNum)-0.5)*2*dblJit,vecAngleMeanAndGainRand,[200],[0.5 0.5 0.5],'.')
			errorbar(vecLoc(2),mean(vecAngleMeanAndGainRand),std(vecAngleMeanAndGainRand)/sqrt(intStimNum),'bx','CapSize',20)
			hold off
			set(gca,'xtick',vecLoc,'xticklabel',{'Gain','Rand'});
			xlim([vecLoc(1)-0.2 vecLoc(end)+0.2]);
			xlabel('Neural axis');
			ylabel('Alignment with mean-axis (cos sim)');
			title(sprintf('T-test, cos-sim mean with gain vs rand: p=%.2e',pCosSim));
			fixfig;
			
			%symmetry around mean or gain axis
			vecLoc = [0.2 0.8 2.2 2.8 4.2 4.8];
			dblJit = 0.05;
			subplot(2,3,4);
			hold on
			errorbar(vecLoc(1),mean(vecReflMuDistToOrth),std(vecReflMuDistToOrth)/sqrt(intStimNum),'kx','CapSize',20)
			errorbar(vecLoc(2),mean(vecReflMuDistToSelf),std(vecReflMuDistToSelf)/sqrt(intStimNum),'bx','CapSize',20)
			errorbar(vecLoc(3),mean(vecReflGainDistToOrth),std(vecReflGainDistToOrth)/sqrt(intStimNum),'kx','CapSize',20)
			errorbar(vecLoc(4),mean(vecReflGainDistToSelf),std(vecReflGainDistToSelf)/sqrt(intStimNum),'bx','CapSize',20)
			errorbar(vecLoc(5),mean(vecReflRandDistToOrth),std(vecReflRandDistToOrth)/sqrt(intStimNum),'kx','CapSize',20)
			errorbar(vecLoc(6),mean(vecReflRandDistToSelf),std(vecReflRandDistToSelf)/sqrt(intStimNum),'bx','CapSize',20)
			hold off
			set(gca,'xtick',vecLoc,'xticklabel',{'Mu/Orth','Mu/Self','Gain/Orth','Gain/Self','Rand/Orth','Rand/Self'});
			xlim([vecLoc(1)-0.2 vecLoc(end)+0.2]);
			xlabel('Reflection over/stim ori');
			ylabel('Prediction error after reflection (\DeltaHz)');
			fixfig;
			
			% symmetry around gain axis
			vecLoc = 1:3;
			subplot(2,3,5);
			[h,pSymMuG]=ttest2(vecMeanSymmetry,vecGainSymmetry);
			[h,pSymMuR]=ttest2(vecMeanSymmetry,vecRandSymmetry);
			[h,pSymGainR]=ttest2(vecRandSymmetry,vecGainSymmetry);
			hold on
			errorbar(1,mean(vecMeanSymmetry),std(vecMeanSymmetry)/sqrt(intStimNum),'bd','CapSize',20)
			errorbar(2,mean(vecGainSymmetry),std(vecGainSymmetry)/sqrt(intStimNum),'ko','CapSize',20)
			errorbar(3,mean(vecRandSymmetry),std(vecRandSymmetry)/sqrt(intStimNum),'x','color',[0.2 0 0.8],'CapSize',20)
			hold off
			ylabel('Manifold symmetry (%)');
			set(gca,'xtick',vecLoc,'xticklabel',{'Pop. mean ax','Pop. gain ax','Rand axis'});
			xlabel('Axis of symmetry');
			xlim([vecLoc(1)-0.2 vecLoc(end)+0.2]);
			title(sprintf('T-test p, Mu-G=%.2e, Mu-R=%.2e, R-G=%.2e',pSymMuG,pSymMuR,pSymGainR));
			fixfig;
			
			%make example figures of the above
			intS1 = 1;
			intS1Adja = 2;
			intS2 = 1+round(intStimNum/2);
			
			dblDistToGainCenter = mean(vecDistFromOrigin);
			dblAlignment=cossim(vecGainAx,ones(size(vecGainAx)));
			dblAng = acos(dblAlignment)/2+deg2rad(45);
			vecGainCenter = [cos(dblAng) sin(dblAng)]*dblDistToGainCenter;
			
			dblAng1 = dblAng + atan((mean(vecDistFromOrthOri(intS1)))/dblDistToGainCenter);
			dblAng1A = dblAng + atan((mean(vecDistFromOrthOri(intS1)) - mean(vecDistFromAdjaOri(intS1))/4)/dblDistToGainCenter);
			dblAng2 = dblAng - atan((mean(vecDistFromOrthOri(intS2)))/dblDistToGainCenter);
			
			vecCenter1 = [cos(dblAng1) sin(dblAng1)]*vecDistFromOrigin(intS1);
			vecCenter1A = [cos(dblAng1A) sin(dblAng1A)]*vecDistFromOrigin(intS1Adja);
			vecCenter2 = [cos(dblAng2) sin(dblAng2)]*vecDistFromOrigin(intS2);
			
			%reflect
			vecSelfR = vecCenter1';
			[dummy,vecXprimeG]=getProjOnLine(vecSelfR,vecGainCenter');
			vecXreflG = 2*vecXprimeG - vecSelfR;
			[dummy,vecXprime]=getProjOnLine(vecSelfR,ones(size(vecSelfR)));
			vecXrefl = 2*vecXprime - vecSelfR;
			
			subplot(2,3,6)%[3 4 7 8]);
			cla;
			hold on
			%mean
			dblMaxLim = ceil((max([vecCenter1 vecCenter2]) + max(vecSdProjGain))/5)*5;
			intOffset = 1;
			plot([0 dblMaxLim],[0 dblMaxLim],'-','color',[0.5 0.5 0.5])
			text(dblMaxLim*0.2+intOffset,dblMaxLim*0.2-intOffset,'Pop. mean axis','FontSize',16,'Color',[0.5 0.5 0.5],'Rotation',45)
			%gain
			plot([0 vecGainCenter(1)*1.4],[0 vecGainCenter(2)*1.4],'-','color','k')
			text(vecGainCenter(1)*0.3+intOffset,vecGainCenter(2)*0.3-intOffset,'Pop. gain axis','FontSize',16,'Color','k','Rotation',rad2deg(dblAng))
			
			%centers
			scatter(vecCenter1(1),vecCenter1(2),'rx');
			scatter(vecCenter2(1),vecCenter2(2),'bx');
			scatter(vecCenter1A(1),vecCenter1A(2),[],[0.5 0 0.5],'x');
			
			
			%covars
			ellipse(vecCenter1(1),vecCenter1(2),mean(vecSdProjGain(intS1)),mean(vecSdProjOrth(intS1)),dblAng1,'Color','r','LineStyle','-');
			ellipse(vecCenter2(1),vecCenter2(2),mean(vecSdProjGain(intS2)),mean(vecSdProjOrth(intS2)),dblAng2,'Color','b','LineStyle','-');
			
			%draw reflections
			plot([vecCenter1(1) vecXreflG(1)],[vecCenter1(2) vecXreflG(2)],'k--');
			plot([vecCenter1(1) vecXrefl(1)],[vecCenter1(2) vecXrefl(2)],'--','color',[0.5 0.5 0.5]);
			scatter(vecXreflG(1),vecXreflG(2),'ro');
			scatter(vecXrefl(1),vecXrefl(2),'rd');
			
			%text
			text(vecCenter1(1)+intOffset,vecCenter1(2)+intOffset,'Ref','FontSize',16,'Color','r')
			text(vecCenter2(1)+intOffset,vecCenter2(2)+intOffset,'Orth','FontSize',16,'Color','b')
			text(vecCenter1A(1)+intOffset,vecCenter1A(2)+intOffset,'Adja','FontSize',16,'Color',[0.5 0 0.5])
			text(vecXreflG(1)+intOffset,vecXreflG(2)-intOffset,'Gain-reflect','FontSize',16,'Color','r')
			text(vecXrefl(1)+intOffset,vecXrefl(2)+intOffset,'Mean-reflect','FontSize',16,'Color','r')
			
			%finish fig
			hold off
			axis equal;
			%xlim([0 dblMaxLim*(2/3)]);
			%ylim([0 dblMaxLim]);
			xlabel('Spiking rate axis X (Hz)')
			ylabel('Spiking rate axis Y (Hz)')
			fixfig;grid off;
			
			if boolSaveFigs
				%% save fig
				export_fig(fullpath(strFigurePathSR,sprintf('QC2_Projections_%s.tif',strTypeRec)));
				export_fig(fullpath(strFigurePathSR,sprintf('QC2_Projections_%s.pdf',strTypeRec)));
			end
		end
		
		%% predict gain from pupil
		%get pupil size per trial
		vecPupilSize = [];
		vecTrialMean = [];
		vecTrialMean_PupilHat = [];
		vecTrialMean_RunningHat = [];
		vecTrialGain = [];
		vecTrialGain_PupilHat = [];
		vecTrialGain_RunningHat = [];
		
		if isfield(sThisRec,'Pupil') && numel(vecPupilStimOn) == intTrialNum
			sPupil = sThisRec.Pupil;
			vecTime = sPupil.vecTime;
			vecVals = sPupil.vecRadius;
			vecEvents = vecPupilStimOn;
			dblDur = 1.5;
			vecPupMov = [0 sqrt(diff(sPupil.vecCenterX).^2 + diff(sPupil.vecCenterY).^2)];
			vecPupilMov = mean(getRespMat(sPupil.vecTime,vecPupMov,vecPupilStimOn-0.5,[0 dblDur]),2);
			vecPupilSize = mean(getRespMat(sPupil.vecTime,sPupil.vecRadius2,vecPupilStimOn-0.5,[0 dblDur]),2);
			vecPupilBlinks = mean(getRespMat(sPupil.vecTime,sPupil.vecBlinks,vecPupilStimOn-0.5,[0 dblDur]),2);
			vecPupilIsEdited = mean(getRespMat(sPupil.vecTime,sPupil.vecIsEdited,vecPupilStimOn-0.5,[0 dblDur]),2);
			[varDataOut,vecUnique,vecCounts]=val2idx(vecPupilSize);
			indRepeatedVals = ismember(varDataOut,find(vecCounts>2));
			indRem = vecPupilBlinks>0.1 | indRepeatedVals;
			
			vecTrialGain = getProjOnLine(matMeanRate,vecGainAx);
			vecTrialMean = getProjOnLine(matMeanRate,ones(size(vecGainAx))./norm(ones(size(vecGainAx))));
			%vecTrialMean = mean(matMeanRate,1)';
			indUseTrials = ~isnan(vecPupilSize) & ~isnan(vecTrialMean) & ~isnan(vecTrialGain) & ~indRem;
			[rGainPupilS,PG,RGL,RGU] = corrcoef(vecTrialGain(indUseTrials),vecPupilSize(indUseTrials));
			[rMeanPupilS,PM,RML,RMU] = corrcoef(vecTrialMean(indUseTrials),vecPupilSize(indUseTrials));
			dblRG = rGainPupilS(1,2);
			dblRG_UB = RGU(1,2)-dblRG;
			dblRG_LB = dblRG-RGL(1,2);
			dblPG = PG(1,2);
			
			dblRM = rMeanPupilS(1,2);
			dblRM_UB = RMU(1,2)-dblRM;
			dblRM_LB = dblRM-RML(1,2);
			dblPM = PM(1,2);
			
			figure;maxfig;
			subplot(2,3,1);
			scatter(vecTrialGain(indUseTrials),vecPupilSize(indUseTrials),'b.');
			xlabel('Trial gain');
			ylabel('Pupil size');
			title(sprintf('%s; Gain prediction from pupil size',strTypeRec),'interpreter','none');
			fixfig;
			
			subplot(2,3,2);
			scatter(vecTrialMean(indUseTrials),vecPupilSize(indUseTrials),'k.');
			xlabel('Trial mean');
			ylabel('Pupil size');
			title('Mean prediction from pupil size');
			fixfig;
			
			subplot(2,3,3);
			hold on
			errorbar(0.2,dblRG,dblRG_LB,dblRG_UB,'bx');
			errorbar(0.8,dblRM,dblRM_LB,dblRM_UB,'kx');
			set(gca,'xtick',[0.2 0.8],'xticklabel',{'Gain','Mean'});
			ylabel('Pupil correlation (Pearson r)');
			xlim([0 1]);
			dblMax = max(get(gca,'ylim'));
			if dblMax > 0,ylim([0 dblMax]);end
			title(sprintf('T-tests vs 0: gain-p=%.2e,mean-p=%.2e',dblPG,dblPM));
			fixfig;
			
			if boolSaveFigs
				%% save fig
				export_fig(fullpath(strFigurePathSR,sprintf('QC2_PupilPredGain_%s.tif',strTypeRec)));
				export_fig(fullpath(strFigurePathSR,sprintf('QC2_PupilPredGain_%s.pdf',strTypeRec)));
			end
		else
			dblRG = nan;
			dblRM = nan;
		end
		
		%% save data
		if boolSaveData
			%save data
			%% predictions
			sPrediction = struct;
			sPrediction.matMeanRate = matMeanRate;
			sPrediction.matTunePredCV = matTunePredCV;
			sPrediction.matMeanPredCV = matMeanPredCV;
			sPrediction.matGainPredCV = matGainPredCV;
			sPrediction.matGain1PredCV = matGain1PredCV;
			
			sPrediction.dblR2_Tune = dblR2_Tune;
			sPrediction.dblR2_Mean = dblR2_Mean;
			sPrediction.dblR2_Gain = dblR2_Gain;
			sPrediction.dblR2_Gain1 = dblR2_Gain1;
			
			sPrediction.vecR2_Tune = vecR2_Tune;
			sPrediction.vecR2_Mean = vecR2_Mean;
			sPrediction.vecR2_Gain = vecR2_Gain;
			sPrediction.vecR2_Gain1 = vecR2_Gain1;
			
			sPrediction.dblPredPerNeuronTune = dblPredPerNeuronTune;
			sPrediction.dblPredPerNeuronMean = dblPredPerNeuronMean;
			sPrediction.dblPredPerNeuronGain = dblPredPerNeuronGain;
			sPrediction.dblPredPerNeuronGain1 = dblPredPerNeuronGain1;
			
			%pupil & running
			sPrediction.vecTrialGain_PupilHat = vecTrialGain_PupilHat;
			sPrediction.vecTrialMean_PupilHat = vecTrialMean_PupilHat;
			sPrediction.vecTrialGain_RunningHat = vecTrialGain_RunningHat;
			sPrediction.vecTrialMean_RunningHat = vecTrialMean_RunningHat;
			sPrediction.vecTrialGain = vecTrialGain;
			sPrediction.vecTrialMean = vecTrialMean;
			
			%% projections
			sProjection = struct;
			sProjection.vecReflMuDistToOrth = vecReflMuDistToOrth;
			sProjection.vecReflMuDistToSelf = vecReflMuDistToSelf;
			sProjection.vecReflGainDistToOrth = vecReflGainDistToOrth;
			sProjection.vecReflGainDistToSelf = vecReflGainDistToSelf;
			sProjection.vecReflRandDistToOrth = vecReflRandDistToOrth;
			sProjection.vecReflRandDistToSelf = vecReflRandDistToSelf;
			
			sProjection.vecGainSymmetry = vecGainSymmetry;
			sProjection.vecMeanSymmetry = vecMeanSymmetry;
			sProjection.vecRandSymmetry = vecRandSymmetry;
			
			sProjection.vecAngleMeanAndGain = vecAngleMeanAndGain;
			sProjection.vecAngleMeanAndGainRand = vecAngleMeanAndGainRand;
			
			sProjection.vecSdProjMean = vecSdProjMean;
			sProjection.vecSdProjGain1 = vecSdProjGain1;
			sProjection.vecSdProjGain = vecSdProjGain;
			sProjection.vecSdProjRand = vecSdProjRand;
			sProjection.vecSdProjOrth = vecSdProjOrth;
			sProjection.vecSdProjAdja = vecSdProjAdja;
			
			sProjection.vecDistFromOrigin = vecDistFromOrigin;
			sProjection.vecDistFromOrthOri = vecDistFromOrthOri;
			sProjection.vecDistFromAdjaOri = vecDistFromAdjaOri;
			
			%% save
			save([strTargetDataPath 'QC2Data' strTypeRec '.mat'],...
				'sPrediction',...
				'sProjection',...
				'dblRM','dblRG',...
				'strRec','strArea','strType');
		end
	end
end
toc
