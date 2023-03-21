%% aim
%{
https://www.nature.com/articles/s41598-021-95037-z
https://www.nature.com/articles/s42003-021-02437-y

q1: can single-neuron variability be explained as population-level gain multiplication?
> estimate tuning curve from real data, then apply trial-by-trial gain to all neurons
>model predicting firing rate as combination of mean-rate multiplied by tuning curve => what do
residuals look like?
A: gain axis and mean-rate axis do not align, but ori tuning is distributed as conic manifold around
gain axis. Using stim-specific gain coupling calculated cross-validated for each neuron, allows pop
response to be predicted with R^2 of 0.72

%see articles:
https://elifesciences.org/articles/8998
https://www.nature.com/articles/nn.3711
https://www.jneurosci.org/content/39/37/7344.abstract
etc

%}
%% define qualifying areas
clear all;
strRunArea = 'Primary visual area';%'posteromedial visual area' 'Primary visual area'
cellUseAreas = {strRunArea};

boolMakeFigs = true;
boolSaveFigs = true;
boolSaveData = true;
boolHome = false;
if boolHome
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
if ~exist('sAggStim','var') || isempty(sAggStim)
	[sAggStim,sAggNeuron,sAggSources]=loadDataNpx('','driftinggrating',strDataPath);
end
sAggNeuron(strcmpi({sAggNeuron.SubjectType},'DBA')) = [];

%% pre-allocate matrices
intAreas = numel(cellUseAreas);
matDecPerf = [];
dblStartT = 0;

%% go through recordings
tic
for intRandomize=1:3
	for intRec=1:numel(sAggStim)
		close all;
		% get matching recording data
		close all;
		strRec = sAggStim(intRec).Exp;
		sThisRec = sAggStim(strcmpi(strRec,{sAggStim(:).Exp}));
		sThisSource = sAggSources(strcmpi(strRec,{sAggSources(:).Exp}));
		
		%prep grating data
		[sUseNeuron,vecStimOnTime,vecStimOffTime,vecOrientation,structStim] = NpxPrepGrating(sAggNeuron,sThisRec,cellUseAreas);
		[vecStimIdx,vecUniqueStims,vecRepNum,cellSelect,vecTrialRepetition] = val2idx(vecOrientation);
		intTrialNum = numel(vecStimOnTime);
		intStimNum = numel(unique(vecOrientation));
		intRepNum = intTrialNum/intStimNum;
		if numel(sUseNeuron) == 0, continue;end
		
		% concatenate stimulus structures
		%if contains(strRec,'20191213_MP3_RunDriftingGratingsR01')
		%	structStim = catstim(sThisRec.cellBlock(1));
		%	intMaxRep = 40;
		intMaxRep = inf;
		
		%% get orientation responses & single-trial population noise
		%get stim vars
		vecTempFreq = cell2vec({structStim.sStimObject(structStim.vecTrialStimTypes).TemporalFrequency})';
		vecPhase = structStim.Phase;
		vecDelayTimeBy = vecPhase./vecTempFreq;
		vecOri180 = mod(vecOrientation,180)*2;
		[vecStimIdx,vecUniqueStims,vecRepNum,cellSelect,vecTrialRepetition] = val2idx(vecOri180);
		indRem=vecTrialRepetition>intMaxRep;
		vecOrientation(indRem) = [];
		vecOri180(indRem) = [];
		vecDelayTimeBy(indRem) = [];
		vecStimOnTime(indRem) = [];
		vecStimOffTime(indRem) = [];
		[vecStimIdx,vecUniqueStims,vecRepNum,cellSelect,vecTrialRepetition] = val2idx(vecOri180);
		intTrialNum = numel(vecStimOnTime);
		intStimNum = numel(vecUniqueStims);
		intRepNum = intTrialNum/intStimNum;
		
		%pupil
		if isfield(structStim,'vecPupilStimOn')
			vecPupilStimOn = structStim.vecPupilStimOn;
			vecPupilStimOff = structStim.vecPupilStimOff;
			vecPupilStimOn(indRem) = [];
			vecPupilStimOff(indRem) = [];
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
		
		
		%change name
		if intRandomize == 1
			strType = 'Real';
		elseif intRandomize == 2
			strType = 'Shuff';
		elseif intRandomize == 3
			strType = 'Poiss';
		end
		strRec = [strType '_' strRec];
		
		%% select area 1
		for intArea=1:intAreas
			strArea = cellUseAreas{intArea};
			indArea1Neurons = contains({sUseNeuron.Area},strArea,'IgnoreCase',true);
			if sum(indArea1Neurons) == 0, continue;end
			
			%% get orientation responses & single-trial population noise
			sArea1Neurons = sUseNeuron(indArea1Neurons);
			cellSpikeTimesRaw = {sArea1Neurons.SpikeTimes};
			[matData,indTuned,cellSpikeTimes,sOut,cellSpikeTimesPerCellPerTrial,vecStimOnStitched,vecNonStat,dblBC,dblMaxDevFrac] = ...
				NpxPrepData(cellSpikeTimesRaw,vecStimOnTime,vecStimOffTime,vecOrientation);
			intTunedN = sum(indTuned);
			intRespN = size(matData,1);
			
			%% define quantiles and remove zero-variance neurons
			%constants
			[vecTrialTypeIdx,vecUnique,vecPriorDistribution,cellSelect,vecRepetition] = val2idx(vecOri180);
			intStimNr = numel(vecUnique);
			dblLambda = 1;%1
			intTypeCV = 2;
			dblUseStartT = 0;
			dblUseMaxDur = dblMaxDur-dblUseStartT;
			intUseMax = inf;
			intReps = mean(vecPriorDistribution);
			
			%remove zero-variance neurons
			matSpikeCounts_pre = cellfun(@(x) sum(x>dblUseStartT),cellSpikeTimesPerCellPerTrial);
			matMeanRate_pre = matSpikeCounts_pre./dblUseMaxDur;
			vecPopRate_pre = sum(matMeanRate_pre,1);
			
			intQuantiles = 5;
			vecStartTrials = round(intReps*linspace(1/intReps,(1+1/intReps),intQuantiles+1));
			vecStartTrials(end)=[];
			intSplitTrialsPerOri = min(floor(cellfun(@sum,cellSelect)/intQuantiles));
			indZeroVarNeurons = false(intRespN,1);
			for intQ=1:intQuantiles
				vecUseTrialsTemp = vecStartTrials(intQ):(vecStartTrials(intQ)+intSplitTrialsPerOri-1);
				for intOri=1:intStimNr
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
			if intRandomize > 1
				for intStim=1:intStimNr
					vecUseT = find(vecTrialTypeIdx==intStim);
					for intN=1:intRespN
						if intRandomize == 2
							%shuffle spikes
							matMeanRate(intN,vecUseT) = matMeanRate(intN,vecUseT(randperm(numel(vecUseT))));
						elseif intRandomize == 3
							%generate spikes
							dblMean = mean(cellfun(@(x) sum(x>dblUseStartT),cellSpikeTimesPerCellPerTrial(intN,vecUseT)));
							vecRates = poissrnd(dblMean,size(vecUseT));
							matMeanRate(intN,vecUseT) = vecRates;
						else
							error
						end
					end
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
				title(sprintf('%s',strRec),'interpreter','none');
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
					export_fig(fullpath(strFigurePathSR,sprintf('QC2_PopRespPrediction_%s.tif',strRec)));
					export_fig(fullpath(strFigurePathSR,sprintf('QC2_PopRespPrediction_%s.pdf',strRec)));
				end
			end
			
			
			%% project
			%for a given cloud of pop responses (n=trials), what fraction of noise aligns with the
			%mean-rate axis? what proportion with f'? how does this compare to a random direction (=spectral decomposition of covariance matrix)?
			
			%run projections
			fprintf('   Running multi-dim analysis [%s]\n',getTime);
			vecGainAx = mean(matTuningCurves,2);
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
					export_fig(fullpath(strFigurePathSR,sprintf('QC2_Projections_%s.tif',strRec)));
					export_fig(fullpath(strFigurePathSR,sprintf('QC2_Projections_%s.pdf',strRec)));
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
				vecPupilSize = mean(getRespMat(sPupil.vecTime,sPupil.vecRadius,vecPupilStimOn-0.5,[0 dblDur]),2);
				vecPupilMov = mean(getRespMat(sPupil.vecTime,vecPupMov,vecPupilStimOn-0.5,[0 dblDur]),2);
				vecTrialGain = getProjOnLine(matMeanRate,vecGainAx);
				vecTrialMean = getProjOnLine(matMeanRate,ones(size(vecGainAx)));
				vecTrialMean2 = mean(matMeanRate,1);
				rGainPupilS = nancorr(vecTrialGain(:),vecPupilSize(:));
				rMeanPupilS = nancorr(vecTrialMean(:),vecPupilSize(:));
				indUseTrials = ~isnan(vecPupilSize) & ~isnan(vecTrialMean) & ~isnan(vecTrialGain);
				
				%run predictions
				[vecTrialGain_PupilHat,dblR2_CV_G,matB_G] = doCrossValidatedDecodingRR(zscore(vecPupilSize(indUseTrials)),zscore(vecTrialGain(indUseTrials)),1,1);
				[vecTrialMean_PupilHat,dblR2_CV_M,matB_M] = doCrossValidatedDecodingRR(zscore(vecPupilSize(indUseTrials)),zscore(vecTrialMean(indUseTrials)),1,1);
				
				[dblR2G,dblSS_totG,dblSS_resG,dblTG,dblPG,dblR2_adjustedG,dblR2_SEG] = getR2(zscore(vecTrialGain(indUseTrials)),vecTrialGain_PupilHat,1);
				[dblR2M,dblSS_totM,dblSS_resM,dblTM,dblPM,dblR2_adjustedM,dblR2_SEM] = getR2(zscore(vecTrialMean(indUseTrials)),vecTrialMean_PupilHat,1);
				
				figure;maxfig;
				subplot(2,3,1);
				scatter(zscore(vecTrialGain(indUseTrials)),vecTrialGain_PupilHat,'b.');
				xlabel('Z-scored trial gain');
				ylabel('CV predicted norm trial gain');
				title(sprintf('%s; Gain prediction from pupil size',strRec));
				fixfig;
				
				subplot(2,3,2);
				scatter(zscore(vecTrialMean(indUseTrials)),vecTrialMean_PupilHat,'k.');
				xlabel('Z-scored trial mean');
				ylabel('CV predicted norm trial mean');
				title('Mean prediction from pupil size');
				fixfig;
				
				subplot(2,3,3);
				hold on
				errorbar(0.2,dblR2G,dblR2_SEG,'bx');
				errorbar(0.8,dblR2M,dblR2_SEM,'kx');
				set(gca,'xtick',[0.2 0.8],'xticklabel',{'Gain','Mean'});
				ylabel('CV predictability (R^2)');
				xlim([0 1]);
				title(sprintf('T-tests vs 0: gain-p=%.2e,mean-p=%.2e',dblPG,dblPM));
				fixfig;
				
				if boolSaveFigs
					%% save fig
					export_fig(fullpath(strFigurePathSR,sprintf('QC2_PupilPredGain_%s.tif',strRec)));
					export_fig(fullpath(strFigurePathSR,sprintf('QC2_PupilPredGain_%s.pdf',strRec)));
				end
			end
			
			%% use locomotion
			if ~isempty(vecRunningTime) && ~isempty(vecRunningSpeed)
				vecTime = vecRunningTime;
				vecVals = vecRunningSpeed;
				vecEvents = vecStimOnTime;
				vecWindow  =[0 dblDur];
				
				vecTrialSpeed = mean(getRespMat(vecTime,vecVals,vecEvents,vecWindow),2);
				
				vecTrialGain = getProjOnLine(matMeanRate,vecGainAx);
				vecTrialMean = getProjOnLine(matMeanRate,ones(size(vecGainAx)));
				vecTrialMean2 = mean(matMeanRate,1);
				rGainPupilS = nancorr(vecTrialGain(:),vecTrialSpeed(:));
				rMeanPupilS = nancorr(vecTrialMean(:),vecTrialSpeed(:));
				indUseTrials = ~isnan(vecTrialSpeed) & ~isnan(vecTrialMean) & ~isnan(vecTrialGain);
				
				%run predictions
				matX = vecTrialSpeed(indUseTrials);
				dblMaxMean = max(abs(mean(matX) ./ range(matX)));
				
				vecTrialGain_RunningHat = doCrossValidatedDecodingRR(zscore(vecTrialSpeed(indUseTrials)),zscore(vecTrialGain(indUseTrials)),1,1);
				vecTrialMean_RunningHat = doCrossValidatedDecodingRR(zscore(vecTrialSpeed(indUseTrials)),zscore(vecTrialMean(indUseTrials)),1,1);
				
				[dblR2G,dblSS_totG,dblSS_resG,dblTG,dblPG,dblR2_adjustedG,dblR2_SEG] = getR2(zscore(vecTrialGain(indUseTrials)),vecTrialGain_RunningHat,1);
				[dblR2M,dblSS_totM,dblSS_resM,dblTM,dblPM,dblR2_adjustedM,dblR2_SEM] = getR2(zscore(vecTrialMean(indUseTrials)),vecTrialMean_RunningHat,1);
				
				figure;maxfig;
				subplot(2,3,1);
				scatter(zscore(vecTrialGain(indUseTrials)),vecTrialGain_RunningHat,'b.');
				xlabel('Z-scored trial gain');
				ylabel('CV predicted norm trial gain');
				title('Gain prediction from pupil size');
				fixfig;
				
				subplot(2,3,2);
				scatter(zscore(vecTrialMean(indUseTrials)),vecTrialMean_RunningHat,'k.');
				xlabel('Z-scored trial mean');
				ylabel('CV predicted norm trial mean');
				title('Mean prediction from pupil size');
				fixfig;
				
				subplot(2,3,3);
				hold on
				errorbar(0.2,dblR2G,dblR2_SEG,'bx');
				errorbar(0.8,dblR2M,dblR2_SEM,'kx');
				set(gca,'xtick',[0.2 0.8],'xticklabel',{'Gain','Mean'});
				ylabel('CV predictability (R^2)');
				xlim([0 1]);
				title(sprintf('T-tests vs 0: gain-p=%.2e,mean-p=%.2e',dblPG,dblPM));
				fixfig;
				
				if boolSaveFigs
					%% save fig
					export_fig(fullpath(strFigurePathSR,sprintf('QC2_RunningPredGain_%s.tif',strRec)));
					export_fig(fullpath(strFigurePathSR,sprintf('QC2_RunningPredGain_%s.pdf',strRec)));
				end
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
				sPrediction.vecTrialMean2 = vecTrialMean2;
				
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
				save([strTargetDataPath 'QC2Data' strRec '.mat'],...
					'sPrediction',...
					'sProjection',...
					'strRec','strArea');
			end
		end
	end
end
toc
