%% aim
%{
https://www.nature.com/articles/s41598-021-95037-z
https://www.nature.com/articles/s42003-021-02437-y

q1: Progression of orientation information over time: how does trial activity evolve? what is the function of the onset peak?
> is decoding better when matched for stimulus phase? => no.

q2: "Spike time and rate coding can be represented within a single model of spiking probability as a
function of time: rate codes are uniform over a certain period tau, while spike time codes are
temporally localized peaks"
> is this true?

q3: Rate codes do not exist; a rate code is simply a subset of spike time codes where the temporal
integration window is very large. But what about multi dim codes? Those are all rate based. Can we
formulate a multidimensional spike-time code? I.e., can we make a taxonomy of neural codes?

q4: How does information evolve over time, is initial peak indeed less tuned? Is pop activity
rhythmic? Are stimuli encoded invariant to brain state? Eg, high arousal, low arousal. Or is
stimulus manifold dynamic over time? Does manifold scale with arousal? => How does manifold depend
on binning size? What is the optimal time window?

%}
%% define qualifying areas
clear all;
cellUseAreas = {...
	'Primary visual area',...
	...'posteromedial visual area',...
	};
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
	[sAggStim,sAggNeuron]=loadDataNpx('','driftinggrating',strDataPath);
end
sAggNeuron(strcmpi({sAggNeuron.SubjectType},'DBA')) = [];

%% pre-allocate matrices
intAreas = numel(cellUseAreas);
matR_PP_All = nan(intAreas,intAreas,numel(sAggStim),2);
matR_OP_All = nan(intAreas,intAreas,numel(sAggStim),2);
matR_OO_All = nan(intAreas,intAreas,numel(sAggStim),2);
matR_PO_All = nan(intAreas,intAreas,numel(sAggStim),2);
mat_xR_All = nan(intAreas,intAreas,numel(sAggStim),2);

%% go through recordings
tic
for intRec=11%1:numel(sAggStim) %19 || weird: 11
	% get matching recording data
	strRec = sAggStim(intRec).Exp;
	sThisRec = sAggStim(strcmpi(strRec,{sAggStim(:).Exp}));
	
	%get AP data
	sLoad = load(sThisRec.File);
	sAP = sLoad.sAP;
	
	%remove stimulus sets that are not 24 stim types
	sThisRec.cellBlock(cellfun(@(x) x.intTrialNum/x.intNumRepeats,sThisRec.cellBlock) ~= 24) = [];
	% concatenate stimulus structures
	intBlock = 1;
	structStim = catstim(sThisRec.cellBlock(intBlock));%(1:end));
	vecStimOnTime = structStim.vecStimOnTime;
	vecStimOffTime = structStim.vecStimOffTime;
	vecStimOnTime1 = structStim.vecStimOnTime;
	vecStimOffTime1 = structStim.vecStimOffTime;
	%vecStimOnTime = structStim.ActOnNI;
	%vecStimOffTime = structStim.ActOffNI;
	vecOrientation = cell2vec({structStim.sStimObject(structStim.vecTrialStimTypes).Orientation})';
	vecTempFreq = cell2vec({structStim.sStimObject(structStim.vecTrialStimTypes).TemporalFrequency})';
	vecPhase = structStim.Phase;
	vecDelayTimeBy = vecPhase./vecTempFreq;
	[vecOriIdx,vecUniqueOris,vecRepNum,cellSelect,vecTrialRepetition] = val2idx(vecOrientation);
	indRem=vecTrialRepetition>min(vecRepNum);
	vecOrientation(indRem) = [];
	vecDelayTimeBy(indRem) = [];
	vecStimOnTime(indRem) = [];
	vecStimOffTime(indRem) = [];
	[vecOriIdx,vecUniqueOris,vecRepNum,cellSelect,vecTrialRepetition] = val2idx(vecOrientation);
	intTrialNum = numel(vecStimOnTime);
	intStimNum = numel(unique(vecOrientation));
	intRepNum = intTrialNum/intStimNum;
	
	%remove neurons from other recs
	indQualifyingNeurons = contains({sAggNeuron.Exp},strRec);
	
	%remove neurons in incorrect areas
	indConsiderNeurons = contains({sAggNeuron.Area},cellUseAreas,'IgnoreCase',true);
	
	%remove bad clusters
	indGoodNeurons = (cell2vec({sAggNeuron.KilosortGood}) == 1) | (cell2vec({sAggNeuron.Contamination}) < 0.1);
	
	%subselect from total
	indUseNeurons = indQualifyingNeurons(:) & indConsiderNeurons(:) & indGoodNeurons(:);
	sUseNeuron = sAggNeuron(indUseNeurons);
	
	% select area 1
	for intArea=1:intAreas
		strArea = cellUseAreas{intArea};
		indArea1Neurons = contains({sUseNeuron.Area},strArea,'IgnoreCase',true);
		if sum(indArea1Neurons) == 0, continue;end
		%% get orientation responses & single-trial population noise
		sArea1Neurons = sUseNeuron(indArea1Neurons);
		intRec
		%% remove untuned cells
		%get data matrix
		cellSpikeTimes = {sArea1Neurons.SpikeTimes};
		dblDur = median(vecStimOffTime-vecStimOnTime);
		matData = getSpikeCounts(cellSpikeTimes,vecStimOnTime,dblDur)./dblDur;
		
		%remove untuned cells
		vecOri180 = mod(vecOrientation,180)*2;
		sOut = getTuningCurves(matData,vecOri180,0);
		dblMinRate = 0.1;
		indTuned = sOut.vecOriAnova<0.05;
		indResp = cellfun(@min,{sArea1Neurons.ZetaP}) < 0.05 & sum(matData,2)'>(size(matData,2)/dblDur)*dblMinRate;
		
		%prep
		vecPrefOri = rad2deg(sOut.matFittedParams(indResp,1))/2;
		vecPrefRad = sOut.matFittedParams(indResp,1);
		cellSpikeTimes(~indResp)=[];
		indTuned(~indResp)=[];
		
		intTunedN = sum(indTuned);
		intNumN = sum(indResp);
		matResp = matData(indResp,:);
		
		vecOri180 = mod(vecOrientation,180)*2;
		[vecOriIdx,vecUnique,vecPriorDistribution,cellSelect,vecRepetition] = val2idx(vecOri180);
		intOriNum = numel(vecUnique);
		intRepNum = min(vecPriorDistribution);
		dblStimDur = roundi(median(vecStimOffTime - vecStimOnTime),1,'ceil');
		dblBinWidth = 5/1000;%/32
		dblPreTime = 10*dblBinWidth;
		dblPostTime = 30*dblBinWidth;
		dblMaxDur = dblStimDur+dblPreTime+dblPostTime;
		%dblMaxDur = dblPreTime+dblPostTime;
		vecBinEdges = 0:dblBinWidth:dblMaxDur;
		vecStimTime = vecBinEdges(2:end)-dblBinWidth/2 - dblPreTime;
		indStimBins = vecStimTime > 0 & vecStimTime < dblStimDur;
		intBinNum = numel(vecBinEdges)-1;
		matBNSR = nan(intBinNum,intNumN,intOriNum,intRepNum);
		matBNT = nan(intBinNum,intNumN,intTrialNum);
		%matBNT_shifted = nan(intBinNum,intNumN,intTrialNum);
		
		for intN=1:intNumN
			vecRepCounter = zeros(1,intOriNum);
			[vecTrialPerSpike,vecTimePerSpike] = getSpikesInTrial(cellSpikeTimes{intN},vecStimOnTime-dblPreTime,dblMaxDur);
			for intTrial=1:intTrialNum
				intTrialOriIdx = vecOriIdx(intTrial);
				vecRepCounter(intTrialOriIdx) = vecRepCounter(intTrialOriIdx) + 1;
				intRep = vecRepCounter(intTrialOriIdx);
				vecSpikeT = vecTimePerSpike(vecTrialPerSpike==intTrial);
				vecSpikeC = histcounts(vecSpikeT,vecBinEdges);
				if any(isnan(vecSpikeC))
					error
				end
				matBNSR(:,intN,intTrialOriIdx,intRep) = vecSpikeC;
				matBNT(:,intN,intTrial) = vecSpikeC;
				%matBNT_shifted(:,intN,intTrial) = vecSpikeC;
				%matBNT_shifted(indStimBins,intN,intTrial) = circshift(matBNT_shifted(indStimBins,intN,intTrial),-round(intBinNum*vecDelayTimeBy(intTrial)));
			end
		end
		
		%% time progression
		vecSpikesPerBin = mean(sum(matBNT,2),3)';
		figure
		subplot(2,3,1)
		plot(vecStimTime,vecSpikesPerBin./dblBinWidth)
		title(sprintf('%s - block %d',strRec,intBlock));
		xlabel('Time after onset (s)');
		ylabel('Binned population spiking rate (Hz)');
		fixfig;
	end
end
toc
