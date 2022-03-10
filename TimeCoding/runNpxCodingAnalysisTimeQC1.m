%% aim
%{
https://www.nature.com/articles/s41598-021-95037-z
https://www.nature.com/articles/s42003-021-02437-y

q1: is coding better in trials with higher or lower firing rates?

q2: can we define peaks in the IFR as population events, and find which cells spike in the beginning
or end? does this ordering differ between orientations?

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
	strFigurePath = 'F:\Data\Results\PopTimeCoding';
else
	strDataPath = 'E:\DataPreProcessed\';
	strFigurePath = 'D:\Data\Results\PopTimeCoding';
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
boolSaveFigs = true;
%% go through recordings
tic
for intRec=19%1:numel(sAggStim)
	% get matching recording data
	strRec = sAggStim(intRec).Exp;
	sThisRec = sAggStim(strcmpi(strRec,{sAggStim(:).Exp}));
	sThisSource = sAggSources(strcmpi(strRec,{sAggSources(:).Exp}));
	
	%remove stimulus sets that are not 24 stim types
	sThisRec.cellBlock(cellfun(@(x) x.intTrialNum/x.intNumRepeats,sThisRec.cellBlock) ~= 24) = [];
	% concatenate stimulus structures
	structStim = catstim(sThisRec.cellBlock(1:2));
	vecStimOnTime = structStim.vecStimOnTime;
	vecStimOffTime = structStim.vecStimOffTime;
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
	intOriNum = numel(unique(vecOrientation));
	intRepNum = intTrialNum/intOriNum;
	
	%remove neurons from other recs
	indQualifyingNeurons = contains({sAggNeuron.Exp},strRec);
	
	%remove neurons in incorrect areas
	indConsiderNeurons = contains({sAggNeuron.Area},cellUseAreas,'IgnoreCase',true);
	
	%remove bad clusters
	indGoodNeurons = (cell2vec({sAggNeuron.KilosortGood}) == 1) | (cell2vec({sAggNeuron.Contamination}) < 0.1);
	
	%subselect from total
	indUseNeurons = indQualifyingNeurons(:) & indConsiderNeurons(:) & indGoodNeurons(:);
	sUseNeuron = sAggNeuron(indUseNeurons);
	
	%% select area 1
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
		matData = getSpikeCounts(cellSpikeTimes,vecStimOnTime,dblDur);
		
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
		
		dblStimDur = roundi(median(vecStimOffTime - vecStimOnTime),1,'ceil');
		dblPreTime = -dblStartT;%0.3;
		dblPostTime = 0;%0.3;
		dblMaxDur = dblStimDur+dblPreTime+dblPostTime;
		dblBinWidth = 0.05;
		vecBinEdges = 0:dblBinWidth:dblMaxDur;
		vecStimTime = vecBinEdges(2:end)-dblBinWidth/2 - dblPreTime;
		indStimBins = vecStimTime > 0 & vecStimTime < dblStimDur;
		intBinNum = numel(vecBinEdges)-1;
		matBNSR = nan(intBinNum,intNumN,intOriNum,intRepNum);
		matBNT = nan(intBinNum,intNumN,intTrialNum);
		%matBNT_shifted = nan(intBinNum,intNumN,intTrialNum);
		
		vecRepCounter = zeros(1,intOriNum);
		%get spikes per trial per neuron
		cellSpikeTimesStitched = cell(intNumN,intTrialNum);
		cellSpikeTimesPerCellPerTrial = cell(intNumN,intTrialNum);
		cellSpikeTimesPerCellPerTrial_S = cell(intNumN,intTrialNum); %single-neuron ISIs, shuffled per trial
		for intN=1:intNumN
			% build pseudo data, stitching stimulus periods
			[vecPseudoSpikeTimes,vecPseudoEventT] = getPseudoSpikeVectors(cellSpikeTimes{intN},vecStimOnTime-dblPreTime,dblMaxDur);
			cellSpikeTimesStitched{intN} = vecPseudoSpikeTimes;
			
			%real
			[vecTrialPerSpike,vecTimePerSpike] = getSpikesInTrial(vecPseudoSpikeTimes,vecPseudoEventT,dblMaxDur);
			for intTrial=1:intTrialNum
				vecSpikeT = vecTimePerSpike(vecTrialPerSpike==intTrial);
				vecISI = diff(vecSpikeT);
				if isempty(vecSpikeT)
					vecGenSpikesS = [];
				else
					vecISIS = diff(vecSpikeT);
					vecGenSpikesS = cumsum([vecSpikeT(1);vecISI(randperm(numel(vecISI)))]);
				end
				
				cellSpikeTimesPerCellPerTrial{intN,intTrial} = vecSpikeT;
				cellSpikeTimesPerCellPerTrial_S{intN,intTrial} = vecGenSpikesS;
			end
		end
		vecStimOnStitched = vecPseudoEventT;
		
		%% decode using first spike delay
		%constants
		[vecTrialTypeIdx,vecUnique,vecPriorDistribution,cellSelect,vecRepetition] = val2idx(vecOri180);
		intStimNr = numel(vecUnique);
		dblLambda = 1;%1
		intTypeCV = 2;
		dblUseStartT = 0;
		dblUseMaxDur = dblMaxDur-dblUseStartT;
		intUseMax = inf;
		intReps = mean(vecPriorDistribution);
		
		%simple "rate code"
		intQuantiles = 5;
		matSpikeCounts = cellfun(@(x) sum(x>dblUseStartT),cellSpikeTimesPerCellPerTrial);
		matMeanRate = matSpikeCounts./dblUseMaxDur;
		vecPopRate = sum(matMeanRate,1);
		intSplitTrialsPerOri = min(floor(cellfun(@sum,cellSelect)/intQuantiles));
		vecPriorDistributionSplit = ones(size(vecPriorDistribution))*intSplitTrialsPerOri;
		vecStartTrials = round(intReps*linspace(1/intReps,(1+1/intReps),intQuantiles+1));
		vecStartTrials(end)=[];
		
		%% plot population mean + sd over neurons, compare with neuron mean+sd over trials 
		intTrials = size(matMeanRate,2);
		matGenRate = poissrnd(repmat(vecNeuronMean,[1 intTrials]));
		
		vecPopMean = mean(matMeanRate,1);
		vecPopSd = std(matMeanRate,[],1);
		vecPopCv = vecPopSd./vecPopMean;
		vecPopFano = (vecPopSd.^2)./vecPopMean;
		
		vecNeuronMean = mean(matMeanRate,2);
		vecNeuronSd = std(matMeanRate,[],2);
		vecNeuronCv = vecNeuronSd./vecNeuronMean;
		vecNeuronFano = (vecNeuronSd.^2)./vecNeuronMean;
		
		vecGenPopMean = mean(matGenRate,1);
		vecGenPopSd = std(matGenRate,[],1);
		vecGenPopCv = vecGenPopSd./vecGenPopMean;
		vecGenPopFano = (vecGenPopSd.^2)./vecGenPopMean;
		
		vecGenNeuronMean = mean(matGenRate,2);
		vecGenNeuronSd = std(matGenRate,[],2);
		vecGenNeuronCv = vecGenNeuronSd./vecGenNeuronMean;
		vecGenNeuronFano = (vecGenNeuronSd.^2)./vecGenNeuronMean;
		
		figure;maxfig
		subplot(3,4,1)
		scatter(vecPopMean,vecPopSd,'.')
		xlabel('Population mean rate (Hz)');
		ylabel('Pop sd (Hz)');fixfig;
		
		subplot(3,4,5)
		scatter(vecPopMean,vecPopCv,'.')
		xlabel('Population mean rate (Hz)');
		ylabel('Pop cv');fixfig;
		
		subplot(3,4,9)
		scatter(vecPopMean,vecPopFano,'.')
		xlabel('Population mean rate (Hz)');
		ylabel('Pop Fano factor');fixfig;
		
		subplot(3,4,2)
		scatter(vecNeuronMean,vecNeuronSd,'.')
		xlabel('Neuron mean rate (Hz)');
		ylabel('Neuron sd (Hz)');fixfig;
		
		subplot(3,4,6)
		scatter(vecNeuronMean,vecNeuronCv,'.')
		xlabel('Neuron mean rate (Hz)');
		ylabel('Neuron cv');fixfig;
		
		subplot(3,4,10)
		scatter(vecNeuronMean,vecNeuronFano,'.')
		xlabel('Neuron mean rate (Hz)');
		ylabel('Neuron Fano factor');fixfig;
		
		subplot(3,4,3)
		scatter(vecGenPopMean,vecGenPopSd,'.')
		xlabel('Gen Population mean rate (Hz)');
		ylabel('Gen Pop sd (Hz)');fixfig;
		
		subplot(3,4,7)
		scatter(vecGenPopMean,vecGenPopCv,'.')
		xlabel('Gen Population mean rate (Hz)');
		ylabel('Gen Pop cv');fixfig;
		
		subplot(3,4,11)
		scatter(vecGenPopMean,vecGenPopFano,'.')
		xlabel('Gen Population mean rate (Hz)');
		ylabel('Gen Pop Fano factor');fixfig;
		
		subplot(3,4,4)
		scatter(vecGenNeuronMean,vecGenNeuronSd,'.')
		xlabel('Gen Neuron mean rate (Hz)');
		ylabel('Gen Neuron sd (Hz)');fixfig;
		
		subplot(3,4,8)
		scatter(vecGenNeuronMean,vecGenNeuronCv,'.')
		xlabel('Gen Neuron mean rate (Hz)');
		ylabel('Gen Neuron cv');fixfig;
		
		subplot(3,4,12)
		scatter(vecGenNeuronMean,vecGenNeuronFano,'.')
		xlabel('Gen Neuron mean rate (Hz)');
		ylabel('Gen Neuron Fano factor');fixfig;
		
		return
		%% pre-allocate
		vecQuantilePerf = nan(1,intQuantiles);
		matQuantileConfusion = nan(intStimNr,intStimNr,intQuantiles);
		
		%split trials into quantiles
		figure;maxfig;
		for intQ=1:intQuantiles
			matUseTrials = nan(intStimNr,intSplitTrialsPerOri);
			vecUseTrials = vecStartTrials(intQ):(vecStartTrials(intQ)+intSplitTrialsPerOri-1);
			for intOri=1:intStimNr
				vecThisOri = find(cellSelect{intOri});
				[vecSorted,vecReorder]=sort(vecPopRate(vecThisOri));
				matUseTrials(intOri,:) = vecThisOri(vecReorder(vecUseTrials));
			end
			
			vecUseTrials = sort(matUseTrials(:));
			matUseRate = matMeanRate(:,vecUseTrials);
			vecOriUse = vecOri180(vecUseTrials);
			
			[dblPerfQ,vecDecodedIndexRateCV,matPosteriorProbabilityRate,dblMeanErrorDegsRate,matConfusion,matWeightsRate] = ...
				doCrossValidatedDecodingLR(matUseRate,vecOriUse,intTypeCV,vecPriorDistributionSplit,dblLambda);
			
			vecQuantilePerf(intQ) = dblPerfQ;
			matQuantileConfusion(:,:,intQ) = matConfusion;
			
			subplot(2,3,intQ)
			imagesc(matConfusion)
			xlabel('Real Orientation')
			ylabel('Decoded Orientation')
			title(sprintf('Quantile %d',intQ))
			fixfig;grid off;
		end
		subplot(2,3,6)
		hold on
		plot([1 numel(vecQuantilePerf)],(1/intStimNr)*[1 1],'--','color',[0.5 0.5 0.5]);
		plot(1:numel(vecQuantilePerf),vecQuantilePerf,'color',lines(1));
		xlabel('Quantile');
		ylabel('CV decoding accuracy');
		fixfig;
		
	end
end
toc
