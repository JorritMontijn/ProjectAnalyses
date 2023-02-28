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
boolHome = true;
if boolHome
	strDataPath = 'F:\Data\Processed\Neuropixels\';
	strFigurePathSR = 'F:\Drive\PopTimeCoding\single_recs';
	strFigurePath = 'F:\Drive\PopTimeCoding\figures\';
	strTargetDataPath = 'F:\Drive\PopTimeCoding\data\';
else
	strDataPath = 'E:\DataPreProcessed\';
	strFigurePath = 'D:\Data\Results\PopTimeCoding';
end


%% select all neurons in LP and drifting grating stimuli
if ~exist('sAggStim','var') || isempty(sAggStim)
	[sAggStim,sAggNeuron]=loadDataNpx('','driftinggrating',strDataPath);
end
sAggNeuron(strcmpi({sAggNeuron.SubjectType},'DBA')) = [];

%% pre-allocate matrices
intAreas = numel(cellUseAreas);
matDecPerf = [];
dblStartT = 0.1;

%% go through recordings
tic
for intRec=1:numel(sAggStim)
	close all;
	% get matching recording data
	strRec = sAggStim(intRec).Exp;
	sThisRec = sAggStim(strcmpi(strRec,{sAggStim(:).Exp}));
	
	%remove stimulus sets that are not 24 stim types
	sThisRec.cellBlock(cellfun(@(x) x.intTrialNum/x.intNumRepeats,sThisRec.cellBlock) ~= 24) = [];
	% concatenate stimulus structures
	structStim = catstim(sThisRec.cellBlock(1:end));
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
		%indTuned = sOut.vecOriAnova<0.05;
		dblMinRate = 0.1;
		indResp = cellfun(@min,{sArea1Neurons.ZetaP}) < 0.05 & sum(matData,2)'>(size(matData,2)/dblDur)*dblMinRate;
		
		%prep
		vecPrefOri = rad2deg(sOut.matFittedParams(indResp,1))/2;
		vecPrefRad = sOut.matFittedParams(indResp,1);
		intTunedN = sum(indResp);
		intNumN = intTunedN;
		cellSpikeTimes(~indResp)=[];
		
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
		cellSpikeTimesPerCellPerTrial_SN = cell(intNumN,intTrialNum); %shuffled ISIs per neuron over all trials
		cellSpikeTimesPerCellPerTrial_SS = cell(intNumN,intTrialNum); %shuffled ISIs per neuron over all trials,single-neuron ISIs, shuffled per trial
		for intN=1:intNumN
			% build pseudo data, stitching stimulus periods
			[vecPseudoSpikeTimes,vecPseudoEventT] = getPseudoSpikeVectors(cellSpikeTimes{intN},vecStimOnTime-dblPreTime,dblMaxDur);
			cellSpikeTimesStitched{intN} = vecPseudoSpikeTimes;
			vecISI_Overall = diff(vecPseudoSpikeTimes);
			vecOverallSpikeT_S = cumsum([vecPseudoSpikeTimes(1);vecISI_Overall(randperm(numel(vecISI_Overall)))]);
			
			%real
			[vecTrialPerSpike,vecTimePerSpike] = getSpikesInTrial(vecPseudoSpikeTimes,vecPseudoEventT,dblMaxDur);
			%shuffled
			[vecTrialPerSpikeS,vecTimePerSpikeS] = getSpikesInTrial(vecOverallSpikeT_S,vecPseudoEventT,dblMaxDur);
			for intTrial=1:intTrialNum
				%real
				vecSpikeT = sort(vecTimePerSpike(vecTrialPerSpike==intTrial));
				vecISI = diff(vecSpikeT);
				if isempty(vecSpikeT)
					vecGenSpikesS = [];
				else
					vecISIS = diff(vecSpikeT);
					vecGenSpikesS = cumsum([vecSpikeT(1);vecISI(randperm(numel(vecISI)))]);
				end
				cellSpikeTimesPerCellPerTrial{intN,intTrial} = vecSpikeT;
				cellSpikeTimesPerCellPerTrial_S{intN,intTrial} = vecGenSpikesS;
				
				%overall shuffle
				vecSpikeT_S = sort(vecTimePerSpikeS(vecTrialPerSpikeS==intTrial));
				vecISI_S = diff(vecSpikeT_S);
				if isempty(vecSpikeT_S)
					vecGenSpikesSS = [];
				else
					vecISIS = diff(vecSpikeT_S);
					vecGenSpikesSS = cumsum([vecSpikeT_S(1);vecISI_S(randperm(numel(vecISI_S)))]);
				end
				cellSpikeTimesPerCellPerTrial_SN{intN,intTrial} = vecSpikeT_S;
				cellSpikeTimesPerCellPerTrial_SS{intN,intTrial} = vecGenSpikesSS;
			end
		end
		vecStimOnStitched = vecPseudoEventT;
		
		%% calc IFRs per trial
		cellIFR_perTrial = cell(intTrialNum,1);
		cellTimeIFR_perTrial = cell(intTrialNum,1);
		cellIFR_perTrial_S = cell(intTrialNum,1);
		cellTimeIFR_perTrial_S = cell(intTrialNum,1);
		cellIFR_perTrial_SN = cell(intTrialNum,1);
		cellTimeIFR_perTrial_SN = cell(intTrialNum,1);
		cellIFR_perTrial_SS = cell(intTrialNum,1);
		cellTimeIFR_perTrial_SS = cell(intTrialNum,1);
		for intTrial=1:intTrialNum
			%real
			vecAllSpikes = sort(cell2vec(cellSpikeTimesPerCellPerTrial(:,intTrial)));
			vecISI0 = [vecAllSpikes(2:end) - vecAllSpikes(1:(end-1)); inf];
			vecAllSpikes(vecISI0==0)=vecAllSpikes(vecISI0==0)-(10^-5)*rand();
			vecAllSpikes = uniquetol(vecAllSpikes,1e-7);
			[vecTimeIFR,vecIFR] = getIFR(vecAllSpikes,0,dblMaxDur,[],[],[],0);
			cellIFR_perTrial{intTrial} = vecIFR;
			cellTimeIFR_perTrial{intTrial} = vecTimeIFR;
			
			%shuffled single-trial, single-neuron ISIs
			vecAllSpikesShuff = sort(cell2vec(cellSpikeTimesPerCellPerTrial_S(:,intTrial)));
			vecISI0 = [vecAllSpikesShuff(2:end) - vecAllSpikesShuff(1:(end-1)); inf];
			vecAllSpikesShuff(vecISI0==0)=vecAllSpikesShuff(vecISI0==0)-(10^-5)*rand();
			vecAllSpikesShuff = uniquetol(vecAllSpikesShuff,1e-7);
			[vecTimeIFRS,vecIFRS] = getIFR(vecAllSpikesShuff,0,dblMaxDur,[],[],[],0);
			cellIFR_perTrial_S{intTrial} = vecIFRS;
			cellTimeIFR_perTrial_S{intTrial} = vecTimeIFRS;
		
			%shuffled overall ISIs
			vecAllSpikes = sort(cell2vec(cellSpikeTimesPerCellPerTrial_SN(:,intTrial)));
			vecISI0 = [vecAllSpikes(2:end) - vecAllSpikes(1:(end-1)); inf];
			vecAllSpikes(vecISI0==0)=vecAllSpikes(vecISI0==0)-(10^-5)*rand();
			vecAllSpikes = uniquetol(vecAllSpikes,1e-7);
			[vecTimeIFR,vecIFR] = getIFR(vecAllSpikes,0,dblMaxDur,[],[],[],0);
			cellIFR_perTrial_SN{intTrial} = vecIFR;
			cellTimeIFR_perTrial_SN{intTrial} = vecTimeIFR;
			
			%shuffled overall ISIs, single-neuron ISIs
			vecAllSpikesShuff = sort(cell2vec(cellSpikeTimesPerCellPerTrial_SS(:,intTrial)));
			vecISI0 = [vecAllSpikesShuff(2:end) - vecAllSpikesShuff(1:(end-1)); inf];
			vecAllSpikesShuff(vecISI0==0)=vecAllSpikesShuff(vecISI0==0)-(10^-5)*rand();
			vecAllSpikesShuff = uniquetol(vecAllSpikesShuff,1e-7);
			[vecTimeIFRS,vecIFRS] = getIFR(vecAllSpikesShuff,0,dblMaxDur,[],[],[],0);
			cellIFR_perTrial_SS{intTrial} = vecIFRS;
			cellTimeIFR_perTrial_SS{intTrial} = vecTimeIFRS;
		end
		
		%% calculate activity during low/high epochs
		%are only cells tuned to the orientation active during low phases? Is this different from high 5%?
		
		%pre-alloc
		intQuantileNum = 5; %5=20%
		
		%run
		for intShuff=[0 1]
			%pre-allocate
			matAggR_temp = nan(3,intTrialNum,intTunedN);%save middle60/lower20/upper20
			for intTrial=1:intTrialNum
				%% define data
				if intShuff == 0
					%get IFRs
					vecTrialIFR = cellIFR_perTrial{intTrial}(2:(end-1));
					vecTrialTimeIFR = cellTimeIFR_perTrial{intTrial}(2:(end-1));
					cellSpikes = cellSpikeTimesPerCellPerTrial(:,intTrial);
				else
					%get IFRs
					vecTrialIFR = cellIFR_perTrial_S{intTrial}(2:(end-1));
					vecTrialTimeIFR = cellTimeIFR_perTrial_S{intTrial}(2:(end-1));
					cellSpikes = cellSpikeTimesPerCellPerTrial_S(:,intTrial);
				end
				
				
				%% divide quantiles
				[vecIFR_sorted,vecReorder] = sort(vecTrialIFR);
				vecTimeIFR_sorted = vecTrialTimeIFR(vecReorder);
				intSamples = numel(vecIFR_sorted);
				intEndLow = round(intSamples/intQuantileNum);
				dblLow_UpperBound = vecIFR_sorted(intEndLow);
				
				intStartHigh = intSamples-intEndLow+1;
				dblHigh_LowerBound = vecIFR_sorted(intStartHigh);
				vecHighLowIdx = ones(1,intSamples);
				vecHighLowIdx(vecTrialIFR<dblLow_UpperBound) = 2; %low
				vecHighLowIdx(vecTrialIFR>dblHigh_LowerBound) = 3; %high
				
				%plot
				if 0
					%% plot
					plot(vecTrialTimeIFR,vecTrialIFR,'k');
					hold on
					scatter(vecTrialTimeIFR(vecHighLowIdx==2),vecTrialIFR(vecHighLowIdx==2),'b.');
					scatter(vecTrialTimeIFR(vecHighLowIdx==3),vecTrialIFR(vecHighLowIdx==3),'r.');
					hold off
				end
				
				%% determine epochs
				vecChanges = find(diff(vecHighLowIdx)~=0);
				intEpochs = numel(vecChanges)+1;
				vecEpochType = nan(1,intEpochs);
				vecEpochStarts = nan(1,intEpochs);
				vecEpochStops = nan(1,intEpochs);
				dblLastEpoch = 0;
				for intEpoch=1:(intEpochs-1)
					intChangeIdx = vecChanges(intEpoch);
					dblNewEpoch = (vecTrialTimeIFR(intChangeIdx) + vecTrialTimeIFR(intChangeIdx+1))/2;
					vecEpochType(intEpoch) = vecHighLowIdx(intChangeIdx);
					vecEpochStarts(intEpoch) = dblLastEpoch;
					vecEpochStops(intEpoch) = dblNewEpoch;
					dblLastEpoch = dblNewEpoch;
				end
				vecEpochType(intEpochs) = vecHighLowIdx(end);
				vecEpochStarts(intEpochs) = dblNewEpoch;
				vecEpochStops(intEpochs) = dblMaxDur;
				
				%% assign epochs
				for intNeuron=1:intTunedN
					vecSpikes = cellSpikes{intNeuron};
					
					%do stuff here
					vecSpikeQ = zeros(size(vecSpikes));
					vecCountsPerType = zeros(1,3);
					for intSpikeIdx=1:numel(vecSpikes)
						intType = vecEpochType(vecEpochStarts < vecSpikes(intSpikeIdx) & vecEpochStops > vecSpikes(intSpikeIdx));
						vecCountsPerType(intType) = vecCountsPerType(intType) + 1;
					end
					
					%save
					matAggR_temp(:,intTrial,intNeuron) = vecCountsPerType; %low
				end
			end
			
			%% save
			if intShuff == 0
				matAggR = matAggR_temp;
			else
				matAggR_shuff = matAggR_temp;
			end
		end
		
		% mean per trial per quantile
		matLowR = squeeze(matAggR(2,:,:));
		matHighR = squeeze(matAggR(3,:,:));
		vecOriLow = vecOri180;
		vecOriHigh = vecOri180;
		
		%% so why is there no effect on the decoder?
		%=> compare greedy decoder of high vs low; greedy in neurons or spikes?
		intUseNeurons = 20;
		intUseNumDecoders = 1;
		for intDecoder=1:intUseNumDecoders
			dblLambda = 1;
			for intSwitchQ=0:1
				if intSwitchQ == 0
					matUseR = matLowR;
					vecUseOri = vecOriLow;
					vecPrior = vecPriorDistributionLow;
					strHighLow = 'low';
				else
					matUseR = matHighR;
					vecUseOri = vecOriHigh;
					vecPrior = vecPriorDistributionHigh;
					strHighLow = 'high';
				end
				intTotN = size(matUseR,2);
				indNeuronsUsed = false(1,intTotN);
				vecGreedyPerf = nan(1,intUseNeurons);
				vecSelectNeuronOrder = nan(1,intUseNeurons);
				for intNeuronsUsed=1:intUseNeurons
					intUnusedNeurons = intTotN-sum(indNeuronsUsed);
					vecUnusedNeurons = find(~indNeuronsUsed);
					vecUsedNeurons = find(indNeuronsUsed);
					vecPerf = nan(1,intUnusedNeurons);
					for intIdxN=1:intUnusedNeurons
						intNeuron = vecUnusedNeurons(intIdxN);
						
						vecUseNeurons = [vecUsedNeurons intNeuron];
						
						%decode
						if intDecoder == 1
							vecPerf(intIdxN) = doCrossValidatedDecodingLR(matUseR(:,vecUseNeurons),vecUseOri,2,vecPrior,dblLambda);
						elseif intDecoder == 2
							vecPerf(intIdxN) = doCrossValidatedDecodingML(matUseR(:,vecUseNeurons),vecUseOri,2,vecPrior);
						elseif intDecoder == 3
							vecPerf(intIdxN) = doCrossValidatedDecodingTM(matUseR(:,vecUseNeurons),vecUseOri,2,vecPrior);
						elseif intDecoder == 4
							vecPerf(intIdxN) = doCrossValidatedDecodingMD(matUseR(:,vecUseNeurons),vecUseOri,2); %add regularization to covar?
						end
					end
					
					%check which one to add
					[dblMaxPerf,intUseIdx]=max(vecPerf);
					intUseN = vecUnusedNeurons(intUseIdx);
					indNeuronsUsed(intUseN) = true;
					vecSelectNeuronOrder(intNeuronsUsed) = intUseN;
					vecGreedyPerf(intNeuronsUsed) = dblMaxPerf;
					fprintf('Decoded %s, neuron %d/%d (%d); current best perf is %.3f; selected neuron %d [%s]\n',strHighLow,intNeuronsUsed,intUseNeurons,intTotN,dblMaxPerf,intUseN,getTime);
				end
				if intSwitchQ == 0
					vecSelectNeuronOrderLow = vecSelectNeuronOrder;
					vecGreedyPerfLow = vecGreedyPerf;
				else
					vecSelectNeuronOrderHigh = vecSelectNeuronOrder;
					vecGreedyPerfHigh = vecGreedyPerf;
				end
			end
			if intDecoder == 1
				vecSelectNeuronOrderLowLR = vecSelectNeuronOrderLow;
				vecGreedyPerfLowLR = vecGreedyPerfLow;
				vecSelectNeuronOrderHighLR = vecSelectNeuronOrderHigh;
				vecGreedyPerfHighLR = vecGreedyPerfHigh;
			elseif intDecoder == 2
				vecSelectNeuronOrderLowML = vecSelectNeuronOrderLow;
				vecGreedyPerfLowML = vecGreedyPerfLow;
				vecSelectNeuronOrderHighML = vecSelectNeuronOrderHigh;
				vecGreedyPerfHighML = vecGreedyPerfHigh;
			elseif intDecoder == 3
				vecSelectNeuronOrderLowTM = vecSelectNeuronOrderLow;
				vecGreedyPerfLowTM = vecGreedyPerfLow;
				vecSelectNeuronOrderHighTM = vecSelectNeuronOrderHigh;
				vecGreedyPerfHighTM = vecGreedyPerfHigh;
			elseif intDecoder == 4
				vecSelectNeuronOrderLowMD = vecSelectNeuronOrderLow;
				vecGreedyPerfLowMD = vecGreedyPerfLow;
				vecSelectNeuronOrderHighMD = vecSelectNeuronOrderHigh;
				vecGreedyPerfHighMD = vecGreedyPerfHigh;
			end
		end
		
		%% plot
		for intUseDec = 1:intUseNumDecoders
			if intUseDec == 1
				strDec = 'LR';
				vecSelectNeuronOrderLow = vecSelectNeuronOrderLowLR;
				vecGreedyPerfLow = vecGreedyPerfLowLR;
				vecSelectNeuronOrderHigh = vecSelectNeuronOrderHighLR;
				vecGreedyPerfHigh = vecGreedyPerfHighLR;
			elseif intUseDec == 2
				strDec = 'ML';
				vecSelectNeuronOrderLow = vecSelectNeuronOrderLowML;
				vecGreedyPerfLow = vecGreedyPerfLowML;
				vecSelectNeuronOrderHigh = vecSelectNeuronOrderHighML;
				vecGreedyPerfHigh = vecGreedyPerfHighML;
			elseif intUseDec == 3
				strDec = 'TM';
				vecSelectNeuronOrderLow = vecSelectNeuronOrderLowTM;
				vecGreedyPerfLow = vecGreedyPerfLowTM;
				vecSelectNeuronOrderHigh = vecSelectNeuronOrderHighTM;
				vecGreedyPerfHigh = vecGreedyPerfHighTM;
			elseif intUseDec == 4
				strDec = 'MD';
				vecSelectNeuronOrderLow = vecSelectNeuronOrderLowMD;
				vecGreedyPerfLow = vecGreedyPerfLowMD;
				vecSelectNeuronOrderHigh = vecSelectNeuronOrderHighMD;
				vecGreedyPerfHigh = vecGreedyPerfHighMD;
			end
			
			figure;maxfig;
			subplot(2,3,1)
			vecX = 1:intUseNeurons;
			dblChanceP = sum((vecPrior./sum(vecPrior)).*vecPrior)./sum(vecPrior);
			hold on
			plot(vecX,vecGreedyPerfLow,'b');
			plot(vecX,vecGreedyPerfHigh,'r');
			plot(vecX([1 end]),dblChanceP.*[1 1],'k--');
			hold off
			title(sprintf('Greedy %s decoder; neurons sorted by performance',strDec));
			legend({'Low quantile','High quantile'},'location','best');
			xlabel('# of neurons used');
			ylabel('Fraction correct');
			fixfig;
			
			subplot(2,3,2)
			hold on
			plot(vecX,cumsum(sum(matLowR(:,vecSelectNeuronOrderLow),1)),'b')
			plot(vecX,cumsum(sum(matHighR(:,vecSelectNeuronOrderHigh),1)),'r')
			hold off
			title('Total # of spikes used');
			xlabel('# of neurons used');
			ylabel('# of spikes used');
			fixfig;
			
			%performance per spike
			subplot(2,3,3)
			hold on
			vecRelPerfLow = (vecGreedyPerfLow-dblChanceP)./cumsum(sum(matLowR(:,vecSelectNeuronOrderLow),1));
			vecRelPerfLow(vecRelPerfLow<0)=nan;
			vecRelPerfHigh = (vecGreedyPerfHigh-dblChanceP)./cumsum(sum(matHighR(:,vecSelectNeuronOrderHigh),1));
			vecRelPerfHigh(vecRelPerfHigh<0)=nan;
			plot(vecX,vecRelPerfLow,'b')
			plot(vecX,vecRelPerfHigh,'r')
			hold off
			title('Efficiency of neural code');
			xlabel('# of neurons used');
			ylabel('Fraction correct/spike used');
			fixfig;
			
			%performance per spike
			subplot(2,3,4)
			vecX_Low = cumsum(sum(matLowR(:,vecSelectNeuronOrderLow),1));
			vecX_High = cumsum(sum(matHighR(:,vecSelectNeuronOrderHigh),1));
			dblChanceP = sum((vecPrior./sum(vecPrior)).*vecPrior)./sum(vecPrior);
			hold on
			plot(vecX_Low,vecGreedyPerfLow,'b');
			plot(vecX_High,vecGreedyPerfHigh,'r');
			plot([vecX_Low(1) vecX_High(end)],dblChanceP.*[1 1],'k--');
			hold off
			title(sprintf('Neurons sorted by performance'));
			legend({'Low quantile','High quantile'},'location','best');
			xlabel('# of spikes used');
			ylabel('Fraction correct');
			fixfig;
			
			%%
			export_fig(fullpath(strFigurePathSR,sprintf('B5_%s_QuantileDecodingT%s_%s.tif',strDec,num2str(dblStartT),strRec)));
			export_fig(fullpath(strFigurePathSR,sprintf('B5_%s_QuantileDecodingT%s_%s.pdf',strDec,num2str(dblStartT),strRec)));
		end
	end
end
toc
