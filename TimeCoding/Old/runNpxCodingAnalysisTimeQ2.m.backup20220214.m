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


%% select all neurons in LP and drifting grating stimuli
if ~exist('sAggStim','var') || isempty(sAggStim)
	[sAggStim,sAggNeuron]=loadDataNpx('','driftinggrating');
end
sAggNeuron(strcmpi({sAggNeuron.SubjectType},'DBA')) = [];

%% pre-allocate matrices
intAreas = numel(cellUseAreas);
matDecPerf = [];

%% go through recordings
tic
for intRec=19%1:numel(sAggStim)
	% get matching recording data
	strRec = sAggStim(intRec).Exp;
	sThisRec = sAggStim(strcmpi(strRec,{sAggStim(:).Exp}));
	
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
		%% get spike times
		cellSpikeTimes = {sArea1Neurons.SpikeTimes};
		intNumN = numel(cellSpikeTimes);
		dblStimDur = roundi(median(vecStimOffTime - vecStimOnTime),1,'ceil');
		dblPreTime = 0;%0.3;
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
		cellSpikeTimesPerCellPerTrial = cell(intNumN,intTrialNum);
		for intN=1:intNumN
			[vecTrialPerSpike,vecTimePerSpike] = getSpikesInTrial(cellSpikeTimes{intN},vecStimOnTime-dblPreTime,dblMaxDur);
			for intTrial=1:intTrialNum
				cellSpikeTimesPerCellPerTrial{intN,intTrial} = vecTimePerSpike(vecTrialPerSpike==intTrial)-dblPreTime;
			end
		end
		%%
		boolPlot = false;
		vecRperTrial = nan(intTrialNum,1);
		vecSperTrial = nan(intTrialNum,1);
		vecHperTrial = nan(intTrialNum,1);
		vecLperTrial = nan(intTrialNum,1);
		vecMperTrial = nan(intTrialNum,1);
		
		%shuffled
		vecSperTrial_S = nan(intTrialNum,1);
		vecHperTrial_S = nan(intTrialNum,1);
		vecLperTrial_S = nan(intTrialNum,1);
		vecMperTrial_S = nan(intTrialNum,1);
		
		for intTrial=1:intTrialNum
			
			%intTrial=1
			vecAllSpikes = sort(cell2vec(cellSpikeTimesPerCellPerTrial(:,intTrial)));
			vecISI0 = abs([vecAllSpikes(1:(end-1)) - vecAllSpikes(2:end); inf]);
			vecAllSpikes(vecISI0==0)=vecAllSpikes(vecISI0==0)-(10^-5)*rand();
			vecAllSpikes = uniquetol(vecAllSpikes,1e-7);
			vecD = diff(vecAllSpikes);
			vecD1 = vecD(1:(end-1));
			vecD2 = vecD(2:end);
			[r,p]=corr(vecD1,vecD2);
			vecRperTrial(intTrial) = r;
			[vecTimeIFR,vecIFR] = getIFR(vecAllSpikes+dblPreTime,0,dblMaxDur,[],[],[],0);
			vecR_sorted = sort(vecIFR);
			dbl5perc = round(numel(vecR_sorted)/20);
			vecHperTrial(intTrial) = max(vecR_sorted);%mean(vecR_sorted((1+end-2*dbl5perc):(end-dbl5perc)));
			vecLperTrial(intTrial) = min(vecR_sorted);%mean(vecR_sorted((1+dbl5perc):(2*dbl5perc)));
			vecHperTrial(intTrial) = mean(vecR_sorted((1+end-1*dbl5perc):(end-0*dbl5perc)));
			vecLperTrial(intTrial) = mean(vecR_sorted((1+0*dbl5perc):(1*dbl5perc)));
			vecMperTrial(intTrial) = mean(vecR_sorted);
			vecSperTrial(intTrial) = std(vecR_sorted);
			
			%ISI
			vecISI1 = abs([vecAllSpikes(1:(end-1)) - vecAllSpikes(2:end); inf]);
			vecISI2 = abs([inf; vecAllSpikes(2:end) - vecAllSpikes(1:(end-1))]);
			vecISI = vecISI1(1:(end-1));%min([vecISI1 vecISI2],[],2);
			vecNSI = min([vecISI1 vecISI2],[],2); %nearest spike interval
			indSelect = true(1,numel(vecNSI));
			
			%shuffle to exponential distribution
			intIters = 1;
			for intIter=1%:intIters
				vecGenSpikes = cumsum([0;vecISI(randperm(numel(vecISI)))])+vecAllSpikes(1);
				[vecTimeIFRS,vecIFRS] = getIFR(vecGenSpikes+dblPreTime,0,dblMaxDur,[],[],[],0);
				
				%error add expected low/high values from shuffling
				vecR_sorted = sort(vecIFRS);
				dbl5perc = round(numel(vecR_sorted)/20);
				vecHperTrial_S(intTrial) = mean(vecR_sorted((1+end-1*dbl5perc):(end-0*dbl5perc)));
				vecLperTrial_S(intTrial) = mean(vecR_sorted((1+0*dbl5perc):(1*dbl5perc)));
				vecMperTrial_S(intTrial) = mean(vecR_sorted);
				vecSperTrial_S(intTrial) = std(vecR_sorted);
			end
			
			if boolPlot
			%%
			if numel(vecIFR)>numel(vecNSI)
				indSelect = cat(2,false,indSelect);
			end
			if numel(vecIFR)>(numel(vecNSI)+1)
				indSelect = cat(2,indSelect,false);
			end
			vecIFR= vecIFR(indSelect);
			
			%shuffle to exponential distribution
			intIters = 100;
			matNSIS = nan(numel(vecNSI),intIters);
			matIFRS = nan(numel(vecNSI),intIters);
			matISIS = nan(numel(vecNSI)-1,intIters);
			for intIter=1:intIters
				vecGenSpikes = cumsum([0;vecISI(randperm(numel(vecISI)))])+vecAllSpikes(1);
				vecISI1 = abs([vecGenSpikes(1:(end-1)) - vecGenSpikes(2:end); inf]);
				vecISI2 = abs([inf; vecGenSpikes(2:end) - vecGenSpikes(1:(end-1))]);
				matNSIS(:,intIter) = min([vecISI1 vecISI2],[],2); %nearest spike interval
				matISIS(:,intIter) = diff(vecGenSpikes); %ISI
				
				[vecTimeIFRS,vecIFRS] = getIFR(vecGenSpikes+dblPreTime,0,dblMaxDur,[],[],[],0);
				matIFRS(:,intIter) = vecIFRS(indSelect);
			end
			
			% plot
			clf;
			subplot(2,3,1)
			plot(vecAllSpikes,vecIFR)
			xlabel('Time after onset (s)');
			ylabel('Instant. firing rate (Hz)');
			title(sprintf('Population spiking rate, trial %d',intTrial));
			fixfig;
			
			subplot(2,3,2)
			dblBinSizeISI = 1/1000;
			vecBinEdgesISI = 0:dblBinSizeISI:0.015;
			vecBinCentersISI = vecBinEdgesISI(2:end)-dblBinSizeISI/2;
			dblLambda = 1./mean(vecISI);
			vecExpPdf = dblLambda.*exp(-dblLambda.*vecBinCentersISI);
			vecCounts = histcounts(vecISI,vecBinEdgesISI);
			hold on
			plot(vecBinCentersISI,vecCounts./sum(vecCounts(:)))
			plot(vecBinCentersISI,vecExpPdf./sum(vecExpPdf(:)))
			hold off
			set(gca,'yscale','log');
			xlabel('Inter-spike interval (s)');
			ylabel('Normalized count (n)');
			legend({'Observed','Exponential'});
			fixfig;
			
			
			subplot(2,3,3)
			scatter(vecIFR(2:end),vecISI,'.');
			hold on
			plot(sort(vecIFR),1./sort(vecIFR))
			hold off
			xlabel('Instant. firing rate (Hz)');
			ylabel('ISI');
			legend({'Observed','Exponential'});
			fixfig;
			
			dblBinSize = 25;
			vecBinEdges = 0:dblBinSize:(max(vecIFR)+dblBinSize);
			vecBinCenters = vecBinEdges(2:end)-dblBinSize/2;
			subplot(2,3,4)
			[vecCounts,vecMeans,vecSDs,cellVals,cellIDs] = makeBins(vecIFR,vecNSI,vecBinEdges);
			errorbar(vecBinCenters,vecMeans,vecSDs./sqrt(vecCounts))
			hold on
			[vecCountsS,vecMeansS,vecSDsS,cellValsS,cellIDsS] = makeBins(matIFRS(:),matNSIS(:),vecBinEdges);
			errorbar(vecBinCenters,vecMeansS,vecSDsS./sqrt(vecCountsS))
			%scatter(matIFR(:),matNSI(:))
			hold off
			xlabel('Instant. firing rate (Hz)');
			ylabel('Nearest-spike interval (s)');
			legend({'Observed','Exponential (shuffled)'});
			fixfig;
			
			subplot(2,3,5)
			dblStep = 50;
			vecBinEdgesIFR = 0:dblStep:1100;
			vecBinCentersIFR = vecBinEdgesIFR(2:end) - dblStep/2;
			[vecCounts,vecMeans,vecSDs,cellVals,cellIDs] = makeBins(vecIFR(2:end),vecISI,vecBinEdgesIFR);
			intMinSpikes = 10;
			indPlot = vecCounts>intMinSpikes;
			hold on
			errorbar(vecBinCentersIFR(indPlot),vecMeans(indPlot),vecSDs(indPlot)./sqrt(vecCounts(indPlot)))
			
			[vecCountsS,vecMeansS,vecSDsS,cellValsS,cellIDsS] = makeBins(flat(matIFRS(2:end,:)),flat(matISIS),vecBinEdgesIFR);
			indPlotS = vecCountsS>intMinSpikes;
			errorbar(vecBinCentersIFR(indPlotS),vecMeansS(indPlotS),vecSDsS(indPlotS)./sqrt(vecCountsS(indPlotS)))
			hold off
			xlabel('Instant. firing rate (Hz)');
			ylabel('Inter-spike interval (s)');
			legend({'Observed','Exponential (shuffled)'});
			fixfig;
			
			subplot(2,3,6)
			vecD = diff(vecAllSpikes);
			vecD1 = vecD(1:(end-1));
			vecD2 = vecD(2:end);
			scatter(vecD1,vecD2);
			[r,p]=corr(vecD1,vecD2);
			dblMaxT = 0.015;
			dblStep = 0.002;
			vecBinD = 0:dblStep:dblMaxT;
			vecBinD_c = vecBinD(2:end)-dblStep/2;
			[vecCounts,vecMeans,vecSDs,cellVals,cellIDs] = makeBins(vecD1,vecD2,vecBinD);
			errorbar(vecBinD_c,vecMeans,vecSDs./sqrt(vecCounts))
			xlabel('ISI spikes i,j (s)');
			ylabel('ISI spikes i+1,j+1');
			title(sprintf('ISI correlation r(d(i,j),d(i+1,j+1)), r=%.3f, p=%.3f',r,p));
			fixfig;
			
			pause
			end
		end
		
		%% compare with decoding
		matResp = cellfun(@numel,cellSpikeTimesPerCellPerTrial);
		intNeurons = size(matResp,1);
		vecPopSparseness = 1-mean(matResp,1).^2 ./ sum((matResp.^2)./intNeurons,1);
		
		dblLambda = 1;
		intTypeCV = 2;
		vecOri180 = mod(vecOrientation,180)*2;
		vecRepNum180 = vecRepNum(1:12)*2;
		intOriNum180 = intOriNum/2;
		[vecTrialTypeIdx,vecUnique,vecPriorDistribution,cellSelect,vecRepetition] = val2idx(vecOri180);
		[dblPerformanceCV,vecDecodedIndexCV,matPosteriorProbability,dblMeanErrorDegs,matConfusion,matWeights] = ...
				doCrossValidatedDecodingLR(matResp,vecOri180,intTypeCV,vecPriorDistribution,dblLambda);
		
		vecConfidence = nan(intTrialNum,1);
		for intTrial=1:intTrialNum
			vecConfidence(intTrial) = matPosteriorProbability(vecTrialTypeIdx(intTrial),intTrial);
		end
		dblActBinW = 50;
		vecActBins = 0:dblActBinW:700;
		vecActBinsC = vecActBins(2:end)-dblActBinW/2;
		
		%%
		figure;maxfig;
		subplot(2,3,1)
		histx(vecMperTrial)
		ylabel('Count (trials)')
		xlabel('Population activity, mean over time (Hz)')
		fixfig;
		
		subplot(2,3,2)
		histx(vecSperTrial)
		ylabel('Count (trials)')
		xlabel('Sd of population activity, mean over time (Hz)')
		fixfig;
		
		subplot(2,3,3)
		[vecCounts,vecMeansV,vecSDsV] = makeBins(vecMperTrial,vecSperTrial,vecActBins);
		indPlotBins = vecCounts>10;
		hold on
		scatter(vecMperTrial,vecSperTrial(:),[],[0.5 0.5 1],'.');
		errorbar(vecActBinsC(indPlotBins),vecMeansV(indPlotBins),vecSDsV(indPlotBins)./sqrt(vecCounts(indPlotBins)),'color',[0 0 1])
		xlabel('Mean population activity (Hz)')
		ylabel('Sd of population activity (Hz)')
		fixfig;
		
		
		subplot(2,3,4)
		vecCVperTrial = vecSperTrial./vecMperTrial;
		[vecCounts,vecMeansCV,vecSDsCV] = makeBins(vecMperTrial,vecCVperTrial,vecActBins);
		hold on
		scatter(vecMperTrial,vecCVperTrial(:),[],[0.5 0.5 1],'.');
		errorbar(vecActBinsC(indPlotBins),vecMeansCV(indPlotBins),vecSDsCV(indPlotBins)./sqrt(vecCounts(indPlotBins)),'color',[0 0 1])
		xlabel('Mean population activity (Hz)')
		ylabel('CV of population activity (Hz)')
		fixfig;
		
		%[r2b,p2b]=corr(vecMperTrial,vecCVperTrial);
		%title(sprintf('Corr(M,Sd)=%.3f, p=%.1e',r2b,p2b));
		
		%subplot(2,3,3)
		%scatter(abs(zscore(vecRperTrial)),vecConfidence)
		[r1,p1]=corr(vecMperTrial,vecConfidence(:))
		[r2,p2]=corr(vecConfidence,vecSperTrial)
		
		%%
		figure
		[rML,pML]=corr(vecMperTrial(:),vecLperTrial(:));
		[rMH,pMH]=corr(vecMperTrial,vecHperTrial(:));
		[rHL,pHL]=corr(vecLperTrial(:),vecHperTrial(:));
		
		subplot(2,3,1)
		[vecCounts,vecMeansL,vecSDsL] = makeBins(vecMperTrial,vecLperTrial,vecActBins);
		[vecCounts,vecMeansH,vecSDsH] = makeBins(vecMperTrial,vecHperTrial,vecActBins);
		[vecCounts_S,vecMeansL_S,vecSDsL_S] = makeBins(vecMperTrial_S,vecLperTrial_S,vecActBins);
		[vecCounts_S,vecMeansH_S,vecSDsH_S] = makeBins(vecMperTrial_S,vecHperTrial_S,vecActBins);hold on
		
		scatter(vecMperTrial,vecLperTrial_S(:),[],[0.7 0.7 0.7],'.');
		scatter(vecMperTrial,vecHperTrial_S(:),[],[0.7 0.7 0.7],'.');
		scatter(vecMperTrial,vecLperTrial(:),[],[0.5 0.5 1],'.');
		scatter(vecMperTrial,vecHperTrial(:),[],[1 0.5 0.5],'.');
		
		errorbar(vecActBinsC(indPlotBins),vecMeansL_S(indPlotBins),vecSDsL_S(indPlotBins)./sqrt(vecCounts_S(indPlotBins)),'color',[0.5 0.5 0.5])
		errorbar(vecActBinsC(indPlotBins),vecMeansH_S(indPlotBins),vecSDsH_S(indPlotBins)./sqrt(vecCounts_S(indPlotBins)),'color',[0.5 0.5 0.5])
		errorbar(vecActBinsC(indPlotBins),vecMeansL(indPlotBins),vecSDsL(indPlotBins)./sqrt(vecCounts(indPlotBins)),'color',[0 0 1])
		errorbar(vecActBinsC(indPlotBins),vecMeansH(indPlotBins),vecSDsH(indPlotBins)./sqrt(vecCounts(indPlotBins)),'color',[1 0 0])
		
		hold off
		fixfig;
		
		subplot(2,3,2)
		[vecCounts,vecMeansL,vecSDsL] = makeBins(vecMperTrial,vecLperTrial(:)./vecMperTrial,vecActBins);
		[vecCounts,vecMeansH,vecSDsH] = makeBins(vecMperTrial,vecHperTrial(:)./vecMperTrial,vecActBins);
		[vecCounts_S,vecMeansL_S,vecSDsL_S] = makeBins(vecMperTrial_S,vecLperTrial_S(:)./vecMperTrial_S,vecActBins);
		[vecCounts_S,vecMeansH_S,vecSDsH_S] = makeBins(vecMperTrial_S,vecHperTrial_S(:)./vecMperTrial_S,vecActBins);hold on
		
		hold on
		scatter(vecMperTrial,vecLperTrial_S(:)./vecMperTrial,[],[0.7 0.7 0.7],'.');
		scatter(vecMperTrial,vecHperTrial_S(:)./vecMperTrial,[],[0.7 0.7 0.7],'.');
		scatter(vecMperTrial,vecLperTrial(:)./vecMperTrial,[],[0.5 0.5 1],'.');
		scatter(vecMperTrial,vecHperTrial(:)./vecMperTrial,[],[1 0.5 0.5],'.');
		
		errorbar(vecActBinsC(indPlotBins),vecMeansL_S(indPlotBins),vecSDsL_S(indPlotBins)./sqrt(vecCounts_S(indPlotBins)),'color',[0.5 0.5 0.5])
		errorbar(vecActBinsC(indPlotBins),vecMeansH_S(indPlotBins),vecSDsH_S(indPlotBins)./sqrt(vecCounts_S(indPlotBins)),'color',[0.5 0.5 0.5])
		errorbar(vecActBinsC(indPlotBins),vecMeansL(indPlotBins),vecSDsL(indPlotBins)./sqrt(vecCounts(indPlotBins)),'color',[0 0 1])
		errorbar(vecActBinsC(indPlotBins),vecMeansH(indPlotBins),vecSDsH(indPlotBins)./sqrt(vecCounts(indPlotBins)),'color',[1 0 0])
		hold off
		fixfig;
		
		subplot(2,3,3)
		scatter(vecLperTrial,vecHperTrial);
		
		subplot(2,3,4)
		r1=corr(vecLperTrial(:),vecPopSparseness(:))
		r2=corr(vecMperTrial(:),vecPopSparseness(:))
		r3=corr(vecHperTrial(:),vecPopSparseness(:))
		scatter(vecMperTrial,vecPopSparseness);
		
		subplot(2,3,5)
		scatter(vecLperTrial,vecPopSparseness);
		
		subplot(2,3,6)
		scatter(vecHperTrial,vecPopSparseness);
		
		return

	end
end
toc
