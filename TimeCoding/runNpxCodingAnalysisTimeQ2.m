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
	strFigurePath = 'F:\Data\Results\PopTimeCoding';
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
		%%
		boolPlot = true;
		vecRperTrial = nan(intTrialNum,1);
		vecSperTrial = nan(intTrialNum,1);
		vecHperTrial = nan(intTrialNum,1);
		vecLperTrial = nan(intTrialNum,1);
		vecMperTrial = nan(intTrialNum,1);
			
		%shuffling ISIs within one trial for each neuron: spike count is identical at trial level
		vecSperTrial_S = nan(intTrialNum,1);
		vecHperTrial_S = nan(intTrialNum,1);
		vecLperTrial_S = nan(intTrialNum,1);
		vecMperTrial_S = nan(intTrialNum,1);
		
		%shuffl
		vecSperTrial_SN = nan(intTrialNum,1);
		vecHperTrial_SN = nan(intTrialNum,1);
		vecLperTrial_SN = nan(intTrialNum,1);
		vecMperTrial_SN = nan(intTrialNum,1);
		
		%shuffl
		vecSperTrial_SS = nan(intTrialNum,1);
		vecHperTrial_SS = nan(intTrialNum,1);
		vecLperTrial_SS = nan(intTrialNum,1);
		vecMperTrial_SS = nan(intTrialNum,1);
		
		for intTrial=1:intTrialNum
			
			%real
			vecAllSpikes = sort(cell2vec(cellSpikeTimesPerCellPerTrial(:,intTrial)));
			vecISI0 = [vecAllSpikes(2:end) - vecAllSpikes(1:(end-1)); inf];
			vecAllSpikes(vecISI0==0)=vecAllSpikes(vecISI0==0)-(10^-5)*rand();
			vecAllSpikes = uniquetol(vecAllSpikes,1e-7);
			intSpikes = numel(vecAllSpikes);
			vecD = diff(vecAllSpikes);
			vecD1 = vecD(1:(end-1));
			vecD2 = vecD(2:end);
			[r,p]=corr(vecD1,vecD2);
			vecRperTrial(intTrial) = r;
			vecIFR = cellIFR_perTrial{intTrial};
			vecTimeIFR = cellTimeIFR_perTrial{intTrial};
			
			vecR_sorted = sort(vecIFR);
			intHighQ = round(numel(vecR_sorted)/20);
			intLowQ = round(numel(vecR_sorted)/20);
			vecHperTrial(intTrial) = mean(vecR_sorted((1+end-1*intHighQ):(end-0*intHighQ)));
			vecLperTrial(intTrial) = mean(vecR_sorted((1+0*intLowQ):(1*intLowQ)));
			vecMperTrial(intTrial) = mean(vecR_sorted);
			vecSperTrial(intTrial) = std(vecR_sorted);
			
			%ISI
			vecISI1 = abs([vecAllSpikes(1:(end-1)) - vecAllSpikes(2:end); inf]);
			vecISI2 = abs([inf; vecAllSpikes(2:end) - vecAllSpikes(1:(end-1))]);
			vecISI = vecISI1(1:(end-1));%min([vecISI1 vecISI2],[],2);
			vecNSI = min([vecISI1 vecISI2],[],2); %nearest spike interval
			
			%shuffled single-neuron ISIs
			vecIFRS = cellIFR_perTrial_S{intTrial};
			vecR_sorted = sort(vecIFRS);
			intHighQ = round(numel(vecR_sorted)/20);
			vecHperTrial_S(intTrial) = mean(vecR_sorted((1+end-1*intHighQ):(end-0*intHighQ)));
			vecLperTrial_S(intTrial) = mean(vecR_sorted((1+0*intHighQ):(1*intHighQ)));
			vecMperTrial_S(intTrial) = mean(vecR_sorted);
			vecSperTrial_S(intTrial) = std(vecR_sorted);
			
			%shuffled overall ISIs
			vecIFRS = cellIFR_perTrial_SN{intTrial};
			vecR_sorted = sort(vecIFRS);
			intHighQ = round(numel(vecR_sorted)/20);
			vecHperTrial_SN(intTrial) = mean(vecR_sorted((1+end-1*intHighQ):(end-0*intHighQ)));
			vecLperTrial_SN(intTrial) = mean(vecR_sorted((1+0*intHighQ):(1*intHighQ)));
			vecMperTrial_SN(intTrial) = mean(vecR_sorted);
			vecSperTrial_SN(intTrial) = std(vecR_sorted);
			
			%shuffled overall ISIs, single-neuron ISIs
			vecIFRS = cellIFR_perTrial_SS{intTrial};
			vecR_sorted = sort(vecIFRS);
			intHighQ = round(numel(vecR_sorted)/20);
			vecHperTrial_SS(intTrial) = mean(vecR_sorted((1+end-1*intHighQ):(end-0*intHighQ)));
			vecLperTrial_SS(intTrial) = mean(vecR_sorted((1+0*intHighQ):(1*intHighQ)));
			vecMperTrial_SS(intTrial) = mean(vecR_sorted);
			vecSperTrial_SS(intTrial) = std(vecR_sorted);
			
			
			if boolPlot
				%%
				%shuffle to exponential distribution
				intIters = 100;
				matNSIS = nan(intSpikes,intIters);
				matIFRS = nan(intSpikes+2,intIters);
				matISIS = nan(intSpikes-1,intIters);
				for intIter=1:intIters
					vecGenSpikes = cumsum([0;vecISI(randperm(numel(vecISI)))])+vecAllSpikes(1);
					vecISI1 = abs([vecGenSpikes(1:(end-1)) - vecGenSpikes(2:end); inf]);
					vecISI2 = abs([inf; vecGenSpikes(2:end) - vecGenSpikes(1:(end-1))]);
					matNSIS(:,intIter) = min([vecISI1 vecISI2],[],2); %nearest spike interval
					matISIS(:,intIter) = diff(vecGenSpikes); %ISI
					
					[vecTimeIFRS,vecIFRS] = getIFR(vecGenSpikes,0,dblMaxDur,[],[],[],0);
					matIFRS(:,intIter) = vecIFRS;
				end
				
				% plot
				figure;maxfig;
				subplot(2,3,1)
				plot(vecTimeIFR+dblStartT,vecIFR)
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
				legend({'Observed','Theory (Exponential)'});
				fixfig;
				
				
				subplot(2,3,3)
				scatter(vecIFR(2:(end-2)),vecISI,'.');
				hold on
				plot(sort(vecIFR),1./sort(vecIFR))
				hold off
				xlabel('Instant. firing rate (Hz)');
				ylabel('ISI');
				legend({'Observed','Theory (Exponential)'});
				fixfig;
				
				dblBinSize = 25;
				vecBinEdges = 0:dblBinSize:(max(vecIFR)+dblBinSize);
				vecBinCenters = vecBinEdges(2:end)-dblBinSize/2;
				subplot(2,3,4)
				[vecCounts,vecMeans,vecSDs,cellVals,cellIDs] = makeBins(vecIFR(2:(end-1)),vecNSI,vecBinEdges);
				errorbar(vecBinCenters,vecMeans,vecSDs./sqrt(vecCounts))
				hold on
				[vecCountsS,vecMeansS,vecSDsS,cellValsS,cellIDsS] = makeBins(flat(matIFRS(2:(end-1),:)),matNSIS(:),vecBinEdges);
				errorbar(vecBinCenters,vecMeansS,vecSDsS./sqrt(vecCountsS))
				%scatter(matIFR(:),matNSI(:))
				hold off
				xlabel('Instant. firing rate (Hz)');
				ylabel('Nearest-spike interval (s)');
				legend({'Observed','Shuffled ISIs'});
				fixfig;
				
				subplot(2,3,5)
				dblStep = 50;
				vecBinEdgesIFR = 0:dblStep:1100;
				vecBinCentersIFR = vecBinEdgesIFR(2:end) - dblStep/2;
				[vecCounts,vecMeans,vecSDs,cellVals,cellIDs] = makeBins(vecIFR(2:(end-2)),vecISI,vecBinEdgesIFR);
				intMinSpikes = 10;
				indPlot = vecCounts>intMinSpikes;
				hold on
				errorbar(vecBinCentersIFR(indPlot),vecMeans(indPlot),vecSDs(indPlot)./sqrt(vecCounts(indPlot)))
				
				[vecCountsS,vecMeansS,vecSDsS,cellValsS,cellIDsS] = makeBins(flat(matIFRS(2:(end-2),:)),flat(matISIS),vecBinEdgesIFR);
				indPlotS = vecCountsS>intMinSpikes;
				errorbar(vecBinCentersIFR(indPlotS),vecMeansS(indPlotS),vecSDsS(indPlotS)./sqrt(vecCountsS(indPlotS)))
				hold off
				xlabel('Instant. firing rate (Hz)');
				ylabel('Inter-spike interval (s)');
				legend({'Observed','Shuffled ISIs'});
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
				
				export_fig(fullpath(strFigurePath,sprintf('B1_ExampleActivityT%s_%sTrial%d.tif',num2str(dblStartT),strRec,intTrial)));
				export_fig(fullpath(strFigurePath,sprintf('B1_ExampleActivityT%s_%sTrial%d.pdf',num2str(dblStartT),strRec,intTrial)));
				boolPlot = false;
			end
		end
		
		%% calc sparseness
		matResp = cellfun(@numel,cellSpikeTimesPerCellPerTrial);
		intNeurons = size(matResp,1);
		vecPopSparseness = 1-mean(matResp,1).^2 ./ sum((matResp.^2)./intNeurons,1);
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
		
		export_fig(fullpath(strFigurePath,sprintf('B2_PopActStatisticsT%s_%s.tif',num2str(dblStartT),strRec)));
		export_fig(fullpath(strFigurePath,sprintf('B2_PopActStatisticsT%s_%s.pdf',num2str(dblStartT),strRec)));
		
		%%
		figure;maxfig;
		[rML,pML]=corr(vecMperTrial(:),vecLperTrial(:));
		[rMH,pMH]=corr(vecMperTrial,vecHperTrial(:));
		[rHL,pHL]=corr(vecLperTrial(:),vecHperTrial(:));
		
		subplot(2,3,1)
		[vecCounts,vecMeansL,vecSDsL] = makeBins(vecMperTrial,vecLperTrial,vecActBins);
		[vecCounts,vecMeansH,vecSDsH] = makeBins(vecMperTrial,vecHperTrial,vecActBins);
		[vecCounts_S,vecMeansL_S,vecSDsL_S] = makeBins(vecMperTrial_S,vecLperTrial_S,vecActBins);
		[vecCounts_S,vecMeansH_S,vecSDsH_S] = makeBins(vecMperTrial_S,vecHperTrial_S,vecActBins);hold on
		
		scatter(cat(1,vecMperTrial_S,vecMperTrial_S),cat(1,vecLperTrial_S(:),vecHperTrial_S(:)),[],[0.7 0.7 0.7],'.');
		scatter(vecMperTrial,vecLperTrial(:),[],[0.5 0.5 1],'.');
		scatter(vecMperTrial,vecHperTrial(:),[],[1 0.5 0.5],'.');
		
		errorbar(vecActBinsC(indPlotBins),vecMeansL_S(indPlotBins),vecSDsL_S(indPlotBins)./sqrt(vecCounts_S(indPlotBins)),'color',[0.5 0.5 0.5])
		errorbar(vecActBinsC(indPlotBins),vecMeansH_S(indPlotBins),vecSDsH_S(indPlotBins)./sqrt(vecCounts_S(indPlotBins)),'color',[0.5 0.5 0.5])
		errorbar(vecActBinsC(indPlotBins),vecMeansL(indPlotBins),vecSDsL(indPlotBins)./sqrt(vecCounts(indPlotBins)),'color',[0 0 1])
		errorbar(vecActBinsC(indPlotBins),vecMeansH(indPlotBins),vecSDsH(indPlotBins)./sqrt(vecCounts(indPlotBins)),'color',[1 0 0])
		
		hold off
		fixfig;
		xlabel('Mean firing rate during trial (Hz)');
		ylabel('Firing rate during upper/lower 5% (Hz)');
		legend({'Shuffled neuron ISIs','Lowest 5%','Highest 5%'},'location','best');
		fixfig;
		
		subplot(2,3,2)
		dblStep = 0.05;
		vecBinE = -0.5:dblStep:1.5;
		vecBinC = vecBinE(2:end)-dblStep/2;
		vecCountsLow = histcounts((vecLperTrial-vecLperTrial_S)./vecLperTrial_S,vecBinE);
		vecCountsHigh = histcounts((vecHperTrial-vecHperTrial_S)./vecHperTrial_S,vecBinE);
		plot(vecBinC*100,vecCountsLow,'color',[0 0 1]);
		hold on
		plot(vecBinC*100,vecCountsHigh,'color',[1 0 0]);
		hold off
		xlabel('% change in firing rate over shuffled pop ISIs');
		legend({'Lowest 5%','Highest 5%'});
		ylabel('Number of trials (count)');
		vecL = (vecLperTrial-vecLperTrial_S)./vecLperTrial_S;
		vecH = (vecHperTrial-vecHperTrial_S)./vecHperTrial_S;
		[h,pL]=ttest(vecL);
		[h,pH]=ttest(vecH);
		title(sprintf('Low, mean=%.1f%%, p=%.1e; high, mean=+%.1f%%, p=%.1e',mean(vecL)*100,pL,mean(vecH)*100,pH));
		fixfig;
		
		
		subplot(2,3,3)
		r1=corr(vecLperTrial(:),vecPopSparseness(:));
		r2=corr(vecMperTrial(:),vecPopSparseness(:));
		[rSpH,pSpH]=corr(vecHperTrial(:),vecPopSparseness(:));
		vecActBinsH = 0:dblActBinW:1700;
		vecActBinsHC = vecActBinsH(2:end)-dblActBinW/2;
		[vecCounts2,vecMeans2,vecSDs2] = makeBins(vecHperTrial,vecPopSparseness,vecActBinsH);
		indPlotBins2 = vecCounts2>10;
		mdl = fitlm(vecHperTrial,vecPopSparseness);
		ci = coefCI(mdl);
		[ypred,yci] = predict(mdl,vecActBinsHC(indPlotBins2)');
		
		hold on;
		scatter(vecHperTrial,vecPopSparseness,[],1-(1-lines(1))*(2/3),'.');
		%errorbar(vecActBinsHC(indPlotBins2),vecMeans2(indPlotBins2),vecSDs2(indPlotBins2)./sqrt(vecCounts2(indPlotBins2)),'color',lines(1));
		plot(vecActBinsHC(indPlotBins2),yci(:,1),'--','color',lines(1));
		plot(vecActBinsHC(indPlotBins2),ypred,'color',lines(1));
		plot(vecActBinsHC(indPlotBins2),yci(:,2),'--','color',lines(1));
		hold off;
		ylabel('Population sparseness per trial');
		xlabel('Max. pop. firing rate per trial (Hz) ');
		title(sprintf('Corr(sparse, max rate); r=%.3f, p=%.1e',rSpH,pSpH));
		fixfig;
		
		%split population in highest 50/lowest 50
		intUseUpperCells = min(sum(matResp>0,1));
		vecHighAct = nan(intTrialNum,1);
		vecLowAct =  nan(intTrialNum,1);
		vecQuantiles = [1/3 1/2 2/3];
		vecQuantileIdx = round(vecQuantiles*intNeurons);
		matQuantileAct = nan(intTrialNum,numel(vecQuantileIdx));
		vecMeanOfActiveCells = nan(intTrialNum,1);
		for intTrial=1:intTrialNum
			vecR = sort(matResp(:,intTrial));
			matQuantileAct(intTrial,:) = vecR(vecQuantileIdx);
			vecMeanOfActiveCells(intTrial) = mean(vecR((end-intUseUpperCells+1):end));
		end
		[vecHsorted,vecReorder]=sort(vecHperTrial);
		
		
		%subplot(2,3,4);
		%cellSpikeTimesPerCellPerTrial_Shuffled{intN,intTrial} = vecTimePerSpikeS(vecTrialPerSpikeS==intTrial)-dblPreTime;
		%scatter(vecHsorted,mean(matResp(:,vecReorder)));
		
		%subplot(2,3,5);
		%scatter(vecHsorted,sum(matResp(:,vecReorder)>0));
		
		%subplot(2,3,6);
		%scatter(vecHsorted,vecMeanOfActiveCells(vecReorder),'.')
		
		%% cell-shuffled
		matRespS = cellfun(@numel,cellSpikeTimesPerCellPerTrial_S);
		intNeuronsS = size(matRespS,1);
		vecPopSparsenessS = 1-mean(matRespS,1).^2 ./ sum((matRespS.^2)./intNeuronsS,1);
		
		[rML,pML]=corr(vecMperTrial_SN(:),vecLperTrial_SN(:));
		[rMH,pMH]=corr(vecMperTrial_SN,vecHperTrial_SN(:));
		[rHL,pHL]=corr(vecLperTrial_SN(:),vecHperTrial_SN(:));
		
		subplot(2,3,4)
		[vecCounts,vecMeansL,vecSDsL] = makeBins(vecMperTrial_SN,vecLperTrial_SN,vecActBins);
		[vecCounts,vecMeansH,vecSDsH] = makeBins(vecMperTrial_SN,vecHperTrial_SN,vecActBins);
		[vecCounts_S,vecMeansL_S,vecSDsL_S] = makeBins(vecMperTrial_SS,vecLperTrial_SS,vecActBins);
		[vecCounts_S,vecMeansH_S,vecSDsH_S] = makeBins(vecMperTrial_SS,vecHperTrial_SS,vecActBins);hold on
		
		scatter(cat(1,vecMperTrial_SS,vecMperTrial_SS),cat(1,vecLperTrial_SS(:),vecHperTrial_SS(:)),[],[0.7 0.7 0.7],'.');
		scatter(vecMperTrial_SN,vecLperTrial_SN(:),[],[0.5 0.5 1],'.');
		scatter(vecMperTrial_SN,vecHperTrial_SN(:),[],[1 0.5 0.5],'.');
		
		errorbar(vecActBinsC(indPlotBins),vecMeansL_S(indPlotBins),vecSDsL_S(indPlotBins)./sqrt(vecCounts_S(indPlotBins)),'color',[0.5 0.5 0.5])
		errorbar(vecActBinsC(indPlotBins),vecMeansH_S(indPlotBins),vecSDsH_S(indPlotBins)./sqrt(vecCounts_S(indPlotBins)),'color',[0.5 0.5 0.5])
		errorbar(vecActBinsC(indPlotBins),vecMeansL(indPlotBins),vecSDsL(indPlotBins)./sqrt(vecCounts(indPlotBins)),'color',[0 0 1])
		errorbar(vecActBinsC(indPlotBins),vecMeansH(indPlotBins),vecSDsH(indPlotBins)./sqrt(vecCounts(indPlotBins)),'color',[1 0 0])
		
		hold off
		fixfig;
		xlabel('Mean firing rate during trial (Hz)');
		ylabel('Firing rate during upper/lower 5% (Hz)');
		legend({'Shuffled pop+neuron ISIs','Lowest 5%','Highest 5%'},'location','best');
		title('Neuron-ISI shuffle control');
		fixfig;
		
		subplot(2,3,5)
		dblStep = 0.05;
		vecBinE = -0.5:dblStep:1.5;
		vecBinC = vecBinE(2:end)-dblStep/2;
		vecCountsLow = histcounts((vecLperTrial_SN-vecLperTrial_SS)./vecLperTrial_SS,vecBinE);
		vecCountsHigh = histcounts((vecHperTrial_SN-vecHperTrial_SS)./vecHperTrial_SS,vecBinE);
		plot(vecBinC*100,vecCountsLow,'color',[0 0 1]);
		hold on
		plot(vecBinC*100,vecCountsHigh,'color',[1 0 0]);
		hold off
		xlabel('% change in firing rate over shuffled pop ISIs');
		legend({'Lowest 5%','Highest 5%'});
		ylabel('Number of trials (count)');
		vecL = (vecLperTrial_SN-vecLperTrial_SS)./vecLperTrial_SS;
		vecH = (vecHperTrial_SN-vecHperTrial_SS)./vecHperTrial_SS;
		[h,pL]=ttest(vecL);
		[h,pH]=ttest(vecH);
		title(sprintf('Low, mean=%.1f%%, p=%.1e; high, mean=+%.1f%%, p=%.1e',mean(vecL)*100,pL,mean(vecH)*100,pH));
		fixfig;
		
		subplot(2,3,6)
		r1=corr(vecLperTrial_SN(:),vecPopSparsenessS(:));
		r2=corr(vecMperTrial_SN(:),vecPopSparsenessS(:));
		[rSpH,pSpH]=corr(vecHperTrial_SN(:),vecPopSparsenessS(:));
		vecActBinsH = 0:dblActBinW:1700;
		vecActBinsHC = vecActBinsH(2:end)-dblActBinW/2;
		[vecCounts2,vecMeans2,vecSDs2] = makeBins(vecHperTrial_SN,vecPopSparsenessS,vecActBinsH);
		indPlotBins2 = vecCounts2>10;
		mdl = fitlm(vecHperTrial_SN,vecPopSparsenessS);
		ci = coefCI(mdl);
		[ypred,yci] = predict(mdl,vecActBinsHC(indPlotBins2)');
		
		hold on;
		scatter(vecHperTrial_SN,vecPopSparsenessS,[],1-(1-lines(1))*(2/3),'.');
		%errorbar(vecActBinsHC(indPlotBins2),vecMeans2(indPlotBins2),vecSDs2(indPlotBins2)./sqrt(vecCounts2(indPlotBins2)),'color',lines(1));
		plot(vecActBinsHC(indPlotBins2),yci(:,1),'--','color',lines(1));
		plot(vecActBinsHC(indPlotBins2),ypred,'color',lines(1));
		plot(vecActBinsHC(indPlotBins2),yci(:,2),'--','color',lines(1));
		hold off;
		ylabel('Population sparseness per trial');
		xlabel('Max. pop. firing rate per trial (Hz) ');
		title(sprintf('Corr(sparse, max rate); r=%.3f, p=%.1e',rSpH,pSpH));
		fixfig;
		
		export_fig(fullpath(strFigurePath,sprintf('B3_QuantileDeviationsT%s_%s.tif',num2str(dblStartT),strRec)));
		export_fig(fullpath(strFigurePath,sprintf('B3_QuantileDeviationsT%s_%s.pdf',num2str(dblStartT),strRec)));
		
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
		
		%% calculate prior & decode with full population
		%low
		intTypeCV = 2;
		dblLambda = 1;
		[vecTrialTypeIdxLow,vecUnique,vecPriorDistributionLow,cellSelect,vecRepetition] = val2idx(vecOriLow);
		[dblPerformanceLow,vecDecodedIndexLow,matPosteriorProbabilityLow,dblMeanErrorDegsLow,matConfusionLow,matWeightsLow] = ...
			doCrossValidatedDecodingLR(matLowR,vecOriLow,intTypeCV,vecPriorDistributionLow,dblLambda);
		vecConfidenceLow = nan(intTrialNum,1);
		for intTrial=1:intTrialNum
			vecConfidenceLow(intTrial) = matPosteriorProbabilityLow(vecTrialTypeIdxLow(intTrial),intTrial);
		end
		
		%high
		[vecTrialTypeIdxHigh,vecUnique,vecPriorDistributionHigh,cellSelect,vecRepetition] = val2idx(vecOriHigh);
		[dblPerformanceHigh,vecDecodedIndexHigh,matPosteriorProbabilityHigh,dblMeanErrorDegsHigh,matConfusionHigh,matWeightsHigh] = ...
			doCrossValidatedDecodingLR(matHighR,vecOriHigh,intTypeCV,vecPriorDistributionHigh,dblLambda);
		vecConfidenceHigh = nan(intTrialNum,1);
		for intTrial=1:intTrialNum
			vecConfidenceHigh(intTrial) = matPosteriorProbabilityHigh(vecTrialTypeIdxHigh(intTrial),intTrial);
		end
		
		% ori tuning is stable across quantiles
		sOutLow = getTuningCurves(matLowR',vecOri180,0);
		vecPrefRadLow = sOutLow.matFittedParams(:,1);
		sOutHigh = getTuningCurves(matHighR',vecOri180,0);
		vecPrefRadHigh = sOutHigh.matFittedParams(:,1);
		%ori tuning is stable
		
		% ori tuning diff within low
		intHalfTrials = floor(intTrialNum/2);
		sOutLow1 = getTuningCurves(matLowR(1:intHalfTrials,:)',vecOri180(1:intHalfTrials),0);
		vecPrefRadLow1 = sOutLow1.matFittedParams(:,1);
		sOutLow2 = getTuningCurves(matLowR((intHalfTrials+1):end,:)',vecOri180((intHalfTrials+1):end),0);
		vecPrefRadLow2 = sOutLow2.matFittedParams(:,1);
		% ori tuning diff within high
		intHalfTrials = floor(intTrialNum/2);
		sOutHigh1 = getTuningCurves(matHighR(1:intHalfTrials,:)',vecOri180(1:intHalfTrials),0);
		vecPrefRadHigh1 = sOutHigh1.matFittedParams(:,1);
		sOutHigh2 = getTuningCurves(matHighR((intHalfTrials+1):end,:)',vecOri180((intHalfTrials+1):end),0);
		vecPrefRadHigh2 = sOutHigh2.matFittedParams(:,1);
		
		%% plot
		figure;maxfig;
		subplot(2,3,1);
		vecEdges = 0:90:360;
		[matCounts,matValMeans,matValSDs,cellVals,cellIDs] = ...
			makeBins2(vecPrefRadLow,vecPrefRadHigh,ones(size(vecPrefRadHigh)),vecEdges,vecEdges);
		%imagesc(matCounts);axis xy
		scatter(rad2deg(vecPrefRadLow),rad2deg(vecPrefRadHigh),'x');
		xlim([-10 370]);
		set(gca,'xtick',0:90:360);
		set(gca,'ytick',0:90:360);
		ylim([-10 370]);
		xlabel('Preferred ori low q')
		ylabel('Preferred ori high q')
		fixfig;
		
		% ori tuning diff within low
		subplot(2,3,2);
		scatter(rad2deg(vecPrefRadLow1),rad2deg(vecPrefRadLow2),'x');
		xlim([-10 370]);
		set(gca,'xtick',0:90:360);
		set(gca,'ytick',0:90:360);
		ylim([-10 370]);
		xlabel('Preferred ori low q, 1st half')
		ylabel('Preferred ori low q, 2nd half')
		fixfig;
		
		
		% ori tuning diff within high
		subplot(2,3,3);
		scatter(rad2deg(vecPrefRadHigh1),rad2deg(vecPrefRadHigh2),'x');
		xlim([-10 370]);
		set(gca,'xtick',0:90:360);
		set(gca,'ytick',0:90:360);
		ylim([-10 370]);
		xlabel('Preferred ori high q, 1st half')
		ylabel('Preferred ori high q, 2nd half')
		fixfig;
		
		% ori tuning diff within high
		subplot(2,3,4);
		scatter(rad2deg(vecPrefRadHigh1),rad2deg(vecPrefRadLow2),'x');
		xlim([-10 370]);
		set(gca,'xtick',0:90:360);
		set(gca,'ytick',0:90:360);
		ylim([-10 370]);
		xlabel('Preferred ori high q, 1st half')
		ylabel('Preferred ori low q, 2nd half')
		fixfig;
		
		vecDiffHL1 = abs(circ_dist(vecPrefRadHigh1,vecPrefRadLow2));
		vecDiffHL2 = abs(circ_dist(vecPrefRadHigh2,vecPrefRadLow1));
		vecDiffHL = rad2deg(vecDiffHL1);%+vecDiffHL2)/2;
		vecDiffLL = rad2deg(abs(circ_dist(vecPrefRadLow1,vecPrefRadLow2)));
		vecDiffHH = rad2deg(abs(circ_dist(vecPrefRadHigh1,vecPrefRadHigh2)));
		
		dblBinS = 22.5;
		vecBinE = 0:dblBinS:180;
		vecBinC = vecBinE(2:end)-dblBinS/2;
		vecCHL = histcounts(vecDiffHL,vecBinE);
		vecCLL = histcounts(vecDiffLL,vecBinE);
		vecCHH = histcounts(vecDiffHH,vecBinE);
		
		[h,pHL_LL]=ttest(vecDiffHL,vecDiffLL);
		[h,pHL_HH]=ttest(vecDiffHL,vecDiffHH);
		[h,pHH_LL]=ttest(vecDiffHH,vecDiffLL);
		subplot(2,3,5)
		hold on
		plot([0.5 3.5],[90 90],'--','color',[0.5 0.5 0.5]);
		errorbar(1,mean(vecDiffHL),std(vecDiffHL)./sqrt(numel(vecDiffHL)),'xk');
		errorbar(2,mean(vecDiffLL),std(vecDiffLL)./sqrt(numel(vecDiffLL)),'xb');
		errorbar(3,mean(vecDiffHH),std(vecDiffHH)./sqrt(numel(vecDiffHH)),'xr');
		hold off
		ylabel('Angular diff. pref. ori.');
		set(gca,'xtick',[1 2 3],'xticklabel',{'High-low','Low-low','High-high'});
		title(sprintf('T-tests: HL-LL,p=%.3f, HL-HH,p=%.3f, HH-LL,p=%.3f',pHL_LL,pHL_HH,pHH_LL));
		fixfig;
		
		%% cross-decode high-low & low-high
		%low weights, low act
		matDataPlusLin = [matLowR'; ones(1,size(matLowR',2))];
		matActivation = matWeightsLow'*matDataPlusLin;
		matPP_LowLow = exp(bsxfun(@minus,matActivation,logsumexp(matActivation,1))); %softmax
		vecConfidenceLowLow = nan(intTrialNum,1);
		for intTrial=1:intTrialNum
			vecConfidenceLowLow(intTrial) = matPP_LowLow(vecTrialTypeIdxLow(intTrial),intTrial);
		end
		
		%high weights, high act
		matDataPlusLin = [matHighR'; ones(1,size(matHighR',2))];
		matActivation = matWeightsHigh'*matDataPlusLin;
		matPP_HighHigh = exp(bsxfun(@minus,matActivation,logsumexp(matActivation,1))); %softmax
		vecConfidenceHighHigh = nan(intTrialNum,1);
		for intTrial=1:intTrialNum
			vecConfidenceHighHigh(intTrial) = matPP_HighHigh(vecTrialTypeIdxHigh(intTrial),intTrial);
		end
		
		%low weights, high act
		matDataPlusLin = [matHighR'; ones(1,size(matHighR',2))];
		matActivation = matWeightsLow'*matDataPlusLin;
		matPP_LowHigh = exp(bsxfun(@minus,matActivation,logsumexp(matActivation,1))); %softmax
		vecConfidenceLowHigh = nan(intTrialNum,1);
		for intTrial=1:intTrialNum
			vecConfidenceLowHigh(intTrial) = matPP_LowHigh(vecTrialTypeIdxHigh(intTrial),intTrial);
		end
		
		%high weights, low act
		matDataPlusLin = [matLowR'; ones(1,size(matLowR',2))];
		matActivation = matWeightsHigh'*matDataPlusLin;
		matPP_HighLow = exp(bsxfun(@minus,matActivation,logsumexp(matActivation,1))); %softmax
		vecConfidenceHighLow = nan(intTrialNum,1);
		for intTrial=1:intTrialNum
			vecConfidenceHighLow(intTrial) = matPP_HighLow(vecTrialTypeIdxLow(intTrial),intTrial);
		end
		
		[h,pPP_HL]=ttest(vecConfidenceHighHigh,vecConfidenceLowLow);
		subplot(2,3,6)
		hold on
		plot([0.5 4.5],[1 1]./numel(vecUnique),'--','color',[0.5 0.5 0.5]);
		errorbar(1,mean(vecConfidenceLowLow),std(vecConfidenceLowLow)./sqrt(numel(vecConfidenceLowLow)),'xb');
		errorbar(2,mean(vecConfidenceLowHigh),std(vecConfidenceLowHigh)./sqrt(numel(vecConfidenceLowHigh)),'x','color',[0.2 0 0.8]);
		errorbar(3,mean(vecConfidenceHighLow),std(vecConfidenceHighLow)./sqrt(numel(vecConfidenceHighLow)),'x','color',[0.8 0 0.2]);
		errorbar(4,mean(vecConfidenceHighHigh),std(vecConfidenceHighHigh)./sqrt(numel(vecConfidenceHighHigh)),'xr');
		hold off
		ylabel('Decoder prob. of correct ori');
		set(gca,'xtick',[1 2 3 4],'xticklabel',{'Low W-low A','Low W-high A','High W-low A','High W-high A'});
		title(sprintf('Train+test on low+low vs hi+hi act, t-test,p=%.2e',pPP_HL));
		xtickangle(15);
		fixfig;
		
		%%
		export_fig(fullpath(strFigurePath,sprintf('B4_LRDecConfidence_Qsplit_OriCodingT%s_%s.tif',num2str(dblStartT),strRec)));
		export_fig(fullpath(strFigurePath,sprintf('B4_LRDecConfidence_Qsplit_OriCodingT%s_%s.pdf',num2str(dblStartT),strRec)));
		
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
			export_fig(fullpath(strFigurePath,sprintf('B5_%s_QuantileDecodingT%s_%s.tif',strDec,num2str(dblStartT),strRec)));
			export_fig(fullpath(strFigurePath,sprintf('B5_%s_QuantileDecodingT%s_%s.pdf',strDec,num2str(dblStartT),strRec)));
		end
	end
end
toc
