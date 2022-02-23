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
strFigurePath = 'F:\Data\Results\PopTimeCoding';

%% select all neurons in LP and drifting grating stimuli
if ~exist('sAggStim','var') || isempty(sAggStim)
	[sAggStim,sAggNeuron]=loadDataNpx('','driftinggrating');
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
		indTuned = sOut.vecFitP<0.05 & sum(matData,2)>(size(matData,2)/dblDur)*dblMinRate;
		
		%prep
		vecPrefOri = rad2deg(sOut.matFittedParams(indTuned,1))/2;
		vecPrefRad = sOut.matFittedParams(indTuned,1);
		intTunedN = sum(indTuned);
		intNumN = intTunedN;
		cellSpikeTimes(~indTuned)=[];
		
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
		cellSpikeTimesPerCellPerTrial_Shuffled = cell(intNumN,intTrialNum);
		for intN=1:intNumN
			% build pseudo data, stitching stimulus periods
			[vecPseudoSpikeTimes,vecPseudoEventT] = getPseudoSpikeVectors(cellSpikeTimes{intN},vecStimOnTime-dblPreTime,dblMaxDur);
			cellSpikeTimesStitched{intN} = vecPseudoSpikeTimes;
			
			vecISI = diff(vecPseudoSpikeTimes);
			vecGenSpikes = cumsum([mean(vecISI);vecISI(randperm(numel(vecISI)))]);
			[vecTrialPerSpikeS,vecTimePerSpikeS] = getSpikesInTrial(vecGenSpikes,vecPseudoEventT,dblMaxDur);
			
			%real
			[vecTrialPerSpike,vecTimePerSpike] = getSpikesInTrial(vecPseudoSpikeTimes,vecPseudoEventT,dblMaxDur);
			for intTrial=1:intTrialNum
				cellSpikeTimesPerCellPerTrial{intN,intTrial} = vecTimePerSpike(vecTrialPerSpike==intTrial);
				cellSpikeTimesPerCellPerTrial_Shuffled{intN,intTrial} = vecTimePerSpikeS(vecTrialPerSpikeS==intTrial);
			end
		end
		vecStimOnStitched = vecPseudoEventT;
		
		%%
		boolPlot = true;
		vecRperTrial = nan(intTrialNum,1);
		vecSperTrial = nan(intTrialNum,1);
		vecHperTrial = nan(intTrialNum,1);
		vecLperTrial = nan(intTrialNum,1);
		vecMperTrial = nan(intTrialNum,1);
		cellIFR_perTrial = cell(intTrialNum,1);
		cellTimeIFR_perTrial = cell(intTrialNum,1);
			
		%neuron shuffled
		vecSperTrial_SN = nan(intTrialNum,1);
		vecHperTrial_SN = nan(intTrialNum,1);
		vecLperTrial_SN = nan(intTrialNum,1);
		vecMperTrial_SN = nan(intTrialNum,1);
		cellIFR_perTrial_Shuffled = cell(intTrialNum,1);
		cellTimeIFR_perTrial_Shuffled = cell(intTrialNum,1);
		
		%pop shuffled
		vecSperTrial_S = nan(intTrialNum,1);
		vecHperTrial_S = nan(intTrialNum,1);
		vecLperTrial_S = nan(intTrialNum,1);
		vecMperTrial_S = nan(intTrialNum,1);
		
		
		% & neuronneuron shuffled
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
			[vecTimeIFR,vecIFR] = getIFR(vecAllSpikes,0,dblMaxDur,[],[],[],0);
			vecR_sorted = sort(vecIFR);
			int5perc = round(numel(vecR_sorted)/20);
			vecHperTrial(intTrial) = mean(vecR_sorted((1+end-1*int5perc):(end-0*int5perc)));
			vecLperTrial(intTrial) = mean(vecR_sorted((1+0*int5perc):(1*int5perc)));
			vecMperTrial(intTrial) = mean(vecR_sorted);
			vecSperTrial(intTrial) = std(vecR_sorted);
			cellIFR_perTrial{intTrial} = vecIFR;
			cellTimeIFR_perTrial{intTrial} = vecTimeIFR;
			
			%ISI
			vecISI1 = abs([vecAllSpikes(1:(end-1)) - vecAllSpikes(2:end); inf]);
			vecISI2 = abs([inf; vecAllSpikes(2:end) - vecAllSpikes(1:(end-1))]);
			vecISI = vecISI1(1:(end-1));%min([vecISI1 vecISI2],[],2);
			vecNSI = min([vecISI1 vecISI2],[],2); %nearest spike interval
			
			%shuffled single-neuron ISIs
			vecAllSpikesShuff = sort(cell2vec(cellSpikeTimesPerCellPerTrial_Shuffled(:,intTrial)));
			vecISI0 = [vecAllSpikesShuff(2:end) - vecAllSpikesShuff(1:(end-1)); inf];
			vecAllSpikesShuff(vecISI0==0)=vecAllSpikesShuff(vecISI0==0)-(10^-5)*rand();
			vecAllSpikesShuff = uniquetol(vecAllSpikesShuff,1e-7);
			vecD = diff(vecAllSpikesShuff);
			vecD1 = vecD(1:(end-1));
			vecD2 = vecD(2:end);
			[r,p]=corr(vecD1,vecD2);
			vecRperTrial(intTrial) = r;
			[vecTimeIFRS,vecIFRS] = getIFR(vecAllSpikesShuff,0,dblMaxDur,[],[],[],0);
			vecR_sortedS = sort(vecIFRS);
			int5perc = round(numel(vecR_sortedS)/20);
			vecHperTrial_SN(intTrial) = mean(vecR_sortedS((1+end-1*int5perc):(end-0*int5perc)));
			vecLperTrial_SN(intTrial) = mean(vecR_sortedS((1+0*int5perc):(1*int5perc)));
			vecMperTrial_SN(intTrial) = mean(vecR_sortedS);
			vecSperTrial_SN(intTrial) = std(vecR_sortedS);
			cellIFR_perTrial_Shuffled{intTrial} = vecIFRS;
			cellTimeIFR_perTrial_Shuffled{intTrial} = vecTimeIFRS;
			vecISIS = diff(vecAllSpikesShuff);
			
			%shuffle to exponential distribution at population level
			intIters = 1;
			for intIter=1%:intIters
				vecGenSpikes = cumsum([0;vecISI(randperm(numel(vecISI)))])+vecAllSpikes(1);
				[vecTimeIFRS,vecIFRS] = getIFR(vecGenSpikes,0,dblMaxDur,[],[],[],0);
				
				%error add expected low/high values from shuffling
				vecR_sorted = sort(vecIFRS);
				int5perc = round(numel(vecR_sorted)/20);
				vecHperTrial_S(intTrial) = mean(vecR_sorted((1+end-1*int5perc):(end-0*int5perc)));
				vecLperTrial_S(intTrial) = mean(vecR_sorted((1+0*int5perc):(1*int5perc)));
				vecMperTrial_S(intTrial) = mean(vecR_sorted);
				vecSperTrial_S(intTrial) = std(vecR_sorted);
				
				%pop & neuron shuffled
				vecGenSpikesS = cumsum([0;vecISIS(randperm(numel(vecISIS)))])+vecAllSpikesShuff(1);
				[vecTimeIFRSS,vecIFRSS] = getIFR(vecGenSpikesS,0,dblMaxDur,[],[],[],0);
				
				%error add expected low/high values from shuffling
				vecRS_sorted = sort(vecIFRSS);
				int5percS = round(numel(vecRS_sorted)/20);
				vecHperTrial_SS(intTrial) = mean(vecRS_sorted((1+end-1*int5percS):(end-0*int5percS)));
				vecLperTrial_SS(intTrial) = mean(vecRS_sorted((1+0*int5percS):(1*int5percS)));
				vecMperTrial_SS(intTrial) = mean(vecRS_sorted);
				vecSperTrial_SS(intTrial) = std(vecRS_sorted);
			end
			
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
				clf;maxfig;
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
		
		scatter(cat(1,vecMperTrial,vecMperTrial),cat(1,vecLperTrial_S(:),vecHperTrial_S(:)),[],[0.7 0.7 0.7],'.');
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
		legend({'Shuffled pop ISIs','Lowest 5%','Highest 5%'},'location','best');
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
		matRespS = cellfun(@numel,cellSpikeTimesPerCellPerTrial_Shuffled);
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
		
		scatter(cat(1,vecMperTrial_SN,vecMperTrial_SN),cat(1,vecLperTrial_SS(:),vecHperTrial_SS(:)),[],[0.7 0.7 0.7],'.');
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
		legend({'Shuffled','Lowest 5%','Highest 5%'},'location','best');
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
		
		%% what is circular variance of preferred oris of cells that are active during low-5% and high 5% epochs?
		%are only cells tuned to the orientation active during low phases? Is this different from high 5%?
		
		%pre-alloc
		intQuantileNum = 5; %5=20%
		
		%run
		for intShuff=[0 1]
			%pre-allocate
			matAggR_temp = nan(3,intTrialNum,intTunedN);%save middle60/lower20/upper20
			vecCircPrecLow_temp = nan(1,intTrialNum);
			vecCircPrecHigh_temp = nan(1,intTrialNum);
			
			for intTrial=1:intTrialNum
				%% define data
				if intShuff == 0
					%get IFRs
					vecTrialIFR = cellIFR_perTrial{intTrial}(2:(end-1));
					vecTrialTimeIFR = cellTimeIFR_perTrial{intTrial}(2:(end-1));
					cellSpikes = cellSpikeTimesPerCellPerTrial(:,intTrial);
				else
					%get IFRs
					vecTrialIFR = cellIFR_perTrial_Shuffled{intTrial}(2:(end-1));
					vecTrialTimeIFR = cellTimeIFR_perTrial_Shuffled{intTrial}(2:(end-1));
					cellSpikes = cellSpikeTimesPerCellPerTrial_Shuffled(:,intTrial);
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
				vecChanges = find(diff(vecHighLowIdx)~=1);
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
				
				%% calc circ prec
				vecSumLow = squeeze(matAggR_temp(2,intTrial,:));
				if sum(vecSumLow) > 1
					dblCircVarLow = circ_var(vecPrefRad(vecSumLow>0), vecSumLow(vecSumLow>0));
				else
					dblCircVarLow = nan;
				end
				
				vecSumHigh = squeeze(matAggR_temp(3,intTrial,:));
				if sum(vecSumHigh) > 1
					dblCircVarHigh = circ_var(vecPrefRad(vecSumHigh>0), vecSumHigh(vecSumHigh>0));
				else
					dblCircVarHigh = nan;
				end
				
				vecCircPrecLow_temp(intTrial) = 1-dblCircVarLow;
				vecCircPrecHigh_temp(intTrial) = 1-dblCircVarHigh;
				
			end
			
			%% save
			if intShuff == 0
				%get IFRs
				vecCircPrecLow = vecCircPrecLow_temp;
				vecCircPrecHigh = vecCircPrecHigh_temp;
				matAggR = matAggR_temp;
			else
				%get IFRs
				vecCircPrecLow_shuff = vecCircPrecLow_temp;
				vecCircPrecHigh_shuff = vecCircPrecHigh_temp;
				matAggR_shuff = matAggR_temp;
			end
		end
		
		%% mean per trial per quantile
		matLowR = squeeze(matAggR(2,:,:));
		matHighR = squeeze(matAggR(3,:,:));
		vecOriLow = vecOri180;
		vecOriHigh = vecOri180;
		
		%%
		figure;maxfig;
		subplot(2,3,1)
		plot([0 1],[0 1],'k--')
		hold on
		scatter(vecCircPrecLow,vecCircPrecLow_shuff,[],lines(1),'.')
		hold off
		title(sprintf('Lowest 20%%, real=%.3f, shuff=%.3f',mean(vecCircPrecLow),mean(vecCircPrecLow_shuff)))
		xlabel('Real orientation selectivity (OSI)');
		ylabel('Shuffled orientation selectivity (OSI)');
		fixfig;
		
		subplot(2,3,2)
		plot([0 1],[0 1],'k--')
		hold on
		scatter(vecCircPrecHigh,vecCircPrecHigh_shuff,[],lines(1),'.')
		hold off
		title(sprintf('Highest 20%%, real=%.3f, shuff=%.3f',mean(vecCircPrecHigh),mean(vecCircPrecHigh_shuff)))
		xlabel('Real orientation selectivity (OSI)');
		ylabel('Shuffled orientation selectivity (OSI)');
		fixfig;
		
		subplot(2,3,3)
		vecEffect2 = vecCircPrecLow-vecCircPrecLow_shuff;
		vecEffect4 = vecCircPrecHigh-vecCircPrecHigh_shuff;
		plot([-1 1]*0.5,0.5*[-1 1],'k--')
		hold on
		scatter(vecEffect2,vecEffect4,[],lines(1),'.')
		hold off
		xlabel('dOSI low q; d(real,shuffled)');
		ylabel('dOSI high q; d(real,shuffled)');
		fixfig;
		
		
		subplot(2,3,4)
		dblStepCV = 0.05;
		vecBinECV = -0.6:dblStepCV:0.6;
		vecBinCCV = vecBinECV(2:end)-dblStepCV/2;
		vecCountsLow = histcounts(vecEffect2,vecBinECV);
		plot(vecBinCCV,vecCountsLow);
		[h,pLow]=ttest(vecEffect2);
		title(sprintf('mean precision low q=%.3f, p=%.1e',mean(vecEffect2),pLow))
		xlabel('dOSI low q; d(real,shuffled)');
		ylabel('Number of trials (count)')
		fixfig;
		
		subplot(2,3,5)
		vecCountsHigh = histcounts(vecEffect4,vecBinECV);
		plot(vecBinCCV,vecCountsHigh);
		[h,pHigh]=ttest(vecEffect4);
		title(sprintf('mean precision high q=%.3f, p=%.1e',mean(vecEffect4),pHigh))
		xlabel('dOSI high q; d(real,shuffled)');
		ylabel('Number of trials (count)')
		fixfig;
		
		subplot(2,3,6)
		vecCountsDiff = histcounts(vecEffect2-vecEffect4,vecBinECV);
		plot(vecBinCCV,vecCountsDiff);
		[h,pDiff]=ttest(vecEffect2-vecEffect4);
		title(sprintf('d(High,Low)=%.3f, p=%.1e',mean(vecEffect2-vecEffect4),pDiff))
		xlabel('Difference dOSI d(low,high)');
		ylabel('Number of trials (count)')
		fixfig;
		
		export_fig(fullpath(strFigurePath,sprintf('B4_QuantileTuningT%s_%s.tif',num2str(dblStartT),strRec)));
		export_fig(fullpath(strFigurePath,sprintf('B4_QuantileTuningT%s_%s.pdf',num2str(dblStartT),strRec)));
		
		%%
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
		
		vecDiffHighLow = vecEffect2-vecEffect4;
		[r,p]=corr(vecEffect2',vecConfidenceLow);
		[r2,p2]=corr(vecEffect4',vecConfidenceHigh);
		
		%% ori tuning is stable across quantiles
		sOut = getTuningCurves(matLowR',vecOri180,0);
		vecPrefRad;
		vecPrefRadLow = sOut.matFittedParams(:,1);
		sOut = getTuningCurves(matHighR',vecOri180,0);
		vecPrefRadHigh = sOut.matFittedParams(:,1);
		%ori tuning is stable
		
		%% so why is there no effect on the decoder?
		%=> compare greedy decoder of high vs low; greedy in neurons or spikes?
		for intDecoder=1:4
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
				indNeuronsUsed = false(1,intTunedN);
				vecGreedyPerf = nan(1,intTunedN);
				vecSelectNeuronOrder = nan(1,intTunedN);
				for intNeuronsUsed=1:intTunedN
					intUnusedNeurons = intTunedN-sum(indNeuronsUsed);
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
					fprintf('Decoded %s, neuron %d/%d; current best perf is %.3f; selected neuron %d [%s]\n',strHighLow,intNeuronsUsed,intTunedN,dblMaxPerf,intUseN,getTime);
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
		for intUseDec = 1:4
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
			vecX = 1:intTunedN;
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
