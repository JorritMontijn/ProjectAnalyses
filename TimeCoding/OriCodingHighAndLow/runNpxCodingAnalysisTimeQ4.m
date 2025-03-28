%% aim
%{
https://www.nature.com/articles/s41598-021-95037-z
https://www.nature.com/articles/s42003-021-02437-y

q4; compare proportion of tuned cells & interneurons during low/high epochs

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
	[sAggStim,sAggNeuron,sAggSources]=loadDataNpx('','driftinggrating',strDataPath);
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
	
	%% is cell an interneuron (fast/narrow spiking) or pyramid (regular/broad spiking)?
	%load waveform
	[sThisRec,sUseNeuron] = loadWaveforms(sThisRec,sUseNeuron);
	
	%calculate waveform properties
	dblSampRateIM = sThisRec.sample_rate;
	dblRecDur = max(cellfun(@max,{sUseNeuron.SpikeTimes})) - min(cellfun(@min,{sUseNeuron.SpikeTimes}));
	vecSpikeRate = cellfun(@numel,{sUseNeuron.SpikeTimes})/dblRecDur;
	matAreaWaveforms = cell2mat({sUseNeuron.Waveform}'); %[cell x sample]
	intNeurons=size(matAreaWaveforms,1);
	vecSpikeDur = nan(1,intNeurons);
	vecSpikePTR = nan(1,intNeurons);
	for intNeuron=1:intNeurons
		%find trough
		[dblTroughVal,intTrough]=min(matAreaWaveforms(intNeuron,:));
		[dblPeakVal,intTroughToPeak]=max(matAreaWaveforms(intNeuron,intTrough:end));
		intPeak = intTrough + intTroughToPeak - 1;
		
		dblTroughTime = intTrough/dblSampRateIM;
		dblTroughToPeakTime = intTroughToPeak/dblSampRateIM;
		dblPeakTime = intPeak/dblSampRateIM;
		vecSpikeDur(intNeuron) = dblTroughToPeakTime;
		vecSpikePTR(intNeuron) = abs(dblPeakVal/dblTroughVal);
	end
	dblPTT = 0.5;
	dblSWT = 0.5/1000;
	vecNarrow = vecSpikeDur < dblSWT & vecSpikePTR > dblPTT; %Ctx BL6
	vecBroad = vecSpikeDur > dblSWT & vecSpikePTR < dblPTT; %Ctx BL6
	vecOther = ~vecNarrow & ~vecBroad;
	vecCol = vecNarrow + vecBroad*2 + vecOther*3;
	%scatter(vecSpikeDur,vecSpikePTR,[],vecCol)
	
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
		intTunedN = sum(indTuned);
		intNumN = sum(indResp);
		cellSpikeTimes(~indResp)=[];
		vecNarrow(~indResp)=[];
		vecBroad(~indResp)=[];
		vecOther(~indResp)=[];
		indTuned(~indResp)=[];
		
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
		
		%% calc IFRs per trial
		cellIFR_perTrial = cell(intTrialNum,1);
		cellTimeIFR_perTrial = cell(intTrialNum,1);
		cellIFR_perTrial_S = cell(intTrialNum,1);
		cellTimeIFR_perTrial_S = cell(intTrialNum,1);
		for intTrial=1:intTrialNum
			%real
			vecAllSpikes = sort(cell2vec(cellSpikeTimesPerCellPerTrial(:,intTrial)));
			vecISI0 = [vecAllSpikes(2:end) - vecAllSpikes(1:(end-1)); inf];
			vecAllSpikes(vecISI0==0)=vecAllSpikes(vecISI0==0)-(10^-5)*rand();
			vecAllSpikes = uniquetol(vecAllSpikes,1e-7);
			[vecTimeIFR,vecIFR] = getIFR(vecAllSpikes,0,dblMaxDur,[],[],[],0);
			cellIFR_perTrial{intTrial} = vecIFR;
			cellTimeIFR_perTrial{intTrial} = vecTimeIFR;
			vecIFR = cellIFR_perTrial{intTrial};
			vecTimeIFR = cellTimeIFR_perTrial{intTrial};
			
			%shuffled single-neuron ISIs
			vecAllSpikesShuff = sort(cell2vec(cellSpikeTimesPerCellPerTrial_S(:,intTrial)));
			vecISI0 = [vecAllSpikesShuff(2:end) - vecAllSpikesShuff(1:(end-1)); inf];
			vecAllSpikesShuff(vecISI0==0)=vecAllSpikesShuff(vecISI0==0)-(10^-5)*rand();
			vecAllSpikesShuff = uniquetol(vecAllSpikesShuff,1e-7);
			[vecTimeIFRS,vecIFRS] = getIFR(vecAllSpikesShuff,0,dblMaxDur,[],[],[],0);
			cellIFR_perTrial_S{intTrial} = vecIFRS;
			cellTimeIFR_perTrial_S{intTrial} = vecTimeIFRS;
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
			int5perc = round(numel(vecR_sorted)/20);
			vecHperTrial(intTrial) = mean(vecR_sorted((1+end-1*int5perc):(end-0*int5perc)));
			vecLperTrial(intTrial) = mean(vecR_sorted((1+0*int5perc):(1*int5perc)));
			vecMperTrial(intTrial) = mean(vecR_sorted);
			vecSperTrial(intTrial) = std(vecR_sorted);
			
			%ISI
			vecISI1 = abs([vecAllSpikes(1:(end-1)) - vecAllSpikes(2:end); inf]);
			vecISI2 = abs([inf; vecAllSpikes(2:end) - vecAllSpikes(1:(end-1))]);
			vecISI = vecISI1(1:(end-1));%min([vecISI1 vecISI2],[],2);
			vecNSI = min([vecISI1 vecISI2],[],2); %nearest spike interval
			
			%shuffled single-neuron ISIs
			vecIFRS = cellIFR_perTrial_S{intTrial};
			vecTimeIFRS = cellTimeIFR_perTrial_S{intTrial};
			
			%error add expected low/high values from shuffling
			vecR_sorted = sort(vecIFRS);
			int5perc = round(numel(vecR_sorted)/20);
			vecHperTrial_S(intTrial) = mean(vecR_sorted((1+end-1*int5perc):(end-0*int5perc)));
			vecLperTrial_S(intTrial) = mean(vecR_sorted((1+0*int5perc):(1*int5perc)));
			vecMperTrial_S(intTrial) = mean(vecR_sorted);
			vecSperTrial_S(intTrial) = std(vecR_sorted);
			
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
				
				export_fig(fullpath(strFigurePath,sprintf('D1_ExampleActivityT%s_%sTrial%d.tif',num2str(dblStartT),strRec,intTrial)));
				export_fig(fullpath(strFigurePath,sprintf('D1_ExampleActivityT%s_%sTrial%d.pdf',num2str(dblStartT),strRec,intTrial)));
				boolPlot = false;
			end
		end
		
		%% wq4; compare proportion of tuned cells & interneurons during low/high epochs
		%are only cells tuned to the orientation active during low phases? Is this different from high 5%?
		
		%pre-alloc
		intQuantileNum = 5; %5=20%
		vecTuning = sOut.vecFitR2(indResp);
		
		%run
		for intShuff=[0 1]
			%pre-allocate
			matAggR_temp = nan(3,intTrialNum,intNumN);%save middle60/lower20/upper20
			vecRatioNBHigh_temp = nan(1,intTrialNum);
			vecRatioNBLow_temp = nan(1,intTrialNum);
			vecTunedLow_temp = nan(1,intTrialNum);
			vecTunedHigh_temp = nan(1,intTrialNum);
			for intTrial=1:intTrialNum
				%% define data
				%fprintf('Running shuff=%d; trial %d/%d [%s]\n',intShuff,intTrial,intTrialNum,getTime);
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
				for intNeuron=1:intNumN
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
				
				%% calc narrow/broad ratio
				vecSumLow = squeeze(matAggR_temp(2,intTrial,:));
				vecSumHigh = squeeze(matAggR_temp(3,intTrial,:));
				
				%save
				vecRatioNBLow_temp(intTrial) = sum(vecNarrow'.*vecSumLow) / (sum(vecNarrow'.*vecSumLow) + sum(vecBroad'.*vecSumLow));
				vecRatioNBHigh_temp(intTrial) = sum(vecNarrow'.*vecSumHigh) / (sum(vecNarrow'.*vecSumHigh) + sum(vecBroad'.*vecSumHigh));
				vecTunedLow_temp(intTrial) = sum(indTuned.*vecSumLow) / (sum(indTuned.*vecSumLow) + sum((~indTuned).*vecSumLow));
				vecTunedHigh_temp(intTrial) = sum(indTuned.*vecSumHigh) / (sum(indTuned.*vecSumHigh) + sum((~indTuned).*vecSumHigh));
			end
			
			%% save
			if intShuff == 0
				%get IFRs
				vecRatioNBLow = vecRatioNBLow_temp;
				vecRatioNBHigh = vecRatioNBHigh_temp;
				vecTunedLow = vecTunedLow_temp;
				vecTunedHigh = vecTunedHigh_temp;
				matAggR = matAggR_temp;
			else
				%get IFRs
				vecRatioNBLow_shuff = vecRatioNBLow_temp;
				vecRatioNBHigh_shuff = vecRatioNBHigh_temp;
				vecTunedLow_shuff = vecTunedLow_temp;
				vecTunedHigh_shuff = vecTunedHigh_temp;
				matAggR_shuff = matAggR_temp;
			end
		end
		
		%% mean per trial per quantile
		matLowR = squeeze(matAggR(2,:,:));
		matHighR = squeeze(matAggR(3,:,:));
		vecOriLow = vecOri180;
		vecOriHigh = vecOri180;
		
		%% tuning R2
		vecColH = [0.8 0 0];
		dblQuantileSize = round(100/intQuantileNum);
		dblMaxVal = 0.08;
		figure;maxfig;
		subplot(2,3,1)
		%plot(dblMaxVal*[0 1],dblMaxVal*[0 1],'k--')
		hold on
		scatter(vecRatioNBLow_shuff,vecRatioNBLow,[],lines(1),'.')
		hold off
		title(sprintf('Lowest %d%%, real=%.3f, shuff=%.3f',dblQuantileSize,mean(vecRatioNBLow),mean(vecRatioNBLow_shuff)))
		ylabel('Real narrow/broad ratio)');
		xlabel('Shuffled narrow/broad ratio');
		fixfig;
		
		subplot(2,3,2)
		%plot(dblMaxVal*[0 1],dblMaxVal*[0 1],'k--')
		hold on
		scatter(vecRatioNBHigh_shuff,vecRatioNBHigh,[],vecColH,'.')
		hold off
		title(sprintf('Highest %d%%, real=%.3f, shuff=%.3f',dblQuantileSize,mean(vecRatioNBHigh),mean(vecRatioNBHigh_shuff)))
		ylabel('Real narrow/broad ratio');
		xlabel('Shuffled narrow/broad ratio');
		fixfig;
		
		subplot(2,3,3)
		vecEffectLow_R2 = (vecRatioNBLow-vecRatioNBLow_shuff)./(vecRatioNBLow + vecRatioNBLow_shuff);
		vecEffectHigh_R2 = (vecRatioNBHigh-vecRatioNBHigh_shuff)./(vecRatioNBHigh + vecRatioNBHigh_shuff);
		%plot(dblMaxVal*[0 1],dblMaxVal*[0 1],'k--')
		hold on
		scatter(vecEffectLow_R2,vecEffectHigh_R2,[],lines(1),'.')
		hold off
		xlabel('dMean narrow/broad ratio low q');
		ylabel('dMean narrow/broad ratio high q');
		fixfig;
		
		
		subplot(2,3,4)
		dblStepCV = 0.05;
		vecBinECV = -0.5:dblStepCV:0.5;
		vecBinCCV = vecBinECV(2:end)-dblStepCV/2;
		vecCountsLow = histcounts(vecEffectLow_R2,vecBinECV);
		plot(vecBinCCV,vecCountsLow);
		[h,pLow]=ttest(vecEffectLow_R2);
		title(sprintf('mean narrow/broad ratio low q=%.3f, p=%.1e',mean(vecEffectLow_R2),pLow))
		xlabel('d(narrow/broad ratio) low q; d(real,shuffled)');
		ylabel('Number of trials (count)')
		fixfig;
		
		subplot(2,3,5)
		vecCountsHigh = histcounts(vecEffectHigh_R2,vecBinECV);
		plot(vecBinCCV,vecCountsHigh,'color',vecColH);
		[h,pHigh]=ttest(vecEffectHigh_R2);
		title(sprintf('mean narrow/broad ratio high q=%.3f, p=%.1e',mean(vecEffectHigh_R2),pHigh))
		xlabel('d(narrow/broad ratio) high q; d(real,shuffled)');
		ylabel('Number of trials (count)')
		fixfig;
		
		subplot(2,3,6)
		vecLocX = [0.4 0.6];
		errorbar(vecLocX(1),mean(vecEffectLow_R2),std(vecEffectLow_R2)./sqrt(intTrialNum),'x','color',lines(1));
		hold on
		errorbar(vecLocX(2),mean(vecEffectHigh_R2),std(vecEffectHigh_R2)./sqrt(intTrialNum),'x','color',vecColH);
		hold off
		[h,pDiff]=ttest(vecEffectLow_R2,vecEffectHigh_R2);
		title(sprintf('Diff, t-test, p=%.2e',pDiff));
		xlim([0.3 0.7]);
		set(gca,'xtick',vecLocX,'xticklabel',{'Low Q','High Q'});
		ylabel('NBR increase over shuffled')
		fixfig;
		
		
		export_fig(fullpath(strFigurePath,sprintf('D2_QuantileNBRT%s_%s.tif',num2str(dblStartT),strRec)));
		export_fig(fullpath(strFigurePath,sprintf('D2_QuantileNBRT%s_%s.pdf',num2str(dblStartT),strRec)));
		
		%% circ var
		dblQuantileSize = round(100/intQuantileNum);
		dblMaxVal = 0.08;
		figure;maxfig;
		subplot(2,3,1)
		%plot(dblMaxVal*[0 1],dblMaxVal*[0 1],'k--')
		hold on
		scatter(vecTunedLow_shuff,vecTunedLow,[],lines(1),'.')
		hold off
		title(sprintf('Lowest %d%%, real=%.3f, shuff=%.3f',dblQuantileSize,mean(vecTunedLow),mean(vecTunedLow_shuff)))
		ylabel('Real tuned fraction');
		xlabel('Shuffled tuned fraction');
		fixfig;
		
		subplot(2,3,2)
		%plot(dblMaxVal*[0 1],dblMaxVal*[0 1],'k--')
		hold on
		scatter(vecTunedHigh_shuff,vecTunedHigh,[],vecColH,'.')
		hold off
		title(sprintf('Highest %d%%, real=%.3f, shuff=%.3f',dblQuantileSize,mean(vecTunedHigh),mean(vecTunedHigh_shuff)))
		ylabel('Real tuned fraction');
		xlabel('Shuffled tuned fraction');
		fixfig;
		
		subplot(2,3,3)
		vecEffectLow_Tuning = (vecTunedLow-vecTunedLow_shuff)./(vecTunedLow+vecTunedLow_shuff);
		vecEffectHigh_Tuning = (vecTunedHigh-vecTunedHigh_shuff)./(vecTunedHigh+vecTunedHigh_shuff);
		%plot(dblMaxVal*[0 1],dblMaxVal*[0 1],'k--')
		hold on
		scatter(vecEffectLow_Tuning,vecEffectHigh_Tuning,[],lines(1),'.')
		hold off
		xlabel('dMean tuned fraction low q');
		ylabel('dMean tuned fraction high q');
		fixfig;
		
		subplot(2,3,4)
		dblStepCV = 0.05;
		vecBinECV = -0.5:dblStepCV:0.5;
		vecBinCCV = vecBinECV(2:end)-dblStepCV/2;
		vecCountsLow = histcounts(vecEffectLow_Tuning,vecBinECV);
		plot(vecBinCCV,vecCountsLow);
		[h,pLow]=ttest(vecEffectLow_Tuning);
		title(sprintf('mean tuned fraction low q=%.3f, p=%.1e',mean(vecEffectLow_Tuning),pLow))
		xlabel('d(tuned fraction) low q; d(real,shuffled)');
		ylabel('Number of trials (count)')
		fixfig;
		
		subplot(2,3,5)
		vecCountsHigh = histcounts(vecEffectHigh_Tuning,vecBinECV);
		plot(vecBinCCV,vecCountsHigh,'color',vecColH);
		[h,pHigh]=ttest(vecEffectHigh_Tuning);
		title(sprintf('mean tuned fraction high q=%.3f, p=%.1e',mean(vecEffectHigh_Tuning),pHigh))
		xlabel('d(tuned fraction) high q; d(real,shuffled)');
		ylabel('Number of trials (count)')
		fixfig;
		
		subplot(2,3,6)
		vecLocX = [0.4 0.6];
		errorbar(vecLocX(1),mean(vecEffectLow_Tuning),std(vecEffectLow_Tuning)./sqrt(intTrialNum),'x','color',lines(1));
		hold on
		errorbar(vecLocX(2),mean(vecEffectHigh_Tuning),std(vecEffectHigh_Tuning)./sqrt(intTrialNum),'x','color',vecColH);
		hold off
		xlim([0.3 0.7]);
		set(gca,'xtick',vecLocX,'xticklabel',{'Low Q','High Q'});
		[h,pDiff]=ttest(vecEffectLow_Tuning,vecEffectHigh_Tuning);
		title(sprintf('Diff, t-test, p=%.2e',pDiff));
		ylabel('Tuned fraction increase over shuffled')
		fixfig;
		
		export_fig(fullpath(strFigurePath,sprintf('D3_QuantileTunedT%s_%s.tif',num2str(dblStartT),strRec)));
		export_fig(fullpath(strFigurePath,sprintf('D3_QuantileTunedT%s_%s.pdf',num2str(dblStartT),strRec)));
		
	end
end
toc
