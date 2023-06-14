%% aim
%{
https://www.nature.com/articles/s41598-021-95037-z
https://www.nature.com/articles/s42003-021-02437-y

q1: Can absence of stimulus information in initial period be explained by neurons responding only to
black/white patches in the first x ms?

q2: How precise are spike times in the natural movie repetitions?


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
		[matMeanRate,cellSpikeTimesReal] = ...
			NpxPrepMovieData(cellSpikeTimesRaw,vecStimOnTime,vecStimOffTime,vecFrameIdx);
		
		%% load prepro T0 data
		% pool spikes from all neurons, but save the time+id per spike, then calculate IFR over all
		% spikes at pop level, detect peaks, and split into trials.
		%get spikes per trial per neuron
		sLoad = load(fullpath(strTargetDataPath,sprintf('T0Data_%s',strRec)));
		
		%% detect peaks
		% filter
		%real
		intLag = round((numel(sLoad.sReal.vecIFR)/numel(vecOrigStimOnTime))/2);
		if (intLag/2) == round(intLag/2)
			intLag = intLag - 1;
		end
		dblThreshZ = 1;
		dblInfluence = 0.5;
		
		%% transform time indices
		%events
		dblStartEpoch = vecOrigStimOnTime(1)-dblStimDur;
		dblEpochDur = vecOrigStimOnTime(end)-vecOrigStimOnTime(1)+dblStimDur;
		dblStopEpoch = dblStartEpoch + dblEpochDur;
		
		%% plot
		sReal = struct;
		sShuff = struct;
		sShuffTid = struct;
		sPoiss = struct;
		cellTypes = {'Real','Poiss','ShuffTid','Shuff','PoissGain'};
		for intType=5
			%which type?
			cellSpikeTimes = cell(1,intNumN);
			strType = cellTypes{intType};
			if strcmp(strType,'Real')
				sSource = sLoad.sReal;
			elseif strcmp(strType,'Poiss')
				sSource = sLoad.sPoiss;
			elseif strcmp(strType,'ShuffTid')
				sSource = sLoad.sShuffTid;
			elseif strcmp(strType,'Shuff')
				sSource = sLoad.sShuff;
			elseif strcmp(strType,'PoissGain')
			
			end
			vecTime = sSource.vecTime;
			vecIFR = sSource.vecIFR;
			vecAllSpikeTime = sSource.vecAllSpikeTime;
			vecAllSpikeNeuron = sSource.vecAllSpikeNeuron;
			intNumN = max(vecAllSpikeNeuron);
			
			%% replace IFR with binned rates?
			if 0
				dblBinDur = (5/1000);
				vecBins=(dblStartEpoch:dblBinDur:dblStopEpoch)';
				vecIFR = flat(histcounts(vecTime,vecBins)./dblBinDur);
				vecTime = vecBins(2:end)-dblBinDur/2;
			end
			
			%% run analyses
			%remove spikes outside epoch & build shuffled array
			cellUseSpikeTimes = cell(1,intNumN);
			for intN=1:intNumN
				vecT = vecAllSpikeTime(vecAllSpikeNeuron==intN);
				indRem = (vecT > dblStopEpoch) | (vecT < dblStartEpoch);
				cellUseSpikeTimes{intN} = vecT(~indRem);
			end
			
			%mean-subtraction
			[signals,avgFilter,stdFilter] = detectpeaks(vecIFR,intLag,dblThreshZ,dblInfluence);
			%vecIFR = (vecIFR - avgFilter)./stdFilter;
			vecNormIFR = (vecIFR - avgFilter)./avgFilter;
			indRem = isinf(vecNormIFR);
			vecNormIFR(indRem) = [];
			vecNormIFR = vecNormIFR + (1e-10)*rand(size(vecNormIFR)); %break identical values
			vecNormTime = vecTime(~indRem);
			
			% peaks
			indStimSpikes = vecNormTime>dblStartEpoch & vecNormTime<dblStopEpoch;
			vecStimIFR_Raw = vecIFR(indStimSpikes);
			vecStimIFR = vecNormIFR(indStimSpikes);
			vecStimTime = vecNormTime(indStimSpikes);
			[vecPeakHeight,vecPeakLocs,w,p] = findpeaks(vecStimIFR);
			
			%threshold peaks
			dblCutOff = 0.65;
			vecStimPeakLocs = vecPeakLocs(vecPeakHeight>dblCutOff);
			
			%retain only stim epoch
			indStimEpoch = vecAllSpikeTime > dblStartEpoch & vecAllSpikeTime < (dblStartEpoch + dblEpochDur);
			vecPopSpikeTime = vecAllSpikeTime(indStimEpoch);
			vecPopNeuronId = vecAllSpikeNeuron(indStimEpoch);
			
			%merged peaks
			matPeakDomain = mergepeaks(vecStimTime,vecStimIFR,vecStimPeakLocs);
			vecMergedPeakLocs = matPeakDomain(:,1);
			
			%get merged pop events
			vecMergedPopEventTimes = vecStimTime(vecMergedPeakLocs);
			intMergedPopEventNum = numel(vecMergedPopEventTimes);
			vecMergedPopEventLocs = nan(1,intMergedPopEventNum);
			vecMergedPopEventLocsIFR = nan(1,intMergedPopEventNum);
			parfor intEvent=1:intMergedPopEventNum
				[dummy,vecMergedPopEventLocs(intEvent)] = min(abs(vecMergedPopEventTimes(intEvent)-vecPopSpikeTime));
				[dummy,vecMergedPopEventLocsIFR(intEvent)] = min(abs(vecMergedPopEventTimes(intEvent)-vecStimTime));
			end
			
			%% what does peak look like? plot PSTH
			dblStep = (1/50);
			dblRange = 10;%0.02
			vecRandMergedPopEventTimes = dblRange*(rand(1,numel(vecMergedPopEventTimes))-0.5)+linspace(min(vecMergedPopEventTimes)+2*dblRange,max(vecMergedPopEventTimes)-2*dblRange,numel(vecMergedPopEventTimes));
			vecEventBins = (-dblRange+dblStep/2):dblStep:(dblRange-dblStep/2);
			vecEventBinsC = vecEventBins(2:end)-dblStep/2;
			
			%merged peaks, IFR
			[matPET,vecEventBinsC] = doPET(vecTime,vecIFR,vecMergedPopEventTimes,vecEventBins);
			
			%merged random peaks, IFR
			matPETR = doPET(vecTime,vecIFR,vecRandMergedPopEventTimes,vecEventBins);
			
			%normalized pop rates
			vecMeanM = nanmean(matPET,1);
			vecSEMM = nanstd(matPET,[],1)./sqrt(size(matPET,1));
			vecMeanMR = nanmean(matPETR,1);
			vecSEMMR = nanstd(matPETR,[],1)./sqrt(size(matPETR,1));
			vecMeanM(vecEventBinsC==0)=nan;
			vecMeanZ = (vecMeanM-mean(matPETR(:)))./std(matPETR(:));
			
			%short timescale
			dblStep_Short = (1/2000);
			dblRange_Short = 0.04;%0.02
			vecEventBins_Short = (-dblRange_Short+dblStep_Short/2):dblStep_Short:(dblRange_Short-dblStep_Short/2);
			
			%merged peaks
			[dummy,dummy,vecEventBinsC_Short,matPETM_Short] = doPEP(cellUseSpikeTimes,vecEventBins_Short,vecMergedPopEventTimes);
			vecEventBinsC_Short = vecEventBinsC_Short*1000;%ms
			matPopRateM_Short = sum(matPETM_Short,3);
			vecMeanM_Short = mean(matPopRateM_Short,1);
			vecSEMM_Short = std(matPopRateM_Short,[],1)./sqrt(size(matPopRateM_Short,1));
			vecMeanM_Short(vecEventBinsC_Short==0)=nan;
			
			%merged random peak times
			[dummy,dummy,dummy,matPETMR_Short] = doPEP(cellUseSpikeTimes,vecEventBins_Short,vecRandMergedPopEventTimes);
			matPopRateMR_Short = sum(matPETMR_Short,3);
			vecMeanMR_Short = mean(matPopRateMR_Short,1);
			vecSEMMR_Short = std(matPopRateMR_Short,[],1)./sqrt(size(matPopRateMR_Short,1));
			
			%{
			%spike count rates
			[dummy,dummy,vecEventBinsCM,matPETM] = doPEP(cellUseSpikeTimes,vecEventBins,vecMergedPopEventTimes);
			vecEventBinsCM = vecEventBinsCM*1000;%ms
			matPopRateM = sum(matPETM,3);
			vecMeanM = mean(matPopRateM,1);
			vecSEMM = std(matPopRateM,[],1)./sqrt(size(matPopRateM,1));
			vecMeanM(vecEventBinsC==0)=nan;
			
			%merged random peak times
			[dummy,dummy,dummy,matPETMR] = doPEP(cellUseSpikeTimes,vecEventBins,vecRandMergedPopEventTimes);
			matPopRateMR = sum(matPETMR,3);
			vecMeanMR = mean(matPopRateMR,1);
			vecSEMMR = std(matPopRateMR,[],1)./sqrt(size(matPopRateMR,1));
			
			%short timescale
			dblStep = (1/100);
			dblRange = 2;%0.02
			vecEventBins = (-dblRange+dblStep/2):dblStep:(dblRange-dblStep/2);
			
			%merged peaks
			[dummy,dummy,vecEventBinsC_Short,matPETM_Short] = doPEP(cellUseSpikeTimes,vecEventBins_Short,vecMergedPopEventTimes);
			vecEventBinsC_Short = vecEventBinsC_Short*1000;%ms
			matPopRateM_Short = sum(matPETM_Short,3);
			vecMeanM_Short = mean(matPopRateM_Short,1);
			vecSEMM_Short = std(matPopRateM_Short,[],1)./sqrt(size(matPopRateM_Short,1));
			vecMeanM_Short(vecEventBinsC_Short==0)=nan;
			
			%merged random peak times
			[dummy,dummy,dummy,matPETMR_Short] = doPEP(cellUseSpikeTimes,vecEventBins_Short,vecRandMergedPopEventTimes);
			matPopRateMR_Short = sum(matPETMR_Short,3);
			vecMeanMR_Short = mean(matPopRateMR_Short,1);
			vecSEMMR_Short = std(matPopRateMR_Short,[],1)./sqrt(size(matPopRateMR_Short,1));
			
			%normalized pop rates
			vecMeanZ = (vecMeanM-mean(matPopRateMR(:)))./std(matPopRateMR(:));
			%}
			%ISIs
			dblStepISI = (1/30000);
			dblRangeISI = 0.04;%0.02
			vecEventBinsISI = 0:dblStepISI:dblRangeISI;
			vecEventBinsISIC = 1000*(vecEventBinsISI(2:end)-median(diff(vecEventBinsISI))/2);
			matISIC = nan(intNumN,numel(vecEventBinsISIC));
			for intNeuron=1:intNumN
				matISIC(intNeuron,:) = histcounts(diff(sort(cellUseSpikeTimes{intNeuron})),vecEventBinsISI);
			end
			
			% bin
			%IFR
			dblBinDur = 1/1000;
			vecBinEdgesIFR = min(vecStimTime):dblBinDur:(max(vecStimTime)+dblBinDur);
			vecBinTime = vecBinEdgesIFR(2:end)-median(diff(vecBinEdgesIFR))/2;
			[vecCounts,vecMeans] = makeBins(vecStimTime,vecStimIFR_Raw,vecBinEdgesIFR);
			%interpolate missing values
			indNan = isnan(vecMeans);
			vecBinnedIFR = interp1(vecBinTime(~indNan),vecMeans(~indNan),vecBinTime);
			
			%count spikes
			vecRate = vecCounts./dblBinDur;
			
			% DFA
			%https://en.wikipedia.org/wiki/Detrended_fluctuation_analysis
			dblBase = 1.2;
			intN = numel(vecBinnedIFR);
			vecIntervals = unique(round(dblBase.^(2:floor(log(intN)/log(dblBase)))));
			dblMin = 40/1000;
			intMin = max(3,dblMin/dblBinDur);
			dblMax = 10;
			intMax = dblMax/dblBinDur;
			vecIntervals(vecIntervals<intMin | vecIntervals>intMax)=[];
			[alphaRate, intervalsRate, fluctsRate] = fastdfa(vecRate', vecIntervals');
			
			%% plot
			figure;maxfig;
			subplot(2,3,1)
			plot(vecEventBinsC,vecMeanMR,'color','k')
			hold on
			plot(vecEventBinsC,vecMeanM,'color',lines(1))
			hold off
			%ylim([0 3000]);
			title(sprintf('Merged peaks; %s',strType),'interpreter','none');
			ylabel('Pop IFR (Hz)');
			xlabel('Time after merged pop event (s)');
			
			subplot(2,3,2)
			hold on
			plot(vecEventBinsC,vecMeanZ,'color',lines(1))
			hold off
			title(sprintf('Rate normalized to random'),'interpreter','none');
			ylabel('Normalized pop IFR (z-score)');
			xlabel('Time after merged pop event (s)');
			
			
			subplot(2,3,5)
			plot(intervalsRate,fluctsRate)
			set(gca,'xscale','log','yscale','log')
			title(sprintf('DFA Binned spiking rate; %s=%.3f',getGreek('alpha'),alphaRate))
			xlabel('Sample window (n)')
			ylabel('Fluctuation (F(n))');
			
			subplot(2,3,3)
			%average
			plot(vecEventBinsISIC,sum(matISIC,1))
			xlabel('ISI (ms)');
			ylabel('# of spikes');
			
			subplot(2,3,4)
			plot(vecEventBinsC_Short,vecMeanMR_Short,'color','k')
			hold on
			plot(vecEventBinsC_Short,vecMeanM_Short,'color',lines(1))
			hold off
			ylabel('Pop binned rate (Hz)');
			xlabel('Time after merged pop event (ms)');
			
			subplot(2,3,6)
			%per neuron
			imagesc(vecEventBinsISIC,1:intNumN,matISIC)
			xlabel('ISI (ms)');
			ylabel('Neuron #');
			
			fixfig;
			
			%% save plot
			export_fig(fullpath(strFigurePathSR,sprintf('T2_PeakTimescale_%s_%s.tif',strRec,strType)));
			export_fig(fullpath(strFigurePathSR,sprintf('T2_PeakTimescale_%s_%s.pdf',strRec,strType)));
			
			%% save data
			if intType == 1
				sReal.vecPeakHeight = vecPeakHeight;
				sReal.vecEventBinsC = vecEventBinsC;
				sReal.vecMeanM = vecMeanM;
			elseif intType == 2
				sPoiss.vecPeakHeight = vecPeakHeight;
				sPoiss.vecEventBinsC = vecEventBinsC;
				sPoiss.vecMeanM = vecMeanM;
			elseif intType == 3
				sShuffTid.vecPeakHeight = vecPeakHeight;
				sShuffTid.vecEventBinsC = vecEventBinsC;
				sShuffTid.vecMeanM = vecMeanM;
			elseif intType == 4
				sShuff.vecPeakHeight = vecPeakHeight;
				sShuff.vecEventBinsC = vecEventBinsC;
				sShuff.vecMeanM = vecMeanM;
			end
			
			fprintf('Finished %s [%s]\n',strType,getTime);
		end
		
		%% comparison
		figure;maxfig;
		subplot(2,3,1)
		vecBins = -0.5:0.02:1.5;
		vecBinsC = vecBins(2:end)-mean(diff(vecBins))/2;
		vecCounts_Real = histcounts(sReal.vecPeakHeight(:),vecBins);
		vecCounts_Shuff = histcounts(sShuff.vecPeakHeight(:),vecBins);
		vecCounts_ShuffTid = histcounts(sShuffTid.vecPeakHeight(:),vecBins);
		vecCounts_Poiss = histcounts(sPoiss.vecPeakHeight(:),vecBins);
		vecNormC_Real = vecCounts_Real./sum(vecCounts_Real);
		vecNormC_Shuff = vecCounts_Shuff./sum(vecCounts_Shuff);
		vecNormC_ShuffTid = vecCounts_ShuffTid./sum(vecCounts_ShuffTid);
		vecNormC_Poiss = vecCounts_Poiss./sum(vecCounts_Poiss);
		plot(vecBinsC,vecNormC_Real);
		hold on
		plot(vecBinsC,vecNormC_Shuff);
		plot(vecBinsC,vecNormC_Poiss);
		plot(vecBinsC,vecNormC_ShuffTid);
		hold off
		xlabel('Peak height')
		ylabel('Normalized count');
		legend({'Real','Shuffled','Poisson','ShuffTid'},'Location','best');
		
		subplot(2,3,2)
		plot(vecBinsC,(vecNormC_Real-vecNormC_Shuff+eps)./(vecNormC_Real+vecNormC_Shuff+eps));
		hold on
		plot(vecBinsC([1 end]),[0 0],'k--');
		xlabel('Peak height')
		ylabel('Real/shuffle ratio');
		
		%plot real minus shufftid
		subplot(2,3,3)
		plot(sReal.vecEventBinsC,sReal.vecMeanM-sShuffTid.vecMeanM,'color',lines(1))
		%ylim([0 3000]);
		title(sprintf('Merged peaks; %s',strType),'interpreter','none');
		ylabel('Pop IFR (Hz)');
		xlabel('Time after merged pop event (s)');
		
		fixfig;
		
		%% save plot
		export_fig(fullpath(strFigurePathSR,sprintf('T2_PeakTimescaleNorm_%s_%s.tif',strRec,strType)));
		export_fig(fullpath(strFigurePathSR,sprintf('T2_PeakTimescaleNorm_%s_%s.pdf',strRec,strType)));
			
		return
	end
end
toc
