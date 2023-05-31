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
boolHome = false;
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
		[matMeanRate,cellSpikeTimes] = ...
			NpxPrepMovieData(cellSpikeTimesRaw,vecStimOnTime,vecStimOffTime,vecFrameIdx);
		
		%% load prepro T0 data
		% pool spikes from all neurons, but save the time+id per spike, then calculate IFR over all
		% spikes at pop level, detect peaks, and split into trials.
		%get spikes per trial per neuron
		sLoad = load(fullpath(strTargetDataPath,sprintf('T0Data_%s',strRec)));
		
		vecAllSpikeNeuron_Real = sLoad.vecAllSpikeNeuron_Real;
		vecAllSpikeTime_Real = sLoad.vecAllSpikeTime_Real;
		vecTime_Real = sLoad.vecTime_Real;
		vecIFR_Real = sLoad.vecIFR_Real;
		
		vecAllSpikeNeuron_Shuff = sLoad.vecAllSpikeNeuron_Shuff;
		vecAllSpikeTime_Shuff = sLoad.vecAllSpikeTime_Shuff;
		vecTime_Shuff = sLoad.vecTime_Shuff;
		vecIFR_Shuff = sLoad.vecIFR_Shuff;
		
		vecAllSpikeNeuron_Poiss = sLoad.vecAllSpikeNeuron_Poiss;
		vecAllSpikeTime_Poiss = sLoad.vecAllSpikeTime_Poiss;
		vecTime_Poiss = sLoad.vecTime_Poiss;
		vecIFR_Poiss = sLoad.vecIFR_Poiss;
		
		dblStartEpoch = sLoad.dblStartEpoch;
		dblEpochDur = sLoad.dblEpochDur;
		dblStopEpoch = dblStartEpoch+dblEpochDur;
		
		%% re-generate poisson
		intNumN = numel(cellSpikeTimes);
		%remove spikes outside epoch
		for intN=1:intNumN
			vecT = cellSpikeTimes{intN};
			indRem = (vecT > dblStopEpoch) | (vecT < dblStartEpoch);
			cellSpikeTimes{intN} = vecT(~indRem);
		end
		
		%generate poisson
		cellSpikeTimes_Poiss = cell(1,intNumN);
		for intN=1:intNumN
			dblT0 = cellSpikeTimes{intN}(1);
			dblTN = cellSpikeTimes{intN}(end);
			dblTotT = dblTN-dblT0;
			%add 3ms refractory period
			dblRt = (3/1000);
			dblLambda = numel(cellSpikeTimes{intN})/(dblTotT-dblRt*numel(cellSpikeTimes{intN}));
			vecISI = dblRt+exprnd(1./dblLambda,[1,round(dblLambda*dblTotT*2)]);
			vecSpikeT = dblT0+cumsum(vecISI)-vecISI(1);
			cellSpikeTimes_Poiss{intN} = vecSpikeT(vecSpikeT<dblTN);
		end
		
		%get spikes per trial per neuron
		intTotS = sum(cellfun(@numel,cellSpikeTimes));
		intTotS_Poiss = sum(cellfun(@numel,cellSpikeTimes_Poiss));
		vecAllSpikeTime_Poiss = nan(1,intTotS_Poiss);
		vecAllSpikeNeuron_Poiss = zeros(1,intTotS_Poiss,'int16');
		intSP = 1;
		for intN=1:intNumN
			%poisson
			intThisSP = numel(cellSpikeTimes_Poiss{intN});
			vecSpikeT_Poiss = cellSpikeTimes_Poiss{intN};
			vecAllSpikeTime_Poiss(intSP:(intSP+intThisSP-1)) = vecSpikeT_Poiss;
			vecAllSpikeNeuron_Poiss(intSP:(intSP+intThisSP-1)) = intN;
			intSP = intSP + intThisSP;
		end
		indRemPoiss = (vecAllSpikeTime_Poiss > dblStopEpoch) | (vecAllSpikeTime_Poiss < dblStartEpoch);
		vecAllSpikeTime_Poiss(indRemPoiss) = [];
		vecAllSpikeNeuron_Poiss(indRemPoiss) = [];
		
		%poisson
		[vecAllSpikeTime_Poiss,vecReorder] = sort(vecAllSpikeTime_Poiss);
		vecAllSpikeNeuron_Poiss = vecAllSpikeNeuron_Poiss(vecReorder);
		vecTime_Poiss = vecAllSpikeTime_Poiss;
		
		%% replace IFR with binned rates
		dblBinDur = (5/1000);
		vecBins=(dblStartEpoch:dblBinDur:dblStopEpoch)';
		vecIFR_Real = flat(histcounts(sLoad.vecTime_Real,vecBins)./dblBinDur);
		vecTime_Real = vecBins(2:end)-dblBinDur/2;
		
		vecIFR_Shuff = flat(histcounts(sLoad.vecTime_Shuff,vecBins)./dblBinDur);
		vecTime_Shuff = vecBins(2:end)-dblBinDur/2;
		
		vecIFR_Poiss = flat(histcounts(vecTime_Poiss,vecBins)./dblBinDur);
		vecTime_Poiss = vecBins(2:end)-dblBinDur/2;
		
		%% detect peaks
		% filter
		%real
		intLag = round((numel(vecIFR_Real)/numel(vecOrigStimOnTime))/2);
		if (intLag/2) == round(intLag/2)
			intLag = intLag - 1;
		end
		dblThreshZ = 1;
		dblInfluence = 0.5;
		[signals,avgFilter_Real,stdFilter_Real] = detectpeaks(vecIFR_Real,intLag,dblThreshZ,dblInfluence);
		%vecIFR_Real = (vecIFR_Real - avgFilter_Real)./stdFilter_Real;
		vecNormIFR_Real = (vecIFR_Real - avgFilter_Real)./avgFilter_Real;
		indRem_Real = isinf(vecNormIFR_Real);
		vecNormIFR_Real(indRem_Real) = [];
		vecNormTime_Real  = vecTime_Real(~indRem_Real);
		
		%shuff
		[signals,avgFilter_Shuff,stdFilter_Shuff] = detectpeaks(vecIFR_Shuff,intLag,dblThreshZ,dblInfluence);
		%vecIFR_Shuff = (vecIFR_Shuff - avgFilter_Shuff)./stdFilter_Shuff;
		vecNormIFR_Shuff = (vecIFR_Shuff - avgFilter_Shuff)./avgFilter_Shuff;
		indRem_Shuff = isinf(vecNormIFR_Shuff);
		vecNormIFR_Shuff(indRem_Shuff) = [];
		vecNormTime_Shuff = vecTime_Shuff(~indRem_Shuff);
		
		%poiss
		[signals,avgFilter_Poiss,stdFilter_Poiss] = detectpeaks(vecIFR_Poiss,intLag,dblThreshZ,dblInfluence);
		%vecIFR_Poiss = (vecIFR_Poiss - avgFilter_Poiss)./stdFilter_Poiss;
		vecNormIFR_Poiss = (vecIFR_Poiss - avgFilter_Poiss)./avgFilter_Poiss;
		indRem_Poiss = isinf(vecNormIFR_Poiss);
		vecNormIFR_Poiss(indRem_Poiss) = [];
		vecNormTime_Poiss = vecTime_Poiss(~indRem_Poiss);
		
		% real peaks
		indStimSpikes_Real = vecNormTime_Real>(vecOrigStimOnTime(1)-dblStimDur) & vecNormTime_Real<(vecOrigStimOnTime(end)+dblStimDur);
		vecStimIFR_RealRaw = vecIFR_Real(indStimSpikes_Real);
		vecStimIFR_Real = vecNormIFR_Real(indStimSpikes_Real);
		vecStimTime_Real = vecNormTime_Real(indStimSpikes_Real);
		
		[vecPeakHeight_Real,vecPeakLocs_Real,w_Real,p_Real] = findpeaks(vecStimIFR_Real);
		
		% shuff peaks
		indStimSpikes_Shuff = vecNormTime_Shuff>(vecOrigStimOnTime(1)-dblStimDur) & vecNormTime_Shuff<(vecOrigStimOnTime(end)+dblStimDur);
		vecStimIFR_Shuff = vecNormIFR_Shuff(indStimSpikes_Shuff);
		vecStimTime_Shuff = vecNormTime_Shuff(indStimSpikes_Shuff);
		
		[vecPeakHeight_Shuff,vecPeakLocs_Shuff,w_Shuff,p_Shuff] = findpeaks(vecStimIFR_Shuff);
		
		% poiss peaks
		indStimSpikes_Poiss = vecNormTime_Poiss>(vecOrigStimOnTime(1)-dblStimDur) & vecNormTime_Poiss<(vecOrigStimOnTime(end)+dblStimDur);
		vecStimIFR_Poiss = vecNormIFR_Poiss(indStimSpikes_Poiss);
		vecStimTime_Poiss = vecNormTime_Poiss(indStimSpikes_Poiss);
		
		[vecPeakHeight_Poiss,vecPeakLocs_Poiss,w_Poiss,p_Poiss] = findpeaks(vecStimIFR_Poiss);
		
		%merge peaks
		%dblCutOff = 0.6;
		dblCutOff = 2.5;
		vecStimPeakLocs_Real = vecPeakLocs_Real(vecPeakHeight_Real>dblCutOff);
		vecStimPeakLocs_Shuff = vecPeakLocs_Shuff(vecPeakHeight_Shuff>dblCutOff);
		vecStimPeakLocs_Poiss = vecPeakLocs_Shuff(vecPeakHeight_Poiss>dblCutOff);
		vecT = vecStimTime_Real;
		vecV = vecStimIFR_Real;
		vecP = vecStimPeakLocs_Real;
		%matPeakDomain = mergepeaks(vecStimTime_Real,vecStimIFR_Real,vecStimPeakLocs_Real);
		
		%% transform time indices
		%events
		dblStartEpoch = vecOrigStimOnTime(1)-dblStimDur;
		dblEpochDur = vecOrigStimOnTime(end)-vecOrigStimOnTime(1)+dblStimDur;
		
		%retain only stim epoch
		indStimEpoch = vecAllSpikeTime_Real > dblStartEpoch & vecAllSpikeTime_Real < (dblStartEpoch + dblEpochDur);
		vecPopSpikeTime_Real = vecAllSpikeTime_Real(indStimEpoch);
		vecPopNeuronId_Real = vecAllSpikeNeuron_Real(indStimEpoch);
		
		%get event times
		vecPopEventTimes_Real = vecStimTime_Real(vecStimPeakLocs_Real);
		vecPopEventTimes_Shuff = vecStimTime_Shuff(vecStimPeakLocs_Shuff);
		vecPopEventTimes_Poiss = vecStimTime_Poiss(vecStimPeakLocs_Poiss);
		intPopEventNum = numel(vecPopEventTimes_Real);
		vecPopEventLocs_Real = nan(1,intPopEventNum);
		vecPopEventLocsIFR_Real = nan(1,intPopEventNum);
		
		for intEvent=1:intPopEventNum
			[dummy,vecPopEventLocs_Real(intEvent)] = min(abs(vecPopEventTimes_Real(intEvent)-vecPopSpikeTime_Real));
			[dummy,vecPopEventLocsIFR_Real(intEvent)] = min(abs(vecPopEventTimes_Real(intEvent)-vecStimTime_Real));
		end
		
		%% plot
		figure;maxfig;
		subplot(2,4,1)
		h0=plot(vecStimTime_Real,vecStimIFR_RealRaw);
		hold on
		scatter(vecOrigStimOnTime,median(vecStimIFR_RealRaw)*ones(size(vecOrigStimOnTime)),'k.');
		title(sprintf('%s; raw IFR',strRec),'interpreter','none');
		xlabel('Time (s)');
		ylabel('IFR (Hz)');
		
		subplot(2,4,2)
		h1=plot(vecStimTime_Real,vecStimIFR_Real);
		hold on
		scatter(vecPopEventTimes_Real,vecStimIFR_Real(vecPopEventLocsIFR_Real),'.');
		scatter(vecOrigStimOnTime,zeros(size(vecOrigStimOnTime)),'k.');
		title(sprintf('%d peaks from filtered deviations',numel(vecPopEventLocsIFR_Real)));
		xlabel('Time (s)');
		ylabel('IFR normalized to local mean');
		
		subplot(2,4,3)
		h2=plot(vecStimTime_Real,vecStimIFR_Real);
		hold on
		%scatter(vecStimTime_Real(matPeakDomain(:,1)),vecStimIFR_Real(matPeakDomain(:,1)),'.');
		%scatter(vecOrigStimOnTime,zeros(size(vecOrigStimOnTime)),'k.');
		%title(sprintf('%d peaks after merging',numel(matPeakDomain(:,1))));
		xlabel('Time (s)');
		ylabel('IFR normalized to local mean');
		
		vecTrialBins = 0:0.25:dblStimDur;
		vecTrialBinC = vecTrialBins(2:end)-mean(diff(vecTrialBins))/2;
		[vecMean,vecSEM,vecWindowBinCenters,matPET] = doPEP(cellSpikeTimes,vecTrialBins,vecOrigStimOnTime,-1);
		
		subplot(2,4,4)
		%plotRaster(vecStimTime_Real(matPeakDomain(:,1)),vecOrigStimOnTime,dblStimDur,inf);
		title('Population events');
		
		subplot(2,4,5)
		vecBins = -0.5:0.02:1.5;
		vecBinsC = vecBins(2:end)-mean(diff(vecBins))/2;
		vecCounts_Real = histcounts(vecPeakHeight_Real(:),vecBins);
		vecCounts_Shuff = histcounts(vecPeakHeight_Shuff(:),vecBins);
		vecCounts_Poiss = histcounts(vecPeakHeight_Poiss(:),vecBins);
		vecNormC_Real = vecCounts_Real./sum(vecCounts_Real);
		vecNormC_Shuff = vecCounts_Shuff./sum(vecCounts_Shuff);
		vecNormC_Poiss = vecCounts_Poiss./sum(vecCounts_Poiss);
		plot(vecBinsC,vecNormC_Real);
		hold on
		plot(vecBinsC,vecNormC_Shuff);
		plot(vecBinsC,vecNormC_Poiss);
		hold off
		xlabel('Peak height')
		ylabel('Normalized count');
		legend({'Real','Shuffled','Poisson'},'Location','best');
		
		subplot(2,4,6)
		plot(vecBinsC,(vecNormC_Real-vecNormC_Shuff+eps)./(vecNormC_Real+vecNormC_Shuff+eps));
		hold on
		plot(vecBinsC([1 end]),[0 0],'k--');
		xlabel('Peak height')
		ylabel('Real/shuffle ratio');
		
		subplot(2,4,7)
		imagesc(vecTrialBinC,[],matPET(:,:,40));
		title('Example single neuron firing rate (#40)');
		xlabel('Time (s)')
		ylabel('Trial #');
		colorbar;
		axis xy
		
		subplot(2,4,8)
		imagesc(vecTrialBinC,[],sum(matPET,3));
		title('Population firing rate');
		xlabel('Time (s)')
		ylabel('Trial #');
		colorbar;
		axis xy
		
		fixfig;
		h0.LineWidth = 1;
		h1.LineWidth = 1;
		h2.LineWidth = 1;
		
		export_fig(fullpath(strFigurePath,sprintf('T0_PeakDetection_%s.tif',strRec)));
		export_fig(fullpath(strFigurePath,sprintf('T0_PeakDetection_%s.pdf',strRec)));
		
		%% what does peak look like? plot PSTH
		%reconstruct spiking arrays
		intNumN = numel(cellSpikeTimes);
		cellSpikeTimes_Real = cell(1,intNumN);
		cellSpikeTimes_Shuff = cell(1,intNumN);
		cellSpikeTimes_Poiss = cell(1,intNumN);
		for intN=1:intNumN
			cellSpikeTimes_Real{intN} = vecAllSpikeTime_Real(vecAllSpikeNeuron_Real==intN);
			cellSpikeTimes_Shuff{intN} = vecAllSpikeTime_Shuff(vecAllSpikeNeuron_Shuff==intN);
			cellSpikeTimes_Poiss{intN} = vecAllSpikeTime_Poiss(vecAllSpikeNeuron_Poiss==intN);
		end
		
		sReal = struct;
		sShuff = struct;
		sPoiss = struct;
		for intType=1:3
			if intType==1
				cellUseSpikeTimes = cellSpikeTimes_Real;
				vecUsePopEventTimes = vecPopEventTimes_Real;
				strType = 'Real spikes';
				%elseif intType == 2
				%	cellUseSpikeTimes = cellSpikeTimes;
				%	vecUsePopEventTimes = vecPopEventTimes_Real;
				%	strType = 'Original spikes';
			elseif intType == 2
				cellUseSpikeTimes = cellSpikeTimes_Poiss;
				vecUsePopEventTimes = vecPopEventTimes_Poiss;
				strType = 'Poisson spikes';
			else
				cellUseSpikeTimes = cellSpikeTimes_Shuff;
				vecUsePopEventTimes = vecPopEventTimes_Shuff;
				strType = 'Shuffled spikes';
			end
			dblStep = (1/30000);
			vecEventBins = (-0.02+dblStep/2):dblStep:(0.02-dblStep/2);
			[dummy,dummy,vecEventBinsC,matPET] = doPEP(cellUseSpikeTimes,vecEventBins,vecUsePopEventTimes);
			vecEventBinsC = vecEventBinsC*1000;%ms
			matPopRate = sum(matPET,3);
			vecMean = mean(matPopRate,1);
			vecSEM = std(matPopRate,[],1)./sqrt(size(matPopRate,1));
			
			%events range from -1ms to +1ms around peak
			dblStep = (1/30000);
			dblStartEvent = -0.002+dblStep/2;
			dblStopEvent = 0.002-dblStep/2;
			vecOrderBins = dblStartEvent:dblStep:dblStopEvent;
			[cellEventPerSpike,cellTimePerSpike] = getSpikesInTrial(cellUseSpikeTimes,vecUsePopEventTimes+dblStartEvent,range(vecOrderBins));
			
			%order cells by average latency
			vecOrderBinsC = vecOrderBins(2:end)-mean(diff(vecOrderBins))/2;
			matSpikes = nan(intNumN,numel(vecOrderBinsC));
			vecLatencyMu = nan(intNumN,1);
			vecLatencySd = nan(intNumN,1);
			for intN=1:intNumN
				vecT = cellTimePerSpike{intN}+dblStartEvent;
				indRemCenter = vecT>-dblStep & vecT<dblStep;
				matSpikes(intN,:) = histcounts(vecT,vecOrderBins);
				vecLatencyMu(intN) = mean(vecT(~indRemCenter));
				vecLatencySd(intN) = std(vecT(~indRemCenter));
			end
			%
			[vecLatencyMu,vecReorder]=sort(vecLatencyMu);
			vecLatencySd = vecLatencySd(vecReorder);
			matSpikes = matSpikes(vecReorder,:);
			cellTimePerSpike = cellTimePerSpike(vecReorder);
			matSpikesNorm = matSpikes./mean(matSpikes,2);
			matSpikesNorm(:,(floor(numel(vecOrderBins)/2))) = 0;
			matSpikesDev = matSpikesNorm - mean(matSpikesNorm,1);
			matSpikesSmooth = imfilt(matSpikesNorm,normpdf(-2:2)./sum(normpdf(-2:2)));
			
			vecMean(vecEventBinsC==0)=nan;
			
			%% plot
			figure;maxfig;
			subplot(2,3,1)
			plot(vecEventBinsC,vecMean)
			ylim([0 2500]);
			title('Pop rate; spikes at t=0 removed');
			ylabel('Pop rate (Hz)');
			xlabel('Time after pop event (ms)');
			
			%subplot(2,3,3)
			%imagesc(vecOrderBinsC*1000,1:intNumN,matSpikesNorm)
			%xlabel('Time after pop event (ms)');
			%ylabel('Neuron #');
			%title(strType)
			
			subplot(2,3,4)
			%neuron rate
			matNeuronRate = squeeze(sum(matPET,1));
			matNeuronRate(vecEventBinsC==0,:)=0;
			%smooth
			matNeuronRate = imfilt(matNeuronRate,normpdf(-10:10,0,5)'./sum(normpdf(-10:10,0,5)));
			matNeuronRate = (matNeuronRate - mean(matNeuronRate,1)) ./ std(matNeuronRate,[],1);
			intNeuronNum = size(matNeuronRate,2);
			vecCoM = nan(1,intNeuronNum);
			for intNeuron=1:intNeuronNum
				vecCoM(intNeuron) = calcCenterOfMass(matNeuronRate(:,intNeuron),1);
			end
			imagesc(vecEventBinsC,1:intNeuronNum,matNeuronRate');
			title('Smoothed and normalized rates around pop events')
			xlabel('Time after pop event (ms)');
			ylabel('Neuron #');
			
			%sum of auto-correlograms
			matACG = nan(size(matNeuronRate));
			parfor intNeuron=1:intNeuronNum
				vecACG = zeros(1,numel(vecEventBinsC));
				vecT = cellUseSpikeTimes{intNeuron};
				for intSpikeIdx=1:numel(vecT)
					vecACG = vecACG + histcounts(vecT-vecT(intSpikeIdx),vecEventBins);
				end
				matACG(:,intNeuron) = vecACG;
			end
			matACG(vecEventBinsC==0,:)=0;
			
			subplot(2,3,2)
			plot(vecEventBinsC,sum(matACG,2))
			ylim([0 2500]);
			title('Sum of auto-correlograms');
			ylabel('Avg rate (Hz)');
			xlabel('Time after spike (ms)');
			
			subplot(2,3,5)
			%smooth
			matPlotACG = imfilt(matACG,normpdf(-10:10,0,5)'./sum(normpdf(-10:10,0,5)));
			matPlotACG = (matPlotACG) ./ std(matPlotACG,[],1);
			imagesc(vecEventBinsC,1:intNeuronNum,matPlotACG');
			title('Smoothed and normalized ACGs')
			xlabel('Time after spike (ms)');
			ylabel('Neuron #');
			
			
			%plot ISIs
			vecEventBinsISI = vecEventBins(vecEventBins>=0);
			vecEventBinsISIC = 1000*(vecEventBinsISI(2:end)-median(diff(vecEventBinsISI))/2);
			matISIC = nan(intNeuronNum,numel(vecEventBinsISIC));
			for intNeuron=1:intNeuronNum
				matISIC(intNeuron,:) = histcounts(diff(sort(cellUseSpikeTimes{intNeuron})),vecEventBinsISI);
			end
			
			subplot(2,3,3)
			%average
			plot(vecEventBinsISIC,sum(matISIC,1))
			xlabel('ISI (ms)');
			ylabel('# of spikes');
			
			subplot(2,3,6)
			%per neuron
			imagesc(vecEventBinsISIC,1:intNeuronNum,matISIC)
			xlabel('ISI (ms)');
			ylabel('Neuron #');
			
			fixfig;
			
			if intType == 1
				sReal.cellEventPerSpike = cellEventPerSpike;
				sReal.cellTimePerSpike = cellTimePerSpike;
				sReal.matMean = matPopRate;
				sReal.vecEventBins = vecEventBins;
				sReal.vecOrderBins = vecOrderBins;
			elseif intType == 2
				sPoiss.cellEventPerSpike = cellEventPerSpike;
				sPoiss.cellTimePerSpike = cellTimePerSpike;
				sPoiss.matMean = matPopRate;
				sPoiss.vecEventBins = vecEventBins;
				sPoiss.vecOrderBins = vecOrderBins;
			elseif intType == 3
				sShuff.cellEventPerSpike = cellEventPerSpike;
				sShuff.cellTimePerSpike = cellTimePerSpike;
				sShuff.matMean = matPopRate;
				sShuff.vecEventBins = vecEventBins;
				sShuff.vecOrderBins = vecOrderBins;
			end
		end
		return
	end
end
toc
