%% aim
%{
https://www.nature.com/articles/s41598-021-95037-z
https://www.nature.com/articles/s42003-021-02437-y

q1: Can absence of stimulus information in initial period be explained by neurons responding only to
black/white patches in the first x ms?

q2: How precise are spike times in the natural movie repetitions?


%}
%% set parameters
cellDataTypes = {'Npx','Sim','ABI','SWN'};%topo, model, allen, nora
intRunDataType = 1;
strRunStim = 'DG';%DG or NM? => superseded to WS by SWN
cellTypes = {'Real','ShuffTid','Uniform'};%, 'UniformTrial', 'ShuffTid'};%, 'PoissGain'}; %PoissGain not done
boolFixSpikeGroupSize = false;
dblRemOnset = 0; %remove onset period in seconds; 0.125 for sim, 0.25 for npx
runHeaderPopTimeCoding;

%% go through recs
tic
for intRec=1:intRecNum %19 || weird: 11
	%% prep ABI or Npx data
	if strcmp(strRunType,'ABI')
		error to be updated
		runRecPrepABI;
		strThisRec = strRec;
	elseif strcmp(strRunType,'Sim')
		%load
		runRecPrepSim;
		
		%edit vars
		strThisRec = strRec;
		strDataPathT0=strDataPathSimT0;
		vecOri180 = mod(vecOrientation,180)*2;
		vecStimIdx = vecOri180;
		
		%get cell props
		%vecNeuronPrefOri = pi-[sLoad.sNeuron.PrefOri];
		vecNeuronType = [sLoad.sNeuron.Types]; %1=pyr,2=interneuron
		
	elseif strcmp(strRunType,'Npx') || strcmp(strRunType,'SWN')
		%prep
		runRecPrepNpx;
		strThisRec = strRec;
		strDataPathT0 = strTargetDataPath;
		
		
		%get cell props
		vecNeuronType = ones(size(indTuned)); %1=pyr,2=interneuron
		%narrow vs broad not done
	else
		error impossible
	end
	
	%% move onset
	%remove first x ms
	vecStimOnTime = vecStimOnTime + dblRemOnset;
	
	%% load prepped data
	strTarget = fullpath(strDataPathT0,[sprintf('T0Data_%s%s%s%s%s',strThisRec,strRunType,strRunStim,'Real') '.mat']);
	if ~exist(strTarget,'file')
		fprintf('Prepped T0 file did not exist for %s; skipping...\n',strThisRec);
		continue;
	end
	
	%check fr
	sSource = load(strTarget);
	cellSpikeTimes = sSource.cellSpikeTimes;
	intNumN = numel(cellSpikeTimes);
	cellSpikeTimesPerCellPerTrial = cell(intNumN,intOrigTrialNum);
	for intN=1:intNumN
		%real
		[vecTrialPerSpike,vecTimePerSpike] = getSpikesInTrial(cellSpikeTimes{intN},vecStimOnTime,dblStimDur);
		for intTrial=1:numel(vecStimOnTime)
			vecSpikeT = vecTimePerSpike(vecTrialPerSpike==intTrial);
			cellSpikeTimesPerCellPerTrial{intN,intTrial} = vecSpikeT;
		end
	end
	matData = cellfun(@numel,cellSpikeTimesPerCellPerTrial)./dblStimDur;
	if mean(sum(matData)) < 90%90 / 50
		fprintf('Avg # of spikes per trial was %.1f for %s; skipping...\n',mean(sum(matData)),strThisRec);
		continue;
	end
	
	
	%% get ori vars
	intTrialNum = numel(vecStimOnTime);
	vecOri180 = mod(vecOrientation,180);
	[vecOriIdx,vecUnique,vecPriorDistribution,cellSelect,vecRepetition] = val2idx(vecOri180);
	intOriNum = numel(vecUnique);
	intRepNum = min(vecPriorDistribution);
	dblStimDur = median(vecStimOffTime - vecStimOnTime);
	
	%types: Real, UniformTrial, ShuffTid, PoissGain
	for intType=1:numel(cellTypes)
		strType = cellTypes{intType};
		
		%% load prepped data and ifr
		sSource = load(fullpath(strDataPathT0,sprintf('T0Data_%s%s%s%s%s',strThisRec,strRunType,strRunStim,strType)));
		vecTime = sSource.vecTime;
		vecIFR = sSource.vecIFR;
		vecAllSpikeTime = sSource.vecAllSpikeTime;
		vecAllSpikeNeuron = sSource.vecAllSpikeNeuron;
		cellSpikeTimes = sSource.cellSpikeTimes;
		if numel(cellSpikeTimes) ~= intNumN
			error('neuron # mismatch!')
		end
		if isempty(vecTime),continue;end
		dblStartEpoch = sSource.dblStartEpoch;
		dblEpochDur = sSource.dblEpochDur;
		dblStopEpoch = dblStartEpoch + dblEpochDur;
		
		%% build trial-neuron cell matrix
		cellSpikeTimesPerCellPerTrial = cell(intNumN,intOrigTrialNum);
		for intN=1:intNumN
			%real
			[vecTrialPerSpike,vecTimePerSpike] = getSpikesInTrial(cellSpikeTimes{intN},vecStimOnTime,dblStimDur);
			for intTrial=1:intOrigTrialNum
				vecSpikeT = vecTimePerSpike(vecTrialPerSpike==intTrial);
				cellSpikeTimesPerCellPerTrial{intN,intTrial} = vecSpikeT;
			end
		end
		matData = cellfun(@numel,cellSpikeTimesPerCellPerTrial)./dblStimDur;
		sTuning = getTuningCurves(matData,vecOri180,0);
		vecNeuronPrefOri = sTuning.matFittedParams(:,1);
		vecNeuronBandwidth = real(sTuning.matBandwidth);
		
		%% check rates
		if boolFixSpikeGroupSize
			intSpikeGroupSize = 10;
			strSGS = ['Fixed' num2str(intSpikeGroupSize)];
		else
			intSpikeGroupSize = ceil(mean(sum(matData))/30); %20
			strSGS = ['Var' num2str(intSpikeGroupSize)];
		end
		
		%% detect peaks
		% filter
		%real
		intLag = round((numel(vecIFR)/numel(vecStimOnTime))/2);
		if (intLag/2) == round(intLag/2)
			intLag = intLag - 1;
		end
		dblThreshZ = 1;
		dblInfluence = 0.5;
		
		%% run analyses
		%add data
		sAllSpike = struct;
		sAllSpike.vecAllSpikeTime = vecAllSpikeTime;
		sAllSpike.vecAllSpikeNeuron = vecAllSpikeNeuron;
		sAllSpike.dblStartEpoch = dblStartEpoch;
		sAllSpike.dblEpochDur = dblEpochDur;
		
		%set params and run pop event detection
		dblCutOff = 0.65;
		sPeakOpts = struct;
		sPeakOpts.intLag = intLag;
		sPeakOpts.dblThreshZ = dblThreshZ;
		sPeakOpts.dblInfluence = dblInfluence;
		[sPopEvents,sMergedPopEvents,vecStimTime,vecStimIFR,vecStimIFR_Raw] = getPopEvents(vecIFR,vecTime,vecStimOnTime,sPeakOpts,sAllSpike,dblCutOff);
		
		%extract
		vecPopEventTimes = [sPopEvents.Time];
		vecPopEventLocs = [sPopEvents.Loc];
		vecMergedPopEventTimes = [sMergedPopEvents.Time];
		vecMergedPopEventLocs = [sMergedPopEvents.Loc];
		vecPeakHeight = vecStimIFR(vecMergedPopEventLocs);
		intHalfGroup = ceil(intSpikeGroupSize/2);
		vecRisingPhaseDeltaIFR = vecStimIFR(vecMergedPopEventLocs) - vecStimIFR(vecMergedPopEventLocs-intHalfGroup);
		vecFallingPhaseDeltaIFR = vecStimIFR(vecMergedPopEventLocs+intHalfGroup) - vecStimIFR(vecMergedPopEventLocs);
		%{
		%% get spike blocks in rising and falling phase of peaks
		intNextEvent = 2;
		dblNextEventT = vecMergedPopEventTimes(1);
		vecPhaseTrialIdx = zeros(1,(numel(vecMergedPopEventTimes)-1)*2);
		matPhaseData = nan(intNumN,(numel(vecMergedPopEventTimes)-1)*2);
		vecPhaseDeltaIFR = nan(1,(numel(vecMergedPopEventTimes)-1)*2);
		vecPrePost = zeros(1,(numel(vecMergedPopEventTimes)-1)*2);
		for intSpike=1:(numel(vecAllSpikeTime)-intSpikeGroupSize)
			if vecAllSpikeTime(intSpike) >= dblNextEventT
				intTrialIdx = sum(vecAllSpikeTime(intSpike)>vecStimOnTime);
				vecAssignTo = [-1 0]+(intNextEvent-1)*2;
				
				%increment event
				intNextEvent = intNextEvent + 1;
				if intNextEvent>numel(vecMergedPopEventTimes)
					break;
				else
					dblNextEventT = vecMergedPopEventTimes(intNextEvent);
				end
				
				%skip if in ITI
				if vecAllSpikeTime(intSpike) > vecStimOnTime(intTrialIdx)+dblStimDur
					continue;
				end
				
				vecPrecedingEntries = [(intSpike-intHalfGroup):(intSpike)]-1;
				vecPostEntries = [intSpike:(intSpike+intHalfGroup)]+1;
				vecPreN = accumarray(vecAllSpikeNeuron(vecPrecedingEntries)',ones(size(vecPrecedingEntries))',[intNumN 1]);
				vecPostN = accumarray(vecAllSpikeNeuron(vecPostEntries)',ones(size(vecPostEntries))',[intNumN 1]);
				matPhaseData(:,vecAssignTo(1)) = vecPreN;
				matPhaseData(:,vecAssignTo(2)) = vecPostN;
				vecPhaseTrialIdx(vecAssignTo) = intTrialIdx;
				vecPhaseDeltaIFR(vecAssignTo(1)) = vecIFR(intSpike) - vecIFR(intSpike-intHalfGroup);
				vecPhaseDeltaIFR(vecAssignTo(2)) = vecIFR(intSpike+intHalfGroup) - vecIFR(intSpike);
				vecPrePost(vecAssignTo) = [1 2];
			end
		end
		
		%remove 0
		indRem=vecPhaseTrialIdx==0;
		vecPhaseTrialIdx(indRem) = [];
		matPhaseData(:,indRem) = [];
		vecPhaseDeltaIFR(indRem) = [];
		vecPrePost(indRem) = [];
		vecPhaseTrialType = vecOriIdx(vecPhaseTrialIdx);
		
		%% do logistic regression
		dblLambda = 1;
		intTypeCV = 0;
		intOriNum = numel(unique(vecOriIdx));
		vecPriorDistribution = accumarray(vecPhaseTrialType(:),1,[intOriNum 1]);
		intVerbose = 1;
		[dblPerformanceCV,vecSpikeGroupDecodedTrial,matPosteriorProbability,dblMeanErrorDegs,matConfusion,matWeights,matAggActivation,matAggWeights,vecRepetition,vecTrialTypeIdx] = ...
			doCrossValidatedDecodingLR(matPhaseData,vecPhaseTrialType,intTypeCV,[],dblLambda,intVerbose);
		
		vecSpikeGroupCorrect = vecSpikeGroupDecodedTrial(:) == vecTrialTypeIdx;
		vecSpikeGroupConfidence = nan(size(vecSpikeGroupCorrect));
		for intStimType=1:size(matPosteriorProbability,1)
			vecSpikeGroupConfidence(vecTrialTypeIdx==intStimType) = matPosteriorProbability(intStimType,vecTrialTypeIdx==intStimType);
		end
		
		%% plot rising/falling phase analysis
		figure;maxfig;
		subplot(2,3,1)
		scatter(vecPhaseDeltaIFR,vecSpikeGroupConfidence,'.');
		[r,p]=corrcoef(vecPhaseDeltaIFR',vecSpikeGroupConfidence');
		
		subplot(2,3,2);
		vecConfRising = vecSpikeGroupConfidence(vecPrePost==1);
		vecConfFalling = vecSpikeGroupConfidence(vecPrePost==2);
		
		errorbar([1 2],[mean(vecConfRising) mean(vecConfFalling)],[std(vecConfRising) std(vecConfFalling)]./sqrt(numel(vecSpikeGroupConfidence)/2))
		%}
		%% plot
		if 1
			figure;maxfig;
			subplot(2,4,1)
			intStartTrial = 9;
			dblStartT = vecOrigStimOnTime(intStartTrial)+3;
			dblStopT = dblStartT + 1;
			intStopTrial = find(vecOrigStimOnTime > dblStopT,1)-1;
			intStartSample = find(vecStimTime > dblStartT,1);
			intStopSample = find(vecStimTime > dblStopT,1);
			vecPlotSamples = intStartSample:intStopSample;
			h0=plot(vecStimTime(vecPlotSamples),vecStimIFR_Raw(vecPlotSamples));
			hold on
			vecSubPopEvT = vecMergedPopEventLocs > intStartSample & vecMergedPopEventLocs < intStopSample;
			vecSubSampleLocs = vecMergedPopEventLocs(vecSubPopEvT);
			scatter(vecMergedPopEventTimes(vecSubPopEvT),vecStimIFR_Raw(vecSubSampleLocs),'o');
			scatter(vecOrigStimOnTime(intStartTrial:intStopTrial),median(vecStimIFR_Raw)*ones([1 intStopTrial-intStartTrial+1]),'k.');
			title(sprintf('%s; %s',strRec,strType),'interpreter','none');
			xlabel('Time (s)');
			ylabel('IFR (Hz)');
			xlim([dblStartT dblStopT]);
			
			subplot(2,4,2)
			h0=plot(vecStimTime,vecStimIFR_Raw);
			hold on
			scatter(vecMergedPopEventTimes,vecStimIFR_Raw(vecMergedPopEventLocs),'.');
			scatter(vecOrigStimOnTime,median(vecStimIFR_Raw)*ones(size(vecOrigStimOnTime)),'k.');
			title(sprintf('%s; %s',strRec,strType),'interpreter','none');
			xlabel('Time (s)');
			ylabel('IFR (Hz)');
			
			subplot(2,4,3)
			h1=plot(vecStimTime,vecStimIFR);
			hold on
			scatter(vecPopEventTimes,vecStimIFR(vecPopEventLocs),'.');
			scatter(vecOrigStimOnTime,zeros(size(vecOrigStimOnTime)),'k.');
			title(sprintf('%d peaks from filtered deviations',numel(vecPopEventLocs)));
			xlabel('Time (s)');
			ylabel('IFR normalized to local mean');
			
			subplot(2,4,4)
			h2=plot(vecStimTime,vecStimIFR);
			hold on
			scatter(vecMergedPopEventTimes,vecStimIFR(vecMergedPopEventLocs),'.');
			scatter(vecOrigStimOnTime,zeros(size(vecOrigStimOnTime)),'k.');
			title(sprintf('%d peaks after merging',numel(vecMergedPopEventLocs)));
			xlabel('Time (s)');
			ylabel('IFR normalized to local mean');
			
			vecTrialBins = 0:0.25:dblStimDur;
			vecTrialBinC = vecTrialBins(2:end)-mean(diff(vecTrialBins))/2;
			[vecMean,vecSEM,vecWindowBinCenters,matPET] = doPEP(cellUseSpikeTimes,vecTrialBins,vecOrigStimOnTime,-1);
			
			subplot(2,4,5)
			vecBinsE = 0:10:1500;
			vecBinsC = vecBinsE(2:end)-diff(vecBinsE(1:2))/2;
			vecC = histcounts(vecStimIFR_Raw,vecBinsE);
			dblMeanIFR = mean(vecStimIFR_Raw);
			dblMedianIFR = median(vecStimIFR_Raw);
			dblSdIFR = std(vecStimIFR_Raw);
			
			plot(vecBinsC,vecC/sum(vecC(:)));
			title(sprintf('IFR distro %s',strType));
			xlabel('IFR')
			ylabel('Norm count');
			
			subplot(2,4,6)
			plotRaster(vecStimTime(vecMergedPopEventLocs),vecOrigStimOnTime,dblStimDur,inf);
			title('Population events');
			
			subplot(2,4,7)
			intPlotN = round(size(matPET,3)/2);
			imagesc(vecTrialBinC,[],matPET(:,:,intPlotN));
			title(sprintf('Example single neuron firing rate (#%d)',intPlotN));
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
			set(gca,'clim',[0 1500]);
			
			fixfig;
			h0.LineWidth = 1;
			h1.LineWidth = 1;
			h2.LineWidth = 1;
			
			export_fig(fullpath(strFigurePathSR,sprintf('T0_PeakDetection_%s_%s.tif',strRec,strType)));
			export_fig(fullpath(strFigurePathSR,sprintf('T0_PeakDetection_%s_%s.pdf',strRec,strType)));
		end
		
		%% what does peak look like? plot PSTH
		%raw peaks
		dblStep = (1/30000);
		vecEventBins = (-0.02+dblStep/2):dblStep:(0.02-dblStep/2);
		[dummy,dummy,vecEventBinsC,matPET] = doPEP(cellUseSpikeTimes,vecEventBins,vecPopEventTimes);
		vecEventBinsC = vecEventBinsC*1000;%ms
		matPopRate = sum(matPET,3);
		vecMean = mean(matPopRate,1);
		vecSEM = std(matPopRate,[],1)./sqrt(size(matPopRate,1));
		vecMean(vecEventBinsC==0)=nan;
		
		%raw random peak times
		vecRandPopEventTimes = linspace(min(vecPopEventTimes)-1,max(vecPopEventTimes)+1,numel(vecPopEventTimes));
		[dummy,dummy,dummy,matPETR] = doPEP(cellUseSpikeTimes,vecEventBins,vecRandPopEventTimes);
		matPopRateR = sum(matPETR,3);
		vecMeanR = mean(matPopRateR,1);
		vecSEMR = std(matPopRateR,[],1)./sqrt(size(matPopRateR,1));
		
		%merged peaks
		[dummy,dummy,vecEventBinsC,matPETM] = doPEP(cellUseSpikeTimes,vecEventBins,vecMergedPopEventTimes);
		matPopRateM = sum(matPETM,3);
		vecMeanM = mean(matPopRateM,1);
		vecSEMM = std(matPopRateM,[],1)./sqrt(size(matPopRateM,1));
		vecMeanM(vecEventBinsC==0)=nan;
		
		%merged random peak times
		vecRandMergedPopEventTimes = linspace(min(vecMergedPopEventTimes)-1,max(vecMergedPopEventTimes)+1,numel(vecMergedPopEventTimes));
		[dummy,dummy,dummy,matPETMR] = doPEP(cellUseSpikeTimes,vecEventBins,vecRandMergedPopEventTimes);
		matPopRateMR = sum(matPETMR,3);
		vecMeanMR = mean(matPopRateMR,1);
		vecSEMMR = std(matPopRateMR,[],1)./sqrt(size(matPopRateMR,1));
		
		%% output
		%raw peaks; events range from -1ms to +1ms around peak
		dblStep = (1/30000);
		dblStartEvent = -0.002+dblStep/2;
		dblStopEvent = 0.002-dblStep/2;
		vecPeakBins = dblStartEvent:dblStep:dblStopEvent;
		%[cellEventPerSpike,cellTimePerSpike] = getSpikesInTrial(cellUseSpikeTimes,vecPopEventTimes+dblStartEvent,range(vecPeakBins));
		
		%merged peaks
		dblStep = (1/30000);
		dblStartEvent = -0.002+dblStep/2;
		dblStopEvent = 0.002-dblStep/2;
		vecPeakBins = dblStartEvent:dblStep:dblStopEvent;
		%[cellMergedEventPerSpike,cellMergedTimePerSpike] = getSpikesInTrial(cellUseSpikeTimes,vecMergedPopEventTimes+dblStartEvent,range(vecPeakBins));
		
		%neuron rate
		matNeuronRate = squeeze(sum(matPET,1));
		matNeuronRate(vecEventBinsC==0,:)=0;
		%smooth
		matNeuronRate = imfilt(matNeuronRate,normpdf(-10:10,0,5)'./sum(normpdf(-10:10,0,5)));
		matNeuronRate = (matNeuronRate - mean(matNeuronRate,1)) ./ std(matNeuronRate,[],1);
		intNeuronNum = size(matNeuronRate,2);
		
		% 			%order cells by average latency
		% 			vecOrderBinsC = vecOrderBins(2:end)-mean(diff(vecOrderBins))/2;
		% 			matSpikes = nan(intNumN,numel(vecOrderBinsC));
		% 			vecLatencyMu = nan(intNumN,1);
		% 			vecLatencySd = nan(intNumN,1);
		% 			for intN=1:intNumN
		% 				vecT = cellTimePerSpike{intN}+dblStartEvent;
		% 				indRemCenter = vecT>-dblStep & vecT<dblStep;
		% 				matSpikes(intN,:) = histcounts(vecT,vecOrderBins);
		% 				vecLatencyMu(intN) = mean(vecT(~indRemCenter));
		% 				vecLatencySd(intN) = std(vecT(~indRemCenter));
		% 			end
		% 			%
		% 			[vecLatencyMu,vecReorder]=sort(vecLatencyMu);
		% 			vecLatencySd = vecLatencySd(vecReorder);
		% 			matSpikes = matSpikes(vecReorder,:);
		% 			cellTimePerSpike = cellTimePerSpike(vecReorder);
		% 			matSpikesNorm = matSpikes./mean(matSpikes,2);
		% 			matSpikesNorm(:,(floor(numel(vecOrderBins)/2))) = 0;
		% 			matSpikesDev = matSpikesNorm - mean(matSpikesNorm,1);
		% 			matSpikesSmooth = imfilt(matSpikesNorm,normpdf(-2:2)./sum(normpdf(-2:2)));
		
		
		%% plot
		figure;maxfig;
		subplot(2,3,1)
		plot(vecEventBinsC,vecMeanR,'color','k')
		hold on
		plot(vecEventBinsC,vecMean,'color',lines(1))
		hold off
		ylim([0 3000]);
		title(sprintf('%s; %s',strRec,strType),'interpreter','none');
		ylabel('Pop rate (Hz)');
		xlabel('Time after pop event (ms)');
		drawnow;
		
		subplot(2,3,4)
		plot(vecEventBinsC,vecMeanMR,'color','k')
		hold on
		plot(vecEventBinsC,vecMeanM,'color',lines(1))
		hold off
		ylim([0 3000]);
		title('Merged peaks','interpreter','none');
		ylabel('Pop rate (Hz)');
		xlabel('Time after merged pop event (ms)');
		
		
		%subplot(2,3,3)
		%imagesc(vecOrderBinsC*1000,1:intNumN,matSpikesNorm)
		%xlabel('Time after pop event (ms)');
		%ylabel('Neuron #');
		%title(strType)
		
		% 			subplot(2,3,4)
		% 			vecCoM = nan(1,intNeuronNum);
		% 			for intNeuron=1:intNeuronNum
		% 				vecCoM(intNeuron) = calcCenterOfMass(matNeuronRate(:,intNeuron),1);
		% 			end
		% 			imagesc(vecEventBinsC,1:intNeuronNum,matNeuronRate');
		% 			title('Smoothed and normalized rates around pop events')
		% 			xlabel('Time after pop event (ms)');
		% 			ylabel('Neuron #');
		
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
		
		%% save plot
		export_fig(fullpath(strFigurePathSR,sprintf('T1_PeakPSTH_%s_%s.tif',strRec,strType)));
		export_fig(fullpath(strFigurePathSR,sprintf('T1_PeakPSTH_%s_%s.pdf',strRec,strType)));
		
		%% save data
		sData = struct;
		sData.matMean = matPopRate;
		sData.vecEventBins = vecEventBins;
		sData.vecOrderBins = vecPeakBins;
		sData.vecPeakHeight = vecPeakHeight;
		sData.vecBinsIFR = vecBinsC;
		sData.vecBinsCount = vecC;
		sData.dblMeanIFR = dblMeanIFR;
		sData.dblMedianIFR = dblMedianIFR;
		sData.dblSdIFR = dblSdIFR;
		if strcmp(strType,'Real')
			sReal = sData;
		elseif strcmp(strType,'Poiss')
			sPoiss = sData;
		elseif strcmp(strType,'ShuffTid')
			sShuffTid = sData;
		elseif strcmp(strType,'Shuff')
			sShuff = sData;
		elseif strcmp(strType,'PoissGain')
			sPoissGain = sData;
		end
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
	hold off
	xlabel('Peak height')
	ylabel('Normalized count');
	legend({'Real','Shuffled','Poisson'},'Location','best');
	
	subplot(2,3,2)
	plot(vecBinsC,(vecNormC_Real-vecNormC_Shuff+eps)./(vecNormC_Real+vecNormC_Shuff+eps));
	hold on
	plot(vecBinsC([1 end]),[0 0],'k--');
	xlabel('Peak height')
	ylabel('Real/shuffle ratio');
	close all;
	
	%% save data
	save(fullpath(strTargetDataPath,sprintf('T1Data_%s',strRec)),...
		'sReal',...
		'sPoiss',...
		'sShuffTid',...
		'sShuff',...
		'sPoissGain',...
		'strRec'...
		);
end

