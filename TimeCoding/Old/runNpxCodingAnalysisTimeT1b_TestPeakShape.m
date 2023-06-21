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
		%get data matrix & remove cells with rates <0.1Hz
		fprintf('Running %s (%d/%d) [%s]\n',strRec,intRec,numel(sAggStim),getTime);
		[matMeanRate,cellRawSpikeTimes] = ...
			NpxPrepMovieData({sArea1Neurons.SpikeTimes},vecStimOnTime,vecStimOffTime,vecFrameIdx);
		%clear other data
		clear sAggNeuron;
		clear sAggStim;
		
		%events
		dblStartEpoch = vecOrigStimOnTime(1)-10;
		dblStopEpoch = vecOrigStimOnTime(1)+dblStimDur+10;
		dblEpochDur = dblStopEpoch - dblStartEpoch;
		intNumN = numel(cellRawSpikeTimes);
		
		%remove spikes outside epoch
		for intN=1:intNumN
			vecT = cellRawSpikeTimes{intN};
			indRem = (vecT > dblStopEpoch) | (vecT < dblStartEpoch);
			cellRawSpikeTimes{intN} = unique(vecT(~indRem));
		end
		
		%% pool spikes from all neurons, but save the time+id per spike, then calculate IFR over all spikes at pop level
		intTotSpikeNum = sum(cellfun(@numel,cellRawSpikeTimes));
		intSpikeType = 2;
		if intSpikeType == 1
			%generate poisson-process spikes
			cellSpikeTimes = cell(1,intNumN);
			for intN=1:intNumN
				dblT0 = cellRawSpikeTimes{intN}(1);
				dblTN = dblStopEpoch;
				dblTotT = dblTN-dblT0;
				%add 3ms refractory period
				dblRt = 0;%(3/1000);
				dblLambda = numel(cellSpikeTimes{intN})/(dblTotT-dblRt*numel(cellSpikeTimes{intN}));
				vecISI = dblRt+exprnd(1./dblLambda,[1,round(dblLambda*dblTotT*2)]);
				vecSpikeT = dblT0+cumsum(vecISI)-vecISI(1);
				vecSpikeT = vecSpikeT(vecSpikeT<dblTN);
				
				%remove times outside epoch
				indRem = (vecSpikeT > dblStopEpoch) | (vecSpikeT < dblStartEpoch);
				cellSpikeTimes{intN} = vecSpikeT(~indRem);
				
				%swap time
				cellSpikeTimes{intN} = sort(abs(cellSpikeTimes{intN}-dblTN)+dblT0);
			end
			
			%get spikes per trial per neuron
			intTotS_Poiss = sum(cellfun(@numel,cellSpikeTimes));
			vecAllSpikeTime = nan(1,intTotS_Poiss);
			vecAllSpikeNeuron = zeros(1,intTotS_Poiss,'int16');
			intSP = 1;
			for intN=1:intNumN
				%poisson
				intThisSP = numel(cellSpikeTimes{intN});
				vecSpikeT_Poiss = cellSpikeTimes{intN};
				vecAllSpikeTime(intSP:(intSP+intThisSP-1)) = vecSpikeT_Poiss;
				vecAllSpikeNeuron(intSP:(intSP+intThisSP-1)) = intN;
				intSP = intSP + intThisSP;
			end
		else
			%random events
			vecAllSpikeTime = (rand(1,intTotSpikeNum)*dblEpochDur)+dblStartEpoch;
			vecAllSpikeNeuron = ones(size(vecAllSpikeTime));
			cellSpikeTimes = {vecAllSpikeTime};
		end
		intNeuronNum = numel(cellSpikeTimes);
		return
		%% remove spikes outside epoch
		indRem = (vecAllSpikeTime > dblStopEpoch) | (vecAllSpikeTime < dblStartEpoch);
		vecAllSpikeNeuron(indRem) = [];
		vecAllSpikeTime(indRem) = [];
		[vecAllSpikeTime,vecReorder] = sort(vecAllSpikeTime);
		vecAllSpikeNeuron = vecAllSpikeNeuron(vecReorder);
		
		% transform time indices
		%events
		dblStartEpoch = vecOrigStimOnTime(1)-dblStimDur;
		dblEpochDur = vecOrigStimOnTime(end)-vecOrigStimOnTime(1)+dblStimDur;
		dblStopEpoch = dblStartEpoch + dblEpochDur;
		
		% replace IFR with binned rates?
		intRateType = 1;
		if intRateType==1
			%binned counts
			dblBinDur = (0.5/1000);
			vecBins=(dblStartEpoch:dblBinDur:dblStopEpoch)';
			vecIFR = flat(histcounts(vecAllSpikeTime,vecBins)./dblBinDur);
			vecTime = vecBins(2:end)-dblBinDur/2;
		else
			%IFR
			[vecTime,vecIFR] = getIFR(vecAllSpikeTime,dblStartEpoch,dblEpochDur,0,[],[],0); %takes about 1 minute
			vecTime = vecTime + dblStartEpoch(1);
		end
		
		% peaks
		[vecPeakHeight,vecPeakLocs,w,p] = findpeaks(vecIFR);
		
		%threshold peaks
		dblCutOff = -inf;%0.65;
		vecStimPeakLocs = vecPeakLocs(vecPeakHeight>dblCutOff);
		
		%get raw pop events
		vecPopEventTimes = vecTime(vecStimPeakLocs);
		
		% what does peak look like? plot PSTH
		%raw peaks
		dblStep = (1/30000);
		vecEventBins = (-0.01+dblStep/2):dblStep:(0.01-dblStep/2);
		[dummy,dummy,vecEventBinsC,matPET] = doPEP(cellSpikeTimes,vecEventBins,vecPopEventTimes);
		vecEventBinsC = vecEventBinsC*1000;%ms
		matPopRate = sum(matPET,3);
		vecMean = mean(matPopRate,1);
		vecSEM = std(matPopRate,[],1)./sqrt(size(matPopRate,1));
		vecMean(vecEventBinsC==0)=nan;
		
		%raw random peak times
		vecRandPopEventTimes = linspace(min(vecPopEventTimes)-1,max(vecPopEventTimes)+1,numel(vecPopEventTimes));
		[dummy,dummy,dummy,matPETR] = doPEP(cellSpikeTimes,vecEventBins,vecRandPopEventTimes);
		matPopRateR = sum(matPETR,3);
		vecMeanR = mean(matPopRateR,1);
		vecSEMR = std(matPopRateR,[],1)./sqrt(size(matPopRateR,1));
		
		% plot
		figure
		plot(vecEventBinsC,vecMeanR,'color','k')
		hold on
		plot(vecEventBinsC,vecMean,'color',lines(1))
		hold off
		ylim([0 3000]);
		title(sprintf('%s',strRec),'interpreter','none');
		ylabel('Pop rate (Hz)');
		xlabel('Time after pop event (ms)');
		drawnow;
		fixfig;
		return
	end
end
toc
