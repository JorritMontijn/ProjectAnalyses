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
		
		dblCutOff = 0.5;%sLoad.dblCutOff;
		vecTime_Real = sLoad.vecTime_Real;
		vecIFR_Real = sLoad.vecIFR_Real;
		vecPeakHeight_Real = sLoad.vecPeakHeight_Real;
		vecPeakHeight_Shuff = sLoad.vecPeakHeight_Shuff;
		vecPeakLocs_Real = sLoad.vecPeakLocs_Real;
		vecPeakLocs_Shuff = sLoad.vecPeakLocs_Shuff;
		vecNormTime_Real = sLoad.vecNormTime_Real;
		vecNormTime_Shuff = sLoad.vecNormTime_Shuff;
		vecAllSpikeNeuron_Real = sLoad.vecAllSpikeNeuron_Real;
		vecAllSpikeTime_Real = sLoad.vecAllSpikeTime_Real;
		
		%events
		dblStartEpoch = vecOrigStimOnTime(1)-10;
		dblEpochDur = vecOrigStimOnTime(end)-vecOrigStimOnTime(1)+dblStimDur+10;
		
		%retain only stim epoch
		indStimEpoch = vecAllSpikeTime_Real > dblStartEpoch & vecAllSpikeTime_Real < (dblStartEpoch + dblEpochDur);
		vecPopSpikeTime_Real = vecAllSpikeTime_Real(indStimEpoch);
		vecPopNeuronId_Real = vecAllSpikeNeuron_Real(indStimEpoch);
		
		%get event times
		vecPopEventTimes_Real = vecNormTime_Real(vecPeakLocs_Real(vecPeakHeight_Real>dblCutOff));
		intPopEventNum = numel(vecPopEventTimes_Real);
		vecPopEventLocs_Real = nan(1,intPopEventNum);
		vecPopEventLocsIFR_Real = nan(1,intPopEventNum);
		
		for intEvent=1:intPopEventNum
			[dummy,vecPopEventLocs_Real(intEvent)] = min(abs(vecPopEventTimes_Real(intEvent)-vecPopSpikeTime_Real));
			[dummy,vecPopEventLocsIFR_Real(intEvent)] = min(abs(vecPopEventTimes_Real(intEvent)-vecTime_Real));
		end
		
		% plot
		figure
		subplot(2,3,1)
		plot(vecTime_Real,vecIFR_Real);
		hold on
		scatter(vecPopEventTimes_Real,vecIFR_Real(vecPopEventLocsIFR_Real));
		scatter(vecOrigStimOnTime,ones(size(vecOrigStimOnTime)),'kx');
		
		vecTrialBins = 0:0.25:dblStimDur;
		[vecMean,vecSEM,vecWindowBinCenters,matPET] = doPEP(cellSpikeTimes,vecTrialBins,vecOrigStimOnTime,-1);
		
		subplot(2,3,2)
		imagesc(mean(matPET,3));
		
		subplot(2,3,4)
		imagesc(matPET(:,:,1));
		subplot(2,3,5)
		imagesc(matPET(:,:,2));
		subplot(2,3,6)
		imagesc(matPET(:,:,3));
		return
	end
end
toc
