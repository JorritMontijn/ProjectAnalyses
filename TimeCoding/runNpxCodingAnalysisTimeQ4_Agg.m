%% aim
%{
is population gain correlated with pupil size?
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
		[matMeanRate,cellSpikeTimes] = ...
			NpxPrepMovieData(cellSpikeTimesRaw,vecStimOnTime,vecStimOffTime,vecFrameIdx);
		
		%events
		dblStartEpoch = vecOrigStimOnTime(1)-10;
		dblStopEpoch = vecOrigStimOnTime(end)+dblStimDur+10;
		dblEpochDur = dblStopEpoch - dblStartEpoch;
		intNumN = numel(cellSpikeTimes);
		
		%% calc gain
		vecGainAx = mean(matMeanRate,2);
		%calc gain
		dblWindow = 1; %secs
		vecGainE = dblStartEpoch:dblWindow:dblStopEpoch;
		vecGainT = vecGainE(2:end) - dblWindow/2;
		matBinnedRate = nan(intNumN,numel(vecGainT));
		for intN=1:intNumN
			matBinnedRate(intN,:)= histcounts(cellSpikeTimes{intN},vecGainE)./dblWindow;
		end
		
		vecGain=getProjOnLine(matBinnedRate,vecGainAx)';
		vecNormGain = vecGain./mean(vecGain(:));
		vecMean=mean(matBinnedRate,1);
		vecNormMean = vecMean./mean(vecMean(:));
		
		%% are fluctuations in the direction of the mean-axis or gain-axis?
		
		%% is population gain correlated with pupil size?
			
		
		%get pupil size
		
		%% save data
		
	end
end
toc
