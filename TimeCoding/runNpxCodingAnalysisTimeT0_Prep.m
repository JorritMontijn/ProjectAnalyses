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
		
		%%
		% pool spikes from all neurons, but save the time+id per spike, then calculate IFR over all
		% spikes at pop level, detect peaks, and split into trials.
		
		%generate poisson-process spikes
		intNumN = size(matMeanRate,1);
		cellSpikeTimes_Poiss = cell(1,intNumN);
		for intN=1:intNumN
			dblT0 = cellSpikeTimes{intN}(1);
			dblTN = cellSpikeTimes{intN}(end);
			dblTotT = dblTN-dblT0;
			dblLambda = numel(cellSpikeTimes{intN})/dblTotT;
			vecISI = exprnd(1./dblLambda,[1,round(dblLambda*dblTotT*2)]);
			vecSpikeT = dblT0+cumsum(vecISI)-vecISI(1);
			cellSpikeTimes_Poiss{intN} = vecSpikeT(vecSpikeT<dblTN);
		end
		
		%get spikes per trial per neuron
		intTotS = sum(cellfun(@numel,cellSpikeTimes));
		intTotS_Poiss = sum(cellfun(@numel,cellSpikeTimes_Poiss));
		vecAllSpikeTime_Real = nan(1,intTotS);
		vecAllSpikeNeuron_Real = zeros(1,intTotS,'int16');
		vecAllSpikeTime_Shuff = nan(1,intTotS);
		vecAllSpikeNeuron_Shuff = zeros(1,intTotS,'int16');
		vecAllSpikeTime_Poiss = nan(1,intTotS_Poiss);
		vecAllSpikeNeuron_Poiss = zeros(1,intTotS_Poiss,'int16');
		intNumN = size(matMeanRate,1);
		intS = 1;
		intSP = 1;
		%events
		dblStartEpoch = vecOrigStimOnTime(1)-10;
		dblEpochDur = vecOrigStimOnTime(end)-vecOrigStimOnTime(1)+dblStimDur+10;
		dblStopEpoch = dblStartEpoch+dblEpochDur;
		for intN=1:intNumN
			%add spikes
			intThisS = numel(cellSpikeTimes{intN});
			vecSpikeT = cellSpikeTimes{intN};
			vecAllSpikeTime_Real(intS:(intS+intThisS-1)) = vecSpikeT;
			vecAllSpikeNeuron_Real(intS:(intS+intThisS-1)) = intN;
			
			%shuffle
			vecISI = diff(sort(vecSpikeT));
			vecShuffSpikes = cumsum([0;vecISI(randperm(numel(vecISI)))])+vecSpikeT(1);
			vecAllSpikeTime_Shuff(intS:(intS+intThisS-1)) = vecShuffSpikes;
			vecAllSpikeNeuron_Shuff(intS:(intS+intThisS-1)) = intN;
			intS = intS + intThisS;
			
			%poisson
			intThisSP = numel(cellSpikeTimes_Poiss{intN});
			vecSpikeT_Poiss = cellSpikeTimes_Poiss{intN};
			vecAllSpikeTime_Poiss(intSP:(intSP+intThisSP-1)) = vecSpikeT_Poiss;
			vecAllSpikeNeuron_Poiss(intSP:(intSP+intThisSP-1)) = intN;
			intSP = intSP + intThisSP;
		end
		
		%remove spikes outside epoch
		indRemReal = (vecAllSpikeTime_Real > dblStopEpoch) | (vecAllSpikeTime_Real < dblStartEpoch);
		indRemShuff = (vecAllSpikeTime_Shuff > dblStopEpoch) | (vecAllSpikeTime_Shuff < dblStartEpoch);
		indRemPoiss = (vecAllSpikeTime_Poiss > dblStopEpoch) | (vecAllSpikeTime_Poiss < dblStartEpoch);
		vecAllSpikeNeuron_Real(indRemReal) = [];
		vecAllSpikeTime_Real(indRemReal) = [];
		vecAllSpikeTime_Shuff(indRemShuff) = [];
		vecAllSpikeNeuron_Shuff(indRemShuff) = [];
		vecAllSpikeTime_Poiss(indRemPoiss) = [];
		vecAllSpikeNeuron_Poiss(indRemPoiss) = [];
	
		%sort
		[vecAllSpikeTime_Real,vecReorder] = sort(vecAllSpikeTime_Real);
		vecAllSpikeNeuron_Real = vecAllSpikeNeuron_Real(vecReorder);
		[vecTime_Real,vecIFR_Real] = getIFR(vecAllSpikeTime_Real,dblStartEpoch,dblEpochDur,0,[],[],0); %takes about 1 minute
		vecTime_Real = vecTime_Real + dblStartEpoch(1);
		
		%shuffled
		[vecAllSpikeTime_Shuff,vecReorder] = sort(vecAllSpikeTime_Shuff);
		vecAllSpikeNeuron_Shuff = vecAllSpikeNeuron_Shuff(vecReorder);
		[vecTime_Shuff,vecIFR_Shuff] = getIFR(vecAllSpikeTime_Shuff,dblStartEpoch,dblEpochDur,0,[],[],0); %takes about 1 minute
		vecTime_Shuff = vecTime_Shuff + dblStartEpoch(1);
		
		%poisson
		[vecAllSpikeTime_Poiss,vecReorder] = sort(vecAllSpikeTime_Poiss);
		vecAllSpikeNeuron_Poiss = vecAllSpikeNeuron_Poiss(vecReorder);
		[vecTime_Poiss,vecIFR_Poiss] = getIFR(vecAllSpikeTime_Poiss,dblStartEpoch,dblEpochDur,0,[],[],0); %takes about 1 minute
		vecTime_Poiss = vecTime_Poiss + dblStartEpoch(1);
		
		%% save intermediate data
		save(fullpath(strTargetDataPath,sprintf('T0Data_%s',strRec)),...
			...%epoch
			'dblStartEpoch',...
			'dblEpochDur',...
			...%real
			'vecAllSpikeTime_Real',...
			'vecAllSpikeNeuron_Real',...
			'vecTime_Real',...
			'vecIFR_Real',...
			...%shuffled
			'vecAllSpikeTime_Shuff',...
			'vecAllSpikeNeuron_Shuff',...
			'vecTime_Shuff',...
			'vecIFR_Shuff',...
			...%shuffled
			'vecAllSpikeTime_Poiss',...
			'vecAllSpikeNeuron_Poiss',...
			'vecTime_Poiss',...
			'vecIFR_Poiss'...
			);
	end
end
toc
