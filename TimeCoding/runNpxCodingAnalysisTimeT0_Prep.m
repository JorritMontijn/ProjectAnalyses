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
		cellSpikeTimesRaw = {sArea1Neurons.SpikeTimes};
		[matMeanRate,cellSpikeTimes] = ...
			NpxPrepMovieData(cellSpikeTimesRaw,vecStimOnTime,vecStimOffTime,vecFrameIdx);
		
		%events
		dblStartEpoch = vecOrigStimOnTime(1)-10;
		dblStopEpoch = vecOrigStimOnTime(end)+dblStimDur+10;
		dblEpochDur = dblStopEpoch - dblStartEpoch;
		intNumN = numel(cellSpikeTimes);
		
		%remove spikes outside epoch
		for intN=1:intNumN
			vecT = cellSpikeTimes{intN};
			indRem = (vecT > dblStopEpoch) | (vecT < dblStartEpoch);
			cellSpikeTimes{intN} = unique(vecT(~indRem));
		end
		
		%% pool spikes from all neurons, but save the time+id per spike, then calculate IFR over all spikes at pop level
		%shuffle trial ids
		cellSpikeTimes_ShuffTid = cell(size(cellSpikeTimes));
		cellSpikeTimesPerCellPerTrial = cell(intNumN,intOrigTrialNum);
		boolDiscardEdges = true;
		for intN=1:intNumN
			%real
			[vecTrialPerSpike,vecTimePerSpike] = getSpikesInTrial(cellSpikeTimes{intN},vecOrigStimOnTime,dblStimDur);
			for intTrial=1:intOrigTrialNum
				vecSpikeT = vecTimePerSpike(vecTrialPerSpike==intTrial);
				cellSpikeTimesPerCellPerTrial{intN,intTrial} = vecSpikeT;
			end
			
			%get beginning & end vectors
			vecStartSpikes = cellSpikeTimes{intN}(cellSpikeTimes{intN}<vecOrigStimOnTime(1));
			vecEndSpikes = cellSpikeTimes{intN}(cellSpikeTimes{intN}>(vecOrigStimOnTime(end)+dblStimDur));
			
			%shuffle trial id
			vecRandTrialIds = randperm(intOrigTrialNum);
			cellShuffTidTrials = cellSpikeTimesPerCellPerTrial(intN,:);
			for intTrial=1:intOrigTrialNum
				dblRandStart = vecOrigStimOnTime(vecRandTrialIds(intTrial));
				cellShuffTidTrials{intTrial} = cellShuffTidTrials{intTrial}+dblRandStart;
			end
			cellSpikeTimes_ShuffTid{intN} = unique(sort([vecStartSpikes;cell2vec(cellShuffTidTrials);vecEndSpikes]));
		end
		
		%generate poisson-process spikes
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
		intTotS_ShuffTid= sum(cellfun(@numel,cellSpikeTimes_ShuffTid));
		vecAllSpikeTime_Real = nan(1,intTotS);
		vecAllSpikeNeuron_Real = zeros(1,intTotS,'int16');
		vecAllSpikeTime_Poiss = nan(1,intTotS_Poiss);
		vecAllSpikeNeuron_Poiss = zeros(1,intTotS_Poiss,'int16');
		vecAllSpikeTime_ShuffTid = nan(1,intTotS_ShuffTid);
		vecAllSpikeNeuron_ShuffTid = zeros(1,intTotS_ShuffTid,'int16');
		vecAllSpikeTime_Shuff = nan(1,intTotS);
		vecAllSpikeNeuron_Shuff = zeros(1,intTotS,'int16');
		intS = 1;
		intSP = 1;
		intSS = 1;
		for intN=1:intNumN
			%add spikes
			intThisS = numel(cellSpikeTimes{intN});
			vecSpikeT = cellSpikeTimes{intN};
			vecAllSpikeTime_Real(intS:(intS+intThisS-1)) = vecSpikeT;
			vecAllSpikeNeuron_Real(intS:(intS+intThisS-1)) = intN;
			
			%poisson
			intThisSP = numel(cellSpikeTimes_Poiss{intN});
			vecSpikeT_Poiss = cellSpikeTimes_Poiss{intN};
			vecAllSpikeTime_Poiss(intSP:(intSP+intThisSP-1)) = vecSpikeT_Poiss;
			vecAllSpikeNeuron_Poiss(intSP:(intSP+intThisSP-1)) = intN;
			intSP = intSP + intThisSP;
			
			%shuffle trial id
			intThisSS = numel(cellSpikeTimes_ShuffTid{intN});
			vecSpikeT_ShuffTid = cellSpikeTimes_ShuffTid{intN};
			vecAllSpikeTime_ShuffTid(intSS:(intSS+intThisSS-1)) = vecSpikeT_ShuffTid;
			vecAllSpikeNeuron_ShuffTid(intSS:(intSS+intThisSS-1)) = intN;
			intSS = intSS + intThisSS;
			
			%shuffle isi
			vecISI = diff(sort(vecSpikeT));
			vecSpikeT_Shuff = cumsum([0;vecISI(randperm(numel(vecISI)))])+vecSpikeT(1);
			vecAllSpikeTime_Shuff(intS:(intS+intThisS-1)) = vecSpikeT_Shuff;
			vecAllSpikeNeuron_Shuff(intS:(intS+intThisS-1)) = intN;
			intS = intS + intThisS;
			
		end
		
		%generate variables
		sReal = struct;
		sPoiss = struct;
		sShuffTid = struct;
		sShuff = struct;
		for intType=1:4
			if intType==1
				vecAllSpikeTime = vecAllSpikeTime_Real;
				vecAllSpikeNeuron = vecAllSpikeNeuron_Real;
				strType = 'Real';
			elseif intType == 2
				vecAllSpikeTime = vecAllSpikeTime_Poiss;
				vecAllSpikeNeuron = vecAllSpikeNeuron_Poiss;
				strType = 'Poiss';
			elseif intType == 3
				vecAllSpikeTime = vecAllSpikeTime_ShuffTid;
				vecAllSpikeNeuron = vecAllSpikeNeuron_ShuffTid;
				strType = 'ShuffTid';
			else
				vecAllSpikeTime = vecAllSpikeTime_Shuff;
				vecAllSpikeNeuron = vecAllSpikeNeuron_Shuff;
				strType = 'Shuff';
			end
			
			%remove spikes outside epoch
			indRem = (vecAllSpikeTime > dblStopEpoch) | (vecAllSpikeTime < dblStartEpoch);
			vecAllSpikeNeuron(indRem) = [];
			vecAllSpikeTime(indRem) = [];
			
			%sort
			[vecAllSpikeTime,vecReorder] = sort(vecAllSpikeTime);
			vecAllSpikeNeuron = vecAllSpikeNeuron(vecReorder);
			[vecTime,vecIFR] = getIFR(vecAllSpikeTime,dblStartEpoch,dblEpochDur,0,[],[],0); %takes about 1 minute
			vecTime = vecTime + dblStartEpoch(1);
			
			%create variables
			eval(['s' strType '.vecAllSpikeTime = vecAllSpikeTime;']);
			eval(['s' strType '.vecAllSpikeNeuron = vecAllSpikeNeuron;']);
			eval(['s' strType '.vecTime = vecTime;']);
			eval(['s' strType '.vecIFR = vecIFR;']);
			
		end
		
		%% save intermediate data
		save(fullpath(strTargetDataPath,sprintf('T0Data_%s',strRec)),...
			...%epoch
			'dblStartEpoch',...
			'dblEpochDur',...
			...%real
			'sReal',...
			...%shuffled isi
			'sShuff',...
			...%shuffled trial id
			'sShuffTid',...
			...%poisson
			'sPoiss'...
			);
	end
end
toc
