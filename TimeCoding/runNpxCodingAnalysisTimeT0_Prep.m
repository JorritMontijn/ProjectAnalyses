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
		[matMeanRate,cellSpikeTimesReal] = ...
			NpxPrepMovieData(cellSpikeTimesRaw,vecStimOnTime,vecStimOffTime,vecFrameIdx);
		
		%events
		dblStartEpoch = vecOrigStimOnTime(1)-10;
		dblStopEpoch = vecOrigStimOnTime(end)+dblStimDur+10;
		dblEpochDur = dblStopEpoch - dblStartEpoch;
		intNumN = numel(cellSpikeTimesReal);
		
		%remove spikes outside epoch
		for intN=1:intNumN
			vecT = cellSpikeTimesReal{intN};
			indRem = (vecT > dblStopEpoch) | (vecT < dblStartEpoch);
			cellSpikeTimesReal{intN} = unique(vecT(~indRem));
		end
		
		%% pool spikes from all neurons, but save the time+id per spike, then calculate IFR over all spikes at pop level
		cellTypes = {'Real','Poiss','ShuffTid','Shuff','PoissGain'};
		for intType=5
			%which type?
			cellSpikeTimes = cell(1,intNumN);
			strType = cellTypes{intType};
			if strcmp(strType,'Real')
				%real data
				cellSpikeTimes = cellSpikeTimesReal;
			elseif strcmp(strType,'Poiss')
				%poisson process neurons
				cellSpikeTimes = cell(1,intNumN);
				for intN=1:intNumN
					dblT0 = cellSpikeTimesReal{intN}(1);
					dblTN = cellSpikeTimesReal{intN}(end);
					dblTotT = dblTN-dblT0;
					%add 3ms refractory period
					dblRt = (3/1000);
					dblLambda = numel(cellSpikeTimesReal{intN})/(dblTotT-dblRt*numel(cellSpikeTimesReal{intN}));
					vecISI = dblRt+exprnd(1./dblLambda,[1,round(dblLambda*dblTotT*2)]);
					vecSpikeT = dblT0+cumsum(vecISI)-vecISI(1);
					cellSpikeTimes{intN} = vecSpikeT(vecSpikeT<dblTN);
				end
				
			elseif strcmp(strType,'ShuffTid')
				%shuffle trial ids
				cellSpikeTimesPerCellPerTrial = cell(intNumN,intOrigTrialNum);
				for intN=1:intNumN
					%real
					[vecTrialPerSpike,vecTimePerSpike] = getSpikesInTrial(cellSpikeTimesReal{intN},vecOrigStimOnTime,dblStimDur);
					for intTrial=1:intOrigTrialNum
						vecSpikeT = vecTimePerSpike(vecTrialPerSpike==intTrial);
						cellSpikeTimesPerCellPerTrial{intN,intTrial} = vecSpikeT;
					end
					
					%get beginning & end vectors
					vecStartSpikes = cellSpikeTimesReal{intN}(cellSpikeTimesReal{intN}<vecOrigStimOnTime(1));
					vecEndSpikes = cellSpikeTimesReal{intN}(cellSpikeTimesReal{intN}>(vecOrigStimOnTime(end)+dblStimDur));
					
					%shuffle trial id
					vecRandTrialIds = randperm(intOrigTrialNum);
					cellShuffTidTrials = cellSpikeTimesPerCellPerTrial(intN,:);
					for intTrial=1:intOrigTrialNum
						dblRandStart = vecOrigStimOnTime(vecRandTrialIds(intTrial));
						cellShuffTidTrials{intTrial} = cellShuffTidTrials{intTrial}+dblRandStart;
					end
					cellSpikeTimes{intN} = unique(sort([vecStartSpikes;cell2vec(cellShuffTidTrials);vecEndSpikes]));
				end
				
			elseif strcmp(strType,'Shuff')
				%shuffled ISI per neuron
				for intN=1:intNumN
					%real
					vecSpikeT = cellSpikeTimesReal{intN};
					vecISI = diff(sort(vecSpikeT));
					vecSpikeT_Shuff = cumsum([0;vecISI(randperm(numel(vecISI)))])+vecSpikeT(1);
					cellSpikeTimes{intN} = vecSpikeT_Shuff;
				end
			elseif strcmp(strType,'PoissGain')
				%poisson process neurons with variable lambda matched to population-gain
				dblWindow = 1; %secs
				vecAllSpikes = cell2vec(cellSpikeTimesReal);
				vecGainE = dblStartEpoch:dblWindow:dblStopEpoch;
				vecGainT = vecGainE(2:end) - dblWindow/2;
				vecR = histcounts(vecAllSpikes,vecGainE);
				vecGain = vecR ./ mean(vecR);
				
				%poisson process
				cellSpikeTimes = cell(1,intNumN);
				for intN=1:intNumN
					dblT0 = cellSpikeTimesReal{intN}(1);
					dblTN = cellSpikeTimesReal{intN}(end);
					dblTotT = dblTN-dblT0;
					%add 3ms refractory period
					dblRt = (3/1000);
					dblLambda = numel(cellSpikeTimesReal{intN})/(dblTotT-dblRt*numel(cellSpikeTimesReal{intN}));
					vecSpikeT = nan(1,2*numel(cellSpikeTimesReal{intN}));
					vecSpikeT(1) = dblT0;
					dblLastSpike = dblT0;
					intSpikeIdx = 1;
					while dblLastSpike < dblTN
						%get current lambda
						intGainIdx = find(vecGainT>dblLastSpike,1);
						if isempty(intGainIdx),intGainIdx=numel(vecGainT);end
						dblCurrLambda = dblLambda*vecGain(intGainIdx);
						dblLastSpike = dblLastSpike + dblRt+exprnd(1./dblCurrLambda);
						intSpikeIdx = intSpikeIdx + 1;
						vecSpikeT(intSpikeIdx) = dblLastSpike;
					end
					cellSpikeTimes{intN} = vecSpikeT(vecSpikeT<dblTN);
				end
			end
			
			%get spikes per trial per neuron
			intTotS = sum(cellfun(@numel,cellSpikeTimes));
			vecAllSpikeTime = nan(1,intTotS);
			vecAllSpikeNeuron = zeros(1,intTotS,'int16');
			intS = 1;
			for intN=1:intNumN
				%add spikes
				intThisS = numel(cellSpikeTimes{intN});
				vecSpikeT = cellSpikeTimes{intN};
				vecAllSpikeTime(intS:(intS+intThisS-1)) = vecSpikeT;
				vecAllSpikeNeuron(intS:(intS+intThisS-1)) = intN;
				intS = intS + intThisS;
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
			
			%% save intermediate data
			save(fullpath(strTargetDataPath,sprintf('T0Data_%s%s',strRec,strType)),...
				...%epoch
				'dblStartEpoch',...
				'dblEpochDur',...
				...%type
				'strType',...
				...%data
				'vecAllSpikeTime',...
				'vecAllSpikeNeuron',...
				'vecTime',...
				'vecIFR'...
				);
		end
	end
end
toc
