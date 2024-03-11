%% aim
%{
https://www.nature.com/articles/s41598-021-95037-z
https://www.nature.com/articles/s42003-021-02437-y

q1: Can absence of stimulus information in initial period be explained by neurons responding only to
black/white patches in the first x ms?

q2: How precise are spike times in the natural movie repetitions?


%}

%% define qualifying data
strRunType = 'Npx'; %ABI or Npx?
strRunStim = 'DG';%DG or NM?

runHeaderPopTimeCoding;
if strcmp(strRunType,'ABI')
	runLoadABI;
else
	runLoadNpx;
end
intAreas = numel(cellUseAreas{1});

%% go through recordings
tic
for intRec=1:intRecNum
	%% prep ABI or Npx data
	if strcmp(strRunType,'ABI')
		runRecPrepABI;
	elseif strcmp(strRunType,'Npx')
		runRecPrepNpx;
	end
	if intNeuronsInArea == 0 || intNeuronNum < 25
		fprintf('Number of neurons is %d for %s: skipping... [%s]\n',intNeuronNum,strRecOrig,getTime);
		continue;
	end
	cellSpikeTimesReal = cellSpikeTimes;
	
	%events
	dblStartEpoch = vecStimOnTime(1)-10;
	dblStopEpoch = vecStimOnTime(end)+dblStimDur+10;
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
	for intType=[5:-1:1]
		%which type?
		cellSpikeTimes = cell(1,intNumN);
		strType = cellTypes{intType};
		fprintf('Running %s - %s (%d/%d) [%s]\n',strRec,strType,intRec,numel(sAggStim),getTime);
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
			cellSpikeTimes = cell(1,intNumN);
			if strcmp(strRunStim,'NM')
				
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
					
			elseif strcmp(strRunStim,'DG')
				%create normal data
				cellSpikeTimesPerCellPerTrial = cell(intNumN,intTrialNum);
				cellStartSpikes = cell(intNumN,1);
				cellEndSpikes = cell(intNumN,1);
				
				%get beginning & end vectors
				for intN=1:intNumN
					%real
					[vecTrialPerSpike,vecTimePerSpike] = getSpikesInTrial(cellSpikeTimesReal{intN},vecStimOnTime,dblStimDur);
					for intTrial=1:intTrialNum
						vecSpikeT = vecTimePerSpike(vecTrialPerSpike==intTrial);
						cellSpikeTimesPerCellPerTrial{intN,intTrial} = vecSpikeT;
					end
					cellStartSpikes{intN} = cellSpikeTimesReal{intN}(cellSpikeTimesReal{intN}<vecStimOnTime(1));
					cellEndSpikes{intN} = cellSpikeTimesReal{intN}(cellSpikeTimesReal{intN}>(vecStimOnTime(end)+dblStimDur));
				end
				
				%shuffle trial ids only within stim type if DG
				cellUseSpikeTimesPerCellPerTrial = cell(size(cellSpikeTimesPerCellPerTrial));
				for intStimType=1:intStimNum
					vecOrigTrials = find(vecStimIdx==intStimType);
					for intN=1:intNumN
						vecShuffTrials = vecOrigTrials(randperm(numel(vecOrigTrials)));
						cellUseSpikeTimesPerCellPerTrial(intN,vecShuffTrials) = cellSpikeTimesPerCellPerTrial(intN,vecOrigTrials);
					end
				end
				
				%compile
				for intN=1:intNumN
					cellSpikeTimes{intN} = unique(sort([cellStartSpikes{intN};...
						cell2vec(cellUseSpikeTimesPerCellPerTrial(intN,:));...
						cellEndSpikes{intN}]));
				end
			else
				error
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
			vecGainE = dblStartEpoch:dblWindow:dblStopEpoch;
			vecGainT = vecGainE(2:end) - dblWindow/2;
			
			matSpikeCounts = zeros(intNumN,numel(vecGainT));
			for intN=1:intNumN
				matSpikeCounts(intN,:) = histcounts(cellSpikeTimesReal{intN},vecGainE);
			end
			vecGainAx = mean(matSpikeCounts,2);
			vecGain=getProjOnLine(matSpikeCounts,vecGainAx)';
			vecGain = vecGain./norm(vecGainAx);
			
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
		vecSpikeTimes = vecAllSpikeTime;
		vecEventStarts = dblStartEpoch;
		dblUseMaxDur = dblEpochDur;
		intSmoothSd = 0;
		dblMinScale = [];
		dblBase = [];
		boolUseParallel = 0;
		[vecTime,vecIFR] = getIFR(vecSpikeTimes,vecEventStarts,dblUseMaxDur,intSmoothSd,dblMinScale,dblBase,boolUseParallel); %takes about 1 minute
		vecTime = vecTime + dblStartEpoch(1);
		
		%% save intermediate data
		save(fullpath(strTargetDataPath,sprintf('T0Data_%s%s%s%s',strRec,strRunType,strRunStim,strType)),...
			...%epoch
			'dblStartEpoch',...
			'dblEpochDur',...
			...%type
			'strType',...
			...%data
			'vecAllSpikeTime',...
			'vecAllSpikeNeuron',...
			'vecTime',...
			'vecIFR',...
			'strRunStim'...
			);
	end
end
toc
