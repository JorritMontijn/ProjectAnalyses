%% aim
%{
https://www.nature.com/articles/s41598-021-95037-z
https://www.nature.com/articles/s42003-021-02437-y

q1: Can absence of stimulus information in initial period be explained by neurons responding only to
black/white patches in the first x ms?

q2: How precise are spike times in the natural movie repetitions?


%}

%% define qualifying data
strRunType = 'Sim'; %Sim or ABI or Npx?
strRunStim = 'DG';%DG or NM?

runHeaderPopTimeCoding;
if strcmp(strRunType,'ABI')
	runLoadABI;
elseif strcmp(strRunType,'Sim')
	intSelectCells = 100; %how many cells to keep?
	intBatchNr=1;
	%vecSubSelectIdx = sort(randperm(1200,intSelectCells));
	vecSubSelectIdx = (intBatchNr-1)*intSelectCells+(1:intSelectCells);
	runLoadSim;
else
	runLoadNpx;
end
intAreas = numel(cellUseAreas{1});
cellSupraGranuInfra = {'Supragranular','Granular','Infragranular',''};

%% go through recordings
tic
for intRec=1:intRecNum
	%% prep ABI or Npx data
	if strcmp(strRunType,'ABI')
		runRecPrepABI;
		strThisRec = strRec;
	elseif strcmp(strRunType,'Sim')
		%not necessary
		
		%edit vars
		strThisRec = [strRec 'N' num2str(intSelectCells) 'B' num2str(intBatchNr)];
		strDataPathT0=strDataPathSimT0;
		
		%get layers
		cellSpikeTimesOrig = cellSpikeTimes;
		vecSupraGranuInfra = 3*ones(size(cellSpikeTimes));
		
	elseif strcmp(strRunType,'Npx')
		%prep
		runRecPrepNpx;
		strThisRec = strRec;
		strDataPathT0 = strTargetDataPath;
		
		%get layers
		cellSpikeTimesOrig = cellSpikeTimes;
		sArea1Neurons = sUseNeuron(indArea1Neurons);
		sRespNeurons = sArea1Neurons(indResp); %FR above 0.1 Hz
		cellAreas = {sRespNeurons.Area};
		vecCorticalLayer = cellfun(@(x) str2double(x(regexp(x,'layer.*')+6)),cellAreas);
		vecDepth = [sRespNeurons.DepthBelowIntersect];
		vecSupraGranuInfra = double(vecCorticalLayer < 4) + 2*double(vecCorticalLayer == 4) + 3*double(vecCorticalLayer > 4);
		
	end

	if intNeuronsInArea == 0% || intNeuronNum < 25
		fprintf('Number of neurons is %d for %s: skipping... [%s]\n',intNeuronNum,strRecOrig,getTime);
		continue;
	end
	
	%% get cortical depth per cell
	%scatter(vecSupraGranuInfra,vecDepth)
	%hold on
	%text(vecSupraGranuInfra,vecDepth,cellAreas)
	for intCortLayer = 4%1:3
		if intCortLayer == 4
		indUseNeurons = indTuned;
		strLayer = '';
		else
		indUseNeurons = vecSupraGranuInfra(:)==intCortLayer & indTuned;
		strLayer = cellSupraGranuInfra{intCortLayer};
		end
	if sum(indUseNeurons) < 5,continue;end
	cellSpikeTimes = cellSpikeTimesOrig(indUseNeurons);
	cellSpikeTimesReal = cellSpikeTimes;
	
	%events
	dblStartEpoch = vecStimOnTime(1)-10;
	dblStopEpoch = vecStimOnTime(end)+dblStimDur+10;
	dblEpochDur = dblStopEpoch - dblStartEpoch;
	intNumN = numel(cellSpikeTimesReal);
	dblTrialDur = median(vecStimOnTime(2:end)-vecStimOnTime(1:(end-1)));
	
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
		fprintf('Running %s - %s (%d/%d) [%s]\n',strRec,strType,intRec,intRecNum,getTime);
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
			%%
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
				cellOtherSpikes = cell(intNumN,1);
				
				%get beginning & end vectors
				for intN=1:intNumN
					%real
					[vecTrialPerSpike,vecTimePerSpike] = getSpikesInTrial(cellSpikeTimesReal{intN},vecStimOnTime);
					for intTrial=1:intTrialNum
						vecSpikeT = vecTimePerSpike(vecTrialPerSpike==intTrial);
						vecSpikeT(vecSpikeT>dblTrialDur)=[];
						cellSpikeTimesPerCellPerTrial{intN,intTrial} = vecSpikeT;
					end
					
					vecSpikes = cellSpikeTimesReal{intN};
					indUsedSpikes = false(size(vecSpikes));
					for intTrial=intTrialNum:-1:1
						%get spikes
						indTrialSpikes = ~indUsedSpikes & (vecSpikes >= vecStimOnTime(intTrial)) & vecSpikes < (vecStimOnTime(intTrial)+ dblTrialDur);
						indUsedSpikes(indTrialSpikes) = true;
						vecTheseSpikes = vecSpikes(indTrialSpikes);
						vecTheseSpikes = vecTheseSpikes - vecStimOnTime(intTrial);
						cellSpikeTimesPerCellPerTrial{intN,intTrial} = vecTheseSpikes;
					end
					cellOtherSpikes{intN} = vecSpikes(~indUsedSpikes);
				end
				
				%shuffle trial ids only within stim type if DG
				vecStimIdx = val2idx(vecStimIdx);
				intStimNum = numel(unique(vecStimIdx));
				cellUseSpikeTimesPerCellPerTrial = cell(size(cellSpikeTimesPerCellPerTrial));
				for intStimType=1:intStimNum
					vecOrigTrials = find(vecStimIdx==intStimType);
					for intN=1:intNumN
						vecShuffTrials = vecOrigTrials(randperm(numel(vecOrigTrials)));
						cellUseSpikeTimesPerCellPerTrial(intN,vecShuffTrials) = cellSpikeTimesPerCellPerTrial(intN,vecOrigTrials);
					end
				end
				
				%re-add trial onsets
				vecOrigTrials = find(vecStimIdx==intStimType);
				for intTrial=1:intTrialNum
					dblOnset = vecStimOnTime(intTrial);
					for intN=1:intNumN
						cellUseSpikeTimesPerCellPerTrial{intN,intTrial} = cellUseSpikeTimesPerCellPerTrial{intN,intTrial} + dblOnset;
					end
				end
				
				%compile
				for intN=1:intNumN
					cellSpikeTimes{intN} = unique(sort([cellOtherSpikes{intN}(:);...
						cell2vec(cellUseSpikeTimesPerCellPerTrial(intN,:))]));
				end
			else
				error
			end
		elseif strcmp(strType,'Shuff')
			%shuffled ISI per neuron
			for intN=1:intNumN
				%real
				vecSpikeT = cellSpikeTimesReal{intN};
				vecISI = diff(sort(vecSpikeT(:)));
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
		save(fullpath(strTargetDataPath,sprintf('T0Data_%s%s%s%s%s',strThisRec,strRunType,strRunStim,strType,strLayer)),...
			...%epoch
			'dblStartEpoch',...
			'dblEpochDur',...
			...%type
			'strType',...
			'strLayer',...
			...%data
			'vecAllSpikeTime',...
			'vecAllSpikeNeuron',...
			'vecTime',...
			'vecIFR',...
			'strRunStim'...
			);
	end
	end
end
toc
