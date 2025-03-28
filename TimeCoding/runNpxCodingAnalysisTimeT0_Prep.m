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
intRunDataType = 3;
strRunStim = 'NS';%DG or NM or NS? => superseded to WS by SWN
cellTypes = {'Real','Poiss','ShuffTid','Shuff','PoissGain','Uniform','RandTid','RandTxClass'};
boolFixSpikeGroupSize = false;
dblRemOnset = 0; %remove onset period in seconds; 0.125 for sim, 0.25 for npx
runHeaderPopTimeCoding;

%% go through recordings
intAreas = numel(cellUseAreas{1});
cellSupraGranuInfra = {'Supragranular','Granular','Infragranular',''};
tic
for intRec=1:intRecNum
	%% prep ABI or Npx data
	fprintf('Running rec %d/%d [%s]\n',intRec,intRecNum,getTime);
	if strcmp(strRunType,'ABI')
		runRecPrepABI;
		strThisRec = strRec;
		cellSpikeTimesOrig = cellSpikeTimes;
		vecSupraGranuInfra = 3*ones(size(cellSpikeTimes));
		
		strDataPathT0 = strTargetDataPath;
		
		%% is cell an interneuron (fast/narrow spiking) or pyramid (regular/broad spiking)?
		fprintf('Loading and processing waveforms... [%s]\n',getTime);
		sThisRec = sRec;
		sUseNeuron = sRec.sNeuron;
		
		%calculate waveform properties
		dblSampRateIM = 3e4;
		dblRecDur = max(cellfun(@max,{sUseNeuron.SpikeTimes})) - min(cellfun(@min,{sUseNeuron.SpikeTimes}));
		vecSpikeRate = cellfun(@numel,{sUseNeuron.SpikeTimes})/dblRecDur;
		matAreaWaveforms = cell2mat({sUseNeuron.Waveform})'; %[cell x sample]
		intNeurons=size(matAreaWaveforms,1);
		vecSpikeDur = nan(1,intNeurons);
		vecSpikePTR = nan(1,intNeurons);
		for intNeuron=1:intNeurons
			%find trough
			[dblTroughVal,intTrough]=min(matAreaWaveforms(intNeuron,:));
			[dblPeakVal,intTroughToPeak]=max(matAreaWaveforms(intNeuron,intTrough:end));
			intPeak = intTrough + intTroughToPeak - 1;
			
			dblTroughTime = intTrough/dblSampRateIM;
			dblTroughToPeakTime = intTroughToPeak/dblSampRateIM;
			dblPeakTime = intPeak/dblSampRateIM;
			vecSpikeDur(intNeuron) = dblTroughToPeakTime;
			vecSpikePTR(intNeuron) = abs(dblPeakVal/dblTroughVal);
		end
		dblRatio_PTT = 0.4;
		dblDur_SWT = 0.4/1000;
		vecNarrow = vecSpikeDur < dblDur_SWT & vecSpikePTR > dblRatio_PTT; %Ctx BL6
		vecBroad = vecSpikeDur > dblDur_SWT;% & vecSpikePTR < dblRatio_PTT; %Ctx BL6
		vecOther = ~vecNarrow & ~vecBroad;
		vecNeuronType = vecNarrow + vecBroad*2 + vecOther*3;
		vecNeuronSpikeDur = vecSpikeDur;
		vecNeuronPeakToTrough = vecSpikePTR;
		matNeuronWaveform = matAreaWaveforms;
		cellNeuronTypes = {'Narrow','Broad','Other'};
		
		%% plot
		vecSpikePTR(vecSpikePTR>1.25)=1.25;
		matCol = [0 0 1;0 1 0;0.5 0.5 0.5];
		figure;hold on
		colormap(matCol);
		scatter(vecSpikePTR(vecNeuronType==1),vecSpikeDur(vecNeuronType==1)*1000,[],matCol(1,:));
		scatter(vecSpikePTR(vecNeuronType==2),vecSpikeDur(vecNeuronType==2)*1000,[],matCol(2,:));
		scatter(vecSpikePTR(vecNeuronType==3),vecSpikeDur(vecNeuronType==3)*1000,[],matCol(3,:));
		xlabel('trough to peak ratio');
		ylabel('trough to peak dur (ms)');
		title(strRec,'interpreter','none');
		legend(cellNeuronTypes,'location','best');
		fixfig;
		
		%%
		export_fig(fullpath(strFigurePathSR,sprintf('T0_NarrowBroadClassification_%s.tif',strThisRec)));
		export_fig(fullpath(strFigurePathSR,sprintf('T0_NarrowBroadClassification_%s.pdf',strThisRec)));
		
		%% close fig
		close;
	elseif strcmp(strRunType,'Sim')
		%load
		runRecPrepSim;
		
		%edit vars
		strThisRec = strRec;
		strDataPathT0=strDataPathSimT0;
		vecOri180 = mod(vecOrientation,180)*2;
		vecStimIdx = vecOri180;
		
		%% get data matrix
		[matData,indTuned,cellSpikeTimes,sOut,cellSpikeTimesPerCellPerTrial,indResp] = ...
			SimPrepData(cellSpikeTimesRaw,vecStimOnTime,vecStimOffTime,vecOrientation);
		intTunedN = sum(indTuned);
		intRespN = size(matData,1);
		intNumN = numel(cellSpikeTimes);
		dblStimDur = roundi(median(vecStimOffTime - vecStimOnTime),1,'ceil');
		dblPreTime = 0;%0.3;
		dblPostTime = 0;%0.3;
		dblMaxDur = dblStimDur+dblPreTime+dblPostTime;
		intTrialNum = numel(vecStimOnTime);
		
		
		% move onset
		vecStimOnTime = vecStimOnTime + dblRemOnset;
		
		%get layers
		cellSpikeTimesOrig = cellSpikeTimes;
		vecSupraGranuInfra = 3*ones(size(cellSpikeTimes));
		
		%inerneuron/pyramid
		vecExcInh = [sLoad.sNeuron.Types];
		vecSwap = [2 1];%to match Npx
		vecNeuronType = vecSwap(vecExcInh);
		vecNeuronSpikeDur = [];
		vecNeuronPeakToTrough = [];
		matNeuronWaveform = [];
		cellNeuronTypes = {'Narrow','Broad','Other'};
		
		
	elseif strcmp(strRunType,'Npx') || strcmp(strRunType,'SWN')
		%prep
		runRecPrepNpx;
		if strcmp(strRunType,'SWN')
			indTuned = true(size(indTuned));
		end
		strThisRec = strRec;
		strDataPathT0 = strTargetDataPath;
		
		
		% move onset
		vecStimOnTime = vecStimOnTime + dblRemOnset;
		
		%get layers
		cellSpikeTimesOrig = cellSpikeTimes;
		sArea1Neurons = sUseNeuron(indArea1Neurons);
		sRespNeurons = sArea1Neurons(indResp); %FR above 0.1 Hz
		cellAreas = {sRespNeurons.Area};
		vecCorticalLayer = cellfun(@(x) str2double(x(regexp(x,'layer.*')+6)),cellAreas);
		vecDepth = [sRespNeurons.DepthBelowIntersect];
		vecSupraGranuInfra = double(vecCorticalLayer < 4) + 2*double(vecCorticalLayer == 4) + 3*double(vecCorticalLayer > 4);
	
		%% is cell an interneuron (fast/narrow spiking) or pyramid (regular/broad spiking)?
		fprintf('Loading and processing waveforms... [%s]\n',getTime);
		%load waveform
		if ~isfield(sUseNeuron,'Waveform')
			[sThisRec,sUseNeuron] = loadWaveforms(sThisRec,sUseNeuron);
		else
			sThisRec.sample_rate = str2double(sAggSources(intRec).sMetaAP.imSampRate);
		end
		
		%calculate waveform properties
		dblSampRateIM = sThisRec.sample_rate;
		dblRecDur = max(cellfun(@max,{sUseNeuron.SpikeTimes})) - min(cellfun(@min,{sUseNeuron.SpikeTimes}));
		vecSpikeRate = cellfun(@numel,{sUseNeuron.SpikeTimes})/dblRecDur;
		matAreaWaveforms = cell2mat({sUseNeuron.Waveform}'); %[cell x sample]
		intNeurons=size(matAreaWaveforms,1);
		vecSpikeDur = nan(1,intNeurons);
		vecSpikePTR = nan(1,intNeurons);
		for intNeuron=1:intNeurons
			%find trough
			[dblTroughVal,intTrough]=min(matAreaWaveforms(intNeuron,:));
			[dblPeakVal,intTroughToPeak]=max(matAreaWaveforms(intNeuron,intTrough:end));
			intPeak = intTrough + intTroughToPeak - 1;
			
			dblTroughTime = intTrough/dblSampRateIM;
			dblTroughToPeakTime = intTroughToPeak/dblSampRateIM;
			dblPeakTime = intPeak/dblSampRateIM;
			vecSpikeDur(intNeuron) = dblTroughToPeakTime;
			vecSpikePTR(intNeuron) = abs(dblPeakVal/dblTroughVal);
		end
		dblRatio_PTT = 0.4;
		dblDur_SWT = 0.4/1000;
		vecNarrow = vecSpikeDur < dblDur_SWT & vecSpikePTR > dblRatio_PTT; %Ctx BL6
		vecBroad = vecSpikeDur > dblDur_SWT;% & vecSpikePTR < dblRatio_PTT; %Ctx BL6
		vecOther = ~vecNarrow & ~vecBroad;
		vecNeuronType = vecNarrow + vecBroad*2 + vecOther*3;
		vecNeuronSpikeDur = vecSpikeDur;
		vecNeuronPeakToTrough = vecSpikePTR;
		matNeuronWaveform = matAreaWaveforms;
		cellNeuronTypes = {'Narrow','Broad','Other'};
		
		%% plot
		vecSpikePTR(vecSpikePTR>1.25)=1.25;
		matCol = [0 0 1;0 1 0;0.5 0.5 0.5];
		figure;hold on
		colormap(matCol);
		scatter(vecSpikePTR(vecNeuronType==1),vecSpikeDur(vecNeuronType==1)*1000,[],matCol(1,:));
		scatter(vecSpikePTR(vecNeuronType==2),vecSpikeDur(vecNeuronType==2)*1000,[],matCol(2,:));
		scatter(vecSpikePTR(vecNeuronType==3),vecSpikeDur(vecNeuronType==3)*1000,[],matCol(3,:));
		xlabel('trough to peak ratio');
		ylabel('trough to peak dur (ms)');
		title(strRec,'interpreter','none');
		legend(cellNeuronTypes,'location','best');
		fixfig;
		%%
		export_fig(fullpath(strFigurePathSR,sprintf('T0_NarrowBroadClassification_%s.tif',strThisRec)));
		export_fig(fullpath(strFigurePathSR,sprintf('T0_NarrowBroadClassification_%s.pdf',strThisRec)));
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
			indUseNeurons = true(size(find(indTuned)));
			strLayer = '';
		else
			indUseNeurons = vecSupraGranuInfra(:)==intCortLayer & true(size(find(indTuned)));
			strLayer = cellSupraGranuInfra{intCortLayer};
		end
		if sum(indUseNeurons) < 5
			fprintf('   Number of cells (%d) is too low; skipping\n...',sum(indUseNeurons));
			continue;
		end
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
			cellSpikeTimesReal{intN} = getUniqueSpikes(vecT(~indRem));
		end
		
		%% pool spikes from all neurons, but save the time+id per spike, then calculate IFR over all spikes at pop level
		for intType=[numel(cellTypes):-1:1]
			%which type?
			cellSpikeTimes = cell(1,intNumN);
			strType = cellTypes{intType};
			fprintf('Running %s - %s (%d/%d) [%s]\n',strRec,strType,intRec,intRecNum,getTime);
			if strcmp(strType,'Real')
				%real data
				cellSpikeTimes = cellSpikeTimesReal;
			elseif strcmp(strType,'Uniform')
				[cellUseSpikeTimesPerCellPerTrial,cellSpikeTimes] = buildUniformSpikes(cellSpikeTimesReal,vecStimOnTime,vecStimIdx,dblStimDur);
			elseif strcmp(strType,'Poiss')
				%poisson process neurons
				cellSpikeTimes = cell(1,intNumN);
				for intN=1:intNumN
					%add 3ms refractory period
					dblRt = (3/1000);
					dblLambda = numel(cellSpikeTimesReal{intN})/(dblEpochDur-dblRt*numel(cellSpikeTimesReal{intN}));
					vecISI = dblRt+exprnd(1./dblLambda,[1,round(dblLambda*dblEpochDur*2)]);
					vecSpikeT = dblStartEpoch+cumsum(vecISI)-vecISI(1);
					cellSpikeTimes{intN} = vecSpikeT(vecSpikeT<dblStopEpoch);
				end
				
			elseif strcmp(strType,'ShuffTid')
				%%
				cellSpikeTimes = cell(1,intNumN);
				if strcmp(strRunStim,'NM')
					error to do
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
					
				elseif strcmp(strRunStim,'DG') || strcmp(strRunStim,'WS') || strcmp(strRunStim,'NS')
					%create data
					[cellUseSpikeTimesPerCellPerTrial,cellSpikeTimes] = buildShuffTidSpikes(cellSpikeTimesReal,vecStimOnTime,vecStimIdx,dblTrialDur);
				
				else
					error
				end
			elseif strcmp(strType,'RandTid')
				%%
				cellSpikeTimes = cell(1,intNumN);
				if strcmp(strRunStim,'NM')
					error
				elseif strcmp(strRunStim,'DG') || strcmp(strRunStim,'WS') || strcmp(strRunStim,'NS')
					%create data
					[cellUseSpikeTimesPerCellPerTrial,cellSpikeTimes] = buildRandTidSpikes(cellSpikeTimesReal,vecStimOnTime,vecStimIdx,dblTrialDur);
				
				else
					error
				end
			elseif strcmp(strType,'RandTxClass')
				%%
				cellSpikeTimes = cell(1,intNumN);
				if strcmp(strRunStim,'NM')
					error
				elseif strcmp(strRunStim,'DG') || strcmp(strRunStim,'WS') || strcmp(strRunStim,'NS')
					%create data
					vecUseStimIdx = ones(size(vecStimIdx)); %collapse all classes
					[cellUseSpikeTimesPerCellPerTrial,cellSpikeTimes] = buildRandTidSpikes(cellSpikeTimesReal,vecStimOnTime,vecUseStimIdx,dblTrialDur);
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
					dblLambda = numel(cellSpikeTimesReal{intN})/(dblEpochDur-dblRt*numel(cellSpikeTimesReal{intN}));
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
			[vecTime,vecIFR] = getIFR(vecSpikeTimes,vecEventStarts,dblUseMaxDur); %takes about 1 minute
			vecTime = vecTime + dblStartEpoch(1);
			
			%% save intermediate data
			save(fullpath(strDataPathT0,sprintf('T0Data_%s%s%s%s%s%s',strThisRec,strRunType,strRunStim,strType,strLayer,strOnset)),...
				...%epoch
				'dblStartEpoch',...
				'dblEpochDur',...
				...%type
				'strType',...
				'strLayer',...
				...%stim
				'vecStimOnTime',...
				'vecStimIdx',...
				'dblStimDur',...
				...%data
				'vecAllSpikeTime',...
				'vecAllSpikeNeuron',...
				'cellSpikeTimes',...
				'vecTime',...
				'vecIFR',...
				'strRunStim',...
				'dblRemOnset',...
				'strRunStim',...
				'vecSupraGranuInfra',...
				...%neuron type
				'vecNeuronType',...
				'vecNeuronSpikeDur',...
				'vecNeuronPeakToTrough',...
				'matNeuronWaveform',...
				'cellNeuronTypes'...
				);
		end
	end
end
toc
