function [cellSpikeTimesPerCellPerTrial,cellSpikeTimes] = buildShuffTidSpikes(cellSpikeTimesReal,vecStimOnTime,vecStimIdx,dblTrialDur)
	%buildShuffTidSpikes Summary of this function goes here
	%   [cellSpikeTimesPerCellPerTrial,cellSpikeTimes] = buildShuffTidSpikes(cellSpikeTimesReal,vecStimOnTime,vecStimIdx,dblTrialDur);
	
	%get beginning & end vectors
	intNumN = numel(cellSpikeTimesReal);
	intTrialNum = numel(vecStimOnTime);
	cellOtherSpikes = cell(intNumN,1);
	cellTrialSpikes = cell(intNumN,1);
	cellTrialTrials = cell(intNumN,1);
	cellShuffTrials = cell(intNumN,1);
	cellSpikeTimes = cellSpikeTimesReal;
	cellSpikeTimesPerCellPerTrial = cell(intNumN,intTrialNum);
	hTic=tic;
	for intN=1:intNumN
		%msg
		if toc(hTic) > 5
			fprintf('  ShuffTid %d/%d [%s]\n',intN,intNumN,getTime);
			hTic=tic;
		end
		
		%real
		vecSpikes = cellSpikeTimesReal{intN};
		vecTrialPerSpike = zeros(size(vecSpikes));
		vecTimePerSpike = zeros(size(vecSpikes));
		indUsedSpikes = false(size(vecSpikes));
		for intTrial=intTrialNum:-1:1
			%get spikes
			indTrialSpikes = ~indUsedSpikes & (vecSpikes >= vecStimOnTime(intTrial)) & vecSpikes < (vecStimOnTime(intTrial)+ dblTrialDur);
			indUsedSpikes(indTrialSpikes) = true;
			vecTheseSpikes = vecSpikes(indTrialSpikes);
			vecTheseSpikes = vecTheseSpikes - vecStimOnTime(intTrial);
			vecTimePerSpike(indTrialSpikes) = vecTheseSpikes;
			vecTrialPerSpike(indTrialSpikes) = intTrial;
		end
		cellOtherSpikes{intN} = vecSpikes(~indUsedSpikes);
		cellTrialSpikes{intN} = vecTimePerSpike;
		cellTrialTrials{intN} = vecTrialPerSpike;
	end
	
	%shuffle trial ids only within stim type if DG
	vecStimIdx = val2idx(vecStimIdx);
	intStimNum = numel(unique(vecStimIdx));
	for intN=1:intNumN
		cellShuffTrials{intN}(cellTrialTrials{intN}==0) = 0;
		for intStimType=1:intStimNum
			vecTrialsOfType = find(vecStimIdx==intStimType);
			indSelectSpikes = ismember(cellTrialTrials{intN},vecTrialsOfType);
			
			%assign to random trials: keep # of spikes in each trial the same; only shuffle
			vecOrigTrials = cellTrialTrials{intN}(indSelectSpikes);
			vecShuffTrials = vecOrigTrials(randperm(numel(vecOrigTrials)));
			cellShuffTrials{intN}(indSelectSpikes) = vecShuffTrials;
		end
	end
	
	%re-add trial onsets
	cellGlobalSpikeTimesPerCellPerTrial = cellSpikeTimesPerCellPerTrial;
	for intTrial=1:intTrialNum
		dblOnset = vecStimOnTime(intTrial);
		for intN=1:intNumN
			%vecSpikeT_orig = cellTrialSpikes{intN}(cellTrialTrials{intN}==intTrial);
			vecSpikeT = cellTrialSpikes{intN}(cellShuffTrials{intN}==intTrial);
			cellSpikeTimesPerCellPerTrial{intN,intTrial} = vecSpikeT;
			cellGlobalSpikeTimesPerCellPerTrial{intN,intTrial} = vecSpikeT + dblOnset;
		end
	end
	
	%compile
	for intN=1:intNumN
		cellSpikeTimes{intN} = sort([cellOtherSpikes{intN}(:);...
			cell2vec(cellGlobalSpikeTimesPerCellPerTrial(intN,:))]);
	end
end

