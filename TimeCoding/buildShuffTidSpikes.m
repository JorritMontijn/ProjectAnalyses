function [cellUseSpikeTimesPerCellPerTrial,cellSpikeTimes] = buildShuffTidSpikes(cellSpikeTimesReal,vecStimOnTime,vecStimIdx,dblTrialDur);
	%buildShuffTidSpikes Summary of this function goes here
	%   [cellUseSpikeTimesPerCellPerTrial,cellSpikeTimes] = buildShuffTidSpikes(cellSpikeTimesReal,vecStimOnTime,vecStimIdx,dblTrialDur);
	
	%get beginning & end vectors
	intNumN = numel(cellSpikeTimesReal);
	intTrialNum = numel(vecStimOnTime);
	cellOtherSpikes = cell(intNumN,1);
	cellSpikeTimesPerCellPerTrial = cell(intNumN,intTrialNum);
	hTic=tic;
	for intN=1:intNumN
		%msg
		if toc(hTic) > 5
			fprintf('  ShuffTid %d/%d [%s]\n',intN,intNumN,getTime);
			hTic=tic;
		end
		
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
	cellGlobalSpikeTimesPerCellPerTrial = cellUseSpikeTimesPerCellPerTrial;
	vecOrigTrials = find(vecStimIdx==intStimType);
	for intTrial=1:intTrialNum
		dblOnset = vecStimOnTime(intTrial);
		for intN=1:intNumN
			cellGlobalSpikeTimesPerCellPerTrial{intN,intTrial} = cellGlobalSpikeTimesPerCellPerTrial{intN,intTrial} + dblOnset;
		end
	end
	
	%compile
	for intN=1:intNumN
		cellSpikeTimes{intN} = getUniqueSpikes(sort([cellOtherSpikes{intN}(:);...
			cell2vec(cellGlobalSpikeTimesPerCellPerTrial(intN,:))]));
	end
end

