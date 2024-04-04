function [cellUseSpikeTimesPerCellPerTrial,cellSpikeTimes] = buildUniformSpikes(cellSpikeTimesReal,vecStimOnTime,vecStimIdx,dblTrialDur)
	%buildShuffTidSpikes Summary of this function goes here
	%   [cellUseSpikeTimesPerCellPerTrial,cellSpikeTimes] = buildShuffTidSpikes(cellSpikeTimesReal,vecStimOnTime,vecStimIdx,dblTrialDur)
	
	
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
			%make spike times uniform in trial
			indReplaceSpikes = vecTrialPerSpike==intTrial;
			vecSpikeT = vecTimePerSpike(indReplaceSpikes);
			%save pre- and post-spikes
			dblStimDur = vecStimOffTime(intTrial) - vecStimOnTime(intTrial);
			indPrepost = (vecSpikeT > (dblTrialDur));
			
			%randomize during-spikes
			vecSpikeT_prepost = vecSpikeT(indPrepost);
			vecSpikeT_unirand = rand(sum(~indPrepost),1)*dblTrialDur;
			vecSpikeT = sort(cat(1,vecSpikeT_unirand,vecSpikeT_prepost));
			error is onset removed?
			
			%assign
			cellSpikeTimesReal{intN}(indReplaceSpikes) = vecSpikeT;
			cellSpikeTimesPerCellPerTrial{intN,intTrial} = vecSpikeT;
		end
		cellSpikeTimesReal{intN} = sort(cellSpikeTimesReal{intN});
	end
end

