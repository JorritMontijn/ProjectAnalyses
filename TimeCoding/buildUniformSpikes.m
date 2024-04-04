function [cellSpikeTimesPerCellPerTrial,cellSpikeTimes] = buildUniformSpikes(cellSpikeTimes,vecStimOnTime,vecStimIdx,dblTrialDur)
	%buildShuffTidSpikes Create randomized data
	%   [cellUseSpikeTimesPerCellPerTrial,cellSpikeTimes] = buildShuffTidSpikes(cellSpikeTimesReal,vecStimOnTime,vecStimIdx,dblTrialDur)
	%Build data such that for each neuron and trial spike times are uniform, but spike  number is
	%matched on neuron and trial basis
	
	%get beginning & end vectors
	intNumN = numel(cellSpikeTimes);
	intTrialNum = numel(vecStimOnTime);
	cellSpikeTimesPerCellPerTrial = cell(intNumN,intTrialNum);
	hTic=tic;
	for intN=1:intNumN
		%msg
		if toc(hTic) > 5
			fprintf('  Creating uniform random spike times per trial for neuron %d/%d [%s]\n',intN,intNumN,getTime);
			hTic=tic;
		end
		
		%real
		[vecTrialPerSpike,vecTimePerSpike] = getSpikesInTrial(cellSpikeTimes{intN},vecStimOnTime);
		for intTrial=1:intTrialNum
			%make spike times uniform in trial
			indReplaceSpikes = vecTrialPerSpike==intTrial;
			vecSpikeT = vecTimePerSpike(indReplaceSpikes);
			%save pre- and post-spikes
			indPrepost = (vecSpikeT > (dblTrialDur));
			
			%randomize during-spikes
			vecSpikeT_prepost = vecSpikeT(indPrepost);
			vecSpikeT_unirand = rand(sum(~indPrepost),1)*dblTrialDur;
			vecSpikeT = sort(cat(1,vecSpikeT_unirand,vecSpikeT_prepost));
			
			%assign
			cellSpikeTimes{intN}(indReplaceSpikes) = vecSpikeT+vecStimOnTime(intTrial);
			cellSpikeTimesPerCellPerTrial{intN,intTrial} = vecSpikeT;
		end
		cellSpikeTimes{intN} = sort(cellSpikeTimes{intN});
	end
end

