function [matMeanRate,indTuned,cellSpikeTimes,sOut,cellSpikeTimesPerCellPerTrial,indResp] = ...
		SimPrepData(cellSpikeTimesIn,vecStimOnTime,vecStimOffTime,vecStimTypes)
	%NpxPrepData Summary of this function goes here
	%   [matMeanRate,indTuned,cellSpikeTimes,sOut,cellSpikeTimesPerCellPerTrial,vecStimOnStitched,vecNonStat,dblBC,dblMaxDevFrac] = ...
	%		NpxPrepData(cellSpikeTimesIn,vecStimOnTime,vecStimOffTime,vecStimTypes)
	
	%% set cut-offs
	dblMinRate = 0.1;
	
	%get dur
	cellSpikeTimes = cellSpikeTimesIn;
	dblDur = median(vecStimOffTime-vecStimOnTime);
	matRawData = getSpikeCounts(cellSpikeTimes,vecStimOnTime,dblDur)./dblDur;
	indResp = sum(matRawData,2)'>(size(matRawData,2)/dblDur)*dblMinRate;
	matMeanRate = matRawData(indResp,:);
	
	%% remove non-responsive cells
	[vecStimIdx,vecUnique,vecPriorDistribution,cellSelect,vecRepetition] = val2idx(vecStimTypes);
	sOut = getTuningCurves(matMeanRate,vecStimTypes,0);
	cellSpikeTimes(~indResp)=[];
	indTuned = sOut.vecOriAnova<0.05;
	intTrialNum = numel(vecStimTypes);
	intRespN = size(matMeanRate,1);
	
	%% assign spike times per trial
	%get spikes per trial per neuron
	cellSpikeTimesPerCellPerTrial = cell(intRespN,intTrialNum);
	hTic=tic;
	for intN=1:intRespN
		if toc(hTic) > 5
			fprintf('   Prepping %d/%d [%s]\n',intN,intRespN,getTime);
			hTic=tic;
		end
		%real
		cellSpikeTimes{intN} = cellSpikeTimes{intN}(:);
		[vecTrialPerSpike,vecTimePerSpike] = getSpikesInTrial(cellSpikeTimes{intN},vecStimOnTime);
		%remove spikes during offset
		indRemSpikes = vecTimePerSpike>dblDur | vecTrialPerSpike == 0;
		vecTimePerSpike(indRemSpikes)=[];
		vecTrialPerSpike(indRemSpikes)=[];
		intSpikeNum = numel(vecTrialPerSpike);
		intNextStart = 1;
		for intTrial=1:intTrialNum
			intThisStart = intNextStart;
			for intNextStart=intNextStart:intSpikeNum
				if vecTrialPerSpike(intNextStart) ~= intTrial
					break
				end
			end
			vecUseSpikes = intThisStart:(intNextStart-1);
			cellSpikeTimesPerCellPerTrial{intN,intTrial} = vecTimePerSpike(vecUseSpikes);
		end
	end
end
