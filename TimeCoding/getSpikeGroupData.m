function [sSpikeGroup,matSpikeGroupData] = getSpikeGroupData(cellUseSpikeTimesPerCellPerTrial,intSpikeGroupSize,vecOriIdx,vecStimOnTime,vecTime,vecIFR)
	%getSpikeGroupData Retrieve data per spike group
	%   [sSpikeGroup,matSpikeGroupData] = getSpikeGroupData(cellUseSpikeTimesPerCellPerTrial,intSpikeGroupSize,vecOriIdx,vecStimOnTime,vecTime,vecIFR)
	
	%assign placefolder
	sSpikeGroup = [];
	
	%throw away spikes not in stimulus period and assign all spikes
	intTotalSpikeNum = sum(sum(cellfun(@numel,cellUseSpikeTimesPerCellPerTrial)));
	[intNumN,intTrialNum] = size(cellUseSpikeTimesPerCellPerTrial);
	matSpikeGroupData = nan(floor(intTotalSpikeNum/intSpikeGroupSize),intNumN);
	vecSpikeGroupStartIdx = nan(floor(intTotalSpikeNum/intSpikeGroupSize),1);
	vecSpikeGroupStopIdx = nan(floor(intTotalSpikeNum/intSpikeGroupSize),1);
	vecSpikeGroupDuration = nan(floor(intTotalSpikeNum/intSpikeGroupSize),1);
	vecSpikeGroupRateChange = nan(floor(intTotalSpikeNum/intSpikeGroupSize),1);
	vecSpikeGroupAvgRate = nan(floor(intTotalSpikeNum/intSpikeGroupSize),1);
	vecSpikeGroupLatency = nan(floor(intTotalSpikeNum/intSpikeGroupSize),1);
	vecSpikeGroupTrialType = nan(floor(intTotalSpikeNum/intSpikeGroupSize),1);
	vecSpikeGroupTrialNumber = nan(floor(intTotalSpikeNum/intSpikeGroupSize),1);
	vecTrialType = nan(floor(intTotalSpikeNum/intSpikeGroupSize),1);
	intSpikeGroupCounter = 0;
	intGlobalSpikeEntry = 1;
	intTotSpikeNumIFR = numel(vecTime);
	for intTrial=1:intTrialNum
		dblTrialStart = vecStimOnTime(intTrial);
		cellSpikesInTrial = cellUseSpikeTimesPerCellPerTrial(:,intTrial);
		intSpikesInTrial = sum(sum(cellfun(@numel,cellSpikesInTrial)));
		vecSpikeTimes = nan(1,intSpikesInTrial);
		vecSpikeNeuron= nan(1,intSpikesInTrial);
		intSpikeC = 1;
		for intN=1:intNumN
			vecSpikeT = cellSpikesInTrial{intN};
			intNumS = numel(vecSpikeT);
			vecAssign = intSpikeC:(intSpikeC+intNumS-1);
			intSpikeC = intSpikeC + intNumS;
			
			vecSpikeTimes(vecAssign) = vecSpikeT;
			vecSpikeNeuron(vecAssign) = intN;
		end
		
		%sort spikes
		[vecSpikeTimes,vecSort]=sort(vecSpikeTimes);
		vecSpikeNeuron = vecSpikeNeuron(vecSort);
		
		%assign to spike group
		intAssignGroups = floor(numel(vecSpikeTimes)/intSpikeGroupSize);
		for intGroup=1:intAssignGroups
			%get data
			intEndSpike = intGroup*intSpikeGroupSize;
			intStartSpike = intEndSpike-intSpikeGroupSize+1;
			vecUseSpikes = intStartSpike:intEndSpike;
			vecT = vecSpikeTimes(vecUseSpikes);
			vecN = vecSpikeNeuron(vecUseSpikes);
			vecCounts = accumarray(vecN(:),1,[intNumN 1]);
			
			%get entry in global time
			dblGlobalTimeStart = vecT(1) + dblTrialStart;
			dblGlobalTimeStop = vecT(end) + dblTrialStart;
			intStartEntry = [];
			intStopEntry = [];
			for intGlobalSpikeEntry=intGlobalSpikeEntry:intTotSpikeNumIFR
				if isempty(intStartEntry) && vecTime(intGlobalSpikeEntry) >= dblGlobalTimeStart
					intStartEntry = intGlobalSpikeEntry;
				end
				if vecTime(intGlobalSpikeEntry) >= dblGlobalTimeStop
					intStopEntry = intGlobalSpikeEntry;
					break;
				end
			end
			
			%assign
			intAssignTo = intGroup + intSpikeGroupCounter;
			matSpikeGroupData(intAssignTo,:) = vecCounts;
			vecSpikeGroupStartIdx(intAssignTo) = intStartEntry;
			vecSpikeGroupStopIdx(intAssignTo) = intStopEntry;
			vecSpikeGroupDuration(intAssignTo) = range(vecT);
			vecSpikeGroupRateChange(intAssignTo) = vecIFR(intStopEntry) - vecIFR(intStartEntry);
			vecSpikeGroupAvgRate(intAssignTo) = mean(vecIFR(intStartEntry:intStopEntry));
			vecSpikeGroupLatency(intAssignTo) = mean(vecT);
			vecSpikeGroupTrialType(intAssignTo) = vecOriIdx(intTrial);
			vecSpikeGroupTrialNumber(intAssignTo) = intTrial;
		end
		intSpikeGroupCounter = intSpikeGroupCounter + intAssignGroups;
	end
	if intSpikeGroupCounter < intTrialNum,return;end
	
	%remove trailing entries
	intSpikeGroupNum = intSpikeGroupCounter;
	matSpikeGroupData = matSpikeGroupData(1:intSpikeGroupNum,:);
	vecSpikeGroupStartIdx = vecSpikeGroupStartIdx(1:intSpikeGroupNum);
	vecSpikeGroupStopIdx = vecSpikeGroupStopIdx(1:intSpikeGroupNum);
	vecSpikeGroupDuration = vecSpikeGroupDuration(1:intSpikeGroupNum);
	vecSpikeGroupRateChange = vecSpikeGroupRateChange(1:intSpikeGroupNum);
	vecSpikeGroupAvgRate = vecSpikeGroupAvgRate(1:intSpikeGroupNum);
	vecSpikeGroupLatency = vecSpikeGroupLatency(1:intSpikeGroupNum);
	vecSpikeGroupTrialType = vecSpikeGroupTrialType(1:intSpikeGroupNum);
	vecSpikeGroupTrialNumber = vecSpikeGroupTrialNumber(1:intSpikeGroupNum);
	
	%% decode all spike groups
	dblLambda = 1;
	intTypeCV = 0;
	
	%do logistic regression
	intOriNum = numel(unique(vecOriIdx));
	vecPriorDistribution = accumarray(vecSpikeGroupTrialType(:),1,[intOriNum 1]);
	intVerbose = 1;
	[dblPerformanceCV,vecSpikeGroupDecodedTrial,matPosteriorProbability,dblMeanErrorDegs,matConfusion,matWeights,matAggActivation,matAggWeights,vecRepetition,vecTrialTypeIdx] = ...
		doCrossValidatedDecodingLR(matSpikeGroupData,vecSpikeGroupTrialType,intTypeCV,[],dblLambda,intVerbose);
	
	vecSpikeGroupCorrect = vecSpikeGroupDecodedTrial(:) == vecTrialTypeIdx;
	vecSpikeGroupConfidence = nan(intSpikeGroupNum,1);
	for intStimType=1:size(matPosteriorProbability,1)
		vecSpikeGroupConfidence(vecTrialTypeIdx==intStimType) = matPosteriorProbability(intStimType,vecTrialTypeIdx==intStimType);
	end
	
	%put in struct
	sSpikeGroup = struct;
	for intSpikeGroup=intSpikeGroupNum:-1:1
		sSpikeGroup(intSpikeGroup).Duration = vecSpikeGroupDuration(intSpikeGroup);
		sSpikeGroup(intSpikeGroup).Latency = vecSpikeGroupLatency(intSpikeGroup);
		sSpikeGroup(intSpikeGroup).TrialType = vecTrialTypeIdx(intSpikeGroup);
		sSpikeGroup(intSpikeGroup).TrialNumber = vecSpikeGroupTrialNumber(intSpikeGroup);
		sSpikeGroup(intSpikeGroup).DecodedType = vecSpikeGroupDecodedTrial(intSpikeGroup);
		sSpikeGroup(intSpikeGroup).Correct = vecSpikeGroupCorrect(intSpikeGroup);
		sSpikeGroup(intSpikeGroup).Confidence = vecSpikeGroupConfidence(intSpikeGroup);
		sSpikeGroup(intSpikeGroup).StartIdx = vecSpikeGroupStartIdx(intSpikeGroup);
		sSpikeGroup(intSpikeGroup).StopIdx = vecSpikeGroupStopIdx(intSpikeGroup);
		sSpikeGroup(intSpikeGroup).RateChange = vecSpikeGroupRateChange(intSpikeGroup);
		sSpikeGroup(intSpikeGroup).AvgRate = vecSpikeGroupAvgRate(intSpikeGroup);
	end
	
end

