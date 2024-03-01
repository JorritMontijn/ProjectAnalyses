function [sPopEvents,sMergedPopEvents,vecStimTime,vecStimIFR] = getPopEvents(vecIFR,vecTime,vecOrigStimOnTime,sPeakOpts,sAllSpike,dblCutOff)
	%getPopEvents Summary of this function goes here
	%   [sPopEvents,sMergedPopEvents,vecStimTime,vecStimIFR] = getPopEvents(vecIFR,vecTime,vecOrigStimOnTime,sPeakOpts,sAllSpike,dblCutOff)
	
	%get params
	intLag = sPeakOpts.intLag;
	dblThreshZ = sPeakOpts.dblThreshZ;
	dblInfluence = sPeakOpts.dblInfluence;
	vecAllSpikeTime = sAllSpike.vecAllSpikeTime;
	vecAllSpikeNeuron = sAllSpike.vecAllSpikeNeuron;
	dblStartEpoch = sAllSpike.dblStartEpoch;
	dblEpochDur = sAllSpike.dblEpochDur;

	%mean-subtraction
	[signals,avgFilter,stdFilter] = detectpeaks(vecIFR,intLag,dblThreshZ,dblInfluence);
	%vecIFR = (vecIFR - avgFilter)./stdFilter;
	vecNormIFR = (vecIFR - avgFilter)./avgFilter;
	indRem = isinf(vecNormIFR);
	vecNormIFR(indRem) = [];
	vecNormIFR = vecNormIFR + (1e-10)*rand(size(vecNormIFR)); %break identical values
	vecNormTime = vecTime(~indRem);
	
	% peaks
	dblStimDur = min(diff(vecOrigStimOnTime));
	indStimSpikes = vecNormTime>(vecOrigStimOnTime(1)-dblStimDur) & vecNormTime<(vecOrigStimOnTime(end)+dblStimDur);
	vecStimIFR_Raw = vecIFR(indStimSpikes);
	vecStimIFR = vecNormIFR(indStimSpikes);
	vecStimTime = vecNormTime(indStimSpikes);
	[vecPeakHeight,vecPeakLocs,w,p] = findpeaks(vecStimIFR);
	
	%remove peaks with identical subsequent values
	%vecPeakLocs(vecStimIFR(vecPeakLocs)==vecStimIFR(vecPeakLocs+1))=[];
	
	%threshold peaks
	vecStimPeakLocs = vecPeakLocs(vecPeakHeight>dblCutOff);
	
	%retain only stim epoch
	indStimEpoch = vecAllSpikeTime > dblStartEpoch & vecAllSpikeTime < (dblStartEpoch + dblEpochDur);
	vecPopSpikeTime = vecAllSpikeTime(indStimEpoch);
	vecPopNeuronId = vecAllSpikeNeuron(indStimEpoch);
	
	%get raw pop events
	vecPopEventTimes = vecStimTime(vecStimPeakLocs);
	intPopEventNum = numel(vecPopEventTimes);
	vecPopEventLocs = nan(1,intPopEventNum);
	vecPopEventLocsIFR = nan(1,intPopEventNum);
	parfor intEvent=1:intPopEventNum
		[dummy,vecPopEventLocs(intEvent)] = min(abs(vecPopEventTimes(intEvent)-vecPopSpikeTime));
		[dummy,vecPopEventLocsIFR(intEvent)] = min(abs(vecPopEventTimes(intEvent)-vecStimTime));
	end
	
	%merged peaks
	matPeakDomain = mergepeaks(vecStimTime,vecStimIFR,vecStimPeakLocs);
	vecMergedPeakLocs = matPeakDomain(:,1);
	
	%get merged pop events
	vecMergedPopEventTimes = vecStimTime(vecMergedPeakLocs);
	intMergedPopEventNum = numel(vecMergedPopEventTimes);
	vecMergedPopEventLocs = nan(1,intMergedPopEventNum);
	vecMergedPopEventLocsIFR = nan(1,intMergedPopEventNum);
	parfor intEvent=1:intMergedPopEventNum
		[dummy,vecMergedPopEventLocs(intEvent)] = min(abs(vecMergedPopEventTimes(intEvent)-vecPopSpikeTime));
		[dummy,vecMergedPopEventLocsIFR(intEvent)] = min(abs(vecMergedPopEventTimes(intEvent)-vecStimTime));
	end
	
	%assign output to structs
	sPopEvents = struct;
	sPopEvents.Time = vecPopEventTimes;
	sPopEvents.Loc = vecPopEventLocs;
	sPopEvents.LocIFR = vecPopEventLocsIFR;
	sMergedPopEvents = struct;
	sMergedPopEvents.Time = vecMergedPopEventTimes;
	sMergedPopEvents.Loc = vecMergedPopEventLocs;
	sMergedPopEvents.LocIFR = vecMergedPopEventLocsIFR;
end

