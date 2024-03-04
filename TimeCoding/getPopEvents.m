function [sPopEvents,sMergedPopEvents,vecStimTime_Filt,vecStimIFR_Filt,vecStimIFR_Raw] = getPopEvents(vecIFR,vecTime,vecOrigStimOnTime,sPeakOpts,sAllSpike,dblCutOff)
	%getPopEvents Summary of this function goes here
	%   [sPopEvents,sMergedPopEvents,vecStimTime,vecStimIFR,vecStimIFR_Raw] = getPopEvents(vecIFR,vecTime,vecOrigStimOnTime,sPeakOpts,sAllSpike,dblCutOff)
	
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
	vecFiltIFR = (vecIFR - avgFilter)./avgFilter;
	indRem = isinf(vecFiltIFR);
	vecIFR(indRem) = [];
	vecFiltIFR(indRem) = [];
	vecFiltIFR = vecFiltIFR + (1e-10)*rand(size(vecFiltIFR)); %break identical values
	vecFiltTime = vecTime(~indRem);
	
	% peaks
	dblStimDur = min(diff(vecOrigStimOnTime));
	indStimSpikes = vecFiltTime>(vecOrigStimOnTime(1)-dblStimDur) & vecFiltTime<(vecOrigStimOnTime(end)+dblStimDur);
	vecStimIFR_Raw = vecIFR(indStimSpikes);
	vecStimIFR_Filt = vecFiltIFR(indStimSpikes);
	vecStimTime_Filt = vecFiltTime(indStimSpikes);
	[vecPeakHeight,vecPeakLocs,w,p] = findpeaks(vecStimIFR_Filt);
	
	%remove peaks with identical subsequent values
	%vecPeakLocs(vecStimIFR(vecPeakLocs)==vecStimIFR(vecPeakLocs+1))=[];
	
	%threshold peaks
	vecCulledPeakLocs = vecPeakLocs(vecPeakHeight>dblCutOff);
	
	%retain only stim epoch
	indStimEpoch = vecAllSpikeTime > dblStartEpoch & vecAllSpikeTime < (dblStartEpoch + dblEpochDur);
	vecPopSpikeTime = vecAllSpikeTime(indStimEpoch);
	vecPopNeuronId = vecAllSpikeNeuron(indStimEpoch);
	
	%get raw pop events
	vecPopEventTimes = vecStimTime_Filt(vecCulledPeakLocs);
	intPopEventNum = numel(vecPopEventTimes);
	vecPopEventLocs = nan(1,intPopEventNum);
	parfor intEvent=1:intPopEventNum
		[dummy,vecPopEventLocs(intEvent)] = min(abs(vecPopEventTimes(intEvent)-vecStimTime_Filt));
	end
	
	%merged peaks
	matPeakDomain = mergepeaks(vecStimTime_Filt,vecStimIFR_Filt,vecCulledPeakLocs);
	vecMergedPeakLocs = matPeakDomain(:,1);
	
	%get merged pop events
	vecMergedPopEventTimes = vecStimTime_Filt(vecMergedPeakLocs);
	intMergedPopEventNum = numel(vecMergedPopEventTimes);
	vecMergedPopEventLocs = nan(1,intMergedPopEventNum);
	parfor intEvent=1:intMergedPopEventNum
		[dummy,vecMergedPopEventLocs(intEvent)] = min(abs(vecMergedPopEventTimes(intEvent)-vecStimTime_Filt));
	end
	
	%assign output to structs
	sPopEvents = struct;
	sPopEvents.Time = vecPopEventTimes;
	sPopEvents.Loc = vecPopEventLocs;
	sMergedPopEvents = struct;
	sMergedPopEvents.Time = vecMergedPopEventTimes;
	sMergedPopEvents.Loc = vecMergedPopEventLocs;
end

