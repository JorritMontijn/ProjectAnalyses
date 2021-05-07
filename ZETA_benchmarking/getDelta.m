function [vecThisDiff] = getDelta(vecSpikeTimes,vecEventStarts,dblUseMaxDur)
	%getDelta Summary of this function goes here
	%   [vecThisDiff] = getDelta(vecSpikeTimes,vecEventStarts,dblUseMaxDur)
	
	vecSpikeT = getSpikeT(vecSpikeTimes,vecEventStarts,dblUseMaxDur);
	%pre-allocate
	vecThisSpikeTimes = unique(getSpikeT(vecSpikeTimes,vecEventStarts,dblUseMaxDur));
	vecThisSpikeFracs = linspace(1/numel(vecThisSpikeTimes),1,numel(vecThisSpikeTimes))';
	vecThisFrac = interp1(vecThisSpikeTimes,vecThisSpikeFracs,vecSpikeT);
	
	%get linear fractions
	vecThisFracLinear = (vecSpikeT./dblUseMaxDur);
	vecThisDiff = vecThisFrac - vecThisFracLinear;
	
end

