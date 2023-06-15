function dblErr = getPupilAlignmentError(dblCorrection)
	%UNTITLED2 Summary of this function goes here
	%   Detailed explanation goes here
	
	%get globals
	global vecTime;
	global vecVals;
	global vecT_E;
	global vecSpikeNum;
	
	
	%bin
	[vecCounts,vecPupilSize] = makeBins(vecTime+dblCorrection,vecVals,vecT_E);
	vecPupilSize = vecPupilSize(:);
	vecSpikeNum = vecSpikeNum(:);
	indUseTrials = ~isnan(vecPupilSize) & ~isnan(vecSpikeNum);
	dblErr = 1-corr(vecPupilSize(:),vecSpikeNum(:));
	
end

