function dblCorrection = getPupilOffset(cellSpikeTimes,x0,sPupil)
	%UNTITLED Summary of this function goes here
	%   Detailed explanation goes here
	
	%set globals
	global vecTime;
	global vecVals;
	global vecT_E;
	global vecSpikeNum;
	
	%pupil data
	vecTime = sPupil.vecTime;
	vecVals = sPupil.vecRadius;
	
	%get mean activity
	vecAllSpikes = cell2vec(cellSpikeTimes);
	dblMin = min(min(vecAllSpikes),min(vecTime))+10;
	dblMax = max(max(vecAllSpikes),max(vecTime))-10;
	vecT_E = dblMin:dblMax;
	vecT_C = vecT_E(2:end)-median(diff(vecT_E));
	vecSpikeNum = histcounts(vecAllSpikes,vecT_E);
	
	%get minimum error
	dblCorrection = fminsearch(@getPupilAlignmentError,x0);
end