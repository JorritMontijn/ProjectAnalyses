function [matResp,indTuned,indResp,cellSpikeTimes,sOut] = NpxPrepData(cellSpikeTimes,vecStimOnTime,vecStimOffTime,vecOrientation)
	%NpxPrepData Summary of this function goes here
	%   [matResp,indTuned,indResp,cellSpikeTimes,sOut] = NpxPrepData(cellSpikeTimes,vecStimOnTime,vecStimOffTime,vecOrientation)
	
	%set cut-off
	dblMinRate = 0.1;
	
	%get dur
	dblDur = median(vecStimOffTime-vecStimOnTime);
	matData = getSpikeCounts(cellSpikeTimes,vecStimOnTime,dblDur)./dblDur;
	indResp = sum(matData,2)'>(size(matData,2)/dblDur)*dblMinRate;
	matResp = matData(indResp,:);
	
	%remove untuned cells
	vecOri180 = mod(vecOrientation,180)*2;
	sOut = getTuningCurves(matResp,vecOri180,0);
	cellSpikeTimes(~indResp)=[];
	indTuned = sOut.vecOriAnova<0.05;
end

