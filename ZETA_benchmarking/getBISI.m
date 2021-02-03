function [dblP] = getBISI(vecSpikeTimes,matEventTimes,dblUseMaxDur,intResampleNum)
	
	%% check inputs and pre-allocate error output
	dblP = 1;
	if numel(vecSpikeTimes) < 3
		return;
	end
	
	%% prep data
	%ensure orientation
	vecSpikeTimes = vecSpikeTimes(:);
	
	%calculate stim/base difference?
	if size(matEventTimes,2) > 2
		matEventTimes = matEventTimes';
	end
	
	%trial dur
	if ~exist('dblUseMaxDur','var') || isempty(dblUseMaxDur)
		dblUseMaxDur = min(diff(matEventTimes(:,1)));
	end
	
	%get boolPlot
	if ~exist('intResampleNum','var') || isempty(intResampleNum)
		intResampleNum = 0;
	end
	
	%% build onset/offset vectors
	vecEventStarts = matEventTimes(:,1);
	
	%% prepare interpolation points
	vecStartOnly = vecEventStarts(:,1);
	vecSpikeT = getSpikeT(vecSpikeTimes,vecStartOnly,dblUseMaxDur);
	
	%% run normal
	%get data
	[vecRealDiff,vecRealFrac,vecRealFracLinear] = ...
		getTempOffset(vecSpikeT,vecSpikeTimes,vecStartOnly,dblUseMaxDur);
	
	%% get ISI distro
	vecUseSpikeTimes = vecSpikeTimes(vecSpikeTimes < (max(vecStartOnly)+3*dblUseMaxDur) & (vecSpikeTimes > (min(vecStartOnly)-3*dblUseMaxDur)));
	if numel(vecRealDiff) < 3 || numel(vecUseSpikeTimes) < 4
		return
	end
	vecISI = diff(vecUseSpikeTimes);
	intSpikes = numel(vecRealDiff);
	matRandDiff = nan(intSpikes,intResampleNum);
	dblT0 = vecUseSpikeTimes(1);
	
	%% run bootstraps
	parfor intResampling=1:intResampleNum
		%% get random subsample
		vecRandISI = vecISI(randperm(numel(vecISI)));
		vecRandSpikeTimes = dblT0 + cumsum(vecRandISI) - vecRandISI(1);
		
		%get temp offset
		vecRandDiff = getTempOffset(vecSpikeT,vecRandSpikeTimes,vecStartOnly,dblUseMaxDur);
		
		%assign data
		matRandDiff(:,intResampling) = vecRandDiff - mean(vecRandDiff);
	end
	
	%% calculate measure of effect size (for equal n, d' equals Cohen's d)
	dblP = 1;
	if numel(vecRealDiff) < 3
		return
	end
	
	%find highest peak and retrieve value
	vecMaxRandD = max(abs(matRandDiff),[],1);
	dblRandMu = mean(vecMaxRandD);
	dblRandVar = var(vecMaxRandD);
	dblPosD= max(abs(vecRealDiff));
	
	%calculate statistical significance using Gumbel distribution
	dblP = getGumbel(dblRandMu,dblRandVar,dblPosD);
end

