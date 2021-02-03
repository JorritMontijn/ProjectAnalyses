function [dblP] = getISItest(vecSpikeTimes,matEventTimes,dblUseMaxDur,intResampleNum,boolUseGumbel)
	
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
	
	%get boolPlot
	if ~exist('boolUseGumbel','var') || isempty(boolUseGumbel)
		boolUseGumbel = true;
	end
	
	%% build onset/offset vectors
	vecEventStarts = matEventTimes(:,1);
	
	%% get PSTH
	vecUseSpikeTimes = vecSpikeTimes(vecSpikeTimes < (max(vecEventStarts)+3*dblUseMaxDur) & (vecSpikeTimes > (min(vecEventStarts)-3*dblUseMaxDur)));
	dblBinSize = 1/1000;
	vecBinEdges = 0:dblBinSize:dblUseMaxDur;
	intBins = numel(vecBinEdges)-1;
	[vecTrialPerSpike,vecTimePerSpike] = getSpikesInTrial(vecUseSpikeTimes,vecEventStarts,dblUseMaxDur);
	intTrials = numel(vecEventStarts);
	matPSTH = nan(intBins,intTrials);
	for intTrial=1:intTrials
		matPSTH(:,intTrial) = histcounts(vecTimePerSpike(vecTrialPerSpike==intTrial),vecBinEdges);
	end
	
	%% get ISI distro
	vecISI = diff(vecUseSpikeTimes);
	matRandPSTH = nan(intBins,intTrials,intResampleNum);
	dblT0 = vecUseSpikeTimes(1);
	parfor intIter = 1:intResampleNum
		%generate spikes
		vecRandISI = vecISI(randperm(numel(vecISI)));
		vecRandSpikeTimes = dblT0 + cumsum(vecRandISI) - vecRandISI(1);
		
		%build PSTH
		[vecTrialPerSpike,vecTimePerSpike] = getSpikesInTrial(vecRandSpikeTimes,vecEventStarts,dblUseMaxDur);
		for intTrial=1:intTrials
			matRandPSTH(:,intTrial,intIter) = histcounts(vecTimePerSpike(vecTrialPerSpike==intTrial),vecBinEdges);
		end
	end
	
	%% calculate significance
	matRandPSTH(isnan(matRandPSTH))=0;
	vecMeanPSTH = mean(matPSTH,2);
	matMeanRandPSTH = squeeze(mean(matRandPSTH,2));
	if ~boolUseGumbel
		vecZ = (vecMeanPSTH-mean(matMeanRandPSTH(:)))./std(matMeanRandPSTH(:));
		dblZ = max(abs(vecZ));
		dblP = normcdf(-dblZ)*2;
	else
		%real peak
		vecRealDiff = vecMeanPSTH - mean(vecMeanPSTH);
		
		%find highest peak and retrieve value
		matRandDiff = matMeanRandPSTH - mean(matMeanRandPSTH,1);
		vecMaxRandD = max(abs(matRandDiff),[],1);
		dblRandMu = mean(vecMaxRandD);
		dblRandVar = var(vecMaxRandD);
		dblPosD= max(abs(vecRealDiff));
		
		%calculate statistical significance using Gumbel distribution
		dblP = getGumbel(dblRandMu,dblRandVar,dblPosD);
	end
end

