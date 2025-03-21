function [dblP_KS,dblP_G,dblIntP_G] = getISItest(vecSpikeTimes,matEventTimes,dblUseMaxDur,intResampleNum)
	%[dblP_KS,dblP_G,dblIntP_G] = getISItest(vecSpikeTimes,matEventTimes,dblUseMaxDur,intResampleNum)
	
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
		intResampleNum = 100;
	end
	
	%% build onset/offset vectors
	vecEventStarts = matEventTimes(:,1);
	
	%% pre-allocate
	dblP_KS = 1;
	dblP_G = 1;
	dblIntP_G = 1;
	
	%% get PSTH
	dblStopHorizon = (max(vecEventStarts)+3*dblUseMaxDur);
	dblStartHorizon = (min(vecEventStarts)-3*dblUseMaxDur);
	vecUseSpikeTimes = vecSpikeTimes(vecSpikeTimes < dblStopHorizon & (vecSpikeTimes > dblStartHorizon));
	if numel(vecUseSpikeTimes)<3,return;end
	dblLambda = numel(vecUseSpikeTimes)/(dblStopHorizon - dblStartHorizon);
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
	matRandPSTH = zeros(intBins,intTrials,intResampleNum);
	dblT0 = vecUseSpikeTimes(1);
	for intIter = 1:intResampleNum
		%generate spikes
		vecRandISI = vecISI(randperm(numel(vecISI)));
		vecRandSpikeTimes = dblT0 + cumsum(vecRandISI) - vecRandISI(1);
		
		%build PSTH
		[vecTrialPerSpike,vecTimePerSpike] = getSpikesInTrial(vecRandSpikeTimes,vecEventStarts,dblUseMaxDur);
		for intTrial=1:intTrials
			matRandPSTH(:,intTrial,intIter) = histcounts(vecTimePerSpike(vecTrialPerSpike==intTrial),vecBinEdges);
		end
	end
	
	
	%% k-s test
	matRandPSTH(isnan(matRandPSTH))=0;
	matRandMeansPSTH = squeeze(mean(matRandPSTH,2));
	vecMeanPSTH = mean(matPSTH,2);
	[h,dblP_KS]=kstest2(vecMeanPSTH,matRandMeansPSTH(:));
	
	%% Gumbel
	if nargout > 1
		%real peak
		vecRealDiff = vecMeanPSTH - mean(vecMeanPSTH);
		
		%find highest peak and retrieve value
		matRandDiff = matRandMeansPSTH - mean(matRandMeansPSTH,1);
		vecMaxRandD = max(abs(matRandDiff),[],1);
		dblRandMu = mean(vecMaxRandD);
		dblRandVar = var(vecMaxRandD);
		dblPosD= max(abs(vecRealDiff));
		
		%calculate statistical significance using Gumbel distribution
		dblP_G = getGumbel(dblRandMu,dblRandVar,dblPosD);
	end
	
	%% cumsum
	if nargout > 1
		vecMeanD = vecMeanPSTH - mean(vecMeanPSTH);
		vecCumSum = cumsum(vecMeanD);
		vecD = vecCumSum - mean(vecCumSum);
		
		matMeanRandD = bsxfun(@minus,matRandMeansPSTH,mean(matRandMeansPSTH,1));
		matCumSum = cumsum(matMeanRandD,1);
		matRandD = bsxfun(@minus,matCumSum,mean(matCumSum,1));
		
		%clf;
		%plot(vecD)
		%hold on
		%plot(matRandD)
		%hold off
		
		%gumbel
		vecIntMaxRandD = max(abs(matRandD),[],1);
		dblIntRandMu = mean(vecIntMaxRandD);
		dblIntRandVar = var(vecIntMaxRandD);
		dblIntPosD= max(abs(vecD));
		
		%calculate statistical significance using Gumbel distribution
		dblIntP_G = getGumbel(dblIntRandMu,dblIntRandVar,dblIntPosD);
	end
end

