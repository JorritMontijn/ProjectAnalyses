function [dblP] = getPoissonTest(vecSpikeTimes,matEventTimes,dblUseMaxDur,intResampleNum)
	
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
	
	%% get PSTH
	dblStopHorizon = (max(vecEventStarts)+3*dblUseMaxDur);
	dblStartHorizon = (min(vecEventStarts)-3*dblUseMaxDur);
	vecUseSpikeTimes = vecSpikeTimes(vecSpikeTimes < dblStopHorizon & (vecSpikeTimes > dblStartHorizon));
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
	
	%% k-s test
	matRandH0 = poissrnd(dblLambda*dblBinSize,[intBins intTrials intResampleNum]);
	vecRandMeansH0 = squeeze(mean(matRandH0,2));	
	vecMeans = mean(matPSTH,2);
	[h,dblP]=kstest2(vecMeans,vecRandMeansH0(:));
end

