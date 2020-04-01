function [dblPerformance,vecDecodedIndexCV,matPosteriorProbability,matActivation] = decodeLR(matWeights,matData,vecTrialTypes)
	%decodeLR Uses LR weight matrix to decode data
	%   [dblPerformance,vecDecodedIndexCV,matPosteriorProbability,matActivation] = decodeLR(matWeights,matData,vecTrialTypes)
	
	%add bias if necessary
	if size(matData,1) == (size(matWeights,1) - 1)
		matDataPlusLin = [matData; ones(1,size(matData,2))];
	else
		matDataPlusLin = matData;
	end
	
	%decode
	vecTrialTypeIdx = label2idx(vecTrialTypes(:));
	matActivation = matWeights'*matDataPlusLin;
	matPosteriorProbability = exp(bsxfun(@minus,matActivation,logsumexp(matActivation,1)))'; %softmax
	[dummy, vecDecodedIndexCV] = max(matPosteriorProbability,[],2);
	dblPerformance=sum(vecDecodedIndexCV==vecTrialTypeIdx)/numel(vecTrialTypeIdx);
end

