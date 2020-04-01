function [dblPercVarInDirOverShuffled,vecCI,dblCohensD] = getVarInDir(matPoints,vecRef,intIters)
	%getVarInDir Returns the variance in the supplied data set (matPoints)
	%along the specified dimension (vecRef), relative to a shuffled,
	%randomized version of that data set. Syntax:
	%   [dblPercVarInDirOverShuffled,vecCI,dblCohensD] = getVarInDir(matPoints,vecRef,intIters)
	
	%get data
	intD=size(matPoints,2);
	intPoints = size(matPoints,1);
	if intPoints < intD
		error([mfilename ':WrongDims'],'Number of dimensions is larger than number of points; please make sure matrix is in form [Trials x Neurons]');
	end
	
	%check inputs
	if nargin < 3 || ~exist('intIters','var') || isempty(intIters)
		intIters = 1000;
	end
	
	%recenter
	vecMu1 = xmean(matPoints,1);
	matRecentered1 = bsxfun(@minus,matPoints,vecMu1);
	vecFprime = vecRef-vecMu1;
	
	%calc cosine similarities original
	vecNumerator = sum(bsxfun(@mtimes,matRecentered1,vecFprime),2);
	vecDenominator = sqrt(sum(matRecentered1.^2,2)) * norm(vecFprime);
	vecCosSim = abs(vecNumerator ./ vecDenominator);
	
	%shuffle
	matShuffled1 = nan(size(matRecentered1));
	vecPercDiffs = nan(1,intIters);
	for intIter=1:intIters
		for intDim=1:size(matRecentered1,2)
			matShuffled1(:,intDim) = matRecentered1(randperm(intPoints),intDim);
		end
		
		%calc cosine similarities shuffled
		vecNumerator = sum(bsxfun(@mtimes,matShuffled1,vecFprime),2);
		vecDenominator = sqrt(sum(matShuffled1.^2,2)) * norm(vecFprime);
		vecCosSimShuffled = abs(vecNumerator ./ vecDenominator);
		
		vecPercDiffs(intIter) = (mean(vecCosSim) ./ mean(vecCosSimShuffled)) * 100;
		
	end
	
	%get statistics
	dblPercVarInDirOverShuffled = mean(vecPercDiffs);
	dblCohensD = (mean(vecPercDiffs)-100)/std(vecPercDiffs);
	vecCI = getCI(vecPercDiffs);
end
