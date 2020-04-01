function [vecProjectedPoints1D,dblRefNorm,matProjectedPointsShuffled,vecProjectedPointsClass2] = getProjOnLine(matPoints,varRef,intIters)
	%UNTITLED Summary of this function goes here
	%   Detailed explanation goes here
	
	%get data
	intD=size(matPoints,2);
	intPoints = size(matPoints,1);
	if intPoints < intD
		error([mfilename ':WrongDims'],'Number of dimensions is larger than number of points; please make sure matrix is in form [Trials x Neurons]');
	end
	
	%check inputs
	if nargin < 3 || ~exist('intIters','var') || isempty(intIters)
		intIters = 1;
	end
	if isvector(varRef)
		vecRef = varRef;
		matPoints2 = [];
	else
		matPoints2 = varRef;
		vecRef = xmean(matPoints2,1);
	end
	
	%recenter
	vecMu1 = xmean(matPoints,1);
	matRecentered1 = bsxfun(@minus,matPoints,vecMu1);
	vecFprime = vecRef-vecMu1;
	
	%calculate projected points
	vecProjectedPoints1D = nan(size(matRecentered1,1),1);
	for intTrial=1:size(matRecentered1,1)
		vecPoint = matRecentered1(intTrial,:);
		vecOrth = ((vecPoint*vecFprime')/(vecFprime*vecFprime')).*vecFprime - vecFprime*100;
		vecProjectedPoints1D(intTrial) = norm(vecOrth);
	end
	vecProjectedPoints1D = vecProjectedPoints1D - norm(vecFprime*100);
	dblRefNorm = norm(vecFprime);
	
	if nargout > 2
		%shuffle
		matShuffled1 = nan(size(matRecentered1));
		matProjectedPointsShuffled = nan(size(matRecentered1,1),intIters);
		for intIter=1:intIters
			for intDim=1:size(matRecentered1,2)
				matShuffled1(:,intDim) = matRecentered1(randperm(intPoints),intDim);
			end
			
			%calculate projected points
			for intTrial=1:size(matShuffled1,1)
				vecPoint = matShuffled1(intTrial,:);
				vecOrth = ((vecPoint*vecFprime')/(vecFprime*vecFprime')).*vecFprime - vecFprime*100;
				matProjectedPointsShuffled(intTrial,intIter) = norm(vecOrth);
			end
		end
		matProjectedPointsShuffled = matProjectedPointsShuffled - norm(vecFprime*100);
	end
	if nargout > 3
		%recenter
		matRecentered2 = bsxfun(@minus,matPoints2,vecMu1);
		
		%calculate projected points
		vecProjectedPointsClass2 = nan(size(matRecentered2,1),1);
		for intTrial=1:size(matRecentered2,1)
			vecPoint = matRecentered2(intTrial,:);
			vecOrth = ((vecPoint*vecFprime')/(vecFprime*vecFprime')).*vecFprime - vecFprime*100;
			vecProjectedPointsClass2(intTrial) = norm(vecOrth);
		end
		vecProjectedPointsClass2 = vecProjectedPointsClass2 - norm(vecFprime*100);
	end
end
