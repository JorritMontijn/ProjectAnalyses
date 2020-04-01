function [matSubspace_PrivateX,matSubspace_PrivateY,matWeights,intDim,vecPerformance,cellSubspPrivX,cellSubspPrivY,cellWeights] = getSpacesPrivateCV(matX1,matY1,matX2,matY2,dblLambda)
	%% prep
	if ~exist('dblLambda','var') || isempty(dblLambda)
		dblLambda = 1;
	end
	
	%% pre-allocate
	cellWeights = cell(1,size(matX1,2));
	
	%% get private space for X
	%check size of private y
	[vecR2PrivX,vecR2PrivX_NonCV,cellSubspPrivX] = getRedRankRegCV(matX1,matX1,matX2,matX2,dblLambda);
	
	%% get private space for Y
	%check size of private y
	[vecR2PrivY,vecR2PrivY_NonCV,cellSubspPrivY] = getRedRankRegCV(matY1,matY1,matY2,matY2,dblLambda);
	
	%% get subspaces and separator
	%get cut off points
	%intDimY = find(vecR2PrivY ./ vecR2PrivY(end) > 0.9,1);
	%intDimX = find(vecR2PrivX ./ vecR2PrivX(end) > 0.9,1);
	%intUseDim = max([intDimY intDimX]);
	intMaxDim = numel(cellSubspPrivX);
	vecPerformance = nan(1,intMaxDim);
	for intUseDim=1:intMaxDim
		%get basis vectors for xy-shared space
		matSubspace_PrivateX = cellSubspPrivX{intUseDim};
		
		%get basis vectors for xy-shared space
		matSubspace_PrivateY = cellSubspPrivY{intUseDim};
		
		%get subspace activity
		maxPrivXX1 = matX1*matSubspace_PrivateX;
		maxPrivXX2 = matX2*matSubspace_PrivateX;
		
		maxPrivXY1 = matY1*matSubspace_PrivateX;
		maxPrivXY2 = matY2*matSubspace_PrivateX;
		
		matAggX = [maxPrivXX1;maxPrivXX2;maxPrivXY1;maxPrivXY2]';
		
		%get subspace activity
		maxPrivYX1 = matX1*matSubspace_PrivateY;
		maxPrivYX2 = matX2*matSubspace_PrivateY;
		
		maxPrivYY1 = matY1*matSubspace_PrivateY;
		maxPrivYY2 = matY2*matSubspace_PrivateY;
		
		matAggY = [maxPrivYX1;maxPrivYX2;maxPrivYY1;maxPrivYY2]';
		
		intTrials = size(matAggX,2);
		vecActX = nan(1,intTrials);
		vecActY = nan(1,intTrials);
		for intTrial=1:intTrials
			vecActX(intTrial) = norm(matAggX(:,intTrial));
			vecActY(intTrial) = norm(matAggY(:,intTrial));
		end
		
		%% get separator
		matAgg = [vecActX;vecActY];
		intT = size(matAgg,2)/2;
		vecAggTrialTypes = cat(1,zeros(intT,1),ones(intT,1));
		vecAggRepVec = cat(2,1:(intT),1:(intT))';
		
		%decode
		[dblPerformance,vecDecodedIndexCV,matPosteriorProbability,matWeights,dblMeanErrorDegs,matConfusion] = ...
			doCrossValidatedDecodingLR(matAgg,vecAggTrialTypes,vecAggRepVec,1/10);
		vecPerformance(intUseDim) = dblPerformance;
		cellWeights{intUseDim} = matWeights;
		%[dblPerformance2,vecDecodedIndexCV2] = decodeLR(matWeights,matAgg,vecAggTrialTypes)
		%vecPerformance2(intUseDim) = dblPerformance2
		
	end
	[d,intDim]=find((vecPerformance - 0.5) > 0.9*(max(vecPerformance) - 0.5),1);
	if isempty(intDim),intDim=1;end
	
	%assign outputs
	matSubspace_PrivateX = cellSubspPrivX{intDim};
	matSubspace_PrivateY = cellSubspPrivY{intDim};
	matWeights = cellWeights{intDim};
	
	%dblPerformance
	
	%dblOptimLambda = 5e3;%HyperOptim('-doCrossValidatedDecodingLR',4,{matAgg,vecAggTrialTypes,1},[],1,1e4,'fminbnd')
	%dblOptimLambda2 = HyperOptim('-doCrossValidatedDecodingLR',4,{matAgg,vecAggTrialTypes,1},dblOptimLambda,[],[],'fminsearch')
	
end