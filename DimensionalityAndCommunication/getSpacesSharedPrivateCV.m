function [vecR2XY,vecR2YX,matSubspace_Shared,matSubspace_PrivateX,matSubspace_PrivateY,matSeparator,vecR2PrivX,vecR2PrivY] = getSpacesSharedPrivateCV(matX1,matY1,matX2,matY2,dblLambda)
	%% prep
	intMaxDim = size(matY1,2);
	
	if ~exist('dblLambda','var') || isempty(dblLambda)
		dblLambda = 1;
	end
	
	%pre-allocate
	cellCommSubspXY = cell(1,intMaxDim);
	cellCommSubspYX = cell(1,intMaxDim);
	vecR2XY = nan(1,intMaxDim);
	vecR2YX = nan(1,intMaxDim);
	vecR2_NonCV = nan(1,intMaxDim);
	
	vecR2PrivX = nan(1,intMaxDim);
	cellSubspPrivX = cell(1,intMaxDim);
		
	vecR2PrivY = nan(1,intMaxDim);
	cellSubspPrivY = cell(1,intMaxDim);
	
	%% pre-calculations
	matCovX1 = cov(matX1);
	matCovX2 = cov(matX2);
	matCovY1 = cov(matY1);
	matCovY2 = cov(matY2);
	%% predictions
	for intDim=1:intMaxDim
		%% prediction X-Y
		%get subspace
		[matUseSubspaceXY, dblMSE, intRankOutT, sSuppOut] = doRidgeRRR(matX1,matY1,intDim,dblLambda);
		vecR2_NonCV(intDim) = sSuppOut.dblR2;
		
		%get bias
		matMu = mean(matY2)' - (matUseSubspaceXY')*(mean(matX2)');
		
		% get prediction
		matY_pred = repmat(matMu',[size(matY2,1) 1]) + matX2*matUseSubspaceXY;
		
		% compute MSE
		matErr = (matY2 - matY_pred).^2;
		
		%get R^2
		vecMu = mean(matY2);
		dblSSRes = sum(sum(matErr));
		dblSSTot = sum(sum(bsxfun(@minus,matY2,vecMu).^2));
		vecR2XY(intDim) = 1 - dblSSRes / dblSSTot;
		%save subspace
		cellCommSubspXY{intDim} = matUseSubspaceXY;
		
		%% prediction Y-X
		%get subspace
		[matUseSubspaceYX, dblMSE, intRankOutT, sSuppOut] = doRidgeRRR(matY1,matX1,intDim,dblLambda);
		vecR2_NonCV(intDim) = sSuppOut.dblR2;
		
		%get bias
		matMu = mean(matX2)' - (matUseSubspaceYX')*(mean(matY2)');
		
		% get prediction
		matY_pred = repmat(matMu',[size(matX2,1) 1]) + matY2*matUseSubspaceYX;
		
		% compute MSE
		matErr = (matX2 - matY_pred).^2;
		
		%get R^2
		vecMu = mean(matX2);
		dblSSRes = sum(sum(matErr));
		dblSSTot = sum(sum(bsxfun(@minus,matX2,vecMu).^2));
		vecR2YX(intDim) = 1 - dblSSRes / dblSSTot;
		%save subspace
		cellCommSubspYX{intDim} = matUseSubspaceYX;
		
	end
	%% get shared comm subspace, and private subspaces
	%choose # of dims
	intDimXY = find(vecR2XY ./ vecR2XY(end) > 0.9,1);
	intDimYX = find(vecR2YX ./ vecR2YX(end) > 0.9,1);
	intUseDims = min([intDimXY intDimYX]);
	intPrivDims = intMaxDim-intUseDims;
	
	%get basis vectors for xy-shared space
	matUseSubspaceXY = cellCommSubspXY{intUseDims};
	M_xy = matCovX1*matUseSubspaceXY;
	[U_xy,S,V_xy] = svd(M_xy);
	
	%get basis vectors for yx-shared space
	matUseSubspaceYX = cellCommSubspYX{intUseDims};
	M_yx = matCovY1*matUseSubspaceYX;
	[U_yx,S,V_yx] = svd(M_yx);
	
	%combine
	matSubspace_Shared = [[U_xy(:,1:intUseDims) U_yx(:,1:intUseDims)] * [V_xy(:,1:intUseDims) V_yx(:,1:intUseDims)]'];
	
	%just a sanity check to see the predictability after removing all shared dimensions
	%{
	%M*Q is zero matrix
	%Q'*Q is identity ([intMaxDim - intDim]-dimensional)
	matX1_Shared = matX1 * matSubspace_Shared;%matSubspaceShared;
	matX2_Shared = matX2 * matSubspace_Shared;%matSubspaceShared;
	
	%ridge regression between reduced space X and Y
	B_ridge_Z = (matX1_Shared' * matX1_Shared + dblLambda*eye(size(matX1_Shared,2))) \ (matX1_Shared' * matY1); %left-divide is same as inverse and multiplication
	
	%test
	matSharedRemX2_pred = matX2_Shared * B_ridge_Z;
	
	% compute MSE
	matErr = (matY2 - matSharedRemX2_pred).^2;
	%matErr = (matY2 - matY_pred).^2;
	
	%get R^2
	vecMu = mean(matY2);
	dblSSRes = sum(sum(matErr));
	dblSSTot = sum(sum(bsxfun(@minus,matY2,vecMu).^2));
	dblR2_shared = 1 - dblSSRes / dblSSTot;
	dblR2_orig = vecR2XY(intUseDims);
	%}
	
	%% get private space for X
	%remove shared dimensions
	MX1 = matSubspace_Shared;%*matCovX1;
	[U,S,V] = svd(MX1);
	QX = V(:,(intUseDims+1):end);
	
	%M*Q is zero matrix
	%Q'*Q is identity ([intMaxDim - intDim]-dimensional)
	matX1_SharedRem = matX1 * QX;
	matX2_SharedRem = matX2 * QX;
	
	%check size of private y
	[vecR2PrivX,vecR2PrivX_NonCV,cellSubspPrivX] = getRedRankRegCV(matX1_SharedRem,matX1_SharedRem,matX2_SharedRem,matX2_SharedRem,dblLambda);
	intDimX = find(vecR2PrivX ./ vecR2PrivX(end) > 0.9,1);
	
	%get basis vectors for xy-shared space
	matSubspace_PrivateX = QX*cellSubspPrivX{intDimX} * QX';
	
	%% get private space for Y
	%remove shared dimensions
	MY1 = matSubspace_Shared;%*matCovY1;
	[U,S,V] = svd(MY1);
	%QY = cat(2,zeros(intMaxDim,intUseDims),V(:,(intUseDims+1):end));
	QY = V(:,(intUseDims+1):end);
	%M*Q is zero matrix
	%Q'*Q is identity ([intMaxDim - intDim]-dimensional)
	matY1_SharedRem = matY1 * QY;
	matY2_SharedRem = matY2 * QY;
	
	%check size of private y
	[vecR2PrivY,vecR2PrivY_NonCV,cellSubspPrivY] = getRedRankRegCV(matY1_SharedRem,matY1_SharedRem,matY2_SharedRem,matY2_SharedRem,dblLambda);
	intDimY = find(vecR2PrivY ./ vecR2PrivY(end) > 0.9,1);
	
	%get basis vectors for xy-shared space
	matSubspace_PrivateY = QY*cellSubspPrivY{intDimY}*QY';
	
	%% get separator
	matAgg = [matX1;matX2;matY1;matY2]';
	intT = size(matAgg,2)/2;
	vecAggTrialTypes = cat(1,zeros(intT,1),ones(intT,1));
	vecAggRepVec = cat(2,1:(intT),1:(intT))';
	
	%decode
	[dblPerformance,vecDecodedIndexCV,matPosteriorProbability,matWeights,dblMeanErrorDegs,matConfusion] = ...
		doCrossValidatedDecodingLR(matAgg,vecAggTrialTypes,vecAggRepVec,1/10);
	matSeparator = matWeights(1:(end-1),:);
	
	[dblPerformance2,vecDecodedIndexCV2] = decodeLR(matWeights,matAgg,vecAggTrialTypes)

	%dblPerformance
	
	%dblOptimLambda = 5e3;%HyperOptim('-doCrossValidatedDecodingLR',4,{matAgg,vecAggTrialTypes,1},[],1,1e4,'fminbnd')
	%dblOptimLambda2 = HyperOptim('-doCrossValidatedDecodingLR',4,{matAgg,vecAggTrialTypes,1},dblOptimLambda,[],[],'fminsearch')
	
end