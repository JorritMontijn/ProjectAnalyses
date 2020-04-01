function [vecR2_PredRem,vecR2,vecR2_NonCV,cellQ,cellB_PredRem] = getRemPredDimsCV(matX1,matY1,matX2,matY2,dblLambda)
	%% prep
	intMaxDim = size(matY1,2);
	
	if ~exist('dblLambda','var') || isempty(dblLambda)
		dblLambda = 0;
	end
	
	%pre-allocate
	cellB_PredRem = cell(1,intMaxDim);
	cellQ = cell(1,intMaxDim);
	vecR2_PredRem = nan(1,intMaxDim);
	vecR2 = nan(1,intMaxDim);
	vecR2_NonCV = nan(1,intMaxDim);
	
	%% pre-calculations
	matCovX1 = cov(matX1);
	
	%% prediction
	for intDim=1:intMaxDim
		%get subspace
		[matUseSubspace, dblMSE, intRankOutT, sSuppOut] = doRidgeRRR(matX1,matY1,intDim,dblLambda);
		vecR2_NonCV(intDim) = sSuppOut.dblR2;
		
		%get bias
		matMu = mean(matY2)' - (matUseSubspace')*(mean(matX2)');
		
		% get prediction
		matY_pred = repmat(matMu',[size(matY2,1) 1]) + matX2*matUseSubspace;
		
		% compute MSE
		matErr = (matY2 - matY_pred).^2;
		
		%get R^2
		vecMu = mean(matY2);
		dblSSRes = sum(sum(matErr));
		dblSSTot = sum(sum(bsxfun(@minus,matY2,vecMu).^2));
		vecR2(intDim) = 1 - dblSSRes / dblSSTot;
		
		
		%remove predictive dimensions
		M = matUseSubspace*matCovX1;
		[U,S,V] = svd(M);
		Q = V(:,(intDim+1):end);
		%M*Q is zero matrix
		%Q'*Q is identity ([intMaxDim - intDim]-dimensional)
		matX1_PredRem = matX1 * Q;
		matX2_PredRem = matX2 * Q;
		
		%ridge regression between reduced space X and Y
		B_ridge_Z = (matX1_PredRem' * matX1_PredRem + dblLambda*eye(size(matX1_PredRem,2))) \ (matX1_PredRem' * matY1); %left-divide is same as inverse and multiplication
		
		%test
		matPredRemY2_pred = matX2_PredRem * B_ridge_Z;
	
		% compute MSE
		matErr = (matY2 - matPredRemY2_pred).^2;
		
		%get R^2
		vecMu = mean(matY2);
		dblSSRes = sum(sum(matErr));
		dblSSTot = sum(sum(bsxfun(@minus,matY2,vecMu).^2));
		vecR2_PredRem(intDim) = 1 - dblSSRes / dblSSTot;
		
		%save
		if nargout > 3
			cellB_PredRem{intDim} = B_ridge_Z;
		end
		%save
		if nargout > 4
			cellQ{intDim} = Q;
		end
		
	end
end