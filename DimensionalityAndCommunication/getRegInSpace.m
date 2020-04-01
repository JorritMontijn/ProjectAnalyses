function vecR2 = getRegInSpace(matX,matSubspace,matY,dblLambda)
	%% prep
	intMaxDim = size(matSubspace,2);
	vecR2 = nan(1,intMaxDim);
	
	if ~exist('dblLambda','var') || isempty(dblLambda)
		dblLambda = 0;
	end
	
	%% prediction
	for intDim=1:intMaxDim
		%select
		matUseSubspace = matSubspace(:,1:intDim)';
		
		%project activity into PC space
		matZ = (matUseSubspace * matX')';
		
		%ridge regression between reduced space X and Y
		B_ridge_Z = (matZ' * matZ + dblLambda*eye(size(matZ,2))) \ (matZ' * matY); %left-divide is same as inverse and multiplication
		Y_pred_Z = matZ * B_ridge_Z;
		
		% compute R^2
		vecMu = mean(matY);
		dblSSRes_Z = sum(sum((matY - Y_pred_Z).^2));
		dblSSTot = sum(sum(bsxfun(@minus,matY,vecMu).^2));
		dblR2_Z = 1 - dblSSRes_Z / dblSSTot;
		vecR2(intDim) = dblR2_Z;
	end
end