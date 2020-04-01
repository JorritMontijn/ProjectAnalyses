function vecR2 = getRegInSpaceCV(matTrainX,matSubspace,matTrainY,matTestX,matTestY,dblLambda)
	%% prep
	intMaxDim = size(matSubspace,2);
	vecR2 = nan(1,intMaxDim);
	
	if ~exist('dblLambda','var') || isempty(dblLambda)
		dblLambda = 0;
	end
	
	%% prediction
	for intDim=1:intMaxDim
		%select
		matUseSubspace = matSubspace(:,1:intDim);
		
		%project activity into PC space
		matTrainZ = matTrainX * matUseSubspace;
		
		%ridge regression between reduced space X and Y
		B_ridge_Z = (matTrainZ' * matTrainZ + dblLambda*eye(size(matTrainZ,2))) \ (matTrainZ' * matTrainY); %left-divide is same as inverse and multiplication
		
		%test
		matTestZ = matTestX * matUseSubspace;
		Y_pred_Z = matTestZ * B_ridge_Z;
		
		% compute R^2
		vecMu = mean(matTestY);
		dblSSRes_Z = sum(sum((matTestY - Y_pred_Z).^2));
		dblSSTot = sum(sum(bsxfun(@minus,matTestY,vecMu).^2));
		dblR2_Z = 1 - dblSSRes_Z / dblSSTot;
		vecR2(intDim) = dblR2_Z;
	end
end