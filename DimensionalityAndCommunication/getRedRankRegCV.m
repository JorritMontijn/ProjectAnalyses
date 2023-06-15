function [vecR2,vecR2_NonCV,cellB_rrr] = getRedRankRegCV(matX_Train,matY_Train,matX_Test,matY_Test,dblLambda)
	%% prep
	intMaxDim = size(matY_Train,2);
	vecR2 = nan(1,intMaxDim);
	vecR2_NonCV = nan(1,intMaxDim);
	
	if ~exist('dblLambda','var') || isempty(dblLambda)
		dblLambda = 0;
	end
	intType = 2;
	
	%pre-allocate
	cellB_rrr = cell(1,intMaxDim);
	
	%% prediction
	for intDim=1:intMaxDim
		if intType == 1
			%get subspace
			[matUseSubspace, dblMSE, intRankOutT, sSuppOut] = doRdRankReg(matX_Train, matY_Train, intDim);
			
			%suppress warning
			[dummy,strWID]=lastwarn;
			if strcmp(strWID,'MATLAB:nearlySingularMatrix')
				warning('off','MATLAB:nearlySingularMatrix')
			end
			
			%project activity into PC space
			matTrainZ = matX_Train*matUseSubspace + repmat(sSuppOut.matMu',[size(matX_Train,1) 1]);
			
			%ridge regression between reduced space X and Y
			B_ridge_Z = (matTrainZ' * matTrainZ + dblLambda*eye(size(matTrainZ,2))) \ (matTrainZ' * matY_Train); %left-divide is same as inverse and multiplication
			
			%test
			matTestZ = matX_Test*matUseSubspace + repmat(sSuppOut.matMu',[size(matX_Test,1) 1]);
			Y_pred_Z = matTestZ * B_ridge_Z;
			
			% compute R^2
			vecMu = mean(matY_Test);
			dblSSRes_Z = sum(sum((matY_Test - Y_pred_Z).^2));
			dblSSTot = sum(sum(bsxfun(@minus,matY_Test,vecMu).^2));
			dblR2_Z = 1 - dblSSRes_Z / dblSSTot;
			vecR2(intDim) = dblR2_Z;
			
			%save
			if nargout > 2
				cellB_rrr{intDim} = B_ridge_Z;
			end
		elseif intType == 2
			%get subspace
			[matUseSubspace, dblMSE, intRankOutT, sSuppOut] = doRidgeRRR(matX_Train,matY_Train,intDim,dblLambda);
			vecR2_NonCV(intDim) = sSuppOut.dblR2;
			
			%get bias
			matMu = mean(matY_Test)' - (matUseSubspace')*(mean(matX_Test)');
			
			% get prediction
			matY_pred = repmat(matMu',[size(matY_Test,1) 1]) + matX_Test*matUseSubspace;
			
			% compute MSE
			matErr = (matY_Test - matY_pred).^2;
			
			%get R^2
			vecMu = mean(matY_Test);
			dblSSRes = sum(sum(matErr));
			dblSSTot = sum(sum(bsxfun(@minus,matY_Test,vecMu).^2));
			vecR2(intDim) = 1 - dblSSRes / dblSSTot;
			
			%save
			if nargout > 2
				cellB_rrr{intDim} = matUseSubspace;
			end
		end
	end
end