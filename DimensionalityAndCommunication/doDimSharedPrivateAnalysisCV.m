function [matAvS,matAvPX,matAvPY,matAggS,matAggPX,matAggPY,matNorms,cellSubspS,cellSubspX,cellSubspY] = doDimSharedPrivateAnalysisCV(matActBin,cellNeuronsX,cellNeuronsY,cellMatX1,cellMatY1,cellMatX2,cellMatY2,dblLambda,boolShuffle)
	%% prep
	if ~exist('dblLambda','var') || isempty(dblLambda)
		dblLambda = 0;
	end
	if ~exist('boolShuffle','var') || isempty(boolShuffle)
		boolShuffle = false;
	end
	
	
	%% pre-allocate
	intBins = size(matActBin,1);
	intTrials = size(matActBin,3);
	intNeurons = size(cellMatX1{1},2);
	intStimTypes = size(cellMatX1,2);
	intResamplings = size(cellMatX1,1);
	
	cellSubspS = cell(intStimTypes,intResamplings);
	cellSubspX = cell(intStimTypes,intResamplings);
	cellSubspY = cell(intStimTypes,intResamplings);
	
	%get dom dim vec
	matNorms = nan(6,intStimTypes,intResamplings);
	matAvS = nan(intBins,intStimTypes,intResamplings);
	matAvPX = nan(intBins,intStimTypes,intResamplings);
	matAvPY = nan(intBins,intStimTypes,intResamplings);
	matAggS = nan(intBins,intTrials,intStimTypes,intResamplings);
	matAggPX = nan(intBins,intTrials,intStimTypes,intResamplings);
	matAggPY = nan(intBins,intTrials,intStimTypes,intResamplings);
	%normalize activity matrix
	matAvAct = mean(matActBin,1);
	vecSd = std(matAvAct,[],3);
	vecMean = mean(matAvAct,3);
	matActBinZ = bsxfun(@rdivide,bsxfun(@minus,matActBin,vecMean),vecSd);
	
	%% run
	for intStimType=1:intStimTypes
		for intResampling=1:intResamplings
			%% run analyses
			%match random ordering of neurons
			[a,vecReorderX] = ismember(min(cellNeuronsX{intResampling}):max(cellNeuronsX{intResampling}),cellNeuronsX{intResampling});
			[a,vecReorderY] = ismember(min(cellNeuronsY{intResampling}):max(cellNeuronsY{intResampling}),cellNeuronsY{intResampling});
			
			% select data
			%{
			matX1 = zscore(cellMatX1{intResampling,intStimType}(:,vecReorderX),[],1);
			matY1 = zscore(cellMatY1{intResampling,intStimType}(:,vecReorderY),[],1);
			matX2 = zscore(cellMatX2{intResampling,intStimType}(:,vecReorderX),[],1);
			matY2 = zscore(cellMatY2{intResampling,intStimType}(:,vecReorderY),[],1);
			%}
			matX1 = cellMatX1{intResampling,intStimType}(:,vecReorderX);
			matY1 = cellMatY1{intResampling,intStimType}(:,vecReorderY);
			matX2 = cellMatX2{intResampling,intStimType}(:,vecReorderX);
			matY2 = cellMatY2{intResampling,intStimType}(:,vecReorderY);
			if (any(range(matX1,1)==0) || any(range(matX2,1)==0) || any(range(matY1,1)==0) || any(range(matY2,1)==0)),continue;end
			%shuffle
			if boolShuffle
				matX2 = matX2(randperm(size(matX2,1)),:);
				matY2 = matY2(randperm(size(matX2,1)),:);
			end
			
			% get shared and private spaces
			[vecR2XY,vecR2YX,matSubspace_Shared,matSubspace_PrivateX,matSubspace_PrivateY,matSeparator,vecR2PrivX,vecR2PrivY] = ...
				getSpacesSharedPrivateCV(matX1,matY1,matX2,matY2,dblLambda);
			
			%assign
			if nargout > 7
				cellSubspS{intStimType,intResampling} = matSubspace_Shared;
				cellSubspX{intStimType,intResampling} = matSubspace_PrivateX;
				cellSubspY{intStimType,intResampling} = matSubspace_PrivateY;
			end
			
			%calculate norms
			dblXS = mean(mean(abs(matX2 * matSubspace_Shared)));
			dblXPX = mean(mean(abs(matX2 * matSubspace_PrivateX)));
			dblXPY = mean(mean(abs(matX2 * matSubspace_PrivateY)));
			dblYS = mean(mean(abs(matY2 * matSubspace_Shared)));
			dblYPX = mean(mean(abs(matY2 * matSubspace_PrivateX)));
			dblYPY = mean(mean(abs(matY2 * matSubspace_PrivateY)));
			matNorms(:,intStimType,intResampling) = ...
				[dblXS dblXPX dblXPY dblYS dblYPX dblYPY];
			
			%calculate temporal progression in spaces from t=0 - t=0.4
			matActT = mean(matActBinZ,3);
			matAvS(:,intStimType,intResampling) = mean(abs(matActT * matSubspace_Shared),2);
			matAvPX(:,intStimType,intResampling) = mean(abs(matActT * matSeparator(:,1)),2);
			matAvPY(:,intStimType,intResampling) = mean(abs(matActT * matSeparator(:,2)),2);
			%[Time x Cell x Trial]
			for intTrial=1:intTrials
				matActT = matActBinZ(:,:,intTrial);
				matAggS(:,intTrial,intStimType,intResampling) = mean(abs(matActT * matSubspace_Shared),2);
				matAggPX(:,intTrial,intStimType,intResampling) = mean(abs(matActT * matSeparator(:,1)),2);
				matAggPY(:,intTrial,intStimType,intResampling) = mean(abs(matActT * matSeparator(:,2)),2);
			end
		end
	end
	
	%% build output
	sOut = struct;
	
	return
	
	%%
	
	px=mean(matAggPX,4);
	py=mean(matAggPY,4);
	
	