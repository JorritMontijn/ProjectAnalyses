function matAlignmentWithX = doDimPrivateAnalysisCV(matActBin,cellNeuronsX,cellNeuronsY,cellMatX1,cellMatY1,cellMatX2,cellMatY2,dblLambda,boolShuffle)
	%% prep
	if ~exist('dblLambda','var') || isempty(dblLambda)
		dblLambda = 0;
	end
	if ~exist('boolShuffle','var') || isempty(boolShuffle)
		boolShuffle = false;
	end
	
	
	%% pre-allocate
	intBins = size(matActBin,1);
	intNeurons = size(matActBin,2);
	intTrials = size(matActBin,3);
	intStimTypes = size(cellMatX1,2);
	intResamplings = size(cellMatX1,1);
	
	%cellWeights = cell(intStimTypes,intResamplings);
	%cellSubspX = cell(intStimTypes,intResamplings);
	%cellSubspY = cell(intStimTypes,intResamplings);
	
	%get dom dim vec
	matAlignmentWithX = nan(intBins,intTrials,intStimTypes,intResamplings);
		
	%% run
	for intStimType=1:intStimTypes
		for intResampling=1:intResamplings
			%% msg
			fprintf('   Running resampling %d/%d [%s]\n',intResampling,intResamplings,getTime);
			
			%% run analyses
			%match random ordering of neurons
			[a,vecReorderX] = ismember(min(cellNeuronsX{intResampling}):max(cellNeuronsX{intResampling}),cellNeuronsX{intResampling});
			[a,vecReorderY] = ismember(min(cellNeuronsY{intResampling}):max(cellNeuronsY{intResampling}),cellNeuronsY{intResampling});
			
			% select data
			matX1 = cellMatX1{intResampling,intStimType}(:,vecReorderX);
			matY1 = cellMatY1{intResampling,intStimType}(:,vecReorderY);
			matX2 = cellMatX2{intResampling,intStimType}(:,vecReorderX);
			matY2 = cellMatY2{intResampling,intStimType}(:,vecReorderY);
			%get mu+sigma
			vecSigX1 = std(matX1,[],1);
			vecSigX2 = std(matX2,[],1);
			vecSigY1 = std(matY1,[],1);
			vecSigY2 = std(matY2,[],1);
			vecSig = mean(cat(1,vecSigX1,vecSigX2,vecSigY1,vecSigY2),1);
			vecMuX1 = mean(matX1,1);
			vecMuX2 = mean(matX2,1);
			vecMuY1 = mean(matY1,1);
			vecMuY2 = mean(matY2,1);
			vecMu = mean(cat(1,vecMuX1,vecMuX2,vecMuY1,vecMuY2),1);
			%pseudo z-score
			matX1Z = (matX1 - vecMu)./vecSig;
			matY1Z = (matX2 - vecMu)./vecSig;
			matX2Z = (matY1 - vecMu)./vecSig;
			matY2Z = (matY2 - vecMu)./vecSig;
			matActBinZ = (matActBin - vecMu)./vecSig; 
			
			if (any(range(matX1,1)==0) || any(range(matX2,1)==0) || any(range(matY1,1)==0) || any(range(matY2,1)==0)),continue;end
			%shuffle
			if boolShuffle
				matX2 = matX2(randperm(size(matX2,1)),:);
				matY2 = matY2(randperm(size(matX2,1)),:);
			end
			
			% get private spaces
			[matSubspace_PrivateX,matSubspace_PrivateY,matWeights,intDim,vecPerformance,cellSubspPrivX,cellSubspPrivY,cellWeights] = ...
				getSpacesPrivateCV(matX1,matY1,matX2,matY2,dblLambda);
			
			%calculate temporal progression
			matAlignmentX = getTempDynamicsInPrivateSpaces(matSubspace_PrivateX,matSubspace_PrivateY,matWeights,matActBin);
			matAlignmentWithX(:,:,intStimType,intResampling) = matAlignmentX';
		end
	end
	