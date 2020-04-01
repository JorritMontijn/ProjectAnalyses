function [matAlignmentWithX,matMagnitudeOfPop] = doDimPrivateAnalysisCV2(matActBin,cellNeuronsX,cellNeuronsY,cellMatX1,cellMatY1,cellMatX2,cellMatY2,cellTrials1,cellTrials2,dblLambda,boolShuffle)
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
	intTrials = size(cellTrials2{1},2);
	intStimTypes = size(cellMatX1,2);
	intResamplings = size(cellMatX1,1);
	
	%cellWeights = cell(intStimTypes,intResamplings);
	%cellSubspX = cell(intStimTypes,intResamplings);
	%cellSubspY = cell(intStimTypes,intResamplings);
	
	%get dom dim vec
	matAlignmentWithX = nan(intBins,intTrials,intStimTypes,intResamplings);
	matMagnitudeOfPop = nan(intBins,intTrials,intStimTypes,intResamplings);
	
	%% check null
	indRemX1 = false(1,intNeurons);
	indRemX2 = false(1,intNeurons);
	indRemY1 = false(1,intNeurons);
	indRemY2 = false(1,intNeurons);
	% run
	for intStimType=1:intStimTypes
		for intResampling=1:intResamplings
			%% msg
			%fprintf('   Running resampling %d/%d [%s]\n',intResampling,intResamplings,getTime);
			
			%% run analyses
			%match random ordering of neurons
			[a,vecReorderX] = ismember(min(cellNeuronsX{intResampling}):max(cellNeuronsX{intResampling}),cellNeuronsX{intResampling});
			[a,vecReorderY] = ismember(min(cellNeuronsY{intResampling}):max(cellNeuronsY{intResampling}),cellNeuronsY{intResampling});
			
			matX1 = cellMatX1{intResampling,intStimType}(:,vecReorderX);
			matY1 = cellMatY1{intResampling,intStimType}(:,vecReorderY);
			matX2 = cellMatX2{intResampling,intStimType}(:,vecReorderX);
			matY2 = cellMatY2{intResampling,intStimType}(:,vecReorderY);
			
			indRemX1 = indRemX1 | range(matX1)==0;
			indRemX2 = indRemX2 | range(matY1)==0;
			indRemY1 = indRemY1 | range(matX2)==0;
			indRemY2 = indRemY2 | range(matY2)==0;
		end
	end
	indRemNeurons = indRemX1 | indRemX2 | indRemY1 | indRemY2;
	matActBin = matActBin(:,~indRemNeurons,:);
			
	%% run
	for intStimType=1:intStimTypes
		for intResampling=1:intResamplings
			%% msg
			%fprintf('   Running resampling %d/%d [%s]\n',intResampling,intResamplings,getTime);
			
			%% run analyses
			%match random ordering of neurons
			[a,vecReorderX] = ismember(min(cellNeuronsX{intResampling}):max(cellNeuronsX{intResampling}),cellNeuronsX{intResampling});
			[a,vecReorderY] = ismember(min(cellNeuronsY{intResampling}):max(cellNeuronsY{intResampling}),cellNeuronsY{intResampling});
			
			% select data
			matX1 = cellMatX1{intResampling,intStimType}(:,vecReorderX);
			matY1 = cellMatY1{intResampling,intStimType}(:,vecReorderY);
			matX2 = cellMatX2{intResampling,intStimType}(:,vecReorderX);
			matY2 = cellMatY2{intResampling,intStimType}(:,vecReorderY);
			%remove neurons
			matX1 = matX1(:,~indRemNeurons);
			matY1 = matY1(:,~indRemNeurons);
			matX2 = matX2(:,~indRemNeurons);
			matY2 = matY2(:,~indRemNeurons);
			
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
			
			%% get data splits
			vecTrainTrials = cellTrials1{intResampling};
			vecTestTrials = cellTrials2{intResampling};
			
			%% get separator
			intRepsPerType = size(matX1,1);
			matTrainData = cat(1,matX1,matY1)';
			%get trial labels
			vecAggTrialTypes = cat(1,zeros(intRepsPerType,1),ones(intRepsPerType,1));
			vecAggRepVec = cat(2,1:(intRepsPerType),1:(intRepsPerType))';
			
			%get weights
			[dblPerformance,vecDecodedIndexCV,matPosteriorProbability,matWeights,dblMeanErrorDegs,matConfusion] = ...
				doCrossValidatedDecodingLR(matTrainData,vecAggTrialTypes,vecAggRepVec,dblLambda);
			
			%% decode bins
			% pre-allocate
			intBins = size(matActBin,1);
			intTrials = numel(vecTestTrials);
			matAlignmentX = nan(intTrials,intBins);
			matMagnitude = nan(intTrials,intBins);
			
			% run
			for intBin=1:intBins
				%get activity matrix
				matThisAct = squeeze(matActBin(intBin,:,vecTestTrials));
				
				%% get separator
				matDataPlusLin = [matThisAct; ones(1,size(matThisAct,2))];
				
				%decode
				matActivation = matWeights'*matDataPlusLin;
				matAlignmentX(:,intBin) = matActivation(1,:)-matActivation(2,:);
				for intTrial=1:size(matThisAct,2)
					matMagnitude(intTrial,intBin) = norm(matThisAct(:,intTrial));
				end
			end
			%calculate temporal progression
			matAlignmentWithX(:,:,intStimType,intResampling) = matAlignmentX';
			matMagnitudeOfPop(:,:,intStimType,intResampling) = matMagnitude';
		end
	end
	