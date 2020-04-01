function [matAlignmentWithX,matMagnitudeOfPop] = doDimPrivateAnalysisCV3(matActBin,matX,matY,dblLambda)
	%% prep
	if ~exist('dblLambda','var') || isempty(dblLambda)
		dblLambda = 0;
	end
	
	%% pre-allocate
	intBins = size(matActBin,1);
	intNeurons = size(matActBin,2);
	intTrials = size(matActBin,3);

	%% get separator
	matX = matX';
	matY = matY';
	intRepsPerType = size(matX,1);
	matTrainData = cat(1,matX,matY)';
	%get trial labels
	vecAggTrialTypes = cat(1,zeros(intRepsPerType,1),ones(intRepsPerType,1));
	vecAggRepVec = cat(2,1:(intRepsPerType),1:(intRepsPerType))';
	
	%get weights
	[dblPerformance,vecDecodedIndexCV,matPosteriorProbability,matWeights,dblMeanErrorDegs,matConfusion] = ...
		doCrossValidatedDecodingLR(matTrainData,vecAggTrialTypes,vecAggRepVec,dblLambda);
	
	%% decode bins
	% pre-allocate
	matAlignmentX = nan(intTrials,intBins);
	matMagnitude = nan(intTrials,intBins);
	
	% run
	for intBin=1:intBins
		%get activity matrix
		matThisAct = squeeze(matActBin(intBin,:,:));
		
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
	matAlignmentWithX = matAlignmentX';
	matMagnitudeOfPop = matMagnitude';
	
	