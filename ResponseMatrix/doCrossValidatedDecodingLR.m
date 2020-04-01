function [dblPerformanceCV,vecDecodedIndexCV,matPosteriorProbabilityCV,matWeights,dblMeanErrorDegs,matConfusion] = doCrossValidatedDecodingLR(matData,vecTrialTypes,dblLambda)
	%UNTITLED9 Summary of this function goes here
	%   Detailed explanation goes here

	%check which kind of cross-validation
	if nargin < 3 || isempty(dblLambda)
		dblLambda = 0.01;
	end
	intVerbose = 0;
	
	%get number of trials
	vecTrialTypes = vecTrialTypes(:);
	intTrials = numel(vecTrialTypes);
	
	%check if matData is [trial x neuron] or [neuron x trial]
	if size(matData,1) == intTrials && size(matData,2) == intTrials
		%number of neurons and trials is the same
		warning([mfilename ':SameNeuronsTrials'],'Number of neurons and trials is identical; please double check the proper orientation of [intNeurons x intTrials]');
	elseif size(matData,1) == intTrials
		%rotate
		matData = matData';
	elseif size(matData,2) == intTrials
		%size is correct
	else
		error([mfilename ':SameNeuronsTrials'],'Size of matData and vecTrialTypes do not match');
	end
	intNeurons = size(matData,1);
	vecUniqueTrialTypes = unique(vecTrialTypes);
	intStimTypes = length(vecUniqueTrialTypes);
	vecTrialTypeIdx = label2idx(vecTrialTypes);
	intReps = intTrials/intStimTypes;
	
	%get weights
	[matWeights, dblLLH] = doMnLogReg(matData,vecTrialTypeIdx,dblLambda);
	
	%get performance
	matDataPlusLin = [matData; ones(1,size(matData,2))];
	matActivation = matWeights'*matDataPlusLin;
	matPosteriorProbabilityCV = exp(bsxfun(@minus,matActivation,logsumexp(matActivation,1))); %softmax
	[dummy, vecDecodedIndexCV] = max(matPosteriorProbabilityCV,[],1);
	
	%decoding accuracy
	dblPerformanceCV = sum(vecDecodedIndexCV(:)==vecTrialTypeIdx(:))/intTrials;
	
	%error
	if nargout > 4
		vecDecodedValuesCV = vecUniqueTrialTypes(vecDecodedIndexCV);
		dblMeanErrorRads = mean(abs(circ_dist(vecDecodedValuesCV,vecTrialTypes)));
		dblMeanErrorDegs = rad2ang(dblMeanErrorRads);
	end
	
	%confusion matrix;
	if nargout > 5
		matConfusion = getFillGrid(zeros(intStimTypes),vecDecodedIndexCV(:),vecTrialTypeIdx(:),ones(intTrials,1))/intReps;
		%imagesc(matConfusion,[0 1]);colormap(hot);axis xy;colorbar;
	end
	
end

