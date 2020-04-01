function [matAlignmentX] = getTempDynamicsInPrivateSpaces(matSubspace_PrivateX,matSubspace_PrivateY,matWeights,matActBin)
	
	%% pre-allocate
	intBins = size(matActBin,1);
	intTrials = size(matActBin,3);
	matAlignmentX = nan(intTrials,intBins);
	
	%% run
	for intBin=1:intBins
		%get activity matrix
		matThisAct = squeeze(matActBin(intBin,:,:))';
			
		%get subspace activity
		maxPrivX = matThisAct*matSubspace_PrivateX;
		maxPrivY = matThisAct*matSubspace_PrivateY;
		
		intTrials = size(maxPrivX,1);
		vecActX = nan(1,intTrials);
		vecActY = nan(1,intTrials);
		for intTrial=1:intTrials
			vecActX(intTrial) = norm(maxPrivX(intTrial,:));
			vecActY(intTrial) = norm(maxPrivY(intTrial,:));
		end
		
		%% get separator
		matAgg = [vecActX;vecActY;ones(1,intTrials)];
		
		%decode
		matActivation = matWeights'*matAgg;
		%matPosteriorProbability = exp(bsxfun(@minus,matActivation,logsumexp(matActivation,1)))'; %softmax
		%vecPostProbX = matPosteriorProbability(:,1);
		%vecAlignmentX(intBin) = mean(vecPostProbX);
		matAlignmentX(:,intBin) = matActivation(1,:)-matActivation(2,:);
	end
end