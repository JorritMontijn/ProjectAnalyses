%% test case with known parameters
% ascertain the difference between Fisher I, mahalanobis distance, and
% predicted fisher information from distributions along F'

%% set parameters
intNeurons = 2; 
dblDistance = 1; %in absolute Euclidian space
dblOffdiag = 0.5; %between 0 (fully uniform; only on-diagonal elements) and 1 (all variability is along 1 (hyper)line only)
dblCosSimFPrime = 1; %is directionality parallel (1) or perpendicular (0) to F'?

%% build known distribution
vecMeanClass1 = zeros(1,intNeurons); %set class 1 to be at origin
vecMeanClass2 = dblDistance*(ones(1,intNeurons)/sqrt(sum(ones(1,intNeurons).^2))); %normalize Euclidian norm to 1, then multiply by requested distance

%build covariance matrix
matSigma = dblOffdiag*ones(intNeurons,intNeurons) + (1-dblOffdiag)*eye(intNeurons,intNeurons)

%% calculate mahalanobis distance
sParams.vecGroupSizes = numel(vecCells);
%sParams.intIters = 1;
matData = matModelRespP(vecCells,:);
[matMahalDists,matMahalDistsShuffled,vecGroupSizes] = getMahalCV(matModelRespP(vecCells,:),vecTrialStimType,sParams);

%% F' prime variance
%get the variance along F' compared to shuffled distribution, and compared
%to a random vector orthogonal to F'
vecF2 = vecFprime;
intN = numel(vecCells);
for intEl = 1:intN
	vecF2(intEl) = -vecFprime(intEl)*sum(vecFprime(~ismember(1:intN,1)).^2);
	if vecF2'*vecFprime < 0,break;end
end
[dblPercVarInDirOverShuffled,vecCIPCIDOS,dblCohensD] = getVarInDir(matR1',vecFprime',intResamplings);
[dblPercVarInDirOverShuffledOrthogonal,vecCIPCIDOSO,dblCohensDOrth] = getVarInDir(matR1',vecF2',intResamplings);

%% F' projection & predicted Fisher I
%get projection onto F'; what is distribution compared to cluster separation?
[vecProjectedPoints1D,dblRefNorm,matProjectedPointsShuffled,vecProjectedPointsClass2] = getProjOnLine(matR1',matR2',intResamplings);
vecIncreaseInDiffCorrs = (xstd(vecProjectedPoints1D,1) ./ xstd(matProjectedPointsShuffled,1)) * 100;
dblIncreaseInDiffCorrs = mean(vecIncreaseInDiffCorrs);
vecIncreaseInDiffCorrCI = getCI(vecIncreaseInDiffCorrs);


%% Fisher I
%get Fisher information
[dblPredA,matPredA,dblDprimeSquaredOffDiagonal,matDprimeSquared,matDprimeSquared_diagonal] = getSeparation(matData,vecTrialTypes,boolLinear)
