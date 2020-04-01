close all;clearvars;
boolReverseType=1;
clearvars -except boolReverseType
%strSimulation = 'xAreaDistributed_C2ExpRet32Col180Ori23Noise5RP_2018-04-30';
%strSimulation = 'xAreaDistributed_C2ExpRet32Col180Ori2Noise5RP_2018-04-30';
%strSimulation = 'xAreaDistributed_C2ExpRet32Col180Ori23Noise1RP_2018-04-30';
%strSimulation = 'xAreaDistributed_C2ExpRet32Col180Ori2Noise1RP_2018-04-30';
strSimulation = 'xAreaExperiment_106r001p26-t_2018-02-21';
%strSimulation = 'xAreaExperiment_107l003p143-t_2018-02-21';
		
boolLoadSpikes = false;
boolLoadLGN = false;
boolLoadSync = false;
boolLoadTrialSplit = true;

%% RUN: #header
strFigDir = 'D:\Data\Results\Block5\';
runAnalysisHeader;

cellIn = strsplit(strSimulation,'_');
strFiller = cellIn{1};
strType = cell2mat(cellIn(2:(end-1)));
strDate = cellIn{end};
strTag = [strType '_' strDate];

%% set globals
global gMatX; %#ok<*TLEV>
global gMatY;
global gIntFolds;

%% set parameters
if ~exist('boolReverseType','var'),boolReverseType = false;end
intUseIndepTrials = 400;
intWithinArea = 3;
boolShowThree = true; %show only pred, rand, full CV
boolDoRegularization = true;
intFoldK = 10;
if ~exist('intUseRepetitionMax','var'),intUseRepetitionMax = 4000;end
intSizeX = 80;
intSizeY = 30;
boolSaveFigs = true;
dblLambdaInfo = 1;
dblLambdaPred = 0;
dblNoiseLevel = 0;
intMaxDimAnal = 30;
intIters = 10;
intRandIters = 5;
vecRank = 1:intSizeY;
cellStrArea = {'V1 from V1','V2 from V2','V2 from V1'};
strPredArea = cellStrArea{intWithinArea};

%% get data
fprintf('Transforming data [%s]\n',getTime);
[vecTrialStimTypeUnsplit,matRespUnsplit,vecTrialStimTypeSplit,matRespSplitNorm,matRespSplit,matSplitTrialIdx] = ...
	prepSplitData(matModelRespTS3,vecTrialStimType);

%% transpose
matRespUnsplit = matRespUnsplit';
matRespSplit = matRespSplit';
vecTrialStimTypeUnsplit2 = vecTrialStimTypeUnsplit;
vecTrialStimTypeSplit2 = vecTrialStimTypeSplit;

%% start
vecUseStimTypes = unique(round(vecTrialStimType));
indUseTrials = ismember(vecTrialStimType,vecUseStimTypes);

%% set default parameters
dblDiffTheta = range(vecStimTypeOris(vecUseStimTypes([1 2])));

%% get data splits
%set parameters
sParamsAnalSplit=struct;
sParamsAnalSplit.intSizeX = intSizeX;
sParamsAnalSplit.intSizeY = intSizeY;
sParamsAnalSplit.intResamplings = intIters;
sParamsAnalSplit.vecCellArea = vecCellArea;
sParamsAnalSplit.intWithinArea = intWithinArea;
sParamsAnalSplit.vecUseStimTypes = vecUseStimTypes;
sParamsAnalSplit.intMaxReps = intUseIndepTrials; %unused
%get splits
[cellMatX,cellNeuronsX,cellMatY,cellNeuronsY,cellTrials] = doDimDataSplits(matRespUnsplit',vecTrialStimTypeUnsplit,sParamsAnalSplit);
fprintf('  Created %dx%d data splits [%s]\n',size(cellNeuronsX),getTime);

matComb = [cellMatX{1,1};cellMatX{1,2}];
vecCombs = [ones(size(cellMatX{1,1},1),1) 2*ones(size(cellMatX{1,2},1),1)];
[~,~,dblInfoRes1]=getSeparation(matComb,vecCombs,false,dblDiffTheta);
fprintf('  Initial information check: %.3f [%s]\n',...
	dblInfoRes1,getTime);

%% compare on same data set
intResamplings = size(cellNeuronsX,1);

%% run through combinations
intStimTypes = 2;
intComparisons = 2;
matCompareTypes = [1 2; 2 1];

%% go through reversals
for boolReverseType=[0 1]
	%reverse
if boolReverseType==1
	strSubspace = 'SplitSubspace';
else
	strSubspace = 'UnsplitSubspace';
end
%% pre-allocate
intMaxRank = numel(vecRank);
intSamplePoints = intResamplings*intComparisons;
%prediction
matAggR2Full_CV = nan(intSamplePoints,intMaxRank);
matAggR2_CV =nan(intSamplePoints,intMaxRank);
matAggR2Rand_CV = nan(intSamplePoints,intMaxRank);

matAggFisherProj = nan(intSamplePoints,intMaxRank);
matAggFisherTotal = nan(intSamplePoints,intMaxRank);
matAggFisherProjRand = nan(intSamplePoints,intMaxRank);

matAggFisherProjT = nan(intSamplePoints,intMaxRank);
matAggFisherTotalT = nan(intSamplePoints,intMaxRank);
matAggFisherProjRandT = nan(intSamplePoints,intMaxRank);

matAggFisherProjSplit = nan(intSamplePoints,intMaxRank);
matAggFisherTotalSplit = nan(intSamplePoints,intMaxRank);
matAggFisherProjRandSplit = nan(intSamplePoints,intMaxRank);

matAggFisherProjTSplit = nan(intSamplePoints,intMaxRank);
matAggFisherTotalTSplit = nan(intSamplePoints,intMaxRank);
matAggFisherProjRandTSplit = nan(intSamplePoints,intMaxRank);

%% run through combinations
for intComparisonIdx=1:intComparisons
	for intResampling=1:intResamplings
		%% start
		intSaveIdx = intResampling+intResamplings*(intComparisonIdx-1);
		fprintf('    Running comp %d/%d, sample %d/%d {%d} [%s]\n',...
			intComparisonIdx,intComparisons,intResampling,intResamplings,intSaveIdx,getTime);
		
		%% get data
		vecNeuronsX = cellNeuronsX{intResampling,1};
		intNeuronsX = intSizeX;
		vecNeuronsY = cellNeuronsY{intResampling,1};
		intNeuronsY = intSizeY;
		
		%% define stim 2
		intStimType = matCompareTypes(intComparisonIdx,1);
		intOtherStim = matCompareTypes(intComparisonIdx,2);
		
		%% select data
		matX1 = cellMatX{intResampling,intStimType};
		matY1 = cellMatY{intResampling,intStimType};
		%if (any(range(matX1,1)==0) || any(range(matY1,1)==0)),continue;end
		vecMeanX1 = mean(matX1);
		vecMeanY1 = mean(matY1);
		matX1N = bsxfun(@minus,matX1,vecMeanX1);
		matY1N = bsxfun(@minus,matY1,vecMeanY1);
		
		%% select other data
		matX2 = cellMatX{intResampling,intOtherStim};
		matY2 = cellMatY{intResampling,intOtherStim};
		matX2N = bsxfun(@minus,matX2,vecMeanX1);
		matY2N = bsxfun(@minus,matY2,vecMeanY1);
		%{
		%% select data
		matX1 = matRespUnsplit(vecTrialStimTypeUnsplit==intStimType,vecNeuronsX);
		matY1 = matRespUnsplit(vecTrialStimTypeUnsplit==intStimType,vecNeuronsY);
		%if (any(range(matX1,1)==0) || any(range(matY1,1)==0)),continue;end
		vecMeanX1 = mean(matX1);
		vecMeanY1 = mean(matY1);
		matX1N = bsxfun(@minus,matX1,vecMeanX1);
		matY1N = bsxfun(@minus,matY1,vecMeanY1);
		
		%% select other data
		matX2 = matRespUnsplit(vecTrialStimTypeUnsplit==intOtherStim,vecNeuronsX);
		matY2 = matRespUnsplit(vecTrialStimTypeUnsplit==intOtherStim,vecNeuronsY);
		matX2N = bsxfun(@minus,matX2,vecMeanX1);
		matY2N = bsxfun(@minus,matY2,vecMeanY1);
		%}
		%get F' and covariance
		vecFprime = ((mean(matX2N,1) - mean(matX1N,1))')/dblDiffTheta;
		matCov = (cov(matX1N) + cov(matX2N))/2;
		%matCov = cov(matX1N);
		
		%get F' and covariance for target
		vecFprimeT = ((mean(matY2N,1) - mean(matY1N,1))')/dblDiffTheta;
		matCovT = (cov(matY1N) + cov(matY2N))/2;
		%matCovT = cov(matY1N);
		
		%% get f' for split data
		vecTrials1 = flat(matSplitTrialIdx(cellTrials{intResampling,intStimType},:));
		matX1Split = matRespSplit(vecTrials1,vecNeuronsX);
		matY1Split = matRespSplit(vecTrials1,vecNeuronsY);
		%if (any(range(matX1,1)==0) || any(range(matY1,1)==0)),continue;end
		vecMeanX1Split = mean(matX1Split);
		vecMeanY1Split = mean(matY1Split);
		matX1NSplit = bsxfun(@minus,matX1Split,vecMeanX1Split);
		matY1NSplit = bsxfun(@minus,matY1Split,vecMeanY1Split);
		
		%% select other data
		vecTrials2 = flat(matSplitTrialIdx(cellTrials{intResampling,intOtherStim},:));
		matX2Split = matRespSplit(vecTrials2,vecNeuronsX);
		matY2Split = matRespSplit(vecTrials2,vecNeuronsY);
		matX2NSplit = bsxfun(@minus,matX2Split,vecMeanX1Split);
		matY2NSplit = bsxfun(@minus,matY2Split,vecMeanY1Split);
		
		%get F' and covariance
		vecFprimeSplit = ((mean(matX2NSplit,1) - mean(matX1NSplit,1))')/dblDiffTheta;
		matCovSplit = (cov(matX1NSplit) + cov(matX2NSplit))/2;
		%matCovSplit = cov(matX1NSplit);
		
		%get F' and covariance for target
		vecFprimeTSplit = ((mean(matY2NSplit,1) - mean(matY1NSplit,1))')/dblDiffTheta;
		matCovTSplit = (cov(matY1NSplit) + cov(matY2NSplit))/2;
		%matCovTSplit = cov(matY1NSplit);
		
		%% total information
		matComb = [matX1N;matX2N];
		vecCombs = [ones(size(matX1N,1),1);2*ones(size(matX2N,1),1)];
		
		matComb2 = [matX1NSplit;matX2NSplit];
		vecCombs2 = [ones(size(matX1NSplit,1),1);2*ones(size(matX2NSplit,1),1)];
		
		[~,~,dblInfoUnsplit]=getSeparation(matComb,vecCombs,false,dblDiffTheta);
		[~,~,dblInfoSplit]=getSeparation(matComb2,vecCombs2,false,dblDiffTheta);
		
		fprintf('    Information unsplit: %.3f; split: %.3f [%s]\n',...
			dblInfoUnsplit,dblInfoSplit,getTime);
		
		%% reverse split/unsplit?
		if boolReverseType
			[matX1N,matX1NSplit] = swap(matX1NSplit,matX1N);
			[matY1N,matY1NSplit] = swap(matY1NSplit,matY1N);
			
			[vecFprime,vecFprimeSplit] = swap(vecFprimeSplit,vecFprime);
			[matCov,matCovSplit] = swap(matCovSplit,matCov);
			[vecFprimeT,vecFprimeTSplit] = swap(vecFprimeTSplit,vecFprimeT);
			[matCovT,matCovTSplit] = swap(matCovTSplit,matCovT);
		end
		
		%% find optimal regularisation parameter for ridge regression
		gMatX = matX1N;
		gMatY = matY1N;
		gIntFolds = intFoldK;
		
		%OLS
		[vecR2_OLS_CV, sOLS, intOutputFlag] = doRidge_CV(matX1N, matY1N, intFoldK, 0);
		matB_ols_sub = sOLS.matB_Full;
		
		if boolDoRegularization
			%find optimal lambda
			dblInitL = 10;
			[dblLambda,dblR2Optim] = fminsearch(@getRidgeWrapper,dblInitL);
			[vecR2_Ridge_CV, sRidge] = doRidge_CV(matX1N, matY1N, intFoldK, dblLambda);
			matB_ridge_sub = sRidge.matB_Full;
			
		else
			%set L to 0
			dblLambda=0;
			dblR2Optim = nan;
			matB_ridge_sub = matB_ols_sub;
			
		end
		
		
		%% figure 4; full prediction versus reduced rank regression, as function of dimensionality
		matR2Full_CV = vecR2_Ridge_CV;
		matR2_CV = nan(numel(vecRank),intFoldK);
		matR2_NonCV = nan(numel(vecRank),intFoldK);
		matR2_Rem = nan(numel(vecRank),intFoldK);
		matR2Rand_CV = nan(numel(vecRank),intRandIters,intFoldK);
		matR2Rand_NonCV = nan(numel(vecRank),intRandIters,intFoldK);
		matR2RandRem = nan(numel(vecRank),intRandIters,intFoldK);
		
		matFisherProj = nan(numel(vecRank),intFoldK);
		matFisherOrth = nan(numel(vecRank),intFoldK);
		matFisherTotal = nan(numel(vecRank),intFoldK);
		matFisherProjRand = nan(numel(vecRank),intRandIters,intFoldK);
		matFisherOrthRand = nan(numel(vecRank),intRandIters,intFoldK);
		matFisherTotalRand = nan(numel(vecRank),intRandIters,intFoldK);
		
		matFisherProjT = nan(numel(vecRank),intFoldK);
		matFisherOrthT = nan(numel(vecRank),intFoldK);
		matFisherTotalT = nan(numel(vecRank),intFoldK);
		matFisherProjRandT = nan(numel(vecRank),intRandIters,intFoldK);
		matFisherOrthRandT = nan(numel(vecRank),intRandIters,intFoldK);
		matFisherTotalRandT = nan(numel(vecRank),intRandIters,intFoldK);
		
		matFisherProjSplit = nan(numel(vecRank),intFoldK);
		matFisherOrthSplit = nan(numel(vecRank),intFoldK);
		matFisherTotalSplit = nan(numel(vecRank),intFoldK);
		matFisherProjRandSplit = nan(numel(vecRank),intRandIters,intFoldK);
		matFisherOrthRandSplit = nan(numel(vecRank),intRandIters,intFoldK);
		matFisherTotalRandSplit = nan(numel(vecRank),intRandIters,intFoldK);
		
		matFisherProjTSplit = nan(numel(vecRank),intFoldK);
		matFisherOrthTSplit = nan(numel(vecRank),intFoldK);
		matFisherTotalTSplit = nan(numel(vecRank),intFoldK);
		matFisherProjRandTSplit = nan(numel(vecRank),intRandIters,intFoldK);
		matFisherOrthRandTSplit = nan(numel(vecRank),intRandIters,intFoldK);
		matFisherTotalRandTSplit = nan(numel(vecRank),intRandIters,intFoldK);
		
	
		vecSubFacProj = nan(1,numel(vecRank));
		vecProdFacProj = nan(1,numel(vecRank));
		vecSubFacOrth = nan(1,numel(vecRank));
		vecProdFacOrth = nan(1,numel(vecRank));
		%%
		intC = 0;
		for intRankInT = vecRank
			intC = intC + 1;
			%% get bias correction factors for Fisher information
			intTrials12 = size(matX1N,1);
			%fprintf('     Running rank %d/%d, #trials=%d [%s]\n',intC,numel(vecRank),intTrials12,getTime);
					
			%proj
			intGroupSizeProj = intRankInT;
			dblSubFacProj =(2*intGroupSizeProj)/(intTrials12*(dblDiffTheta.^2));
			dblProdFacRawProj = ((2*intTrials12-intGroupSizeProj-3)/(2*intTrials12-2));
			vecSubFacProj(intC) = dblSubFacProj;
			vecProdFacProj(intC) = dblProdFacRawProj;
			
			%orth
			intGroupSizeOrth = intNeuronsX-intRankInT;
			dblSubFacOrth =(2*intGroupSizeOrth)/(intTrials12*(dblDiffTheta.^2));
			dblProdFacRawOrth = ((2*intTrials12-intGroupSizeOrth-3)/(2*intTrials12-2));
			vecSubFacOrth(intC) = dblSubFacOrth;
			vecProdFacOrth(intC) = dblProdFacRawOrth;
			
			%orth target
			intGroupSizeOrthT = intNeuronsY-intRankInT;
			dblSubFacOrthT =(2*intGroupSizeOrthT)/(intTrials12*(dblDiffTheta.^2));
			dblProdFacRawOrthT = ((2*intTrials12-intGroupSizeOrthT-3)/(2*intTrials12-2));
			%vecSubFacOrthT(intC) = dblSubFacOrthT;
			%vecProdFacOrthT(intC) = dblProdFacRawOrthT;
			
			%total
			intGroupSizeTotal = intNeuronsX;
			dblSubFacTotal =(2*intGroupSizeTotal)/(intTrials12*(dblDiffTheta.^2));
			dblProdFacRawTotal = ((2*intTrials12-intGroupSizeTotal-3)/(2*intTrials12-2));
			
			%total target
			intGroupSizeTotalT = intNeuronsY;
			dblSubFacTotalT =(2*intGroupSizeTotalT)/(intTrials12*(dblDiffTheta.^2));
			dblProdFacRawTotalT = ((2*intTrials12-intGroupSizeTotalT-3)/(2*intTrials12-2));
			
			%% perform noise prediction using reduced rank regression
			%[matW, dblMSE, intRankOutT, sSuppOut] = doRdRankReg(matX1N, matY1N, 'rank', intRankInT);
			[matW, dblMSE, intRankOutT, sRRR] = doRRR_CV(matX1N, matY1N, intRankInT, intFoldK, dblLambda, matB_ridge_sub, matCov);
			matW_Basis = orth(matW);
			
			matR2_CV(intC,:) = sRRR.vecR2_CV;
			matR2_Rem(intC,:) = sRRR.vecR2_Rem;
			matR2_NonCV(intC,:) = sRRR.vecR2_NonCV;
			for intRandIter=1:intRandIters
				%% perform noise prediction with random subspaces
				matW_RandPred = normrnd(0,1,[intNeuronsX,intRankInT]);
				dblRankB = rank(matW_Basis);
				dblRankW = rank(matW_RandPred);
				if dblRankB ~= dblRankW
					warning([mfilename ':RankInconsistent'],'Rank of projection basis and random matrix do not match');
					fprintf('Ranks; Proj: %d, Rand: %d; Type %d, Resample %d, Rank %d [%s]\n',...
						dblRankB,dblRankW,intStimType,intResampling,intRankInT,getTime);
				end
				
				%calculate reduced matrix
				matReducedX1 = matX1N*matW_RandPred;
				
				%ridge reg
				[vecR2_CV, sRand] = doRidge_CV(matReducedX1, matY1N, intFoldK, dblLambda);
				matR2Rand_CV(intC,intRandIter,:) = vecR2_CV;
				matR2Rand_NonCV(intC,intRandIter,:) = sRand.vecR2_NonCV;
				cellSaveRandB{intRandIter} = sRand.cellB;
				cellSaveRandW{intRandIter} = matW_RandPred;
				
				%% perform prediction when removing random dimensions
				matR2RandRem(intC,intRandIter,:) = 0;
			end
			
			%% go through folds to get information in different subspaces
			for intFold=1:intFoldK
				%% get subspace
				matW = sRRR.matW_Train(:,:,intFold);
				matW_Basis = orth(matW); %basis in source
				matW_BasisT = orth(matW'); %basis in target
						
				%% source pop
				%get Fisher information in reduced space
				[dblFisherProj, dblFisherOrth, dblFisherTotal] = getFisherInSubspace(vecFprime, matCov, matW_Basis);
				matFisherProj(intC,intFold) = dblFisherProj*dblProdFacRawProj -dblSubFacProj;
				matFisherOrth(intC,intFold) = dblFisherOrth*dblProdFacRawOrth-dblSubFacOrth;
				matFisherTotal(intC,intFold) = dblFisherTotal*dblProdFacRawTotal-dblSubFacTotal;
						
				%get Fisher information selecting random dimensions
				for intRandIter=1:intRandIters
					matB = cellSaveRandB{intRandIter}{intFold};
					matW_Rand = orth(cellSaveRandW{intRandIter} * matB);
					
					dblRankB = rank(matW_Basis);
					dblRankW = rank(matW_Rand);
					if dblRankB ~= dblRankW
						warning([mfilename ':RankInconsistent'],'Rank of projection basis and random matrix do not match');
						fprintf('Ranks; Proj: %d, Rand: %d; Type %d, Resample %d, Rank %d [%s]\n',...
							dblRankB,dblRankW,intStimType,intResampling,intRankInT,getTime);
					end
					
					[dblFisherProjRand, dblFisherOrthRand, dblFisherTotalRand] = getFisherInSubspace(vecFprime, matCov, matW_Rand);
					matFisherProjRand(intC,intRandIter,intFold) = dblFisherProjRand*dblProdFacRawProj -dblSubFacProj;
					matFisherOrthRand(intC,intRandIter,intFold) = dblFisherOrthRand*dblProdFacRawOrth-dblSubFacOrth;
					matFisherTotalRand(intC,intRandIter,intFold) = dblFisherTotalRand*dblProdFacRawTotal-dblSubFacTotal;
				end
				
				%% target pop
				%get Fisher information in reduced space
				[dblFisherProjT, dblFisherOrthT, dblFisherTotalT] = getFisherInSubspace(vecFprimeT, matCovT, matW_BasisT);
				matFisherProjT(intC,intFold) = dblFisherProjT*dblProdFacRawProj -dblSubFacProj;
				matFisherOrthT(intC,intFold) = dblFisherOrthT*dblProdFacRawOrthT-dblSubFacOrthT;
				matFisherTotalT(intC,intFold) = dblFisherTotalT*dblProdFacRawTotalT-dblSubFacTotalT;
				
				%get Fisher information selecting random dimensions
				for intRandIter=1:intRandIters
					matB = cellSaveRandB{intRandIter}{intFold};
					matW_RandT = orth((cellSaveRandW{intRandIter} * matB)');
					
					dblRankB = rank(matW_BasisT);
					dblRankW = rank(matW_RandT);
					if dblRankB ~= dblRankW
						warning([mfilename ':RankInconsistent'],'Rank of projection basis and random matrix do not match');
						fprintf('Ranks; Proj: %d, Rand: %d; Type %d, Resample %d, Rank %d [%s]\n',...
							dblRankB,dblRankW,intStimType,intResampling,intRankInT,getTime);
					end
					[dblFisherProjRandT, dblFisherOrthRandT, dblFisherTotalRandT] = getFisherInSubspace(vecFprimeT, matCovT, matW_RandT);
					matFisherProjRandT(intC,intRandIter,intFold) = dblFisherProjRandT*dblProdFacRawProj -dblSubFacProj;
					matFisherOrthRandT(intC,intRandIter,intFold) = dblFisherOrthRandT*dblProdFacRawOrthT-dblSubFacOrthT;
					matFisherTotalRandT(intC,intRandIter,intFold) = dblFisherTotalRandT*dblProdFacRawTotalT-dblSubFacTotalT;
				end
			end
			
			%% go through folds to get information in different subspaces for split data
			for intFold=1:intFoldK
				intTrials12 = size(matX1NSplit,1);
				
				%proj
				intGroupSizeProj = intRankInT;
				dblSubFacProj =(2*intGroupSizeProj)/(intTrials12*(dblDiffTheta.^2));
				dblProdFacRawProj = ((2*intTrials12-intGroupSizeProj-3)/(2*intTrials12-2));
				
				%orth
				intGroupSizeOrth = intNeuronsX-intRankInT;
				dblSubFacOrth =(2*intGroupSizeOrth)/(intTrials12*(dblDiffTheta.^2));
				dblProdFacRawOrth = ((2*intTrials12-intGroupSizeOrth-3)/(2*intTrials12-2));
				
				%orth target
				intGroupSizeOrthT = intNeuronsY-intRankInT;
				dblSubFacOrthT =(2*intGroupSizeOrthT)/(intTrials12*(dblDiffTheta.^2));
				dblProdFacRawOrthT = ((2*intTrials12-intGroupSizeOrthT-3)/(2*intTrials12-2));
				
				%total
				intGroupSizeTotal = intNeuronsX;
				dblSubFacTotal =(2*intGroupSizeTotal)/(intTrials12*(dblDiffTheta.^2));
				dblProdFacRawTotal = ((2*intTrials12-intGroupSizeTotal-3)/(2*intTrials12-2));
				
				%total target
				intGroupSizeTotalT = intNeuronsY;
				dblSubFacTotalT =(2*intGroupSizeTotalT)/(intTrials12*(dblDiffTheta.^2));
				dblProdFacRawTotalT = ((2*intTrials12-intGroupSizeTotalT-3)/(2*intTrials12-2));
				
				%% get subspace
				matW = sRRR.matW_Train(:,:,intFold);
				matW_Basis = orth(matW); %basis in source
				matW_BasisT = orth(matW'); %basis in target
				
				%% source pop
				%get Fisher information in reduced space
				[dblFisherProj, dblFisherOrth, dblFisherTotal] = getFisherInSubspace(vecFprimeSplit, matCovSplit, matW_Basis);
				matFisherProjSplit(intC,intFold) = dblFisherProj*dblProdFacRawProj -dblSubFacProj;
				matFisherOrthSplit(intC,intFold) = dblFisherOrth*dblProdFacRawOrth-dblSubFacOrth;
				matFisherTotalSplit(intC,intFold) = dblFisherTotal*dblProdFacRawTotal-dblSubFacTotal;
				
				%get Fisher information selecting random dimensions
				for intRandIter=1:intRandIters
					matB = cellSaveRandB{intRandIter}{intFold};
					matW_Rand = orth(cellSaveRandW{intRandIter} * matB);
					
					dblRankB = rank(matW_Basis);
					dblRankW = rank(matW_Rand);
					if dblRankB ~= dblRankW
						warning([mfilename ':RankInconsistent'],'Rank of projection basis and random matrix do not match');
						fprintf('Ranks; Proj: %d, Rand: %d; Type %d, Resample %d, Rank %d [%s]\n',...
							dblRankB,dblRankW,intStimType,intResampling,intRankInT,getTime);
					end
					
					[dblFisherProjRand, dblFisherOrthRand, dblFisherTotalRand] = getFisherInSubspace(vecFprimeSplit, matCovSplit, matW_Rand);
					matFisherProjRandSplit(intC,intRandIter,intFold) = dblFisherProjRand*dblProdFacRawProj -dblSubFacProj;
					matFisherOrthRandSplit(intC,intRandIter,intFold) = dblFisherOrthRand*dblProdFacRawOrth-dblSubFacOrth;
					matFisherTotalRandSplit(intC,intRandIter,intFold) = dblFisherTotalRand*dblProdFacRawTotal-dblSubFacTotal;
				end
				
				%% target pop
				%get Fisher information in reduced space
				[dblFisherProjT, dblFisherOrthT, dblFisherTotalT] = getFisherInSubspace(vecFprimeTSplit, matCovTSplit, matW_BasisT);
				matFisherProjTSplit(intC,intFold) = dblFisherProjT*dblProdFacRawProj -dblSubFacProj;
				matFisherOrthTSplit(intC,intFold) = dblFisherOrthT*dblProdFacRawOrthT-dblSubFacOrthT;
				matFisherTotalTSplit(intC,intFold) = dblFisherTotalT*dblProdFacRawTotalT-dblSubFacTotalT;
				
				%get Fisher information selecting random dimensions
				for intRandIter=1:intRandIters
					matB = cellSaveRandB{intRandIter}{intFold};
					matW_RandT = orth((cellSaveRandW{intRandIter} * matB)');
					
					dblRankB = rank(matW_BasisT);
					dblRankW = rank(matW_RandT);
					if dblRankB ~= dblRankW
						warning([mfilename ':RankInconsistent'],'Rank of projection basis and random matrix do not match');
						fprintf('Ranks; Proj: %d, Rand: %d; Type %d, Resample %d, Rank %d [%s]\n',...
							dblRankB,dblRankW,intStimType,intResampling,intRankInT,getTime);
					end
					[dblFisherProjRandT, dblFisherOrthRandT, dblFisherTotalRandT] = getFisherInSubspace(vecFprimeTSplit, matCovTSplit, matW_RandT);
					matFisherProjRandTSplit(intC,intRandIter,intFold) = dblFisherProjRandT*dblProdFacRawProj -dblSubFacProj;
					matFisherOrthRandTSplit(intC,intRandIter,intFold) = dblFisherOrthRandT*dblProdFacRawOrthT-dblSubFacOrthT;
					matFisherTotalRandTSplit(intC,intRandIter,intFold) = dblFisherTotalRandT*dblProdFacRawTotalT-dblSubFacTotalT;
				end
			end
		end
		
		%% put data in aggregates
		if any(isnan(matR2Full_CV)),error;end
		matAggR2Full_CV(intSaveIdx,:) = mean(matR2Full_CV(:));
		matAggR2_CV(intSaveIdx,:) = xmean(matR2_CV,2);
		matAggR2Rand_CV(intSaveIdx,:) = xmean(xmean(matR2Rand_CV,2),3);
		
		matAggFisherProj(intSaveIdx,:) = xmean(matFisherProj,2);
		matAggFisherTotal(intSaveIdx,:) = xmean(matFisherTotal,2);
		matAggFisherProjRand(intSaveIdx,:) = xmean(xmean(matFisherProjRand,2),3);
		
		matAggFisherProjT(intSaveIdx,:) = xmean(matFisherProjT,2);
		matAggFisherTotalT(intSaveIdx,:) = xmean(matFisherTotalT,2);
		matAggFisherProjRandT(intSaveIdx,:) = xmean(xmean(matFisherProjRandT,2),3);
		
		matAggFisherProjSplit(intSaveIdx,:) = xmean(matFisherProjSplit,2);
		matAggFisherTotalSplit(intSaveIdx,:) = xmean(matFisherTotalSplit,2);
		matAggFisherProjRandSplit(intSaveIdx,:) = xmean(xmean(matFisherProjRandSplit,2),3);
		
		matAggFisherProjTSplit(intSaveIdx,:) = xmean(matFisherProjTSplit,2);
		matAggFisherTotalTSplit(intSaveIdx,:) = xmean(matFisherTotalTSplit,2);
		matAggFisherProjRandTSplit(intSaveIdx,:) = xmean(xmean(matFisherProjRandTSplit,2),3);
		
		
	end
end
%%
figure
%full screen
drawnow;
jFig = get(handle(gcf), 'JavaFrame');
jFig.setMaximized(true);
figure(gcf);
drawnow;
intTrialsUnsplit = size(matX1N,1);
subplot(2,3,1)
hold on
errorbar(vecRank,xmean(matAggR2Full_CV,1),xstd(matAggR2Full_CV,1)/sqrt(intSamplePoints),'k-');
errorbar(vecRank,xmean(matAggR2_CV,1),xstd(matAggR2_CV,1)/sqrt(intSamplePoints),'b-');
errorbar(vecRank,xmean(matAggR2Rand_CV,1),xstd(matAggR2Rand_CV,1)/sqrt(intSamplePoints),'r-');

hold off
title(sprintf('Noise prediction in %s [%s]',strPredArea,strSubspace),'Interpreter','none');
ylabel('Noise predictability (R^2)');
xlabel('Dimensionality of subspace');
xlim([min(vecRank)-1 max(vecRank+1)]);
if max(get(gca,'ylim')) < 0, dblMinY = min(get(gca,'ylim'));else dblMinY=0;end
ylim([dblMinY max(get(gca,'ylim'))]);
fixfig;
h=legend({'Full model CV','Predictive subsp CV','Random subsp CV','Pred subsp train','Random Train','Orth to pred subsp'});
set(h,'Location','Best');

subplot(2,3,2)
%matAggFullMinusPred = matAggFisherDimDepTotalAvg - matAggFisherDimDepProjAvg;
hold on
errorbar(vecRank,xmean(matAggFisherTotal,1),xstd(matAggFisherTotal,1)/sqrt(intSamplePoints),'k-');
errorbar(vecRank,xmean(matAggFisherProj,1),xstd(matAggFisherProj,1)/sqrt(intSamplePoints),'b-');
errorbar(vecRank,xmean(matAggFisherProjRand,1),xstd(matAggFisherProjRand,1)/sqrt(intSamplePoints),'r-');

hold off
%h2=legend({'Total','Orth to Rand Subsp','Orth to Pred Subsp','Pred Subspace','Rand Subspace','(Full) - (Pred Subsp)'});
%h2=legend({'Total','Orth to Rand Subsp','Orth to Pred Subsp','Pred Subspace','Rand Subspace'});
%set(h2,'Location','west');
title(sprintf('Info in source subspace [N=%d,T=%d] (+/-SEM)',intSizeX,intTrialsUnsplit));
ylabel('Fisher information (d''^2/d{\theta}^2)');
xlabel('Dimensionality of subspace');
xlim([min(vecRank)-1 max(vecRank+1)]);
if max(get(gca,'ylim')) < 0, dblMinY = min(get(gca,'ylim'));else dblMinY=0;end
ylim([dblMinY max(get(gca,'ylim'))]);
fixfig;


subplot(2,3,3)
hold on
errorbar(vecRank,xmean(matAggFisherTotalT,1),xstd(matAggFisherTotalT,1)/sqrt(intSamplePoints),'k-');
errorbar(vecRank,xmean(matAggFisherProjT,1),xstd(matAggFisherProjT,1)/sqrt(intSamplePoints),'b-');
errorbar(vecRank,xmean(matAggFisherProjRandT,1),xstd(matAggFisherProjRandT,1)/sqrt(intSamplePoints),'r-');

hold off
%h2=legend({'Total','Orth to Rand Subsp','Orth to Pred Subsp','Pred Subspace','Rand Subspace','(Full) - (Pred Subsp)'});
%h2=legend({'Total','Orth to Rand Subsp','Orth to Pred Subsp','Pred Subspace','Rand Subspace'});
%set(h2,'Location','west');
title(sprintf('Info in target subspace [N=%d,T=%d] (+/-SEM)',intSizeY,intTrialsUnsplit));
ylabel('Fisher information (d''^2/d{\theta}^2)');
xlabel('Dimensionality of subspace');
xlim([min(vecRank)-1 max(vecRank+1)]);
if max(get(gca,'ylim')) < 0, dblMinY = min(get(gca,'ylim'));else dblMinY=0;end
ylim([dblMinY max(get(gca,'ylim'))]);
fixfig;

subplot(2,3,4)
title(strTag,'Interpreter','none');
axis off;
fixfig
intTrialsSplit = size(matX1NSplit,1);
subplot(2,3,5)
%matAggFullMinusPred = matAggFisherDimDepTotalAvg - matAggFisherDimDepProjAvg;
hold on
errorbar(vecRank,xmean(matAggFisherTotalSplit,1),xstd(matAggFisherTotalSplit,1)/sqrt(intSamplePoints),'k-');
errorbar(vecRank,xmean(matAggFisherProjSplit,1),xstd(matAggFisherProjSplit,1)/sqrt(intSamplePoints),'b-');
errorbar(vecRank,xmean(matAggFisherProjRandSplit,1),xstd(matAggFisherProjRandSplit,1)/sqrt(intSamplePoints),'r-');

hold off
%h2=legend({'Total','Orth to Rand Subsp','Orth to Pred Subsp','Pred Subspace','Rand Subspace','(Full) - (Pred Subsp)'});
%h2=legend({'Total','Orth to Rand Subsp','Orth to Pred Subsp','Pred Subspace','Rand Subspace'});
%set(h2,'Location','west');
title(sprintf('Info in source subspace [N=%d,T=%d] (+/-SEM)',intSizeX,intTrialsSplit));
ylabel('Fisher information (d''^2/d{\theta}^2)');
xlabel('Dimensionality of subspace');
xlim([min(vecRank)-1 max(vecRank+1)]);
if max(get(gca,'ylim')) < 0, dblMinY = min(get(gca,'ylim'));else dblMinY=0;end
ylim([dblMinY max(get(gca,'ylim'))]);
fixfig;


subplot(2,3,6)
hold on
errorbar(vecRank,xmean(matAggFisherTotalTSplit,1),xstd(matAggFisherTotalTSplit,1)/sqrt(intSamplePoints),'k-');
errorbar(vecRank,xmean(matAggFisherProjTSplit,1),xstd(matAggFisherProjTSplit,1)/sqrt(intSamplePoints),'b-');
errorbar(vecRank,xmean(matAggFisherProjRandTSplit,1),xstd(matAggFisherProjRandTSplit,1)/sqrt(intSamplePoints),'r-');

hold off
%h2=legend({'Total','Orth to Rand Subsp','Orth to Pred Subsp','Pred Subspace','Rand Subspace','(Full) - (Pred Subsp)'});
%h2=legend({'Total','Orth to Rand Subsp','Orth to Pred Subsp','Pred Subspace','Rand Subspace'});
%set(h2,'Location','west');
title(sprintf('Info in target subspace [N=%d,T=%d] (+/-SEM)',intSizeY,intTrialsSplit));
ylabel('Fisher information (d''^2/d{\theta}^2)');
xlabel('Dimensionality of subspace');
xlim([min(vecRank)-1 max(vecRank+1)]);
if max(get(gca,'ylim')) < 0, dblMinY = min(get(gca,'ylim'));else dblMinY=0;end
ylim([dblMinY max(get(gca,'ylim'))]);
fixfig;

%save figure
if boolSaveFigs
	%figure(hFigDD);
	drawnow;
	jFig = get(handle(gcf), 'JavaFrame');
	jFig.setMaximized(true);
	figure(gcf);
	drawnow;
	export_fig([strFigDir  'Block5_' strSubspace '_Area' num2str(intWithinArea) '_T' num2str(intUseIndepTrials) '_' strTag '.tif']);
	export_fig([strFigDir  'Block5_' strSubspace '_Area' num2str(intWithinArea) '_T' num2str(intUseIndepTrials) '_' strTag '.pdf']);
end
end
