function [matR2_hh2h,matR2_hm2h] = doDimCompareSpaceCV_HitsMiss(matRespHits,matRespMiss,intResamplings,dblLambda,boolShuffle)
	%% prep
	if ~exist('dblLambda','var') || isempty(dblLambda)
		dblLambda = 0;
	end
	if ~exist('boolShuffle','var') || isempty(boolShuffle)
		boolShuffle = false;
	end
	%% analyze
	%B1_hit(B1_hit)->B1_hit and compare with B1_hit(B1_miss)->B1_hit
	%B1_hit(B1_hit)->B1_hit: getRegInSpaceCV(B1_hit_train,PC(B1_hit_train),B1_hit_train,B1_hit_test,B1_hit_test)
	%B1_hit(B1_miss)->B1_hit: getRegInSpaceCV(B1_hit_train,PC(B1_miss_train),B1_hit_train,B1_hit_test,B1_hit_test)
	
	%REPEAT FOR:
	%B2_hit(B2_hit)->B2_hit: getRegInSpaceCV(B2_hit_train,PC(B2_hit_train),B1_hit_train,B1_hit_test,B1_hit_test)
	%B2_hit(B2_miss)->B2_hit: getRegInSpaceCV(B2_hit_train,PC(B2_miss_train),B2_hit_train,B2_hit_test,B2_hit_test)

	%% split data
	%split trials for hits into train/test
	matRespHits = matRespHits';
	matRespMiss = matRespMiss';
	
	intNeurons = size(matRespHits,2);
	intTrialsHits = size(matRespHits,1);
	intTrialsMiss = size(matRespMiss,1);
	
	intSplitSizeH = floor(intTrialsHits/2);
	
	cellHitRespTrain = cell(1,intResamplings);
	cellHitRespTest = cell(1,intResamplings);
	for intResampling=1:intResamplings
		vecRandPermTrials = randperm(intTrialsHits);
		cellHitRespTrain{intResampling} = matRespHits(vecRandPermTrials((intSplitSizeH+1):end),:);
		cellHitRespTest{intResampling} = matRespHits(vecRandPermTrials(1:intSplitSizeH),:);
	end
		
	%% get principal components
	%for hits
	[matPCsHits,vecLambdas]=eig(cov(matRespHits),'vector');
	[vecLambdas,vecReorder] = sort(vecLambdas,'descend');
	matPCsHits = matPCsHits(:,vecReorder);
	
	%for misses
	[matPCsMiss,vecLambdas]=eig(cov(matRespMiss),'vector');
	[vecLambdas,vecReorder] = sort(vecLambdas,'descend');
	matPCsMiss = matPCsMiss(:,vecReorder);
		
	%% pre-allocate
	%get dom dim vec
	matR2_hh2h = nan(intResamplings,intNeurons);
	matR2_hm2h = nan(intResamplings,intNeurons);
	
	%% run
	for intResampling=1:intResamplings
		%% run analyses
		% select data
		matTrain = zscore(cellHitRespTrain{intResampling},[],1);
		matTest = zscore(cellHitRespTest{intResampling},[],1);
		
		%shuffle
		if boolShuffle
			matTest = matTest(randperm(size(matTest,1)),:);%is this ori correct?
		end
		
		% X(X) -> X (interne prediction B1: PCA/FA)
		vecR2_hh2h = getRegInSpaceCV(matTrain,matPCsHits,matTrain,matTest,matTest,dblLambda);
		matR2_hh2h(intResampling,:) = vecR2_hh2h;
		
		%Y(X) -> X
		vecR2_hm2h = getRegInSpaceCV(matTrain,matPCsMiss,matTrain,matTest,matTest,dblLambda);
		matR2_hm2h(intResampling,:) = vecR2_hm2h;
	end
	
end
	