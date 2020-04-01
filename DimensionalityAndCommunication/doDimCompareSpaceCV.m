function sOut = doDimCompareSpaceCV(cellMatX1,cellMatY1,cellMatX2,cellMatY2,intMaxDim,dblLambda,boolShuffle)
	%% prep
	if ~exist('dblLambda','var') || isempty(dblLambda)
		dblLambda = 0;
	end
	if ~exist('boolShuffle','var') || isempty(boolShuffle)
		boolShuffle = false;
	end
	
	%% pre-allocate
	intNeurons = size(cellMatX1{1},2);
	intStimTypes = size(cellMatX1,2);
	intResamplings = size(cellMatX1,1);
	%get dom dim vec
	%baselines
	matR2_xx2x = nan(intResamplings,intStimTypes,intMaxDim);
	matR2_yy2y = nan(intResamplings,intStimTypes,intMaxDim);
	matR2_xRRR2y = nan(intResamplings,intStimTypes,intMaxDim);
	matR2_yRRR2x = nan(intResamplings,intStimTypes,intMaxDim);
	%tests, x-y
	matR2_xx2y = nan(intResamplings,intStimTypes,intMaxDim);
	matR2_xy2y = nan(intResamplings,intStimTypes,intMaxDim);
	%controls, y-x
	matR2_yy2x = nan(intResamplings,intStimTypes,intMaxDim);
	matR2_yx2x = nan(intResamplings,intStimTypes,intMaxDim);
	%% run
	for intStimType=1:intStimTypes
		for intResampling=1:intResamplings
			%% run analyses
			% select data
			matX1 = zscore(cellMatX1{intResampling,intStimType},[],1);
			matY1 = zscore(cellMatY1{intResampling,intStimType},[],1);
			matX2 = zscore(cellMatX2{intResampling,intStimType},[],1);
			matY2 = zscore(cellMatY2{intResampling,intStimType},[],1);
			if (any(range(matX1,1)==0) || any(range(matX2,1)==0) || any(range(matY1,1)==0) || any(range(matY2,1)==0)),continue;end
			%shuffle
			if boolShuffle
				matX2 = matX2(randperm(size(matX2,1)),:);
				matY2 = matY2(randperm(size(matX2,1)),:);
			end
			
			% get principal components  for X
			%[matPCs2,score,vecLambdas2,tsquared,explained,mu] = pca(matX); %slow, but gives same answer
			[matPCsX1,vecLambdas]=eig(cov(matX1),'vector');
			[vecLambdas,vecReorder] = sort(vecLambdas,'descend');
			matLambdas = diag(vecLambdas);
			matPCsX1 = matPCsX1(:,vecReorder);
			
			% X(X) -> X (interne prediction B1: PCA/FA)
			vecR2_xx2x = getRegInSpaceCV(matX1,matPCsX1,matX1,matX2,matX2,dblLambda);
			matR2_xx2x(intResampling,intStimType,:) = vecR2_xx2x;
			
			% X(X) -> Y (hoeveel overlapt B2 met interne variance in B1?)
			vecR2_xx2y = getRegInSpaceCV(matX1,matPCsX1,matY1,matX2,matY2,dblLambda);
			matR2_xx2y(intResampling,intStimType,:) = vecR2_xx2y;
			
			% X(RRR) -> Y  (greedy regression: wat is de grootte van de B1->B2 subspace?
			%[vergelijk met (2): zijn predictive en internal dimensions hetzelfde?])
			[vecR2_xRRR2y,vecR2_NonCV_xRRR2y,cellB_xRRR2y] = getRedRankRegCV(matX1,matY1,matX2,matY2,dblLambda);
			matR2_xRRR2y(intResampling,intStimType,:) = vecR2_xRRR2y;
			
			% get principal components  for Y
			%[matPCs2,score,vecLambdas2,tsquared,explained,mu] = pca(matX); %slow, but gives same answer
			[matPCsY1,vecLambdas]=eig(cov(matY1),'vector');
			[vecLambdas,vecReorder] = sort(vecLambdas,'descend');
			matLambdas = diag(vecLambdas);
			matPCsY1 = matPCsY1(:,vecReorder);
			
			% Y(Y) -> Y (interne prediction B2: PCA/FA)
			vecR2_yy2y = getRegInSpaceCV(matY1,matPCsY1,matY1,matY2,matY2,dblLambda);
			matR2_yy2y(intResampling,intStimType,:) = vecR2_yy2y;
			
			% X(Y) -> Y  (inverse regression: zijn interne dimensies van B2 hetzelfde als B1
			%[vergelijk met (2) en (4)]?)
			%Als je dezelfde R^2 curve krijgt uit B1(B1-PCs)->B2 als uit
			%B1(B2-PCs)->B2, dan weet je zeker dat B1 en B2 in identieke
			%subspaces liggen, met als enige caveat dat 1 van de 2 groter
			%kan zijn, maar dat in dat geval iig alle hogere dimensies een
			%kleinere PC moeten hebben dan de dimensies die beide subspaces
			%delen. Als de curves niet identiek zijn, dan weet je dat de
			%spaces niet hetzelfde zijn, en kun je afleiden uit de
			%groeisnelheid welke space groter is.
			vecR2_xy2y = getRegInSpaceCV(matX1,matPCsY1,matY1,matX2,matY2,dblLambda);
			matR2_xy2y(intResampling,intStimType,:) = vecR2_xy2y;
			
			%% controls
			%Y(Y) -> X
			vecR2_yy2x = getRegInSpaceCV(matY1,matPCsY1,matX1,matY2,matX2,dblLambda);
			matR2_yy2x(intResampling,intStimType,:) = vecR2_yy2x;
			
			%Y(RRR) -> X
			[vecR2_yRRR2x,vecR2_NonCV_yRRR2x,cellB_yRRR2x] = getRedRankRegCV(matY1,matX1,matY2,matX2,dblLambda);
			matR2_yRRR2x(intResampling,intStimType,:) = vecR2_yRRR2x;
			
			%Y(X) -> X
			vecR2_yx2x = getRegInSpaceCV(matY1,matPCsX1,matX1,matY2,matX2,dblLambda);
			matR2_yx2x(intResampling,intStimType,:) = vecR2_yx2x;
		end
	end
	
	%% build output
	sOut = struct;
	sOut.matR2_xx2x = matR2_xx2x;
	sOut.matR2_yy2y = matR2_yy2y;
	sOut.matR2_xRRR2y = matR2_xRRR2y;
	sOut.matR2_yRRR2x = matR2_yRRR2x;
	%tests
	sOut.matR2_xx2y = matR2_xx2y;
	sOut.matR2_xy2y = matR2_xy2y;
	%controls
	sOut.matR2_yy2x = matR2_yy2x;
	sOut.matR2_yx2x = matR2_yx2x;
	