%% initialize
boolLoad = false;
boolSaveFigs = false;
if ~exist('matModelResp','var')
	clearvars;%close all;
	intLoadSim = 5; %% SET WHICH SIM TO LOAD HERE
	boolLoad = true;
end
if intLoadSim == 1 && boolLoad
	strSimulation = 'SimpleLine2017-01-18';
elseif intLoadSim == 2 && boolLoad
	strSimulation = 'Line2017-01-17';
elseif intLoadSim == 3 && boolLoad
	strSimulation = 'SimpleSquareGrating2017-01-19';
elseif intLoadSim == 4 && boolLoad
	strSimulation = 'SquareGrating2017-01-18';
	%Cells: 1440
	%Cort Conn: 403200
	%LGN Conn: 32256*2
	%Trials: 500 / Types: 2 / Reps: 250
	%Oris: [42.5 47.5]
elseif intLoadSim == 5 && boolLoad
	strSimulation = 'SquareGrating2017-01-25';
	%Cells: 1440
	%Cort Conn: 403200
	%LGN Conn: 32256*2
	%Trials: 3000 / Types: 2 / Reps: 1500
	%Oris: [42.5 47.5]
elseif intLoadSim == 6 && boolLoad
	strSimulation = 'SquareGrating2017-01-26';
	%Cells: 1440
	%Cort Conn: 403200
	%LGN Conn: 32256*2
	%Trials: 3600 / Types: 12 / Reps: 300
	%Oris: [0    15    30    45    60    75    90   105   120   135   150   165]
end

%% RUN: #header
strFigDir = 'D:\Data\Results\V1_LIFmodel\figures\';
if boolLoad
	runModelHeader;
end

%% calculate information present in data
matModelRespP = matModelResp;
intNeurons = size(matModelRespP,1);
vecGroupSizes = [2.^[0:10]];
%vecGroupSizes = [2.^[0:6]];
dblLambda = 0;
intIters = 5;
intGroups = numel(vecGroupSizes);

%pre-allocate
matI_LogReg_bc = nan(intIters,intGroups);
matI_Direct_bc = nan(intIters,intGroups);
matI_LogReg_bc_shuff = nan(intIters,intGroups);
matI_Direct_bc_shuff = nan(intIters,intGroups);

for intGroupSizeIdx=1:intGroups;
	intGroupSize = vecGroupSizes(intGroupSizeIdx)
	if intGroupSize > 100
		intIters = min([20 intIters]);
	end
	
	for intIter=1:intIters
		vecNeurons = randperm(intNeurons,intGroupSize);
		
		%select only classes 1 and 2
		matData = matModelRespP(vecNeurons,:);
		indClasses12 = (vecTrialOriIdx==1 | vecTrialOriIdx==2);
		matData12 = matData(:,indClasses12);
		vecTrialOriIdx12 = label2idx(vecTrialOriIdx(indClasses12));
		
		matThisD1 = matData12(:,vecTrialOriIdx12==1);
		matThisD2 = matData12(:,vecTrialOriIdx12==2);
		intTrials = size(matThisD1,2);
		
		%get correction factors
		dblDiffTheta = range(vecTrialOris(indClasses12));
		dblSubFac =(2*intGroupSize)/(intTrials*(dblDiffTheta.^2));
		dblProdFacRaw = ((2*intTrials-intGroupSize-3)/(2*intTrials-2));
		
		%get logistic regression output
		[vecWeightsLogReg, dblLLH] = doBinLogReg([matThisD1 matThisD2], [zeros(1,size(matThisD1,2)) ones(1,size(matThisD2,2))], dblLambda);
		vecClass1 = vecWeightsLogReg'*[matThisD1;ones(1,size(matThisD1,2))];
		vecClass2 = vecWeightsLogReg'*[matThisD2;ones(1,size(matThisD2,2))];
		dblDprimeLogReg = getdprime2(vecClass1,vecClass2);
		
		%get direct output
		[dblPredA,matPredA,dblD2,dblD2mat,dblD2_diag] = getSeparation([matThisD1 matThisD2],[zeros(1,size(matThisD1,2)) ones(1,size(matThisD2,2))],0);
		
		%save
		matI_LogReg_bc(intIter,intGroupSizeIdx) = (dblDprimeLogReg.^2)*dblProdFacRaw-dblSubFac;
		matI_Direct_bc(intIter,intGroupSizeIdx) = dblD2*dblProdFacRaw-dblSubFac;
		
		
		%% shuffled
		%shuffle
		matData12_shuff = nan(size(matData12));
		for intClass = [1 2];
			vecThisClass = find(vecTrialOriIdx12==intClass);
			for intThisNeuron=1:size(matData12_shuff,1)
				matData12_shuff(intThisNeuron,vecThisClass) = matData12(intThisNeuron,circshift(vecThisClass,intThisNeuron-1,2));
			end
		end
		matThisD1_shuff = matData12_shuff(:,vecTrialOriIdx12==1);
		matThisD2_shuff = matData12_shuff(:,vecTrialOriIdx12==2);
		
		%get logistic regression output
		[vecWeightsLogReg_shuff, dblLLH] = doBinLogReg([matThisD1_shuff matThisD2_shuff], [zeros(1,size(matThisD1_shuff,2)) ones(1,size(matThisD2_shuff,2))], dblLambda);
		vecClass1_shuff = vecWeightsLogReg_shuff'*[matThisD1_shuff;ones(1,size(matThisD1_shuff,2))];
		vecClass2_shuff = vecWeightsLogReg_shuff'*[matThisD2_shuff;ones(1,size(matThisD2_shuff,2))];
		dblDprimeLogReg_shuff = getdprime2(vecClass1_shuff,vecClass2_shuff);
		
		%calculate directly
		[dblPredA_shuff,matPredA_shuff,dblD2_shuff] = getSeparation([matThisD1_shuff matThisD2_shuff],[zeros(1,size(matThisD1_shuff,2)) ones(1,size(matThisD2_shuff,2))],0);
		
		%save
		matI_LogReg_bc_shuff(intIter,intGroupSizeIdx) = (dblDprimeLogReg_shuff.^2)*dblProdFacRaw-dblSubFac;
		matI_Direct_bc_shuff(intIter,intGroupSizeIdx) = dblD2_shuff*dblProdFacRaw-dblSubFac;
		
	end
end
clearvars matModelRespP;

% plot
%get mean+sd
vecI_LogRegM = nanmean(matI_LogReg_bc,1);
vecI_LogRegS = nanstd(matI_LogReg_bc,[],1);
vecI_Direct_bcM = nanmean(matI_Direct_bc,1);
vecI_Direct_bcS = nanstd(matI_Direct_bc,[],1);

vecI_LogReg_shuffM = nanmean(matI_LogReg_bc_shuff,1);
vecI_LogReg_shuffS = nanstd(matI_LogReg_bc_shuff,[],1);
vecI_Direct_bc_shuffM = nanmean(matI_Direct_bc_shuff,1);
vecI_Direct_bc_shuffS = nanstd(matI_Direct_bc_shuff,[],1);


%bias corrected
%vecI_shuff_bc %direct
figure
errorbar(vecGroupSizes,vecI_LogRegM,vecI_LogRegS,'b-');
hold on
errorbar(vecGroupSizes,vecI_LogReg_shuffM,vecI_LogReg_shuffS,'r-');
errorbar(vecGroupSizes,vecI_Direct_bcM,vecI_Direct_bcS,'b--');
errorbar(vecGroupSizes,vecI_Direct_bc_shuffM,vecI_Direct_bc_shuffS,'r--');
hold off
xlabel('Sample size (neurons)')
ylabel('Fisher information (d''^2)')
title('Blue=raw, red=shuffled; solid=Log reg; dashed=direct bc')
fixfig

%% perform V1-V1 prediction
%% summary
%data:
%88-159 (mean: 112.8) neurons in V1
%24-37 (mean: 29.4) neurons in V2
%stim:
%8 orientations, 22.5 degs, 300-400 repetitions (1.28 stim dur, 1.5 ITI)
%10 spike bins of 100 ms, starting 160 ms after stim onset
%subtract mean response per bin per neuron for residuals, exclude neurons <0.5 spikes / sec
%V1-V2:
%mean-matched sample of same number of neurons [15-31] from V1 as recorded in V2
%repeat mean matching 25 times and average across repeats
%% full regression model prediction to estimate maximal V1-V1 prediction
%Y = X * B
%B_ols = inv(X' * X) * X' * Y
%B_ridge = inv(X' * X + dblLambda * eye) * X' * Y
%largest lamdba within 1 SEM of best performance with 10=fold nested CV (see ref 36)

%get data
indTrials = vecTrialOriIdx==1;
intSamples = sum(indTrials);
intSizeX = 1000;
intSizeY = 500;
matModelRespRR = matModelResp;
matMiniSet = zscore(matModelRespRR(1:intSizeX,indTrials),[],2);

%set inputs
dblLambda = 1000;
X = matMiniSet((intSizeY+1):end,:)'; %source population (predictor)
Y = matMiniSet(1:intSizeY,:)'; %target population to be predicted

%perform ordinary least-squares (OLS) or ridge regression
B_ols = (X' * X) \ (X' * Y); %left-divide is same as inverse and multiplication
B_ridge = (X' * X + dblLambda*eye(intSizeY)) \ (X' * Y); %left-divide is same as inverse and multiplication

%predict responses
Y_pred_ridge = X * B_ridge;
Y_pred_ols = X * B_ols;

% compute R^2
vecMu = mean(Y);
dblSSRes_ridge = sum(sum((Y - Y_pred_ridge).^2));
dblSSRes_ols = sum(sum((Y - Y_pred_ols).^2));
dblSSTot = sum(sum(bsxfun(@minus,Y,vecMu).^2));
dblR2_ridge = 1 - dblSSRes_ridge / dblSSTot;
dblR2_ols = 1 - dblSSRes_ols / dblSSTot;

		
%% reduced rank regression to investigate dimensionality of predictive subspace
%B_rrr = B_ols * V * V'
%V contains top m principal components from Y = X * B
%Y = X * B_rrr = X * B_ols * V * V' = X * ^B * V'
%where ^B = B * V
%{
Req: 1 <= t <= s <= r
s = neurons in V2
n = trials
r = neurons in V1
t = rank of predictor model (projection matrix)

Y[s*n] = mu[s*n] + C[s*r] * X[r*n] + err[s*n]
C = A[s*t] * B[t*r]
%}
indTrials = vecTrialOriIdx==1;
intSamples = sum(indTrials);
intSizeX = 1000;
intSizeY = 500;
matModelRespRR = matModelResp;
matMiniSet = zscore(matModelRespRR(1:intSizeX,indTrials),[],2);
matX = matMiniSet((intSizeY+1):end,:)';
matY = matMiniSet(1:intSizeY,:)';

vecRank = 2:2:intSizeY;
vecR2 = zeros(size(vecRank));
intC = 0;
for intRankInT = vecRank
	intRankInT
	[matC, dblMSE, intRankOutT, sSuppOut] = doRdRankReg(matX, matY, 'rank', intRankInT);
	intC = intC + 1;
	vecR2(intC) = sSuppOut.dblR2;
end
plot(vecRank,vecR2)

%% factor analysis to investigate dimensionality of population dynamics subspace (max var, "dominant")
%select one stimulus type
indTrials = vecTrialOriIdx==1;
matDataFA = zscore(matModelRespRR(:,indTrials),[],2)';
intSamples = sum(indTrials);

%select training and test samples
vecTrainSamples = randperm(intSamples,round(0.9*intSamples));
vecTestSamples = 1:intSamples;
vecTestSamples(vecTrainSamples) = [];

%select neurons
intSizePop = 30;
vecNeurons = randperm(intNeurons,intSizePop);

%split test/train for selected population
matDataFA_All =  matDataFA(:,vecNeurons);
matDataFA_Train = matDataFA(vecTrainSamples,vecNeurons);
matDataFA_Test = matDataFA(vecTestSamples,vecNeurons);

%set variables for FA
intMaxRank = 30;
intRankInT = 2;
vecRank = 1:1:intMaxRank;
vecLogLikeTest = nan(size(vecRank));
vecLogLikeTrain = nan(size(vecRank));
vecTrace = nan(size(vecRank));
intC = 0;
intMaxC = length(vecRank);
intMaxRank = max(vecRank);
for intC=1:intMaxC
	intRankInT = vecRank(intC);
	fprintf('Now at rank %d/%d (%d/%d) [%s]\n',intRankInT,intMaxRank,intC,intMaxC,getTime);drawnow;
	
	%version 1
	[L,Psi,matRotation,stats,F] = factoran(matDataFA_Train,intRankInT);
	F = F';
	X_FA = matDataFA_Train';
	vecMu = xmean(X_FA,2);
	X_mu = repmat(vecMu,[1 size(X_FA,2)]);
	Psi = diag(Psi);
	CV_shared = L*L';
	CV_tot = CV_shared + Psi;
	
	mLH = X_FA -X_mu;
	mRH = L*F;
	mEps = mLH-mRH;
	
	%get likelihood of observing data
	y=mvnpdf(matDataFA_Test,vecMu',CV_tot);
	
	vecLogLikeTest(intC) = mean(log(y));
	vecLogLikeTrain(intC) = stats.loglike;
	vecTrace(intC) = trace(CV_shared);
end

% plot
figure
subplot(2,2,1)
plot(vecRank,vecTrace)
xlabel('Dominant dimensions (FA factors)');
ylabel('Tr(LL^T)');
fixfig

subplot(2,2,2)
plot(vecRank,vecLogLikeTest)
xlabel('Dominant dimensions (FA factors)');
ylabel('Test: Log(P(Test data | Training [LL^T + \Psi] )');
fixfig

subplot(2,2,3)
plot(vecRank,vecLogLikeTrain)
xlabel('Dominant dimensions (FA factors)');
ylabel('Train: Log(P(Training data | Training [LL^T + \Psi] )');
fixfig

if boolSaveFigs
	export_fig([strFigDir strSimulation '_PopComplexity_FA.tif']);
	export_fig([strFigDir strSimulation '_PopComplexity_FA.pdf']);
end

% SELECT m THAT MAXIMIZES CV LOG-LIKELIHOOD
[dblLL,intIdx] = max(vecLogLikeTest);
intM = vecRank(intIdx);

% DEFINE d AS NUMBER OF DIMS REQUIRED FOR EXPLAINING 95% OF LL'
%get full FA model
[L,Psi,matRotation,stats,F] = factoran(matDataFA_All,intM);
Psi = diag(Psi);
CV_shared = L*L';

%select dimensionality that explains >95% of LL'
boolConverged = false;
intOptDim = 0;
while ~boolConverged
	intOptDim = intOptDim + 1;
	[V,D] = eigs(CV_shared,intOptDim);
	dblExplainedSharedVar = trace(V'*CV_shared*V) / trace(CV_shared);
	if dblExplainedSharedVar > 0.95,boolConverged = true;end
end

%calculate percent shared variance per neuron
vecSharedV = nan(1,intSizePop);
vecPercSharedV = nan(1,intSizePop);
for intNeuron=1:intSizePop
	dblSharedV = L(intNeuron,:) * L(intNeuron,:)';
	dblPercSharedV = dblSharedV / (dblSharedV + Psi(intNeuron,intNeuron));
	
	vecSharedV(intNeuron) = dblSharedV;
	vecPercSharedV(intNeuron) = dblPercSharedV;
end

%calculate dominant modes (including both shared and private variance)
[V_full,D_full] = eig(CV_shared);
[D_full_desc,vecReorder]=sort(diag(D_full),'descend');
V_full_desc = V_full(:,vecReorder);
dblSumEigs = sum(D_full_desc(1:intOptDim));
vecSharedVarExplainedByMode = D_full_desc/dblSumEigs;

%calculate dominant modes for mode i and neuron n (including both shared and private variance)
matSharedVarModeI_NeuronN = nan(intSizePop,intSizePop);
for intMode=1:intSizePop
	for intNeuron=1:intSizePop
		dbl_u_i = V_full_desc(intNeuron,intMode);
		
		matSharedVarModeI_NeuronN(intNeuron,intMode) = (D_full_desc(intMode) * (dbl_u_i.^2)) / ( (L(intNeuron,:) * L(intNeuron,:)') + Psi(intNeuron,intNeuron) );
	end
end
vecSharedVarNeurons = sum(matSharedVarModeI_NeuronN,2);

%% check uniqueness of predictive subspace by iteratively removing predictive dimensions and redoing predictions [should quickly drop]
intRankB_hat = 5; %how many predictive dimensions to remove

%get OLS regression (full model)
matSigma = cov(matX);
B_ols_sub = (matX' * matX) \ (matX' * matY); %left-divide is same as inverse and multiplication

%reduced rank regression
[matC, dblMSE, intRankOutT, sSuppOut] = doRdRankReg(matX, matY, 'rank', intRankB_hat); %perform reduced rank regression
matV_rr = sSuppOut.matV_rr; %get reduced rank V

%get B-hat
B_hat = B_ols_sub * matV_rr; %define B-hat as the top predictive dimensions, using B and reduced rank version of V
matM = B_hat' * matSigma; %combine B-hat with covariance matrix

%get orthogonal basis Q for predictive nullspace
[U,S,V] = svd(matM,0); %singular value decomposition of matrix M
Q = V(:,(intRankB_hat+1):end); %orthonormal basis for uncorrelated (non-predictive) subspace

% project source activity onto uncorrelated subspace
X_hat = matX * Q; %uncorrelated subspace of X

%ridge regression between X_hat and Y
B_ridge_hat = (X_hat' * X_hat + dblLambda*eye(size(X_hat,2))) \ (X_hat' * matY); %left-divide is same as inverse and multiplication
Y_pred_ridge_hat = X_hat * B_ridge_hat;

% compute R^2
vecMu = mean(matY);
dblSSRes_ridge_hat = sum(sum((matY - Y_pred_ridge_hat).^2));
dblSSTot = sum(sum(bsxfun(@minus,matY,vecMu).^2));
dblR2_ridge_hat = 1 - dblSSRes_ridge_hat / dblSSTot;



%% check alignment of population dynamics subspace (max var, "dominant") with predictive subspace by iteratively removing pop. dyn. subsp. and redoing predictions [should not drop as quickly]
%calculate dominant modes (including both shared and private variance) for optimal model [intOptDim]
[L,Psi,matRotation,stats,F] = factoran(matDataFA_All,intOptDim);
Psi = diag(Psi);
CV_shared = L*L';

%SVD and eigen decomposition to select top N dominant dimensions
[U_l,S_l,V_l] = svd(L,0);
matZ = (S_l * V_l' * L' * inv(L*L' + Psi) * matDataFA_All')';

%ridge regression between reduced space X and Y
B_ridge_Z = (matZ' * matZ + dblLambda*eye(size(matZ,2))) \ (matZ' * matY); %left-divide is same as inverse and multiplication
Y_pred_Z = matZ * B_ridge_Z;

% compute R^2
vecMu = mean(matY);
dblSSRes_Z = sum(sum((matY - Y_pred_Z).^2));
dblSSTot = sum(sum(bsxfun(@minus,matY,vecMu).^2));
dblR2_Z = 1 - dblSSRes_Z / dblSSTot;





