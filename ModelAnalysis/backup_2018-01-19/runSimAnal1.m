%% initialize
clearvars;
boolLoad = false;
boolSaveFigs = true;
if ~exist('matModelResp','var')
	clearvars -except bool*;%close all;
	intLoadSim = 11; %% SET WHICH SIM TO LOAD HERE; 41 vs 42; 2 vs. 5
	boolLoad = true;
end
if intLoadSim == 11 && boolLoad
	strSimulation = 'xAreaDistributed_OriFull_2017-06-15';
	
elseif intLoadSim == 12 && boolLoad
	strSimulation = 'OldRerunSquareGrating_736803.458055_2017-04-18';
	
	
	
elseif intLoadSim == 21 && boolLoad
	strSimulation = 'SimpleSquareGrating2017-04-07';
	
elseif intLoadSim == 22 && boolLoad
	strSimulation = 'SquareGrating_736792.187807_2017-04-07';
	
	
elseif intLoadSim == 31 && boolLoad
	strSimulation = 'xAreaDistributed_491429_2017-04-21';
elseif intLoadSim == 32 && boolLoad
	strSimulation = 'xAreaDistributed_498422_2017-04-28';
elseif intLoadSim == 33 && boolLoad
	strSimulation = 'xAreaDistributed_491429_2017-05-02';
	
elseif intLoadSim == 41 && boolLoad
	strSimulation = 'xAreaDistributed_491426_2017-04-03';
	
elseif intLoadSim == 42 && boolLoad
	strSimulation = 'xAreaDistributed_493425_2017-04-05';
	
elseif intLoadSim == 43 && boolLoad
	strSimulation = 'SquareGrating2017-04-05';
	
elseif intLoadSim == 5 && boolLoad
	strSimulation = 'SquareGrating2017-01-25';
	%Cells: 1440
	%Cort Conn: 403200
	%LGN Conn: 32256*2
	%Trials: 3000 / Types: 2 / Reps: 1500
	%Oris: [42.5 47.5]
	
elseif intLoadSim == 51 && boolLoad
	strSimulation = 'OldRerun150SquareGrating_736796.588630_2017-04-11';
	%original connectivity structure
elseif intLoadSim == 52 && boolLoad
	strSimulation = 'OldRerun150FixedSquareGrating_736796.589038_2017-04-11';
	
elseif intLoadSim == 6 && boolLoad
	strSimulation = 'SquareGrating2017-01-26';
	%Cells: 1440
	%Cort Conn: 403200
	%LGN Conn: 32256*2
	%Trials: 3600 / Types: 12 / Reps: 300
	%Oris: [0    15    30    45    60    75    90   105   120   135   150   165]
	%elseif intLoadSim == 7 && boolLoad
	%	strSimulation = 'SquareGrating2017-03-15';
	%Cells: 240
	%Cort Conn: 16800
	%LGN Conn: 5376*2
	%Trials: 20 / Types: 2 / Reps: 10
	%Oris: [0    15    30    45    60    75    90   105   120   135   150   165]
	%elseif intLoadSim == 8 && boolLoad
	%strSimulation = 'xAreaLine2017-03-21';
	%	strSimulation = 'xAreaSquareGrating2017-03-22';
	%elseif intLoadSim == 9 && boolLoad
	%strSimulation = 'xAreaLine2017-03-21';
	%	strSimulation = 'SquareGrating_736781.722426_2017-03-27';
end

%% RUN: #header
strFigDir = 'D:\Data\Results\V1_LIFmodel\figures\';
if boolLoad
	runModelHeader;
end

%% parameters
cellC = {'000','010','100'};
intUseStimType = 3;
dblCutOff = 0.80;
dblLambda = 0;
dblNoiseLevel = 1;
strSimulation = ['Contrast' cellC{intUseStimType}];

%% add noise to improve numerical stability
matModelResp = double(matModelResp);
intNeurons = size(matModelResp,1);
intTrials = size(matModelResp,2);
vecNeuronSD = xstd(matModelResp,2);
matNoise = normrnd(zeros(intNeurons,intTrials),repmat(vecNeuronSD*dblNoiseLevel,[1 intTrials]));
matModelRespP = matModelResp + matNoise;

%% calculate information present in data
if numel(unique(vecTrialStimType)) > 1
	vecUseStimTypes = find(vecTrialContrasts==100);
	intNeurons = size(matModelRespP,1);
	intTrials = size(matModelRespP,2);
	%vecGroupSizes = [2.^[0:3]];
	vecGroupSizes = [2.^[0:20]];
	vecGroupSizes(vecGroupSizes>min([intNeurons intTrials]))=[];
	vecGroupSizes((end-1):end) = [];
	vecGroupSizes(8:end) = [];
	dblLambda = 0;
	intIters = 100;
	intGroups = numel(vecGroupSizes);
	
	%pre-allocate
	matI_LogReg_bc = nan(intIters,intGroups);
	matI_Direct_bc = nan(intIters,intGroups);
	matI_LogReg_bc_shuff = nan(intIters,intGroups);
	matI_Direct_bc_shuff = nan(intIters,intGroups);
	
	for intGroupSizeIdx=1:intGroups;
		intGroupSize = vecGroupSizes(intGroupSizeIdx);
		if intGroupSize > 100
			intIters = min([20 intIters]);
		end
		
		for intIter=1:intIters
			vecNeurons = randperm(intNeurons,intGroupSize);
			
			%select only classes 1 and 2
			matData = matModelRespP(vecNeurons,:);
			indClasses12 = (vecTrialStimType==vecUseStimTypes(1) | vecTrialStimType==vecUseStimTypes(2));
			matData12 = matData(:,indClasses12);
			vecTrialStimIdx12 = label2idx(vecTrialStimType(indClasses12));
			
			matThisD1 = matData12(:,vecTrialStimIdx12==1);
			matThisD2 = matData12(:,vecTrialStimIdx12==2);
			intTrials = size(matThisD1,2);
			
			%get correction factors
			if exist('vecOrientations','var') && numel(vecOrientations) > 1
				dblDiffTheta = range(vecOrientations([1 2]));
			else
				dblDiffTheta = 1;
			end
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
				vecThisClass = find(vecTrialStimIdx12==intClass);
				for intThisNeuron=1:size(matData12_shuff,1)
					matData12_shuff(intThisNeuron,vecThisClass) = matData12(intThisNeuron,circshift(vecThisClass,intThisNeuron-1,2));
				end
			end
			matThisD1_shuff = matData12_shuff(:,vecTrialStimIdx12==1);
			matThisD2_shuff = matData12_shuff(:,vecTrialStimIdx12==2);
			
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
	%clearvars matModelRespP;
	
	%% plot
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
	errorbar(vecGroupSizes+1,vecI_Direct_bcM,vecI_Direct_bcS,'b--');
	errorbar(vecGroupSizes+1,vecI_Direct_bc_shuffM,vecI_Direct_bc_shuffS,'r--');
	hold off
	xlabel('Sample size (neurons)')
	ylabel('Fisher information (d''^2)')
	title('Blue=raw, red=shuffled; solid=Log reg; dashed=direct bc')
	fixfig
	
	%save figure
	if boolSaveFigs
		export_fig([strFigDir  'Information_' strSimulation '.tif']);
		export_fig([strFigDir  'Information_' strSimulation '.pdf']);
	end
end
%% perform V1-V1 prediction
%% summary
%data:
% 5 sessions
%88-159 (mean: 112.8) neurons in V1
%24-37 (mean: 29.4) neurons in V2
%stim:
%8 orientations, 22.5 degs, 300-400 repetitions (1.28 stim dur, 1.5 ITI)
%10 spike bins of 100 ms, starting 160 ms after stim onset
%subtract mean response per bin per neuron for residuals, exclude neurons <0.5 spikes / sec
%V1-V2:
%mean-matched sample of same number of neurons [15-31] from V1 as recorded in V2
%repeat mean matching 25 times and average across repeats


%% make figures
%% define data for random selection
%get data
indTrials = vecTrialStimType==intUseStimType;
intSamples = sum(indTrials);
intSizeX = 110;%110
intSizeY = 30;%30
intResamplings = 10;
matModelRespZ = zscore(matModelRespP(:,indTrials),[],2);
%intNeuronsV1 = size(matModelRespZ,1);
if exist('intCellsV1','var')
	intNeuronsV1 = intCellsV1;
else
	intNeuronsV1 = intCortexCells;
end
%build array for resamplings
cellMatX = cell(1,intResamplings);
cellNeuronsX = cell(1,intResamplings);
cellMatY = cell(1,intResamplings);
cellNeuronsY = cell(1,intResamplings);
for intResampling=1:intResamplings
	%% select data
	vecNeuronsX = randperm(intNeuronsV1,intSizeX);
	vecOthers = find(~ismember(1:intNeuronsV1,vecNeuronsX));
	vecNeuronsY = vecOthers(randperm(length(vecOthers),intSizeY));
	
	%get data
	matX = matModelRespZ(vecNeuronsX,:)'; %source population (predictor)
	matY = matModelRespZ(vecNeuronsY,:)'; %target population to be predicted
	
	%save data
	cellMatX{intResampling} = matX;
	cellNeuronsX{intResampling} = vecNeuronsX;
	cellMatY{intResampling} = matY;
	cellNeuronsY{intResampling} = vecNeuronsY;
end

%% single neuron prediction
vecPredictionsR2_SingleNeurons = nan(1,intNeuronsV1);
for intNeuron=1:intNeuronsV1
	%% select data
	vecNeuronsX = intNeuron;
	vecOthers = find(~ismember(1:intNeuronsV1,vecNeuronsX));
	vecNeuronsY = vecOthers(randperm(length(vecOthers),intSizeY));
	
	%get data
	matX = matModelRespZ(vecNeuronsX,:)'; %source population (predictor)
	matY = matModelRespZ(vecNeuronsY,:)'; %target population to be predicted
	
	%% figure 2b; predictive performance whole population/single neuron
	%perform ordinary least-squares (OLS) or ridge regression
	B_ridge = (matX' * matX + dblLambda) \ (matX' * matY); %left-divide is same as inverse and multiplication
	
	%predict responses
	matY_pred_ridge = matX * B_ridge;
	
	% compute R^2
	vecMu = mean(matY);
	dblSSRes_ridge = sum(sum((matY - matY_pred_ridge).^2));
	dblSSTot = sum(sum(bsxfun(@minus,matY,vecMu).^2));
	dblR2_ridge = 1 - dblSSRes_ridge / dblSSTot;
	
	%save data
	vecPredictionsR2_SingleNeurons(intNeuron) = dblR2_ridge;
end

%% full regression model prediction to estimate maximal V1-V1 prediction
vecPredictionsR2 = nan(1,intResamplings);
for intResampling=1:intResamplings
	%% select data
	matX = cellMatX{intResampling};
	matY = cellMatY{intResampling};
	
	%% figure 2b; predictive performance whole population/single neuron
	%perform ordinary least-squares (OLS) or ridge regression
	B_ols = (matX' * matX) \ (matX' * matY); %left-divide is same as inverse and multiplication
	B_ridge = (matX' * matX + dblLambda*eye(intSizeX)) \ (matX' * matY); %left-divide is same as inverse and multiplication
	
	%predict responses
	matY_pred_ridge = matX * B_ridge;
	matY_pred_ols = matX * B_ols;
	
	% compute R^2
	vecMu = mean(matY);
	dblSSRes_ridge = sum(sum((matY - matY_pred_ridge).^2));
	dblSSRes_ols = sum(sum((matY - matY_pred_ols).^2));
	dblSSTot = sum(sum(bsxfun(@minus,matY,vecMu).^2));
	dblR2_ridge = 1 - dblSSRes_ridge / dblSSTot;
	dblR2_ols = 1 - dblSSRes_ols / dblSSTot;
	
	%save data
	vecPredictionsR2(intResampling) = dblR2_ridge;
end
%% plot
figure
histx(vecPredictionsR2);
vecLimY = get(gca,'ylim');
hold on
plot(mean(vecPredictionsR2_SingleNeurons)*[1 1],vecLimY,'b--');
text(mean(vecPredictionsR2_SingleNeurons)-0.04,max(vecLimY)/3,'Single neuron perf.','Color','blue','FontSize',14,'Rotation',90);
hold off
fixfig;
xlim([0 1])
xlabel('Predictive performance')
ylabel('Counts')
title('Semedo Fig 2')

%save figure
if boolSaveFigs
	export_fig([strFigDir  'SemedoFig2_PredPerf' strSimulation '.tif']);
	export_fig([strFigDir  'SemedoFig2_PredPerf' strSimulation '.pdf']);
end

%% reduced rank regression to investigate dimensionality of predictive subspace
vecRank = 1:intSizeY;
matPredDimDepR2 = nan(length(vecRank),intResamplings);
for intResampling=1:intResamplings
	%% select data
	matX = cellMatX{intResampling};
	matY = cellMatY{intResampling};
	
	%% figure 4; full prediction versus reduced rank regression, as function of dimensionality
	vecR2 = nan(size(vecRank));
	intC = 0;
	for intRankInT = vecRank
		[matC, dblMSE, intRankOutT, sSuppOut] = doRdRankReg(matX, matY, 'rank', intRankInT);
		intC = intC + 1;
		vecR2(intC) = sSuppOut.dblR2;
	end
	%put in matrix
	matPredDimDepR2(:,intResampling) = vecR2;
end
vecPredictiveDimensions = sum(bsxfun(@rdivide,matPredDimDepR2,max(matPredDimDepR2,[],1))<dblCutOff,1)+1;

%% plot
figure
errorbar(vecRank,xmean(matPredDimDepR2,2),xstd(matPredDimDepR2,2))
xlim([0 max(vecRank)]);
hold on
scatter(0.5,mean(vecPredictionsR2),80,'xb');
text(1,mean(vecPredictionsR2)+0.02,'Full prediction','Color','blue','FontSize',14,'Rotation',45);
hold off
vecLimY = get(gca,'ylim');
ylim([0 vecLimY(end)]);
fixfig
xlabel('Number of predictive dimensions');
ylabel('Predictive performance');
title('Semedo Fig 4')

%save figure
if boolSaveFigs
	export_fig([strFigDir  'SemedoFig4_PredDim' strSimulation '.tif']);
	export_fig([strFigDir  'SemedoFig4_PredDim' strSimulation '.pdf']);
end

%% factor analysis to investigate dimensionality of population dynamics subspace (max var, "dominant")
intMaxRank = round(intSizeY*(2/3));
vecRankDom = 1:intMaxRank;
vecTargetPopDimensionality = nan(1,intResamplings);
matTargetPopLogLikeTest = nan(length(vecRankDom),intResamplings);
for intResampling=1:intResamplings
	%% select data
	%matX = cellMatX{intResampling};
	matY = cellMatY{intResampling};
	
	%% figure 5; dimensionality of source and target populations (FA)
	%split matY in 90-10 train-test pairs
	vecTrainSamples = randperm(intSamples,round(0.9*intSamples));
	vecTestSamples = 1:intSamples;
	vecTestSamples(vecTrainSamples) = [];
	
	%split test/train for selected population
	matDataFA_All =  matY(:,:);
	matDataFA_Train = matY(vecTrainSamples,:);
	matDataFA_Test = matY(vecTestSamples,:);
	
	%set variables for FA
	vecLogLikeTest = nan(size(vecRankDom));
	vecLogLikeTrain = nan(size(vecRankDom));
	vecTrace = nan(size(vecRankDom));
	intMaxC = length(vecRankDom);
	fprintf('Starting resampling %d/%d [%s]\n',intResampling,intResamplings,getTime);drawnow;
	for intC=1:intMaxC
		intRankInT = vecRankDom(intC);
		
		%version 1
		%warning('off');
		[L, mu, Psi, llh] = factorAnal(matDataFA_Train,intRankInT);
		%[L2,Psi2,matRotation,stats,F] = factoran(matDataFA_Train,intRankInT);
		%F = F';
		X_FA = matDataFA_Train';
		vecMu = xmean(X_FA,2);
		X_mu = repmat(vecMu,[1 size(X_FA,2)]);
		Psi = diag(Psi);
		CV_shared = L*L';
		CV_tot = CV_shared + Psi;
		
		%get likelihood of observing data
		y=mvnpdf(matDataFA_Test,vecMu',CV_tot);
		
		vecLogLikeTest(intC) = mean(log(y));
		%vecLogLikeTrain(intC) = stats.loglike;
		vecLogLikeTrain(intC) = -(llh(end)-llh(end-1));
		vecTrace(intC) = trace(CV_shared);
	end
	matTargetPopLogLikeTest(:,intResampling) = vecLogLikeTest;
	
	% SELECT m THAT MAXIMIZES CV LOG-LIKELIHOOD
	[dblLL,intIdx] = max(vecLogLikeTest);
	intM = vecRankDom(intIdx);
	
	% DEFINE d AS NUMBER OF DIMS REQUIRED FOR EXPLAINING 95% OF LL'
	%get full FA model
	%[L,Psi,matRotation,stats,F] = factoran(matDataFA_All,intM);
	[L, mu, Psi] = factorAnal(matDataFA_All,intM);
	Psi = diag(Psi);
	CV_shared = L*L';
	
	%select dimensionality that explains >95% of LL'
	boolConverged = false;
	intOptDim = 0;
	while ~boolConverged
		intOptDim = intOptDim + 1;
		[V,D] = eigs(CV_shared,intOptDim);
		dblExplainedSharedVar = trace(V'*CV_shared*V) / trace(CV_shared);
		if dblExplainedSharedVar > dblCutOff,boolConverged = true;end
	end
	vecTargetPopDimensionality(intResampling) = intOptDim;
end

%% plot
figure
subplot(2,3,1)
errorbar(vecRankDom,xmean(matTargetPopLogLikeTest,2),xstd(matTargetPopLogLikeTest,2))
xlabel('Target population dimensionality');
ylabel('Log likelihood of CV test data');
fixfig

subplot(2,3,2)
plot([0 10],[0 10],'k--');
hold on
scatter(vecTargetPopDimensionality+0.1*(rand(size(vecTargetPopDimensionality))-0.5),vecPredictiveDimensions+0.1*(rand(size(vecPredictiveDimensions))-0.5))
hold off
%xlim([0 10])
%ylim([0 10])
xlabel('Target population dimensionality')
ylabel('Number of predictive dimensions')
fixfig
title('Semedo Fig 5')

%full screen
drawnow;
jFig = get(handle(gcf), 'JavaFrame');
jFig.setMaximized(true);
figure(gcf);
drawnow;

%save figure
if boolSaveFigs
	export_fig([strFigDir  'SemedoFig5_DomDim' strSimulation '.tif']);
	export_fig([strFigDir  'SemedoFig5_DomDim' strSimulation '.pdf']);
end

%% check uniqueness of predictive subspace by iteratively removing predictive dimensions and redoing predictions [should quickly drop]
vecRemovePredictiveDimensions = 0:max(vecPredictiveDimensions);
matR2RemPredDim = nan(length(vecRemovePredictiveDimensions),intResamplings);
for intResampling=1:intResamplings
	%% select data
	matX = cellMatX{intResampling};
	matY = cellMatY{intResampling};
	
	%% figure 6; prediction when removing predictive dimensions
	for intRankC = 1:length(vecRemovePredictiveDimensions)
		intRankB_hat = vecRemovePredictiveDimensions(intRankC); %how many predictive dimensions to remove
		if intRankB_hat == 0
			matR2RemPredDim(intRankC,intResampling) = vecPredictionsR2(intResampling);
		else
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
			matR2RemPredDim(intRankC,intResampling) = dblR2_ridge_hat;
		end
	end
end
%% plot
figure
subplot(2,3,1)
errorbar(vecRemovePredictiveDimensions,xmean(matR2RemPredDim,2),xstd(matR2RemPredDim,2))
xlim([min(vecRemovePredictiveDimensions)-0.1 max(vecRemovePredictiveDimensions)+0.1]);
vecLimY = get(gca,'ylim');
ylim([0 vecLimY(end)]);
xlabel('Number of predictive dimensions removed');
ylabel('Predictive performance (R^2)')
fixfig

subplot(2,3,2)
vecR2AllRem = getMatVals(matR2RemPredDim,vecPredictiveDimensions,1:intResamplings);
histx(vecR2AllRem)
xlim([0 1]);
xlabel('Performance after removing all pred dims');
ylabel('Counts')
fixfig
title('Semedo Fig 6')

%full screen
drawnow;
jFig = get(handle(gcf), 'JavaFrame');
jFig.setMaximized(true);
figure(gcf);
drawnow;

%save figure
if boolSaveFigs
	export_fig([strFigDir  'SemedoFig6_PredRem' strSimulation '.tif']);
	export_fig([strFigDir  'SemedoFig6_PredRem' strSimulation '.pdf']);
end

%% check alignment of population dynamics subspace (max var, "dominant") with predictive subspace by iteratively removing pop. dyn. subsp. and redoing predictions [should not drop as quickly]
vecUseDomDims = 1:10;
matDomDimDepR2 = nan(length(vecUseDomDims),intResamplings);
for intResampling=1:intResamplings
	%% select data
	matX = cellMatX{intResampling};
	matY = cellMatY{intResampling};
	
	%% figure 7; prediction comparing inclusion of predictive or dominant dimensions
	fprintf('Starting resampling %d/%d [%s]\n',intResampling,intResamplings,getTime);drawnow;
	for intDDC=1:length(vecUseDomDims)
		intDomDim = vecUseDomDims(intDDC);
		%calculate dominant modes (including both shared and private variance) for optimal model [intOptDim]
		%try
			%[L,Psi,matRotation,stats,F] = factoran(matX,intDomDim);
			[L, mu, Psi] = factorAnal(matX,intDomDim);
		%catch
		%	fprintf('Oops %d|%d [%s]\n',intResampling,intDDC,getTime);
		%	continue;
		%end
		
		
		Psi = diag(Psi);
		CV_shared = L*L';
		
		%SVD and eigen decomposition to select top N dominant dimensions
		[U_l,S_l,V_l] = svd(L,0);
		matZ = (((S_l * V_l' * L') / (L*L' + Psi)) * matX')';
		
		%ridge regression between reduced space X and Y
		B_ridge_Z = (matZ' * matZ + dblLambda*eye(size(matZ,2))) \ (matZ' * matY); %left-divide is same as inverse and multiplication
		Y_pred_Z = matZ * B_ridge_Z;
		
		% compute R^2
		vecMu = mean(matY);
		dblSSRes_Z = sum(sum((matY - Y_pred_Z).^2));
		dblSSTot = sum(sum(bsxfun(@minus,matY,vecMu).^2));
		dblR2_Z = 1 - dblSSRes_Z / dblSSTot;
		matDomDimDepR2(intDDC,intResampling) = dblR2_Z;
	end
end

%% plot
%performance as function of dimensions
figure
errorbar(vecUseDomDims,nanmean(matPredDimDepR2(vecUseDomDims,:),2),nanstd(matPredDimDepR2(vecUseDomDims,:),[],2),'x-');
hold on
errorbar(vecUseDomDims,nanmean(matDomDimDepR2(vecUseDomDims,:),2),nanstd(matDomDimDepR2(vecUseDomDims,:),[],2),'bo-');
hold off
xlim([0 max(vecUseDomDims)]);
vecLimY = get(gca,'ylim');
ylim([0 vecLimY(end)]);
fixfig
xlabel('Number of dimensions');
ylabel('Predictive performance');
legend('Predictive','Dominant','Location','Best')
title('Semedo Fig 7')

%save figure
if boolSaveFigs
	export_fig([strFigDir  'SemedoFig7_PredDom' strSimulation '.tif']);
	export_fig([strFigDir  'SemedoFig7_PredDom' strSimulation '.pdf']);
end

%dominant dimensions as a function of predictive dimensions

