%% initialize
boolLoad = true;
boolSaveFigs = false;
if ~exist('matModelResp','var')
	clearvars -except bool*;%close all;
	intLoadSim = 8; %% SET WHICH SIM TO LOAD HERE
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
elseif intLoadSim == 7 && boolLoad
	strSimulation = 'SquareGrating2017-03-15';
	%Cells: 240
	%Cort Conn: 16800
	%LGN Conn: 5376*2
	%Trials: 20 / Types: 2 / Reps: 10
	%Oris: [0    15    30    45    60    75    90   105   120   135   150   165]
elseif intLoadSim == 8 && boolLoad
	strSimulation = 'xAreaSquareGrating2017-03-20';
		
end

%% RUN: #header
strFigDir = 'D:\Data\Results\V1_LIFmodel\figures\';
if boolLoad
	runModelHeader;
end

%%
%{
matData = matModelResp(1:20,:);
[dblPerformanceCV,vecDecodedIndexCV,matPosteriorProbabilityCV,matWeights] = doCrossValidatedDecodingLR(matData,vecTrialTypes,0.1);
dblPerformanceCV

vecW1 = matWeights(:,1);
vecdMu12 = [nanmean(matData(:,vecTrialOriIdx==1),2); 1] - [nanmean(matData(:,vecTrialOriIdx==2),2); 1];

vecdMu12'*matWeights(:,2)

matDataPlusLin = [matData; ones(1,size(matData,2))];
matActivation = matWeights'*matDataPlusLin;
	
vecDeltaMu

matCov1 = cov(matData(:,vecTrialOriIdx==1)');
%}
%% reproduce fig. 4 from Series et al.
%matModelRespP = matModelResp(vecCellTypes==1,:);
matModelRespP = matModelResp;
intNeurons = size(matModelRespP,1);
vecGroupSizes = [2.^[0:10]];
%vecGroupSizes = [2.^[0:8]];
dblLambda = 0;
intIters = 1;
intGroups = numel(vecGroupSizes);

%pre-allocate
matI_LogReg_bc = nan(intIters,intGroups);
matI_Direct_bc = nan(intIters,intGroups);
matI_Direct = nan(intIters,intGroups);
matPerf = nan(intIters,intGroups);
matPerfPred = nan(intIters,intGroups);
		
matI_LogReg_bc_shuff = nan(intIters,intGroups);
matI_Direct_bc_shuff = nan(intIters,intGroups);
matI_Direct_shuff = nan(intIters,intGroups);
matPerf_shuff = nan(intIters,intGroups);
matPerfPred_shuff = nan(intIters,intGroups);

for intGroupSizeIdx=1:intGroups;
	intGroupSize = vecGroupSizes(intGroupSizeIdx)
	if intGroupSize > 100
		intIters = min([20 intIters]);
	end
	%vecNeurons = linspace(1,intNeurons,intGroupSize+1);
	%vecNeurons = ceil(vecNeurons(1:(end-1)));
	%vecNeurons = 1:intGroupSize;
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
		y = ~round(exp(-log1pexp([vecClass1 vecClass2])));
		dblPerfLogReg = sum(y==[zeros(1,numel(vecClass1)) ones(1,numel(vecClass2))])/numel(y);
		

		%get direct output
		[dblPredA,matPredA,dblD2,dblD2mat,dblD2_diag] = getSeparation([matThisD1 matThisD2],[zeros(1,size(matThisD1,2)) ones(1,size(matThisD2,2))],0);
		
		%save
		matI_LogReg_bc(intIter,intGroupSizeIdx) = (dblDprimeLogReg.^2);
		matI_Direct_bc(intIter,intGroupSizeIdx) = dblD2*dblProdFacRaw-dblSubFac;
		matI_Direct(intIter,intGroupSizeIdx) = dblD2;
		matPerf(intIter,intGroupSizeIdx) = dblPerfLogReg;
		matPerfPred(intIter,intGroupSizeIdx) = dblPredA;
		
		
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
		y = ~round(exp(-log1pexp([vecClass1_shuff vecClass2_shuff])));
		dblPerfLogReg_shuff = sum(y==[zeros(1,numel(vecClass1)) ones(1,numel(vecClass2))])/numel(y);
		
		%calculate directly
		[dblPredA_shuff,matPredA_shuff,dblD2_shuff] = getSeparation([matThisD1_shuff matThisD2_shuff],[zeros(1,size(matThisD1_shuff,2)) ones(1,size(matThisD2_shuff,2))],0);
		
		%save
		matI_LogReg_bc_shuff(intIter,intGroupSizeIdx) = (dblDprimeLogReg_shuff.^2);
		matI_Direct_bc_shuff(intIter,intGroupSizeIdx) = dblD2_shuff*dblProdFacRaw-dblSubFac;
		matI_Direct_shuff(intIter,intGroupSizeIdx) = dblD2_shuff;
		matPerf_shuff(intIter,intGroupSizeIdx) = dblPerfLogReg_shuff;
		matPerfPred_shuff(intIter,intGroupSizeIdx) = dblPredA_shuff;
		
	end
end
%% plot
% information
%logistic regression non-shuffled, CV 50/50 (cross)
%vecI

%logistic regression shuffled, CV 50/50 (cross)
%vecI_shuff

%get mean+sd
vecI_LogRegM = nanmean(matI_LogReg_bc,1);
vecI_LogRegS = nanstd(matI_LogReg_bc,[],1);
vecI_Direct_bcM = nanmean(matI_Direct_bc,1);
vecI_Direct_bcS = nanstd(matI_Direct_bc,[],1);
vecI_DirectM = nanmean(matI_Direct,1);
vecI_DirectS = nanstd(matI_Direct,[],1);
vecPerfM = nanmean(matPerf,1);
vecPerfS = nanstd(matPerf,[],1);
vecPerfPredM = nanmean(matPerfPred,1);
vecPerfPredS = nanstd(matPerfPred,[],1);

vecI_LogReg_shuffM = nanmean(matI_LogReg_bc_shuff,1);
vecI_LogReg_shuffS = nanstd(matI_LogReg_bc_shuff,[],1);
vecI_Direct_bc_shuffM = nanmean(matI_Direct_bc_shuff,1);
vecI_Direct_bc_shuffS = nanstd(matI_Direct_bc_shuff,[],1);
vecI_Direct_shuffM = nanmean(matI_Direct_shuff,1);
vecI_Direct_shuffS = nanstd(matI_Direct_shuff,[],1);
vecPerf_shuffM = nanmean(matPerf_shuff,1);
vecPerf_shuffS = nanstd(matPerf_shuff,[],1);
vecPerfPred_shuffM = nanmean(matPerfPred_shuff,1);
vecPerfPred_shuffS = nanstd(matPerfPred_shuff,[],1);

%bias corrected
%vecI_shuff_bc %direct
figure
subplot(2,2,1)
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

subplot(2,2,2)
errorbar(vecGroupSizes,vecI_DirectM,vecI_DirectS,'b-');
hold on
errorbar(vecGroupSizes,vecI_Direct_shuffM,vecI_Direct_shuffS,'r-');
errorbar(vecGroupSizes,vecI_Direct_bcM,vecI_Direct_bcS,'b--');
errorbar(vecGroupSizes,vecI_Direct_bc_shuffM,vecI_Direct_bc_shuffS,'r--');
hold off
xlabel('Sample size (neurons)')
ylabel('Fisher information (d''^2)')
title('Direct, blue=raw, red=shuffled; dashed=bias corr')
fixfig

% decoding accuracy
%vecPerf
%vecPerfPred


%vecPerf_shuff
%vecPerfPred_shuff
%vecPerfPred_shuff_bc

subplot(2,2,3)
errorbar(vecGroupSizes,vecPerfM,vecPerfS,'b-');
hold on
errorbar(vecGroupSizes,vecPerfPredM,vecPerfPredS,'b--');
hold off
xlabel('Sample size (neurons)')
ylabel('Decoding accuracy')
set(gca,'xscale','log')
title('Solid=real, dashed=predicted')
fixfig

subplot(2,2,4)
errorbar(vecGroupSizes,vecPerf_shuffM,vecPerf_shuffS,'r-');
hold on
errorbar(vecGroupSizes,vecPerfPred_shuffM,vecPerfPred_shuffS,'r--');
hold off
xlabel('Sample size (neurons)')
ylabel('Decoding accuracy')
set(gca,'xscale','log')
title('Solid=real, dashed=predicted')
fixfig
%{
subplot(2,3,6)
plot(vecGroupSizes,vecPerfPred_shuff_bc,'xm--');
hold off
xlabel('Sample size')
ylabel('Decoding accuracy')
set(gca,'xscale','log')
title('Solid=real, dashed=predicted')
fixfig
%}
%%

%save figure
drawnow;
jFig = get(handle(gcf), 'JavaFrame');
jFig.setMaximized(true);
figure(gcf);
drawnow;
export_fig(['D:\Data\Results\V1_LIFmodel\figures\' strSimulation '_InformationAndDecoding.tif']);
export_fig(['D:\Data\Results\V1_LIFmodel\figures\' strSimulation '_InformationAndDecoding.pdf']);
return

%% get relationship predicted/actual
intIters=10;
vecSizes = [2.^(0:8)];
intSizes = numel(vecSizes);
matDprimeSquared = nan(intIters,intSizes);
matPredictedAcc = nan(intIters,intSizes);
matActualAcc = nan(intIters,intSizes);
intNeuronCounter = 0;
indOriPair = vecTrialOriIdx==1 | vecTrialOriIdx==2;
matModelRespPair = matModelResp(:,indOriPair);
vecOriPair = vecTrialOriIdx(indOriPair);
for intSize=1:intSizes
	intGroupSize = vecSizes(intSize);
	fprintf('Running group size %d/%d (%d) [%s]\n',intSize,intSizes,intGroupSize,getTime);
	parfor intIter=1:intIters
		vecNeurons = randperm(intNeurons,intGroupSize);
		
		[dblPredA,matPredA,dblDprimeSquaredOffDiagonal] = getSeparation(matModelResp(vecNeurons,:),vecTrialOriIdx);
		matDprimeSquared(intIter,intSize) = dblDprimeSquaredOffDiagonal;
		matPredictedAcc(intIter,intSize) = dblPredA(end,1);
		matActualAcc(intIter,intSize) = doCrossValidatedDecodingLR(matModelResp(vecNeurons,:),vecTrialOriIdx,0.01);
	end
end

%

plot([0 1],[0 1],'k--')
hold on
scatter(mean(matActualAcc,1),mean(matPredictedAcc,1))
hold off
%xlim([0.5 1])
%ylim([0.5 1])

%% plot
clf
subplot(2,2,1);
plot(sort(vecDprimeSquared));
xlabel('Number of neurons');
ylabel('D''^2');
title('Greedy classifier');
xlim([0 intNeurons]);

subplot(2,2,2);
plot(sort(vecPredictedAcc));
xlabel('Number of neurons');
ylabel('Predicted decoding accuracy');
xlim([0 intNeurons]);
ylim([0 1]);
title('Decoding accuracy predicted from d''^2' )

subplot(2,2,3);
plot(vecPerformance);
xlabel('Number of neurons');
ylabel('Real decoding accuracy');
xlim([0 intNeurons]);
ylim([0 1]);
title('Real decoding accuracy, neurons sorted by d''^2 contribution' )

subplot(2,2,4)
title('Model')
axes off
return

export_fig(['D:\Data\Results\V1_LIFmodel\' strSimulation '_greedyDecoding2.tif']);
export_fig(['D:\Data\Results\V1_LIFmodel\' strSimulation '_greedyDecoding2.pdf']);



%% get neuron group
vecUseNeurons = [1 2 3];
matDataPoints = matModelResp(vecUseNeurons,:)';

%calc orth/para
[vecSD_Orth,vecSD_Para] = calcOrthPara(matDataPoints,vecTrialOriIndex);
[vecShuffSD_Orth,vecShuffSD_Para] = calcOrthPara(matDataPoints,vecTrialOriIndex,true);
clc
dblOrthRel = mean(vecSD_Orth)/mean(vecShuffSD_Orth)
dblParaRel = mean(vecSD_Para)/mean(vecShuffSD_Para)

%transform matrix
intStimTypes = numel(vecOrientations);
intReps = intTrials/intStimTypes;
matStimResp = nan(intNeurons,intStimTypes,intReps);
for intStimType=1:intStimTypes
	matStimResp(:,intStimType,:) = matModelResp(:,vecTrialOriIdx==intStimType);
end

mat3Resp = matStimResp(vecUseNeurons,:,:);
plotTube(mat3Resp)

%% investigate orientation coding per spike

