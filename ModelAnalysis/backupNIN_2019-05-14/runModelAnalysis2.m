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

%% decode
intIters = 1;
intMaxSize = 30;
matPerfIndep = nan(intMaxSize,intIters);
matPerfMD = nan(intMaxSize,intIters);

parfor intSize=2:intMaxSize
	intSize
	for intIter = 1:intIters
		vecUseNeurons = randperm(intNeurons,intSize);
		matThisData = matModelResp(vecUseNeurons,:);
		%dblPerformanceIndep_NCV = doCrossValidatedDecodingML(matThisData,vecTrialOriIndex,0)
		%dblPerformanceIndep_LOO = doCrossValidatedDecodingML(matThisData,vecTrialOriIndex,1)
		dblPerformanceIndep_LRO = doCrossValidatedDecodingML(matThisData,vecTrialOriIndex,2);
		matPerfIndep(intSize,intIter) = dblPerformanceIndep_LRO;
		
		%dblPerformanceMultD_NCV = doCrossValidatedDecodingMD(matThisData,vecTrialOriIndex,0)
		%dblPerformanceMultD_LOO = doCrossValidatedDecodingMD(matThisData,vecTrialOriIndex,1)
		dblPerformanceMultD_LRO = doCrossValidatedDecodingMD(matThisData,vecTrialOriIndex,2);
		matPerfMD(intSize,intIter) = dblPerformanceMultD_LRO;
		intIter
	end
end
return
dblLambda = 0.01;
dblPerformanceLR = doCrossValidatedDecodingLR(matData,vecTrialOriIndex,dblLambda)


[dblPredA,matPredA,dblDprime2_Indep,matDprimeSquared] = getSeparation(matData,vecTrialOriIdx,1);
dblDprime2_Indep

[dblPredA_MD,matPredA_MD,dblDprime2_MultD,matDprimeSquared_MD] = getSeparation(matData,vecTrialOriIdx,0);
dblDprime2_MultD
%%


return
%% greedy separator
vecDprimeSquared = nan(1,intNeurons);
vecPredictedAcc = nan(1,intNeurons);
vecNeuronsUsedD2 = false(1,intNeurons);
vecNeuronsLeftD2 = true(1,intNeurons);
vecNeuronOrderD2 = nan(1,intNeurons);
intNeuronCounter = 0;

while any(vecNeuronsLeftD2)
	%select neurons
	intNeuronCounter = intNeuronCounter + 1;
	if intNeuronCounter <= intNeurons
		vecRunNeurons = find(vecNeuronsLeftD2);
		vecThisD2 = nan(1,numel(vecRunNeurons));
		vecThisPredA = nan(1,numel(vecRunNeurons));
		vecOtherNeurons = find(vecNeuronsUsedD2);
		intRunN = numel(vecRunNeurons);
		parfor intCounterN=1:intRunN
			intNeuron = vecRunNeurons(intCounterN);
			[dblPredA,matPredA,dblDprimeSquaredOffDiagonal] = getSeparation(matModelResp([vecOtherNeurons intNeuron],:),vecTrialOriIdx);
			vecThisD2(intCounterN) = dblDprimeSquaredOffDiagonal;
			vecThisPredA(intCounterN) = dblPredA;
			%if mod(intCounterN,1000) == 0,fprintf('Processing neuron %d (%d/%d) [%s]\n',intNeuron,intCounterN,intRunN,getTime);end
		end
		[dblD2,intN] = max(vecThisD2);
		intThisNeuron = vecRunNeurons(intN);
		vecNeuronsLeftD2(intThisNeuron) = false;
		vecNeuronsUsedD2(intThisNeuron) = true;
		vecDprimeSquared(intThisNeuron) = dblD2;
		vecPredictedAcc(intThisNeuron) = vecThisPredA(intN);
		vecNeuronOrderD2(intThisNeuron) = intNeuronCounter;
		
		%message
		fprintf('Processed size %d/%d; neuron %d had largest contribution with %.3f; current group: [%s] [%s]\n',sum(vecNeuronsUsedD2),numel(vecNeuronsUsedD2),intThisNeuron,dblD2,num2str(find(vecNeuronsUsedD2)),getTime);
	else
		vecNeuronsLeftD2 = false;
	end
end

%% greedy decoder
fprintf('Starting final run [%s]\n',getTime);
[vecSortedDprime,vecSortedRunNeurons]=sort(vecDprimeSquared,'descend');
vecPerformance = nan(size(vecDprimeSquared));
parfor intGroupSize=1:length(vecSortedRunNeurons)
	vecNeurons = vecSortedRunNeurons(1:intGroupSize);
	
	%msg
	vecPerformance(intGroupSize) = doCrossValidatedDecodingML(matModelResp(vecNeurons,:),vecTrialOriIdx,1);
	fprintf('Processed group %d/%d, performance was: %.3f (%s)\n',intGroupSize,intNeurons,vecPerformance(intGroupSize),getTime);
end


%save
sGreedySim = struct;
sGreedySim.vecNeuronsLeftD2 = vecNeuronsLeftD2;
sGreedySim.vecNeuronsUsedD2 = vecNeuronsUsedD2;
sGreedySim.vecDprimeSquared = vecDprimeSquared;
sGreedySim.vecNeuronOrderD2 = vecNeuronOrderD2;
sGreedySim.vecPerformance = vecPerformance;
sGreedySim.vecPredictedAcc = vecPredictedAcc;

strDir = 'D:\Data\Results\V1_LIFmodel\';
save([strDir strSimulation '_GreedyDecoder' getDate],'sGreedySim');

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

