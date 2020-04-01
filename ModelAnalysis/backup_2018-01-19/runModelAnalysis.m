%% initialize
%clearvars;
boolLoad = false;
boolSaveFigs = true;
if ~exist('matModelResp','var')
	clearvars -except bool*;%close all;
	intLoadSim = 42; %% SET WHICH SIM TO LOAD HERE; 41 vs 42; 2 vs. 5
	boolLoad = true;
end
if intLoadSim == 11 && boolLoad
	
elseif intLoadSim == 42 && boolLoad
	strSimulation = 'xAreaDistributed_Ori2NewConn_2017-09-18'; %new connectivity
	
end
%% RUN: #header
strFigDir = 'D:\Data\Results\V1_LIFmodel\figures\';
if boolLoad
	runModelHeader;
end


%% #continue
intTrials = length(vecTrialStimType);
intStimTypes = numel(unique(vecTrialStimType));
cellSelect = cell(1,intStimTypes);
for intStimType=1:intStimTypes
	cellSelect{intStimType} = vecTrialStimType==intStimType;
end
%{
structOut = calcStimCorrsRespMat(matModelResp,cellSelect);

figure
subplot(2,2,1)
histx(structOut.matNoiseCorrs(:))

subplot(2,2,2)
imagesc(structOut.matNoiseCorrs,[-1 1]);colormap(redblue);colorbar
%}

%% fano factor
matModelResp = double(matModelResp);
dblBin = 0.2;
vecBins = (dblBin/2):dblBin:4;
figure
subplot(2,3,1)
vecFanos = var(matModelResp,0,2)./mean(matModelResp,2);
hist(vecFanos,vecBins)
xlim([0 4]);
xlabel('Fano Factor')
ylabel('Number of neurons')
title('All neurons')

subplot(2,3,2)
hist(vecFanos(vecCellTypes==1),vecBins)
xlim([0 4]);
xlabel('Fano Factor')
ylabel('Number of pyramidal neurons')
title('Only excitatory')

subplot(2,3,3)
hist(vecFanos(vecCellTypes==2),vecBins)
xlim([0 4]);
xlabel('Fano Factor')
ylabel('Number of inhibitory neurons')
title('Only inhibitory')

subplot(2,3,4)
scatter(log(xmean(matModelResp,2)),log(xstd(matModelResp,2)));
title('All neurons; Each point is one neuron')
xlabel('log(mean)')
ylabel('log(std)')

subplot(2,3,5)
scatter(log(xmean(matModelResp(vecCellTypes==1,:),2)),log(xstd(matModelResp(vecCellTypes==1,:),2)));
title('Only excitatory; Each point is one neuron')
xlabel('log(mean)')
ylabel('log(std)')

subplot(2,3,6)
scatter(log(xmean(matModelResp(vecCellTypes==2,:),2)),log(xstd(matModelResp(vecCellTypes==2,:),2)));
title('Only inhibitory; Each point is one neuron')
xlabel('log(mean)')
ylabel('log(std)')

if boolSaveFigs
	export_fig(['D:\Data\Results\V1_LIFmodel\figures\' strSimulation '_FanoFactors.tif']);
	export_fig(['D:\Data\Results\V1_LIFmodel\figures\' strSimulation '_FanoFactors.pdf']);
end

%% fano factors on shorter time scales
dblBinSize = 0.5/3; %50 ms bins
dblStimDur = vecStimStopSecs(1)-vecStimStartSecs(1);
intBins = dblStimDur/dblBinSize;
intNeurons = numel(cellSpikeTimesCortex);

matSpikeCountsTriEpoch = nan(intNeurons,intTrials,intBins);
parfor intTrial=1:intTrials
	dblTrialStartT = vecStimStartSecs(intTrial);
	for intBin=1:intBins
		dblStartT = dblTrialStartT+(intBin-1)*dblBinSize;
		dblStopT = dblStartT+dblBinSize;
		vecCellSize = size(cellSpikeTimesCortex);
		matSpikeCountsTriEpoch(:,intTrial,intBin) = cellfun(@sum,...
			cellfun(@and,...
				cellfun(@gt,cellSpikeTimesCortex,cellfill(dblStartT,vecCellSize),'UniformOutput',false),...
				cellfun(@lt,cellSpikeTimesCortex,cellfill(dblStopT,vecCellSize),'UniformOutput',false),...
			'UniformOutput',false)...
		);
	
	end
	if mod(intTrial,10) == 0
		fprintf('Now at trial %d/%d [%s]\n',intTrial,intTrials,getTime);
	end
end


%%
intOris = numel(vecTrialOris);
matMeanResp = nan(sum(vecCellTypes==1),intOris,intBins);
matVarResp = nan(sum(vecCellTypes==1),intOris,intBins);
for intOriIdx=1:intOris
	matMeanResp(:,intOriIdx,:) = xmean(matSpikeCountsTriEpoch(vecCellTypes==1,vecTrialOriIdx==intOriIdx,:),2);
	matVarResp(:,intOriIdx,:) = xvar(matSpikeCountsTriEpoch(vecCellTypes==1,vecTrialOriIdx==intOriIdx,:),2);
end
figure
for i=1:intBins
	matMu = matMeanResp(:,:,i);
	matVar = (matVarResp(:,:,i));
	
	dblPlotBinSize = 0.1;
	vecBinX = 0:dblPlotBinSize:6;
	vecBinY = 0:dblPlotBinSize:2;
	vecMuBins = min(max(round((matMu(:)-vecBinX(1))/dblPlotBinSize),1),numel(vecBinX));
	vecVarBins = min(max(round((matVar(:)-vecBinY(1))/dblPlotBinSize),1),numel(vecBinY));
	
	
	matProbDens = getFillGrid(zeros(numel(vecBinY),numel(vecBinX)),vecVarBins,vecMuBins);
	
	subplot(2,3,i)
	imagesc(vecBinX,vecBinY,log(matProbDens));
	axis xy
	xlabel('Mean (spike count)')
	ylabel('Var (spike count)')
	title(sprintf('Time after onset: %.3fs - %.3fs' ,dblBinSize*(i-1),dblBinSize*(i)))
end

% make pairwise comparisons epochs
matMu1 = matMeanResp(:,:,1);
matVar1 = matVarResp(:,:,1);
vecFano1 = matVar1(:)./matMu1(:);

matMu2 = matMeanResp(:,:,2);
matVar2 = matVarResp(:,:,2);
vecFano2 = matVar2(:)./matMu2(:);

matMu3 = matMeanResp(:,:,3);
matVar3 = matVarResp(:,:,3);
vecFano3 = matVar3(:)./matMu3(:);


vecLim = [0 2];
subplot(2,3,4)
plot(vecLim,vecLim,'k--')
hold on
scatter(vecFano1,vecFano2)
hold off
xlim(vecLim)
ylim(vecLim)
xlabel('Fano factor 0.000s - 0.167s');
ylabel('Fano factor 0.167s - 0.333s');

subplot(2,3,5)
plot(vecLim,vecLim,'k--')
hold on
scatter(vecFano2,vecFano3)
hold off
xlim(vecLim)
ylim(vecLim)
xlabel('Fano factor 0.167s - 0.333s');
ylabel('Fano factor 0.333s - 0.500s');

subplot(2,3,6)
plot(vecLim,vecLim,'k--')
hold on
scatter(vecFano1,vecFano3)
hold off
xlim(vecLim)
ylim(vecLim)
xlabel('Fano factor 0.000s - 0.167s');
ylabel('Fano factor 0.333s - 0.500s');
drawnow;
if boolSaveFigs
export_fig(['D:\Results\V1SubPopulations\simulations\' strSimulation '_FanoFactors2.tif']);
export_fig(['D:\Results\V1SubPopulations\simulations\' strSimulation '_FanoFactors2.pdf']);
end

%% select subset of trials as also presented during experiment
if max(vecTrialOris) < (2*pi)
	vecTrialOriDegs = rad2ang(vecTrialOris);
else
	vecTrialOriDegs = vecTrialOris;
end
vecOrisExperiment = [0:45:(180-1)];
intOrisSub = numel(vecOrisExperiment);
indSelect = ismember(vecTrialOriDegs,vecOrisExperiment);
if all(~indSelect)
	vecOrisExperiment = [42.5 47.5];
	intOrisSub = numel(vecOrisExperiment);
	indSelect = ismember(vecTrialOriDegs,vecOrisExperiment);
end
vecTrialOrisSub = vecTrialOris(indSelect);
matRespSub = matModelResp(:,indSelect);
[dblPerformanceML,vecDecodedIndexML,matPosteriorProbabilityML,dblMeanErrorDegsML,matConfusionML] = doCrossValidatedDecodingML(matRespSub',vecTrialOrisSub,1);

%% select subset of repetitions
%{
intUseReps = 8;
intNumRepsOrig = numel(vecTrialOrisSub)/numel(unique(vecTrialOrisSub));
vecUseTrials = nan(1,intUseReps*intOrisSub);
for intOri=1:intOrisSub
	vecTempSelect = randperm(intNumRepsOrig,intUseReps);
	vecTempThisOri = find(vecTrialOrisSub == ang2rad(vecOrisExperiment(intOri)));
	vecKeep = vecTempThisOri(vecTempSelect);
	intStart = sum(~isnan(vecUseTrials));
	vecUseTrials((intStart+1):(intStart+intUseReps)) = vecKeep;
end

%get subset of subset
vecTheseOris = vecTrialOrisSub(vecUseTrials);
matTheseResps = matRespSub(:,vecUseTrials);

[dblPerformanceML,vecDecodedIndexML,matPosteriorProbabilityML,dblMeanErrorDegsML,matConfusionML] = doCrossValidatedDecodingML(matTheseResps',vecTheseOris,1);
%}

%% subsample neurons
intIters=1000;
vecUseNeurons = 2.^(0:10);
matPerformance = nan(numel(vecUseNeurons),intIters);
for intUseNeuronIdx=1:numel(vecUseNeurons)
	intUseNeurons = vecUseNeurons(intUseNeuronIdx);
	fprintf('Running neuron group size %d/%d (%d) [%s]\n',intUseNeuronIdx,numel(vecUseNeurons),intUseNeurons,getTime);
	for intIter=1:intIters
		% select subset of repetitions
		intUseReps = 8;
		intNumRepsOrig = numel(vecTrialOrisSub)/intOrisSub;
		vecUseTrials = nan(1,intUseReps*intOrisSub);
		for intOri=1:intOrisSub
			vecTempSelect = randperm(intNumRepsOrig,intUseReps);
			vecTempThisOri = find(vecTrialOrisSub == ang2rad(vecOrisExperiment(intOri)));
			vecKeep = vecTempThisOri(vecTempSelect);
			intStart = sum(~isnan(vecUseTrials));
			vecUseTrials((intStart+1):(intStart+intUseReps)) = vecKeep;
		end
		
		%get subset of subset
		vecTheseOris = vecTrialOrisSub(vecUseTrials);
		matTheseResps = matRespSub(:,vecUseTrials);
		
		%get random neurons
		vecSelectNeurons = randperm(intNeurons,intUseNeurons);
		[dblPerformanceML] = doCrossValidatedDecodingML(matTheseResps(vecSelectNeurons,:),vecTheseOris,1);
		matPerformance(intUseNeuronIdx,intIter) = dblPerformanceML;
	end
end

%%
subplot(2,2,2)
vecLimX = [min(log2(vecUseNeurons))-0.1 max(log2(vecUseNeurons))+0.1];
plot(vecLimX,(1/intOrisSub)*[1 1],'k--');
hold on
errorbar(log2(vecUseNeurons),mean(matPerformance,2),std(matPerformance,[],2)/sqrt(intIters))
set(gca,'xtick',log2(vecUseNeurons),'xticklabel',vecUseNeurons)
ylim([0 1]);
hold off
xlim(vecLimX);
xlabel('Number of neurons used in decoding')
ylabel('Orientation decoding accuracy')
title('Model')

%%
if boolSaveFigs
export_fig('Orientation_decoding_comparison_model_experiment.tif')
export_fig('Orientation_decoding_comparison_model_experiment.pdf')
end

%% decode inhibitory/excitatory
%exc
matRespExc = matRespSub(vecCellTypes==1,:);
intNeuronsExc = size(matRespExc,1);
intIters=1000;
vecUseNeurons = 2.^(0:7);
matPerformanceExc = nan(numel(vecUseNeurons),intIters);
for intUseNeuronIdx=1:numel(vecUseNeurons)
	intUseNeurons = vecUseNeurons(intUseNeuronIdx);
	fprintf('Running neuron group size %d/%d (%d) [%s]\n',intUseNeuronIdx,numel(vecUseNeurons),intUseNeurons,getTime);
	for intIter=1:intIters
		% select subset of repetitions
		intUseReps = 8;
		intNumRepsOrig = numel(vecTrialOrisSub)/intOrisSub;
		vecUseTrials = nan(1,intUseReps*intOrisSub);
		for intOri=1:intOrisSub
			vecTempSelect = randperm(intNumRepsOrig,intUseReps);
			vecTempThisOri = find(vecTrialOrisSub == ang2rad(vecOrisExperiment(intOri)));
			vecKeep = vecTempThisOri(vecTempSelect);
			intStart = sum(~isnan(vecUseTrials));
			vecUseTrials((intStart+1):(intStart+intUseReps)) = vecKeep;
		end
		
		%get subset of subset
		vecTheseOris = vecTrialOrisSub(vecUseTrials);
		matTheseResps = matRespExc(:,vecUseTrials);
		
		%get random neurons
		vecSelectNeurons = randperm(intNeuronsExc,intUseNeurons);
		[dblPerformanceML,vecDecodedIndexML,matPosteriorProbabilityML,dblMeanErrorDegsML,matConfusionML] = doCrossValidatedDecodingML(matTheseResps(vecSelectNeurons,:),vecTheseOris,1);
		matPerformanceExc(intUseNeuronIdx,intIter) = dblPerformanceML;
	end
end

%inhib
matRespInh = matRespSub(vecCellTypes==2,:);
intNeuronsInh = size(matRespInh,1);
intIters=1000;
vecUseNeurons = 2.^(0:7);
matPerformanceInh = nan(numel(vecUseNeurons),intIters);
for intUseNeuronIdx=1:numel(vecUseNeurons)
	intUseNeurons = vecUseNeurons(intUseNeuronIdx);
	fprintf('Running neuron group size %d/%d (%d) [%s]\n',intUseNeuronIdx,numel(vecUseNeurons),intUseNeurons,getTime);
	for intIter=1:intIters
		% select subset of repetitions
		intUseReps = 8;
		intNumRepsOrig = numel(vecTrialOrisSub)/intOrisSub;
		vecUseTrials = nan(1,intUseReps*intOrisSub);
		for intOri=1:intOrisSub
			vecTempSelect = randperm(intNumRepsOrig,intUseReps);
			vecTempThisOri = find(vecTrialOrisSub == ang2rad(vecOrisExperiment(intOri)));
			vecKeep = vecTempThisOri(vecTempSelect);
			intStart = sum(~isnan(vecUseTrials));
			vecUseTrials((intStart+1):(intStart+intUseReps)) = vecKeep;
		end
		
		%get subset of subset
		vecTheseOris = vecTrialOrisSub(vecUseTrials);
		matTheseResps = matRespInh(:,vecUseTrials);
		
		%get random neurons
		vecSelectNeurons = randperm(intNeuronsInh,intUseNeurons);
		[dblPerformanceML,vecDecodedIndexML,matPosteriorProbabilityML,dblMeanErrorDegsML,matConfusionML] = doCrossValidatedDecodingML(matTheseResps(vecSelectNeurons,:),vecTheseOris,1);
		matPerformanceInh(intUseNeuronIdx,intIter) = dblPerformanceML;
	end
end

%plot
figure
subplot(2,2,1)
vecLimX = [min(log2(vecUseNeurons))-0.1 max(log2(vecUseNeurons))+0.1];
plot(vecLimX,(1/intOrisSub)*[1 1],'k--');
hold on
errorbar(log2(vecUseNeurons),mean(matPerformanceInh,2),std(matPerformanceInh,[],2)/sqrt(intIters),'r')
errorbar(log2(vecUseNeurons),mean(matPerformanceExc,2),std(matPerformanceExc,[],2)/sqrt(intIters),'b')
set(gca,'xtick',log2(vecUseNeurons),'xticklabel',vecUseNeurons)
ylim([0 1]);
hold off
xlim(vecLimX);
xlabel('Number of neurons used in decoding')
ylabel('Orientation decoding accuracy')
title('Inhibitory (red) vs Excitatory (blue); 4 oris, 8 reps')
legend({'Chance','Inhibitory','Excitatory'});
grid on

%
vecOris = unique(vecTrialOris);
intOris = numel(vecOris);
%exc
matRespExcFull = matModelResp(vecCellTypes==1,:);
intNeuronsExc = size(matRespExcFull,1);
intIters=100;
vecUseNeurons = 2.^(0:7);
matPerformanceExcFull = nan(numel(vecUseNeurons),intIters);
for intUseNeuronIdx=1:numel(vecUseNeurons)
	intUseNeurons = vecUseNeurons(intUseNeuronIdx);
	fprintf('Running neuron group size %d/%d (%d) [%s]\n',intUseNeuronIdx,numel(vecUseNeurons),intUseNeurons,getTime);
	for intIter=1:intIters
		% select subset of repetitions
		intUseReps = 8;
		intNumRepsOrig = numel(vecTrialOris)/intOris;
		vecUseTrials = nan(1,intUseReps*intOris);
		for intOri=1:intOris
			vecTempSelect = randperm(intNumRepsOrig,intUseReps);
			vecTempThisOri = find(vecTrialOris == vecOris(intOri));
			vecKeep = vecTempThisOri(vecTempSelect);
			intStart = sum(~isnan(vecUseTrials));
			vecUseTrials((intStart+1):(intStart+intUseReps)) = vecKeep;
		end
		
		%get subset of subset
		vecTheseOris = vecTrialOris(vecUseTrials);
		
		%get random neurons
		vecSelectNeurons = randperm(intNeuronsExc,intUseNeurons);
		[dblPerformanceML,vecDecodedIndexML,matPosteriorProbabilityML,dblMeanErrorDegsML,matConfusionML_Exc] = doCrossValidatedDecodingML(matRespExcFull(vecSelectNeurons,vecUseTrials),vecTheseOris,1);
		matPerformanceExcFull(intUseNeuronIdx,intIter) = dblPerformanceML;
	end
end

%inh
matRespInhFull = matModelResp(vecCellTypes==2,:);
intNeuronsInh = size(matRespInhFull,1);
intIters=100;
vecUseNeurons = 2.^(0:7);
matPerformanceInhFull = nan(numel(vecUseNeurons),intIters);
for intUseNeuronIdx=1:numel(vecUseNeurons)
	intUseNeurons = vecUseNeurons(intUseNeuronIdx);
	fprintf('Running neuron group size %d/%d (%d) [%s]\n',intUseNeuronIdx,numel(vecUseNeurons),intUseNeurons,getTime);
	for intIter=1:intIters
		% select subset of repetitions
		intUseReps = 8;
		intNumRepsOrig = numel(vecTrialOris)/intOris;
		vecUseTrials = nan(1,intUseReps*intOris);
		for intOri=1:intOris
			vecTempSelect = randperm(intNumRepsOrig,intUseReps);
			vecTempThisOri = find(vecTrialOris == vecOris(intOri));
			vecKeep = vecTempThisOri(vecTempSelect);
			intStart = sum(~isnan(vecUseTrials));
			vecUseTrials((intStart+1):(intStart+intUseReps)) = vecKeep;
		end
		
		%get subset of subset
		vecTheseOris = vecTrialOris(vecUseTrials);
		
		%get random neurons
		vecSelectNeurons = randperm(intNeuronsInh,intUseNeurons);
		[dblPerformanceML,vecDecodedIndexML,matPosteriorProbabilityML,dblMeanErrorDegsML,matConfusionML_Inh] = doCrossValidatedDecodingML(matRespInhFull(vecSelectNeurons,vecUseTrials),vecTheseOris,1);
		matPerformanceInhFull(intUseNeuronIdx,intIter) = dblPerformanceML;
	end
end

%plot
subplot(2,2,2)
vecLimX = [min(log2(vecUseNeurons))-0.1 max(log2(vecUseNeurons))+0.1];
plot(vecLimX,(1/intOris)*[1 1],'k--');
hold on
errorbar(log2(vecUseNeurons),mean(matPerformanceInhFull,2),std(matPerformanceInhFull,[],2)/sqrt(intIters),'r')
errorbar(log2(vecUseNeurons),mean(matPerformanceExcFull,2),std(matPerformanceExcFull,[],2)/sqrt(intIters),'b')
set(gca,'xtick',log2(vecUseNeurons),'xticklabel',vecUseNeurons)
ylim([0 1]);
hold off
xlim(vecLimX);
xlabel('Number of neurons used in decoding')
ylabel('Orientation decoding accuracy')
title('Inhibitory (red) vs Excitatory (blue); 12 oris, 8 reps')
legend({'Chance','Inhibitory','Excitatory'});
grid on

%%
vecOris = unique(vecTrialOris);
intOris = numel(vecOris);
%exc
matRespExcFull = matModelResp(vecCellTypes==1,:);
intNeuronsExc = size(matRespExcFull,1);
intIters=10;
vecUseNeurons = 2.^(0:7);
matPerformanceExcFullMax = nan(numel(vecUseNeurons),intIters);
for intUseNeuronIdx=1:numel(vecUseNeurons)
	intUseNeurons = vecUseNeurons(intUseNeuronIdx);
	fprintf('Running neuron group size %d/%d (%d) [%s]\n',intUseNeuronIdx,numel(vecUseNeurons),intUseNeurons,getTime);
	for intIter=1:intIters
		%get random neurons
		vecSelectNeurons = randperm(intNeuronsExc,intUseNeurons);
		[dblPerformanceML,vecDecodedIndexML,matPosteriorProbabilityML,dblMeanErrorDegsML,matConfusionML_Exc] = doCrossValidatedDecodingML(matRespExcFull(vecSelectNeurons,vecUseTrials),vecTheseOris,1);
		matPerformanceExcFullMax(intUseNeuronIdx,intIter) = dblPerformanceML;
	end
end

%inh
matRespInhFull = matModelResp(vecCellTypes==2,:);
intNeuronsInh = size(matRespInhFull,1);
intIters=10;
vecUseNeurons = 2.^(0:7);
matPerformanceInhFullMax = nan(numel(vecUseNeurons),intIters);
for intUseNeuronIdx=1:numel(vecUseNeurons)
	intUseNeurons = vecUseNeurons(intUseNeuronIdx);
	fprintf('Running neuron group size %d/%d (%d) [%s]\n',intUseNeuronIdx,numel(vecUseNeurons),intUseNeurons,getTime);
	for intIter=1:intIters
		%get random neurons
		vecSelectNeurons = randperm(intNeuronsInh,intUseNeurons);
		[dblPerformanceML,vecDecodedIndexML,matPosteriorProbabilityML,dblMeanErrorDegsML,matConfusionML_Inh] = doCrossValidatedDecodingML(matRespInhFull(vecSelectNeurons,:),vecTrialOris,1);
		matPerformanceInhFullMax(intUseNeuronIdx,intIter) = dblPerformanceML;
	end
end

%plot
subplot(2,2,3)
vecLimX = [min(log2(vecUseNeurons))-0.1 max(log2(vecUseNeurons))+0.1];
plot(vecLimX,(1/intOris)*[1 1],'k--');
hold on
errorbar(log2(vecUseNeurons),mean(matPerformanceInhFullMax,2),std(matPerformanceInhFullMax,[],2)/sqrt(intIters),'r')
errorbar(log2(vecUseNeurons),mean(matPerformanceExcFullMax,2),std(matPerformanceExcFullMax,[],2)/sqrt(intIters),'b')
set(gca,'xtick',log2(vecUseNeurons),'xticklabel',vecUseNeurons)
ylim([0 1]);
hold off
xlim(vecLimX);
xlabel('Number of neurons used in decoding')
ylabel('Orientation decoding accuracy')
title('Inhibitory (red) vs Excitatory (blue); 12 oris, 100 reps')
legend({'Chance','Inhibitory','Excitatory'});
grid on
drawnow

%%
vecOris = unique(vecTrialOris);
intOris = numel(vecOris);
%exc
matRespExcFull = matModelResp(vecCellTypes==1,:);
intNeuronsExc = size(matRespExcFull,1);
intIters=100;
intUseNeurons = 4;
vecUseRepetitions = 2.^(2:6);
matPerformanceExcReps = nan(numel(vecUseRepetitions),intIters);
for intUseRepIdx=1:numel(vecUseRepetitions)
	intUseReps = vecUseRepetitions(intUseRepIdx);
	fprintf('Running repetition size %d/%d (%d) [%s]\n',intUseRepIdx,numel(vecUseRepetitions),intUseReps,getTime);
	for intIter=1:intIters
		% select subset of repetitions
		intNumRepsOrig = numel(vecTrialOris)/intOris;
		vecUseTrials = nan(1,intUseReps*intOris);
		for intOri=1:intOris
			vecTempSelect = randperm(intNumRepsOrig,intUseReps);
			vecTempThisOri = find(vecTrialOris == vecOris(intOri));
			vecKeep = vecTempThisOri(vecTempSelect);
			intStart = sum(~isnan(vecUseTrials));
			vecUseTrials((intStart+1):(intStart+intUseReps)) = vecKeep;
		end
		
		%get subset of subset
		vecTheseOris = vecTrialOris(vecUseTrials);
		
		%get random neurons
		vecSelectNeurons = randperm(intNeuronsExc,intUseNeurons);
		[dblPerformanceML,vecDecodedIndexML,matPosteriorProbabilityML,dblMeanErrorDegsML,matConfusionML_Exc] = doCrossValidatedDecodingML(matRespExcFull(vecSelectNeurons,vecUseTrials),vecTheseOris,1);
		matPerformanceExcReps(intUseRepIdx,intIter) = dblPerformanceML;
	end
end


%inh
matRespInhFull = matModelResp(vecCellTypes==2,:);
intNeuronsInh = size(matRespInhFull,1);
intIters=100;
intUseNeurons = 4;
vecUseRepetitions = 2.^(2:6);
matPerformanceInhReps = nan(numel(vecUseRepetitions),intIters);
for intUseRepIdx=1:numel(vecUseRepetitions)
	intUseReps = vecUseRepetitions(intUseRepIdx);
	fprintf('Running repetition size %d/%d (%d) [%s]\n',intUseRepIdx,numel(vecUseRepetitions),intUseReps,getTime);
	for intIter=1:intIters
		% select subset of repetitions
		intNumRepsOrig = numel(vecTrialOris)/intOris;
		vecUseTrials = nan(1,intUseReps*intOris);
		for intOri=1:intOris
			vecTempSelect = randperm(intNumRepsOrig,intUseReps);
			vecTempThisOri = find(vecTrialOris == vecOris(intOri));
			vecKeep = vecTempThisOri(vecTempSelect);
			intStart = sum(~isnan(vecUseTrials));
			vecUseTrials((intStart+1):(intStart+intUseReps)) = vecKeep;
		end
		
		%get subset of subset
		vecTheseOris = vecTrialOris(vecUseTrials);
		
		%get random neurons
		vecSelectNeurons = randperm(intNeuronsInh,intUseNeurons);
		[dblPerformanceML,vecDecodedIndexML,matPosteriorProbabilityML,dblMeanErrorDegsML,matConfusionML_Inh] = doCrossValidatedDecodingML(matRespInhFull(vecSelectNeurons,vecUseTrials),vecTheseOris,1);
		matPerformanceInhReps(intUseRepIdx,intIter) = dblPerformanceML;
	end
end
%%
%plot
subplot(2,2,4)
cla reset
vecLimX = [min(log2(vecUseRepetitions))-0.1 max(log2(vecUseRepetitions))+0.1];
%plot(vecLimX,(1/intOris)*[1 1],'k--');
hold on
errorbar(log2(vecUseRepetitions),mean(matPerformanceInhReps,2),std(matPerformanceInhReps,[],2)/sqrt(intIters),'r')
errorbar(log2(vecUseRepetitions),mean(matPerformanceExcReps,2),std(matPerformanceExcReps,[],2)/sqrt(intIters),'b')
set(gca,'xtick',log2(vecUseRepetitions),'xticklabel',vecUseRepetitions)
%ylim([0 1]);
hold off
xlim(vecLimX);
xlabel('Number of repetitions used in decoding')
ylabel('Orientation decoding accuracy')
title('Inhibitory (red) vs Excitatory (blue); 12 oris, 4 neurons')
legend({'Chance','Inhibitory','Excitatory'});
grid on
drawnow

%%
if boolSaveFigs
strDir = 'D:\Data\Results\V1_LIFmodel\';
export_fig([strDir strSimulation '_DecodingInhExc.tif']);
export_fig([strDir strSimulation '_DecodingInhExc.pdf']);
end
