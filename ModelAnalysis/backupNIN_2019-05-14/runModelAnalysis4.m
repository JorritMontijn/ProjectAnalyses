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
	strSimulation = 'xAreaSquareGrating2017-03-21';
		
end

%% RUN: #header
strFigDir = 'D:\Data\Results\V1_LIFmodel\figures\';
if boolLoad
	runModelHeader;
end


%% continue
matModelResp = matModelResp./dblStimDur;
close all
figure
imagesc(matModelResp(vecCellTypes==1,:));colorbar
imagesc(zscore(matModelResp(vecCellTypes==1,:),[],2));colorbar

%% get spike counts for small time bins
dblBinSize = dblStimDur/20;%0.5/3; %50 ms bins
intBins = dblStimDur/dblBinSize;
intTrials = numel(vecTrialOris);
intNeurons = numel(cellSpikeTimesCortex);
if ~exist('matSpikeCounts','var')
	matSpikeCounts = nan(intNeurons,intTrials,intBins);
	parfor intTrial=1:intTrials
		dblTrialStartT = vecStimStartSecs(intTrial);
		for intBin=1:intBins
			dblStartT = dblTrialStartT+(intBin-1)*dblBinSize;
			dblStopT = dblStartT+dblBinSize;
			
			matSpikeCounts(:,intTrial,intBin) = cellfun(@sum,...
				cellfun(@and,cellfun(@gt,cellSpikeTimesCortex,cellfill(dblStartT,vecCellSize),'UniformOutput',false),...
				cellfun(@lt,cellSpikeTimesCortex,cellfill(dblStopT,vecCellSize),'UniformOutput',false),...
				'UniformOutput',false));
		end
		if mod(intTrial,10) == 0
			fprintf('Now at trial %d/%d [%s]\n',intTrial,intTrials,getTime);
		end
	end
	
	%save
	save(['D:\Data\Results\V1_LIFmodel\Simulation_' strSimulation '_prepro.mat'],'sData','matModelResp','matSpikeCounts','dblBinSize','intBins');
end

%put in 2D matrix
intT=intTrials * intBins;
matActivity = nan(intNeurons,intT);
for intTrial=1:intTrials
	matActivity(:,(1+(intTrial-1)*intBins):(intTrial*intBins)) = matSpikeCounts(:,intTrial,:);
end
 [vecHeterogeneity,vecActivity] = calcMatRespHeteroGen(matActivity);

 vecCorrX=nan(1,intBins*2);
 for intBin=0:(intBins*2-1);
	 vecCorrX(intBin+1)=corr(vecActivity(1:(end-intBin)),vecActivity((intBin+1):end));
 end

 figure
 subplot(2,2,1)
 plot((0:(intBins*2-1))*dblBinSize,vecCorrX)
 
 vecSpikePerSecPerNeuron = (vecActivity/dblBinSize);
 subplot(2,2,2)
 plot(sort(vecSpikePerSecPerNeuron))
 
figure
imagesc(matActivity)


%% pca
[coeff,score,latent,tsquared,explained,mu] = pca(matModelResp');
coefforth = inv(diag(std(matModelResp')))* coeff;

figure
subplot(2,2,1)
imagesc(coeff,max(abs(coeff(:)))*[-1 1]);
colormap(redblue)
colorbar

subplot(2,2,2)
plot(cumsum(explained))

subplot(2,2,3)
imagesc(coefforth,max(abs(coefforth(:)))*[-1 1]);
colormap(redblue)
colorbar


matXCentered = score * coeff';

%select only classes 1 and 2
vecNeurons = 1:intNeurons;
matData = matModelResp(vecNeurons,:);
indClasses12 = (vecTrialOriIdx==1 | vecTrialOriIdx==2);
matData12 = matData(:,indClasses12);
vecTrialOriIdx12 = label2idx(vecTrialOriIdx(indClasses12));

matThisD1 = matData12(:,vecTrialOriIdx12==1);
matThisD2 = matData12(:,vecTrialOriIdx12==2);
intTrials = size(matThisD1,2);

[vecWeightsLogReg, vecLLH] = logitBin([matThisD1 matThisD2], [zeros(1,size(matThisD1,2)) ones(1,size(matThisD2,2))], 0.001);
vecClass1 = vecWeightsLogReg'*[matThisD1;ones(1,size(matThisD1,2))];
vecClass2 = vecWeightsLogReg'*[matThisD2;ones(1,size(matThisD2,2))];
dblDprimeLogReg = getdprime2(vecClass1,vecClass2);
y = ~round(exp(-log1pexp([vecClass1 vecClass2])));
dblPerfLogReg = sum(y==[zeros(1,numel(vecClass1)) ones(1,numel(vecClass2))])/numel(y);

%get princomps
vecStimPred = nan(1,intNeurons);
for vecSelectPCs=1:intNeurons
matDataPC = matXCentered(:,vecSelectPCs)';
matDataPC12 = matDataPC(:,indClasses12);
matThisPCD1 = matDataPC12(:,vecTrialOriIdx12==1);
matThisPCD2 = matDataPC12(:,vecTrialOriIdx12==2);

[vecWeightsLogRegPC, vecLLH] = logitBin([matThisPCD1 matThisPCD2], [zeros(1,size(matThisPCD1,2)) ones(1,size(matThisPCD2,2))], 0.001);
vecClass1PC = vecWeightsLogRegPC'*[matThisPCD1;ones(1,size(matThisPCD1,2))];
vecClass2PC = vecWeightsLogRegPC'*[matThisPCD2;ones(1,size(matThisPCD2,2))];
dblDprimeLogRegPC = getdprime2(vecClass1PC,vecClass2PC);
y = ~round(exp(-log1pexp([vecClass1PC vecClass2PC])));
dblPerfLogRegPC = sum(y==[zeros(1,numel(vecClass1PC)) ones(1,numel(vecClass2PC))])/numel(y);

vecStimPred(vecSelectPCs) = dblDprimeLogRegPC;
end
[a,b]=sort(vecStimPred,'descend')

vecSelectPCs = b(1)
matDataPC = matXCentered(:,vecSelectPCs)';
matDataPC12 = matDataPC(:,indClasses12);
matThisPCD1 = matDataPC12(:,vecTrialOriIdx12==1);
matThisPCD2 = matDataPC12(:,vecTrialOriIdx12==2);

[vecWeightsLogRegPC, vecLLH] = logitBin([matThisPCD1 matThisPCD2], [zeros(1,size(matThisPCD1,2)) ones(1,size(matThisPCD2,2))], 0.001);
vecClass1PC = vecWeightsLogRegPC'*[matThisPCD1;ones(1,size(matThisPCD1,2))];
vecClass2PC = vecWeightsLogRegPC'*[matThisPCD2;ones(1,size(matThisPCD2,2))];
dblDprimeLogRegPC = getdprime2(vecClass1PC,vecClass2PC);
y = ~round(exp(-log1pexp([vecClass1PC vecClass2PC])));
dblPerfLogRegPC = sum(y==[zeros(1,numel(vecClass1PC)) ones(1,numel(vecClass2PC))])/numel(y);




%F-prime based
vecNeurons = 1:intNeurons;
matData = matModelResp(vecNeurons,:);
matShuffledData = getShuffledData(matData',vecTrialOriIdx,2)';

indClasses12 = (vecTrialOriIdx==1 | vecTrialOriIdx==2);
matData12 = matShuffledData(:,indClasses12);
vecTrialOriIdx12 = label2idx(vecTrialOriIdx(indClasses12));

matThisD1 = matData12(:,vecTrialOriIdx12==1);
matThisD2 = matData12(:,vecTrialOriIdx12==2);
intTrials = size(matThisD1,2);

vecFprime = (xmean(matThisD1,2)-xmean(matThisD2,2));
vecSep = vecFprime' * matData12;

[n1, xout1] = histx(vecSep(vecTrialOriIdx12==1));
[n2, xout2] = histx(vecSep(vecTrialOriIdx12==2));

subplot(2,2,1)
hold on
plot(xout1,n1,'b')
plot(xout2,n2,'r')
hold off

[vecWeightsLogReg, vecLLH] = logitBin([matThisD1 matThisD2], [zeros(1,size(matThisD1,2)) ones(1,size(matThisD2,2))], 0.001);
vecClass1 = vecWeightsLogReg'*[matThisD1;ones(1,size(matThisD1,2))];
vecClass2 = vecWeightsLogReg'*[matThisD2;ones(1,size(matThisD2,2))];
dblDprimeLogReg = getdprime2(vecClass1,vecClass2);
y = ~round(exp(-log1pexp([vecClass1 vecClass2])));
dblPerfLogReg = sum(y==[zeros(1,numel(vecClass1)) ones(1,numel(vecClass2))])/numel(y);



[vecWeightsLogRegF, vecLLH] = logitBin(vecSep, vecTrialOriIdx12, 0.001);
vecClass1F = vecWeightsLogRegF'*[vecSep(vecTrialOriIdx12==1);ones(1,size(vecSep(vecTrialOriIdx12==1),2))];
vecClass2F = vecWeightsLogRegF'*[vecSep(vecTrialOriIdx12==2);ones(1,size(vecSep(vecTrialOriIdx12==2),2))];
dblDprimeLogRegF = getdprime2(vecClass1F,vecClass2F);
y = ~round(exp(-log1pexp([vecClass1F vecClass2F])));
dblPerfLogRegF = sum(y==[zeros(1,numel(vecClass1F)) ones(1,numel(vecClass2F))])/numel(y);

