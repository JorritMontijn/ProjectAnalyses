%% description
%{
Sliding window analysis
B1(B1)->B2 ... => B2(B1)->B2
%}
%% clear
clear all;

%% set variables
%On which timestamp to align as t=0
strAlignOn = 'Change';
intFullIters = 10;
intResamplings = 10;
boolShuffle = false;
dblLambda = 1/10;
boolSavePlots = true;
strFigDir = 'D:\Data\Results\BumpsMatthijs\';

% parameters for window of interest:
dblStartBaseT = -0.4;
dblStopBaseT = -0.2;
dblStartEp1T = 0;
dblStopEp1T = 0.2;
dblStartEp2T = 0.2;
dblStopEp2T = 0.4;

%get neuronal activity from t=0 - t=0.4
dblBinSizeSecs = 100/1000;%25/1000;
vecBinEdges = [0 dblBinSizeSecs];
dblBinOffsetStep = (dblBinSizeSecs/10);
vecBinOffsets = -2:(dblBinSizeSecs/10):2;
vecBinStartT = vecBinOffsets;
intBins = numel(vecBinStartT);
vecRemoveSessions = [3 6 9]; %due to bad performance

%% load data
fprintf('Loading data... [%s]\n',getTime);
[Data] = MOL_GetData('D:\Data\Processed\MatthijsOudeLohuis\','CHDET',{'ChangeDetectionConflict'},{'2003' '2009' '2010' '2011' '2012' '2013'},[],{'sessionData' 'trialData' 'spikeData'});
sessionData     = Data.sessionData;
trialData       = Data.trialData;
spikeData       = Data.spikeData;

%% Remove last 20 trials:
trialData = MOL_RemoveLastnTrials(trialData,20);

%% Filter good neurons:
spikeData = MOL_filterNeurons(sessionData,trialData,spikeData);

%% split data by session
[cellSessions,ia,vecSessionIdx] = unique(spikeData.session_ID);

%% split by area
[cellAreas,ia,vecAreaIdx] = unique(spikeData.area);

%% Get significant first and second bump responses:
%[sessionData,trialData,spikeData] = MOL_calc_signresp(sessionData,trialData,spikeData);

%% make populations
[cellPopID,ia,vecPopulationIdx] = unique(vecSessionIdx*10 + vecAreaIdx);

%% define V1 populations
indV1=vecAreaIdx==find(strcmp(cellAreas,'V1'));
vecPopsV1 = unique(vecPopulationIdx(indV1));
intPops = numel(vecPopsV1);

%% pre-allocate output
vecCellsPerPop = nan(1,intPops);
matSubspaceChange = nan(intBins,intPops,intFullIters);

%save output
cellAggMA = cell(intPops,intFullIters);
cellAggMM = cell(intPops,intFullIters);

%  analyze
dblTimeConversion = 10^6;
for intPopIdx=1:intPops
	%% prepare, retrieve and transform data
	%retrieve cells
	intPopulation = vecPopsV1(intPopIdx);
	indUseCells = vecPopulationIdx==intPopulation;
	intCells = sum(indUseCells);
	intCutOffCellNr = 16;
	vecCellsPerPop(intPopIdx) = intCells;
	if intCells < intCutOffCellNr
		fprintf('Pop %d, Number of cells is %d, which is under %d, skipping... [%s]\n',intPopIdx,intCells,intCutOffCellNr,getTime);
		continue;
	end
	
	%get spikes and transform to seconds
	cellSpikes = cellfun(@rdivide,spikeData.ts(indUseCells),cellfill(dblTimeConversion,[intCells 1]),'uniformoutput',false);
	strSesID = unique(spikeData.session_ID(indUseCells));
	if numel(strSesID)~=1,error([mfilename ':SessionNotUnique'],'Something went wrong with population selection...');end
	
	%get trial data
	indSesTrials = strcmp(trialData.session_ID,strSesID);
	intSesTrials = sum(indSesTrials);
	
	vecTrialType = trialData.trialType(indSesTrials);
	vecTrialVisOriPre = trialData.visualOriPreChange(indSesTrials);
	vecTrialVisOriChange = trialData.visualOriChange(indSesTrials);
	vecTrialAudFreqPre = trialData.audioFreqPreChange(indSesTrials);
	vecTrialAudFreqChange = trialData.audioFreqChange(indSesTrials);
	vecTrialCorrect = trialData.correctResponse(indSesTrials);
	vecTrialRespSide = trialData.responseSide(indSesTrials);
	
	vecTrialLickSecs = cellfun(@rdivide,trialData.lickTime(indSesTrials),cellfill(dblTimeConversion,[intSesTrials 1]),'uniformoutput',false);
	
	vecTrialStartSecs = trialData.trialStart(indSesTrials)./dblTimeConversion;
	vecTrialChangeSecs = trialData.stimChange(indSesTrials)./dblTimeConversion;
	vecTrialEndSecs = trialData.trialEnd(indSesTrials)./dblTimeConversion;
	vecTrialStartITISecs = trialData.itiStart(indSesTrials)./dblTimeConversion;
	vecTrialEndITISecs = trialData.itiEnd(indSesTrials)./dblTimeConversion;
	vecTrialRespSecs = trialData.responseLatency(indSesTrials)./dblTimeConversion;
	
	%sub-select trials that are large change, visual, and correct
	%indKeepTrials = abs(vecTrialVisOriChange) == 90;
	indKeepTrials = abs(vecTrialAudFreqChange) == 0 & abs(vecTrialVisOriChange) > 0;
	%indKeepTrials = abs(vecTrialVisOriChange) > 0;
	dblCutOff = 0.4;
	dblFracCorr = sum(vecTrialCorrect(indKeepTrials))/sum(indKeepTrials);
	if ismember(intPopIdx,vecRemoveSessions)
		fprintf('Pop %d, Correct response proportion is %.3f, which is under %.3f, skipping... [%s]\n',intPopIdx,dblFracCorr,dblCutOff,getTime);
		continue;
	end
	fprintf('Rec%d; %d vis trials, %d av trials, %d aud trials\n',intPopIdx,sum(indKeepTrials),sum(abs(vecTrialAudFreqChange) > 0 & abs(vecTrialVisOriChange) > 0),sum(abs(vecTrialAudFreqChange) > 0 & abs(vecTrialVisOriChange) == 0));
	continue
	
	%message
	fprintf('Analyzing pop %d/%d; %d neurons, %d trials [%s]\n',intPopIdx,intPops,intCells,sum(indKeepTrials),getTime);
	
	%define 1st vs 2nd bump windows
	if strcmp(strAlignOn,'Change')
		vecT0 = vecTrialChangeSecs(indKeepTrials);
	end
	%build vectors
	intTrials = numel(vecT0);
	vecStartBaseT = vecT0 + dblStartBaseT;
	vecStopBaseT = vecT0 + dblStopBaseT;
	vecStartEp1T = vecT0 + dblStartEp1T;
	vecStopEp1T = vecT0 + dblStopEp1T;
	vecStartEp2T = vecT0 + dblStartEp2T;
	vecStopEp2T = vecT0 + dblStopEp2T;
	vecEdges = sort(cat(1,vecStartBaseT,vecStopBaseT,vecStartEp1T,vecStopEp1T,vecStartEp2T,vecStopEp2T));
	vecBaseEdges = sort(cat(1,vecStartBaseT,vecStopBaseT));
	vecEp1Edges = sort(cat(1,vecStartEp1T,vecStopEp1T));
	vecEp2Edges = sort(cat(1,vecStartEp2T,vecStopEp2T));
	%calculate response matrix
	matRespBase = nan(intCells,intTrials);
	matRespEp1 = nan(intCells,intTrials);
	matRespEp2 = nan(intCells,intTrials);
	%pre-alloc spiking vectors
	cellTrialPerSpike = cell(intCells,1);
	cellTimePerSpike = cell(intCells,1);
	%bin spikes
	matActBin = zeros(intBins,intCells,intTrials);
	for intNeuron=1:intCells
		vecBaseR = histcounts(cellSpikes{intNeuron},vecBaseEdges);
		matRespBase(intNeuron,:) = vecBaseR(1:2:end);
		vecEp1R = histcounts(cellSpikes{intNeuron},vecEp1Edges);
		matRespEp1(intNeuron,:) = vecEp1R(1:2:end);
		vecEp2R = histcounts(cellSpikes{intNeuron},vecEp2Edges);
		matRespEp2(intNeuron,:) = vecEp2R(1:2:end);
		[vecTrialPerSpike,vecTimePerSpike] = getSpikesInTrial(cellSpikes{intNeuron},vecT0+vecBinOffsets(1));
		cellTrialPerSpike{intNeuron} = vecTrialPerSpike;
		cellTimePerSpike{intNeuron} = vecTimePerSpike+vecBinOffsets(1);
		%bin spikes
		for intTrial=1:intTrials
			vecSpikeTimes = vecTimePerSpike(vecTrialPerSpike==intTrial)+vecBinOffsets(1);
			for intBin=1:intBins
				matActBin(intBin,intNeuron,intTrial) = histcounts(vecSpikeTimes,vecBinEdges+vecBinOffsets(intBin));
			end
		end
	end
	
	%% prep data
	matData = cat(1,matRespEp1,matRespEp2);
	vecCellArea = cat(1,ones(intCells,1),ones(intCells,1)*2);
	vecTrialStimType = ones(1,intTrials);
	intStimTypes = numel(unique(vecTrialStimType));
	
	%% go through iters
	fprintf('    Iter (tot: %d) : \n',intFullIters);
	for intIter=1:intFullIters
		%% msg
		fprintf('\b[%d]\n',intIter);
		
		%% get data splits
		%set parameters
		dblCutOff = 0.9;
		sParams=struct;
		intUseSize = intCells;%16%floor(intCells/2);
		sParams.intSizeX = intUseSize;%floor(intCells/2);
		sParams.intSizeY = intUseSize;%floor(intCells/2);
		sParams.intResamplings = intResamplings;
		sParams.vecCellArea = vecCellArea;
		sParams.intWithinArea = 3;%0, across; 1, within 1; 2, within 2
		sParams.vecUseStimTypes = 1;
		cellStr = {'B1','B2','B1-B2'};
		strC = cellStr{sParams.intWithinArea};
		
		%get splits
		[cellMatX1,cellMatX2,cellNeuronsX,cellMatY1,cellMatY2,cellNeuronsY,cellTrials12,cellTrials1,cellTrials2] = ...
			doDimDataSplitsCV(matData,vecTrialStimType,sParams);
		
		%% analyze
		intMaxDim = sParams.intSizeY;
		[matAlignmentWithX,matMagnitudeOfPop] = doDimPrivateAnalysisCV2(matActBin,cellNeuronsX,cellNeuronsY,cellMatX1,cellMatY1,cellMatX2,cellMatY2,cellTrials1,cellTrials2,dblLambda,boolShuffle);
		matMeanAlignment = mean(mean(matAlignmentWithX,4),3);
		matMeanMagnitude = mean(mean(matMagnitudeOfPop,4),3);
		vecNormMeanAlignment = mean(matMeanAlignment,2);
		%vecSdAlignment = std(matMeanAlignment,[],2);
		
		%% save output
		matSubspaceChange(:,intPopIdx,intIter) = vecNormMeanAlignment; %1D
		cellAggMA{intPopIdx,intIter} = matMeanAlignment;
		cellAggMM{intPopIdx,intIter} = matMeanMagnitude;
	end
end
warning('on','MATLAB:nearlySingularMatrix')

%% prep data
%alignment
cellAggMA(all(cellfun(@isempty,cellAggMA),2),:) = [];
cellMeanMA = cellfun(@nanmean,cellAggMA,cellfill(2,size(cellAggMA)),'Uniformoutput',false);
cellMeanMA=cellfun(@transpose,reshape(cellMeanMA,[size(cellMeanMA,1) 1 size(cellMeanMA,2)]),'Uniformoutput',false);
matNeuralDynamics = nanmean(cell2mat(cellMeanMA),3);
matNeuralDynamics(all(isnan(matNeuralDynamics),2),:) = [];
matNeuralDynamics = mean(matNeuralDynamics,3);
%magnitude
cellAggMM(all(cellfun(@isempty,cellAggMM),2),:) = [];
cellMeanMM = cellfun(@nanmean,cellAggMM,cellfill(2,size(cellAggMM)),'Uniformoutput',false);
cellMeanMM=cellfun(@transpose,reshape(cellMeanMM,[size(cellMeanMM,1) 1 size(cellMeanMM,2)]),'Uniformoutput',false);
matNeuralMagnitude = nanmean(cell2mat(cellMeanMM),3);
matNeuralMagnitude(all(isnan(matNeuralMagnitude),2),:) = [];
matNeuralMagnitude = mean(matNeuralMagnitude,3);

%get start/stop bins
intStartBin = find(vecBinStartT>=0,1);
intStopBin = find(vecBinStartT>=0.2,1);

%mean/sd alignments
vecMeanAlignment = mean(matNeuralDynamics,1);
vecSdAlignment = std(matNeuralDynamics,[],1);
vecLinAlign = linspace(vecMeanAlignment(intStartBin),vecMeanAlignment(intStopBin),intStopBin-intStartBin+1);
intRecordings = size(matNeuralDynamics,1);

%mean/sd magnitudes
vecMeanMagnitude = mean(matNeuralMagnitude,1);
vecSdMagnitude = std(matNeuralMagnitude,[],1);

%normalize to [-1 1]
matNormDynamics = bsxfun(@minus,matNeuralDynamics,matNeuralDynamics(:,intStopBin));
matNormDynamics = bsxfun(@rdivide,matNormDynamics,matNormDynamics(:,intStartBin));
matNormDynamics = (matNormDynamics - 0.5)*2;
%mean/sd normalized alignments
vecNormMeanAlignment = mean(matNormDynamics,1);
vecNormSdAlignment = std(matNormDynamics,[],1);
vecNormLinAlign = linspace(vecNormMeanAlignment(intStartBin),vecNormMeanAlignment(intStopBin),intStopBin-intStartBin+1);


%prep plot
%close all;
figure
subplot(2,2,1)
hold on
%plot(vecBinStartT(intStartBin:intStopBin),vecLinAlign,'Color',[0.5 0.5 0.5]);
errorfill(vecBinStartT,vecMeanAlignment,vecSdAlignment/sqrt(intRecordings),'color',lines(1))
hold off
title(sprintf('A) %.1f sliding window; x-val is window start',dblBinSizeSecs))
ylabel('Alignment along B1-B2 axis');
xlabel('Time after stim change (s)');
fixfig;

%color limits
vecLimC = [min(vecBinStartT) max(vecBinStartT)];
subplot(2,2,2)
colormap(redbluepurple)
cline(vecMeanAlignment,vecMeanMagnitude,vecBinStartT);
h=colorbar;
h.Ticks = [0 0.5 1];
h.TickLabels = linspace(vecLimC(1),vecLimC(2),numel(h.Ticks));
xlabel('Alignment along B1-B2 axis');
ylabel('Activation magnitude (vector norm)');
ylabel(h, 'Time after stim change (s)')
xlim([-1.5 1.5])
title(sprintf('B) Neural space trajectory'))
fixfig;

subplot(2,2,3)
hold on
%plot(vecBinStartT(intStartBin:intStopBin),vecNormLinAlign,'Color',[0.5 0.5 0.5]);
errorfill(vecBinStartT,vecNormMeanAlignment,vecNormSdAlignment/sqrt(intRecordings),'color',lines(1))
hold off
title('C) Normalized per recording to [-1 -1]')
ylabel('Alignment along B1-B2 axis');
xlabel('Time after stim change (s)');
fixfig;

subplot(2,2,4)
colormap(redbluepurple);
cline(vecNormMeanAlignment,vecMeanMagnitude,vecBinStartT);
h=colorbar;
h.Ticks = [0 0.5 1];
h.TickLabels = linspace(vecLimC(1),vecLimC(2),numel(h.Ticks));
xlabel('Normalized B1-B2 alignment');
ylabel('Activation magnitude (vector norm)');
ylabel(h, 'Time after stim change (s)')
xlim([-1.5 1.5])
title(sprintf('D) Neural space trajectory (normalized)'))
fixfig;
maxfig();

if boolSavePlots
	strFigFile = ['DimCom3ST2_Window' sprintf('%.1f',dblBinSizeSecs) 'ProgressionNeuralDynamicsSummary' getDate];
	export_fig([strFigDir strFigFile '.tif']);
	export_fig([strFigDir strFigFile '.pdf']);
end



