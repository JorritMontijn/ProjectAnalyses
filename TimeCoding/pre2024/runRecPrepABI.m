% get matching recording data
sRec = sAggABI(intRec);
strRecOrig = sRec.Exp;
strRec = strrep(strRecOrig,'ecephys_session_','ABI');
fprintf('Preparing %d/%d: %s, area %s [%s]\n',intRec,numel(vecUseRec),strRecOrig,strArea,getTime);

strField = ['s' strRunStim];
if strcmp(strRunStim,'DG') && isfield(sRec.structStimAgg,strField)
	%% DG
	% concatenate stimulus structures
	structStim = sRec.structStimAgg.sDG; %sDG, sNM, sNS
	
	%generate stimulus vectors
	vecStimLabels = mod(structStim.orientation(:)',360);
	vecStimOnTime = structStim.startT(:)';
	vecStimOffTime = structStim.stopT(:)';
	%remove nans
	indRem = isnan(vecStimLabels);
	vecStimLabels(indRem) = [];
	vecStimOnTime(indRem) = [];
	vecStimOffTime(indRem) = [];
	
	
elseif strcmp(strRunStim,'NM') && isfield(sRec.structStimAgg,strField)
	%% NM
	% concatenate stimulus structures
	structStim = sRec.structStimAgg.sNM; %sDG, sNM, sNS
	vecOrigStartIdx = [1; 1+find(diff(structStim.frame)<0)];
	vecOrigStimOnTime = flat(structStim.startT(vecOrigStartIdx))';
	dblDur = median(diff(vecOrigStimOnTime));
	vecOrigStimOffTime = vecOrigStimOnTime+dblDur;
	intFramesInMovie = max(structStim.frame)+1;
	dblBinAvg = 10;
	intUseBinsInMovie = intFramesInMovie/dblBinAvg;
	dblBinRate = round(intUseBinsInMovie/dblDur);
	dblBinDur = 1/dblBinRate;
	vecBinEdges = linspace(0,dblBinDur*intUseBinsInMovie,intUseBinsInMovie+1);
	
	%generate fake stimulus vectors
	vecUniqueStims = 1:intUseBinsInMovie;
	vecFrameIdx = repmat(vecUniqueStims,[1 numel(vecOrigStimOnTime)]);
	vecStimLabels = (vecFrameIdx/max(vecFrameIdx))*180;
	vecStimOnTime = flat(vecBinEdges(1:(end-1))' + vecOrigStimOnTime)';
	vecStimOffTime = vecStimOnTime + dblBinDur;
end
%general prep
[vecStimIdx,vecUniqueStims,vecRepNum,cellSelect,vecTrialRepetition] = val2idx(vecStimLabels);
intRepNum = min(vecRepNum);
indRem2 = vecTrialRepetition>intRepNum;
vecStimLabels(indRem2) = [];
vecStimOnTime(indRem2) = [];
vecStimOffTime(indRem2) = [];
[vecStimIdx,vecUniqueStims,vecRepNum,cellSelect,vecTrialRepetition] = val2idx(vecStimLabels);
intTrialNum = numel(vecStimOnTime);
intStimNum = numel(vecUniqueStims);

%% remove untuned cells
%get data matrix
indArea1Neurons = contains(sRec.cellAreas,strArea,'IgnoreCase',true);
intNeuronsInArea = sum(indArea1Neurons);
if intNeuronsInArea==0,return;end
cellSpikeTimes = sRec.cellSpikes(indArea1Neurons);
dblDur = median(vecStimOffTime-vecStimOnTime);
matData = getSpikeCounts(cellSpikeTimes,vecStimOnTime,dblDur);
[matRespNSR,vecStimTypes,vecUnique] = getStimulusResponses(matData,vecStimIdx);
matAvgR = mean(matRespNSR,3);
matTrialR = squeeze(mean(matRespNSR,2));
vecP_A = nan(1,size(matData,1));
vecP_Z = nan(1,size(matData,1));
for intN=1:size(matData,1)
	vecP_A(intN) = anova1(matData(intN,:),vecStimIdx,'off');
	vecP_Z(intN) = zetatest(cellSpikeTimes{intN},vecStimOnTime,[],[],0);
end
vecP = min(vecP_Z,vecP_A);

%remove untuned cells
sOut = getTuningCurves(matData,vecStimLabels,0);
dblMinCount = 100;
indTuned = vecP<0.05;
indResp = sum(matData,2)'>dblMinCount;

%prep
cellSpikeTimes(~indResp)=[];
indTuned(~indResp)=[];
intRespN = sum(indResp);

dblStimDur = roundi(median(vecStimOffTime - vecStimOnTime),1,'ceil');
dblPreTime = -dblStartT;%0.3;
dblPostTime = 0;%0.3;
dblMaxDur = dblStimDur+dblPreTime+dblPostTime;
dblBinWidth = 0.05;
vecBinEdges = 0:dblBinWidth:dblMaxDur;
vecStimTime = vecBinEdges(2:end)-dblBinWidth/2 - dblPreTime;
indStimBins = vecStimTime > 0 & vecStimTime < dblStimDur;
intBinNum = numel(vecBinEdges)-1;
matBNSR = nan(intBinNum,intRespN,intStimNum,intRepNum);
matBNT = nan(intBinNum,intRespN,intTrialNum);
%matBNT_shifted = nan(intBinNum,intRespN,intTrialNum);

%% check non-stationarity
vecRepCounter = zeros(1,intStimNum);
%get spikes per trial per neuron
cellSpikeTimesPerCellPerTrial = cell(intRespN,intTrialNum);
vecNonStat = nan(1,intRespN);
vecViolIdx1ms = nan(1,intRespN);
vecViolIdx2ms = nan(1,intRespN);
boolDiscardEdges = true;
for intN=1:intRespN
	% build pseudo data, stitching stimulus periods
	[vecPseudoSpikeTimes,vecPseudoEventT] = getPseudoSpikeVectors(cellSpikeTimes{intN},vecStimOnTime-dblPreTime,dblMaxDur,boolDiscardEdges);
	
	%real
	[vecTrialPerSpike,vecTimePerSpike] = getSpikesInTrial(vecPseudoSpikeTimes,vecPseudoEventT,dblMaxDur);
	for intTrial=1:intTrialNum
		vecSpikeT = vecTimePerSpike(vecTrialPerSpike==intTrial);
		cellSpikeTimesPerCellPerTrial{intN,intTrial} = vecSpikeT;
	end
	
	%calc non-stationarity
	vecSortedSpikeTimes = sort(vecPseudoSpikeTimes,'ascend') - min(vecPseudoSpikeTimes);
	dblAUC = sum(vecSortedSpikeTimes);
	dblLinAUC = (max(vecSortedSpikeTimes) * numel(vecSortedSpikeTimes) ) / 2;
	vecNonStat(intN) = (dblAUC - dblLinAUC) / dblLinAUC;
end
vecStimOnStitched = vecPseudoEventT;
matDataZ = zscore(log(1+matData),[],2);
vecMeanZ = mean(matDataZ,1);
vecFilt = normpdf(-4:4,0,1)/sum(normpdf(-4:4,0,1));
vecFiltM = imfilt(vecMeanZ,vecFilt);

%calc metrics
[h,pKS,ksstat,cv] = kstest(vecFiltM);
[BF, dblBC] = bimodalitycoeff(vecFiltM);
dblMaxDevFrac = max(abs(vecFiltM));

if boolSaveFigs
	figure;maxfig;
	subplot(1,2,1);
	imagesc(matDataZ);
	fixfig;grid off;
	subplot(1,2,2);
	plot(vecFiltM);
	title(sprintf('%s: K-S test,p =%.1e; BC=%.3f; Dev =%.3f',strRecOrig,pKS,dblBC,dblMaxDevFrac),'interpreter','none');
	drawnow;
	fixfig;
	
	%% save fig
	export_fig(fullpath(strFigurePathSR,sprintf('BimoCheck%s_%s.tif',strRunStim,strRecOrig)));
	export_fig(fullpath(strFigurePathSR,sprintf('BimoCheck%s_%s.pdf',strRunStim,strRecOrig)));
end

%if dblBC > dblBimoThreshold || dblMaxDevFrac > dblDevThreshold
%	fprintf('Population drift is bimodal (BC=%.3f) for %s: skipping... [%s]\n',dblBC,strRec,getTime);
%	continue;
%end

%% define quantiles and remove zero-variance neurons
%constants
[vecTrialTypeIdx,vecUnique,vecPriorDistribution,cellSelect,vecRepetition] = val2idx(vecStimIdx);
intStimNum = numel(vecUnique);
dblLambda = 1;%1
intTypeCV = 2;
dblUseStartT = 0;
dblUseMaxDur = dblMaxDur-dblUseStartT;
intUseMax = inf;
intRepNum = mean(vecPriorDistribution);

%remove zero-variance neurons
matSpikeCounts_pre = cellfun(@(x) sum(x>dblUseStartT),cellSpikeTimesPerCellPerTrial);
matMeanRate_pre = matSpikeCounts_pre./dblUseMaxDur;
vecPopRate_pre = sum(matMeanRate_pre,1);

intQuantiles = 5;
vecStartTrials = round(intRepNum*linspace(1/intRepNum,(1+1/intRepNum),intQuantiles+1));
vecStartTrials(end)=[];
intSplitTrialsPerOri = min(floor(cellfun(@sum,cellSelect)/intQuantiles));
indZeroVarNeurons = false(intRespN,1);
for intQ=1:intQuantiles
	vecUseTrialsTemp = vecStartTrials(intQ):(vecStartTrials(intQ)+intSplitTrialsPerOri-1);
	for intStim=1:intStimNum
		vecThisStim = find(cellSelect{intStim});
		[vecSorted,vecReorder]=sort(vecPopRate_pre(vecThisStim));
		vecQualifyingTrials = vecThisStim(vecReorder(vecUseTrialsTemp));
		indZeroVarNeurons = indZeroVarNeurons | (var(matMeanRate_pre(:,vecQualifyingTrials),[],2) == 0);
	end
end
indZeroVarNeurons = false(size(indZeroVarNeurons));
vecUseNeurons = find(~indZeroVarNeurons);
vecRemNeurons = find(indZeroVarNeurons);
cellSpikeTimesPerCellPerTrial(vecRemNeurons,:) = [];
intNeuronNum = numel(vecUseNeurons);
intTrialsPerQ = intSplitTrialsPerOri*intStimNum;

%simple "rate code"
matSpikeCounts = cellfun(@(x) sum(x>dblUseStartT),cellSpikeTimesPerCellPerTrial);
matMeanRate = matSpikeCounts./dblUseMaxDur;
