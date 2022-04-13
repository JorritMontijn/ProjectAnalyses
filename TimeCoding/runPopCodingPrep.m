%% remove untuned cells
%get data matrix
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
	export_fig(fullpath(strFigurePath,sprintf('BimoCheck%s_%s.tif',strRunStim,strRecOrig)));
	export_fig(fullpath(strFigurePath,sprintf('BimoCheck%s_%s.pdf',strRunStim,strRecOrig)));
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
intReps = mean(vecPriorDistribution);

%remove zero-variance neurons
matSpikeCounts_pre = cellfun(@(x) sum(x>dblUseStartT),cellSpikeTimesPerCellPerTrial);
matMeanRate_pre = matSpikeCounts_pre./dblUseMaxDur;
vecPopRate_pre = sum(matMeanRate_pre,1);

intQuantiles = 5;
vecStartTrials = round(intReps*linspace(1/intReps,(1+1/intReps),intQuantiles+1));
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

