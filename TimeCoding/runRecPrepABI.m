% get matching recording data
sRec = sAggABI(intRec);
strRecOrig = sRec.Exp;
strRec = strrep(strRecOrig,'ecephys_session_','ABI');
fprintf('Preparing %d: %s, area %s [%s]\n',intRec,strRecOrig,strArea,getTime);

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
elseif strcmp(strRunStim,'NS') && isfield(sRec.structStimAgg,strField)
	%% NS; natural scenes with task
	% concatenate stimulus structures
	structStim = sRec.structStimAgg.sNS; %sDG, sNM, sNS
	
	
	%remove omitted stimuli
	cellIms = cellfun(@getFlankedBy,structStim.image_name,cellfill('im',size(structStim.image_name)),cellfill('_',size(structStim.image_name)),...
		'UniformOutput',false);
	vecIms = str2double(cellIms);
	structStim.image_nr = vecIms;
	indRem = isnan(vecIms);
	cellFields = fieldnames(structStim);
	for i=1:numel(cellFields)
		varVals = structStim.(cellFields{i});
		structStim.(cellFields{i}) = varVals(~indRem);
	end
	
	%transform names
	vecStimLabels = structStim.image_nr;
	vecStimOnTime = structStim.start_time(:)';
	vecStimOffTime = structStim.stop_time(:)';
	
	vecStimIsChange = structStim.is_change(:)';
	vecStimIsRewarded = structStim.rewarded(:)';
	vecStimTrialId = structStim.trials_id(:)';
	
	%remove nans
	indRem = isnan(vecStimLabels);
	vecStimLabels(indRem) = [];
	vecStimOnTime(indRem) = [];
	vecStimOffTime(indRem) = [];
	
end

%check if trials structure exists
if isfield(sRec.structStimAgg,'sTrials')
	sTrials = sRec.structStimAgg.sTrials;
end

%general prep
if ~strcmp(strRunStim,'NS')
	%equalize reps
	[vecStimIdx,vecUniqueStims,vecRepNum,cellSelect,vecTrialRepetition] = val2idx(vecStimLabels);
	intRepNum = min(vecRepNum);
	indRem2 = vecTrialRepetition>intRepNum;
	vecStimLabels(indRem2) = [];
	vecStimOnTime(indRem2) = [];
	vecStimOffTime(indRem2) = [];
end
[vecStimIdx,vecUniqueStims,vecRepNum,cellSelect,vecTrialRepetition] = val2idx(vecStimLabels);
intTrialNum = numel(vecStimOnTime);
intStimNum = numel(vecUniqueStims);


%% get neurons in this area
indArea1Neurons = strcmp(sRec.cellAreas,strArea);
intNeuronsInArea = sum(indArea1Neurons);
intNeuronNum = intNeuronsInArea;
if intNeuronsInArea==0,return;end

%% prep data
dblMinRate = 0.1;

%get dur
cellSpikes = {sRec.sNeuron.SpikeTimes};
cellSpikeTimes = cellSpikes(indArea1Neurons);
dblDur = median(vecStimOffTime-vecStimOnTime);
matRawData = getSpikeCounts(cellSpikeTimes,vecStimOnTime,dblDur)./dblDur;
indResp = sum(matRawData,2)'>(size(matRawData,2)/dblDur)*dblMinRate;
matData = matRawData;%(indResp,:);

%% check non-stationarity
%get spikes per trial per neuron
cellSpikeTimesPerCellPerTrial = cell(numel(cellSpikeTimes),intTrialNum);
vecNonStat = nan(1,numel(cellSpikeTimes));
boolDiscardEdges = true;
for intN=1:numel(cellSpikeTimes)
	% build pseudo data, stitching stimulus periods
	[vecPseudoSpikeTimes,vecPseudoEventT] = getPseudoSpikeVectors(cellSpikeTimes{intN},vecStimOnTime,dblDur,boolDiscardEdges);
	
	%real
	[vecTrialPerSpike,vecTimePerSpike] = getSpikesInTrial(vecPseudoSpikeTimes,vecPseudoEventT,dblDur);
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

dblStimDur = roundi(median(vecStimOffTime - vecStimOnTime),1,'ceil');
dblPreTime = 0;%0.3;
dblPostTime = 0;%0.3;
dblMaxDur = dblStimDur+dblPreTime+dblPostTime;

boolSaveFigs = true;
if boolSaveFigs
	figure;maxfig;
	subplot(1,2,1);
	imagesc(matDataZ);
	ylabel(sprintf('Neuron in %s',strArea))
	xlabel(sprintf('Stimulus # (%s)',strRunStim))
	title(sprintf('Norm act in %s',strRecOrig),'interpreter','none')
	
	fixfig;grid off;
	subplot(1,2,2);
	plot(vecFiltM);
	xlabel(sprintf('Stimulus # (%s)',strRunStim))
	ylabel(sprintf('Avg norm act (z-score)'))
	title(sprintf('%s: K-S test,p =%.1e; BC=%.3f; Dev =%.3f',strRecOrig,pKS,dblBC,dblMaxDevFrac),'interpreter','none');
	drawnow;
	fixfig;
	
	%% save fig
	export_fig(fullpath(strFigurePathSR,sprintf('BimoCheck%s_%s_%s.tif',strRunStim,strArea,strRecOrig)));
	export_fig(fullpath(strFigurePathSR,sprintf('BimoCheck%s_%s_%s.pdf',strRunStim,strArea,strRecOrig)));
	
	%delete fig
	close;
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
intRepNum = min(vecPriorDistribution);

%remove zero-variance neurons
matSpikeCounts_pre = cellfun(@(x) sum(x>dblUseStartT),cellSpikeTimesPerCellPerTrial);
matMeanRate_pre = matSpikeCounts_pre./dblUseMaxDur;
vecPopRate_pre = sum(matMeanRate_pre,1);

intQuantiles = 5;
vecStartTrials = round(intRepNum*linspace(1/intRepNum,(1+1/intRepNum),intQuantiles+1));
vecStartTrials(end)=[];
intSplitTrialsPerOri = min(floor(cellfun(@sum,cellSelect)/intQuantiles));
indZeroVarNeurons = false(size(matSpikeCounts_pre,1),1);
for intQ=1:intQuantiles
	vecUseTrialsTemp = vecStartTrials(intQ):(vecStartTrials(intQ)+intSplitTrialsPerOri-1);
	for intStim=1:intStimNum
		vecThisStim = find(cellSelect{intStim});
		[vecSorted,vecReorder]=sort(vecPopRate_pre(vecThisStim));
		vecQualifyingTrials = vecThisStim(vecReorder(vecUseTrialsTemp));
		indZeroVarNeurons = indZeroVarNeurons | (var(matMeanRate_pre(:,vecQualifyingTrials),[],2) == 0);
	end
end
%% remove non-responsive cells
indTuned = ~indZeroVarNeurons(:) & indResp(:);
vecUseNeurons = find(indTuned);
vecRemNeurons = find(~indTuned);
cellSpikeTimes = cellSpikeTimes(vecUseNeurons);
cellSpikeTimesPerCellPerTrial(vecRemNeurons,:) = [];
intNeuronNum = numel(vecUseNeurons);
intTrialsPerQ = intSplitTrialsPerOri*intStimNum;

%simple "rate code"
matSpikeCounts = cellfun(@(x) sum(x>dblUseStartT),cellSpikeTimesPerCellPerTrial);
matMeanRate = matSpikeCounts./dblUseMaxDur;
