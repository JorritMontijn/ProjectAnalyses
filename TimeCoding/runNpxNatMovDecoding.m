%% aim
%{
https://www.nature.com/articles/s41598-021-95037-z
https://www.nature.com/articles/s42003-021-02437-y

q1: can single-neuron variability be explained as population-level gain multiplication?
> estimate tuning curve from real data, then apply trial-by-trial gain to all neurons
>model predicting firing rate as combination of mean-rate multiplied by tuning curve => what do
residuals look like?
A: gain axis and mean-rate axis do not align, but ori tuning is distributed as conic manifold around
gain axis. Using stim-specific gain coupling calculated cross-validated for each neuron, allows pop
response to be predicted with R^2 of 0.72

%see articles:
https://elifesciences.org/articles/8998
https://www.nature.com/articles/nn.3711
https://www.jneurosci.org/content/39/37/7344.abstract
etc

%}
%% define qualifying areas
clearvars -except sAggStim sAggNeuron sAggSources;
%cellUseAreas = {'hippocampal formation','Dentate gyrus','Field CA1','Field CA2','Field CA3','subiculum'};
%cellUseAreas = {'Dentate gyrus','Field CA1','Field CA2','Field CA3'};
%cellUseAreas = {'posteromedial'};
cellUseAreas = {'anteromedial'};
	
intRandomize = 1; %1=real data, 2=shuffled, 3=generated
boolSaveFigs = true;
boolHome = true;
if boolHome
	strDataPath = 'F:\Data\Processed\Neuropixels\';
	strFigurePath = 'F:\Data\Results\NatMovDec\figures\';
	strTargetDataPath = 'F:\Data\Results\NatMovDec\data\';
else
	strDataPath = 'E:\DataPreProcessed\';
	strFigurePath = 'D:\Data\Results\NatMovDec\figures\';
	strTargetDataPath = 'D:\Data\Results\NatMovDec\data\';
end

%% select all neurons in LP and drifting grating stimuli
%strUseRec = '20191211_MP3_RunDriftingGratingsR01_g0_t0';
%strUseRec = '20191216_MP3_RunNaturalMovieR01_g0_t0';
if ~exist('sAggStim','var') || isempty(sAggStim)
	[sAggStim,sAggNeuron,sAggSources]=loadDataNpx('','natural',strDataPath);
end
vecUseRec = find(contains({sAggStim.Exp},'MP'));
if isempty(vecUseRec)
	fprintf('No recordings found\n');
	return
end
if ~isfolder(strFigurePath),mkdir(strFigurePath);end
if ~isfolder(strTargetDataPath),mkdir(strTargetDataPath);end

%% pre-allocate matrices
matDecPerf = [];
dblStartT = 0;

%% go through recordings
tic
for intRec=vecUseRec
	% get matching recording data
	strRec = sAggStim(intRec).Exp;
	%if ~strcmp(strRec,'20191213_MP3_RunDriftingGratingsR01_g0_t0'),continue;end
	%if ~strcmp(strRec,'20200116_MP4_RunDriftingGratingsR01_g0_t0'),continue;end
	if ~strcmp(strRec,'20191212_MP3_RunNaturalMovieR01_g0_t0'),continue;end
	
	
	sThisRec = sAggStim(strcmpi(strRec,{sAggStim(:).Exp}));
	sThisSource = sAggSources(strcmpi(strRec,{sAggSources(:).Exp}));
	
	%remove stimulus sets that are not 100 reps
	sThisRec.cellBlock(cellfun(@(x) x.intNumRepeats,sThisRec.cellBlock) < 100) = [];
	
	% concatenate stimulus structures
	structStim = catstim(sThisRec.cellBlock(1));
	vecOrigStimOnTime = structStim.vecStimOnTime;
	dblDur = median(diff(vecOrigStimOnTime));
	vecOrigStimOffTime = vecOrigStimOnTime+dblDur;
	intFramesInMovie = 500;
	dblBinAvg = 25;
	intUseBinsInMovie = round(intFramesInMovie/dblBinAvg);
	dblBinRate = roundi(intUseBinsInMovie/dblDur,2);
	dblBinDur = 1/dblBinRate;
	vecBinEdges = linspace(0,dblBinDur*intUseBinsInMovie,intUseBinsInMovie+1);
	
	%generate fake stimulus vectors
	vecUniqueStims = 1:intUseBinsInMovie;
	vecFrameIdx = repmat(vecUniqueStims,[1 numel(vecOrigStimOnTime)]);
	vecStimOnTime = flat(vecBinEdges(1:(end-1))' + vecOrigStimOnTime)';
	vecStimOffTime = vecStimOnTime + dblBinDur;
	
	[vecStimIdx,vecUniqueStims,vecRepNum,cellSelect,vecTrialRepetition] = val2idx(vecFrameIdx);
	intTrialNum = numel(vecStimOnTime);
	intStimNum = numel(vecUniqueStims);
	intRepNum = intTrialNum/intStimNum;
	%remove neurons from other recs
	indQualifyingNeurons = contains({sAggNeuron.Exp},strRec);
	sTheseNeurons = sAggNeuron(indQualifyingNeurons);
	
	%remove neurons in incorrect areas
	indConsiderNeurons = contains({sAggNeuron.Area},cellUseAreas,'IgnoreCase',true);
	
	%remove bad clusters
	indGoodNeurons = (cell2vec({sAggNeuron.KilosortGood}) == 1) | (cell2vec({sAggNeuron.Contamination}) < 0.1);
	
	%subselect from total
	indUseNeurons = indQualifyingNeurons(:) & indConsiderNeurons(:) & indGoodNeurons(:);
	sUseNeuron = sAggNeuron(indUseNeurons);
	if numel(sUseNeuron) < 10,continue;end
	
	%get data matrix
	cellSpikeTimes = {sUseNeuron.SpikeTimes};
	dblDur = median(vecStimOffTime-vecStimOnTime);
	matDataC = getSpikeCounts(cellSpikeTimes,vecStimOnTime,dblDur);
	matDataL = getSpikeCounts(cellSpikeTimes,vecStimOnTime-0.5,dblDur);
	matDataR = getSpikeCounts(cellSpikeTimes,vecStimOnTime+0.5,dblDur);
	matDataLL = getSpikeCounts(cellSpikeTimes,vecStimOnTime-1,dblDur);
	matDataRR = getSpikeCounts(cellSpikeTimes,vecStimOnTime+1,dblDur);
	matData = (matDataC + matDataL + matDataR)./3;
	[matRespNSR,vecStimTypes,vecUnique] = getStimulusResponses(matData,vecStimIdx);
	vecNonStat = cell2vec({sUseNeuron.NonStationarity});
	matAvgR = mean(matRespNSR,3);
	matTrialR = squeeze(mean(matRespNSR,2));
	
	%remove untuned cells
	sOut = getTuningCurves(matData,(vecFrameIdx/max(vecFrameIdx))*180,0);
	dblMinRate = 0.1;
	indTuned = sOut.vecOriAnova<0.05;
	indResp = abs(vecNonStat) < 0.3;% & indTuned;% & sum(matData,2)>(size(matData,2)/dblDur)*dblMinRate;
	if sum(indResp) < 10,continue;end
	
	%% decode
	vecTrialTypes = vecStimIdx;
	intTypeCV = 2;
	vecPriorDistribution = vecRepNum;
	dblLambda = 1;
	[dblPerformanceCV,vecDecodedIndexCV,matPosteriorProbability,dblMeanErrorDegs,matConfusion,matWeights,matAggActivation,matAggWeights,vecRepetition] = ...
		doCrossValidatedDecodingLR(matData(indResp,:),vecTrialTypes,intTypeCV,vecPriorDistribution,dblLambda);
	
	pBinom=myBinomTest(dblPerformanceCV*intTrialNum,intTrialNum,1/intStimNum);
	%%
	matFilt = normpdf(-2:2,0,0.5)' * normpdf(-2:2,0,0.5);
	matConfusion2 = imfilt(matConfusion,matFilt./sum(matFilt(:)));
	vecT = vecBinEdges(2:end) - mean(diff(vecBinEdges))/2;
	
	matConfusion3 = matConfusion(1:2:end,:) + matConfusion(2:2:end,:);
	matConfusion3 = matConfusion3(:,1:2:end) + matConfusion3(:,2:2:end);
	
	intReps = mean(vecRepNum);
	matConfusionPlot = matConfusion2;
	figure
	imagesc(vecT,vecT,matConfusionPlot./intReps,[0 20]./intReps);
	colormap(hot);
	vecTransitions = ([0 3 10 13 20]./20);
	hold on
	strCol = 'b';
	matBlockConf = nan(size(matConfusionPlot));
	for intTransition=1:(numel(vecTransitions)-1)
		for intTransition2=1:(numel(vecTransitions)-1)
			vecUseTrans = vecTransitions*size(matConfusionPlot,1);
			vecSelectX = round(vecUseTrans(intTransition2):vecUseTrans(intTransition2+1));
			vecSelectY = round(vecUseTrans(intTransition):vecUseTrans(intTransition+1));
			vecSelectX(vecSelectX==0)=[];
			vecSelectY(vecSelectY==0)=[];
			matBlockConf(vecSelectY,vecSelectX) = mean(flat(matConfusionPlot(vecSelectY,vecSelectX)));
		end
		dblT0 = vecTransitions(intTransition)*dblBinDur*intUseBinsInMovie;
		dblT1 = vecTransitions(intTransition+1)*dblBinDur*intUseBinsInMovie;
		plot([dblT1 dblT1],[dblT0 dblT1],strCol);
		plot([dblT0 dblT1],[dblT1 dblT1],strCol);
		plot([dblT1 dblT1],[dblT0 dblT1],strCol);
		plot([dblT0 dblT1],[dblT1 dblT1],strCol);
		
		plot([dblT0 dblT0],[dblT0 dblT1],strCol);
		plot([dblT0 dblT1],[dblT0 dblT0],strCol);
		plot([dblT0 dblT0],[dblT0 dblT1],strCol);
		plot([dblT0 dblT1],[dblT0 dblT0],strCol);
	end
	hold off
	colorbar;
	axis image;
	xlabel('Real movie time (s)');
	ylabel('Decoded movie time (s)');
	fixfig;
	grid off;
	title(sprintf('%s\n%s - %d%%, N=%d,p=%.2e',strRec,cellUseAreas{1},round((sum(diag(matConfusionPlot))./sum(matConfusionPlot(:)))*100),sum(indResp),pBinom),'interpreter','none');
	drawnow;
	strFigName = [cellUseAreas{1} '_' strRec];
	export_fig([strFigurePath strFigName '.tif']);
	export_fig([strFigurePath strFigName '.pdf']);
	
	%% save
	strTargetFile = [strTargetDataPath sprintf('NatMovDec_%s',strRec)];
	save(strTargetFile,...
		'dblPerformanceCV','vecDecodedIndexCV','matPosteriorProbability','dblMeanErrorDegs','matConfusion','matWeights','matAggActivation',...
		'matAggWeights','vecRepetition');
	fprintf('Saved data to %s [%s]\n',strTargetFile,getTime);
end
toc
