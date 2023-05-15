%% aim
%{
https://www.nature.com/articles/s41598-021-95037-z
https://www.nature.com/articles/s42003-021-02437-y

q1: Can absence of stimulus information in initial period be explained by neurons responding only to
black/white patches in the first x ms? 

q2: How precise are spike times in the natural movie repetitions?


%}
%% define qualifying areas
clear all;
cellUseAreas = {...
	'Primary visual area',...
	...'posteromedial visual area',...
	};
boolHome = false;
if isfolder('F:\Drive\PopTimeCoding') && isfolder('F:\Data\Processed\Neuropixels\')
	strDataPath = 'F:\Data\Processed\Neuropixels\';
	strFigurePathSR = 'F:\Drive\PopTimeCoding\single_recs';
	strFigurePath = 'F:\Drive\PopTimeCoding\figures\';
	strTargetDataPath = 'F:\Drive\PopTimeCoding\data\';
else
	strDataPath = 'E:\DataPreProcessed\';
	strFigurePathSR = 'C:\Drive\PopTimeCoding\single_recs';
	strFigurePath = 'C:\Drive\PopTimeCoding\figures\';
	strTargetDataPath = 'C:\Drive\PopTimeCoding\data\';
end


%% select all neurons in LP and drifting grating stimuli
if ~exist('sAggStim','var') || isempty(sAggStim) || isempty(sAggNeuron)
	[sAggStim,sAggNeuron]=loadDataNpx('','natural',strDataPath);
end
indRemDBA = strcmpi({sAggNeuron.SubjectType},'DBA');
fprintf('Removing %d cells of DBA animals; %d remaining [%s]\n',sum(indRemDBA),sum(~indRemDBA),getTime);
sAggNeuron(indRemDBA) = [];
	
%% pre-allocate matrices
intAreas = numel(cellUseAreas);

%% go through recordings
tic
for intRec=1:numel(sAggStim) %19 || weird: 11
	% get matching recording data
	strRec = sAggStim(intRec).Exp;
	sThisRec = sAggStim(strcmpi(strRec,{sAggStim(:).Exp}));
	
	%prep nm data
	[sUseNeuron,vecStimOnTime,vecStimOffTime,vecStimIdx] = NpxPrepMovies(sAggNeuron,sThisRec,cellUseAreas);
	[vecFrameIdx,vecUniqueFrames,vecRepNum,cellSelect,vecTrialRepetition] = val2idx(vecStimIdx);
	intTrialNum = numel(vecStimOnTime);
	intStimNum = numel(unique(vecStimIdx));
	intRepNum = intTrialNum/intStimNum;
	vecOrigStimOnTime = sThisRec.cellBlock{1}.vecStimOnTime;
	dblStimDur = min(diff(vecOrigStimOnTime));
	intOrigTrialNum = numel(vecOrigStimOnTime);
	
	%% select area 1
	for intArea=1:intAreas
		strArea = cellUseAreas{intArea};
		indArea1Neurons = contains({sUseNeuron.Area},strArea,'IgnoreCase',true);
		if sum(indArea1Neurons) == 0, continue;end
		
		%% get orientation responses & single-trial population noise
		sArea1Neurons = sUseNeuron(indArea1Neurons);
		
		%% prep data
		%get data matrix
		cellSpikeTimesRaw = {sArea1Neurons.SpikeTimes};
		[matData,cellSpikeTimes] = NpxPrepMovieData(cellSpikeTimesRaw,vecStimOnTime,vecStimOffTime,vecFrameIdx);
		intNumN = size(matData,1);
		dblBinDur = median(vecStimOffTime-vecStimOnTime);
	
		
		%% run analysis
		%intCombEntries = 10;
		%matDataReduced = zeros(size(matData,1),size(matData,2)/intCombEntries);
		%for intRedIdx=1:intCombEntries
		%	matDataReduced = matDataReduced + matData(:,intRedIdx:intCombEntries:end);
		%end
		matDataZ = zscore(matData,[],2);
		[matRespNSRZ,vecStimTypes,vecUnique] = getStimulusResponses(matDataZ.*dblBinDur,vecFrameIdx);
		
		[matRespNSR,vecStimTypes,vecUnique] = getStimulusResponses(matData.*dblBinDur,vecFrameIdx);
		
		vecBinEdges = 0:dblBinDur:(dblStimDur+0.05);
		intBinNum = numel(vecBinEdges)-1;
		matRate = nan(intNumN,intBinNum);
		matMinDt = nan(intNumN,intBinNum);
		matMinDt_Raw = nan(intNumN,intBinNum);
		for intNeuron=1:intNumN
			intNeuron
			vecSpikeTimes = cellSpikeTimes{intNeuron};
			[vecTrialPerSpike,vecTimePerSpike] = getSpikesInTrial(vecSpikeTimes,vecOrigStimOnTime,dblStimDur);
			
			%the expected distance to the nearest neighbour for a point on the unit circle with N-1
			%other points is pi/n
			%https://math.stackexchange.com/questions/2931257/the-expected-distance-to-the-nearest-neighbor-when-n-points-are-placed-randomly
			
			%make times circular
			vecTimePerSpikePhase = (vecTimePerSpike/dblStimDur)*2*pi;
			vecMinDt = nan(size(vecTimePerSpikePhase));
			vecMinDt_Norm = nan(size(vecTimePerSpikePhase));
			for intTrial=1:intRepNum
				indTrialSpikes = vecTrialPerSpike==intTrial;
				matDist = circ_dist(vecTimePerSpikePhase(indTrialSpikes)',vecTimePerSpikePhase(~indTrialSpikes));
				vecMinDt(indTrialSpikes) = min(abs(matDist));
				vecMinDt_Norm(indTrialSpikes) = min(abs(matDist))/(pi/(sum(~indTrialSpikes)+1));
			end
			vecMinDt = (vecMinDt/(2*pi))*dblStimDur;

			[vecCounts,vecMeans,vecSDs,cellVals,cellIDs] = makeBins(vecTimePerSpike,vecMinDt_Norm,vecBinEdges);
			[vecCounts,vecMeans_Raw,vecSDs,cellVals,cellIDs] = makeBins(vecTimePerSpike,vecMinDt,vecBinEdges);

			matRate(intNeuron,:) = (vecCounts./dblBinDur)./intRepNum;
			matMinDt(intNeuron,:) = vecMeans;
			matMinDt_Raw(intNeuron,:) = vecMeans_Raw;
		end
		
		%% plot
		subplot(2,3,1)
		imagesc(matRate)
		colorbar;
		subplot(2,3,2)
		imagesc(matMinDt,[0 2])
		colorbar;
		subplot(2,3,3)
		imagesc(matMinDt_Raw)
		colorbar;
		
		subplot(2,3,4)
		imagesc(log(matRate))
		colorbar;
		subplot(2,3,5)
		imagesc(log(matMinDt),[-2 2])
		colormap(redblue)
		colorbar;
		subplot(2,3,6)
		imagesc(log(matMinDt_Raw))
		colorbar;
		
		figure
		subplot(2,3,1)
		imagesc(mean(matRespNSR,3))
		
		subplot(2,3,2)
		imagesc(matRate)
		
		subplot(2,3,3)
		scatter(log(matRate(:)),log(matMinDt(:)))
		
		subplot(2,3,4)
		imagesc(mean(matRespNSRZ,3))
		
		subplot(2,3,5)
		imagesc(zscore(matRate,[],2));
		
		subplot(2,3,6)
		scatter(matRate(1,:),matMinDt(1,:))
		return
	end
	close all;
end
toc
