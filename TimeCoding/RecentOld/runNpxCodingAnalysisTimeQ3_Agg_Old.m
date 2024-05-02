%% aim
%{
https://www.nature.com/articles/s41598-021-95037-z
https://www.nature.com/articles/s42003-021-02437-y

q1: Progression of orientation information over time: how does trial activity evolve? what is the function of the onset peak?
> is decoding better when matched for stimulus phase? => no.

q2: "Spike time and rate coding can be represented within a single model of spiking probability as a
function of time: rate codes are uniform over a certain period tau, while spike time codes are
temporally localized peaks"
> is this true?

q3: Rate codes do not exist; a rate code is simply a subset of spike time codes where the temporal
integration window is very large. But what about multi dim codes? Those are all rate based. Can we
formulate a multidimensional spike-time code? I.e., can we make a taxonomy of neural codes?

q4: How does information evolve over time, is initial peak indeed less tuned? Is pop activity
rhythmic? Are stimuli encoded invariant to brain state? Eg, high arousal, low arousal. Or is
stimulus manifold dynamic over time? Does manifold scale with arousal? => How does manifold depend
on binning size? What is the optimal time window?

%}
%% define qualifying areas
clear all;
cellUseAreas = {...
	'Primary visual area',...
	...'posteromedial visual area',...
	};
boolHome = false;
if boolHome
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
	[sAggStim,sAggNeuron]=loadDataNpx('','driftinggrating',strDataPath);
end
indRemDBA = strcmpi({sAggNeuron.SubjectType},'DBA');
fprintf('Removing %d cells of DBA animals; %d remaining [%s]\n',sum(indRemDBA),sum(~indRemDBA),getTime);
sAggNeuron(indRemDBA) = [];

%% pre-allocate matrices
intAreas = numel(cellUseAreas);
matDecPerf = [];
dblStartT = 0.1;

%% go through recordings
tic
for intRec=1:numel(sAggStim)
	clearvars -except sAggStim sAggNeuron intRec dblStartT intAreas cellUseAreas strDataPath strFigurePathSR strFigurePath strTargetDataPath
	close all;
	% get matching recording data
	strRec = sAggStim(intRec).Exp;
	sThisRec = sAggStim(strcmpi(strRec,{sAggStim(:).Exp}));
	
	%prep grating data
	[sUseNeuron,vecStimOnTime,vecStimOffTime,vecOrientation] = NpxPrepGrating(sAggNeuron,sThisRec,cellUseAreas);
	[vecOriIdx,vecUniqueOris,vecRepNum,cellSelect,vecTrialRepetition] = val2idx(vecOrientation);
	intTrialNum = numel(vecStimOnTime);
	intStimNum = numel(unique(vecOrientation));
	intRepNum = intTrialNum/intStimNum;
	numel(sUseNeuron)
	
	%% select area 1
	for intArea=1:intAreas
		strArea = cellUseAreas{intArea};
		indArea1Neurons = contains({sUseNeuron.Area},strArea,'IgnoreCase',true);
		if sum(indArea1Neurons) == 0, continue;end
		%% get orientation responses & single-trial population noise
		sArea1Neurons = sUseNeuron(indArea1Neurons);
		intRec
		
		%% prep data
		%get data matrix
		cellSpikeTimesRaw = {sArea1Neurons.SpikeTimes};
		[matData,indTuned,cellSpikeTimes,sOut,cellSpikeTimesPerCellPerTrial,vecStimOnStitched,vecNonStat,dblBC,dblMaxDevFrac] = ...
			NpxPrepData(cellSpikeTimesRaw,vecStimOnTime,vecStimOffTime,vecOrientation);
		intTunedN = sum(indTuned);
		intNumN = size(matData,1);
		
		%% is pop activity more stable over trials as euclidian distance to origin or as mean?
		%mean is L1 norm divided by # of neurons, so scaling by mean magnitude makes them
		%equaivalent. If population activity showed a constant mean over trials, it would live on a
		%hyperplane orthogonal to the diagonal. If, instead, the population activity shows a
		%constant Euclidian distance to the origin, this means it lives on a hypersphere. L1 and L2
		%norms have intrinsically different magnitudes, however: L1 norms are larger than L2 norms
		%(except in the case where all neurons but one are silent), so a fair comparison requires
		%comparing the variability of L1 and L2 norms after normalizing for the average magnitude.
		
		%split data into reps per stim type
		[vecOriIdx,vecUniqueOris,vecRepNum,cellSelect,vecTrialRepetition] = val2idx(vecOrientation);
		[matRespNSR,vecStimTypes,vecUnique] = getStimulusResponses(matData,vecOrientation);
		intTypeNum = numel(vecUnique);
		matSumPopSignal = sum(matRespNSR,1);
		matDistR_L1 = nan(mean(vecRepNum),intTypeNum);
		matDistR_L2 = nan(mean(vecRepNum),intTypeNum);
		
		%per stim type
		vecRelVarMean = nan(1,intTypeNum);
		vecRelVarNorm = nan(1,intTypeNum);
		for intStimType=1:intTypeNum
			vecDistR_L1 = squeeze(matSumPopSignal(1,intStimType,:)); %hyperplane: mean distance
			vecDistR_L2 = nan(size(vecDistR_L1)); %hypersphere: euclidian distance
			matR = squeeze(matRespNSR(:,intStimType,:));
			for intRep=1:size(matR,2)
				vecDistR_L2(intRep) = norm(matR(:,intRep));
			end
			matDistR_L1(:,intStimType) = vecDistR_L1;
			matDistR_L2(:,intStimType) = vecDistR_L2;
			
			vecRelVarMean(intStimType) = std(vecDistR_L1)/mean(vecDistR_L1);
			vecRelVarNorm(intStimType) = std(vecDistR_L2)/mean(vecDistR_L2);
		end
		
		vecPercDiff = 100*((vecRelVarNorm./vecRelVarMean)-1);
		vecBins = -20:2.5:20;
		vecBinsC = vecBins(2:end)-median(diff(vecBins))/2;
		vecCounts =histcounts(vecPercDiff,vecBins);
		
		dblRelVarMean = std(matDistR_L1(:))/mean(matDistR_L1(:));
		dblRelVarNorm = std(matDistR_L2(:))/mean(matDistR_L2(:));
		
		matDistR_L1_Norm = matDistR_L1 ./ mean(matDistR_L1(:));
		matDistR_L2_Norm = matDistR_L2 ./ mean(matDistR_L2(:));
		
		if 0
			%plot
			figure;maxfig;
			subplot(2,3,1)
			bar(vecBinsC,vecCounts,'hist');
			xlabel('Variability of pop act. as Euclidian vs mean distance ((L2/L1) - 1) (%)');
			title('Trial-to-trial variability of pop act per stim type');
			ylabel('Number of stimulus types');
			fixfig;
			
			%overall
			subplot(2,3,2)
			%xlim([0 25]);
			hold on
			for intStimType=1:intTypeNum
				h=scatter(vecUniqueOris(intStimType)*ones(size(matDistR_L1(:,intStimType))),matDistR_L1(:,intStimType));
			end
			ylim([0 max(get(gca,'ylim'))]);
			xlabel('Stimulus type (direction in degrees)');
			title('L1 norm per trial');
			ylabel('L1 distance of population activity (Hz)');
			fixfig;
			
			subplot(2,3,3)
			hold on
			for intStimType=1:intTypeNum
				h=scatter(vecUniqueOris(intStimType)*ones(size(matDistR_L1(:,intStimType))),matDistR_L2(:,intStimType));
			end
			ylim([0 max(get(gca,'ylim'))]);
			xlabel('Stimulus type (direction in degrees)');
			title('L2 norm per trial');
			ylabel('L2 distance of population activity (Hz)');
			fixfig;
			
			subplot(2,3,4)
			scatter(matDistR_L1_Norm(:),matDistR_L2_Norm(:));
			hold on
			plot([0 max(get(gca,'xlim'))],[0 max(get(gca,'ylim'))],'k--');
			title('Normalized by averaged magnitude');
			xlabel('L1-norm of trial pop act');
			ylabel('L2-norm of trial pop act');
			fixfig;
			
			%{
		subplot(2,3,5)
		vecDiffNorm = 100*((matDistR_L2_Norm(:)./matDistR_L1_Norm(:))-1);
		vecBins = -50:2:50;
		vecBinsC = vecBins(2:end)-median(diff(vecBins))/2;
		vecCounts =histcounts(vecDiffNorm,vecBins);
		hold on
		bar(vecBinsC,vecCounts,'hist');
		xlabel('Variability of pop act. as Euclidian vs mean distance ((L2/L1) - 1) (%)');
		title('Trial-to-trial variability of pop act per trial');
		ylabel('Number of trials');
		fixfig;
			%}
			
			subplot(2,3,5)
			title(strRec,'interpreter','none');
			
			% save fig
			export_fig(fullpath(strFigurePathSR,sprintf('C1_VarofL1vsL2_%s.tif',strRec)));
			export_fig(fullpath(strFigurePathSR,sprintf('C1_VarofL1vsL2_%s.pdf',strRec)));
		end
		%% save data
		save(fullpath(strTargetDataPath,sprintf('Q3Data_%s',strRec)),...
			'dblStartT',...
			'strRec',...
			'vecRelVarMean',...
			'vecRelVarNorm',...
			'dblRelVarMean',...
			'dblRelVarNorm'...
			);
	end
end
toc
