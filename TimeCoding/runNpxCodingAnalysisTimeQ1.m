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
	strFigurePath = 'F:\Data\Results\PopTimeCoding';
else
	strDataPath = 'E:\DataPreProcessed\';
	strFigurePath = 'D:\Data\Results\PopTimeCoding';
end

%% select all neurons in LP and drifting grating stimuli
if ~exist('sAggStim','var') || isempty(sAggStim)
	[sAggStim,sAggNeuron]=loadDataNpx('','driftinggrating',strDataPath);
end
sAggNeuron(strcmpi({sAggNeuron.SubjectType},'DBA')) = [];

%% pre-allocate matrices
intAreas = numel(cellUseAreas);
matR_PP_All = nan(intAreas,intAreas,numel(sAggStim),2);
matR_OP_All = nan(intAreas,intAreas,numel(sAggStim),2);
matR_OO_All = nan(intAreas,intAreas,numel(sAggStim),2);
matR_PO_All = nan(intAreas,intAreas,numel(sAggStim),2);
mat_xR_All = nan(intAreas,intAreas,numel(sAggStim),2);
matDecPerf = [];
matDecConf = [];

%% go through recordings
tic
for intRec=19%1:numel(sAggStim)
	% get matching recording data
	strRec = sAggStim(intRec).Exp;
	sThisRec = sAggStim(strcmpi(strRec,{sAggStim(:).Exp}));
	
	%remove stimulus sets that are not 24 stim types
	sThisRec.cellBlock(cellfun(@(x) x.intTrialNum/x.intNumRepeats,sThisRec.cellBlock) ~= 24) = [];
	% concatenate stimulus structures
	structStim = catstim(sThisRec.cellBlock);
	vecStimOnTime = structStim.vecStimOnTime;
	vecStimOffTime = structStim.vecStimOffTime;
	vecOrientation = cell2vec({structStim.sStimObject(structStim.vecTrialStimTypes).Orientation})';
	vecTempFreq = cell2vec({structStim.sStimObject(structStim.vecTrialStimTypes).TemporalFrequency})';
	vecPhase = structStim.Phase;
	vecDelayTimeBy = vecPhase./vecTempFreq;
	[vecOriIdx,vecUniqueOris,vecRepNum,cellSelect,vecTrialRepetition]  = val2idx(vecOrientation);
	indRem=vecTrialRepetition>min(vecRepNum);
	vecOrientation(indRem) = [];
	vecDelayTimeBy(indRem) = [];
	vecStimOnTime(indRem) = [];
	vecStimOffTime(indRem) = [];
	[vecOriIdx,vecUniqueOris,vecRepNum,cellSelect,vecTrialRepetition]  = val2idx(vecOrientation);
	intTrialNum = numel(vecStimOnTime);
	intOriNum = numel(unique(vecOrientation));
	intRepNum = intTrialNum/intOriNum;
	
	%remove neurons from other recs
	indQualifyingNeurons = contains({sAggNeuron.Exp},strRec);
	
	%remove neurons in incorrect areas
	indConsiderNeurons = contains({sAggNeuron.Area},cellUseAreas,'IgnoreCase',true);
	
	%remove bad clusters
	indGoodNeurons = (cell2vec({sAggNeuron.KilosortGood}) == 1) | (cell2vec({sAggNeuron.Contamination}) < 0.1);
	
	%subselect from total
	indUseNeurons = indQualifyingNeurons(:) & indConsiderNeurons(:) & indGoodNeurons(:);
	sUseNeuron = sAggNeuron(indUseNeurons);
	
	%% select area 1
	for intArea=1:intAreas
		strArea = cellUseAreas{intArea};
		indArea1Neurons = contains({sUseNeuron.Area},strArea,'IgnoreCase',true);
		if sum(indArea1Neurons) == 0, continue;end
		%% get orientation responses & single-trial population noise
		sArea1Neurons = sUseNeuron(indArea1Neurons);
		intRec
		%% get spike times
		cellSpikeTimes = {sArea1Neurons.SpikeTimes};
		intNumN = numel(cellSpikeTimes);
		dblStimDur = roundi(median(vecStimOffTime - vecStimOnTime),1,'ceil');
		dblPreTime = 0.3;
		dblPostTime = 0.3;
		dblMaxDur = dblStimDur+dblPreTime+dblPostTime;
		dblBinWidth = 0.1;
		vecBinEdges = 0:dblBinWidth:dblMaxDur;
		vecStimTime = vecBinEdges(2:end)-dblBinWidth/2 - dblPreTime;
		indStimBins = vecStimTime > 0 & vecStimTime < dblStimDur;
		intBinNum = numel(vecBinEdges)-1;
		matBNSR = nan(intBinNum,intNumN,intOriNum,intRepNum);
		matBNT = nan(intBinNum,intNumN,intTrialNum);
		%matBNT_shifted = nan(intBinNum,intNumN,intTrialNum);
		
		vecRepCounter = zeros(1,intOriNum);
		for intN=1:intNumN
			[vecTrialPerSpike,vecTimePerSpike] = getSpikesInTrial(cellSpikeTimes{intN},vecStimOnTime-dblPreTime,dblMaxDur);
			for intTrial=1:intTrialNum
				if intN==1
					intTrialOriIdx = vecOriIdx(intTrial);
					vecRepCounter(intTrialOriIdx) = vecRepCounter(intTrialOriIdx) + 1;
					intRep = vecRepCounter(intTrialOriIdx);
				end
				vecSpikeT = vecTimePerSpike(vecTrialPerSpike==intTrial);
				vecSpikeC = histcounts(vecSpikeT,vecBinEdges);
				matBNSR(:,intN,intTrialOriIdx,intRep) = vecSpikeC;
				matBNT(:,intN,intTrial) = vecSpikeC;
				%matBNT_shifted(:,intN,intTrial) = vecSpikeC;
				%matBNT_shifted(indStimBins,intN,intTrial) = circshift(matBNT_shifted(indStimBins,intN,intTrial),-round(intBinNum*vecDelayTimeBy(intTrial)));
			end
		end
		
		%% time progression
		matDecPerf(:,end+1)=nan;
		matDecConf(:,end+1)=nan;
		dblLambda = 1;
		intTypeCV = 2;
		vecOri180 = mod(vecOrientation,180)*2;
		vecRepNum180 = vecRepNum(1:12)*2;
		intOriNum180 = intOriNum/2;
		matDecConfusion = nan(intOriNum180,intOriNum180,intBinNum);
		vecSpikesPerBin = mean(sum(matBNT,2),3)';
		matAcrossTimeDecoder = nan(intBinNum,intBinNum);
		[vecTrialTypeIdx,vecUnique,vecPriorDistribution,cellSelect,vecRepetition] = val2idx(vecOri180);
		for intBinIdx=1:intBinNum
			intBinIdx
			[dblPerformanceCV,vecDecodedIndexCV,matPosteriorProbability,dblMeanErrorDegs,matConfusion,matWeights] = ...
				doCrossValidatedDecodingLR(squeeze(matBNT(intBinIdx,:,:)),vecOri180,intTypeCV,vecPriorDistribution,dblLambda);
			vecDecErr(intBinIdx) = dblMeanErrorDegs;
			matDecPerf(intBinIdx,end) = dblPerformanceCV;
			vecConf = nan(size(vecDecodedIndexCV));
			for intTrial=1:numel(vecTrialTypeIdx)
				vecConf(intTrial) = matPosteriorProbability(vecTrialTypeIdx(intTrial),intTrial);
			end
			matDecConf(intBinIdx,end) = mean(vecConf);
			matDecConfusion(:,:,intBinIdx) = matConfusion;
			
			%% apply on all bins
			for intTestBinIdx=1:intBinNum
				if intTestBinIdx == intBinIdx
					matAcrossTimeDecoder(intBinIdx,intTestBinIdx) = dblPerformanceCV;
					continue;
				end
				%get performance
				matTestData = squeeze(matBNT(intTestBinIdx,:,:));
				matDataPlusLin = [matTestData; ones(1,size(matTestData,2))];
				matActivation = matWeights'*matDataPlusLin;
				matPosteriorProbability = exp(bsxfun(@minus,matActivation,logsumexp(matActivation,1))); %softmax
				intTrials = size(matTestData,2);
				vecPriorCopy = vecPriorDistribution;
				% normal decoding or with prior distro?
				if isempty(vecPriorCopy)
					%calculate output
					[dummy, vecDecodedIndexCV] = max(matPosteriorProbability,[],1);
				else
					%% loop through trials and assign next most certain trial
					vecDecodedIndexCV = nan(intTrials,1);
					indAssignedTrials = false(intTrials,1);
					matTempProbs = matPosteriorProbability;
					for intTrial=1:intTrials
						%check if we're done
						if sum(vecPriorCopy==0)==(numel(vecPriorCopy)-1)
							vecDecodedIndexCV(~indAssignedTrials) = find(vecPriorCopy>0);
							break;
						end
						
						%remove trials of type that has been chosen max number
						matTempProbs(vecPriorCopy==0,:) = nan;
						matTempProbs(:,indAssignedTrials) = nan;
						
						%calculate probability of remaining trials and types
						[vecTempProbs,vecTempDecodedIndexCV]=max(matTempProbs,[],1);
						%get 2nd most likely stim per trial
						matDist2 = matTempProbs;
						for intT2=1:intTrials
							matDist2(vecTempDecodedIndexCV(intT2),intT2) = nan;
						end
						[vecTempProbs2,vecTempDecodedIndexCV2]=max(matDist2,[],1);
						
						%use trial with largest difference between most likely and 2nd most likely stimulus
						vecMaxDiff = abs(vecTempProbs2 - vecTempProbs);
						%assign trial
						[dummy,intAssignTrial]=max(vecMaxDiff);
						intAssignType = vecTempDecodedIndexCV(intAssignTrial);
						if vecPriorCopy(intAssignType) == 0
							intAssignType = vecTempDecodedIndexCV2(intAssignTrial);
						end
						vecDecodedIndexCV(intAssignTrial) = intAssignType;
						indAssignedTrials(intAssignTrial) = true;
						vecPriorCopy(intAssignType) = vecPriorCopy(intAssignType) - 1;
						%fprintf('assigned %d to %d; %s\n',intAssignType,intAssignTrial,sprintf('%d ',vecPriorCopy))
						%pause
					end
				end
				dblPerformanceX=sum(vecDecodedIndexCV(:) == vecTrialTypeIdx(:))/length(vecDecodedIndexCV);
				matAcrossTimeDecoder(intBinIdx,intTestBinIdx) = dblPerformanceX;
			end
		end
		%%
		%matDecPerf = matDecConf;
		
		dblChance = 1/numel(vecPriorDistribution);
		figure;maxfig;
		subplot(2,3,1)
		plot(vecStimTime,matDecPerf(:,end)');
		hold on
		plot([vecStimTime(1) vecStimTime(end)],[dblChance dblChance],'--','color',[0.5 0.5 0.5]);
		hold off
		title(sprintf('Dec perf; %s',strArea))
		xlabel('Time after onset (s)');
		ylabel('Fraction correct decoded');
		fixfig;
		
		subplot(2,3,2)
		plot(vecStimTime,vecSpikesPerBin./dblBinWidth)
		title('Spikes per bin')
		xlabel('Time after onset (s)');
		ylabel('Binned population spiking rate (Hz)');
		fixfig;
		
		subplot(2,3,3)
		plot(vecStimTime,(matDecPerf(:,end)./dblChance)'./(vecSpikesPerBin./dblBinWidth))
		title('Dec perf / spike')
		xlabel('Time after onset (s)');
		ylabel('Performance/spike');
		fixfig;
		
		hS=subplot(2,3,4);
		cMap=colormap(hS,blueredblue);
		hB=colorbar;
		set(gca,'clim',[min(vecStimTime) max(vecStimTime)]);
		h=cline([vecSpikesPerBin(:); vecSpikesPerBin(1)]./dblBinWidth,[matDecPerf(:,end); matDecPerf(1,end)],[vecStimTime(:); vecStimTime(1)]);
		set(h,'LineWidth',2);
		hold on
		plot([min(vecSpikesPerBin) max(vecSpikesPerBin)]./dblBinWidth,[dblChance dblChance],'--','color',[0.5 0.5 0.5]);
		hold off
		hB.Label.String = 'Time (s)';
		xlabel('Binned population spiking rate (Hz)');
		ylabel('Decoding performance');
		fixfig;
		
		hS2=subplot(2,3,5);
		colormap(hS2,parula);
		imagesc(vecStimTime,vecStimTime,matAcrossTimeDecoder);
		hB2=colorbar;
		hB2.Label.String = 'Decoding performance';
		axis xy
		xlabel('Training bin');
		ylabel('Testing bin');
		fixfig;grid off
		
		export_fig(fullpath(strFigurePath,sprintf('A1_PopActDynamics_%s.tif',strRec)));
		export_fig(fullpath(strFigurePath,sprintf('A1_PopActDynamics_%s.pdf',strRec)));
		
	end
end
toc
