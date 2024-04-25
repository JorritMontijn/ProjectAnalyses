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
	[sAggStim,sAggNeuron]=loadDataNpx('','driftinggrating',strDataPath);
end
indRemDBA = strcmpi({sAggNeuron.SubjectType},'DBA');
fprintf('Removing %d cells of DBA animals; %d remaining [%s]\n',sum(indRemDBA),sum(~indRemDBA),getTime);
sAggNeuron(indRemDBA) = [];

%% pre-allocate matrices
cellTypes = {'Real', 'UniformTrial', 'ShuffTid', 'PoissGain'};
intNumTypes = numel(cellTypes);
intAreas = numel(cellUseAreas);
matR_PP_All = nan(intAreas,intAreas,numel(sAggStim),2,intNumTypes);
matR_OP_All = nan(intAreas,intAreas,numel(sAggStim),2,intNumTypes);
matR_OO_All = nan(intAreas,intAreas,numel(sAggStim),2,intNumTypes);
matR_PO_All = nan(intAreas,intAreas,numel(sAggStim),2,intNumTypes);
mat_xR_All = nan(intAreas,intAreas,numel(sAggStim),2,intNumTypes);

%% go through recordings
tic
for intRec=1:numel(sAggStim) %19 || weird: 11
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
		
		%% prep data
		%get data matrix
		cellSpikeTimesRaw = {sArea1Neurons.SpikeTimes};
		[matData,indTuned,cellSpikeTimes,sOut,cellSpikeTimesPerCellPerTrial,vecStimOnStitched,vecNonStat,dblBC,dblMaxDevFrac] = ...
			NpxPrepData(cellSpikeTimesRaw,vecStimOnTime,vecStimOffTime,vecOrientation);
		intTunedN = sum(indTuned);
		intNumN = size(matData,1);
		
		%get ori vars
		vecOri180 = mod(vecOrientation,180)*2;
		[vecOriIdx,vecUnique,vecPriorDistribution,cellSelect,vecRepetition] = val2idx(vecOri180);
		intOriNum = numel(vecUnique);
		intRepNum = min(vecPriorDistribution);
		dblStimDur = roundi(median(vecStimOffTime - vecStimOnTime),1,'ceil');
		dblBinWidth = 5/1000;%/32
		dblPreTime = 0.3;%10*dblBinWidth;
		dblPostTime = 0.2;%30*dblBinWidth;
		dblMaxDur = dblStimDur+dblPreTime+dblPostTime;
		%dblMaxDur = dblPreTime+dblPostTime;
		vecBinEdges = 0:dblBinWidth:dblMaxDur;
		vecStimTime = vecBinEdges(2:end)-dblBinWidth/2 - dblPreTime;
		indStimBins = vecStimTime > 0 & vecStimTime < dblStimDur;
		intBinNum = numel(vecBinEdges)-1;
		matBNSR = nan(intBinNum,intNumN,intOriNum,intRepNum);
		matBNT = nan(intBinNum,intNumN,intTrialNum);
		%matBNT_shifted = nan(intBinNum,intNumN,intTrialNum);
		
		%types: Real, UniformTrial, ShuffTid, PoissGain
		for intType=1:numel(cellTypes)
			matData = getSpikeCounts(cellSpikeTimes,vecStimOnTime,dblStimDur);
			for intN=1:intNumN
				vecRepCounter = zeros(1,intOriNum);
				[vecTrialPerSpike,vecTimePerSpike] = getSpikesInTrial(cellSpikeTimes{intN},vecStimOnTime-dblPreTime,dblMaxDur);
				for intTrial=1:intTrialNum
					intTrialOriIdx = vecOriIdx(intTrial);
					vecRepCounter(intTrialOriIdx) = vecRepCounter(intTrialOriIdx) + 1;
					intRep = vecRepCounter(intTrialOriIdx);
					vecSpikeT = vecTimePerSpike(vecTrialPerSpike==intTrial);
					vecSpikeC = histcounts(vecSpikeT,vecBinEdges);
					if any(isnan(vecSpikeC))
						error
					end
					matBNSR(:,intN,intTrialOriIdx,intRep) = vecSpikeC;
					matBNT(:,intN,intTrial) = vecSpikeC;
					%matBNT_shifted(:,intN,intTrial) = vecSpikeC;
					%matBNT_shifted(indStimBins,intN,intTrial) = circshift(matBNT_shifted(indStimBins,intN,intTrial),-round(intBinNum*vecDelayTimeBy(intTrial)));
				end
				
				%plot
				if 0
					clf;
					subplot(4,1,2:4)
					matR = squeeze(matBNT(:,intN,:))'./dblBinWidth;
					imagesc(vecStimTime,1:intTrialNum,matR);
					xlabel('Time (s)');
					ylabel('Trial #');
					title(sprintf('%s - N%d',strRec,intN),'interpreter','none');
					%colormap(redwhite)
					
					subplot(4,1,1)
					plot(vecStimTime,mean(matR,1));
					xlabel('Time (s)');
					ylabel('Mean rate (spike count)');
					vecRate1 = matData(intN,:);
					vecRate2 = mean(matR(:,indStimBins),2);
					title(sprintf('%s - N%d; mean rate: %.3f - %.3f',strRec,intN,mean(vecRate1),mean(vecRate2)),'interpreter','none');
					pause
				end
				
				%close;
			end
			
			%% time progression
			matDecActPerTrial = nan(intBinNum,intTrialNum);
			matDecIdxPerTrial = nan(intBinNum,intTrialNum);
			matDecRealIdxPerTrial = nan(intBinNum,intTrialNum);
			matDecConfPerTrial = nan(intBinNum,intTrialNum);
			vecDecCorr = nan(intBinNum,1);
			vecDecConf = nan(intBinNum,1);
			dblLambda = 1;
			intTypeCV = 2;
			vecOri180 = mod(vecOrientation,180)*2;
			vecRepNum180 = vecRepNum(1:12)*2;
			matDecConfusion = nan(intOriNum,intOriNum,intBinNum);
			
			% the first bin of vecSpikesPerBin is not equal to last bin, primarily because of
			% concatenating two sets of 480 trials; the post-period of trial #480 and pre-period of trial
			% #481 are therefore completely different. Secondly, the ITI has jitter; it is closer to
			% 1.55s than 1.50s - therefore inducing a 0.05s shift if setting the window duration to 1.5s
			vecSpikesPerBin = mean(sum(matBNT,2),3)';
			matAcrossTimeDecoder = nan(intBinNum,intBinNum);
			
			[vecTrialTypeIdx,vecUnique,vecPriorDistribution,cellSelect,vecRepetition] = val2idx(vecOri180);
			for intBinIdx=1:intBinNum
				fprintf('%s (%d/%d): Now at bin %d/%d [%s]\n',strRec,intRec,numel(sAggStim),intBinIdx,intBinNum,getTime);
				%% rewrite this to use across-bin decoding
				matTrainData = squeeze(matBNT(intBinIdx,:,:));
				
				%mvn doesn't work because it cannot handle zero-variance predictors
				%[dblPerformanceCV,vecDecodedIndexCV,matPosteriorProbability,dblMeanErrorDegs,matConfusion] = ...
				%	doCrossValidatedDecoding(matTrainData,vecOri180,intTypeCV,vecPriorDistribution,dblLambda);
				
				%do old logistic regression
				[dblPerformanceCV,vecDecodedIndexCV,matPosteriorProbability,dblMeanErrorDegs,matConfusion,matWeights] = ...
					doCrossValidatedDecodingLR(matTrainData,vecOri180,intTypeCV,vecPriorDistribution,dblLambda);
				
				vecDecCorr(intBinIdx) = dblPerformanceCV;
				vecConf = nan(size(vecDecodedIndexCV));
				for intTrial=1:numel(vecTrialTypeIdx)
					vecConf(intTrial) = matPosteriorProbability(vecTrialTypeIdx(intTrial),intTrial);
				end
				vecDecConf(intBinIdx) = mean(vecConf);
				matDecConfusion(:,:,intBinIdx) = matConfusion;
				
				%trial vectors
				matDecActPerTrial(intBinIdx,:)=squeeze(sum(matBNT(intBinIdx,:,:),2));
				matDecIdxPerTrial(intBinIdx,:)=vecDecodedIndexCV;
				matDecRealIdxPerTrial(intBinIdx,:)=vecTrialTypeIdx;
				matDecConfPerTrial(intBinIdx,:)=vecConf;
				
				%% apply on all bins
				for intTestBinIdx=1:intBinNum
					if intTestBinIdx == intBinIdx
						matAcrossTimeDecoder(intBinIdx,intTestBinIdx) =  mean(vecConf);
						continue;
					end
					matTestData = squeeze(matBNT(intTestBinIdx,:,:));
					
					% %do cross-bin decoding without sample splitting, as CV is automatic
					% matTestData = squeeze(matBNT(intTestBinIdx,:,:));
					% matCrossPosterior = doMvnDec(matTrainData,vecTrialTypeIdx,matTestData,dblLambda);
					% vecPriorCopy = vecPriorDistribution;
					
					%do cross-bin decoding without sample splitting, as CV is automatic
					matDataPlusLin = [matTestData; ones(1,size(matTestData,2))];
					matActivation = matWeights'*matDataPlusLin;
					matCrossPosterior = exp(bsxfun(@minus,matActivation,logsumexp(matActivation,1))); %softmax
					%matCrossPosterior = doMvnDec(matTrainData,vecTrialTypeIdx,matTestData,dblLambda);
					
					% normal decoding or with prior distro?
					vecPriorCopy = vecPriorDistribution;
					if isempty(vecPriorCopy)
						%calculate output
						[dummy, vecDecodedIndexCV] = max(matCrossPosterior,[],1);
					else
						vecDecodedIndexCV = doDecClassify(matCrossPosterior,vecPriorCopy);
					end
					dblPerformanceX=sum(vecDecodedIndexCV(:) == vecTrialTypeIdx(:))/length(vecDecodedIndexCV);
					matAcrossTimeDecoder(intBinIdx,intTestBinIdx) = dblPerformanceX;
				end
			end
			
			%%
			%matDecPerf = matDecConf;
			vecDecPerf = vecDecConf;%vecDecCorr
			
			dblChance = 1/numel(vecPriorDistribution);
			figure;maxfig;
			subplot(2,3,1)
			plot(vecStimTime,vecDecPerf(:,end)');
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
			plot(vecStimTime,(vecDecPerf(:,end)./dblChance)'./(vecSpikesPerBin./dblBinWidth))
			title('Dec perf / spike')
			xlabel('Time after onset (s)');
			ylabel('Performance/spike');
			fixfig;
			
			hS=subplot(2,3,4);
			cMap=colormap(hS,circcol);
			hB=colorbar;
			set(gca,'clim',[min(vecStimTime) max(vecStimTime)]);
			h=cline([vecSpikesPerBin(:); vecSpikesPerBin(1)]./dblBinWidth,[vecDecPerf(:,end); vecDecPerf(1,end)],[vecStimTime(:); vecStimTime(1)]);
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
			
			matNormAct = zscore(matDecActPerTrial,[],2);
			matCorrect = matDecIdxPerTrial == matDecRealIdxPerTrial;
			%matDecConfPerTrial(intBinIdx,:)=vecConf;
			
			%%
			export_fig(fullpath(strFigurePathSR,sprintf('A1_PopActDynamics_%s.tif',strRec)));
			export_fig(fullpath(strFigurePathSR,sprintf('A1_PopActDynamics_%s.pdf',strRec)));
			
		end
			%% save data
			save(fullpath(strTargetDataPath,sprintf('Q1Data_%s',strRec)),...
				'vecBinEdges',...
				'vecOri180',...
				'matAcrossTimeDecoder',...
				'vecSpikesPerBin',...
				'dblBinWidth',...
				'vecDecConf',...
				'vecDecCorr',...
				'vecStimTime');
		
	end
	close all;
end
toc
