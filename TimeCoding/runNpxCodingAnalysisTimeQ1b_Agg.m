%% aim
%{
compare mean pop activity and decoding of real, within-trial-randomized spike times, and shufftid data

%}
%% define qualifying areas
clear all;
runHeaderPopTimeCoding;
cellUseAreas = {...
	'Primary visual area',...
	...'posteromedial visual area',...
	};

%% select all neurons in LP and drifting grating stimuli
if ~exist('sAggStim','var') || isempty(sAggStim) || isempty(sAggNeuron)
	[sAggStim,sAggNeuron]=loadDataNpx('','driftinggrating',strDataPath);
end
indRemDBA = strcmpi({sAggNeuron.SubjectType},'DBA');
fprintf('Removing %d cells of DBA animals; %d remaining [%s]\n',sum(indRemDBA),sum(~indRemDBA),getTime);
sAggNeuron(indRemDBA) = [];

%% pre-allocate matrices
cellTypes = {'Real', 'UniformTrial', 'ShuffTid'};%, 'PoissGain'}; %PoissGain not done
intNumTypes = numel(cellTypes);
intAreas = numel(cellUseAreas);

%% go through recordings
tic
for intRec=1:numel(sAggStim) %19 || weird: 11
	% get matching recording data
	strRec = sAggStim(intRec).Exp;
	sThisRec = sAggStim(strcmpi(strRec,{sAggStim(:).Exp}));
	
	%prep grating data
	[sUseNeuron,vecStimOnTime,vecStimOffTime,vecOrientation] = NpxPrepGrating(sAggNeuron,sThisRec,cellUseAreas);
	if isempty(sUseNeuron),continue;end
	[vecOriIdx,vecUniqueOris,vecRepNum,cellSelect,vecTrialRepetition] = val2idx(vecOrientation);
	intTrialNum = numel(vecStimOnTime);
	intStimNum = numel(unique(vecOrientation));
	intRepNum = intTrialNum/intStimNum;
	
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
			strType = cellTypes{intType};
			for intN=1:intNumN
				vecRepCounter = zeros(1,intOriNum);
				[vecTrialPerSpike,vecTimePerSpike] = getSpikesInTrial(cellSpikeTimes{intN},vecStimOnTime-dblPreTime,dblMaxDur);
				vecRandTrials = randperm(intTrialNum);
				for intTrial=1:intTrialNum
					if strcmp(strType,'Real')
						%do nothing
						vecSpikeT = vecTimePerSpike(vecTrialPerSpike==intTrial);
					elseif strcmp(strType,'UniformTrial')
						%make spike times uniform in trial
						vecSpikeT = vecTimePerSpike(vecTrialPerSpike==intTrial);
						vecSpikeT = sort(rand(size(vecSpikeT))*range(vecSpikeT)+min(vecSpikeT));
					elseif strcmp(strType,'ShuffTid')
						%shuffle trial ids for each neuron independently
						vecSpikeT = vecTimePerSpike(vecTrialPerSpike==vecRandTrials(intTrial));
					elseif strcmp(strType,'PoissGain')
						error not done yet
					end
					
					intTrialOriIdx = vecOriIdx(intTrial);
					vecRepCounter(intTrialOriIdx) = vecRepCounter(intTrialOriIdx) + 1;
					intRep = vecRepCounter(intTrialOriIdx);
					vecSpikeC = histcounts(vecSpikeT,vecBinEdges);
					if any(isnan(vecSpikeC))
						error
					end
					matBNSR(:,intN,intTrialOriIdx,intRep) = vecSpikeC;
					matBNT(:,intN,intTrial) = vecSpikeC;
				end
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
				intVerbose = 0;
				[dblPerformanceCV,vecDecodedIndexCV,matPosteriorProbability,dblMeanErrorDegs,matConfusion,matWeights] = ...
					doCrossValidatedDecodingLR(matTrainData,vecOri180,intTypeCV,vecPriorDistribution,dblLambda,intVerbose);
				
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
			title(sprintf('Dec perf; %s',strType))
			xlabel('Time after onset (s)');
			ylabel('Fraction correct decoded');
			fixfig;
			
			subplot(2,3,2)
			plot(vecStimTime,vecSpikesPerBin./dblBinWidth)
			title('Spikes per bin')
			xlabel('Time after onset (s)');
			ylabel('Binned population spiking rate (Hz)');
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
			export_fig(fullpath(strFigurePathSR,sprintf('A1b_PopActDynamics_%s_%s.tif',strRec,strType)));
			export_fig(fullpath(strFigurePathSR,sprintf('A1b_PopActDynamics_%s_%s.pdf',strRec,strType)));
			
			%% save data
			save(fullpath(strTargetDataPath,sprintf('Q1bData_%s_%s',strRec,strType)),...
				'vecBinEdges',...
				'vecOri180',...
				'matAcrossTimeDecoder',...
				'vecSpikesPerBin',...
				'dblBinWidth',...
				'vecDecConf',...
				'vecDecCorr',...
				'vecStimTime');
		end
	end
	close all;
end
toc
