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
intAreas = numel(cellUseAreas);
matR_PP_All = nan(intAreas,intAreas,numel(sAggStim),2);
matR_OP_All = nan(intAreas,intAreas,numel(sAggStim),2);
matR_OO_All = nan(intAreas,intAreas,numel(sAggStim),2);
matR_PO_All = nan(intAreas,intAreas,numel(sAggStim),2);
mat_xR_All = nan(intAreas,intAreas,numel(sAggStim),2);

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
		
		%% make example figure; decode & make quantiles
		% plot example single-neuron response to drifting gratings
		vecOriP = sOut.vecOriAnova;
		[a,intIdx]=min(vecOriP);
		matNeuronR = bsxfun(@rdivide,squeeze(matBNSR(:,intIdx,:,:)),diff(vecBinEdges(:)));
		vecNeuronTrialR = matData(intIdx,:);
		[matRespNSR,vecStimTypes] = getStimulusResponses(matData,vecTrialTypeIdx);
		matMeanResp = mean(matRespNSR,3);
		matSdResp = std(matRespNSR,[],3);
		vecMeanR = matMeanResp(intIdx,:);
		vecSemR = matSdResp(intIdx,:)/sqrt(intRepNum);
		[a,intPrefIdx]=max(vecMeanR);
		vecUseStims = [-1:1];
		intAggRepNum = numel(vecUseStims)*intRepNum;
		vecPrefIdx =  modx(intPrefIdx + vecUseStims,numel(vecMeanR));
		intOrthIdx = modx(intPrefIdx-numel(vecMeanR)/2,numel(vecMeanR));
		intObliIdx = modx(intPrefIdx+1,numel(vecMeanR));
		vecOrthIdx =  modx(intOrthIdx + vecUseStims,numel(vecMeanR));
		matPrefResp = reshape(matNeuronR(:,vecPrefIdx,:),[intBinNum intAggRepNum]);
		matOrthResp = reshape(matNeuronR(:,vecOrthIdx,:),[intBinNum intAggRepNum]);
		%matOrthResp = squeeze(matNeuronR(:,vecOrthIdx,:));
		
		vecColOrth = [0 0 0];
		vecColPref = lines(1);
		vecColObli = [0 0.3 0];%[0.8 0 0.8];
		vecColBckg = [0.6 0.6 0.6];
		figure;maxfig;
		h=subplot(2,3,1);
		hold on
		errorbar(vecStimTime,mean(matOrthResp,2),std(matOrthResp,[],2)/sqrt(intAggRepNum),'color',vecColOrth);
		errorbar(vecStimTime,mean(matPrefResp,2),std(matPrefResp,[],2)/sqrt(intAggRepNum),'color',vecColPref);
		hold off
		legend({'Orth','Pref'});
		xlabel('Time after stim onset (s)');
		ylabel('Spiking rate (Hz)');
		fixfig;
		
		subplot(2,3,2);
		vecPlotOris = vecUnique/2;
		errorbar([vecPlotOris 180],[vecMeanR vecMeanR(1)],[vecSemR vecSemR(1)],'color',vecColBckg);
		hold on
		errorbar(vecPlotOris(intObliIdx),vecMeanR(intObliIdx),vecSemR(intObliIdx),'color',vecColObli);
		hold off
		xlabel('Stimulus orientation (degs)')
		ylabel('Spiking rate (Hz)');
		xlim([-5 185]);
		set(gca,'xtick',0:45:180);
		ylim(gca,get(h,'ylim'));
		fixfig;grid off;
		
		subplot(2,3,3);
		dblSpikeStep = 2;
		vecSpikeE = 0:dblSpikeStep:20;
		vecSpikeC = vecSpikeE(2:end)-dblSpikeStep/2;
		matCounts = nan(intOriNum,numel(vecSpikeC));
		for intOri=1:intOriNum
			matCounts(intOri,:) = histcounts(vecNeuronTrialR(ismember(vecTrialTypeIdx,intOri)),vecSpikeE);
		end
		vecOrthR = vecNeuronTrialR(ismember(vecTrialTypeIdx,intOrthIdx));
		vecObliR = vecNeuronTrialR(ismember(vecTrialTypeIdx,intObliIdx));
		vecPrefR = vecNeuronTrialR(ismember(vecTrialTypeIdx,intPrefIdx));
		vecNotObliR = vecNeuronTrialR(~ismember(vecTrialTypeIdx,intObliIdx));
		dblObliqueDegs = vecUnique(intObliIdx)/2;
		vecCountsPref = histcounts(vecPrefR,vecSpikeE);
		vecCountsOrth = histcounts(vecOrthR,vecSpikeE);
		vecCountsObli = histcounts(vecObliR,vecSpikeE);
		vecCountsNotObli =  histcounts(vecNotObliR,vecSpikeE);
		hold on
		plot(vecSpikeC,matCounts(~ismember(1:intOriNum,intObliIdx),:)./sum(matCounts(~ismember(1:intOriNum,intObliIdx),:),2),'color',vecColBckg);
		%plot(vecSpikeC,vecCountsOrth/sum(vecCountsOrth),'color',vecColBckg);
		%plot(vecSpikeC,vecCountsNotObli/sum(vecCountsNotObli),'color',vecColBckg);
		plot(vecSpikeC,vecCountsObli/sum(vecCountsObli),'color',vecColObli);
		hold off
		xlabel('Spiking rate (Hz)');
		ylabel('Probability (# of trials)');
		%legend({sprintf('Not %d degs',dblObliqueDegs),sprintf('%d degs',dblObliqueDegs)});
		fixfig;
		
		
		subplot(2,3,4);
		vecRelProb = vecCountsObli./(max(matCounts(~ismember(1:intOriNum,intObliIdx),:),[],1));
		colormap(cat(1,[0 0 0],vecColObli));
		hold on
		plot(vecSpikeC([1 end]),[1 1],'-','color',vecColBckg);
		plot(vecSpikeC,vecRelProb,'color',vecColObli);
		%cline(vecSpikeC,vecRelProb,ones(size(vecRelProb)),double(vecRelProb>1)+1,true);
		%plot(mean(vecObliR)*[1 1],[0 max(get(gca,'ylim'))-0.01],'--','color',vecColObli);
		hold off
		xlabel('Spiking rate (Hz)');
		ylabel(sprintf('Ratio stim=%d/(most likely other stim)',dblObliqueDegs));
		fixfig;grid off
		
		%{
		[dblPerformanceCV,vecDecodedIndexCV,matPosteriorProbability,dblMeanErrorDegs,matConfusion,matWeights] = ...
			doCrossValidatedDecodingLR(matResp,vecOri180,intTypeCV,[],0);
		vecConfidence = nan(size(vecDecodedIndexCV));
		for intTrial=1:numel(vecTrialTypeIdx)
			vecConfidence(intTrial) = matPosteriorProbability(vecTrialTypeIdx(intTrial),intTrial);
		end
		
		subplot(2,3,5)
		colormap(redwhite(1024))
		imagesc(vecUnique/2,vecUnique/2,matConfusion);
		xlabel('Real orientation');
		ylabel('Decoded orientation');
		set(gca,'xtick',0:45:179);
		set(gca,'ytick',0:45:179);
		h=colorbar;
		ylabel(h,'# of trials');
		fixfig;grid off;
		
		%%
		% split trials into quantiles (per ori)
		vecPopRate = sum(matResp,1);
		intQuantiles = 3;
		intSplitTrialsPerOri = min(floor(cellfun(@sum,cellSelect)/intQuantiles));
		vecPriorDistributionSplit = ones(size(vecPriorDistribution))*intSplitTrialsPerOri;
		vecStartTrials = round(intRepNum*linspace(1/intRepNum,(1+1/intRepNum),intQuantiles+1));
		vecStartTrials(end)=[];
		vecTrialQuantile = zeros(1,intTrials);
		for intQ=1:intQuantiles
			matUseTrials = nan(intOriNum,intSplitTrialsPerOri);
			vecUseTrialsTemp = vecStartTrials(intQ):(vecStartTrials(intQ)+intSplitTrialsPerOri-1);
			for intOri=1:intOriNum
				vecThisOri = find(cellSelect{intOri});
				[vecSorted,vecReorder]=sort(vecPopRate(vecThisOri));
				matUseTrials(intOri,:) = vecThisOri(vecReorder(vecUseTrialsTemp));
			end
			vecUseTrials = sort(matUseTrials(:));
			vecTrialQuantile(vecUseTrials) = intQ;
		end
		
		%assign per quantile
		intTrialsPerQ = intSplitTrialsPerOri*intOriNum;
		matConfPerQ = nan(intTrialsPerQ,intQuantiles);
		matCorrPerQ = nan(intTrialsPerQ,intQuantiles);
		vecCorr = vecDecodedIndexCV(:) == vecTrialTypeIdx(:);
		for intQ=1:intQuantiles
			%divide into quantiles
			vecUseT = vecStartTrials(intQ):(vecStartTrials(intQ)+intTrialsPerQ-1);
			matCorrPerQ(:,intQ) = vecCorr(vecTrialQuantile==intQ);
			matConfPerQ(:,intQ) = vecConfidence(vecTrialQuantile==intQ);
		end
		
		dblAlphaEquivOfSd = normcdf(1)-normcdf(-1);
		[phat,pci] = binofit(sum(matCorrPerQ,1),size(matCorrPerQ,1),dblAlphaEquivOfSd);
		
		h2=subplot(2,3,6);
		mapC = redbluepurple(intQuantiles);
		hold on;
		%plot([1 intQuantiles],(1/intOriNum180)*[1 1],'--','color',vecColBckg);
		for intQ=1:intQuantiles
			errorbar(intQ,phat(intQ),phat(intQ)-pci(intQ,1)',phat(intQ)-pci(intQ,2)','x','color',mapC(intQ,:));
		end
		hold off;
		xlim([0.5 intQuantiles+0.5]);
		set(gca,'xtick',[1 ceil(intQuantiles/2) intQuantiles],'xticklabel',{'Low','Mid','High'});
		xlabel('Pop. act. quantile');
		ylabel('Decoding accuracy');
		fixfig;grid off;
		%}
		%%
		export_fig(fullpath(strFigurePathSR,sprintf('A2_ExampleNeuralCodes_%s.tif',strRec)));
		export_fig(fullpath(strFigurePathSR,sprintf('A2_ExampleNeuralCodes_%s.pdf',strRec)));
	
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
