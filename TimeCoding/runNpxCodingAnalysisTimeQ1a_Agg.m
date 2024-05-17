%% aim
%{
decode over time
%}

%% set parameters
cellDataTypes = {'Npx','Sim','ABI','SWN'};%topo, model, allen, nora
intRunDataType = 1;
strRunStim = 'DG';%DG or NM? => superseded to WS by SWN
cellTypes = {'Real', 'Uniform', 'ShuffTid', 'PoissGain'};
boolFixSpikeGroupSize = false;
dblRemOnset = 0; %remove onset period in seconds; 0.125 for sim, 0.25 for npx
runHeaderPopTimeCoding;

%% pre-allocate matrices
intNumTypes = numel(cellTypes);
intAreas = numel(cellUseAreas);
matR_PP_All = nan(intAreas,intAreas,numel(sAggStim),2,intNumTypes);
matR_OP_All = nan(intAreas,intAreas,numel(sAggStim),2,intNumTypes);
matR_OO_All = nan(intAreas,intAreas,numel(sAggStim),2,intNumTypes);
matR_PO_All = nan(intAreas,intAreas,numel(sAggStim),2,intNumTypes);
mat_xR_All = nan(intAreas,intAreas,numel(sAggStim),2,intNumTypes);

%% go through recs
tic
for intRec=6%1:intRecNum %19 || weird: 11
	%% prep ABI or Npx data
	if strcmp(strRunType,'ABI')
		error to be updated
		runRecPrepABI;
		strThisRec = strRec;
	elseif strcmp(strRunType,'Sim')
		%load
		runRecPrepSim;
		
		%edit vars
		strThisRec = strRec;
		strDataPathT0=strDataPathSimT0;
		vecOri180 = mod(vecOrientation,180)*2;
		vecStimIdx = vecOri180;
		
		%get cell props
		%vecNeuronPrefOri = pi-[sLoad.sNeuron.PrefOri];
		vecNeuronType = [sLoad.sNeuron.Types]; %1=pyr,2=interneuron
		%swap types to match npx
		vecSwap = [2 1];
		vecNeuronType = vecSwap(vecNeuronType);
		intUseMaxRep = 40;
		indRemTrials = vecTrialRepetition>intUseMaxRep;
		vecOri180(indRemTrials) = [];
		vecStimIdx(indRemTrials) = [];
		vecStimOnTime(indRemTrials) = [];
		vecStimOffTime(indRemTrials) = [];
		vecOrientation(indRemTrials) = [];
		dblMinHz = 0;
	elseif strcmp(strRunType,'Npx') || strcmp(strRunType,'SWN')
		%prep
		runRecPrepNpx;
		strThisRec = strRec;
		strDataPathT0 = strTargetDataPath;
		dblMinHz = 90;
		
		
		%% layers
		sArea1Neurons = sUseNeuron(indArea1Neurons);
		cellAreas = {sArea1Neurons.Area};
		vecCorticalLayer = cellfun(@(x) str2double(x(regexp(x,'layer.*')+6)),cellAreas);
		vecDepth = [sArea1Neurons.DepthBelowIntersect];
		vecSupraGranuInfra = double(vecCorticalLayer < 4) + 2*double(vecCorticalLayer == 4) + 3*double(vecCorticalLayer > 4);
		vecSupraGranuInfra = vecSupraGranuInfra(indResp);
	else
		error impossible
	end
	
	
	%% move onset
	%remove first x ms
	vecStimOnTime = vecStimOnTime + dblRemOnset;
	
	%% load prepped data
	strTarget = fullpath(strDataPathT0,[sprintf('T0Data_%s%s%s%s%s',strThisRec,strRunType,strRunStim,'Real') '.mat']);
	if ~exist(strTarget,'file')
		fprintf('Prepped T0 file did not exist for %s; skipping...\n',strThisRec);
		continue;
	else
		fprintf('Running %s... [%s]\n',strRec,getTime);
	end
	
	%check fr
	sSource = load(strTarget);
	cellSpikeTimes = sSource.cellSpikeTimes;
	intNumN = numel(cellSpikeTimes);
	cellSpikeTimesPerCellPerTrial = cell(intNumN,intTrialNum);
	for intN=1:intNumN
		%real
		[vecTrialPerSpike,vecTimePerSpike] = getSpikesInTrial(cellSpikeTimes{intN},vecStimOnTime,dblStimDur);
		for intTrial=1:numel(vecStimOnTime)
			vecSpikeT = vecTimePerSpike(vecTrialPerSpike==intTrial);
			cellSpikeTimesPerCellPerTrial{intN,intTrial} = vecSpikeT;
		end
	end
	matData = cellfun(@numel,cellSpikeTimesPerCellPerTrial)./dblStimDur;
	if mean(sum(matData)) < dblMinHz%90 / 50
		fprintf('Avg # of spikes per trial was %.1f for %s; skipping...\n',mean(sum(matData)),strThisRec);
		continue;
	end
	sTuning = getTuningCurves(matData,vecOri180,false);
	
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
		
		%% load prepped data and ifr
		sSource = load(fullpath(strDataPathT0,sprintf('T0Data_%s%s%s%s%s',strThisRec,strRunType,strRunStim,strType)));
		vecTime = sSource.vecTime;
		vecIFR = sSource.vecIFR;
		vecAllSpikeTime = sSource.vecAllSpikeTime;
		vecAllSpikeNeuron = sSource.vecAllSpikeNeuron;
		cellSpikeTimes = sSource.cellSpikeTimes;
		vecStimOnTime = sSource.vecStimOnTime;
		if numel(cellSpikeTimes) ~= intNumN
			error('neuron # mismatch!')
		end
		if isempty(vecTime),continue;end
		
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
		if 1
			%% example neuron
			dblPreTime= 0.3;
			intNeuron=25; %rec=6,n=25=sustained/39=onset+bursty
			dblBinWidth = 0.0125;
			dblPlotDur = dblMaxDur;
			vecBinEdges = 0:dblBinWidth:dblPlotDur;
			[duumy,intOriIdx]=max(sTuning.matMeanResp(intNeuron,:));
			vecBinC = vecBinEdges(2:end)-diff(vecBinEdges(1:2))/2-dblPreTime;
			matAct = squeeze(matBNSR(:,intNeuron,intOriIdx,:));
			vecStimT = vecStimOnTime(vecStimIdx==intOriIdx)-dblPreTime;
			vecTuneParams = sTuning.matFittedParams(intNeuron,:);
			matCol=lines(2);
			
			figure;maxfig;
			h1=subplot(2,3,1);
			
			%get spike times in trials
			[vecTrialPerSpike,vecTimePerSpike] = getSpikesInTrial(cellSpikeTimes{intNeuron},vecStimT,dblPlotDur);
			hold all;
			for intTrial=1:numel(vecStimT)
				vecTimes = vecTimePerSpike(vecTrialPerSpike==intTrial);
				vecTimes(vecTimes>dblPlotDur)=[];
				line([vecTimes(:)';vecTimes(:)'],intRepNum+[intTrial*ones(1,numel(vecTimes))-0.5;intTrial*ones(1,numel(vecTimes))+0.5],'Color',matCol(1,:),'LineWidth',1.5);
			end
			
			h2=subplot(2,3,2)
			[vecMean,vecSEM,vecWindowBinCenters,matPET] = doPEP(cellSpikeTimes{intNeuron},vecBinEdges,vecStimT);
			%errorbar(vecBinC,mean(matPET,1),std(matPET,[],1)./sqrt(intRepNum))
			plot(vecBinC,mean(matPET,1));%,std(matPET,[],1)./sqrt(intRepNum))
			
			h3=subplot(2,3,3)
			[vecTime,vecIFR,sIFR] = getIFR(cellSpikeTimes{intNeuron},vecStimT,1.5);
			plot(vecTime-dblPreTime,vecIFR)
			
			%orth
			intOrthIdx = modx(intOriIdx+round(intOriNum/2),intOriNum);
			vecStimT = vecStimOnTime(vecStimIdx==intOrthIdx)-dblPreTime;
			axes(h1);hold on;
			%get spike times in trials
			[vecTrialPerSpike,vecTimePerSpike] = getSpikesInTrial(cellSpikeTimes{intNeuron},vecStimT,dblPlotDur);
			for intTrial=1:numel(vecStimT)
				vecTimes = vecTimePerSpike(vecTrialPerSpike==intTrial);
				vecTimes(vecTimes>dblPlotDur)=[];
				line([vecTimes(:)';vecTimes(:)'],[intTrial*ones(1,numel(vecTimes))-0.5;intTrial*ones(1,numel(vecTimes))+0.5],'Color',matCol(2,:),'LineWidth',1.5);
			end
			hold off
			vecTicks = [0 dblPreTime dblStimDur+dblPreTime dblPlotDur];
			set(gca,'xtick',vecTicks,'xticklabels',vecTicks-dblPreTime);
			ylim([0 intRepNum*2+1]);
			xlim([0 dblPlotDur]);
			xlabel('Time (s)');
			ylabel('Trial #');
			title('blue=pref,orange=orth');
			
			axes(h2);hold on;
			[vecMean,vecSEM,vecWindowBinCenters,matPET] = doPEP(cellSpikeTimes{intNeuron},vecBinEdges,vecStimT);
			%errorbar(vecBinC,mean(matPET,1),std(matPET,[],1)./sqrt(intRepNum))
			plot(vecBinC,mean(matPET,1),'color',matCol(2,:));%,std(matPET,[],1)./sqrt(intRepNum))
			xlim([-dblPreTime dblMaxDur-dblPreTime]);
			set(gca,'xtick',vecTicks-dblPreTime);
			ylim([0 35]);
			xlabel('Time (s)');
			ylabel('Firing rate (Hz)');
			title(sprintf('%s - N%d',strRec,intNeuron),'interpreter','none');
			
			axes(h3);hold on;
			[vecTime,vecIFR,sIFR] = getIFR(cellSpikeTimes{intNeuron},vecStimT,1.5);
			plot(vecTime-dblPreTime,vecIFR,'color',matCol(2,:))
			ylim([0 35]);
			xlim([-dblPreTime dblMaxDur-dblPreTime]);
			set(gca,'xtick',vecTicks-dblPreTime);
			ylim([0 35]);
			xlabel('Time (s)');
			ylabel('IFR (Hz)');
			fixfig;
			
			%subplot(2,3,4)
			%imagesc(vecBinC,1:intRepNum,matAct');axis xy
			
			%subplot(2,3,5)
			%imagesc(vecBinC,1:intRepNum,matPET);axis xy
			
			if 0
				%% save
				export_fig(fullpath(strFigurePathSR,sprintf('Q1a_ExampleNeuron_%s%sN%d.tif',strRec,strType,intNeuron)));
				export_fig(fullpath(strFigurePathSR,sprintf('Q1a_ExampleNeuron_%s%sN%d.pdf',strRec,strType,intNeuron)));
		
			end
		end
		return
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
		intVerbose = 0;
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
		title(sprintf('%s - %s',strRec,strType))
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
		export_fig(fullpath(strFigurePathSR,sprintf('A1_PopActDynamics_%s%s.tif',strRec,strType)));
		export_fig(fullpath(strFigurePathSR,sprintf('A1_PopActDynamics_%s%s.pdf',strRec,strType)));
		
		
		%% save data
		save(fullpath(strTargetDataPath,sprintf('Q1Data_%s%s',strRec,strType)),...
			'vecBinEdges',...
			'vecOri180',...
			'matAcrossTimeDecoder',...
			'vecSpikesPerBin',...
			'dblBinWidth',...
			'vecDecConf',...
			'vecDecCorr',...
			'vecStimTime',...
			'matBNT',...
			'strRec',...
			'strType');
	end
end
close all;
toc
