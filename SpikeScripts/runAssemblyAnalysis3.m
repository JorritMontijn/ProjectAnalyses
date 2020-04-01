%load data
for intMouse=5%1:8
	close all
	drawnow;
	clearvars -except vecFrac intMouse
	
	%use neuropil subtraction?
	boolUseNeuropilSubtraction = false;
	boolOnlyTuned = false;
	
	%get block data
	if ~exist('cellMultiSes','var')
		if intMouse == 1
			strSes = '20140207';
		elseif intMouse == 2
			strSes = '20140314';
		elseif intMouse == 3
			strSes = '20140425';
		elseif intMouse == 4
			strSes = '20140507';
		elseif intMouse == 5
			strSes = '20140530';
		elseif intMouse == 6
			strSes = '20140604';
		elseif intMouse == 7
			strSes = '20140711';
		elseif intMouse == 8
			strSes = '20140715';
		end
		if boolUseNeuropilSubtraction
			strSes = ['NPS' strSes];
		end
		if boolOnlyTuned
			strSes = [strSes 'OT'];
		end
		fprintf('Loading pre-processed data for %s [%s]\n',strSes,getTime);
		load(['E:\UvA_Backup\Data\Results\spikeAnalysis\dataPreProAssemblies' strSes '.mat']);
	end
	
	clear sLoad sSesAggregate ses
	sParams.strFigDir = ['E:\UvA_Backup\Data\Results\spikeAnalysis' filesep strSes filesep];
	sParams.boolSavePlots = false;
	sParams.boolSaveData = false;
	strOldDir = cd(sParams.strFigDir);
	
	%change name for no split
	%#ok<*AGROW>
	strSes = ['AA' strSes]; 
	
	for intPopulation = 1:numel(cellMultiSes)
		%% GET RESPONSES AFTER COMBINING BLOCKS OF SAME POPULATION INTO ONE
		fprintf('Starting %s pop %d [%s]\n',strSes,intPopulation,getTime);
		
		%define stimulus duration parameter
		dblStimSecs = 2;
		
		% remove trials with reaction time <100ms
		dblRemSecs = 0.1;
		indTooFast = (cellMultiSes{intPopulation}.structStim.vecTrialRespSecs-cellMultiSes{intPopulation}.structStim.SecsOn)<dblRemSecs & cellMultiSes{intPopulation}.structStim.vecTrialResponse == 1;
		fprintf('Removed %d trials with responses <%dms\n',sum(indTooFast),round(dblRemSecs*1000));
		%structStim.vecTrialResponse(indTooFast) = 0;
		cellMultiSes{intPopulation}.structStim = remel(cellMultiSes{intPopulation}.structStim,~indTooFast);
		%structStim.FrameOff = structStim.FrameOn+1;
		
		%get stimulus variables
		sTypesC = getStimulusTypes(cellMultiSes{intPopulation}.structStim,{'Contrast'});
		cellSelectC = getSelectionVectors(cellMultiSes{intPopulation}.structStim,sTypesC);
		intContrasts = length(cellSelectC);
		
		sTypesO = getStimulusTypes(cellMultiSes{intPopulation}.structStim,{'Orientation'});
		cellSelectO = getSelectionVectors(cellMultiSes{intPopulation}.structStim,sTypesO);
		intOrientations = length(cellSelectO);
		
		sTypesCO = getStimulusTypes(cellMultiSes{intPopulation}.structStim,{'Orientation','Contrast'});
		cellSelectCO = getSelectionVectors(cellMultiSes{intPopulation}.structStim,sTypesCO);
		
		% get stimulus data
		vecOrientations = unique(cellMultiSes{intPopulation}.structStim.Orientation);
		intTrials = length(cellMultiSes{intPopulation}.structStim.Orientation);
		intOris = length(vecOrientations);
		[vecStimOris,dummy]=find(bsxfun(@eq,cellMultiSes{intPopulation}.structStim.Orientation,vecOrientations'));
		vecStimOris = vecStimOris';
		
		[vecStimContrasts,dummy]=find(bsxfun(@eq,cellMultiSes{intPopulation}.structStim.Contrast,unique(cellMultiSes{intPopulation}.structStim.Contrast)'));
		vecStimContrasts = vecStimContrasts';
		
		indStimResp = cellMultiSes{intPopulation}.structStim.vecTrialRespSecs- cellMultiSes{intPopulation}.structStim.SecsOn<3;
		
		
		%get behavioral variables
		indMiss = cellMultiSes{intPopulation}.structStim.vecTrialResponse == 0;
		indHit = cellMultiSes{intPopulation}.structStim.vecTrialResponse == 1;
		vecRTs = cellMultiSes{intPopulation}.structStim.vecTrialRespSecs -  cellMultiSes{intPopulation}.structStim.SecsOn;
		vecRTs(indMiss) = nan;
		indFast = false(size(indMiss));
		indSlow = false(size(indMiss));
		for intContrast=unique(vecStimContrasts)
			indTheseContrastTrials = vecStimContrasts==intContrast;
			indFast = (vecRTs <= nanmedian(vecRTs(indTheseContrastTrials)) & indTheseContrastTrials) | indFast;
			indSlow = (vecRTs > nanmedian(vecRTs(indTheseContrastTrials)) & indTheseContrastTrials) | indSlow;
		end
		
		%get assembly variables
		matAssemblyActivity = sAssemblies{intPopulation}.matAssemblyActivity; %[assemblies x frames]
		%matRespAssemblies = sAssemblies{intPopulation}.matRespAssemblies; 
		matReorderedAssemblyCorrs=sAssemblies{intPopulation}.matAssemblyCorrelations(sAssemblies{intPopulation}.vecReorder,sAssemblies{intPopulation}.vecReorder);
		matM = getBlockMeans(matReorderedAssemblyCorrs,sAssemblies{intPopulation}.vecAssemblyIdentity);
		
		%get assembly activity during trials
		dblSamplingFreq = cellMultiSes{intPopulation}.structStim.FrameOff(end)/cellMultiSes{intPopulation}.structStim.SecsOff(end);
		sStimTemp = cellMultiSes{intPopulation}.structStim;
		sStimTemp.FrameOff = sStimTemp.FrameOn + round(dblSamplingFreq*dblStimSecs);
		matRespAssemblies = getNeuronResponseRespMat(matAssemblyActivity,sStimTemp);%[assemblies x trials]
		
		%get response matrix
		[matTrialResponse,cellSelectContrasts] = getTrialResponseData(cellMultiSes{intPopulation},sStimTemp);
		
		%get pre-resp assemblies
		structParamsAss.intStopOffset = 0;
		structParamsAss.intStartOffset = round(-1*cellMultiSes{intPopulation}.samplingFreq);
		matPreRespAssemblies = getNeuronResponseRespMat(matAssemblyActivity,sStimTemp,[],[],structParamsAss);%[assemblies x trials]

		%detect spikes
		intNeurons = numel(cellMultiSes{intPopulation}.neuron);
		dblSamplingFreq = cellMultiSes{intPopulation}.structStim.FrameOff(end)/cellMultiSes{intPopulation}.structStim.SecsOff(end);
		dblTau = 0.5;
		intTotDur = length(cellMultiSes{intPopulation}.neuron(1).dFoF);
		dblTotDurSecs = intTotDur/dblSamplingFreq;
		if isfield(sAssemblies{intPopulation},'matSpikeCountsHD')
			matSpikeCountsHD = sAssemblies{intPopulation}.matSpikeCountsHD;
		else
			matSpikeCountsHD = nan(intNeurons,intTotDur);
			for intNeuron=1:intNeurons
				%msg
				fprintf('Performing detection for neuron %d/%d...\n',intNeuron,intNeurons)
				ptrTime = tic;
				
				%AE detection
				[apFrames, apSpikes, vecSpikes, expFit] = doDetectSpikes(cellMultiSes{intPopulation}.neuron(intNeuron).dFoF,dblSamplingFreq,dblTau);
				intTotSpikes = sum(apSpikes);
				
				%assign to neuron
				cellMultiSes{intPopulation}.neuron(intNeuron).apFrames = apFrames;
				cellMultiSes{intPopulation}.neuron(intNeuron).apSpikes = apSpikes;
				cellMultiSes{intPopulation}.neuron(intNeuron).vecSpikes = vecSpikes;
				cellMultiSes{intPopulation}.neuron(intNeuron).expFit = expFit;
				matSpikeCountsHD(intNeuron,:) = vecSpikes;
				
				%msg
				fprintf('\b	Done! Took %.1f seconds; %d transients; mean rate is %.2f Hz\n',toc(ptrTime),intTotSpikes,intTotSpikes/dblTotDurSecs)
			end
		end
		matSpikeCounts = getNeuronResponseRespMat(matSpikeCountsHD,sStimTemp);
		matSpikeCounts = round(matSpikeCounts*dblSamplingFreq);
		
		%% create miss/slow/fast selection vectors
		intFrames = size(matSpikeCountsHD,2);
		%miss
		vecMissOn = cellMultiSes{intPopulation}.structStim.FrameOn(indMiss);
		vecMissOff = vecMissOn + round(dblSamplingFreq*dblStimSecs);
		indMissHD = false(1,intFrames);
		for intTrial=1:length(vecMissOn)
			indMissHD(vecMissOn(intTrial):vecMissOff(intTrial)) = true;
		end
		vecMissHD = find(indMissHD);
		%slow
		vecSlowOn = cellMultiSes{intPopulation}.structStim.FrameOn(indSlow);
		vecSlowOff = vecSlowOn + round(dblSamplingFreq*dblStimSecs);
		indSlowHD = false(1,intFrames);
		for intTrial=1:length(vecSlowOn)
			indSlowHD(vecSlowOn(intTrial):vecSlowOff(intTrial)) = true;
		end
		vecSlowHD = find(indSlowHD);
		%fast
		vecFastOn = cellMultiSes{intPopulation}.structStim.FrameOn(indFast);
		vecFastOff = vecFastOn + round(dblSamplingFreq*dblStimSecs);
		indFastHD = false(1,intFrames);
		for intTrial=1:length(vecFastOn)
			indFastHD(vecFastOn(intTrial):vecFastOff(intTrial)) = true;
		end
		vecFastHD = find(indFastHD);
		
		%% get members of assemblies
		intNeurons = size(sAssemblies{intPopulation}.matAssemblies,1);
		matAssemblyMembers = sAssemblies{intPopulation}.matAssemblies > 0;
		intAssemblies = max(sAssemblies{intPopulation}.vecAssemblies);
		matAssemblyCoreMembers = false(intAssemblies,intNeurons);
		for intAssembly=1:intAssemblies
			%get random memberships
			intRandIters = 1000;
			indAssembly = sAssemblies{intPopulation}.vecAssemblies == intAssembly;
			matRandMembership = nan(intNeurons,intRandIters);
			for intIter=1:intRandIters
				indRandomAssembly = indAssembly(randperm(length(indAssembly)));
				matRandMembership(:,intIter) = sum(matAssemblyMembers(:,indRandomAssembly),2)/sum(indAssembly);
			end
			dblAlpha = 0.01;
			dblCorrectedSDs = qnorm(1-(dblAlpha/2)/intNeurons);
			vecUpperBound = mean(matRandMembership,2)+std(matRandMembership,[],2)*dblCorrectedSDs;
			vecData = sum(matAssemblyMembers(:,indAssembly),2)/sum(indAssembly);
			indHigher = vecData>vecUpperBound;
			
			%save output
			matAssemblyCoreMembers(intAssembly,:) = indHigher;
		end
		%matAssemblyCoreMembers = true(size(matAssemblyCoreMembers));
		%get assembly types
		vecCoreMemberSize = sum(matAssemblyCoreMembers,2);
		vecAssemblyAutoCorr = diag(matM);
		indRecurringClusters = zscore(vecAssemblyAutoCorr)>-2;
		[dummy,intNonRecurringCluster]=min(zscore(vecAssemblyAutoCorr));
		indRecurringClusters(intNonRecurringCluster) = false;
		indMultiNeuronClusters = vecCoreMemberSize>1;
		indRealAssemblies = indRecurringClusters & indMultiNeuronClusters;
		%indRealAssemblies = true(size(indRealAssemblies));
		
		%% get assembly data
		vecAssemblyStarts = sAssemblies{intPopulation}.vecAssemblyStarts;
		vecAssemblyStops = sAssemblies{intPopulation}.vecAssemblyStops;
		vecAssemblies = sAssemblies{intPopulation}.vecAssemblies';
		vecRealAssemblies = ismember(vecAssemblies,find(indRealAssemblies));
		vecNonRecurringAssemblies = vecAssemblies == intNonRecurringCluster;
		vecAssemblyAEs = zeros(size(vecAssemblyStarts));
		dblSecsCCG = (max(vecAssemblyStops-vecAssemblyStarts)/dblSamplingFreq);
		vecPlotCCG = (-ceil(dblSamplingFreq*dblSecsCCG):ceil(dblSamplingFreq*dblSecsCCG));
		intMiddle = ceil(dblSamplingFreq)+1;
		
		%% perform decoding; are stimuli more accurately represented during population events?
		%vecAssemblyStarts
		%vecAssemblyStops
		indStimTrials = ~cellSelectContrasts{1};
		vecOriStimTrials = vecStimOris(indStimTrials);
		vecContrastStimTrials = vecStimContrasts(indStimTrials);
		matRespRealAssemblies = matRespAssemblies(indRealAssemblies,:);
		matPreRespRealAssemblies  = matPreRespAssemblies(indRealAssemblies,:);
		matRespNonRecurring =  matRespAssemblies(~indRealAssemblies,:);
		matPreRespNRAssemblies = matPreRespAssemblies(~indRealAssemblies,:);
		intStimTrials = length(vecOriStimTrials);
		
		
		cellSelectSubO = cell(1,4);
		for i=1:length(cellSelectO)
			cellSelectSubO{i} = cellSelectO{i}(indStimTrials);
		end
		
		
		%decode orientation with cross-validated template matching algorithm
		structStimDecode = remel(sStimTemp,indStimTrials);
		sMetaOut_dFoF = doCrossValidationRespMat(matTrialResponse(:,indStimTrials),structStimDecode);
		dblPerformanceML_dFoF = sum(sMetaOut_dFoF.vecDecodedStimType==sMetaOut_dFoF.vecStimType)/sum(indStimTrials);
		sMetaOut_AEs = doCrossValidationRespMat(matSpikeCounts(:,indStimTrials),structStimDecode);
		dblPerformanceML_AEs = sum(sMetaOut_AEs.vecDecodedStimType==sMetaOut_AEs.vecStimType)/sum(indStimTrials);
		[phat,dblCI_MLdFoF] = binofit(dblPerformanceML_dFoF*intStimTrials,intStimTrials);
		[phat,dblCI_MLAEs] = binofit(dblPerformanceML_AEs*intStimTrials,intStimTrials);
		

		% decode orientation with template matching
		[dblPerformance_dFoF,vecDecodedIndexCV,matDists] = doCrossValidatedTemplateMatching(matTrialResponse(:,indStimTrials),vecOriStimTrials);
		[dblPerformance_AEs,vecDecodedIndexCV_AEs_TM,matDists] = doCrossValidatedTemplateMatching(matSpikeCounts(:,indStimTrials),vecOriStimTrials);
		[phat,dblCI_dFoF] = binofit(dblPerformance_dFoF*intStimTrials,intStimTrials);
		[phat,dblCI_AEs] = binofit(dblPerformance_AEs*intStimTrials,intStimTrials);
		
		%get neuron response properties
		structStimCorrsNAE = calcStimCorrsRespMat(matSpikeCounts(:,indStimTrials),cellSelectSubO);
		sTuningNAE = calcTuningRespMat(matSpikeCounts(:,indStimTrials),cellSelectSubO,vecOrientations);
		
		%get assembly response properties
		structStimCorrsAss = calcStimCorrsRespMat(matRespRealAssemblies(:,indStimTrials),cellSelectSubO);
		sTuningAss = calcTuningRespMat(matRespRealAssemblies(:,indStimTrials),cellSelectSubO,vecOrientations);
		
		%plot decoding
		figure
		%TM
		subplot(2,2,1)
		vecDecodingPerformance = [dblPerformance_dFoF dblPerformance_AEs];
		errorbar([0.3 0.7],vecDecodingPerformance,[dblCI_dFoF(1) dblCI_AEs(1)]-vecDecodingPerformance,[dblCI_dFoF(2) dblCI_AEs(2)]-vecDecodingPerformance,'xb')
		xlim([0 1])
		ylim([0 1])
		hold on
		plot(get(gca,'xlim'),[1 1]/length(cellSelectSubO),'k--')
		set(gca,'XTick',[0.3 0.7],'xticklabel',{'Neuron dF/F0','Neuron AE'})
		ylabel('Orientation decoding accuracy (TM)')
		
		%ML
		subplot(2,2,2)
		vecDecodingPerformance = [dblPerformanceML_dFoF dblPerformanceML_AEs];
		errorbar([0.3 0.7],vecDecodingPerformance,[dblCI_MLdFoF(1) dblCI_MLAEs(1)]-vecDecodingPerformance,[dblCI_MLdFoF(2) dblCI_MLAEs(2)]-vecDecodingPerformance,'xb')
		xlim([0 1])
		ylim([0 1])
		hold on
		plot(get(gca,'xlim'),[1 1]/length(cellSelectSubO),'k--')
		set(gca,'XTick',[0.3 0.7],'xticklabel',{'Neuron dF/F0','Neuron AE'})
		ylabel('Orientation decoding accuracy (TM)')
		
		
		%% collect ordering of all neurons in each population event during miss/slow/fast trials, then check if ordering is more consistent during fast than slow trial assembly occurrences
		cellNeuronSequencesAssemblies = cell(intAssemblies,3);
		cellNeuronSequences = cellfill(nan(intNeurons,length(vecAssemblies)),[1 4]);
		intMissCounterOverall = 0;
		intSlowCounterOverall = 0;
		intFastCounterOverall = 0;
		intAllCounterOverall = 0;
		vecRunAssemblies = 1:intAssemblies;%find(indRealAssemblies)';%1:intAssemblies
		for intAssembly=vecRunAssemblies
			indOccurrences = vecAssemblies == intAssembly;
			vecStarts = vecAssemblyStarts(indOccurrences);
			vecStops = vecAssemblyStops(indOccurrences);
			intOccurrences = numel(vecStarts);
			
			%pre-allocate these matrices
			cellNeuronSequencesAssemblies{intAssembly,1} = nan(intNeurons,intOccurrences);
			cellNeuronSequencesAssemblies{intAssembly,2} = nan(intNeurons,intOccurrences);
			cellNeuronSequencesAssemblies{intAssembly,3} = nan(intNeurons,intOccurrences);
			intMissCounterAssembly = 0;
			intSlowCounterAssembly = 0;
			intFastCounterAssembly = 0;
			
			%get sequences
			for intOccurrence=1:intOccurrences
				%get spike times for this occurrence
				intStart = vecStarts(intOccurrence);
				intStop = vecStops(intOccurrence)-1;
				matAEs = matSpikeCountsHD(:,intStart:intStop);
				
				%get mean and std
				matVals = repmat(1:size(matAEs,2),[size(matAEs,1) 1]);
				vecVals = getValByIdx(matVals,matAEs);
				dblMean = mean(vecVals);
				dblSD = std(vecVals);
				
				%check if during trial & get relative timing to center of mass of this occurrence
				if any(ismember(intStart:intStop,vecMissHD))
					intMissCounterAssembly = intMissCounterAssembly + 1;
					intMissCounterOverall = intMissCounterOverall + 1;
					
					vecThisMuCC = (sum((matAEs.*repmat(1:size(matAEs,2),[size(matAEs,1) 1])),2) ./ sum(matAEs,2) - dblMean)/dblSD;
					cellNeuronSequencesAssemblies{intAssembly,1}(:,intMissCounterAssembly) = vecThisMuCC;
					cellNeuronSequences{1}(:,intMissCounterOverall) = vecThisMuCC;
				end
				if any(ismember(intStart:intStop,vecSlowHD))
					intSlowCounterAssembly = intSlowCounterAssembly + 1;
					intSlowCounterOverall = intSlowCounterOverall + 1;
					
					vecThisMuCC = (sum((matAEs.*repmat(1:size(matAEs,2),[size(matAEs,1) 1])),2) ./ sum(matAEs,2) - dblMean)/dblSD;
					cellNeuronSequencesAssemblies{intAssembly,2}(:,intSlowCounterAssembly) = vecThisMuCC;
					cellNeuronSequences{2}(:,intSlowCounterOverall) = vecThisMuCC;
				end
				if any(ismember(intStart:intStop,vecFastHD))
					intFastCounterAssembly = intFastCounterAssembly + 1;
					intFastCounterOverall = intFastCounterOverall + 1;
					
					vecThisMuCC = (sum((matAEs.*repmat(1:size(matAEs,2),[size(matAEs,1) 1])),2) ./ sum(matAEs,2) - dblMean)/dblSD;
					cellNeuronSequencesAssemblies{intAssembly,3}(:,intFastCounterAssembly) = vecThisMuCC;
					cellNeuronSequences{3}(:,intFastCounterOverall) = vecThisMuCC;
				end
				
				%get overall sd
				vecThisMuCC = (sum((matAEs.*repmat(1:size(matAEs,2),[size(matAEs,1) 1])),2) ./ sum(matAEs,2) - dblMean)/dblSD;
				cellNeuronSequences{4}(:,intOccurrence) = vecThisMuCC;
				intAllCounterOverall = intAllCounterOverall + 1;
			end
			
			%remove pre-allocated nans
			if intMissCounterAssembly == 0,intMissCounterAssembly = 1;end
			if intSlowCounterAssembly == 0,intSlowCounterAssembly = 1;end
			if intFastCounterAssembly == 0,intFastCounterAssembly = 1;end
			cellNeuronSequencesAssemblies{intAssembly,1}(:,(intMissCounterAssembly+1):end) = [];
			cellNeuronSequencesAssemblies{intAssembly,2}(:,(intSlowCounterAssembly+1):end) = [];
			cellNeuronSequencesAssemblies{intAssembly,3}(:,(intFastCounterAssembly+1):end) = [];
		end
		%remove pre-allocated nans
		cellNeuronSequences{1}(:,(intMissCounterOverall+1):end) = [];
		cellNeuronSequences{2}(:,(intSlowCounterOverall+1):end) = [];
		cellNeuronSequences{3}(:,(intFastCounterOverall+1):end) = [];
		cellNeuronSequences{4}(:,(intAllCounterOverall+1):end) = [];
		
		%calculate mean+sd
		vecMeanLatencyMiss = nanmean(cellNeuronSequences{1},2);
		vecMeanLatencySlow = nanmean(cellNeuronSequences{2},2);
		vecMeanLatencyFast = nanmean(cellNeuronSequences{3},2);
		vecMeanLatencyOverall = nanmean(cellNeuronSequences{4},2);
		vecSDLatencyMiss = nanstd(cellNeuronSequences{1},[],2);
		vecSDLatencySlow = nanstd(cellNeuronSequences{2},[],2);
		vecSDLatencyFast = nanstd(cellNeuronSequences{3},[],2);
		vecSDLatencyOverall = nanstd(cellNeuronSequences{4},[],2);
		
		
		%select neurons with sufficient data points
		intMin = 5;
		indSelect1 = (sum(~isnan(cellNeuronSequences{1}),2)>intMin);
		indSelect2 = (sum(~isnan(cellNeuronSequences{2}),2)>intMin);
		indSelect3 = (sum(~isnan(cellNeuronSequences{3}),2)>intMin);
		indSelect4 = (sum(~isnan(cellNeuronSequences{4}),2)>intMin);
		
		vecMeanLatencyMiss(~indSelect1) = [];
		vecMeanLatencySlow(~indSelect2) = [];
		vecMeanLatencyFast(~indSelect3) = [];
		vecMeanLatencyOverall(~indSelect4) = [];
		vecSDLatencyMiss(~indSelect1) = [];
		vecSDLatencySlow(~indSelect2) = [];
		vecSDLatencyFast(~indSelect3) = [];
		vecSDLatencyOverall(~indSelect4) = [];
		
		
		%{
		%plot
		figure
		subplot(2,2,1)
		scatter(vecMeanLatencyMiss,vecSDLatencyMiss,'r')
		
		subplot(2,2,2)
		scatter(vecMeanLatencySlow,vecSDLatencySlow,'m')
		
		subplot(2,2,3)
		scatter(vecMeanLatencyFast,vecSDLatencyFast,'g')
		
		subplot(2,2,4)
		scatter(vecMeanLatencyFast,vecSDLatencyFast,'k')
		%}
		
		%plot
		hSequenceStability = figure;
		
		matLatenciesMiss = cellNeuronSequences{1};
		vecLatenciesMiss = nanmean(matLatenciesMiss,2);
		matLatenciesHit = [cellNeuronSequences{2} cellNeuronSequences{3}];
		vecLatenciesHit = nanmean(matLatenciesHit,2);
		vecLatencies = mean(cat(2,vecLatenciesMiss,vecLatenciesHit),2);
		[dummy,vecReorder] = sort(vecLatencies);
		
		subplot(2,2,1);
		
		vecLimMiss = 2*(nanmean(matLatenciesMiss(:))+nanstd(matLatenciesMiss(:)))*[-1 1];
		matColormap = colormap('redbluepurple');
		implot(matLatenciesMiss(vecReorder,:),vecLimMiss,matColormap,3);
		%colorbar;
		title('Neuronal response latency relative to pop. CoM [Miss]');
		ylabel('Neuron')
		xlabel('Population event during miss trial')
		
		subplot(2,2,2);
		
		
		vecLimHit = 2*(nanmean(matLatenciesHit(:))+nanstd(matLatenciesHit(:)))*[-1 1];
		matColormap = colormap('redbluepurple');
		implot(matLatenciesHit(vecReorder,:),vecLimHit,matColormap,3);
		%colorbar;
		title('Neuronal response latency relative to pop. CoM [Hit]');
		ylabel('Neuron')
		xlabel('Population event during hit trial')
		
		
		subplot(2,2,3);
		errorbar(1,nanmean(vecSDLatencyMiss),nanstd(vecSDLatencyMiss)/sqrt(sum(~isnan(vecSDLatencyMiss))),'xr');
		hold on;
		errorbar(2,nanmean(vecSDLatencySlow),nanstd(vecSDLatencySlow)/sqrt(sum(~isnan(vecSDLatencySlow))),'xm');
		errorbar(3,nanmean(vecSDLatencyFast),nanstd(vecSDLatencyFast)/sqrt(sum(~isnan(vecSDLatencyFast))),'xg');
		hold off;
		%ylim([0 max(get(gca,'ylim'))]);
		set(gca,'xtick',[1 2 3],'xticklabel',{'Miss','Slow','Fast'});
		ylabel('Temporal sequence variability (s)');
		title('Mean consistency of activation order in population events across neurons');
		
		subplot(2,2,4);
		errorbar(1,nanmean(vecMeanLatencyMiss),nanstd(vecMeanLatencyMiss)/sqrt(sum(~isnan(vecMeanLatencyMiss))),'xr');
		hold on;
		errorbar(2,nanmean(vecMeanLatencySlow),nanstd(vecMeanLatencySlow)/sqrt(sum(~isnan(vecMeanLatencySlow))),'xm');
		errorbar(3,nanmean(vecMeanLatencyFast),nanstd(vecMeanLatencyFast)/sqrt(sum(~isnan(vecMeanLatencyFast))),'xg');
		hold off;
		%ylim([0 max(get(gca,'ylim'))]);
		set(gca,'xtick',[1 2 3],'xticklabel',{'Miss','Slow','Fast'});
		ylabel('Temporal sequence variability (s)');
		title('Mean consistency of activation order in population events across neurons');
		
		%save figure
		if sParams.boolSavePlots
			drawnow;
			jFig = get(handle(gcf), 'JavaFrame');
			jFig.setMaximized(true);
			figure(gcf);
			drawnow;
			export_fig(sprintf('%spop%dtemporal_sequence_stability.tif',strSes,intPopulation));
			export_fig(sprintf('%spop%dtemporal_sequence_stability.pdf',strSes,intPopulation));
		end
		
		%save data
		cellSaveTemporalSequenceStability{intPopulation} = [nanmean(vecSDLatencyMiss) nan nanmean(vecSDLatencySlow) nanmean(vecSDLatencyFast) nanmean(vecSDLatencyOverall)]; %correlation temporal order within assemblies with temporal order overall
		
		%% get assembly location data
		matLocationsXY = nan(intNeurons,2);
		for intNeuron = 1:intNeurons
			matLocationsXY(intNeuron,1) = cellMultiSes{intPopulation}.neuron(intNeuron).x;
			matLocationsXY(intNeuron,2) = cellMultiSes{intPopulation}.neuron(intNeuron).y;
		end
		dblPix2Micron = cellMultiSes{intPopulation}.xml.dblActualVoxelSizeX/1000;
		matLocationsXY = matLocationsXY*dblPix2Micron; %transform pixels to microns
		matDiffX = abs(bsxfun(@minus,matLocationsXY(:,1),matLocationsXY(:,1)'));
		matDiffY = abs(bsxfun(@minus,matLocationsXY(:,2),matLocationsXY(:,2)'));
		matDist = sqrt(matDiffX.^2+matDiffY.^2);
		
		%compare mean distance of core assembly members to random group of same size
		intIters = 1000;
		vecMeanDist = nan(1,intAssemblies);
		matShuffledDist = nan(intIters,intAssemblies);
		vecPotentialMembers = find(any(matAssemblyCoreMembers,1));
		intPotentialMembers = length(vecPotentialMembers);
		for intAssembly=1:intAssemblies
			%get members
			vecMembers = matAssemblyCoreMembers(intAssembly,:);
			intMembers = sum(vecMembers);
			if intMembers < 2,continue;end
			matSelect = tril(true(intMembers),-1);
			
			%get real dPA
			matRealDist = matDist(vecMembers,vecMembers);
			vecMeanDist(intAssembly) = mean(matRealDist(matSelect));
			
			%get shuffled dPAs
			for intIter=1:intIters
				vecShuffledMembers = vecPotentialMembers(randperm(intPotentialMembers,intMembers));
				matShuffledDistTemp = matDist(vecShuffledMembers,vecShuffledMembers);
				matShuffledDist(intIter,intAssembly) = mean(matShuffledDistTemp(matSelect));
			end
		end
		
		%calc data
		dblAlpha = 0.05;
		intLowerThreshold = round(intIters*(dblAlpha/2));
		intUpperThreshold = round(intIters*(1-dblAlpha/2));
		matShuffledSorted = sort(matShuffledDist,1);
		vecUpper = matShuffledSorted(intUpperThreshold,:);
		vecLower = matShuffledSorted(intLowerThreshold,:);
		vecMean = mean(matShuffledSorted,1);
		
		% plot
		figure
		%subplot(2,2,1);
		errorbar(1:intAssemblies,vecMean,vecLower-vecMean,vecUpper-vecMean,'xb');
		hold on
		scatter(1:intAssemblies,vecMeanDist,'xr');
		hold off
		xlabel('Assembly');
		ylabel('Mean distance between members (micron)');
		
		if sParams.boolSavePlots
			drawnow;
			%jFig = get(handle(gcf), 'JavaFrame');
			%jFig.setMaximized(true);
			%figure(gcf);
			%drawnow;
			export_fig(sprintf('%spop%dassembly_member_distances.tif',strSes,intPopulation));
			export_fig(sprintf('%spop%dassembly_member_distances.pdf',strSes,intPopulation));
		end
		
		% save data
		cellSaveAssemblyDistances{intPopulation,1} = vecMeanDist;
		cellSaveAssemblyDistances{intPopulation,2} = matShuffledDist;
		
		
		%% OSI neurons vs OSI assemblies (distribution comparison)
		%get assembly OSIs
		indStimC = ~cellSelectC{1};%cellSelectC{6} | cellSelectC{5};
		matRespAssStim = matRespAssemblies(:,indStimC);
		cellSelectStimO = cell(size(cellSelectO));
		for intOri=1:length(cellSelectStimO)
			cellSelectStimO{intOri} = cellSelectO{intOri}(indStimC);
		end
		sTuning = calcTuningRespMat(matRespAssStim,cellSelectStimO,vecOrientations);
		vecAssOSI = sTuning.vecOSI;
		
		%get neuron OSIs
		matRespNeuronStim = matSpikeCounts(:,indStimC); %vecPotentialMembers
		sTuning = calcTuningRespMat(matRespNeuronStim,cellSelectStimO,vecOrientations);
		vecNeuronOSI = sTuning.vecOSI;
		
		% save data
		cellSaveOSIs{intPopulation,1} = vecAssOSI;
		cellSaveOSIs{intPopulation,2} = vecNeuronOSI;
	end
	
	%% save data structures
	if sParams.boolSaveData
		vecClock = fix(clock);
		strDate = num2str(vecClock(1:3));
		strFile = ['dataAssAnalThree' strSes '_' strrep(strDate,'    ','_')];
		save(['D:\Data\Results\spikeAnalysis\' strFile],'cellSave*','-v7.3');
		cd(strOldDir);
	end
end