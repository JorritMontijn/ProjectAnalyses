%load data
for intMouse=1:8
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
		indMissStim = indMiss & vecStimContrasts > 1;
		indCR = indMiss & vecStimContrasts == 1;
		
		indHits = cellMultiSes{intPopulation}.structStim.vecTrialResponse == 1;
		indHitsStim = indHits & vecStimContrasts > 1;
		indFA = indHits & vecStimContrasts == 1;
		
		vecRTs = cellMultiSes{intPopulation}.structStim.vecTrialRespSecs -  cellMultiSes{intPopulation}.structStim.SecsOn;
		vecRTs(indMiss) = nan;
		indFast = false(size(indMiss));
		indSlow = false(size(indMiss));
		for intContrast=unique(vecStimContrasts)
			indTheseContrastTrials = vecStimContrasts==intContrast;
			indFast = (vecRTs <= nanmedian(vecRTs(indTheseContrastTrials)) & indTheseContrastTrials) | indFast;
			indSlow = (vecRTs > nanmedian(vecRTs(indTheseContrastTrials)) & indTheseContrastTrials) | indSlow;
		end
		vecRespType = indMissStim * 1 + indHitsStim * 2 + indCR * 3 + indFA * 4;
		
		%get assembly variables
		matReorderedAssemblyCorrs=sAssemblies{intPopulation}.matAssemblyCorrelations(sAssemblies{intPopulation}.vecReorder,sAssemblies{intPopulation}.vecReorder);
		matM = getBlockMeans(matReorderedAssemblyCorrs,sAssemblies{intPopulation}.vecAssemblyIdentity);
		
		%assign random durations to miss trials
		sStimTemp = cellMultiSes{intPopulation}.structStim;
		%intMissDur = round(3 * cellMultiSes{intPopulation}.samplingFreq);
		%sStimTemp.FrameOff(indMiss) = sStimTemp.FrameOn(indMiss) + intMissDur;
		
		%get response matrix
		[matTrialResponse,cellSelectContrasts,vecTrialDur] = getTrialResponseData(cellMultiSes{intPopulation},sStimTemp);
		
		%detect spikes
		dblSamplingFreq = cellMultiSes{intPopulation}.structStim.FrameOff(end)/cellMultiSes{intPopulation}.structStim.SecsOff(end);
		matSpikeCountsHD = sAssemblies{intPopulation}.matSpikeCountsHD;
		matSpikeCounts = getNeuronResponseRespMat(matSpikeCountsHD,sStimTemp);
		matSpikeCounts = round(matSpikeCounts*dblSamplingFreq);
		
		%% create miss/slow/fast selection vectors
		intFrames = size(matSpikeCountsHD,2);
		%miss
		vecMissOn = sStimTemp.FrameOn(indMiss);
		vecMissOff = sStimTemp.FrameOff(indMiss);
		indMissHD = false(1,intFrames);
		for intTrial=1:length(vecMissOn)
			indMissHD(vecMissOn(intTrial):vecMissOff(intTrial)) = true;
		end
		vecMissHD = find(indMissHD);
		%hits
		vecHitsOn = sStimTemp.FrameOn(indHits);
		vecHitsOff = sStimTemp.FrameOff(indHits);
		indHitsHD = false(1,intFrames);
		for intTrial=1:length(vecHitsOn)
			indHitsHD(vecHitsOn(intTrial):vecHitsOff(intTrial)) = true;
		end
		vecHitsHD = find(indHitsHD);
		%base
		vecBaseOff = sStimTemp.FrameOn;
		vecBaseOn = sStimTemp.FrameOn - round(5*cellMultiSes{intPopulation}.samplingFreq);
		indBaseHD = false(1,intFrames);
		for intTrial=1:length(vecBaseOn)
			indBaseHD(vecBaseOn(intTrial):vecBaseOff(intTrial)) = true;
		end
		vecBaseHD = find(indBaseHD);
		
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
		%%{
		%% collect ordering of all neurons in each population event during miss/slow/fast trials, then check if ordering is more consistent during fast than slow trial assembly occurrences
		cellNeuronSequencesAssemblies = cell(intAssemblies,3);
		cellNeuronSequences = cellfill(nan(intNeurons,length(vecAssemblies)),[1 4]);
		intMissCounterOverall = 0;
		intHitsCounterOverall = 0;
		intBaseCounterOverall = 0;
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
			intHitsCounterAssembly = 0;
			intBaseCounterAssembly = 0;
			
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
				if any(ismember(intStart:intStop,vecHitsHD))
					intHitsCounterAssembly = intHitsCounterAssembly + 1;
					intHitsCounterOverall = intHitsCounterOverall + 1;
					
					vecThisMuCC = (sum((matAEs.*repmat(1:size(matAEs,2),[size(matAEs,1) 1])),2) ./ sum(matAEs,2) - dblMean)/dblSD;
					cellNeuronSequencesAssemblies{intAssembly,2}(:,intHitsCounterAssembly) = vecThisMuCC;
					cellNeuronSequences{2}(:,intHitsCounterOverall) = vecThisMuCC;
				end
				if any(ismember(intStart:intStop,vecBaseHD))
					intBaseCounterAssembly = intBaseCounterAssembly + 1;
					intBaseCounterOverall = intBaseCounterOverall + 1;
					
					vecThisMuCC = (sum((matAEs.*repmat(1:size(matAEs,2),[size(matAEs,1) 1])),2) ./ sum(matAEs,2) - dblMean)/dblSD;
					cellNeuronSequencesAssemblies{intAssembly,3}(:,intBaseCounterAssembly) = vecThisMuCC;
					cellNeuronSequences{3}(:,intBaseCounterOverall) = vecThisMuCC;
				end
				
				%get overall sd
				vecThisMuCC = (sum((matAEs.*repmat(1:size(matAEs,2),[size(matAEs,1) 1])),2) ./ sum(matAEs,2) - dblMean)/dblSD;
				cellNeuronSequences{4}(:,intOccurrence) = vecThisMuCC;
				intAllCounterOverall = intAllCounterOverall + 1;
			end
			
			%remove pre-allocated nans
			if intMissCounterAssembly == 0,intMissCounterAssembly = 1;end
			if intHitsCounterAssembly == 0,intHitsCounterAssembly = 1;end
			if intBaseCounterAssembly == 0,intBaseCounterAssembly = 1;end
			cellNeuronSequencesAssemblies{intAssembly,1}(:,(intMissCounterAssembly+1):end) = [];
			cellNeuronSequencesAssemblies{intAssembly,2}(:,(intHitsCounterAssembly+1):end) = [];
			cellNeuronSequencesAssemblies{intAssembly,3}(:,(intBaseCounterAssembly+1):end) = [];
		end
		%remove pre-allocated nans
		cellNeuronSequences{1}(:,(intMissCounterOverall+1):end) = [];
		cellNeuronSequences{2}(:,(intHitsCounterOverall+1):end) = [];
		cellNeuronSequences{3}(:,(intBaseCounterOverall+1):end) = [];
		cellNeuronSequences{4}(:,(intAllCounterOverall+1):end) = [];
		
		%calculate mean+sd
		vecMeanLatencyMiss = nanmean(cellNeuronSequences{1},2);
		vecMeanLatencyHits = nanmean(cellNeuronSequences{2},2);
		vecMeanLatencyBase = nanmean(cellNeuronSequences{3},2);
		vecMeanLatencyOverall = nanmean(cellNeuronSequences{4},2);
		
		vecSDLatencyMiss = nanstd(cellNeuronSequences{1},[],2);
		vecSDLatencyHits = nanstd(cellNeuronSequences{2},[],2);
		vecSDLatencyBase = nanstd(cellNeuronSequences{3},[],2);
		vecSDLatencyOverall = nanstd(cellNeuronSequences{4},[],2);
		
		%select neurons with sufficient data points
		intMin = 7; %7
		indSelect1 = (sum(~isnan(cellNeuronSequences{1}),2)>intMin);
		indSelect2 = (sum(~isnan(cellNeuronSequences{2}),2)>intMin);
		indSelect3 = (sum(~isnan(cellNeuronSequences{3}),2)>intMin);
		indSelect4 = (sum(~isnan(cellNeuronSequences{4}),2)>intMin);
		
		vecMeanLatencyMiss(~indSelect1) = [];
		vecMeanLatencyHits(~indSelect2) = [];
		vecMeanLatencyBase(~indSelect3) = [];
		vecSDLatencyMiss(~indSelect1) = [];
		vecSDLatencyHits(~indSelect2) = [];
		vecSDLatencyBase(~indSelect3) = [];
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
		
		matLatenciesMiss = cellNeuronSequences{1};
		vecLatenciesMiss = nanmean(matLatenciesMiss,2);
		matLatenciesHits = cellNeuronSequences{2};
		vecLatenciesHits = nanmean(matLatenciesHits,2);
		vecLatencies = mean(cat(2,vecLatenciesMiss,vecLatenciesHits),2);
		[dummy,vecReorder] = sort(vecLatencies);
		%{
		intStep = 0.5;
		vecBins = (-2+intStep/2):intStep:(2-intStep/2);
		figure
		subplot(4,2,1)
		vecNs = [19 38 67 104];
		for intN=1:4
			intNeuron=vecNs(intN);
			subplot(4,2,((intN-1)*2)+1)
			vecMissC = hist(matLatenciesMiss(intNeuron,:),vecBins);
			bar(vecBins,vecMissC/sum(vecMissC),'hist');
			ylim([0 0.6])
			
			subplot(4,2,((intN-1)*2)+2)
			vecHitC = hist(matLatenciesHits(intNeuron,:),vecBins);
			bar(vecBins,vecHitC/sum(vecHitC),'hist');
			ylim([0 0.6])
			
			title(sprintf('neuron %d; sd miss: %.3f; sd hit: %.3f',intNeuron,nanstd(matLatenciesMiss(intNeuron,:)),nanstd(matLatenciesHits(intNeuron,:))))
		end
		%}
		
		hSequenceStability = figure;
		subplot(2,2,1);
		
		vecLimMiss = 2*(nanmean(matLatenciesMiss(:))+nanstd(matLatenciesMiss(:)))*[-1 1];
		matColormap = colormap('redbluepurple');
		implot(matLatenciesMiss(vecReorder,:),vecLimMiss,matColormap,3);
		%colorbar;
		title('Neuronal response latency relative to pop. CoM [Miss]');
		ylabel('Neuron')
		xlabel('Population event during miss trial')
		
		subplot(2,2,2);
		
		
		vecLimHit = 2*(nanmean(matLatenciesHits(:))+nanstd(matLatenciesHits(:)))*[-1 1];
		matColormap = colormap('redbluepurple');
		implot(matLatenciesHits(vecReorder,:),vecLimHit,matColormap,3);
		%colorbar;
		title('Neuronal response latency relative to pop. CoM [Hits]');
		ylabel('Neuron')
		xlabel('Population event during hit trial')
		
		
		subplot(2,2,3);
		errorbar(1,nanmean(vecSDLatencyMiss),nanstd(vecSDLatencyMiss)/sqrt(sum(~isnan(vecSDLatencyMiss))),'xr');
		hold on;
		errorbar(2,nanmean(vecSDLatencyHits),nanstd(vecSDLatencyHits)/sqrt(sum(~isnan(vecSDLatencyHits))),'xg');
		errorbar(3,nanmean(vecSDLatencyBase),nanstd(vecSDLatencyBase)/sqrt(sum(~isnan(vecSDLatencyBase))),'xb');
		hold off;
		%ylim([0 max(get(gca,'ylim'))]);
		set(gca,'xtick',[1 2 3],'xticklabel',{'Miss','Hits','Base'});
		ylabel('Temporal sequence variability (s)');
		title('Mean consistency of activation order in population events across neurons');
		[h,p]=ttest2(vecSDLatencyHits,vecSDLatencyMiss);
		
		subplot(2,2,4);
		errorbar(1,nanmean(vecMeanLatencyMiss),nanstd(vecMeanLatencyMiss)/sqrt(sum(~isnan(vecMeanLatencyMiss))),'xr');
		hold on;
		errorbar(2,nanmean(vecMeanLatencyHits),nanstd(vecMeanLatencyHits)/sqrt(sum(~isnan(vecMeanLatencyHits))),'xg');
		errorbar(3,nanmean(vecMeanLatencyBase),nanstd(vecMeanLatencyBase)/sqrt(sum(~isnan(vecMeanLatencyBase))),'xb');
		hold off;
		%ylim([0 max(get(gca,'ylim'))]);
		set(gca,'xtick',[1 2 3],'xticklabel',{'Miss','Hits','Base'});
		ylabel('Mean latency (s)');
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
		close;
		
		%save data
		cellSaveTemporalSequenceStability{intPopulation} = [nanmean(vecSDLatencyMiss) nanmean(vecSDLatencyHits) nanmean(vecSDLatencyBase)]; %correlation temporal order within assemblies with temporal order overall
		
		%% get latency matrix [neuron x PE]
		intAssOccs = numel(vecAssemblyStarts);
		matLatencies = nan(intNeurons,intAssOccs);
		vecLatencySDs = nan(1,intAssOccs);
		vecAssDuringTrial = nan(1,intAssOccs);
		vecAssPreTrial = nan(1,intAssOccs);
		intPreFrames = round(2*dblSamplingFreq);
		for intOccurrence=1:intAssOccs
			%get spike times for this occurrence
			if ~indMultiNeuronClusters(vecAssemblies(intOccurrence)),continue;end
			intStart = vecAssemblyStarts(intOccurrence);
			intStop = vecAssemblyStops(intOccurrence)-1;
			matAEs = matSpikeCountsHD(:,intStart:intStop);
			
			%get mean and std
			matVals = repmat(1:size(matAEs,2),[size(matAEs,1) 1]);
			vecVals = getValByIdx(matVals,matAEs);
			dblMean = mean(vecVals);
			dblSD = std(vecVals);
			
			%get overall sd
			matLatencies(:,intOccurrence) = (sum((matAEs.*repmat(1:size(matAEs,2),[size(matAEs,1) 1])),2) ./ sum(matAEs,2) - dblMean)/dblSD;
			vecLatencySDs(intOccurrence) = dblSD;
			
			%get trial
			vecDuringTrials = find((sStimTemp.FrameOn > intStart & sStimTemp.FrameOn < intStop) |...
				(sStimTemp.FrameOff > intStart & sStimTemp.FrameOff < intStop));
			if ~isempty(vecDuringTrials)
				vecAssDuringTrial(intOccurrence) = vecDuringTrials;
			end
			
			%get pre-trial
			vecPreTrials = find((sStimTemp.FrameOn-intPreFrames > intStart & sStimTemp.FrameOn-intPreFrames < intStop) |...
				(sStimTemp.FrameOff-intPreFrames > intStart & sStimTemp.FrameOff-intPreFrames < intStop));
			if ~isempty(vecPreTrials)
				vecAssPreTrial(intOccurrence) = vecPreTrials;
			end
		end
		%matLatencies = (matLatencies / cellMultiSes{intPopulation}.samplingFreq)*1000;
		
		%get correlation with overall
		[h,vecSignificantOffset]=ttest(matLatencies');
		vecMeanOffset = nanmean(matLatencies,2);
		vecCorr = nan(1,intAssOccs);
		for intOccurrence=1:intAssOccs
			vecCorr(intOccurrence) = corr(vecMeanOffset,matLatencies(:,intOccurrence),'rows','pairwise','type','Pearson');
		end
		
		%% get variability per trial
		vecConsistencyPerTrial = nan(1,intTrials);
		vecPreConsistencyPerTrial = nan(1,intTrials);
		for intTrial=1:intTrials
			%get consistency
			vecOccurrences = find(vecAssDuringTrial==intTrial);
			if ~isempty(vecOccurrences)
				vecConsistencyPerTrial(intTrial) = nanmean(vecCorr(vecOccurrences));
			end
			
			%get pre-consistency
			vecPreOccurrences = find(vecAssPreTrial==intTrial);
			if ~isempty(vecPreOccurrences)
				vecPreConsistencyPerTrial(intTrial) = nanmean(vecCorr(vecPreOccurrences));
			end
		end
		%get consistency per response type
		vecHitConsistency = vecConsistencyPerTrial(indHits);
		vecMissConsistency = vecConsistencyPerTrial(indMiss);
		vecSlowConsistency = vecConsistencyPerTrial(indSlow);
		vecFastConsistency = vecConsistencyPerTrial(indFast);
		vecPreConsistency = vecPreConsistencyPerTrial;
		
		%save data resp corr
		cellSaveConsistencyPerRespType{intPopulation,1} = [nanmean(vecHitConsistency) nanmean(vecMissConsistency) nanmean(vecSlowConsistency) nanmean(vecFastConsistency) nanmean(vecPreConsistency)];
		
		%get consistencies per stimulus trial, split by highest/lowest 50%
		indStimTrials = ~cellSelectContrasts{1};
		vecOriStimTrials = vecStimOris(indStimTrials);
		intStimTrials = length(vecOriStimTrials);
		
		vecConsistencyPerStimTrial = vecConsistencyPerTrial(indStimTrials);
		vecLowConsistencyTrials = vecConsistencyPerStimTrial < nanmedian(vecConsistencyPerStimTrial);
		vecHighConsistencyTrials = vecConsistencyPerStimTrial >= nanmedian(vecConsistencyPerStimTrial);
		
		% decode orientation with template matching
		[dblPerformance_dFoF,vecDecodedIndexCV,matDists] = doCrossValidatedTemplateMatching(matTrialResponse(:,indStimTrials),vecOriStimTrials);
		[dblPerformance_AEs,vecDecodedIndexCV_AEs_TM,matDists] = doCrossValidatedTemplateMatching(matSpikeCounts(:,indStimTrials),vecOriStimTrials);
		[phat,dblCI_dFoF] = binofit(dblPerformance_dFoF*intStimTrials,intStimTrials);
		[phat,dblCI_AEs] = binofit(dblPerformance_AEs*intStimTrials,intStimTrials);
		
		%transform to correct
		vecDecodedCorrect = vecDecodedIndexCV_AEs_TM == vecOriStimTrials;
		
		%get decoding accuracy for lowest and highest 50% consistency trials
		vecDecAccLowConsist = vecDecodedCorrect(vecLowConsistencyTrials);
		vecDecAccHighConsist = vecDecodedCorrect(vecHighConsistencyTrials);
		
		%get response types during stimulus trials
		indHitStim = indHits(indStimTrials);
		indMissStim = indMiss(indStimTrials);
		
		%save data
		cellSaveConsistencyPerTrial{intPopulation,1} = vecDecAccLowConsist;
		cellSaveConsistencyPerTrial{intPopulation,2} = vecDecAccHighConsist;
		cellSaveConsistencyPerTrial{intPopulation,3} = vecDecodedCorrect(indMissStim);
		cellSaveConsistencyPerTrial{intPopulation,4} = vecDecodedCorrect(indHitStim);
		%}
		%%{
		%post-stimulus, 1 second time window
		%sStimTemp = getRandomMissDurations(sStimTemp,1);
		intWindowSize = round(1 * dblSamplingFreq);
		[vecDecodedIndexCV,vecRealIndex] = doCrossValidatedTemplateMatchingAcrossTime(matSpikeCountsHD,sStimTemp.FrameOn(indStimTrials),sStimTemp.FrameOff(indStimTrials),vecOriStimTrials,intWindowSize,false);
		
		%get centered on stim onset
		sEvents = [];
		sEvents.handleFig = [];
		sEvents.vecOn = sStimTemp.FrameOn(indStimTrials)+round(intWindowSize/2);
		sEvents.vecWindow = [-200 200];
		[cellHandles,sOut] = doPEP(sEvents,vecDecodedIndexCV==vecRealIndex);
		vecTimeOn = sOut.vecLineX;
		vecAccOn = sOut.vecLineY;
		
		%get centered on stim offset
		sEvents = [];
		sEvents.handleFig = [];
		sEvents.vecOn = sStimTemp.FrameOff(indStimTrials)+round(intWindowSize/2);
		sEvents.vecWindow = [-200 200];
		[cellHandles,sOut] = doPEP(sEvents,vecDecodedIndexCV==vecRealIndex);
		vecTimeOff = sOut.vecLineX;
		vecAccOff = sOut.vecLineY;
		
		%save data
		cellSaveDecodingAT{intPopulation,1} = vecTimeOn;
		cellSaveDecodingAT{intPopulation,2} = vecAccOn;
		cellSaveDecodingAT{intPopulation,3} = vecTimeOff;
		cellSaveDecodingAT{intPopulation,4} = vecAccOff;
		%}
		
		%% PE duration 1 frame: assembly occurrence probability drop-off (PSTH relative to occurrences, or to stimulus onset); self vs. other
		%get data
		%%{
		matAssemblyActivity = zeros(size(sAssemblies{intPopulation}.matAssemblyActivity)); %[assemblies x frames]
		for intAssembly=1:intAssemblies
			matAssemblyActivity(intAssembly,sAssemblies{intPopulation}.vecAssemblyStarts(sAssemblies{intPopulation}.vecAssemblies==intAssembly)) = 1;
		end
		dblSamplingFreq = 25.4;
		dblNumSecsPlot = 60;
		dblStimSecs = 2;
		
		%prep
		vecPSTH_Frames1 = -round(dblNumSecsPlot*dblSamplingFreq):round(dblNumSecsPlot*dblSamplingFreq);
		cellPSTH_self = cellfill(zeros(size(vecPSTH_Frames1)),[1 intAssemblies]);
		cellPSTH_others = cellfill(zeros(size(vecPSTH_Frames1)),[1 intAssemblies]);
		cellPSTH_Stim = cellfill(zeros(size(vecPSTH_Frames1)),[1 intAssemblies]);
		cellPSTH_Offset = cellfill(zeros(size(vecPSTH_Frames1)),[1 intAssemblies]);
		vecLaggingRatio = nan(1,intAssemblies);
		vecStimDrivenRatio = nan(1,intAssemblies);
		for intAssembly=find(indRealAssemblies)'
			vecOccurrences = matAssemblyActivity(intAssembly,:);
			indOtherAssemblies = true(1,intAssemblies);
			indOtherAssemblies(intAssembly) = false;
			vecOtherOccurrences = sum(matAssemblyActivity(indOtherAssemblies,:),1);
			for intOccFrame=find(logical(vecOccurrences))
				vecGetFrames = vecPSTH_Frames1+intOccFrame-1;
				vecAssgnFrames = 1:length(vecGetFrames);
				vecAssgnFrames(vecGetFrames<1 | vecGetFrames>length(vecOccurrences)) = [];
				vecGetFrames(vecGetFrames<1 | vecGetFrames>length(vecOccurrences)) = [];
				cellPSTH_self{intAssembly}(vecAssgnFrames) = cellPSTH_self{intAssembly}(vecAssgnFrames) + vecOccurrences(vecGetFrames);
				cellPSTH_others{intAssembly}(vecAssgnFrames) = cellPSTH_others{intAssembly}(vecAssgnFrames) + vecOtherOccurrences(vecGetFrames);
			end
			
			vecStimOn = cellMultiSes{intPopulation}.structStim.FrameOn;
			for intStimFrame=vecStimOn
				vecGetFrames = vecPSTH_Frames1+intStimFrame-1;
				vecAssgnFrames = 1:length(vecGetFrames);
				vecAssgnFrames(vecGetFrames<1 | vecGetFrames>length(vecOccurrences)) = [];
				vecGetFrames(vecGetFrames<1 | vecGetFrames>length(vecOccurrences)) = [];
				cellPSTH_Stim{intAssembly}(vecAssgnFrames) = cellPSTH_Stim{intAssembly}(vecAssgnFrames) + vecOccurrences(vecGetFrames);
			end
			
			vecStimOff = cellMultiSes{intPopulation}.structStim.FrameOff;
			for intStimFrame=vecStimOff(indHits)
				vecGetFrames = vecPSTH_Frames1+intStimFrame-1;
				vecAssgnFrames = 1:length(vecGetFrames);
				vecAssgnFrames(vecGetFrames<1 | vecGetFrames>length(vecOccurrences)) = [];
				vecGetFrames(vecGetFrames<1 | vecGetFrames>length(vecOccurrences)) = [];
				cellPSTH_Offset{intAssembly}(vecAssgnFrames) = cellPSTH_Offset{intAssembly}(vecAssgnFrames) + vecOccurrences(vecGetFrames);
			end
			
			%get data
			vecSelf = cellPSTH_self{intAssembly}/sum(vecOccurrences);
			vecOthers = cellPSTH_others{intAssembly}/sum(vecOccurrences);
			dblCorrFactorOthers = sum(vecOtherOccurrences)/length(vecOtherOccurrences);
			dblCorrFactorSelf = sum(vecOccurrences)/length(vecOccurrences);
			[vecPowerSelf,vecF] = pwelch(vecSelf/dblCorrFactorSelf,[],[],length(vecSelf),dblSamplingFreq,'onesided');
			[vecPowerOthers,vecFOthers] = pwelch(vecOthers/dblCorrFactorOthers,[],[],length(vecOthers),dblSamplingFreq,'onesided');
			
			intMid = ceil(length(vecOthers)/2);
			vecPre = (intMid-round(dblSamplingFreq)):(intMid-1);
			vecPost = (intMid+1):(intMid+round(dblSamplingFreq));
			vecData = vecOthers/dblCorrFactorOthers;
			vecPreData = vecData(vecPre);
			vecPostData = vecData(vecPost);
			
			dblLaggingRatio = mean(vecPreData)/mean(vecPostData);
			vecLaggingRatio(intAssembly) = dblLaggingRatio;
			
			vecSelfStim = cellPSTH_Stim{intAssembly};
			vecSelfOffset = cellPSTH_Offset{intAssembly};
			intStartStim = ceil(length(vecSelfStim)/2);
			intStimDur = round(dblSamplingFreq*dblStimSecs);
			vecStimDur = intStartStim:(intStartStim+intStimDur-1);
			dblBaseAct = mean(vecSelfStim(1:(intStartStim-1)));
			dblBaseActOffset = mean(vecSelfOffset(1:(intStartStim-1)));
			dblMeanStim = mean(vecSelfStim(vecStimDur));
			vecStimDrivenRatio(intAssembly) = dblMeanStim/dblBaseAct;
			
			%save figure
			if sParams.boolSavePlots
				figure
				subplot(2,2,1)
				plot([-dblNumSecsPlot dblNumSecsPlot],[1 1],'k--');
				hold on
				plot(vecPSTH_Frames1/dblSamplingFreq,vecOthers/dblCorrFactorOthers,'b')
				plot(vecPSTH_Frames1/dblSamplingFreq,vecSelf/dblCorrFactorSelf,'r')
				
				hold off;
				title(sprintf('Assembly %d; Blue=others, red=self',intAssembly));
				set(gca,'yscale','log');
				xlim([-dblNumSecsPlot dblNumSecsPlot])
				ylim([0.5 10]);
				xlabel('Time (s)')
				ylabel('Occ. prob. relative to random');
				
				subplot(2,2,2)
				pwelch(vecSelf/dblCorrFactorSelf,[],[],length(vecSelf),dblSamplingFreq,'onesided');
				title('Power spectrum self-occurrence (Welch)')
				
		%{
				subplot(2,2,3)
				errorbar(1:2,[mean(vecPreData) mean(vecPostData)],[std(vecPreData)/sqrt(dblSamplingFreq) std(vecPostData)/sqrt(dblSamplingFreq)],'x');
				title(sprintf('Lagging ratio, pre/post: %.3f',dblLaggingRatio));
				ylabel('Occurrence probability');
				set(gca,'xtick',[1 2],'xticklabel',{'Pre-occurrence','Post-occurrence'});
		%}
				subplot(2,2,3)
				hold on
				plot(vecPSTH_Frames1/dblSamplingFreq,vecSelfOffset/dblBaseAct,'r')
				hold off;
				title(sprintf('centered on stim offset'));
				%set(gca,'yscale','log');
				xlim([-8 8])
				xlabel('Time (s)')
				ylabel('Occ. prob. relative to random');
				
				
				
				subplot(2,2,4)
				hold on
				plot(vecPSTH_Frames1/dblSamplingFreq,vecSelfStim/dblBaseAct,'r')
				hold off;
				title(sprintf('centered on stim onset'));
				%set(gca,'yscale','log');
				xlim([-5 10])
				xlabel('Time (s)')
				ylabel('Occ. prob. relative to random');
				
				drawnow;
				jFig = get(handle(gcf), 'JavaFrame');
				jFig.setMaximized(true);
				figure(gcf);
				drawnow;
				export_fig(sprintf('%spop%d_ass%doccurrence_probability.tif',strSes,intPopulation,intAssembly));
				export_fig(sprintf('%spop%d_ass%doccurrence_probability.pdf',strSes,intPopulation,intAssembly));
			end
			
			%assign data
			cellSaveAssemblyRecurrenceMinDur{intPopulation,intAssembly,1} = vecF;
			cellSaveAssemblyRecurrenceMinDur{intPopulation,intAssembly,2} = vecPowerSelf;
			cellSaveAssemblyRecurrenceMinDur{intPopulation,intAssembly,3} = vecPowerOthers;
			
			%occurrence probability centered on event
			cellSaveAssemblyRecurrenceMinDur{intPopulation,intAssembly,4} = vecPSTH_Frames1/dblSamplingFreq;
			cellSaveAssemblyRecurrenceMinDur{intPopulation,intAssembly,5} = vecSelf/dblCorrFactorSelf;
			cellSaveAssemblyRecurrenceMinDur{intPopulation,intAssembly,6} = vecOthers/dblCorrFactorOthers;
			
			%occurrence probability centered on stim presentation
			cellSaveAssemblyRecurrenceMinDur{intPopulation,intAssembly,7} = vecSelfStim;
			cellSaveAssemblyRecurrenceMinDur{intPopulation,intAssembly,8} = dblBaseAct;
			cellSaveAssemblyRecurrenceMinDur{intPopulation,intAssembly,9} = vecSelfOffset;
			cellSaveAssemblyRecurrenceMinDur{intPopulation,intAssembly,10} = dblBaseActOffset;
		end
		%}
		
		
		%% HMM analysis
		%{
		%set parameters
		intStates = 5;
		matAssemblyActivity = sAssemblies{intPopulation}.matAssemblyActivity; %[assemblies x frames]
		intT = size(matAssemblyActivity,2);
		dblSamplingFreq = 25.4;
		vecSequence = sum(matAssemblyActivity.*repmat((1:intAssemblies)',[1 intT]),1)+1;
		
		%calculate transition probability matrix
		indDiag = tril(true(intStates)) & triu(true(intStates));
		matT = zeros(intStates,intStates);
		matT(indDiag) = 0.99 + 0.009*rand(1,intStates);
		vecDiagT = matT(indDiag);
		matT = repmat((1-matT(indDiag))/(intStates-1),[1 intStates]);
		matT(indDiag) = vecDiagT; %[from-state x to-state]
		
		%calculate emission probability matrix
		matE = repmat(sum(matAssemblyActivity,2)'/intT,[intStates 1]); %[states x symbols]
		matE = [1-sum(matE,2) matE]; %add 0-symbol
		
		
		%estimate
		vecStates = zeros(1,intT);
		for intTrial=find(indStimTrials)
			vecStates(sStimTemp.FrameOn(intTrial):sStimTemp.FrameOff(intTrial)) = vecStimOris(intTrial);
		end
		vecStates = vecStates+1;
		
		%get transition + emission estimates
		[matT,matE] = hmmestimate(vecSequence,vecStates);
	
		%perform fitting
		[ESTTR,ESTEMIT] = hmmtrain(vecSequence,matT,matE,'tolerance',0.0001,'maxiterations',100,'verbose',true);
		
		%decode
		matDecodedStateProb = hmmdecode(vecSequence,ESTTR,ESTEMIT);
		[vecProb,vecDecState]=max(matDecodedStateProb,[],1);
		
		%{
		%plot
		sEvents = [];
		sEvents.vecOn = sStimTemp.FrameOn;
		vecTrace = vecStates==vecDecState;
		vecTrace(vecStates==1) = 0;
		doPEP(sEvents,vecTrace);
		%}
		
		%trial based
		vecDecodedOris = zeros(1,intTrials);
		cellDecode = cell(2,6);
		for intTrial=1:intTrials
			vecFrames = sStimTemp.FrameOn(intTrial):sStimTemp.FrameOff(intTrial);
			vecProbs = sum(matDecodedStateProb(:,vecFrames),2);
			[dummy,intOri]=max(vecProbs(2:5));
			vecDecodedOris(intTrial) = intOri;
			
			%get contrast / resp
			intR = indHits(intTrial)+1;
			intC = vecStimContrasts(intTrial);
			cellDecode{intR,intC} = [cellDecode{intR,intC} intOri==vecStimOris(intTrial)];
		end
		cellSaveOriRespDecoding{intPopulation} = cellDecode;
		cellSaveOverallOriRespDecoding{intPopulation} = cellfun(@mean,cellDecode);
		
		%miss
		dblAlpha = normcdf(1)-normcdf(-1);
		[phatM,pciM] = binofit(cellfun(@sum,cellDecode(1,:)),cellfun(@numel,cellDecode(1,:)),dblAlpha);
		vecCorrFacM = sqrt(cellfun(@numel,cellDecode(1,:)));
		
		%hit
		[phatH,pciH] = binofit(cellfun(@sum,cellDecode(2,:)),cellfun(@numel,cellDecode(2,:)),dblAlpha);
		vecCorrFacH = sqrt(cellfun(@numel,cellDecode(2,:)));
		if sParams.boolSavePlots
			%plot
			figure
			vecContrasts = [0 0.5 2 8 32 100];
			vecPlotC = vecContrasts;
			vecPlotC(1) = 0.2;
			errorfill(vecPlotC,phatM,(pciM(:,2)'-phatM)./vecCorrFacM,(phatM-pciM(:,1)')./vecCorrFacM,[0.7 0 0],[1 0.5 0.5])
			hold on
			errorfill(vecPlotC,phatH,(pciH(:,2)'-phatH)./vecCorrFacH,(phatH-pciH(:,1)')./vecCorrFacH,[0 0.7 0],[0.5 1 0.5])
			hold off
			set(gca,'XScale','log')
			set(gca,'Xtick',vecPlotC,'xticklabel',vecContrasts)
			ylabel('Mean dF/F0')
			xlabel('Stimulus contrast (%)')
			
			%save plot
			drawnow;
			export_fig(sprintf('%spop%d_HMM_oriresp_decoding.tif',strSes,intPopulation));
			export_fig(sprintf('%spop%d_HMM_oriresp_decoding.pdf',strSes,intPopulation));
		end
		%}
		
		
		%% check presence ===================> ALSO FA/CR!!!
		%{
		matAssemblyActivity = sAssemblies{intPopulation}.matAssemblyActivity; %[assemblies x frames]
		intT = size(matAssemblyActivity,2);
		dblSamplingFreq = 25.4;
		vecSequence = sum(matAssemblyActivity.*repmat((1:intAssemblies)',[1 intT]),1);
		vecSequence(vecSequence==0)=nan;
		
		%check sequences for different stimuli
		%trial based
		matAssemblyPresenceTrials = nan(intAssemblies,intTrials);
		for intTrial=1:intTrials
			matAssemblyPresenceTrials(:,intTrial) = ismember(1:intAssembly,vecSequence(sStimTemp.FrameOn(intTrial):sStimTemp.FrameOff(intTrial)));
		end
		
		% check if identity is correlated with stimulus detection
		intHits = sum(indHits);
		intMiss = sum(indMiss);
		vecHitPresRatio = sum(matAssemblyPresenceTrials(:,indHits),2)/intHits;
		vecMissPresRatio = sum(matAssemblyPresenceTrials(:,indMiss),2)/intMiss;
		
		%perform for shuffled
		intIters = 1000;
		matHitPresRatioShuf = nan(intAssemblies,intIters);
		matMissPresRatioShuf = nan(intAssemblies,intIters);
		for intIter=1:intIters
			%shuffle assembly identity
			matAssemblyActivityS = zeros(size(sAssemblies{intPopulation}.matAssemblyActivity)); %[assemblies x frames]
			intOccs = numel(sAssemblies{intPopulation}.vecAssemblyStarts);
			vecAssembliesS = sAssemblies{intPopulation}.vecAssemblies(randperm(intOccs));
			for intAssemblyOcc=1:intOccs
				matAssemblyActivityS(vecAssembliesS(intAssemblyOcc),sAssemblies{intPopulation}.vecAssemblyStarts(intAssemblyOcc):sAssemblies{intPopulation}.vecAssemblyStops(intAssemblyOcc)) = 1;
			end
			vecSequenceS = sum(matAssemblyActivityS.*repmat((1:intAssemblies)',[1 intT]),1);
			vecSequenceS(vecSequenceS==0)=nan;
			
			%check sequences for different stimuli
			%trial based
			matAssemblyPresenceTrialsS = nan(intAssemblies,intTrials);
			for intTrial=1:intTrials
				matAssemblyPresenceTrialsS(:,intTrial) = ismember(1:intAssembly,vecSequenceS(sStimTemp.FrameOn(intTrial):sStimTemp.FrameOff(intTrial)));
			end
			
			% check if identity is correlated with stimulus detection
			matHitPresRatioShuf(:,intIter) = sum(matAssemblyPresenceTrialsS(:,indHits),2)/intHits;
			matMissPresRatioShuf (:,intIter)= sum(matAssemblyPresenceTrialsS(:,indMiss),2)/intMiss;
		end
		
		%transform to z-scores
		vecHitPresRatioZ = (vecHitPresRatio-mean(matHitPresRatioShuf,2))./std(matHitPresRatioShuf,[],2);
		vecMissPresRatioZ = (vecMissPresRatio-mean(matMissPresRatioShuf,2))./std(matMissPresRatioShuf,[],2);
		
		%save data
		cellSavePresRatioZ{intPopulation,1} = vecHitPresRatioZ;
		cellSavePresRatioZ{intPopulation,2} = vecMissPresRatioZ;
		%}
		%% analyze sequences
		matAssemblyActivity = sAssemblies{intPopulation}.matAssemblyActivity; %[assemblies x frames]
		intT = size(matAssemblyActivity,2);
		dblSamplingFreq = 25.4;
		vecSequence = sum(matAssemblyActivity.*repmat((1:intAssemblies)',[1 intT]),1);
		vecSequence(vecSequence==0)=nan;
		
		%check sequences for different stimuli
		%trial based
		cellSequence = cell(1,intTrials);
		for intTrial=1:intTrials
			cellSequence{intTrial} = vecSequence(sStimTemp.FrameOn(intTrial):sStimTemp.FrameOff(intTrial));
		end

		%get sequence similarity
		matSeqSimOrder = nan(intTrials,intTrials);
		matSeqSimOverlap = nan(intTrials,intTrials);
		for intTrial1=1:intTrials
			vecSequence1 = cellSequence{intTrial1};
			vecShortSequence1 = vecSequence1([true (diff(vecSequence1)~=0)] & ~isnan(vecSequence1));
			for intTrial2=(intTrial1+1):intTrials
				vecSequence2 = cellSequence{intTrial2};
				
				%short sequence: order
				vecShortSequence2 = vecSequence2([true (diff(vecSequence2)~=0)] & ~isnan(vecSequence2));
				intMaxShort = min([numel(vecShortSequence1) numel(vecShortSequence2)]);
				if intMaxShort > 0
					matSeqSimOrder(intTrial1,intTrial2) = sum(vecShortSequence1(1:intMaxShort) == vecShortSequence2(1:intMaxShort))/intMaxShort;
					matSeqSimOrder(intTrial2,intTrial1) = matSeqSimOrder(intTrial1,intTrial2);
				end
				
				%long sequence: overlap
				intMaxLong = min([numel(vecSequence1) numel(vecSequence2)]);
				if intMaxLong > 0
					matSeqSimOverlap(intTrial1,intTrial2) = sum(vecSequence1(1:intMaxLong) == vecSequence2(1:intMaxLong))/sum(~isnan(vecSequence1(1:intMaxLong)) & ~isnan(vecSequence1(1:intMaxLong)));
					matSeqSimOverlap(intTrial2,intTrial1) = matSeqSimOverlap(intTrial1,intTrial2);
				end
			end
		end
		
		%get sorting index
		[dummy,vecReorder] = sort(vecStimOris);
		matOriSimOrder = getBlockMeans(matSeqSimOrder,vecStimOris);
		matOriSimOverlap = getBlockMeans(matSeqSimOverlap,vecStimOris);
		
		%hit/miss/CR/FA
		matRespSimOrder = getBlockMeans(matSeqSimOrder,vecRespType);
		matRespSimOverlap = getBlockMeans(matSeqSimOverlap,vecRespType);
		
		
		%save data
		cellSaveEnsembleSequences{intPopulation,1} = matOriSimOrder;
		cellSaveEnsembleSequences{intPopulation,2} = matOriSimOverlap;
		cellSaveEnsembleSequences{intPopulation,3} = matRespSimOrder;
		cellSaveEnsembleSequences{intPopulation,4} = matRespSimOverlap;
		
		
		
		%% repeat analysis; shuffle ensemble identities
		matAssemblyActivityS = zeros(size(sAssemblies{intPopulation}.matAssemblyActivity)); %[assemblies x frames]
		%shuffle assembly identity
		intOccs = numel(sAssemblies{intPopulation}.vecAssemblyStarts);
		vecAssembliesS = sAssemblies{intPopulation}.vecAssemblies(randperm(intOccs));
		for intAssemblyOcc=1:intOccs
			matAssemblyActivityS(vecAssembliesS(intAssemblyOcc),sAssemblies{intPopulation}.vecAssemblyStarts(intAssemblyOcc):sAssemblies{intPopulation}.vecAssemblyStops(intAssemblyOcc)) = 1;
		end
		vecSequenceS = sum(matAssemblyActivityS.*repmat((1:intAssemblies)',[1 intT]),1);
		vecSequenceS(vecSequenceS==0)=nan;
		
		%check sequences for different stimuli
		%trial based
		cellSequenceS = cell(1,intTrials);
		for intTrial=1:intTrials
			cellSequenceS{intTrial} = vecSequenceS(sStimTemp.FrameOn(intTrial):sStimTemp.FrameOff(intTrial));
		end
		
		%get sequence similarity
		matSeqSimOrderS = nan(intTrials,intTrials);
		matSeqSimOverlapS = nan(intTrials,intTrials);
		for intTrial1=1:intTrials
			vecSequence1 = cellSequenceS{intTrial1};
			vecShortSequence1 = vecSequence1([true (diff(vecSequence1)~=0)] & ~isnan(vecSequence1));
			for intTrial2=(intTrial1+1):intTrials
				vecSequence2 = cellSequenceS{intTrial2};
				
				%short sequence: order
				vecShortSequence2 = vecSequence2([true (diff(vecSequence2)~=0)] & ~isnan(vecSequence2));
				intMaxShort = min([numel(vecShortSequence1) numel(vecShortSequence2)]);
				if intMaxShort > 0
					matSeqSimOrderS(intTrial1,intTrial2) = sum(vecShortSequence1(1:intMaxShort) == vecShortSequence2(1:intMaxShort))/intMaxShort;
					matSeqSimOrderS(intTrial2,intTrial1) = matSeqSimOrderS(intTrial1,intTrial2);
				end
				
				%long sequence: overlap
				intMaxLong = min([numel(vecSequence1) numel(vecSequence2)]);
				if intMaxLong > 0
					matSeqSimOverlapS(intTrial1,intTrial2) = sum(vecSequence1(1:intMaxLong) == vecSequence2(1:intMaxLong))/sum(~isnan(vecSequence1(1:intMaxLong)) & ~isnan(vecSequence1(1:intMaxLong)));
					matSeqSimOverlapS(intTrial2,intTrial1) = matSeqSimOverlapS(intTrial1,intTrial2);
				end
			end
		end
		
		%get sorting index
		[dummy,vecReorder] = sort(vecStimOris);
		matOriSimOrderS = getBlockMeans(matSeqSimOrderS,vecStimOris);
		matOriSimOverlapS = getBlockMeans(matSeqSimOverlapS,vecStimOris);
		
		%hit/miss
		matRespSimOrderS = getBlockMeans(matSeqSimOrderS,vecRespType);
		matRespSimOverlapS = getBlockMeans(matSeqSimOverlapS,vecRespType);
		
		
		%save data
		cellSaveEnsembleSequences{intPopulation,5} = matOriSimOrderS;
		cellSaveEnsembleSequences{intPopulation,6} = matOriSimOverlapS;
		cellSaveEnsembleSequences{intPopulation,7} = matRespSimOrderS;
		cellSaveEnsembleSequences{intPopulation,8} = matRespSimOverlapS;
		
		%indCR = indMiss & vecStimContrasts == 1;
		%indFA = indHits & vecStimContrasts == 1;
		
		%plot response corr
		vecSeqMissesShuffled = matSeqSimOrderS(indMissStim,indMissStim);
		vecSeqMissesShuffled = vecSeqMissesShuffled(~isnan(vecSeqMissesShuffled));
		vecSeqHitsShuffled = matSeqSimOrderS(indHitsStim,indHitsStim);
		vecSeqHitsShuffled = vecSeqHitsShuffled(~isnan(vecSeqHitsShuffled));
		vecSeqCRShuffled = matSeqSimOrderS(indCR,indCR);
		vecSeqCRShuffled = vecSeqCRShuffled(~isnan(vecSeqCRShuffled));
		vecSeqFAShuffled = matSeqSimOrderS(indFA,indFA);
		vecSeqFAShuffled = vecSeqFAShuffled(~isnan(vecSeqFAShuffled));
		
		vecSeqMisses = matSeqSimOrder(indMissStim,indMissStim);
		vecSeqMisses = vecSeqMisses(~isnan(vecSeqMisses));
		vecSeqHits = matSeqSimOrder(indHitsStim,indHitsStim);
		vecSeqHits = vecSeqHits(~isnan(vecSeqHits));
		vecSeqCRs = matSeqSimOrder(indCR,indCR);
		vecSeqCRs = vecSeqCRs(~isnan(vecSeqCRs));
		vecSeqFAs = matSeqSimOrder(indFA,indFA);
		vecSeqFAs = vecSeqFAs(~isnan(vecSeqFAs));
		
		vecMeans = [mean(vecSeqMissesShuffled) mean(vecSeqMisses) mean(vecSeqHitsShuffled) mean(vecSeqHits)];
		vecSEMs = [std(vecSeqMissesShuffled)/sqrt(numel(vecSeqMissesShuffled)) std(vecSeqMisses)/sqrt(numel(vecSeqMisses)) std(vecSeqHitsShuffled)/sqrt(numel(vecSeqHitsShuffled)) std(vecSeqHits)/sqrt(numel(vecSeqHits))];
		
		vecMeans2 = [mean(vecSeqCRShuffled) mean(vecSeqCRs) mean(vecSeqFAShuffled) mean(vecSeqFAs)];
		vecSEMs2 = [std(vecSeqCRShuffled)/sqrt(numel(vecSeqCRShuffled)) std(vecSeqCRs)/sqrt(numel(vecSeqCRs)) std(vecSeqFAShuffled)/sqrt(numel(vecSeqFAShuffled)) std(vecSeqFAs)/sqrt(numel(vecSeqFAs))];
		
		%t-tests
		[h,pMissShuffled]=ttest2(vecSeqMissesShuffled,vecSeqMisses);
		[h,pHitShuffled]=ttest2(vecSeqHitsShuffled,vecSeqHits);
		[h,pMissHit]=ttest2(vecSeqHits,vecSeqMisses);
		
		figure
		errorbar(1:8,[vecMeans vecMeans2],[vecSEMs vecSEMs2],'xb');
		set(gca,'xtick',1:8,'xticklabel',{'Miss shuf','Miss','Hit Shuf','Hit','CR shuf','CR','FA Shuf','FA'});
		title(sprintf('T-tests, miss-shuf,p=%.3f; hit-shuf,p=%.3f;hit-miss,p=%.3f',pMissShuffled,pHitShuffled,pMissHit))
		ylim([0.1 0.3+eps])
		ylabel('Sequence similarity')
		drawnow;
		export_fig(sprintf('%spop%d_ensemble_sequence_corr_behav.tif',strSes,intPopulation));
		export_fig(sprintf('%spop%d_ensemble_sequence_corr_behav.pdf',strSes,intPopulation));
		
		%{
		%% plot examples
		%hits: [2 100]
		%misses: [21 22]
		close all
		figure
		hold on
		cellSeq = [];
		intS=0;
		cMap=nanjet(intAssemblies+1);
		for intSeq = [2 100 21 22]
			intS=intS+1;
			vecSeq = cellSequence{intSeq};
			vecSeq([false diff(vecSeq)~=0] & ~isnan([nan vecSeq(1:(end-1))])) = nan;
			cellSeq{intS} = vecSeq;
			
			vecSeq(isnan(vecSeq)) = 12;
			h = cline((1:length(vecSeq))/dblSamplingFreq,intS*ones(size(vecSeq)),[],vecSeq);
			
			set(h,'EdgeColor','flat','linewidth',50)
		end
		h = cline((1:11)/dblSamplingFreq,-1*ones(1,11),[],1:11);
		ylim([0 5])
		xlim([0 3])
		colormap(nanjet(intAssemblies))
		colorbar
		xlabel('Time after stimulus onset (s)')
		ylabel('Trials: [2 100 21 22]')
		title('Hits: 2 100, misses: 21 22')
		
		%get sequence similarity
		matSeqSimOrderExamples = nan(4,4);
		matSeqSimOverlapExamples = nan(4,4);
		for intTrial1=1:4
			vecSequence1 = cellSeq{intTrial1};
			vecShortSequence1 = vecSequence1([true (diff(vecSequence1)~=0)] & ~isnan(vecSequence1));
			for intTrial2=(intTrial1+1):4
				vecSequence2 = cellSeq{intTrial2};
				
				%short sequence: order
				vecShortSequence2 = vecSequence2([true (diff(vecSequence2)~=0)] & ~isnan(vecSequence2));
				intMaxShort = min([numel(vecShortSequence1) numel(vecShortSequence2)]);
				if intMaxShort > 0
					matSeqSimOrderExamples(intTrial1,intTrial2) = sum(vecShortSequence1(1:intMaxShort) == vecShortSequence2(1:intMaxShort))/intMaxShort;
					matSeqSimOrderExamples(intTrial2,intTrial1) = matSeqSimOrderExamples(intTrial1,intTrial2);
				end
				
				%long sequence: overlap
				intMaxLong = min([numel(vecSequence1) numel(vecSequence2)]);
				if intMaxLong > 0
					matSeqSimOverlapExamples(intTrial1,intTrial2) = sum(vecSequence1(1:intMaxLong) == vecSequence2(1:intMaxLong))/sum(~isnan(vecSequence1(1:intMaxLong)) & ~isnan(vecSequence1(1:intMaxLong)));
					matSeqSimOverlapExamples(intTrial2,intTrial1) = matSeqSimOverlapExamples(intTrial1,intTrial2);
				end
			end
		end
		%}
		
		%% plot
		if sParams.boolSavePlots
			%plot
			figure
			subplot(2,2,1)
			imagesc(matSeqSimOrder(1:50,1:50),[0 1]);colormap(hot(256))
			xlabel('Trial #')
			ylabel('Trial #')
			title('Cluster sequence similarity [0-1]')
			
			subplot(2,2,2)
			imagesc(matSeqSimOverlap(1:50,1:50),[0 1]);colormap(hot(256))
			xlabel('Trial #')
			ylabel('Trial #')
			title('Cluster presence across time similarity [0-1]')
			
			
			subplot(2,2,3)
			imagesc(matOriSimOrder,[0 1]);colormap(hot(256))
			xlabel('Stimulus orientation')
			ylabel('Stimulus orientation')
			title('Mean cluster sequence similarity [0-0.2]')
			
			subplot(2,2,4)
			imagesc(matOriSimOverlap,[0 1]);colormap(hot(256))
			xlabel('Stimulus orientation')
			ylabel('Stimulus orientation')
			title('Mean cluster presence across time similarity [0-0.1]')
			
			drawnow;
			jFig = get(handle(gcf), 'JavaFrame');
			jFig.setMaximized(true);
			figure(gcf);
			drawnow;
			export_fig(sprintf('%spop%d_ensemble_sequences.tif',strSes,intPopulation));
			export_fig(sprintf('%spop%d_ensemble_sequences.pdf',strSes,intPopulation));
		end
	end
	
	%% save data structures
	if sParams.boolSaveData
		vecClock = fix(clock);
		strDate = num2str(vecClock(1:3));
		strFile = ['dataAssAnalFive' strSes '_' strrep(strDate,'    ','_')];
		save(['D:\Data\Results\spikeAnalysis\' strFile],'cellSave*','-v7.3');
		cd(strOldDir);
	end
end

