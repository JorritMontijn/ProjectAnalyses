%load data
for intMouse=5
	close all
	drawnow;
	clearvars -except intMouse
	
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
	sParams.strFigDir = ['D:\Data\Results\spikeAnalysis' filesep strSes filesep];
	sParams.boolSavePlots = false;
	sParams.boolSaveData = false;
	%strOldDir = cd(sParams.strFigDir);
	
	%change name for no split
	%#ok<*AGROW>
	strSes = ['AA' strSes];
	
	for intPopulation = 1:numel(cellMultiSes)
		%% GET RESPONSES AFTER COMBINING BLOCKS OF SAME POPULATION INTO ONE
		fprintf('Starting %s pop %d [%s]\n',strSes,intPopulation,getTime);
		
		%define stimulus duration parameter
		dblStimSecs = 2;
		
		% remove trials with reaction time <150ms
		dblRemSecs = 0.15;
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
		
		%% population CCG
		%pop-act based CCGs
		intFrames = size(matSpikeCountsHD,2);
		matSpikeCountsHD(:,end+1) = 0; %#ok<SAGROW>
		matCCG = zeros(intNeurons,length(vecPlotCCG));
		for intNeuron=1:intNeurons
			vecAEs = matSpikeCountsHD(intNeuron,:);
			vecSpikeTimes = [];
			while any(vecAEs>0)
				vecTheseTimes = find(vecAEs>0);
				vecSpikeTimes = [vecSpikeTimes vecTheseTimes];
				vecAEs(vecTheseTimes) = vecAEs(vecTheseTimes) - 1;
			end
			
			for intSpike=1:length(vecSpikeTimes)
				vecFrames = vecPlotCCG + vecSpikeTimes(intSpike);
				vecFrames(vecFrames<1) = intFrames+1;
				vecFrames(vecFrames>intFrames) = intFrames+1;
				matThisCCG = matSpikeCountsHD(:,vecFrames);
				matThisCCG(intNeuron,:)=0;
				matCCG = matCCG + matThisCCG;
			end
		end
		matSpikeCountsHD(:,end) = [];
		%matCCG(end+1,:) = sum(matCCG,1);
		
		%get mucc
		vecMuCC = (sum((matCCG.*repmat(1:size(matCCG,2),[size(matCCG,1) 1])),2) ./ sum(matCCG,2) - ceil(length(vecPlotCCG)/2))/dblSamplingFreq;
		[dummy,vecReorder] = sort(vecMuCC,'ascend');
		matCCGnorm = bsxfun(@rdivide,matCCG,sum(matCCG,2));
		vecFilt = normpdf(-2:2,0,1);
		vecFilt = vecFilt/sum(vecFilt);
		matCCGnormBlur = conv2(matCCGnorm,vecFilt,'same');
		
		%save figure
		figure
		imagesc(vecPlotCCG/dblSamplingFreq,[],matCCGnorm(vecReorder,:));colormap('hot');
		if sParams.boolSavePlots
			drawnow;
			export_fig(sprintf('%spop%dpopulation_muCC.tif',strSes,intPopulation));
			export_fig(sprintf('%spop%dpopulation_muCC.pdf',strSes,intPopulation));
		end
		
		%% collect ordering of all neurons in each population event during miss/slow/fast trials, then check if ordering is more consistent during fast than slow trial assembly occurrences
		cellNeuronSequencesAssemblies = cell(intAssemblies,3);
		cellNeuronSequences = cellfill(nan(intNeurons,length(vecAssemblies)),[1 4]);
		intMissCounterOverall = 0;
		intSlowCounterOverall = 0;
		intFastCounterOverall = 0;
		for intAssembly=1:intAssemblies
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
				intStop = vecStops(intOccurrence);
				matAEs = matSpikeCountsHD(:,intStart:intStop);
				
				%check if during trial & get relative timing to center of mass of this occurrence
				if any(ismember(intStart:intStop,vecMissHD))
					intMissCounterAssembly = intMissCounterAssembly + 1;
					intMissCounterOverall = intMissCounterOverall + 1;
					
					vecThisMuCC = (sum((matAEs.*repmat(1:size(matAEs,2),[size(matAEs,1) 1])),2) ./ sum(matAEs,2) - centerOfMass(matAEs,2))/dblSamplingFreq;
					cellNeuronSequencesAssemblies{intAssembly,1}(:,intMissCounterAssembly) = vecThisMuCC;
					cellNeuronSequences{1}(:,intMissCounterOverall) = vecThisMuCC;
				end
				if any(ismember(intStart:intStop,vecSlowHD))
					intSlowCounterAssembly = intSlowCounterAssembly + 1;
					intSlowCounterOverall = intSlowCounterOverall + 1;
					
					vecThisMuCC = (sum((matAEs.*repmat(1:size(matAEs,2),[size(matAEs,1) 1])),2) ./ sum(matAEs,2) - centerOfMass(matAEs,2))/dblSamplingFreq;
					cellNeuronSequencesAssemblies{intAssembly,2}(:,intSlowCounterAssembly) = vecThisMuCC;
					cellNeuronSequences{2}(:,intSlowCounterOverall) = vecThisMuCC;
				end
				if any(ismember(intStart:intStop,vecFastHD))
					intFastCounterAssembly = intFastCounterAssembly + 1;
					intFastCounterOverall = intFastCounterOverall + 1;
					
					vecThisMuCC = (sum((matAEs.*repmat(1:size(matAEs,2),[size(matAEs,1) 1])),2) ./ sum(matAEs,2) - centerOfMass(matAEs,2))/dblSamplingFreq;
					cellNeuronSequencesAssemblies{intAssembly,3}(:,intFastCounterAssembly) = vecThisMuCC;
					cellNeuronSequences{3}(:,intFastCounterOverall) = vecThisMuCC;
				end
				
				%get overall sd
				vecThisMuCC = (sum((matAEs.*repmat(1:size(matAEs,2),[size(matAEs,1) 1])),2) ./ sum(matAEs,2) - centerOfMass(matAEs,2))/dblSamplingFreq;
				cellNeuronSequences{4}(:,intOccurrence) = vecThisMuCC;
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
		cellNeuronSequences{4}(:,(intFastCounterOverall+1):end) = [];
		
		vecMeanLatencyMiss = nanmean(cellNeuronSequences{1},2);
		vecMeanLatencySlow = nanmean(cellNeuronSequences{2},2);
		vecMeanLatencyFast = nanmean(cellNeuronSequences{3},2);
		vecMeanLatencyHit = nanmean(cat(2,cellNeuronSequences{2},cellNeuronSequences{3}),2);
		
		vecSDLatencyMiss = nanstd(cellNeuronSequences{1},[],2);
		vecSDLatencySlow = nanstd(cellNeuronSequences{2},[],2);
		vecSDLatencyFast = nanstd(cellNeuronSequences{3},[],2);
		vecSDLatencyHits = nanstd(cat(2,cellNeuronSequences{2},cellNeuronSequences{3}),[],2);
		vecSDLatencyOverall = nanstd(cellNeuronSequences{4},[],2);
		
		%plot
		hSequenceStability = figure;
		
		subplot(2,4,1);
		matLatenciesSlow = cellNeuronSequences{2};
		matLatenciesSlowRelative = bsxfun(@minus,matLatenciesSlow,vecMeanLatencySlow);
		vecLimSlow = 2*(nanmean(matLatenciesSlow(:))+nanstd(matLatenciesSlow(:)))*[-1 1];
		matColormap = colormap('redbluepurple');
		implot(matLatenciesSlow,vecLimSlow,matColormap,3);
		%colorbar;
		title('Neuronal response latency relative to pop. CoM [Slow]');
		ylabel('Neuron')
		xlabel('Population event during slow trial')
		
		subplot(2,4,2);
		matLatenciesFast = cellNeuronSequences{3};
		matLatenciesFastRelative = bsxfun(@minus,matLatenciesFast,vecMeanLatencyFast);
		vecLimFast = 2*(nanmean(matLatenciesFast(:))+nanstd(matLatenciesFast(:)))*[-1 1];
		implot(matLatenciesFast,vecLimFast,matColormap,3);
		%colorbar;
		title('Neuronal response latency relative to pop. CoM [Fast]');
		ylabel('Neuron')
		xlabel('Population event during fast trial')
		
		subplot(2,4,3);
		matLatenciesMiss = cellNeuronSequences{1};
		matLatenciesMissRelative = bsxfun(@minus,matLatenciesMiss,vecMeanLatencyMiss);
		vecLimSlow = 2*(nanmean(matLatenciesMiss(:))+nanstd(matLatenciesMiss(:)))*[-1 1];
		matColormap = colormap('redbluepurple');
		implot(matLatenciesMiss,vecLimSlow,matColormap,3);
		%colorbar;
		title('Neuronal response latency relative to pop. CoM [Miss]');
		ylabel('Neuron')
		xlabel('Population event during miss trial')
		
		subplot(2,4,4);
		matLatenciesHit = cat(2,cellNeuronSequences{2},cellNeuronSequences{3});
		matLatenciesHitRelative = bsxfun(@minus,matLatenciesHit,vecMeanLatencyHit);
		vecLimFast = 2*(nanmean(matLatenciesHit(:))+nanstd(matLatenciesHit(:)))*[-1 1];
		implot(matLatenciesHit,vecLimFast,matColormap,3);
		%colorbar;
		title('Neuronal response latency relative to pop. CoM [Hit]');
		ylabel('Neuron')
		xlabel('Population event during hit trial')
		
		subplot(2,2,3);
		errorbar(1,nanmean(vecSDLatencySlow),nanstd(vecSDLatencySlow)/sqrt(sum(~isnan(vecSDLatencySlow))),'xm');
		hold on;
		errorbar(2,nanmean(vecSDLatencyFast),nanstd(vecSDLatencyFast)/sqrt(sum(~isnan(vecSDLatencyFast))),'xg');
		hold off;
		ylim([0 max(get(gca,'ylim'))]);
		set(gca,'xtick',[1 2],'xticklabel',{'Slow','Fast'});
		ylabel('Temporal sequence variability (s)');
		title('Mean consistency of activation order in population events across neurons');
		
		subplot(2,2,4);
		errorbar(1,nanmean(vecSDLatencyMiss),nanstd(vecSDLatencyMiss)/sqrt(sum(~isnan(vecSDLatencyMiss))),'xm');
		hold on;
		errorbar(2,nanmean(vecSDLatencyHits),nanstd(vecSDLatencyHits)/sqrt(sum(~isnan(vecSDLatencyHits))),'xg');
		hold off;
		ylim([0 max(get(gca,'ylim'))]);
		set(gca,'xtick',[1 2],'xticklabel',{'Miss','Hit'});
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
		cellSaveTemporalSequenceStability{intPopulation} = [nanmean(vecSDLatencyMiss) nanmean(vecSDLatencyHits) nanmean(vecSDLatencySlow) nanmean(vecSDLatencyFast) nanmean(vecSDLatencyOverall)]; %correlation temporal order within assemblies with temporal order overall
		
		%% get neuronal CCGs in assemblies
		%{
		intFrames = size(matSpikeCountsHD,2);
		matSpikeCountsHD(:,end+1) = 0; %#ok<SAGROW>
		for intAssembly=find(indRealAssemblies)'
			indOccurrences = vecAssemblies == intAssembly;
			vecStarts = vecAssemblyStarts(indOccurrences);
			vecStops = vecAssemblyStops(indOccurrences);
			intCoreMembers = sum(matAssemblyCoreMembers(intAssembly,:));
			intC1 = 0;
			figure
			title(sprintf('Cluster %d',intAssembly));
			for intNeuron1=find(matAssemblyCoreMembers(intAssembly,:))%1:(intNeurons-1)
				indMember1 = matAssemblyMembers(intNeuron1,:) & indOccurrences;
				intC1 = intC1 + 1;
				intC2 = 0;
				for intNeuron2=find(matAssemblyCoreMembers(intAssembly,:))%(intNeuron1+1):intNeurons
					indMember2 = matAssemblyMembers(intNeuron2,:) & indOccurrences;
					intC2 = intC2 + 1;
					
					vecCCG = zeros(size(vecPlotCCG));
					intSizeCCG = length(vecCCG);
					for intOccurrence=find(indMember1 & indMember2)
						intStart = vecAssemblyStarts(intOccurrence);
						intStop = vecAssemblyStops(intOccurrence);
						vecFrames = intStart:intStop;
						vecFrames(vecFrames<1) = intFrames+1;
						vecFrames(vecFrames>intFrames) = intFrames+1;
						intLength = intStop-intStart+1;
						vecAct1 = matSpikeCountsHD(intNeuron1,vecFrames);
						vecAct2 = matSpikeCountsHD(intNeuron2,vecFrames);
						
						%put in ccg
						while any(vecAct1>0)
							intLoc = find(vecAct1>0,1,'first');
							intShiftPosStart = intMiddle - intLoc + 1;
							vecLoc = intShiftPosStart:(intShiftPosStart+(intLength-1));
							indLoc = ~(vecLoc>intSizeCCG | vecLoc < 1);
							vecCCG(vecLoc(indLoc)) = vecCCG(vecLoc(indLoc)) + vecAct2(indLoc);
							vecAct1(intLoc) = vecAct1(intLoc) - 1;
						end
					end
					subplot(intCoreMembers,intCoreMembers,(intCoreMembers*(intC1-1)+intC2))
					plot(vecPlotCCG,vecCCG)
					drawnow;
				end
			end
			pause
		end
		matSpikeCountsHD(:,end) = [];
		%}
		
		%% population cross correlograms of neurons in assemblies => relative to only spikes from significant members
		%pre-allocate outputs
		matRealCorrs = nan(3,intAssemblies); %mean, upper, lower
		matMuCC_Assembly = nan(intNeurons,intAssemblies);
		intIters = 2;
		matShuffledCorrs = nan(intIters,intAssemblies); %mean, upper, lower
		
		%perform analysis for real data
		for intAssembly=1:intAssemblies
			indOccurrences = vecAssemblies == intAssembly;
			vecStarts = vecAssemblyStarts(indOccurrences);
			vecStops = vecAssemblyStops(indOccurrences);
			vecCoreMembers = find(matAssemblyCoreMembers(intAssembly,:));
			intCoreMembers = sum(matAssemblyCoreMembers(intAssembly,:));
			intFrames = size(matSpikeCountsHD,2);
			
			%get ccg
			matSpikeCountsHD(:,end+1) = 0; %#ok<SAGROW>
			matCCG_assembly = zeros(intNeurons,length(vecPlotCCG));
			for intNeuron=1:intNeurons%vecCoreMembers
				%get spike times for this neuron
				vecAEs = matSpikeCountsHD(intNeuron,:);
				vecSpikeTimes = [];
				while any(vecAEs>0)
					vecTheseTimes = find(vecAEs>0);
					vecSpikeTimes = [vecSpikeTimes vecTheseTimes];
					vecAEs(vecTheseTimes) = vecAEs(vecTheseTimes) - 1;
				end
				
				%get assembly occurrences with this neuron
				vecAssemblySpikes = [];
				indMember = matAssemblyMembers(intNeuron,:) & indOccurrences;
				for intOccurrence=1:length(vecStarts)
					intStart = vecAssemblyStarts(intOccurrence);
					intStop = vecAssemblyStops(intOccurrence);
					vecTheseTimes = vecSpikeTimes(vecSpikeTimes >= intStart & vecSpikeTimes <= intStop);
					vecAssemblySpikes = [vecAssemblySpikes vecTheseTimes];
				end
				
				for intSpike=1:length(vecAssemblySpikes)
					vecFrames = vecPlotCCG + vecAssemblySpikes(intSpike);
					vecFrames(vecFrames<1) = intFrames+1;
					vecFrames(vecFrames>intFrames) = intFrames+1;
					matThisCCG = matSpikeCountsHD(:,vecFrames);
					matThisCCG(intNeuron,:)=0;
					matCCG_assembly = matCCG_assembly + matThisCCG;
				end
			end
			matSpikeCountsHD(:,end) = [];
			%matCCG(end+1,:) = sum(matCCG,1);
			
			%calculate MuCC's and correlations
			vecMuCC = (sum((matCCG.*repmat(1:size(matCCG,2),[size(matCCG,1) 1])),2) ./ sum(matCCG,2) - ceil(length(vecPlotCCG)/2))/dblSamplingFreq;
			vecMuCC_assembly = (sum((matCCG_assembly.*repmat(1:size(matCCG_assembly,2),[size(matCCG_assembly,1) 1])),2) ./ sum(matCCG_assembly,2) - ceil(length(vecPlotCCG)/2))/dblSamplingFreq;
			indRemNeurons = isnan(vecMuCC_assembly); %remove neurons that by chance have no spikes
			[R,P,RLO,RUP]=corrcoef(vecMuCC(~indRemNeurons),vecMuCC_assembly(~indRemNeurons));
			matRealCorrs(:,intAssembly) = [R(1,2); RUP(1,2); RLO(1,2)];
			matMuCC_Assembly(:,intAssembly) = vecMuCC_assembly;
			
			%plot assembly CCG
			%{
			indRemNeurons = isnan(vecMuCC_assembly);
			vecMuCC(indRemNeurons) = [];
			vecMuCC_assembly(indRemNeurons) = [];
			matCCG_assembly(indRemNeurons,:) = [];
			
			dblR = R(2,1);
			dblLowerR = RLO(2,1);
			dblHigherR = RUP(2,1);
			[dummy,vecReorder] = sort(vecMuCC_assembly,'ascend');
			matCCGnorm = bsxfun(@rdivide,matCCG_assembly,sum(matCCG_assembly,2));
			vecFilt = normpdf(-2:2,0,1);
			vecFilt = vecFilt/sum(vecFilt);
			matCCGnormBlur = conv2(matCCGnorm,vecFilt,'same');
			imagesc(vecPlotCCG/dblSamplingFreq,[],matCCGnorm(vecReorder,:));colormap('hot');
			title(sprintf('Cluster %d, temporal order correlation, r=%.3f [95%% CI: %.3f - %.3f]',intAssembly,dblR,dblLowerR,dblHigherR));
			pause
			%}
		end
		
		
		%% make picture
		%{
		[vecMeanSort,vecSort] = sort(vecMuCC);
		matMeanSortedAssemblies = matMuCC_Assembly(vecSort,:);
		
		subplot(2,6,1)
		scatter(vecMeanSort,1:numel(vecMeanSort),[],vecMeanSort,'filled')
		caxis([-0.4 0.4])
		xlim([-0.7 0.7])
		colormap('redbluepurple')
		axis ij
		
		for intAssembly=1:intAssemblies
			subplot(2,6,intAssembly+1)
			vecCoreMembers = find(matAssemblyCoreMembers(intAssembly,:));
			vecThisAssembly = matMeanSortedAssemblies(:,intAssembly);
			
			scatter(vecThisAssembly,1:numel(vecMeanSort),[],vecThisAssembly,'filled');
			%hold on
			%scatter(matMeanSortedAssemblies(vecCoreMembers,intAssembly),vecCoreMembers,[],matMeanSortedAssemblies(vecCoreMembers,intAssembly),'filled');
			%hold off
			dblR = corrcoef(vecMeanSort(~isnan(vecThisAssembly)),vecThisAssembly(~isnan(vecThisAssembly)));
			dblR = dblR(1,2);
			caxis([-0.4 0.4])
			xlim([-0.7 0.7])
			colormap('redbluepurple')
			axis ij
			title(sprintf('corr: %.3f',dblR))
		end
		
		%save figure
		if sParams.boolSavePlots
			drawnow;
			jFig = get(handle(gcf), 'JavaFrame');
			jFig.setMaximized(true);
			figure(gcf);
			drawnow;
			export_fig(sprintf('%spop%dtemporal_sequence_example.tif',strSes,intPopulation));
			export_fig(sprintf('%spop%dtemporal_sequence_example.pdf',strSes,intPopulation));
		end
		%}
		
		
		
		%% compare with shuffled, is corr higher than shuffled?
		%get all spike times and counts per neuron
		vecSpikesPerNeuron = sum(matSpikeCountsHD>0,2);
		intTotalSpikeNr = sum(vecSpikesPerNeuron);
		vecAllSpikeTimes = nan(1,intTotalSpikeNr);
		vecAllSpikeNumbers = nan(1,intTotalSpikeNr);
		vecLoc = cumsum(vecSpikesPerNeuron)';
		vecLoc = [0 vecLoc];
		for intNeuron=1:intNeurons
			vecSpikeLocations = find(matSpikeCountsHD(intNeuron,:)>0);
			vecIndex = (vecLoc(intNeuron)+1):vecLoc(intNeuron+1);
			vecAllSpikeTimes(vecIndex) = vecSpikeLocations;
			vecAllSpikeNumbers(vecIndex) = matSpikeCountsHD(intNeuron,vecSpikeLocations);
		end
		
		%do shuffling
		for intIter=1:3%intIters
			%randomly distribute spike times across neurons (keeping counts intact)
			vecShuffleIndex = randperm(intTotalSpikeNr);
			vecShuffledSpikeTimes = vecAllSpikeTimes(vecShuffleIndex);
			vecShuffledSpikeNumber = vecAllSpikeNumbers(vecShuffleIndex);
			matShuffledSpikeCountsHD = zeros(size(matSpikeCountsHD));
			for intNeuron=1:intNeurons
				vecIndex = (vecLoc(intNeuron)+1):vecLoc(intNeuron+1);
				vecTheseShuffledTimes = vecShuffledSpikeTimes(vecIndex);
				vecTheseShuffledNumber = vecShuffledSpikeNumber(vecIndex);
				matShuffledSpikeCountsHD(intNeuron,vecTheseShuffledTimes) = vecTheseShuffledNumber;
			end
			
			%calculate Assembly MuCC
			for intAssembly=1:intAssemblies
				indOccurrences = vecAssemblies == intAssembly;
				vecStarts = vecAssemblyStarts(indOccurrences);
				vecStops = vecAssemblyStops(indOccurrences);
				vecCoreMembers = find(matAssemblyCoreMembers(intAssembly,:));
				intCoreMembers = sum(matAssemblyCoreMembers(intAssembly,:));
				intFrames = size(matShuffledSpikeCountsHD,2);
				
				%get ccg
				matShuffledSpikeCountsHD(:,end+1) = 0; %#ok<SAGROW>
				matCCG_assembly = zeros(intNeurons,length(vecPlotCCG));
				for intNeuron=1:intNeurons%vecCoreMembers
					%get spike times for this neuron
					vecAEs = matShuffledSpikeCountsHD(intNeuron,:);
					vecSpikeTimes = [];
					while any(vecAEs>0)
						vecTheseTimes = find(vecAEs>0);
						vecSpikeTimes = [vecSpikeTimes vecTheseTimes];
						vecAEs(vecTheseTimes) = vecAEs(vecTheseTimes) - 1;
					end
					
					%get assembly occurrences with this neuron
					vecAssemblySpikes = [];
					indMember = matAssemblyMembers(intNeuron,:) & indOccurrences;
					for intOccurrence=1:length(vecStarts)
						intStart = vecAssemblyStarts(intOccurrence);
						intStop = vecAssemblyStops(intOccurrence);
						vecTheseTimes = vecSpikeTimes(vecSpikeTimes >= intStart & vecSpikeTimes <= intStop);
						vecAssemblySpikes = [vecAssemblySpikes vecTheseTimes];
					end
					
					for intSpike=1:length(vecAssemblySpikes)
						vecFrames = vecPlotCCG + vecAssemblySpikes(intSpike);
						vecFrames(vecFrames<1) = intFrames+1;
						vecFrames(vecFrames>intFrames) = intFrames+1;
						matThisCCG = matShuffledSpikeCountsHD(:,vecFrames);
						matThisCCG(intNeuron,:)=0;
						matCCG_assembly = matCCG_assembly + matThisCCG;
					end
				end
				matShuffledSpikeCountsHD(:,end) = [];
				%matCCG(end+1,:) = sum(matCCG,1);
				
				%calculate MuCC's and correlations
				vecMuCC = (sum((matCCG.*repmat(1:size(matCCG,2),[size(matCCG,1) 1])),2) ./ sum(matCCG,2) - ceil(length(vecPlotCCG)/2))/dblSamplingFreq;
				vecMuCC_assembly = (sum((matCCG_assembly.*repmat(1:size(matCCG_assembly,2),[size(matCCG_assembly,1) 1])),2) ./ sum(matCCG_assembly,2) - ceil(length(vecPlotCCG)/2))/dblSamplingFreq;
				indRemNeurons = isnan(vecMuCC_assembly); %remove neurons that by chance have no spikes
				if sum(~indRemNeurons)>0
					vecMuCC(indRemNeurons) = [];
					vecMuCC_assembly(indRemNeurons) = [];
					dblCorr=corr(vecMuCC,vecMuCC_assembly);
					matShuffledCorrs(intIter,intAssembly) = dblCorr;
				end
				%fprintf('Correlation for asssembly %d/%d [iter %d/%d], normal MuCC-Shuffled MuCC: %.3f [%s]\n',intAssembly,intAssemblies,intIter,intIters,dblCorr,getTime);
			end
			fprintf('Correlation for assemblies computed; iter %d/%d [%s]\n',intIter,intIters,getTime);
		end
		
		%% plot
		%calculate 95% CI of shuffled distro's
		vecPlotRealAssemblies = find(indRealAssemblies)';
		vecPlotAssemblies = find(true(size(indRealAssemblies)))';
		matSC = matShuffledCorrs(:,vecPlotAssemblies);
		vecShuffledMean = nanmean(matShuffledCorrs,1);
		vecShuffledSD = nanstd(matShuffledCorrs,[],1);
		
		%get z-scores
		vecRealMean = matRealCorrs(1,vecPlotAssemblies);
		vecZ = (vecRealMean - vecShuffledMean) ./ vecShuffledSD;
		vecP = 1-erf(vecZ/sqrt(2));
		vecP_corr = bonf_holm(vecP);
		
		
		%plot correlations
		hTemporalOrdering = figure;
		subplot(2,2,1)
		scatter(1:length(vecPlotAssemblies),vecRealMean,'bx','linewidth',2);
		hold on
		errorbar(1:length(vecPlotAssemblies),vecShuffledMean,2*vecShuffledSD,2*vecShuffledSD,'rx');
		hold off
		title(['p-values: ' sprintf('%.6f;',vecP)]);
		ylabel('Pearson correlation of temporal order of neuronal activation (muCC)');
		set(gca,'xtick', 1:length(vecPlotAssemblies),'xticklabel',vecPlotAssemblies);
		xlabel('Cluster number')
		ylim([-1-eps 1+eps]*max(get(gca,'ylim')));
		legend('Real data','Shuffled data','Location','best')
		
		%redefine real assemblies ot those with also significant ordering
		%indSig = false(size(indRealAssemblies));
		%indSig(indRealAssemblies) = vecP'<0.05;
		%indRealAssemblies = indSig;
		
		%% get miss/slow/fast muCC correlations
		%pre-allocate outputs
		matBehavCorrs = nan(3,intAssemblies); %miss, slow, fast
		matBehavMuCC_Assembly = nan(intNeurons,intAssemblies,3); %miss, slow, fast
		
		%perform analysis for real data
		for intAssembly=1:intAssemblies
			%msg
			fprintf('Calculating assembly %d consistency for miss/slow/fast trials [%s]\n',intAssembly,getTime);
			
			%get occurrences
			indOccurrences = vecAssemblies == intAssembly;
			vecStarts = vecAssemblyStarts(indOccurrences);
			vecStops = vecAssemblyStops(indOccurrences);
			
			vecCoreMembers = find(matAssemblyCoreMembers(intAssembly,:));
			intCoreMembers = sum(matAssemblyCoreMembers(intAssembly,:));
			intFrames = size(matSpikeCountsHD,2);
			
			%get ccg
			matSpikeCountsHD(:,end+1) = 0; %#ok<SAGROW>
			matCCG_assemblyMiss = zeros(intNeurons,length(vecPlotCCG));
			matCCG_assemblySlow = zeros(intNeurons,length(vecPlotCCG));
			matCCG_assemblyFast = zeros(intNeurons,length(vecPlotCCG));
			for intNeuron=vecCoreMembers
				%get spike times for this neuron
				vecAEs = matSpikeCountsHD(intNeuron,:);
				vecSpikeTimes = [];
				while any(vecAEs>0)
					vecTheseTimes = find(vecAEs>0);
					vecSpikeTimes = [vecSpikeTimes vecTheseTimes];
					vecAEs(vecTheseTimes) = vecAEs(vecTheseTimes) - 1;
				end
				
				%get assembly occurrences with this neuron
				vecAssemblySpikesMiss = [];
				vecAssemblySpikesSlow = [];
				vecAssemblySpikesFast = [];
				indMember = matAssemblyMembers(intNeuron,:) & indOccurrences;
				for intOccurrence=1:length(vecStarts)
					intStart = vecAssemblyStarts(intOccurrence);
					intStop = vecAssemblyStops(intOccurrence);
					vecTheseTimes = vecSpikeTimes(vecSpikeTimes >= intStart & vecSpikeTimes <= intStop);
					if any(ismember(intStart:intStop,vecMissHD))
						vecAssemblySpikesMiss = [vecAssemblySpikesMiss vecTheseTimes];
					end
					if any(ismember(intStart:intStop,vecSlowHD))
						vecAssemblySpikesSlow = [vecAssemblySpikesSlow vecTheseTimes];
					end
					if any(ismember(intStart:intStop,vecFastHD))
						vecAssemblySpikesFast = [vecAssemblySpikesFast vecTheseTimes];
					end
				end
				
				%miss CCG
				for intSpike=1:length(vecAssemblySpikesMiss)
					vecFrames = vecPlotCCG + vecAssemblySpikesMiss(intSpike);
					vecFrames(vecFrames<1) = intFrames+1;
					vecFrames(vecFrames>intFrames) = intFrames+1;
					matThisCCG = matSpikeCountsHD(:,vecFrames);
					matThisCCG(intNeuron,:)=0;
					matCCG_assemblyMiss = matCCG_assemblyMiss + matThisCCG;
				end
				
				%slow CCG
				for intSpike=1:length(vecAssemblySpikesSlow)
					vecFrames = vecPlotCCG + vecAssemblySpikesSlow(intSpike);
					vecFrames(vecFrames<1) = intFrames+1;
					vecFrames(vecFrames>intFrames) = intFrames+1;
					matThisCCG = matSpikeCountsHD(:,vecFrames);
					matThisCCG(intNeuron,:)=0;
					matCCG_assemblySlow = matCCG_assemblySlow + matThisCCG;
				end
				
				%fast CCG
				for intSpike=1:length(vecAssemblySpikesFast)
					vecFrames = vecPlotCCG + vecAssemblySpikesFast(intSpike);
					vecFrames(vecFrames<1) = intFrames+1;
					vecFrames(vecFrames>intFrames) = intFrames+1;
					matThisCCG = matSpikeCountsHD(:,vecFrames);
					matThisCCG(intNeuron,:)=0;
					matCCG_assemblyFast = matCCG_assemblyFast + matThisCCG;
				end
			end
			matSpikeCountsHD(:,end) = [];
			%matCCG(end+1,:) = sum(matCCG,1);
			
			%calculate MuCC's and correlations for miss
			vecMuCC_assemblyMiss = (sum((matCCG_assemblyMiss.*repmat(1:size(matCCG_assemblyMiss,2),[size(matCCG_assemblyMiss,1) 1])),2) ./ sum(matCCG_assemblyMiss,2) - ceil(length(vecPlotCCG)/2))/dblSamplingFreq;
			indRemNeurons = isnan(vecMuCC_assemblyMiss); %remove neurons that by chance have no spikes
			if sum(~indRemNeurons) == 0,matBehavCorrs(1,intAssembly)=nan;
			else matBehavCorrs(1,intAssembly) = corr(matMuCC_Assembly(~indRemNeurons,intAssembly),vecMuCC_assemblyMiss(~indRemNeurons));end %miss, slow, fast
			matBehavMuCC_Assembly(:,intAssembly,1) = vecMuCC_assemblyMiss; %miss, slow, fast
			
			%calculate MuCC's and correlations for slow
			vecMuCC_assemblySlow = (sum((matCCG_assemblySlow.*repmat(1:size(matCCG_assemblySlow,2),[size(matCCG_assemblySlow,1) 1])),2) ./ sum(matCCG_assemblySlow,2) - ceil(length(vecPlotCCG)/2))/dblSamplingFreq;
			indRemNeurons = isnan(vecMuCC_assemblySlow); %remove neurons that by chance have no spikes
			if sum(~indRemNeurons) == 0,matBehavCorrs(1,intAssembly)=nan;
			else matBehavCorrs(2,intAssembly) = corr(matMuCC_Assembly(~indRemNeurons,intAssembly),vecMuCC_assemblySlow(~indRemNeurons));end %miss, slow, fast
			matBehavMuCC_Assembly(:,intAssembly,2) = vecMuCC_assemblySlow; %miss, slow, fast
			
			%calculate MuCC's and correlations for fast
			vecMuCC_assemblyFast = (sum((matCCG_assemblyFast.*repmat(1:size(matCCG_assemblyFast,2),[size(matCCG_assemblyFast,1) 1])),2) ./ sum(matCCG_assemblyFast,2) - ceil(length(vecPlotCCG)/2))/dblSamplingFreq;
			indRemNeurons = isnan(vecMuCC_assemblyFast); %remove neurons that by chance have no spikes
			if sum(~indRemNeurons) == 0,matBehavCorrs(1,intAssembly)=nan;
			else matBehavCorrs(3,intAssembly) = corr(matMuCC_Assembly(~indRemNeurons,intAssembly),vecMuCC_assemblyFast(~indRemNeurons));end %miss, slow, fast
			matBehavMuCC_Assembly(:,intAssembly,3) = vecMuCC_assemblyFast; %miss, slow, fast
			
			%plot assembly CCG
			%{
			indRemNeurons = isnan(vecMuCC_assembly);
			vecMuCC(indRemNeurons) = [];
			vecMuCC_assembly(indRemNeurons) = [];
			matCCG_assembly(indRemNeurons,:) = [];
			
			dblR = R(2,1);
			dblLowerR = RLO(2,1);
			dblHigherR = RUP(2,1);
			[dummy,vecReorder] = sort(vecMuCC_assembly,'ascend');
			matCCGnorm = bsxfun(@rdivide,matCCG_assembly,sum(matCCG_assembly,2));
			vecFilt = normpdf(-2:2,0,1);
			vecFilt = vecFilt/sum(vecFilt);
			matCCGnormBlur = conv2(matCCGnorm,vecFilt,'same');
			imagesc(vecPlotCCG/dblSamplingFreq,[],matCCGnorm(vecReorder,:));colormap('hot');
			title(sprintf('Cluster %d, temporal order correlation, r=%.3f [95%% CI: %.3f - %.3f]',intAssembly,dblR,dblLowerR,dblHigherR));
			pause
			%}
		end
		
		%% plot miss/slow/fast muCC correlations
		vecMissAssemblyConsistencies = matBehavCorrs(1,:);
		vecSlowAssemblyConsistencies = matBehavCorrs(2,:);
		vecFastAssemblyConsistencies = matBehavCorrs(3,:);
		
		vecMuCC = (sum((matCCG.*repmat(1:size(matCCG,2),[size(matCCG,1) 1])),2) ./ sum(matCCG,2) - ceil(length(vecPlotCCG)/2))/dblSamplingFreq;
		[vecMuCC_All,vecReorder] = sort(vecMuCC,'ascend');
		matAssemblyMuCC_All = matMuCC_Assembly(vecReorder,:);
		matAssemblyMuCCMiss = matBehavMuCC_Assembly(vecReorder,:,1);
		matAssemblyMuCCSlow = matBehavMuCC_Assembly(vecReorder,:,2);
		matAssemblyMuCCFast = matBehavMuCC_Assembly(vecReorder,:,3);
		
		dblCorrOverallMiss = corr(matAssemblyMuCC_All(:),matAssemblyMuCCMiss(:),'rows','pairwise');
		
		figure(hTemporalOrdering)
		
		subplot(2,2,2);
		[h,dblP]=ttest(vecSlowAssemblyConsistencies,vecFastAssemblyConsistencies);
		errorbar(1,nanmean(vecSlowAssemblyConsistencies),nanstd(vecSlowAssemblyConsistencies)/sqrt(intAssemblies),'mx');
		hold on
		errorbar(2,nanmean(vecFastAssemblyConsistencies),nanstd(vecFastAssemblyConsistencies)/sqrt(intAssemblies),'gx');
		hold off
		ylim([-1-eps 1+eps]*max(get(gca,'ylim')));
		ylabel('Mean temporal consistency of assemblies (r)');
		xlabel('Behavioral response speed');
		set(gca,'xtick',[1 2],'xticklabel',{'Slowest 50%','Fastest 50%'});
		title(sprintf('Paired t-test slow-miss, p=%.3f, n=%d assemblies',dblP,intAssemblies));
		
		
		subplot(2,2,3);
		hold on
		for intAssembly=1:size(matAssemblyMuCC_All,2);
			scatter(matAssemblyMuCC_All(:,intAssembly),matAssemblyMuCCSlow(:,intAssembly),'mx');
		end
		hold off
		ylim([-1 1]);
		xlim([-1 1]);
		dblCorrOverallSlow = corr(matAssemblyMuCC_All(:),matAssemblyMuCCSlow(:),'rows','pairwise');
		title(sprintf('Correlation of temporal ordering in assemblies for slow trials with overall, r=%.3f',dblCorrOverallSlow));
		ylabel('Temporal order (MuCC) slow trials (s)')
		xlabel('Temporal order (MuCC) whole data set (s)')
		
		
		subplot(2,2,4);
		hold on
		for intAssembly=1:size(matAssemblyMuCC_All,2);
			scatter(matAssemblyMuCC_All(:,intAssembly),matAssemblyMuCCFast(:,intAssembly),'gx');
		end
		hold off
		ylim([-1 1]);
		xlim([-1 1]);
		dblCorrOverallFast = corr(matAssemblyMuCC_All(:),matAssemblyMuCCFast(:),'rows','pairwise');
		title(sprintf('Correlation of temporal ordering in assemblies for fast trials with overall, r=%.3f',dblCorrOverallFast));
		ylabel('Temporal order (MuCC) slow trials (s)')
		xlabel('Temporal order (MuCC) whole data set (s)')
		
		%save figure
		if sParams.boolSavePlots
			drawnow;
			jFig = get(handle(gcf), 'JavaFrame');
			jFig.setMaximized(true);
			figure(gcf);
			drawnow;
			export_fig(sprintf('%spop%dtemporal_sequence_assemblies.tif',strSes,intPopulation));
			export_fig(sprintf('%spop%dtemporal_sequence_assemblies.pdf',strSes,intPopulation));
		end
		
		
		%% save data
		cellSaveTemporalOrder{intPopulation,1} = [mean(vecRealMean(indRealAssemblies)) mean(vecRealMean(~indRealAssemblies))]; %correlation temporal order within assemblies with temporal order overall
		cellSaveTemporalOrder{intPopulation,2} = vecSlowAssemblyConsistencies; %correlation temporal order of assembly sequence during slow trials with overall assembly sequence
		cellSaveTemporalOrder{intPopulation,3} = vecFastAssemblyConsistencies; %correlation temporal order of assembly sequence during fast trials with overall assembly sequence
		cellSaveTemporalOrder{intPopulation,4} = [dblCorrOverallSlow dblCorrOverallFast]; %overall correlations slow/fast
		
		%% assembly occurrence probability drop-off (PSTH relative to occurrences, or to stimulus onset); self vs. other
		close all
		dblNumSecsPlot = 60;
		vecPSTH_Frames1 = -round(dblNumSecsPlot*dblSamplingFreq):round(dblNumSecsPlot*dblSamplingFreq);
		cellPSTH_self = cellfill(zeros(size(vecPSTH_Frames1)),[1 intAssemblies]);
		cellPSTH_others = cellfill(zeros(size(vecPSTH_Frames1)),[1 intAssemblies]);
		cellPSTH_Stim = cellfill(zeros(size(vecPSTH_Frames1)),[1 intAssemblies]);
		vecLaggingRatio = nan(1,intAssemblies);
		vecStimDrivenRatio = nan(1,intAssemblies);
		for intAssembly=1:intAssemblies
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
			intStartStim = ceil(length(vecSelfStim)/2);
			intStimDur = round(dblSamplingFreq*dblStimSecs);
			vecStimDur = intStartStim:(intStartStim+intStimDur-1);
			dblBaseAct = mean(vecSelfStim(1:(intStartStim-1)));
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
				xlabel('Time (s)')
				ylabel('Occ. prob. relative to random');
				
				subplot(2,2,2)
				pwelch(vecSelf/dblCorrFactorSelf,[],[],length(vecSelf),dblSamplingFreq,'onesided');
				title('Power spectrum self-occurrence (Welch)')
				
				subplot(2,2,3)
				errorbar(1:2,[mean(vecPreData) mean(vecPostData)],[std(vecPreData)/sqrt(dblSamplingFreq) std(vecPostData)/sqrt(dblSamplingFreq)],'x');
				title(sprintf('Lagging ratio, pre/post: %.3f',dblLaggingRatio));
				ylabel('Occurrence probability');
				set(gca,'xtick',[1 2],'xticklabel',{'Pre-occurrence','Post-occurrence'});
				
				
				subplot(2,2,4)
				plot([-10 10],dblBaseAct*[1 1],'k--');
				hold on
				plot(vecPSTH_Frames1/dblSamplingFreq,vecSelfStim,'r')
				hold off;
				title(sprintf('centered on stim onset'));
				%set(gca,'yscale','log');
				xlim([-10 10])
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
			cellSaveAssemblyRecurrence{intPopulation,intAssembly,1} = vecF;
			cellSaveAssemblyRecurrence{intPopulation,intAssembly,2} = vecPowerSelf;
			cellSaveAssemblyRecurrence{intPopulation,intAssembly,3} = vecPowerOthers;
			
			%occurrence probability centered on event
			cellSaveAssemblyRecurrence{intPopulation,intAssembly,4} = vecPSTH_Frames1/dblSamplingFreq;
			cellSaveAssemblyRecurrence{intPopulation,intAssembly,5} = vecSelf/dblCorrFactorSelf;
			cellSaveAssemblyRecurrence{intPopulation,intAssembly,6} = vecOthers/dblCorrFactorOthers;
			
			%occurrence probability centered on stim presentation
			cellSaveAssemblyRecurrence{intPopulation,intAssembly,7} = vecSelfStim;
			cellSaveAssemblyRecurrence{intPopulation,intAssembly,8} = dblBaseAct;
			
		end
		
		%% compare OSI within assemblies vs across assemblies; are assemblies more orientation tuned than random?
		%get assembly OSIs
		indStimC = ~cellSelectC{1};%cellSelectC{6} | cellSelectC{5};
		matRespAssStim = matRespAssemblies(:,indStimC);
		cellSelectStimO = cell(size(cellSelectO));
		for intOri=1:length(cellSelectStimO)
			cellSelectStimO{intOri} = cellSelectO{intOri}(indStimC);
		end
		sTuning = calcTuningRespMat(matRespAssStim,cellSelectStimO,vecOrientations);
		vecAssOSI = sTuning.vecOSI;
		
		%get OSIs of neurons
		matSpikeCountsStim = matSpikeCounts(:,indStimC);
		sTuning = calcTuningRespMat(matSpikeCountsStim,cellSelectStimO,vecOrientations);
		vecNeuronPA = 2*ang2rad(sTuning.vecPrefAngle);
		mat_dPA = abs(circ_dist2(vecNeuronPA,vecNeuronPA));
		
		%compare mean dPA of core assembly members to random group of same size & get OSI
		intIters = 1000;
		matAssemblyMean_dPA_Shuffled = nan(intIters,intAssemblies);
		vecAssemblyMean_dPA = nan(1,intAssemblies);
		matAssembly_OSI_Shuffled = nan(intIters,intAssemblies);
		vecAssembly_OSI = nan(1,intAssemblies);
		for intAssembly=1:intAssemblies
			%get members
			vecMembers = matAssemblyCoreMembers(intAssembly,:);
			intMembers = sum(vecMembers);
			if intMembers < 2,continue;end
			matSelect = tril(true(intMembers),-1);
			
			%get real dPA
			matReal_dPA = mat_dPA(vecMembers,vecMembers);
			vecAssemblyMean_dPA(intAssembly) = mean(matReal_dPA(matSelect));
			
			%get real OSI
			vecAssResp = sum(matSpikeCountsStim(vecMembers,:),1);
			sTuning = calcTuningRespMat(vecAssResp,cellSelectStimO,vecOrientations);
			vecAssembly_OSI(intAssembly) = sTuning.vecOSI;
			
			%get shuffled dPAs & OSIs
			for intIter=1:intIters
				vecShuffledMembers = randperm(intNeurons,intMembers);
				matShuffled_dPA = mat_dPA(vecShuffledMembers,vecShuffledMembers);
				matAssemblyMean_dPA_Shuffled(intIter,intAssembly) = mean(matShuffled_dPA(matSelect));
				
				%calc OSI
				vecAssRespShuf = sum(matSpikeCountsStim(vecShuffledMembers,:),1);
				sTuning = calcTuningRespMat(vecAssRespShuf,cellSelectStimO,vecOrientations);
				matAssembly_OSI_Shuffled(intIter,intAssembly) = sTuning.vecOSI;
			end
		end
		
		%calc data
		dblAlpha = 0.01;
		intLowerThreshold = round(intIters*(dblAlpha/2));
		intUpperThreshold = round(intIters*(1-dblAlpha/2));
		
		matShuffledSorted = sort(matAssemblyMean_dPA_Shuffled,1);
		vecUpper = matShuffledSorted(intUpperThreshold,:);
		vecLower = matShuffledSorted(intLowerThreshold,:);
		vecMean = mean(matShuffledSorted,1);
		
		matShuffledSortedOSI = sort(matAssembly_OSI_Shuffled,1);
		vecUpperOSI = matShuffledSortedOSI(intUpperThreshold,:);
		vecLowerOSI = matShuffledSortedOSI(intLowerThreshold,:);
		vecMeanOSI = mean(matShuffledSortedOSI,1);
		
		
		%% slow/fast/hit/miss correlations of individual assemblies
		vecHitMissRatioPresence = nan(1,intAssemblies);
		vecSlowFastRatioPresence = nan(1,intAssemblies);
		for intAssembly=1:intAssemblies
			dblMissPresence = sum(matAssemblyActivity(intAssembly,indMissHD))/sum(indMissHD);
			dblFastPresence = sum(matAssemblyActivity(intAssembly,indFastHD))/sum(indFastHD);
			dblSlowPresence = sum(matAssemblyActivity(intAssembly,indSlowHD))/sum(indSlowHD);
			dblHitPresence = sum(matAssemblyActivity(intAssembly,indSlowHD|indFastHD))/sum(indSlowHD|indFastHD);
			
			vecHitMissRatioPresence(intAssembly) = dblHitPresence/dblMissPresence;
			vecSlowFastRatioPresence(intAssembly) = dblFastPresence/dblSlowPresence;
		end
		
		figure
		subplot(2,2,1);
		errorbar(1:intAssemblies,vecMean,vecLower-vecMean,vecUpper-vecMean,'xb');
		hold on
		scatter(1:intAssemblies,vecAssemblyMean_dPA,'xr');
		hold off
		xlabel('Assembly');
		ylabel('dPA');
		
		subplot(2,2,2);
		scatter(vecAssemblyMean_dPA(indRealAssemblies),vecSlowFastRatioPresence(indRealAssemblies));
		[dblCorr,dblP] = corr(vecAssemblyMean_dPA(indRealAssemblies)',vecSlowFastRatioPresence(indRealAssemblies)');
		xlabel('Mean dPA');
		ylabel('Fast/Slow presence ratio');
		title(sprintf('Correlation; r=%.3f; p=%.3f',dblCorr,dblP));
		
		subplot(2,2,3);
		errorbar(1:intAssemblies,vecMeanOSI,vecLowerOSI-vecMeanOSI,vecUpperOSI-vecMeanOSI,'xb');
		hold on
		scatter(1:intAssemblies,vecAssembly_OSI,'xr');
		hold off
		xlabel('Assembly');
		ylabel('OSI');
		
		
		subplot(2,2,4);
		scatter(vecAssemblyMean_dPA(indRealAssemblies),vecLaggingRatio(indRealAssemblies));
		[dblCorr,dblP] = corr(vecAssemblyMean_dPA(indRealAssemblies)',vecLaggingRatio(indRealAssemblies)');
		xlabel('Mean dPA');
		ylabel('Lagging ratio');
		title(sprintf('Correlation; r=%.3f; p=%.3f',dblCorr,dblP));
		
		if sParams.boolSavePlots
			drawnow;
			jFig = get(handle(gcf), 'JavaFrame');
			jFig.setMaximized(true);
			figure(gcf);
			drawnow;
			export_fig(sprintf('%spop%dorientation_tuning_assemblies.tif',strSes,intPopulation));
			export_fig(sprintf('%spop%dorientation_tuning_assemblies.pdf',strSes,intPopulation));
		end
		
		% other properties
		figure
		subplot(2,2,1);
		scatter(vecStimDrivenRatio(indRealAssemblies),vecAssOSI(indRealAssemblies));
		[dblCorr,dblP] = corr(vecStimDrivenRatio(indRealAssemblies)',vecAssOSI(indRealAssemblies)');
		xlabel('Stim-driven ratio');
		ylabel('OSI');
		title(sprintf('Correlation; r=%.3f; p=%.3f',dblCorr,dblP));
		
		subplot(2,2,2);
		scatter(vecLaggingRatio(indRealAssemblies),vecSlowFastRatioPresence(indRealAssemblies));
		[dblCorr,dblP] = corr(vecLaggingRatio(indRealAssemblies)',vecSlowFastRatioPresence(indRealAssemblies)');
		xlabel('Lagging ratio');
		ylabel('Fast/Slow presence ratio');
		title(sprintf('Correlation; r=%.3f; p=%.3f',dblCorr,dblP));
		
		subplot(2,2,3);
		scatter(vecLaggingRatio(indRealAssemblies),vecStimDrivenRatio(indRealAssemblies));
		[dblCorr,dblP] = corr(vecLaggingRatio(indRealAssemblies)',vecStimDrivenRatio(indRealAssemblies)');
		xlabel('Lagging ratio');
		ylabel('Stim-driven ratio');
		title(sprintf('Correlation; r=%.3f; p=%.3f',dblCorr,dblP));
		
		subplot(2,2,4);
		scatter(vecSlowFastRatioPresence(indRealAssemblies),vecStimDrivenRatio(indRealAssemblies));
		[dblCorr,dblP] = corr(vecStimDrivenRatio(indRealAssemblies)',vecSlowFastRatioPresence(indRealAssemblies)');
		xlabel('Fast/Slow presence ratio');
		ylabel('Stim-driven ratio');
		title(sprintf('Correlation; r=%.3f; p=%.3f',dblCorr,dblP));
		
		if sParams.boolSavePlots
			drawnow;
			jFig = get(handle(gcf), 'JavaFrame');
			jFig.setMaximized(true);
			figure(gcf);
			drawnow;
			export_fig(sprintf('%spop%dassembly_property_correlations.tif',strSes,intPopulation));
			export_fig(sprintf('%spop%dassembly_property_correlations.pdf',strSes,intPopulation));
		end
		
		%save data
		cellSaveAssemblyProperties{intPopulation,1} = vecAssemblyMean_dPA;
		cellSaveAssemblyProperties{intPopulation,2} = matAssemblyMean_dPA_Shuffled;
		cellSaveAssemblyProperties{intPopulation,3} = vecLaggingRatio;
		cellSaveAssemblyProperties{intPopulation,4} = vecHitMissRatioPresence;
		cellSaveAssemblyProperties{intPopulation,5} = vecSlowFastRatioPresence;
		cellSaveAssemblyProperties{intPopulation,6} = vecAssOSI;
		cellSaveAssemblyProperties{intPopulation,7} = vecStimDrivenRatio;
		cellSaveAssemblyProperties{intPopulation,8} = double(indRealAssemblies');
		cellSaveAssemblyProperties{intPopulation,9} = vecAssembly_OSI;
		cellSaveAssemblyProperties{intPopulation,10} = matAssembly_OSI_Shuffled;
		
	end
	
	%% save data structures
	if sParams.boolSaveData
		vecClock = fix(clock);
		strDate = num2str(vecClock(1:3));
		strFile = ['dataAssAnalTO' strSes '_' strrep(strDate,'    ','_')];
		save(['D:\Data\Results\spikeAnalysis\' strFile],'cellSave*','-v7.3');
		cd(strOldDir);
	end
end