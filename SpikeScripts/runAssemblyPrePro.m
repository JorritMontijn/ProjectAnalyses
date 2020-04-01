%pre-process data
for intMouse=1:8
	close all
	drawnow;
	clearvars -except intMouse
	
	%use neuropil subtraction?
	boolUseNeuropilSubtraction = false;
	boolExcludeLocomotor = false;
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
		fprintf('Loading pre-processed data for %s [%s]\n',strSes,getTime);
		load(['D:\Data\Results\spikeAnalysis\dataPreProAggregate' strSes '.mat']);
	end
	if boolOnlyTuned
		strSes = [strSes 'OT'];
	end
	
	%vecBlock = [1 1 1 2 2 2 2 2]; %define whether neurons are in same population or not
	vecBlockTypes = unique(vecBlock);
	intNumBlocks = length(vecBlockTypes);
	%vecNeuronNum = zeros(1,intNumBlocks);
	%cellKeepList = cell(1,intNumBlocks);
	%#ok<*ASGLU>
	%#ok<*AGROW>
	clear sLoad sSesAggregate ses
	sParams.strFigDir = ['D:\Data\Results\spikeAnalysis' filesep strSes filesep];
	sParams.boolSavePlots = true;
	sParams.boolSaveData = true;
	strOldDir = cd(sParams.strFigDir);
	
	%change name for no split
	strSes = ['AE' strSes];
	
	for intPopulation = vecBlockTypes
		%% GET RESPONSES AFTER COMBINING BLOCKS OF SAME POPULATION INTO ONE
		fprintf('Starting %s pop %d [%s]\n',strSes,intPopulation,getTime);
		
		%get neuronal tuning
		if intMouse==8
			%remove last trial
			%intLastFrame = cellMultiSes{1}.structStim.FrameOn(end-47);
			vecRem = true(size(cellMultiSes{1}.structStim.Orientation));
			vecRem((end-47):end) = false;
			cellMultiSes{1}.structStim = remel(cellMultiSes{1}.structStim,vecRem);
			%remove data
			%for intNeuron=1:numel(cellMultiSes{intPopulation}.neuron)
			%	cellMultiSes{intPopulation}.neuron(intNeuron).dFoF(intLastFrame:end) = [];
			%end
			%cellMultiSes{1} = doRecalcdFoF(cellMultiSes{1},3);
		end
		
		%reperform spike detection
		dblSamplingFreq = cellMultiSes{intPopulation}.structStim.FrameOff(end)/cellMultiSes{intPopulation}.structStim.SecsOff(end);
		sStimTemp = cellMultiSes{intPopulation}.structStim;
		dblStimSecs = 2;
		sStimTemp.FrameOff = sStimTemp.FrameOn + round(dblSamplingFreq*dblStimSecs);
		intNeurons = numel(cellMultiSes{intPopulation}.neuron);
		cellMultiSes{intPopulation}.samplingFreq = dblSamplingFreq;
		dblTau = 0.5;
		intTotDur = length(cellMultiSes{intPopulation}.neuron(1).dFoF);
		dblTotDurSecs = intTotDur/dblSamplingFreq;
		matSpikeCountsHD = nan(intNeurons,intTotDur);
		for intNeuron=1:intNeurons
			%msg
			fprintf('Performing detection for neuron %d/%d...\n',intNeuron,intNeurons)
			ptrTime = tic;

			%AE detection
			[apFrames, apSpikes, vecSpikes, expFit] = doDetectSpikes(cellMultiSes{intPopulation}.neuron(intNeuron).dFoF,dblSamplingFreq,dblTau,10000);
			intTotSpikes = sum(apSpikes);
			
			%assign to neuron
			cellMultiSes{intPopulation}.neuron(intNeuron).apFrames = apFrames;
			cellMultiSes{intPopulation}.neuron(intNeuron).apSpikes = apSpikes;
			cellMultiSes{intPopulation}.neuron(intNeuron).vecSpikes = vecSpikes;
			cellMultiSes{intPopulation}.neuron(intNeuron).expFit = expFit;
			matSpikeCountsHD(intNeuron,:) = vecSpikes;
			
			%msg
			fprintf('\b	Done! Took %.1f seconds; %d transients; mean rate is %.2f Hz\n',toc(ptrTime),intTotSpikes,intTotSpikes/dblTotDurSecs);
			pause(eps);
		end
		matSpikeCounts = getNeuronResponseRespMat(matSpikeCountsHD,sStimTemp);
		matSpikeCounts = matSpikeCounts*round(dblSamplingFreq);
		
		%[indTuned,vecNeuronPrefStim,vecOSI,cellAggregate] = getTunedStimDetectionNeurons(cellSes(vecBlock==intPopulation));
		[indTuned,vecNeuronPrefStim,vecOSI,cellAggregate] = getTunedStimDetectionNeurons(cellMultiSes{intPopulation},matSpikeCounts);
		if ~boolOnlyTuned
			indTuned = true(size(indTuned)); %include all neurons
		end
		structStim = cellMultiSes{intPopulation}.structStim;
		matSpikeCounts(~indTuned,:) = [];
		matSpikeCountsHD(~indTuned,:) = [];
		
		
		% remove trials with reaction time <100ms
		dblRemSecs = 0.1;
		indTooFast = (structStim.vecTrialRespSecs-structStim.SecsOn)<dblRemSecs & structStim.vecTrialResponse == 1;
		fprintf('Removed %d trials with responses <%dms\n',sum(indTooFast),round(dblRemSecs*1000));
		%structStim.vecTrialResponse(indTooFast) = 0;
		structStim = remel(structStim,~indTooFast);
		%structStim.FrameOff = structStim.FrameOn+1;
		
		%remove trials that were too slow
		dblRemSecs = 3;
		indTooSlow = (structStim.vecTrialRespSecs-structStim.SecsOn)>dblRemSecs;
		fprintf('Removed %d trials with responses >%dms\n',sum(indTooSlow),round(dblRemSecs*1000));
		structStim.vecTrialResponse(indTooSlow) = 0;
		structStim.vecTrialRespSecs(indTooSlow) = nan;
		%structStim.FrameOff = structStim.FrameOn+1;
		
		%remove locomotion trials
		if boolExcludeLocomotor
			indLoco = false(1,numel(structStim.Orientation));
			for intTrial=1:numel(structStim.Orientation)
				vecMoveTimes = structStim.cellMoveSecs{intTrial};
				indLoco(intTrial) = sum(vecMoveTimes > structStim.SecsOn(intTrial) & vecMoveTimes < structStim.SecsOff(intTrial)) > 0;
			end
			structStim = remel(structStim,~indLoco);
			fprintf('Removed %d trials with locomotion\n',sum(indLoco));
		end
		
		%take opposite directions as the same
		structStim.Orientation = mod(structStim.Orientation,180);
		vecOrientations = unique(structStim.Orientation);
		cellMultiSes{intPopulation}.structStim = structStim;
		
		%get orientation-based trial selection vectors
		sTypesOri = getStimulusTypes(structStim,{'Orientation'});
		cellSelectOri = getSelectionVectors(structStim,sTypesOri);
		
		%remove non-tuned neurons
		cellMultiSes{intPopulation}.neuron(~indTuned) = [];
		vecNeuronPrefStim(~indTuned) = [];
		intNeurons = numel(cellMultiSes{intPopulation}.neuron);
		intTrials = length(structStim.Orientation);
		intOris = length(vecOrientations);
		
		%%
		%{
		search for frames where multiple neurons are
		active; or where other neurons are active in the subsequent frame;
		keep adding frames until no neurons are active at a certain point;
		then call this an assembly event and save the event with which
		neurons are active at which point
		
		matNeuronsInAssembly = [];
		%}
		
		fprintf('Retrieving population spiking data for assembly formation analysis...   [%s]\n',getTime);
		
		%get number of active neurons per frame
		vecNrActive=sum(matSpikeCountsHD>0,1);
		
		%select which frames are potentially interesting
		vecPotentialAssemblyFrames = vecNrActive > 0;
		
		%remove frames where single neurons fire single spikes in temporal isolation
		vecSingleSpikes = vecNrActive(1:(end-2))==0 & vecNrActive(2:(end-1))<2 & vecNrActive(2:(end-1))>0 & vecNrActive(3:end)==0;
		vecPotentialAssemblyFrames([false vecSingleSpikes false]) = 0;
		
		%get frames where multiple neurons are firing in temporal isolation
		vecTransientAssemblies = vecNrActive(1:(end-2))==0 & vecNrActive(2:(end-1))>2 & vecNrActive(3:end)==0;
		vecPotentialAssemblyFrames([false vecTransientAssemblies false]) = 1;
		
		%define start frames of assemblies
		vecAssemblyStarts = [false vecPotentialAssemblyFrames(1:(end-1))==0 & vecPotentialAssemblyFrames(2:(end))>0];
		vecAssemblyStops = [false vecPotentialAssemblyFrames(1:(end-1))>0 & vecPotentialAssemblyFrames(2:end)==0];
		
		%combine
		vecAssemblyStarts = find(vecAssemblyStarts);
		vecAssemblyStops = find(vecAssemblyStops);
		if vecAssemblyStarts(1) > vecAssemblyStops(1)
			vecAssemblyStarts = [1 vecAssemblyStarts];
		end
		[dummy,vecReorder] = sort(vecAssemblyStarts,'ascend');
		if length(vecAssemblyStops) < length(vecAssemblyStarts)
			vecAssemblyStops(end+1) = length(vecPotentialAssemblyFrames);
		end
		vecAssemblyStarts = vecAssemblyStarts(vecReorder);
		vecAssemblyStops = vecAssemblyStops(vecReorder);
		
		%get assemblies from spike matrix + assembly stop/starts
		matAssemblies = getAssemblies(matSpikeCountsHD,vecAssemblyStarts,vecAssemblyStops)>0;
		
		%remove assemblies with 2 or fewer neurons
		indRealAssemblies = sum(matAssemblies,1)>=5 & sum(matAssemblies,1) < intNeurons;
		matAssemblies = matAssemblies(:,indRealAssemblies);
		vecAssemblyStarts = vecAssemblyStarts(indRealAssemblies);
		vecAssemblyStops = vecAssemblyStops(indRealAssemblies);
		
		%get consistencies
		matNeuronCorrelations = corr(matAssemblies');
		matAssemblyCorrelations = corr(matAssemblies);
		
		%cluster assemblies
		matSelect = tril(true(size(matAssemblyCorrelations)),0) & triu(true(size(matAssemblyCorrelations)),0);
		matDistAN = 1-matAssemblyCorrelations;
		matDistAN(matSelect) = 0;
		[intNumberOfAssemblies,vecSilhouetteDistances] = doFastClustering(matDistAN,min([100 round(sqrt(length(matDistAN))*2)]));
		fprintf('Optimal number of assemblies: %d. Proceeding with assembly assignment... [%s]\n',intNumberOfAssemblies,getTime);
		
		%re-cluster with detected number of assemblies
		matLinkageA = linkage(matDistAN,'ward');
		vecA = cluster(matLinkageA,'maxclust',intNumberOfAssemblies);
		[vecAssemblyIdentity,vecReorderA] = sort(vecA,'ascend');
		
		%return
		intT = length(0:(1/dblSamplingFreq):(cellMultiSes{intPopulation}.structStim.SecsOff(end)+5));
		matAssemblyActivity = zeros(intNumberOfAssemblies,intT);
		matAssemblyActivity = putAssemblies(matAssemblyActivity,vecA,vecAssemblyStarts,vecAssemblyStops);

		if sParams.boolSavePlots
			fprintf('Assemblies assigned. Creating figure and dendrogram... [%s]\n',getTime);
		
			figure
			subplot(2,2,1)
			matLinkageAN = linkage(matDistAN,'ward');
			
			[H,T,vecClusterReorderAN] = dendrogram(matLinkageAN,0);
			drawnow;
			
			subplot(2,2,2)
			vecFiltSil = conv(vecSilhouetteDistances,normpdf(-2:2,0,0.8),'same');
			plot(vecSilhouetteDistances);
			xlabel('Number of clusters');
			ylabel('Silhouette distance');
			title(sprintf('Number of clusters that is optimal; %d [val=%.3f]',intNumberOfAssemblies,vecSilhouetteDistances(intNumberOfAssemblies)));
			
			subplot(2,2,3)
			imagesc(matAssemblyCorrelations(vecReorderA,vecReorderA),[-1 1]);colormap('redblue');freezeColors;colorbar;cbfreeze;
			title('Assembly correlations, ordered by cluster')
			
			subplot(2,4,7)
			imagesc(vecAssemblyIdentity);colormap('jet');freezeColors;colorbar;cbfreeze;
			title('Cluster membership')
			
			subplot(2,4,8)
			matM = getBlockMeans(matAssemblyCorrelations(vecReorderA,vecReorderA),vecAssemblyIdentity);
			vecAssemblyAutoCorr = diag(matM);
			scatter(1:length(vecAssemblyAutoCorr),vecAssemblyAutoCorr');
			ylim([0 1])
			title('Assembly consistency')
			
			drawnow;
			jFig = get(handle(gcf), 'JavaFrame');
			jFig.setMaximized(true);
			figure(gcf);
			drawnow;
			strFig = sprintf('%sagg_assembly_correlations_pop%d_raw',strSes,intPopulation);
			export_fig([sParams.strFigDir strFig '.tif']);
			export_fig([sParams.strFigDir strFig '.pdf']);
		end
		
		%% put in output
		matRespAssemblies = getNeuronResponseRespMat(matAssemblyActivity,cellMultiSes{intPopulation}.structStim);
		sAssemblies{intPopulation}.matAssemblies = matAssemblies;
		sAssemblies{intPopulation}.vecAssemblies = vecA;
		sAssemblies{intPopulation}.vecAssemblyStarts = vecAssemblyStarts;
		sAssemblies{intPopulation}.vecAssemblyStops = vecAssemblyStops;
		sAssemblies{intPopulation}.matAssemblyActivity = matAssemblyActivity;
		sAssemblies{intPopulation}.matAssemblyCorrelations = matAssemblyCorrelations;
		sAssemblies{intPopulation}.vecAssemblyIdentity = vecAssemblyIdentity;
		sAssemblies{intPopulation}.vecReorder = vecReorderA;
		sAssemblies{intPopulation}.matRespAssemblies = matRespAssemblies;
		sAssemblies{intPopulation}.matSpikeCountsHD = matSpikeCountsHD;
	end
	
	%% save data
	fprintf('Saving data for mouse %d (%s)... [%s]\n',intMouse,strSes,getTime);
	strPath = 'D:\Data\Results\spikeAnalysis\';
	strFile = ['dataPreProAssemblies' strSes(3:end) ];
	save([strPath strFile],'cellMultiSes','sAssemblies','-v7.3');
	fprintf('Assembly preprocessing completed! [%s]\n',getTime);
end