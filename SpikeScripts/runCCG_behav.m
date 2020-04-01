%pre-process data
for intMouse=1:8
	close all
	drawnow;
	clearvars -except intMouse
	%get block data
	if ~exist('cellMultiSes','var')
		if intMouse == 1
			strSes = '20140207';
		elseif intMouse == 2
			strSes = '20140314';
		elseif intMouse == 3
			strSes = '20140425';
		elseif intMouse == -1
			%return; %exclude?; bad behavior, weird signals
			strSes = '20140430';
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
		fprintf('Loading pre-processed data for %s [%s]\n',strSes,getTime);
		load(['D:\Data\Results\spikeAnalysis\dataPreProAggregate' strSes '.mat']);
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
		
		%% pre-pro to spikes
		fprintf('Transforming dF/F0 data to activation events [%s]\n',getTime);
		[cellSpikes,structStimSpikes] = getSpikeArrays(cellMultiSes{intPopulation});
		cellMultiSes{intPopulation}.structStim = structStimSpikes;
		fprintf('Transformation complete; calculating preferred stimulus orientations on spikes [%s]\n',getTime);
		matSpikeCounts = getSpikeCounts(cellSpikes,structStimSpikes.TimeOn,structStimSpikes.TimeOff);
		
		%take opposite directions as the same
		cellMultiSes{intPopulation}.structStim.Orientation = mod(cellMultiSes{intPopulation}.structStim.Orientation,180);
		vecOrientations = unique(cellMultiSes{intPopulation}.structStim.Orientation);
		intNeurons = numel(cellMultiSes{intPopulation}.neuron);
		
		%get neuronal responses per trial
		sTypesC = getStimulusTypes(cellMultiSes{intPopulation}.structStim,{'Contrast'});
		cellSelectContrasts = getSelectionVectors(cellMultiSes{intPopulation}.structStim,sTypesC);
		intContrasts = length(cellSelectContrasts);
		
		%get orientation-based trial selection vectors
		sTypesOri = getStimulusTypes(cellMultiSes{intPopulation},{'Orientation'});
		cellSelectOri = getSelectionVectors(cellMultiSes{intPopulation}.structStim,sTypesOri);
		
		%make normalized activity matrix per contrast
		matTrialNormSpikes = zeros(size(matSpikeCounts));
		for intContrast=1:intContrasts
			matTrialNormSpikes(:,cellSelectContrasts{intContrast}) = imnorm(matSpikeCounts(:,cellSelectContrasts{intContrast}),1);
		end
		
		%calculate preferred stimulus orientation for all neurons
		vecOriPrefContrasts = 4:6; % 8% - 100%
		matPref = nan(length(vecOriPrefContrasts),intNeurons);
		intCounter = 0;
		for intContrast = vecOriPrefContrasts
			intCounter = intCounter + 1;
			matResp = matTrialNormSpikes(:,cellSelectContrasts{intContrast});
			structStimC{intContrast} = cellMultiSes{intPopulation}.structStim;
			cellFields = fieldnames(cellMultiSes{intPopulation}.structStim);
			for intField=1:length(cellFields)
				strField = cellFields{intField};
				structStimC{intContrast}.(strField) = structStimC{intContrast}.(strField)(cellSelectContrasts{intContrast});
			end
			cellSelect = getSelectionVectors(structStimC{intContrast},sTypesOri);
			sTuning{intContrast} = calcTuningRespMat(matResp,cellSelect,vecOrientations);
			matPref(intCounter,:) = sTuning{intContrast}.vecPrefIndex;
		end
		vecNeuronPrefStim = nan(1,intNeurons);
		for intOri=1:length(vecOrientations);
			vecNeuronPrefStim(sum(matPref == intOri,1) > 1) = intOri;
		end
		%vecNeuronPrefStim(isnan(vecNeuronPrefStim)) = randi(length(vecOrientations),[1 sum(isnan(vecNeuronPrefStim))]);
		
		%remove non-tuned neurons
		cellMultiSes{intPopulation}.neuron(isnan(vecNeuronPrefStim)) = [];
		cellSpikes(isnan(vecNeuronPrefStim)) = [];
		vecNeuronPrefStim(isnan(vecNeuronPrefStim)) = [];
		intNeurons = numel(cellMultiSes{intPopulation}.neuron);
		intTrials = length(cellMultiSes{intPopulation}.structStim.Orientation);
		
		fprintf('Pref stim calculation complete; decoding...   [%s]\n',getTime);
		
		%% do analysis
		% decoding
		
		%TO DO: decoding over contrasts
		sDec_Spikes = doCrossValidationRespMat(matSpikeCounts,cellMultiSes{intPopulation}.structStim);
		sDec_dFoF = doCrossValidation(cellMultiSes{intPopulation});
		dblTotSpikesCorr = sum(sDec_Spikes.vecDecodedStimType == sDec_Spikes.vecStimType);
		dblTotdFoFCorr = sum(sDec_dFoF.vecDecodedStimType == sDec_dFoF.vecStimType);
		
		[phat_S,pci_S] = binofit(dblTotSpikesCorr,length(sDec_Spikes.vecDecodedStimType));
		[phat_D,pci_D] = binofit(dblTotdFoFCorr,length(sDec_dFoF.vecDecodedStimType));
		
		figure
		subplot(2,2,1)
		
		vecX = [0.8 1.2];
		errorbar(vecX(1),phat_S,phat_S-pci_S(1),phat_S-pci_S(2),'ro');
		hold on;
		errorbar(vecX(2),phat_D,phat_D-pci_D(1),phat_D-pci_D(2),'bo');
		plot(get(gca,'XLim'),[0.25 0.25],'k--');
		hold off;
		set(gca,'XTick',vecX,'XTickLabel',{'Activation Events','dF/F0'});
		vecLimY = get(gca,'YLim');
		vecLimY(1) = 0;
		ylim(vecLimY);
		ylabel('Decoding accuracy');
		title('Decoding performance (mean +/- 95% CI)');
		drawnow;
		
		
		%{
		%pop-act based CCGs
		dblStep = 1/cellMultiSes{intPopulation}.samplingFreq;
		vecCCG_bins = -1:dblStep:1;
		matCCG = zeros(intNeurons,length(vecCCG_bins));
		intNeurons = length(cellSpikes);
		
		vecCCG_bins = -1:dblStep:1;
		for intNeuron=1:intNeurons
			fprintf('Now at neuron %d/%d\n',intNeuron,intNeurons);
			vecSpikeTimes = cellSpikes{intNeuron};
			for intSpike=1:length(vecSpikeTimes)
				vecStart = vecCCG_bins + vecSpikeTimes(intSpike);
				matThisCCG = getSpikeCounts(cellSpikes,vecStart,dblStep);
				matThisCCG(intNeuron,:)=0;
				matCCG = matCCG + matThisCCG;
			end
		end
		%matCCG(end+1,:) = sum(matCCG,1);
		
		%plot mean pop activity
		%{
figure
subplot(2,2,1)
if exist('matSpikeCountsHiRes','var')
	vecPopAct = sum(matSpikeCountsHiRes,1);
	vecBlur = normpdf(-4:4,0,2);
	vecBlur = vecBlur./sum(vecBlur);
	vecPopAct = conv(vecPopAct,vecBlur,'same');
	plot(vecPopAct)
	xlim([0 dblTotTime/dblStep])
	title('Population activity (filtered event rate)')
end
		%}
		
		%try taking only events where pop act is above shuffled spiking prob
		subplot(2,2,2)
		vecBlur = 1;
		matCCG2 = conv2(matCCG,vecBlur,'valid');
		matCCG2 = imnorm(matCCG2,1);
		imagesc(matCCG2);
		colormap('hot')
		title('Cross-correlations normalized to max=1, sorted by random')
		
		subplot(2,2,3)
		matCCG3 = matCCG2./repmat(sum(matCCG2,2),[1 size(matCCG2,2)]);
		
		vecMuCC = (sum((matCCG.*repmat(1:size(matCCG,2),[size(matCCG,1) 1])),2) ./ sum(matCCG,2) - ceil(length(vecCCG_bins)/2)) * dblStep;
		[dummy,vecReorder] = sort(vecMuCC,'ascend');
		%[dummy,vecReorder] = sort(sum(cumsum(matCCG3,2)>0.5,2),'descend');
		
		imagesc(cumsum(matCCG3(vecReorder,:),2))
		colormap('hot')
		title('Cumulative cross-correlations normalized to sum=1, sorted by Mu_c_c')
		
		subplot(2,2,4)
		imagesc(matCCG2(vecReorder,:))
		colormap('hot')
		title('Cross-correlations normalized to max=1, sorted by Mu_c_c')
		
		if sParams.boolSavePlots
			drawnow;
			jFig = get(handle(gcf), 'JavaFrame');
			jFig.setMaximized(true);
			figure(gcf);
			drawnow;
			strFig = sprintf('%sagg_crosscorrelations_decoding_pop%d_raw',strSes,intPopulation);
			export_fig([strFig '.tif']);
			export_fig([strFig '.pdf']);
		end
		
		
		figure
		subplot(2,2,1)
		dblMaxWithoutAuto = max(max(matCCG3));
		imagesc(matCCG3(vecReorder,:),[0 dblMaxWithoutAuto]);
		colormap('hot');
		freezeColors;
		title('Sorted by Mu_cc')
		subplot(2,2,2)
		matCorr = corr(matCCG3(vecReorder,:)');
		imagesc(matCorr,[-1 1]);
		colormap('redblue');
		freezeColors;
		title('Sorted by Mu_c_c')
		
		%sort by mean corr
		[dummy,vecReorder2] = sort(mean(matCorr),'descend');
		
		subplot(2,2,3)
		imagesc(matCCG3(vecReorder(vecReorder2),:),[0 dblMaxWithoutAuto]);
		colormap('hot');
		freezeColors;
		title('Sorted by mean cross-correlation')
		
		subplot(2,2,4)
		imagesc(matCorr(vecReorder2,vecReorder2),[-1 1]);
		colormap('redblue');
		title('Sorted by mean cross-correlation')
		freezeColors;
		
		if sParams.boolSavePlots
			drawnow;
			jFig = get(handle(gcf), 'JavaFrame');
			jFig.setMaximized(true);
			figure(gcf);
			drawnow;
			strFig = sprintf('%sagg_crosscorrelations_Mu_cc_pop%d_raw',strSes,intPopulation);
			export_fig([strFig '.tif']);
			export_fig([strFig '.pdf']);
		end
		
		%% single-neuron stuff
		%single-neuron CCGs
		vecCCG_bins = -dblStep:dblStep:dblStep;
		intNeurons = length(cellSpikes);
		matSynchrony = nan(intNeurons,intNeurons);
		intCombinations = ((intNeurons^2)-intNeurons)/2;
		matPairwiseSynchrony = nan(intCombinations,2);
		matPairwisePairCombo = nan(intCombinations,2);
		intPairCounter = 0;
		for intNeuron1=1:intNeurons
			%get spike times
			vecSpikeTimesNeuron1 = cellSpikes{intNeuron1};
			fprintf('Now at neuron %d/%d\n',intNeuron1,intNeurons);
			
			%calculate CCG for this neuron with all others
			for intNeuron2=1:intNeurons
				vecCCG = zeros(1,length(vecCCG_bins));
				vecSpikeTimesNeuron2 = cellSpikes{intNeuron2};
				for intSpike=1:length(vecSpikeTimesNeuron1)
					vecStart = vecCCG_bins + vecSpikeTimesNeuron1(intSpike);
					vecCCG = vecCCG + getSpikeCounts(vecSpikeTimesNeuron2,vecStart,dblStep);
				end
				
				%normalize synchronicity value by number of spikes of neuron 2
				dblSynchrony = sum(vecCCG)/(length(vecSpikeTimesNeuron2));
				
				%put in output
				matSynchrony(intNeuron1,intNeuron2) = dblSynchrony;
				
				if intNeuron1 == intNeuron2
					%skip
				else
					%get locations
					vecSortedNeuronPair = sort([intNeuron1 intNeuron2],'ascend');
					intLowerNeuron = vecSortedNeuronPair(1);
					intHigherNeuron = vecSortedNeuronPair(2);
					if intNeuron1 == intLowerNeuron
						intLoc = 1;
					elseif intNeuron2 == intLowerNeuron
						intLoc = 2;
					end
					intPairLocation = sum((intNeurons-1):-1:(intNeurons - intLowerNeuron + 1)) + intHigherNeuron - intLowerNeuron;
					matPairwisePairCombo(intPairLocation,:) = [intLowerNeuron; intHigherNeuron];
					
					%get data & assign
					matPairwiseSynchrony(intPairLocation,intLoc) = dblSynchrony;
				end
			end
		end
		
		figure
		subplot(2,2,1)
		matSelect = tril(true(size(matSynchrony)),0) & triu(true(size(matSynchrony)),0);
		matSynchrony(matSelect) = 0;
		imagesc(matSynchrony)
		colormap('hot')
		colorbar
		title('Synchrony')
		
		%reorder by synchrony
		[dummy,vecReorderS] = sort(sum(matSynchrony,1)+sum(matSynchrony,2)','ascend');
		subplot(2,2,2)
		matReordered = matSynchrony(vecReorderS,vecReorderS);
		imagesc(matReordered)
		colormap('hot')
		colorbar
		title('Synchrony, ordered by mean synchony')
		
		%try clustering with ward algorithm, then compute silhouette
		%distance for several depths (up to ~5), then take max of distances
		%with those clusters and cluster means to k-means cluster with
		%those cluster centroids
		matDist = 1-matSynchrony;
		matDist(matSelect) = 0;
		
		
		matLinkage = linkage(matDist,'ward');
		
		subplot(4,2,5)
		[H,T,vecClusterReorder] = dendrogram(matLinkage,0,'colorthreshold','default');
		drawnow;
		
		subplot(4,2,7)
		[intNumberOfClustersThatIsOptimal,vecSilhouetteDistances] = doClustering(matDist,10);
		
		plot(vecSilhouetteDistances);
		xlabel('Number of clusters');
		ylabel('Silhouette distance');
		title(sprintf('Number of clusters that is optimal; %d [val=%.3f]',intNumberOfClustersThatIsOptimal,vecSilhouetteDistances(intNumberOfClustersThatIsOptimal)));
		
		subplot(2,2,4)
		imagesc(matSynchrony(vecClusterReorder,vecClusterReorder));
		colormap('hot')
		colorbar
		title('Synchrony, ordered by cluster')
		
		
		%save plot
		if sParams.boolSavePlots
			drawnow;
			jFig = get(handle(gcf), 'JavaFrame');
			jFig.setMaximized(true);
			figure(gcf);
			drawnow;
			strFig = sprintf('%sagg_synchrony_clusters_pop%d_raw',strSes,intPopulation);
			export_fig([strFig '.tif']);
			export_fig([strFig '.pdf']);
		end
		%close all;drawnow;
		%}
		
		
		%%
		%{
		search for frames where multiple neurons are
		active; or where other neurons are active in the subsequent frame;
		keep adding frames until no neurons are active at a certain point;
		then call this an assembly event and save the event with which
		neurons are active at which point
		
		matNeuronsInAssembly = [];
		%}
		
		fprintf('Retrieving population spiking data for assembly formation analysis. This might take a while...   [%s]\n',getTime);
		
		% get mean population activity
		dblStep = 1/cellMultiSes{intPopulation}.samplingFreq;
		dblTotTime = cellMultiSes{intPopulation}.structStim.TimeOff(end)+5;
		matSpikeCountsHiRes = getSpikeCounts(cellSpikes,0:dblStep:dblTotTime,dblStep);
		fprintf('Population spiking data retrieved [%s]\n',getTime);
		
		%get number of active neurons per frame
		vecNrActive=sum(matSpikeCountsHiRes>0,1);
		
		%select which frames are potentially interesting
		vecPotentialAssemblyFrames = vecNrActive > 0;
		
		%remove frames where single neurons fire single spikes in temporal isolation
		vecSingleSpikes = vecNrActive(1:(end-2))==0 & vecNrActive(2:(end-1))==1 & vecNrActive(3:end)==0;
		vecPotentialAssemblyFrames([false vecSingleSpikes false]) = 0;
		
		%get frames where multiple neurons are firing in temporal isolation
		vecTransientAssemblies = vecNrActive(1:(end-2))==0 & vecNrActive(2:(end-1))>1 & vecNrActive(3:end)==0;
		vecPotentialAssemblyFrames([false vecTransientAssemblies false]) = 0;
		
		%define start frames of assemblies
		vecAssemblyStarts = vecPotentialAssemblyFrames(1:(end-2))==0 & vecPotentialAssemblyFrames(2:(end-1))>0 & vecPotentialAssemblyFrames(3:end)>0;
		vecAssemblyStops = vecPotentialAssemblyFrames(1:(end-2))>0 & vecPotentialAssemblyFrames(2:(end-1))>0 & vecPotentialAssemblyFrames(3:end)==0;
		
		%combine
		vecAssemblyStarts = find(vecAssemblyStarts);
		vecAssemblyStops = find(vecAssemblyStops);
		[dummy,vecReorder] = sort(vecAssemblyStarts,'ascend');
		vecAssemblyStarts = vecAssemblyStarts(vecReorder);
		vecAssemblyStops = vecAssemblyStops(vecReorder);
		
		%get assemblies from spike matrix + assembly stop/starts
		matAssemblies = getAssemblies(matSpikeCountsHiRes,vecAssemblyStarts,vecAssemblyStops);
		
		
		subplot(2,2,3)
		matNeuronCorrelations = corr(matAssemblies');
		imagesc(matNeuronCorrelations,[-1 1]);colormap('redblue');freezeColors;colorbar;cbfreeze;
		
		subplot(2,2,4)
		matAssemblyCorrelations = corr(matAssemblies);
		imagesc(matAssemblyCorrelations,[-1 1]);colormap('redblue');freezeColors;colorbar;cbfreeze;
		
		subplot(2,1,1)
		imagesc(matAssemblies>1);
		colormap('hot');
		drawnow;
		freezeColors;
		cbfreeze;
		
		%save plot
		if sParams.boolSavePlots
			drawnow;
			jFig = get(handle(gcf), 'JavaFrame');
			jFig.setMaximized(true);
			figure(gcf);
			drawnow;
			strFig = sprintf('%sagg_assembly_correlations_pop%d_raw',strSes,intPopulation);
			export_fig([strFig '.tif']);
			export_fig([strFig '.pdf']);
		end
		
		%cluster assemblies
		matSelect = tril(true(size(matAssemblyCorrelations)),0) & triu(true(size(matAssemblyCorrelations)),0);
		matDistAN = 1-matAssemblyCorrelations;
		matDistAN(matSelect) = 0;
		matLinkageAN = linkage(matDistAN,'ward');
		
		figure
		subplot(4,2,1)
		[H,T,vecClusterReorderAN] = dendrogram(matLinkageAN,0);
		drawnow;
		
		subplot(4,2,3)
		[intNumberOfAssemblies,vecSilhouetteDistances] = doFastClustering(matDistAN,round(sqrt(length(matDistAN))*2));
		
		plot(vecSilhouetteDistances);
		xlabel('Number of clusters');
		ylabel('Silhouette distance');
		title(sprintf('Number of clusters that is optimal; %d [val=%.3f]',intNumberOfAssemblies,vecSilhouetteDistances(intNumberOfAssemblies)));
		
		subplot(2,2,2)
		imagesc(matAssemblyCorrelations(vecClusterReorderAN,vecClusterReorderAN),[-1 1]);colormap('redblue');freezeColors;colorbar;cbfreeze;
		title('Assembly membership, ordered by cluster')
		
		
		
		%cluster neurons
		matSelect = tril(true(size(matNeuronCorrelations)),0) & triu(true(size(matNeuronCorrelations)),0);
		matDistNA = 1-matNeuronCorrelations;
		matDistNA(matSelect) = 0;
		matLinkageNA = linkage(matDistNA,'ward');
		
		
		subplot(4,2,5)
		[H,T,vecClusterReorderNA] = dendrogram(matLinkageNA,0);
		drawnow;
		
		subplot(4,2,7)
		[intNumberOfClustersThatIsOptimalNeurons,vecSilhouetteDistances] = doFastClustering(matDistNA,round(sqrt(length(matDistNA))*2));
		
		plot(vecSilhouetteDistances);
		xlabel('Number of clusters');
		ylabel('Silhouette distance');
		title(sprintf('Number of clusters that is optimal; %d [val=%.3f]',intNumberOfClustersThatIsOptimalNeurons,vecSilhouetteDistances(intNumberOfClustersThatIsOptimalNeurons)));
		
		subplot(2,2,4)
		imagesc(matNeuronCorrelations(vecClusterReorderNA,vecClusterReorderNA),[-1 1]);colormap('redblue');freezeColors;colorbar;cbfreeze;
		title('Assembly membership, ordered by cluster')
		
		%save plot
		if sParams.boolSavePlots
			drawnow;
			jFig = get(handle(gcf), 'JavaFrame');
			jFig.setMaximized(true);
			figure(gcf);
			drawnow;
			strFig = sprintf('%sagg_assembly_clusters_pop%d_raw',strSes,intPopulation);
			export_fig([strFig '.tif']);
			export_fig([strFig '.pdf']);
		end
		
		%re-cluster with detected number of assemblies
		matLinkageA = linkage(matDistAN,'ward');
		vecA = cluster(matLinkageA,'maxclust',intNumberOfAssemblies);
		[vecAssemblyIdentity,vecReorderA] = sort(vecA,'ascend');
		figure
		subplot(2,2,1)
		matReorderedAssemblyCorrs = matAssemblyCorrelations(vecReorderA,vecReorderA);
		imagesc(matReorderedAssemblyCorrs,[-1 1]);colormap('redblue');freezeColors;colorbar;cbfreeze;
		title('Assembly membership, ordered by cluster')
		matM = getBlockMeans(matReorderedAssemblyCorrs,vecAssemblyIdentity);
		vecRealClusters = abs(zscore(diag(matM)))<3;
		
		vecReorderedAssemblyStarts = vecAssemblyStarts(vecReorderA);
		vecReorderedAssemblyStops = vecAssemblyStops(vecReorderA);
		
		subplot(2,2,2)
		scatter(1:length(vecReorderedAssemblyStarts),vecReorderedAssemblyStarts)
		
		%save plot
		if sParams.boolSavePlots
			drawnow;
			jFig = get(handle(gcf), 'JavaFrame');
			jFig.setMaximized(true);
			figure(gcf);
			drawnow;
			strFig = sprintf('%sagg_assembly identities_pop%d_raw',strSes,intPopulation);
			export_fig([strFig '.tif']);
			export_fig([strFig '.pdf']);
		end
		
		%return
		intT = size(matSpikeCountsHiRes,2);
		matAssemblyActivity = zeros(intNumberOfAssemblies,intT);
		matAssemblyActivity = putAssemblies(matAssemblyActivity,vecA,vecAssemblyStarts,vecAssemblyStops);
		%figure
		%[dummy,matID] = meshgrid(1:size(matAssemblyActivity,2),1:size(matAssemblyActivity,1));
		%matPlot = matID.*matAssemblyActivity;
		%imagesc(matPlot);colormap('jet');nancolorbar(matPlot,[1 max(matPlot(:))]);
		
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
		sAssemblies{intPopulation}.sDec_Spikes = sDec_Spikes;
		sAssemblies{intPopulation}.sDec_dFoF = sDec_dFoF;

		%% possible analyses
		%decoding; which give most information?
		
		%sDec_Assemblies = doCrossValidationRespMat(matRespAssemblies,cellMultiSes{intPopulation}.structStim);
		
		%number of assemblies vs. heterogeneity
		
		%assembly occurrence over time relative to stimulus [PSTH] for
		%miss, fast, slow
		
		%relative occurrence of assemblies during miss/fast/slow
		
		%consistency of assemblies themselves during miss/fast/slow
		
		%size of assemblies during miss/fast/slow
		
		%orientation tuning of assemblies
		
		%contrast tuning of assemblies
		
		%correlation assemblies with pupil size
	end
	
	%% save data
	strPath = 'D:\Data\Results\spikeAnalysis\';
	strFile = ['dataPreProAssemblies' cellMultiSes{1}.mouse];
	save([strPath strFile],'cellMultiSes','sAssemblies','-v7.3');
end