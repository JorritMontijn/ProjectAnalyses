%load data
vecPEs = zeros(1,8);
vecRPEs = zeros(1,8);
vecNPEs = zeros(1,8);
vecTotDur = zeros(1,8);
vecActCorr = []
for intMouse=1:8
	close all
	drawnow;
	clearvars -except intMouse vecPEs vecActCorr
	
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
		load(['D:\Data\Results\spikeAnalysis\dataPreProAssemblies' strSes '.mat']);
	end
	
	clear sLoad sSesAggregate ses
	sParams.strFigDir = ['D:\Data\Results\spikeAnalysis' filesep strSes filesep];
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
		structParamsAss.intStartOffset = round(-dblStimSecs*cellMultiSes{intPopulation}.samplingFreq);
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
		
		vecActCorr(end+1) = corr(matSpikeCounts(:),matTrialResponse(:))
		return
		vecPEs(intMouse) = vecPEs(intMouse) + numel(sAssemblies{intPopulation}.vecAssemblies)
		continue
		
		%% get members of assemblies
		intNeurons = size(sAssemblies{intPopulation}.matAssemblies,1);
		matAssemblyMembers = sAssemblies{intPopulation}.matAssemblies > 0;
		intAssemblies = max(sAssemblies{intPopulation}.vecAssemblies);
		hAssemblyDetails=figure;
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
			subplot(2,2,1)
			imagesc(matRandMembership);colormap('hot');freezeColors;colorbar;cbfreeze;
			xlabel('Shuffle iteration');
			ylabel('Neuron ID');
			title(sprintf('Assembly %d membership probability',intAssembly));
			
			%mean(matCorr(tril(true(size(matCorr)),-1)))
			subplot(2,2,3)
			dblAlpha = 0.01;
			dblCorrectedSDs = qnorm(1-(dblAlpha/2)/intNeurons);
			matCorr = corr(matAssemblyMembers(:,indAssembly));
			imagesc(matCorr,[-1 1]);colormap(redblue);freezeColors;colorbar;cbfreeze;
			xlabel('Assembly occurrence');
			ylabel('Assembly occurrence');
			title('Assembly consistency (r neuron membership)');
			
			%get actual memberships
			subplot(2,2,2)
			imagesc(matAssemblyMembers(:,indAssembly));colormap('bw');freezeColors;
			xlabel('Assembly occurrence');
			ylabel('Neuron ID');
			title('Assembly membership');
			
			
			
			subplot(2,2,4)
			vecLowerBound = mean(matRandMembership,2)-std(matRandMembership,[],2)*dblCorrectedSDs;
			vecUpperBound = mean(matRandMembership,2)+std(matRandMembership,[],2)*dblCorrectedSDs;
			vecData = sum(matAssemblyMembers(:,indAssembly),2)/sum(indAssembly);
			
			indLower = vecData<vecLowerBound;
			vecLower = zeros(size(vecData));
			vecLower(indLower) = vecData(indLower);
			
			indHigher = vecData>vecUpperBound;
			vecHigher = zeros(size(vecData));
			vecHigher(indHigher) = vecData(indHigher);
			
			vecNeither = zeros(size(vecData));
			vecNeither(~indHigher & ~indLower) = vecData(~indHigher & ~indLower);
			
			bar(vecNeither,'k');
			hold on
			bar(vecHigher,'g');
			bar(vecLower,'r');
			errorbar(1:intNeurons,mean(matRandMembership,2),std(matRandMembership,[],2)*dblCorrectedSDs,'xb');
			scatter(find(indLower),0.95*ones(size(find(indLower))),'r*');
			scatter(find(indHigher),0.95*ones(size(find(indHigher))),'g*');
			hold off
			ylim([0 1]);
			xlim([0 intNeurons+1]);
			xlabel('Neuron ID');
			ylabel('Membership probability');
			title('Significant assembly memberships');
			
			%save output
			matAssemblyCoreMembers(intAssembly,:) = indHigher;
			
			if sParams.boolSavePlots
				drawnow;
				jFig = get(handle(gcf), 'JavaFrame');
				jFig.setMaximized(true);
				figure(gcf);
				drawnow;
				export_fig(sprintf('%spop%dexample_assembly%d.tif',strSes,intPopulation,intAssembly));
				export_fig(sprintf('%spop%dexample_assembly%d.pdf',strSes,intPopulation,intAssembly));
			end
			
			%pause
			%sum(matAssemblyMembers(:,indAssembly),2)
			
		end
		
		%% plot assembly summary
		hCheckAssemblies = figure;
		subplot(3,2,1)
		imagesc(abs(matTrialResponse-max(matTrialResponse(:))),[0.1 max(matTrialResponse(:))]);
		colormap('hot');colorbar;freezeColors;cbfreeze;
		title('dF/F0 activity whole recording')
		xlabel('Trial')
		ylabel('Neuron')
		
		subplot(3,2,3)
		imagesc(abs(matSpikeCounts-max(matSpikeCounts(:))),[2 10]);
		colormap('hot');colorbar;freezeColors;cbfreeze;
		title('Activation events whole recording')
		xlabel('Trial')
		ylabel('Neuron')
		
		subplot(3,2,5)
		imagesc(abs(matRespAssemblies-max(matRespAssemblies(:))));
		colormap('hot');colorbar;freezeColors;cbfreeze;
		title('Assembly activity whole recording')
		xlabel('Trial')
		ylabel('Assembly')
		
		
		subplot(2,2,2)
		imagesc(matReorderedAssemblyCorrs,[-1 1]);
		colormap('redblue');colorbar;freezeColors;cbfreeze;
		title('Population events ordered by assembly')
		xlabel('Population Event')
		ylabel('Population Event')
		
		subplot(2,2,4)
		imagesc(matM,[-1 1]);
		colormap('redblue');colorbar;freezeColors;cbfreeze;
		title('Assembly consistency')
		xlabel('Recurring Assembly')
		ylabel('Recurring Assembly')
		
		if sParams.boolSavePlots
			drawnow;
			jFig = get(handle(gcf), 'JavaFrame');
			jFig.setMaximized(true);
			figure(gcf);
			drawnow;
			export_fig(sprintf('%spop%dassemblies_summary.tif',strSes,intPopulation));
			export_fig(sprintf('%spop%dassemblies_summary.pdf',strSes,intPopulation));
		end
		
		%% basic assembly properties
		hAssemblyProperties = figure;
		subplot(2,2,1);
		imagesc(1-matAssemblyCoreMembers');colormap(hot);freezeColors;
		ylabel('Neuron ID');
		xlabel('Assembly ID')
		title('Significant core member');
		
		subplot(2,2,3);
		imagesc(corr(matAssemblyCoreMembers'),[-1 1]);colormap(redblue);freezeColors;colorbar;cbfreeze;
		ylabel('Assembly ID');
		xlabel('Assembly ID')
		title('Core member cross-correlations');
		
		subplot(3,2,2);
		vecCoreMemberSize = sum(matAssemblyCoreMembers,2);
		vecAssemblyAutoCorr = diag(matM);
		indRecurringClusters = zscore(vecAssemblyAutoCorr)>-2;
		[dummy,intNonRecurringCluster]=min(zscore(vecAssemblyAutoCorr));
		indRecurringClusters(intNonRecurringCluster) = false;
		indMultiNeuronClusters = vecCoreMemberSize>1;
		indRealAssemblies = indRecurringClusters & indMultiNeuronClusters;
		bar(find(indMultiNeuronClusters),vecCoreMemberSize(indMultiNeuronClusters),'g');
		hold on
		bar(find(~indMultiNeuronClusters),vecCoreMemberSize(~indMultiNeuronClusters),'r');
		hold off
		title('Multi-neuron clusters');
		xlabel('Assembly ID');
		ylabel('Number of neuronal core members');
		
		subplot(3,2,4)
		bar(find(indRecurringClusters),vecAssemblyAutoCorr(indRecurringClusters),'g');
		hold on
		bar(find(~indRecurringClusters),vecAssemblyAutoCorr(~indRecurringClusters),'r');
		hold off
		title('Consistently recurring clusters');
		xlabel('Assembly ID');
		ylabel('Cluster consistency (r)');
		
		subplot(3,2,6)
		bar(indRealAssemblies,'g');
		title('Real Assemblies');
		xlabel('Assembly ID');
		
		if sParams.boolSavePlots
			drawnow;
			jFig = get(handle(gcf), 'JavaFrame');
			jFig.setMaximized(true);
			figure(gcf);
			drawnow;
			export_fig(sprintf('%spop%dassembly_properties.tif',strSes,intPopulation));
			export_fig(sprintf('%spop%dassembly_properties.pdf',strSes,intPopulation));
		end
		
		%% assembly correlates
		%decoding; remove 0% contrast trials
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
		

		% decode orientation with template matching
		[dblPerformance_dFoF,vecDecodedIndexCV,matDists] = doCrossValidatedTemplateMatching(matTrialResponse(:,indStimTrials),vecOriStimTrials);
		[dblPerformance_AEs,vecDecodedIndexCV,matDists] = doCrossValidatedTemplateMatching(matSpikeCounts(:,indStimTrials),vecOriStimTrials);
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
		subplot(2,2,1)
		vecDecodingPerformance = [dblPerformance_dFoF dblPerformance_AEs];
		errorbar([0.3 0.7],vecDecodingPerformance,[dblCI_dFoF(1) dblCI_AEs(1)]-vecDecodingPerformance,[dblCI_dFoF(2) dblCI_AEs(2)]-vecDecodingPerformance,'xb')
		xlim([0 1])
		ylim([0 1])
		hold on
		plot(get(gca,'xlim'),[1 1]/length(cellSelectSubO),'k--')
		set(gca,'XTick',[0.3 0.7],'xticklabel',{'Neuron dF/F0','Neuron AE'})
		ylabel('Orientation decoding accuracy (TM)')
		
		%% get contrast repsonse functions
		vecContrasts = [0 0.5 2 8 32 100];
		cellHit_dFoF = cell(1,intContrasts);
		cellMiss_dFoF = cell(1,intContrasts);
		cellHit_AE = cell(1,intContrasts);
		cellMiss_AE = cell(1,intContrasts);
		cellHit_PE = cell(1,intContrasts);
		cellMiss_PE = cell(1,intContrasts);
		
		%make normalized activity matrix per contrast
		for intContrast=1:intContrasts
			cellHit_dFoF{intContrast} = nan(1,sum(cellSelectContrasts{intContrast}));
			cellMiss_dFoF{intContrast} = nan(1,sum(cellSelectContrasts{intContrast}));
			cellHit_AE{intContrast} = nan(1,sum(cellSelectContrasts{intContrast}));
			cellMiss_AE{intContrast} = nan(1,sum(cellSelectContrasts{intContrast}));
			cellHit_PE{intContrast} = nan(1,sum(cellSelectContrasts{intContrast}));
			cellMiss_PE{intContrast} = nan(1,sum(cellSelectContrasts{intContrast}));
		end
		
		%loop through trials
		for intTrial = 1:intTrials
			intContrast = vecStimContrasts(intTrial);
			vec_dFoF = matTrialResponse(:,intTrial);
			vec_AE = matSpikeCounts(:,intTrial);
			dblPE = sum(matRespAssemblies(:,intTrial),1);
			if indStimResp(intTrial)
				cellHit_dFoF{intContrast}(find(isnan(cellHit_dFoF{intContrast}),1)) = mean(vec_dFoF);
				cellHit_AE{intContrast}(find(isnan(cellHit_AE{intContrast}),1)) = mean(vec_AE);
				cellHit_PE{intContrast}(find(isnan(cellHit_PE{intContrast}),1)) = dblPE;
			else
				cellMiss_dFoF{intContrast}(find(isnan(cellMiss_dFoF{intContrast}),1)) = mean(vec_dFoF);
				cellMiss_AE{intContrast}(find(isnan(cellMiss_AE{intContrast}),1)) = mean(vec_AE);
				cellMiss_PE{intContrast}(find(isnan(cellMiss_PE{intContrast}),1)) = dblPE;
			end
		end
		
		%remove trailing nans
		for intContrast=1:length(cellHit_dFoF)
			cellHit_dFoF{intContrast} = cellHit_dFoF{intContrast}(~isnan(cellHit_dFoF{intContrast}));
			cellMiss_dFoF{intContrast} = cellMiss_dFoF{intContrast}(~isnan(cellMiss_dFoF{intContrast}));
			cellHit_AE{intContrast} = cellHit_AE{intContrast}(~isnan(cellHit_AE{intContrast}));
			cellMiss_AE{intContrast} = cellMiss_AE{intContrast}(~isnan(cellMiss_AE{intContrast}));
			cellHit_PE{intContrast} = cellHit_PE{intContrast}(~isnan(cellHit_PE{intContrast}));
			cellMiss_PE{intContrast} = cellMiss_PE{intContrast}(~isnan(cellMiss_PE{intContrast}));
		end
		vecMeanHit_dFoF = cellfun(@mean,cellHit_dFoF);
		vecSEHit_dFoF = cellfun(@std,cellHit_dFoF)./cellfun(@numel,cellHit_dFoF);
		vecMeanMiss_dFoF = cellfun(@mean,cellMiss_dFoF);
		vecSEMiss_dFoF = cellfun(@std,cellMiss_dFoF)./cellfun(@numel,cellMiss_dFoF);
		
		vecMeanHit_AE = cellfun(@mean,cellHit_AE);
		vecSEHit_AE = cellfun(@std,cellHit_AE)./cellfun(@numel,cellHit_AE);
		vecMeanMiss_AE = cellfun(@mean,cellMiss_AE);
		vecSEMiss_AE = cellfun(@std,cellMiss_AE)./cellfun(@numel,cellMiss_AE);
		
		vecMeanHit_PE = cellfun(@mean,cellHit_PE);
		vecSEHit_PE = cellfun(@std,cellHit_PE)./cellfun(@numel,cellHit_PE);
		vecMeanMiss_PE = cellfun(@mean,cellMiss_PE);
		vecSEMiss_PE = cellfun(@std,cellMiss_PE)./cellfun(@numel,cellMiss_PE);
		
		%plot sigmoid dF/F
		subplot(2,2,2)
		vecPlotC = vecContrasts;
		vecPlotC(1) = 0.2;
		hold on
		errorfill(vecPlotC,vecMeanMiss_dFoF,vecSEMiss_dFoF,[1 0 0],[1 0.5 0.5]);
		errorfill(vecPlotC,vecMeanHit_dFoF,vecSEHit_dFoF,[0 1 0],[0.5 1 0.5]);
		hold off
		set(gca,'XScale','log')
		set(gca,'Xtick',vecPlotC,'xticklabel',vecContrasts)
		ylabel('Mean dF/F0')
		xlabel('Stimulus contrast (%)')
		
		%plot sigmoid AEs
		subplot(2,2,3)
		vecPlotC = vecContrasts;
		vecPlotC(1) = 0.2;
		hold on
		errorfill(vecPlotC,vecMeanMiss_AE,vecSEMiss_AE,[1 0 0],[1 0.5 0.5]);
		errorfill(vecPlotC,vecMeanHit_AE,vecSEHit_AE,[0 1 0],[0.5 1 0.5]);
		hold off
		set(gca,'XScale','log')
		set(gca,'Xtick',vecPlotC,'xticklabel',vecContrasts)
		ylabel('Mean AEs')
		xlabel('Stimulus contrast (%)')
		
		%plot sigmoid PEs
		subplot(2,2,4)
		vecPlotC = vecContrasts;
		vecPlotC(1) = 0.2;
		hold on
		errorfill(vecPlotC,vecMeanMiss_PE,vecSEMiss_PE,[1 0 0],[1 0.5 0.5]);
		errorfill(vecPlotC,vecMeanHit_PE,vecSEHit_PE,[0 1 0],[0.5 1 0.5]);
		hold off
		set(gca,'XScale','log')
		set(gca,'Xtick',vecPlotC,'xticklabel',vecContrasts)
		ylabel('Mean PEs')
		xlabel('Stimulus contrast (%)')
		
		
		% save fig
		if sParams.boolSavePlots
			drawnow;
			jFig = get(handle(gcf), 'JavaFrame');
			jFig.setMaximized(true);
			figure(gcf);
			drawnow;
			export_fig(sprintf('%spop%dcontrast_responses.tif',strSes,intPopulation));
			export_fig(sprintf('%spop%dcontrast_responses.pdf',strSes,intPopulation));
		end
		
		%save data
		cellSaveContAnal{1,intPopulation} = vecDecodingPerformance; %dFoF AE Ass
		
		cellSaveContAnal{2,intPopulation} = [vecMeanHit_dFoF;vecMeanMiss_dFoF]; 
		cellSaveContAnal{3,intPopulation} = [vecMeanHit_AE;vecMeanMiss_AE]; 
		cellSaveContAnal{4,intPopulation} = [vecMeanHit_PE;vecMeanMiss_PE]; 
		
		%% make sigmoid plots; presence, consistency of assemblies for hit/miss as function of contrast, baseline occurrence preceding hit, miss, slow and fast, and non-recurring cluster
		hNonRecAssemblyCorrelates = figure;
		hRealAssemblyCorrelates = figure;
		hNonRecAssemblyCorrelatesSlowFast = figure;
		hRealAssemblyCorrelatesSlowFast = figure;
		
		%% presence
		vecContrasts = [0 0.5 2 8 32 100];
		cellHitRealAssAct = cell(1,intContrasts);
		cellMissRealAssAct = cell(1,intContrasts);
		cellHitNRAssAct = cell(1,intContrasts);
		cellMissNRAssAct = cell(1,intContrasts);
		cellFastRealAssAct = cell(1,intContrasts);
		cellSlowRealAssAct = cell(1,intContrasts);
		cellFastNRAssAct = cell(1,intContrasts);
		cellSlowNRAssAct = cell(1,intContrasts);
		
		%make normalized activity matrix per contrast
		for intContrast=1:intContrasts
			cellHitRealAssAct{intContrast} = nan(1,sum(cellSelectContrasts{intContrast}));
			cellMissRealAssAct{intContrast} = nan(1,sum(cellSelectContrasts{intContrast}));
			cellHitNRAssAct{intContrast} = nan(1,sum(cellSelectContrasts{intContrast}));
			cellMissNRAssAct{intContrast} = nan(1,sum(cellSelectContrasts{intContrast}));
			cellFastRealAssAct{intContrast} = nan(1,sum(cellSelectContrasts{intContrast}));
			cellSlowRealAssAct{intContrast} = nan(1,sum(cellSelectContrasts{intContrast}));
			cellFastNRAssAct{intContrast} = nan(1,sum(cellSelectContrasts{intContrast}));
			cellSlowNRAssAct{intContrast} = nan(1,sum(cellSelectContrasts{intContrast}));
		end
		
		%loop through trials
		for intTrial = 1:intTrials
			intContrast = vecStimContrasts(intTrial);
			vecActReal = matRespRealAssemblies(:,intTrial);
			vecActNR = matRespAssemblies(~indRealAssemblies,intTrial);
			if indStimResp(intTrial)
				cellHitRealAssAct{intContrast}(find(isnan(cellHitRealAssAct{intContrast}),1)) = mean(vecActReal);
				cellHitNRAssAct{intContrast}(find(isnan(cellHitNRAssAct{intContrast}),1)) = mean(vecActNR);
				%fast/slow?
				if indFast(intTrial)
					cellFastRealAssAct{intContrast}(find(isnan(cellHitRealAssAct{intContrast}),1)) = mean(vecActReal);
					cellFastNRAssAct{intContrast}(find(isnan(cellHitNRAssAct{intContrast}),1)) = mean(vecActNR);
				elseif indSlow(intTrial)
					cellSlowRealAssAct{intContrast}(find(isnan(cellHitRealAssAct{intContrast}),1)) = mean(vecActReal);
					cellSlowNRAssAct{intContrast}(find(isnan(cellHitNRAssAct{intContrast}),1)) = mean(vecActNR);
				end
			else
				cellMissRealAssAct{intContrast}(find(isnan(cellMissRealAssAct{intContrast}),1)) = mean(vecActReal);
				cellMissNRAssAct{intContrast}(find(isnan(cellMissNRAssAct{intContrast}),1)) = mean(vecActNR);
			end
		end
		
		%remove trailing nans
		for intContrast=1:length(cellHitRealAssAct)
			cellHitRealAssAct{intContrast} = cellHitRealAssAct{intContrast}(~isnan(cellHitRealAssAct{intContrast}));
			cellMissRealAssAct{intContrast} = cellMissRealAssAct{intContrast}(~isnan(cellMissRealAssAct{intContrast}));
			cellHitNRAssAct{intContrast} = cellHitNRAssAct{intContrast}(~isnan(cellHitNRAssAct{intContrast}));
			cellMissNRAssAct{intContrast} = cellMissNRAssAct{intContrast}(~isnan(cellMissNRAssAct{intContrast}));
			cellFastRealAssAct{intContrast} = cellFastRealAssAct{intContrast}(~isnan(cellFastRealAssAct{intContrast}));
			cellSlowRealAssAct{intContrast} = cellSlowRealAssAct{intContrast}(~isnan(cellSlowRealAssAct{intContrast}));
			cellFastNRAssAct{intContrast} = cellFastNRAssAct{intContrast}(~isnan(cellFastNRAssAct{intContrast}));
			cellSlowNRAssAct{intContrast} = cellSlowNRAssAct{intContrast}(~isnan(cellSlowNRAssAct{intContrast}));
		end
		vecMeanHit = cellfun(@mean,cellHitRealAssAct);
		vecSEHit = cellfun(@std,cellHitRealAssAct)./cellfun(@numel,cellHitRealAssAct);
		vecMeanMiss = cellfun(@mean,cellMissRealAssAct);
		vecSEMiss = cellfun(@std,cellMissRealAssAct)./cellfun(@numel,cellMissRealAssAct);
		
		vecMeanHitNR = cellfun(@mean,cellHitNRAssAct);
		vecSEHitNR = cellfun(@std,cellHitNRAssAct)./cellfun(@numel,cellHitNRAssAct);
		vecMeanMissNR = cellfun(@mean,cellMissNRAssAct);
		vecSEMissNR = cellfun(@std,cellMissNRAssAct)./cellfun(@numel,cellMissNRAssAct);
		
		vecMeanSlow = cellfun(@mean,cellSlowRealAssAct);
		vecSESlow = cellfun(@std,cellSlowRealAssAct)./cellfun(@numel,cellSlowRealAssAct);
		vecMeanFast = cellfun(@mean,cellFastRealAssAct);
		vecSEFast = cellfun(@std,cellFastRealAssAct)./cellfun(@numel,cellFastRealAssAct);
		
		vecMeanSlowNR = cellfun(@mean,cellSlowNRAssAct);
		vecSESlowNR = cellfun(@std,cellSlowNRAssAct)./cellfun(@numel,cellSlowNRAssAct);
		vecMeanFastNR = cellfun(@mean,cellFastNRAssAct);
		vecSEFastNR = cellfun(@std,cellFastNRAssAct)./cellfun(@numel,cellFastNRAssAct);
		
		%get pre-trial assembly activity
		vecMeanPreHitRealAssAct = mean(matPreRespRealAssemblies(:,indStimResp),1);
		vecMeanPreMissRealAssAct = mean(matPreRespRealAssemblies(:,~indStimResp),1);
		vecMeanPreHitNRAssAct = mean(matPreRespNRAssemblies(:,indStimResp),1);
		vecMeanPreMissNRAssAct = mean(matPreRespNRAssemblies(:,~indStimResp),1);
		vecMeanPreFastRealAssAct = mean(matPreRespRealAssemblies(:,indFast),1);
		vecMeanPreSlowRealAssAct = mean(matPreRespRealAssemblies(:,indSlow),1);
		vecMeanPreFastNRAssAct = mean(matPreRespNRAssemblies(:,indFast),1);
		vecMeanPreSlowNRAssAct = mean(matPreRespNRAssemblies(:,indSlow),1);
		
		
		%plot recurring assemblies
		figure(hRealAssemblyCorrelates)
		subplot(2,3,2)
		vecPlotC = vecContrasts;
		vecPlotC(1) = 0.2;
		errorbar([vecPlotC(1) vecPlotC(end)],mean(vecMeanPreHitRealAssAct)*[1 1],std(vecMeanPreHitRealAssAct)*[1 1]/sqrt(sum(indStimResp)),'g')
		hold on
		errorbar([vecPlotC(1) vecPlotC(end)],mean(vecMeanPreMissRealAssAct)*[1 1],std(vecMeanPreMissRealAssAct)*[1 1]/sqrt(sum(~indStimResp)),'r')
		errorfill(vecPlotC,vecMeanMiss,vecSEMiss,[1 0 0],[1 0.5 0.5]);
		errorfill(vecPlotC,vecMeanHit,vecSEHit,[0 1 0],[0.5 1 0.5]);
		hold off
		set(gca,'XScale','log')
		set(gca,'Xtick',vecPlotC,'xticklabel',vecContrasts)
		ylabel('Mean assembly presence')
		xlabel('Stimulus contrast (%)')
		title('Recurring')
		
		%plot non-recurring assemblies
		figure(hNonRecAssemblyCorrelates)
		subplot(2,3,2)
		errorbar([vecPlotC(1) vecPlotC(end)],mean(vecMeanPreHitNRAssAct)*[1 1],std(vecMeanPreHitNRAssAct)*[1 1]/sqrt(sum(indStimResp)),'g')
		hold on
		errorbar([vecPlotC(1) vecPlotC(end)],mean(vecMeanPreMissNRAssAct)*[1 1],std(vecMeanPreMissNRAssAct)*[1 1]/sqrt(sum(~indStimResp)),'r')
		errorfill(vecPlotC,vecMeanMissNR,vecSEMissNR,[1 0 0],[1 0.5 0.5]);
		errorfill(vecPlotC,vecMeanHitNR,vecSEHitNR,[0 1 0],[0.5 1 0.5]);
		hold off
		set(gca,'XScale','log')
		set(gca,'Xtick',vecPlotC,'xticklabel',vecContrasts)
		ylabel('Mean assembly presence')
		xlabel('Stimulus contrast (%)')
		title('Non-recurring')
		
		%plot recurring assemblies
		figure(hRealAssemblyCorrelatesSlowFast)
		subplot(2,3,2)
		errorbar([vecPlotC(1) vecPlotC(end)],mean(vecMeanPreFastRealAssAct)*[1 1],std(vecMeanPreFastRealAssAct)*[1 1]/sqrt(sum(indStimResp)),'g')
		hold on
		errorbar([vecPlotC(1) vecPlotC(end)],mean(vecMeanPreSlowRealAssAct)*[1 1],std(vecMeanPreSlowRealAssAct)*[1 1]/sqrt(sum(~indStimResp)),'m')
		errorfill(vecPlotC,vecMeanSlow,vecSESlow,[1 0 0],[1 0.5 0.5]);
		errorfill(vecPlotC,vecMeanFast,vecSEFast,[0 1 0],[0.5 1 0.5]);
		hold off
		set(gca,'XScale','log')
		set(gca,'Xtick',vecPlotC,'xticklabel',vecContrasts)
		ylabel('Mean assembly presence')
		xlabel('Stimulus contrast (%)')
		title('Recurring')
		
		%plot non-recurring assemblies
		figure(hNonRecAssemblyCorrelatesSlowFast)
		subplot(2,3,2)
		errorbar([vecPlotC(1) vecPlotC(end)],mean(vecMeanPreFastNRAssAct)*[1 1],std(vecMeanPreFastNRAssAct)*[1 1]/sqrt(sum(indStimResp)),'g')
		hold on
		errorbar([vecPlotC(1) vecPlotC(end)],mean(vecMeanPreSlowNRAssAct)*[1 1],std(vecMeanPreSlowNRAssAct)*[1 1]/sqrt(sum(~indStimResp)),'m')
		errorfill(vecPlotC,vecMeanSlowNR,vecSESlowNR,[0.5 0 0.5],[1 0.5 1]);
		errorfill(vecPlotC,vecMeanFastNR,vecSEFastNR,[0 1 0],[0.5 1 0.5]);
		hold off
		set(gca,'XScale','log')
		set(gca,'Xtick',vecPlotC,'xticklabel',vecContrasts)
		ylabel('Mean assembly presence')
		xlabel('Stimulus contrast (%)')
		title('Non-recurring')
		
		%% occurrence, size, number of AEs per occurrence, duration per occurrence
		%get all frames during stimulus period
		vecAssemblyStarts = sAssemblies{intPopulation}.vecAssemblyStarts;
		vecAssemblyStops = sAssemblies{intPopulation}.vecAssemblyStops;
		vecAssemblyDuration = (vecAssemblyStops-vecAssemblyStarts)/dblSamplingFreq;
		vecAssemblies = sAssemblies{intPopulation}.vecAssemblies';
		vecRealAssemblies = ismember(vecAssemblies,find(indRealAssemblies));
		vecNonRecurringAssemblies = ismember(vecAssemblies,find(~indRealAssemblies));
		vecAssemblyAEs = zeros(size(vecAssemblyStarts));
		for intPE=1:length(vecAssemblyStarts)
			intStart = vecAssemblyStarts(intPE);
			intStop = vecAssemblyStops(intPE);
			vecAssemblyAEs(intPE) = sum(sum(matSpikeCountsHD(:,intStart:intStop)));
		end
		
		%pre-allocate selection index vectors
		cellRespTypeNames = {'Hit','Miss','Fast','Slow'};
		cellRecurTypeNames = {'Real','NonRec'};
		cellPrePostTypeNames = {'','Pre'};
		for intRespType=1:4
			for intRecurType=1:2
				for intPrePostType=1:2
					eval([sprintf('indOccurrences%s%s%s',cellPrePostTypeNames{intPrePostType},cellRespTypeNames{intRespType},cellRecurTypeNames{intRecurType}) ' = false(size(vecAssemblies));']);
				end
			end
		end
		
		%build selection index vectors
		matAssemblyOccurrences = nan(intTrials,4);
		for intTrial=1:intTrials
			%get selection occurrences
			intOn = sStimTemp.FrameOn(intTrial);
			intOff = sStimTemp.FrameOff(intTrial);
			indTheseRealAssemblies = (vecAssemblyStops > intOn & vecAssemblyStarts < intOff) & vecRealAssemblies;
			indTheseRealPreAssemblies = (vecAssemblyStops > intOn-round(dblStimSecs*dblSamplingFreq) & vecAssemblyStarts < intOn) & vecRealAssemblies;
			indTheseNRAssemblies = (vecAssemblyStops > intOn & vecAssemblyStarts < intOff) & vecNonRecurringAssemblies;
			indTheseNRPreAssemblies = (vecAssemblyStops > intOn-round(dblStimSecs*dblSamplingFreq) & vecAssemblyStarts < intOn) & vecNonRecurringAssemblies;
			
			%assign data
			matAssemblyOccurrences(intTrial,1) = sum(indTheseRealAssemblies); %recur during
			matAssemblyOccurrences(intTrial,2) = sum(indTheseRealPreAssemblies); %recur pre
			matAssemblyOccurrences(intTrial,3) = sum(indTheseNRAssemblies); %non-recur during
			matAssemblyOccurrences(intTrial,4) = sum(indTheseNRPreAssemblies); %non-recur pre
			
			%assign assembly occurrences to selection vectors
			if indHit(intTrial)
				indOccurrencesHitReal(indTheseRealAssemblies) = true; %#ok<SAGROW>
				indOccurrencesPreHitReal(indTheseRealPreAssemblies) = true;%#ok<SAGROW>
				indOccurrencesHitNonRec(indTheseNRAssemblies) = true;%#ok<SAGROW>
				indOccurrencesPreHitNonRec(indTheseNRPreAssemblies) = true;%#ok<SAGROW>
				if indSlow(intTrial)
					indOccurrencesSlowReal(indTheseRealAssemblies) = true; %#ok<SAGROW>
					indOccurrencesPreSlowReal(indTheseRealPreAssemblies) = true;%#ok<SAGROW>
					indOccurrencesSlowNonRec(indTheseNRAssemblies) = true;%#ok<SAGROW>
					indOccurrencesPreSlowNonRec(indTheseNRPreAssemblies) = true;%#ok<SAGROW>
				elseif indFast(intTrial)
					indOccurrencesFastReal(indTheseRealAssemblies) = true; %#ok<SAGROW>
					indOccurrencesPreFastReal(indTheseRealPreAssemblies) = true;%#ok<SAGROW>
					indOccurrencesFastNonRec(indTheseNRAssemblies) = true;%#ok<SAGROW>
					indOccurrencesPreFastNonRec(indTheseNRPreAssemblies) = true;%#ok<SAGROW>
				else
					intTrial
				end
			elseif indMiss(intTrial)
				indOccurrencesMissReal(indTheseRealAssemblies) = true; %#ok<SAGROW>
				indOccurrencesPreMissReal(indTheseRealPreAssemblies) = true;%#ok<SAGROW>
				indOccurrencesMissNonRec(indTheseNRAssemblies) = true;%#ok<SAGROW>
				indOccurrencesPreMissNonRec(indTheseNRPreAssemblies) = true;%#ok<SAGROW>
			else
				intTrial
			end
		end
		
		%{
		%% SANITY CHECK: plot occurrences
		
		%plot to check
		figure
		vecOccHitMean = nan(1,6);
		vecOccHitErr = nan(1,6);
		vecOccMissMean = nan(1,6);
		vecOccMissErr = nan(1,6);
		vecOccHitNum = nan(1,6);
		vecOccTotNum = nan(1,6);
		for intC=1:6
			vecOccC_Hit = matAssemblyOccurrences(cellSelectContrasts{intC}&indHit,1);
			vecOccC_Miss =  matAssemblyOccurrences(cellSelectContrasts{intC}&indMiss,1);
			vecOccHitMean(intC) = mean(vecOccC_Hit);
			vecOccHitErr(intC) = std(vecOccC_Hit)/sqrt(numel(vecOccC_Hit));
			vecOccMissMean(intC) = mean(vecOccC_Miss);
			vecOccMissErr(intC) = std(vecOccC_Miss)/sqrt(numel(vecOccC_Miss));
			vecOccHitNum(intC) = sum(vecOccC_Hit);
			vecOccTotNum(intC) = sum(vecOccC_Hit) + sum(vecOccC_Miss);
		end
		errorbar(1:6,vecOccHitMean,vecOccHitErr,'gx')
		hold on
		errorbar(1:6,vecOccMissMean,vecOccMissErr,'rx')
		hold off
		%}
		
		%% pre-allocate output variables
		for intRespType=1:4
			for intRecurType=1:2
				for intPrePostType=1:2
					eval([sprintf('cellAssembliesPerContrast%s%s%s',cellPrePostTypeNames{intPrePostType},...
						cellRespTypeNames{intRespType},cellRecurTypeNames{intRecurType}) ' = cell(1,intContrasts);']);
					eval([sprintf('cellAssemblyMembers%s%s%s',cellPrePostTypeNames{intPrePostType},...
						cellRespTypeNames{intRespType},cellRecurTypeNames{intRecurType}) ' = cell(intAssemblies,intContrasts);']);
					eval([sprintf('cellAssemblyDurPerContrast%s%s%s',cellPrePostTypeNames{intPrePostType},...
						cellRespTypeNames{intRespType},cellRecurTypeNames{intRecurType}) ' = cell(1,intContrasts);']);
					eval([sprintf('cellAssemblyAEsPerContrast%s%s%s',cellPrePostTypeNames{intPrePostType},...
						cellRespTypeNames{intRespType},cellRecurTypeNames{intRecurType}) ' = cell(1,intContrasts);']);
				end
			end
		end
		
		% get occurrences, duration & AEs per contrast
		for intContrast=1:intContrasts
			indSelectTrialsC = cellSelectContrasts{intContrast};
			indOccurrencesThisContrastDuring = false(size(vecAssemblies));
			indOccurrencesThisContrastPreceding = false(size(vecAssemblies));
			
			for intTrial=find(indSelectTrialsC)
				intOn = sStimTemp.FrameOn(intTrial);
				intOff = sStimTemp.FrameOff(intTrial);
				indOccurrencesThisContrastPreceding(vecAssemblyStops > (intOn-round(dblStimSecs*dblSamplingFreq)) & vecAssemblyStarts < (intOff-round(dblStimSecs*dblSamplingFreq))) = true;
				indOccurrencesThisContrastDuring(vecAssemblyStops > intOn & vecAssemblyStarts < intOff) = true;
			end
			
			for intPrePostType=1:2
				if intPrePostType==1
					indOccurrencesThisContrast = indOccurrencesThisContrastDuring;
				else
					indOccurrencesThisContrast = indOccurrencesThisContrastPreceding;
				end
				for intRespType=1:4
					for intRecurType=1:2
						
						%get selection vector for this type
						strSelVecName = sprintf('indOccurrences%s%s%s',cellPrePostTypeNames{intPrePostType},cellRespTypeNames{intRespType},cellRecurTypeNames{intRecurType});
						
						%get output names for this type
						strVarName1 = sprintf('cellAssembliesPerContrast%s%s%s{intContrast}',cellPrePostTypeNames{intPrePostType},...
							cellRespTypeNames{intRespType},cellRecurTypeNames{intRecurType});
						strVarName2 = sprintf('cellAssemblyDurPerContrast%s%s%s{intContrast}',cellPrePostTypeNames{intPrePostType},...
							cellRespTypeNames{intRespType},cellRecurTypeNames{intRecurType});
						strVarName3 = sprintf('cellAssemblyAEsPerContrast%s%s%s{intContrast}',cellPrePostTypeNames{intPrePostType},...
							cellRespTypeNames{intRespType},cellRecurTypeNames{intRecurType});
						
						%execute commands
						%eval(['sum(indOccurrencesThisContrast & ' strSelVecName ')'])
						eval([strVarName1 ' = [' strVarName1 ' matAssemblyMembers(:,indOccurrencesThisContrast & ' strSelVecName ')];']);
						eval([strVarName2 ' = [' strVarName2 ' vecAssemblyDuration(indOccurrencesThisContrast & ' strSelVecName ')];']);
						eval([strVarName3 ' = [' strVarName3 ' vecAssemblyAEs(indOccurrencesThisContrast & ' strSelVecName ')];']);
						
						%get member name
						strVarName4 = sprintf('cellAssemblyMembers%s%s%s{intAssembly,intContrast}',cellPrePostTypeNames{intPrePostType},...
							cellRespTypeNames{intRespType},cellRecurTypeNames{intRecurType});
					
						for intAssembly=1:intAssemblies
							%get specific assembly occurrences
							indSelectOccurrences = vecAssemblies == intAssembly;
							
							%execute commands
							eval([strVarName4 ' = matAssemblyMembers(:,indSelectOccurrences & indOccurrencesThisContrast & ' strSelVecName ');']);
						end
					end
				end
			end
		end
		
		%% remove inconsistent allocations of recurring/non-recurring assemblies
		for intRespType=1:4
			for intRecurType=1:2
				for intPrePostType=1:2
					%get variable name
					strVarName = sprintf('cellAssemblyMembers%s%s%s',cellPrePostTypeNames{intPrePostType},...
							cellRespTypeNames{intRespType},cellRecurTypeNames{intRecurType});
					
					%create selection string
					if intRecurType == 1
						strSelect = '(~indRealAssemblies,:)';
					else
						strSelect = '(indRealAssemblies,:)';
					end
					
					%execute
					eval([strVarName strSelect ' = [];']);
				end
			end
		end

		%% create plotting variables and perform plotting
		for intRespPlotType=1:2
			for intRecurType=1:2
				%shorten name
				strRecur = cellRecurTypeNames{intRecurType};
				
				%get target figure
				if intRespPlotType == 1 %hit/miss
					strGood = cellRespTypeNames{1}; %hit
					strBad = cellRespTypeNames{2}; %miss
					if intRecurType == 1 %recurring
						figure(hRealAssemblyCorrelates);
					elseif intRecurType == 2 %non-recurring
						figure(hNonRecAssemblyCorrelates);
					end
				elseif intRespPlotType == 2 %fast/slow
					strGood = cellRespTypeNames{3}; %fast
					strBad = cellRespTypeNames{4}; %slow
					if intRecurType == 1 %recurring
						figure(hRealAssemblyCorrelatesSlowFast);
					elseif intRecurType == 2 %non-recurring
						figure(hNonRecAssemblyCorrelatesSlowFast);
					end
				end
				
				%get target variable names
				%members
				strNameMembersDuringGood = sprintf('cellAssemblyMembers%s%s',strGood,strRecur);
				strNameMembersDuringBad = sprintf('cellAssemblyMembers%s%s',strBad,strRecur);
				strNameMembersPreGood = sprintf('cellAssemblyMembersPre%s%s',strGood,strRecur);
				strNameMembersPreBad = sprintf('cellAssemblyMembersPre%s%s',strBad,strRecur);
				
				%assemblies per contrast
				strNameAPCDuringGood = sprintf('cellAssembliesPerContrast%s%s',strGood,strRecur);
				strNameAPCDuringBad = sprintf('cellAssembliesPerContrast%s%s',strBad,strRecur);
				strNameAPCPreGood = sprintf('cellAssembliesPerContrastPre%s%s',strGood,strRecur);
				strNameAPCPreBad = sprintf('cellAssembliesPerContrastPre%s%s',strBad,strRecur);
				
				%duration per contrast
				strNameDPCDuringGood = sprintf('cellAssemblyDurPerContrast%s%s',strGood,strRecur);
				strNameDPCDuringBad = sprintf('cellAssemblyDurPerContrast%s%s',strBad,strRecur);
				strNameDPCPreGood = sprintf('cellAssemblyDurPerContrastPre%s%s',strGood,strRecur);
				strNameDPCPreBad = sprintf('cellAssemblyDurPerContrastPre%s%s',strBad,strRecur);

				%AEs per contrast
				strNameAEsPCDuringGood = sprintf('cellAssemblyAEsPerContrast%s%s',strGood,strRecur);
				strNameAEsPCDuringBad = sprintf('cellAssemblyAEsPerContrast%s%s',strBad,strRecur);
				strNameAEsPCPreGood = sprintf('cellAssemblyAEsPerContrastPre%s%s',strGood,strRecur);
				strNameAEsPCPreBad = sprintf('cellAssemblyAEsPerContrastPre%s%s',strBad,strRecur);

				%mean size pre-trial
				cellSizePreGood = cellfun(@sum,eval(strNameMembersPreGood),'UniformOutput',false);
				vecPreGoodSizes = cell2mat(cellSizePreGood(:)');
				dblMeanSizePreGood = mean(vecPreGoodSizes);
				dblSESizePreGood = std(vecPreGoodSizes)/sqrt(length(vecPreGoodSizes));
				
				cellSizePreBad = cellfun(@sum,eval(strNameMembersPreBad),'UniformOutput',false);
				vecPreBadSizes = cell2mat(cellSizePreBad(:)');
				dblMeanSizePreBad = mean(vecPreBadSizes);
				dblSESizePreBad = std(vecPreBadSizes)/sqrt(length(vecPreBadSizes));
				
				%calculate size
				matMeanSizePerAssembly_Good = cellfun(@mean,cellfun(@sum,eval(strNameMembersDuringGood),'UniformOutput',false));
				matMeanSizePerAssembly_Bad = cellfun(@mean,cellfun(@sum,eval(strNameMembersDuringBad),'UniformOutput',false));
				matMeanSizePerAssembly_PreGood = cellfun(@mean,cellfun(@sum,eval(strNameMembersPreGood),'UniformOutput',false));
				matMeanSizePerAssembly_PreBad = cellfun(@mean,cellfun(@sum,eval(strNameMembersPreBad),'UniformOutput',false));
				
				%mean number of assembly occurrences per trial per contrast
				intTrialsPerContrast = sum(vecStimContrasts==1);
				vecTrialOccurrencesGood = cellfun(@size,eval(strNameAPCDuringGood),cellfill(2,size(eval(strNameAPCDuringGood))));
				vecTrialOccurrencesBad = cellfun(@size,eval(strNameAPCDuringBad),cellfill(2,size(eval(strNameAPCDuringBad))));
				vecTrialOccurrencesPreGood = cellfun(@size,eval(strNameAPCPreGood),cellfill(2,size(eval(strNameAPCPreGood))));
				vecTrialOccurrencesPreBad = cellfun(@size,eval(strNameAPCPreBad),cellfill(2,size(eval(strNameAPCPreBad))));
				
				%mean size of assemblies per trial per contrast
				cellTrialSizeGood = cellfun(@sum,eval(strNameAPCDuringGood),'UniformOutput',false);
				vecTrialMeanSizeGood = cellfun(@mean,cellTrialSizeGood);
				vecTrialSESizeGood = cellfun(@std,cellTrialSizeGood)./sqrt(vecTrialOccurrencesGood);
				cellTrialSizeBad = cellfun(@sum,eval(strNameAPCDuringBad),'UniformOutput',false);
				vecTrialMeanSizeBad = cellfun(@mean,cellTrialSizeBad);
				vecTrialSESizeBad = cellfun(@std,cellTrialSizeBad)./sqrt(vecTrialOccurrencesBad);
				
				%get duration
				vecMeanDurBad = cellfun(@mean,eval(strNameDPCDuringBad));
				vecSEDurBad = cellfun(@std,eval(strNameDPCDuringBad))./sqrt(cellfun(@numel,eval(strNameDPCDuringBad)));

				vecMeanDurGood = cellfun(@mean,eval(strNameDPCDuringGood));
				vecSEDurGood = cellfun(@std,eval(strNameDPCDuringGood))./sqrt(cellfun(@numel,eval(strNameDPCDuringGood)));
		
				dblMeanPreDurBad = mean(cell2mat(eval(strNameDPCPreBad)));
				dblSEPreDurBad = std(cell2mat(eval(strNameDPCPreBad)))/numel(cell2mat(eval(strNameDPCPreBad)));
				
				dblMeanPreDurGood = mean(cell2mat(eval(strNameDPCPreGood)));
				dblSEPreDurGood = std(cell2mat(eval(strNameDPCPreGood)))/numel(cell2mat(eval(strNameDPCPreGood)));
				
				%get AEs
				vecMeanAEsBad = cellfun(@mean,eval(strNameAEsPCDuringBad));
				vecSEAEsBad = cellfun(@std,eval(strNameAEsPCDuringBad))./sqrt(cellfun(@numel,eval(strNameAEsPCDuringBad)));

				vecMeanAEsGood = cellfun(@mean,eval(strNameAEsPCDuringGood));
				vecSEAEsGood = cellfun(@std,eval(strNameAEsPCDuringGood))./sqrt(cellfun(@numel,eval(strNameAEsPCDuringGood)));
		
				dblMeanPreAEsBad = mean(cell2mat(eval(strNameAEsPCPreBad)));
				dblSEPreAEsBad = std(cell2mat(eval(strNameAEsPCPreBad)))/numel(cell2mat(eval(strNameAEsPCPreBad)));
				
				dblMeanPreAEsGood = mean(cell2mat(eval(strNameAEsPCPreGood)));
				dblSEPreAEsGood = std(cell2mat(eval(strNameAEsPCPreGood)))/numel(cell2mat(eval(strNameAEsPCPreGood)));
				
				% do plotting
				%plot occurrence
				subplot(2,3,3)
				vecNrTrialsGood = nan(1,6);
				vecNrTrialsBad = nan(1,6);
				for intC=1:6
					vecNrTrialsGood(intC) = sum(cellSelectContrasts{intC}&eval(['ind' strGood]));
					vecNrTrialsBad(intC) = sum(cellSelectContrasts{intC}&eval(['ind' strBad]));
				end
				plot([vecPlotC(1) vecPlotC(end)],mean(vecTrialOccurrencesPreGood)*[1 1]/mean(vecNrTrialsGood),'color',[0 1 0],'linestyle','--');
				hold on
				plot([vecPlotC(1) vecPlotC(end)],mean(vecTrialOccurrencesPreBad)*[1 1]/mean(vecNrTrialsBad),'color',[1 0 0],'linestyle','--');
				plot(vecPlotC,vecTrialOccurrencesBad./vecNrTrialsBad,'Color',[1 0 0],'LineWidth',1.5);
				plot(vecPlotC,vecTrialOccurrencesGood./vecNrTrialsGood,'Color',[0 1 0],'LineWidth',1.5);
				hold off
				set(gca,'XScale','log')
				set(gca,'Xtick',vecPlotC,'xticklabel',vecContrasts)
				ylabel('Mean assembly occurrences per trial')
				xlabel('Stimulus contrast (%)')
				title([strRecur ' ' strGood '/' strBad]);
				
				%save data
				cellSaveAssAnal{2,intRecurType,intRespPlotType,1,intPopulation} = [vecTrialOccurrencesGood./vecNrTrialsGood;vecTrialOccurrencesBad./vecNrTrialsBad]; %normal
				cellSaveAssAnal{2,intRecurType,intRespPlotType,2,intPopulation} = [mean(vecTrialOccurrencesPreGood)/mean(vecNrTrialsGood);mean(vecTrialOccurrencesPreBad)/mean(vecNrTrialsBad)]; %pre
		
				%plot size
				subplot(2,3,4)
				errorbar([vecPlotC(1) vecPlotC(end)],dblMeanSizePreGood*[1 1],dblSESizePreGood*[1 1],'g')
				hold on
				errorbar([vecPlotC(1) vecPlotC(end)],dblMeanSizePreBad*[1 1],dblSESizePreBad*[1 1],'r')
				errorfill(vecPlotC,vecTrialMeanSizeBad,vecTrialSESizeBad,[1 0 0],[1 0.5 0.5]);
				errorfill(vecPlotC,vecTrialMeanSizeGood,vecTrialSESizeGood,[0 1 0],[0.5 1 0.5]);
				hold off
				set(gca,'XScale','log')
				set(gca,'Xtick',vecPlotC,'xticklabel',vecContrasts)
				ylabel('Mean assembly size (# of neurons)')
				xlabel('Stimulus contrast (%)')
				ylim([0 max(get(gca,'ylim'))]);
				title([strRecur ' ' strGood '/' strBad]);
				%save data
				cellSaveAssAnal{3,intRecurType,intRespPlotType,1,intPopulation} = [vecTrialMeanSizeGood;vecTrialMeanSizeBad]; %normal
				cellSaveAssAnal{3,intRecurType,intRespPlotType,2,intPopulation} = [dblMeanSizePreGood;dblMeanSizePreBad]; %pre
		
				% plot duration
				subplot(2,3,5)
				errorbar([vecPlotC(1) vecPlotC(end)],dblMeanPreDurGood*[1 1],dblSEPreDurGood*[1 1],'g')
				hold on
				errorbar([vecPlotC(1) vecPlotC(end)],dblMeanPreDurBad*[1 1],dblSEPreDurBad*[1 1],'r')
				errorfill(vecPlotC,vecMeanDurBad,vecSEDurBad,[1 0 0],[1 0.5 0.5]);
				errorfill(vecPlotC,vecMeanDurGood,vecSEDurGood,[0 1 0],[0.5 1 0.5]);
				hold off
				set(gca,'XScale','log')
				set(gca,'Xtick',vecPlotC,'xticklabel',vecContrasts)
				ylabel('Mean duration per assembly occurrence (s)')
				xlabel('Stimulus contrast (%)')
				title([strRecur ' ' strGood '/' strBad]);
				%save data
				cellSaveAssAnal{4,intRecurType,intRespPlotType,1,intPopulation} = [vecMeanDurGood;vecMeanDurBad]; %normal
				cellSaveAssAnal{4,intRecurType,intRespPlotType,2,intPopulation} = [dblMeanPreDurGood;dblMeanPreDurBad]; %pre
		
				%plot AEs
				subplot(2,3,6)
				errorbar([vecPlotC(1) vecPlotC(end)],dblMeanPreAEsGood*[1 1],dblSEPreAEsGood*[1 1],'g')
				hold on
				errorbar([vecPlotC(1) vecPlotC(end)],dblMeanPreAEsBad*[1 1],dblSEPreAEsBad*[1 1],'r')
				errorfill(vecPlotC,vecMeanAEsBad,vecSEAEsBad,[1 0 0],[1 0.5 0.5]);
				errorfill(vecPlotC,vecMeanAEsGood,vecSEAEsGood,[0 1 0],[0.5 1 0.5]);
				hold off
				set(gca,'XScale','log')
				set(gca,'Xtick',vecPlotC,'xticklabel',vecContrasts)
				ylabel('Mean AEs per assembly occurrence')
				xlabel('Stimulus contrast (%)')
				ylim([0 max(get(gca,'ylim'))]);
				title([strRecur ' ' strGood '/' strBad]);
				cellSaveAssAnal{5,intRecurType,intRespPlotType,1,intPopulation} = [vecMeanAEsGood;vecMeanAEsBad]; %normal
				cellSaveAssAnal{5,intRecurType,intRespPlotType,2,intPopulation} = [dblMeanPreAEsGood;dblMeanPreAEsBad]; %pre
			end
		end
		
		%% save data
		cellSaveAssAnal{1,1,1,1,intPopulation} = [vecMeanHit;vecMeanMiss]; %hit-miss assembly presence (time) during trials
		cellSaveAssAnal{1,2,1,1,intPopulation} = [vecMeanHitNR;vecMeanMissNR]; %hit-miss assembly presence (time) during trials
		cellSaveAssAnal{1,1,2,1,intPopulation} = [vecMeanFast;vecMeanSlow]; %hit-miss assembly presence (time) during trials
		cellSaveAssAnal{1,2,2,1,intPopulation} = [vecMeanFastNR;vecMeanSlowNR]; %hit-miss assembly presence (time) during trials
		cellSaveAssAnal{1,1,1,2,intPopulation} = [mean(vecMeanPreHitRealAssAct);mean(vecMeanPreMissRealAssAct)]; %hit-miss assembly presence (time) preceding trials
		cellSaveAssAnal{1,2,1,2,intPopulation} = [mean(vecMeanPreHitNRAssAct);mean(vecMeanPreMissNRAssAct)]; %hit-miss assembly presence (time) during trials
		cellSaveAssAnal{1,1,2,2,intPopulation} = [mean(vecMeanPreFastRealAssAct);mean(vecMeanPreSlowRealAssAct)]; %hit-miss assembly presence (time) during trials
		cellSaveAssAnal{1,2,2,2,intPopulation} = [mean(vecMeanPreFastNRAssAct);mean(vecMeanPreSlowNRAssAct)]; %hit-miss assembly presence (time) during trials

		%% save fig
		if sParams.boolSavePlots
			figure(hRealAssemblyCorrelates);
			drawnow;
			jFig = get(handle(gcf), 'JavaFrame');
			jFig.setMaximized(true);
			figure(gcf);
			drawnow;
			export_fig(sprintf('%spop%dassembly_correlates.tif',strSes,intPopulation));
			export_fig(sprintf('%spop%dassembly_correlates.pdf',strSes,intPopulation));
			
			figure(hNonRecAssemblyCorrelates);
			drawnow;
			jFig = get(handle(gcf), 'JavaFrame');
			jFig.setMaximized(true);
			figure(gcf);
			drawnow;
			export_fig(sprintf('%spop%dNRassembly_correlates.tif',strSes,intPopulation));
			export_fig(sprintf('%spop%dNRassembly_correlates.pdf',strSes,intPopulation));
			
			figure(hRealAssemblyCorrelatesSlowFast);
			drawnow;
			jFig = get(handle(gcf), 'JavaFrame');
			jFig.setMaximized(true);
			figure(gcf);
			drawnow;
			export_fig(sprintf('%spop%dassembly_correlates_slowfast.tif',strSes,intPopulation));
			export_fig(sprintf('%spop%dassembly_correlates_slowfast.pdf',strSes,intPopulation));
			
			figure(hNonRecAssemblyCorrelatesSlowFast);
			drawnow;
			jFig = get(handle(gcf), 'JavaFrame');
			jFig.setMaximized(true);
			figure(gcf);
			drawnow;
			export_fig(sprintf('%spop%dNRassembly_correlates_slowfast.tif',strSes,intPopulation));
			export_fig(sprintf('%spop%dNRassembly_correlates_slowfast.pdf',strSes,intPopulation));
			
		end
		
		%% hit-correlation
		vecRTs = cellMultiSes{intPopulation}.structStim.vecTrialRespSecs - cellMultiSes{intPopulation}.structStim.SecsOn; %get RTs for hit trials and selected contrasts
		indSelectRespTrials = ~isnan(vecRTs);
		vecRT = vecRTs(indSelectRespTrials);
		[vecRTsSorted,vecRTSortIndex] = sort(vecRT,'ascend');
		
		
		%% are there specific assemblies correlated with reaction time?
		binVector = 0.15:0.25:3;
		
		%dF/F0
		figure
		subplot(2,2,1)
		vecAct_dFoF = mean(matTrialResponse,1);
		vecActHit = vecAct_dFoF(indSelectRespTrials);
		vecActRT = vecActHit(vecRTSortIndex);
		[nVec,meanVec,stdVec,cellVals,cellIDs] = makeBins(vecRTsSorted,vecActRT,binVector);
		vecPlot = (binVector(1:end-1)+binVector(2:end))/2;
		sStats=regstats(vecActRT,vecRTsSorted,'linear',{'beta','rsquare','tstat'});
		bar(vecPlot,meanVec,1,'w','EdgeColor',[0 0 0])
		hold on
		errorbar(vecPlot,meanVec,stdVec./sqrt(nVec),'xk');
		scatter(vecRTsSorted,vecActRT);
		vecPlotX = [0.15 3];
		plot(vecPlotX,polyval(sStats.beta([2 1]),vecPlotX),'r');
		hold off
		title(sprintf('Lin reg: slope=%.3f; p=%.3f; intercept=%.3f; p=%.3f; R^2=%.3f',...
				sStats.beta(2),sStats.tstat.pval(2),sStats.beta(1),sStats.tstat.pval(1),sStats.rsquare))
		xlabel('Reaction Time (s)');
		ylabel('Mean neuronal dF/F0');
		
		%spike-detected
		subplot(2,2,2)
		vecAE = mean(matSpikeCounts,1);
		vecAEHit = vecAE(indSelectRespTrials);
		vecAERT = vecAEHit(vecRTSortIndex);
		[nVec,meanVec,stdVec,cellVals,cellIDs] = makeBins(vecRTsSorted,vecAERT,binVector);
		vecPlot = (binVector(1:end-1)+binVector(2:end))/2;
		sStats=regstats(vecAERT,vecRTsSorted,'linear',{'beta','rsquare','tstat'});
		bar(vecPlot,meanVec,1,'w','EdgeColor',[0 0 0])
		hold on
		errorbar(vecPlot,meanVec,stdVec./sqrt(nVec),'xk');
		scatter(vecRTsSorted,vecAERT);
		vecPlotX = [0.15 3];
		plot(vecPlotX,polyval(sStats.beta([2 1]),vecPlotX),'r');
		hold off
		title(sprintf('Lin reg: slope=%.3f; p=%.3f; intercept=%.3f; p=%.3f; R^2=%.3f',...
				sStats.beta(2),sStats.tstat.pval(2),sStats.beta(1),sStats.tstat.pval(1),sStats.rsquare))
		xlabel('Reaction Time (s)');
		ylabel('Mean neuronal Activation Events');
		
		%assembly presence
		subplot(2,2,3)
		vecAssPresence = mean(matRespRealAssemblies,1);%vecRespNonRecurring;%mean(matRespRealAssemblies,1);
		vecAssHit = vecAssPresence(indSelectRespTrials);
		vecAssRT = vecAssHit(vecRTSortIndex);
		[nVec,meanVec,stdVec,cellVals,cellIDs] = makeBins(vecRTsSorted,vecAssRT,binVector);
		vecPlot = (binVector(1:end-1)+binVector(2:end))/2;
		sStats=regstats(vecAssRT,vecRTsSorted,'linear',{'beta','rsquare','tstat'});
		bar(vecPlot,meanVec,1,'w','EdgeColor',[0 0 0])
		hold on
		errorbar(vecPlot,meanVec,stdVec./sqrt(nVec),'xk');
		scatter(vecRTsSorted,vecAssRT);
		vecPlotX = [0.15 3];
		plot(vecPlotX,polyval(sStats.beta([2 1]),vecPlotX),'r');
		hold off
		title(sprintf('Lin reg: slope=%.3f; p=%.3f; intercept=%.3f; p=%.3f; R^2=%.3f',...
				sStats.beta(2),sStats.tstat.pval(2),sStats.beta(1),sStats.tstat.pval(1),sStats.rsquare))
		xlabel('Reaction Time (s)');
		ylabel('Mean recurring assembly presence');
		title(sprintf('Recurring; %s pop%d',strSes,intPopulation));
		
		subplot(2,2,4)
		title(sprintf('Non-recurring; %s pop%d',strSes,intPopulation));
		vecNRAssPresence = mean(matRespNonRecurring,1);
		vecNRAssHit = vecNRAssPresence(indSelectRespTrials);
		vecNRAssRT = vecNRAssHit(vecRTSortIndex);
		[nVec,meanVec,stdVec,cellVals,cellIDs] = makeBins(vecRTsSorted,vecAssRT,binVector);
		vecPlot = (binVector(1:end-1)+binVector(2:end))/2;
		sStats=regstats(vecAssRT,vecRTsSorted,'linear',{'beta','rsquare','tstat'});
		bar(vecPlot,meanVec,1,'w','EdgeColor',[0 0 0])
		hold on
		errorbar(vecPlot,meanVec,stdVec./sqrt(nVec),'xk');
		scatter(vecRTsSorted,vecNRAssRT);
		vecPlotX = [0.15 3];
		plot(vecPlotX,polyval(sStats.beta([2 1]),vecPlotX),'r');
		hold off
		title(sprintf('Lin reg: slope=%.3f; p=%.3f; intercept=%.3f; p=%.3f; R^2=%.3f',...
				sStats.beta(2),sStats.tstat.pval(2),sStats.beta(1),sStats.tstat.pval(1),sStats.rsquare))
		xlabel('Reaction Time (s)');
		ylabel('Mean non-recurring assembly presence');
		
		
		%save data
		cellSaveRTDependency{1,intPopulation} = vecRTsSorted;
		cellSaveRTDependency{2,intPopulation} = vecActRT;
		cellSaveRTDependency{3,intPopulation} = vecAERT;
		cellSaveRTDependency{4,intPopulation} = vecAssRT;
		cellSaveRTDependency{5,intPopulation} = vecNRAssRT;
		
		%save fig
		if sParams.boolSavePlots
			drawnow;
			jFig = get(handle(gcf), 'JavaFrame');
			jFig.setMaximized(true);
			figure(gcf);
			drawnow;
			export_fig(sprintf('%spop%dassembly_RT.tif',strSes,intPopulation));
			export_fig(sprintf('%spop%dassembly_RT.pdf',strSes,intPopulation));
		end
		
		%relative occurrence of assemblies during miss/fast/slow
		
		%consistency of assemblies themselves during miss/fast/slow
		
		%size of assemblies during miss/fast/slow
		
		
		%contrast tuning of assemblies
		
		%correlation assemblies with pupil size
		
	end
	
	%% save data structures
	if sParams.boolSaveData
		vecClock = fix(clock);
		strDate = num2str(vecClock(1:3));
		strFile = ['dataAssAnal' strSes '_' strrep(strDate,'    ','_')];
		save(['D:\Data\Results\spikeAnalysis\' strFile],'cellSave*','-v7.3');
		cd(strOldDir);
	end
end