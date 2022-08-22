%function vecPredictorAssemblyZscores = doPlotActivityDistributions(cellContrastDataAD,sMetaDataAD,sPC,sParams)

%if ~exist('sParams','var'),sParams = struct;end
%if isfield(sParams,'intStartSes'),intStartSes=sParams.intStartSes;else intStartSes=1;end
%if isfield(sParams,'intStopSes'),intStopSes=sParams.intStopSes;else intStopSes=length(cellContrastDataAD);end
close all
for intMouse=[4]
	
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
		load(['D:\Data\Results\stimdetection\dataPreProAggregate' strSes '.mat']);
	end
	
	%vecBlock = [1 1 1 2 2 2 2 2]; %define whether neurons are in same population or not
	vecBlockTypes = unique(vecBlock);
	intNumBlocks = length(vecBlockTypes);
	%vecNeuronNum = zeros(1,intNumBlocks);
	%cellKeepList = cell(1,intNumBlocks);
	%#ok<*ASGLU>
	%#ok<*AGROW>
	clear sLoad sSesAggregate ses
	sParams.strFigDir = ['D:\Data\Results\stimdetection' filesep strSes filesep];
	sParams.boolSavePlots = true;
	sParams.boolSaveData = false;
	strOldDir = cd(sParams.strFigDir);
	
	%change name for no split
	strSes = ['NS_Extra_' strSes];
	
	for intPopulation = vecBlockTypes
		%get neuronal tuning
		if intMouse==8
			%remove last trial
			vecRem = true(size(cellMultiSes{1}.structStim.Orientation));
			vecRem((end-47):end) = false;
			cellMultiSes{1}.structStim = remel(cellMultiSes{1}.structStim,vecRem);
			%recalc dfof
			%cellMultiSes{1} = doRecalcdFoF(cellMultiSes{1},3);
		end
		
		%[indTuned,vecNeuronPrefStim,vecOSI,cellAggregate] = getTunedStimDetectionNeurons(cellSes(vecBlock==intPopulation));
		[indTuned,vecNeuronPrefStim,vecOSI,cellAggregate] = getTunedStimDetectionNeurons(cellMultiSes{intPopulation});
		structStim = cellMultiSes{intPopulation}.structStim;
		structStim.vecTrialRepetition = ceil((1:(length(structStim.Orientation)))/48);
		
		% remove trials with reaction time <100ms
		dblRemSecs = 0.15;
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
		%structStim.FrameOff = structStim.FrameOn + 10;
		
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
		intOris = length(vecOrientations);
		
		%% get pre-formatted data variables
		%get neuronal responses per trial
		[matTrialResponse,cellSelectContrasts] = getTrialResponseData(cellMultiSes{intPopulation},structStim);
		intContrasts = length(cellSelectContrasts);
		
		% get stimulus data
		[vecStimOris,dummy]=find(bsxfun(@eq,structStim.Orientation,vecOrientations'));
		vecStimOris = vecStimOris';
		
		[vecStimContrasts,dummy]=find(bsxfun(@eq,structStim.Contrast,unique(structStim.Contrast)'));
		vecStimContrasts = vecStimContrasts';
		vecTrialRepetition = structStim.vecTrialRepetition;
		
		vecRTs = cellMultiSes{intPopulation}.structStim.vecTrialRespSecs - cellMultiSes{intPopulation}.structStim.SecsOn; %get RTs for hit trials and selected contrasts
		indSelectRespTrials = ~isnan(vecRTs);
		intTrials = length(indSelectRespTrials);
		
		%% pre-analyses
		%normalize
		matRespNormPerNeuron = zscore(matTrialResponse,[],2);
		matSelect = tril(true(intNeurons,intNeurons),-1);
		vecHeterogeneity = nan(intTrials,1);
		vecHetRawAct = nan(intTrials,1);
		vecHetMultiDim = nan(intTrials,1);
		vecActPrefPop = nan(intTrials,1);
		vecActWholePop = nan(intTrials,1);
		vecZActAll = nan(intTrials,1);
		for intTrial=1:intTrials
			vecZ = matRespNormPerNeuron(:,intTrial);
			matZ1 = repmat(vecZ,[1 intNeurons]);
			matZ2 = repmat(vecZ',[intNeurons 1]);
			matDist = abs(matZ1 - matZ2);
			vecHeterogeneity(intTrial) = mean(matDist(matSelect));
			%vecHetMultiDim(intTrial) = sqrt(sum(abs(vecZ-mean(vecZ))))/length(vecZ);
			vecHetMultiDim(intTrial) = sqrt(sum((vecZ-mean(vecZ)).^2))/length(vecZ);
			
			vecAct = matTrialResponse(:,intTrial);
			matA1 = repmat(vecAct,[1 intNeurons]);
			matA2 = repmat(vecAct',[intNeurons 1]);
			matDistA = abs(matA1 - matA2);
			vecHetRawAct(intTrial) = mean(matDistA(matSelect));
			
			%get pref pop
			%intTrialOri = [];
			vecPrefNeurons = vecNeuronPrefStim == vecStimOris(intTrial);
			vecActPrefPop(intTrial) = mean(vecAct(vecPrefNeurons));
			vecActWholePop(intTrial) = mean(vecAct);
			vecZActAll(intTrial) = mean(matRespNormPerNeuron(:,intTrial));
			
		end
		
		%% analyze orthogonal vs parallel properties
		figure
		subplot(2,2,1)
		%scatter3(vecZActAll,vecHeterogeneity,vecHetMultiDim)
		%xlabel('Z-scored dF/F0');
		%ylabel('Heterogeneity');
		%zlabel('Multi-dimensional heterogeneity');
		
		vecMeanActResp = nan(1,6);
		vecSEActResp = nan(1,6);
		vecMeanHResp = nan(1,6);
		vecSEHResp = nan(1,6);
		vecMeanActNoResp = nan(1,6);
		vecSEActNoResp = nan(1,6);
		vecMeanHNoResp = nan(1,6);
		vecSEHNoResp = nan(1,6);
		vecHitMissSeparabilityH = nan(1,6);
		vecHitMissSeparabilityA = nan(1,6);
		for intC=1:6
			vecHet = vecHeterogeneity;%vecHetMultiDim;%vecHeterogeneity
			vecAct = vecActWholePop;%vecZActAll;%vecActPrefPop
			subplot(2,3,intC)
			indC = cellSelectContrasts{intC};
			h= scatter(vecAct(indC)',vecHet(indC)','CDATA',1-double(indSelectRespTrials(indC)));
			colormap(redgreen);
			xlabel('Mean activity (dF/F0)');
			ylabel('Heterogeneity');
			
			%get data
			vecRespA = vecAct(indC & indSelectRespTrials);
			vecNoRespA = vecAct(indC & ~indSelectRespTrials);
			vecRespH = vecHet(indC & indSelectRespTrials);
			vecNoRespH = vecHet(indC & ~indSelectRespTrials);
			
			%mean+sem
			vecMeanActResp(intC) = mean(vecRespA);
			vecSEActResp(intC) = std(vecRespA)/sqrt(length(vecRespA));
			vecMeanHResp(intC) =  mean(vecRespH);
			vecSEHResp(intC) = std(vecRespH)/sqrt(length(vecRespH));
			vecMeanActNoResp(intC) =  mean(vecNoRespA);
			vecSEActNoResp(intC) = std(vecNoRespA)/sqrt(length(vecNoRespA));
			vecMeanHNoResp(intC) =  mean(vecNoRespH);
			vecSEHNoResp(intC) = std(vecNoRespH)/sqrt(length(vecNoRespH));
			
			%separability
			indTheseTrialsResp = indSelectRespTrials(indC);
			
			%het
			vecTheseTrialsHet =  vecHet(indC);
			if numel(unique(indTheseTrialsResp))==1
				AUC=1;
				h= scatter(vecAct(indC)',vecHet(indC)');
			else
				[X,Y,T,AUC] = perfcurve(indTheseTrialsResp',vecTheseTrialsHet,true);
			end
			vecHitMissSeparabilityH(intC) = AUC;
			
			%df/f
			vecTheseTrialsAct =  vecAct(indC);
			if numel(unique(indTheseTrialsResp))==1,AUC=1;else
				[X,Y,T,AUC] = perfcurve(indTheseTrialsResp',vecTheseTrialsAct,true);
			end
			vecHitMissSeparabilityA(intC) = AUC;
			
			title(sprintf('Separability; Act: %.3f; Het: %.3f',vecHitMissSeparabilityA(intC),vecHitMissSeparabilityH(intC)));
			
			%% =======> calculate mean hit/miss for different measures and save to output variable
			%vecHeterogeneity = nan(intTrials,1);
			%vecHetRawAct = nan(intTrials,1);
			%vecHetMultiDim = nan(intTrials,1);
			%vecActPrefPop = nan(intTrials,1);
			%vecActWholePop = nan(intTrials,1);
			%vecZActAll = nan(intTrials,1);

		end
		%continue
		cellSaveSeparability(intPopulation,1) = mean(vecHitMissSeparabilityH(2:5));
		cellSaveSeparability(intPopulation,2) = mean(vecHitMissSeparabilityA(2:5));
		
		%{
		figure
		subplot(2,2,1)
		vecPlotC=[0.2 0.5 2 8 32 100];
		errorfill(vecPlotC,vecMeanActNoResp,vecSEActNoResp,[1 0 0],[1 0.5 0.5]);
		hold on
		errorfill(vecPlotC,vecMeanActResp,vecSEActResp,[0 1 0],[0.5 1 0.5]);
		hold off
		ylabel('dF/F0')
		xlabel('Stimulus contrast (%)')
		set(gca,'xscale','log');
		
		subplot(2,2,2)
		vecPlotC=[0.2 0.5 2 8 32 100];
		errorfill(vecPlotC,vecMeanHNoResp,vecSEHNoResp,[1 0 0],[1 0.5 0.5]);
		hold on
		errorfill(vecPlotC,vecMeanHResp,vecSEHResp,[0 1 0],[0.5 1 0.5]);
		hold off
		ylabel('Heterogeneity')
		xlabel('Stimulus contrast (%)')
		set(gca,'xscale','log');
		%}
		
		
		%% test for multidimensional non-uniformity
		%calculate distribution of pairwise inter-point distances in
		%multidimensional space
		
		%type (norm or raw)
		intType = 1;
		indTestTrials = vecStimContrasts ~= 1 & vecStimContrasts ~= 6;
		indTestRespTrials = indSelectRespTrials(indTestTrials);
		if intType == 1
			matTestResps = matTrialResponse(:,indTestTrials);
			%binning vector
			dblStep = 0.02; %dF/F: 0.02
			vecBins = 0:dblStep:1.4;
			strLabelX = 'Distance (dF/F0)';
		else
			matTestResps = matRespNormPerNeuron(:,indTestTrials);
			%binning vector
			dblStep = 0.4; %dF/F: 0.02
			vecBins = 0:dblStep:30;
			strLabelX = 'Distance (z-score units)';
		end
		intTestTrials = size(matTestResps,2);
		matPairwiseDistances = nan(intTestTrials,intTestTrials);
		for intTrial=1:intTestTrials
			matPairwiseDistances(intTrial,:) = sqrt(sum(bsxfun(@minus,matTestResps,matTestResps(:,intTrial)).^2,1));
		end
		matSelectTrials = tril(true(intTestTrials,intTestTrials),-1);
		matSelectHitTrials = matSelectTrials;
		matSelectHitTrials(~indTestRespTrials,~indTestRespTrials) = false;
		matSelectMissTrials = matSelectTrials;
		matSelectMissTrials(indTestRespTrials,indTestRespTrials) = false;
		
		figure
		%subplot(2,2,1)
		%imagesc(matPairwiseDistances);colormap(hot);colorbar;
		%title('Pairwise (trials) Euclidian distances');
		%xlabel('Trial')
		%ylabel('Trial')
		
		subplot(2,2,1);
		vecDistroAll = hist(matPairwiseDistances(matSelectTrials),vecBins);
		vecDistroAll = vecDistroAll/sum(vecDistroAll);
		vecDistroHit = hist(matPairwiseDistances(matSelectHitTrials),vecBins);
		vecDistroHit = vecDistroHit/sum(vecDistroHit);
		vecDistroMiss = hist(matPairwiseDistances(matSelectMissTrials),vecBins);
		vecDistroMiss = vecDistroMiss/sum(vecDistroMiss);
		
		
		stairs(vecBins,vecDistroAll,'Color','k');
		hold on
		stairs(vecBins,vecDistroHit,'Color','g');
		stairs(vecBins,vecDistroMiss,'Color','r');
		hold off
		legend({'all','hit','miss'});
		xlabel(strLabelX);
		ylabel('Normalized count of trial-pairs');
		title(sprintf('Mean distances, all=%.3f, hit=%.3f, miss=%.3f',mean(matPairwiseDistances(matSelectTrials)),mean(matPairwiseDistances(matSelectHitTrials)),mean(matPairwiseDistances(matSelectMissTrials))));
		
		%for each point (i.e., trial), reflect/mirror that neuron through
		%the multidimensional minimal-heterogeneity axis (i.e., mean
		%response gradient)
		matMirrorDistroAll = nan(intTestTrials,length(vecBins));
		matMirrorDistroHit = nan(intTestTrials,length(vecBins));
		matMirrorDistroMiss = nan(intTestTrials,length(vecBins));
		for intTrial = 1:intTestTrials;
			vecPoint = matTestResps(:,intTrial);
			dblMean = mean(vecPoint);
			vecRel = vecPoint - dblMean;
			vecMirrored = -vecRel + dblMean;
			matMirrorResps = matTestResps;
			matMirrorResps(:,intTrial) = vecMirrored;
			
			%recalculate the distribution of pairwise inter-point distances in
			%multidimensional space
			matPairwiseDistancesMirrored = matPairwiseDistances;
			vecMirrorDistances = sqrt(sum(bsxfun(@minus,matMirrorResps,matMirrorResps(:,intTrial)).^2,1));
			matPairwiseDistancesMirrored(intTrial,:) = vecMirrorDistances;
			matPairwiseDistancesMirrored(:,intTrial) = vecMirrorDistances;
			
			vecMirrorDistroAll = hist(matPairwiseDistancesMirrored(matSelectTrials),vecBins);
			vecMirrorDistroAll = vecMirrorDistroAll/sum(vecMirrorDistroAll);
			vecMirrorDistroHit = hist(matPairwiseDistancesMirrored(matSelectHitTrials),vecBins);
			vecMirrorDistroHit = vecMirrorDistroHit/sum(vecMirrorDistroHit);
			vecMirrorDistroMiss = hist(matPairwiseDistancesMirrored(matSelectMissTrials),vecBins);
			vecMirrorDistroMiss = vecMirrorDistroMiss/sum(vecMirrorDistroMiss);
			
			matMirrorDistroAll(intTrial,:) = vecMirrorDistroAll-vecDistroAll;
			matMirrorDistroHit(intTrial,:) = vecMirrorDistroHit-vecDistroHit;
			matMirrorDistroMiss(intTrial,:) = vecMirrorDistroMiss-vecDistroMiss;
			
		end
		%compare this distribution to the original; if changes occur, the
		%population response is asymmetricaly distributed around the axis,
		%and specific multidimensional response clusters that correspond to
		%hit-trials are present. if no changes occur, the population
		%response is symmetric, and heterogeneity is the main factor
		%determining stimulus detection
		
		
		vecMirrorDistroAll = mean(matMirrorDistroAll,1);
		%vecMirrorDistroAll = vecMirrorDistroAll/sum(vecMirrorDistroAll);
		vecMirrorDistroHit = mean(matMirrorDistroHit,1);
		%vecMirrorDistroHit = vecMirrorDistroHit/sum(vecMirrorDistroHit);
		vecMirrorDistroMiss = mean(matMirrorDistroMiss,1);
		%vecMirrorDistroMiss = vecMirrorDistroMiss/sum(vecMirrorDistroMiss);
		
		subplot(2,2,2);
		%stairs(vecBins,vecDistroAll,'Color','k');
		hold on
		stairs(vecBins,vecMirrorDistroAll,'Color','b');
		hold off
		xlabel(strLabelX);
		ylabel('Normalized count of trial-pairs');
		title(sprintf('Mean distance All trials; Mirrored-Raw'));%, data=%.3f, mirrored=%.3f',mean(matPairwiseDistances(matSelectTrials)),mean(matPairwiseDistances(matSelectHitTrials)),mean(matPairwiseDistances(matSelectMissTrials))));
		
		subplot(2,2,3);
		%stairs(vecBins,vecDistroHit,'Color','k');
		hold on
		stairs(vecBins,vecMirrorDistroHit,'Color','b');
		hold off
		xlabel(strLabelX);
		ylabel('Normalized count of trial-pairs');
		title(sprintf('Mean distance Hit trials; Mirrored-Raw'));%, data=%.3f, mirrored=%.3f',mean(matPairwiseDistances(matSelectTrials)),mean(matPairwiseDistances(matSelectHitTrials)),mean(matPairwiseDistances(matSelectMissTrials))));
		
		subplot(2,2,4);
		%stairs(vecBins,vecDistroMiss,'Color','k');
		hold on
		stairs(vecBins,vecMirrorDistroMiss,'Color','b');
		hold off
		xlabel(strLabelX);
		ylabel('Normalized count of trial-pairs');
		title(sprintf('Mean distance Miss trials; Mirrored-Raw'));%, data=%.3f, mirrored=%.3f',mean(matPairwiseDistances(matSelectTrials)),mean(matPairwiseDistances(matSelectHitTrials)),mean(matPairwiseDistances(matSelectMissTrials))));
		
		%save data
		cellSaveMultiDimDistanceDistros(1,:,intPopulation) = vecDistroAll;
		cellSaveMultiDimDistanceDistros(2,:,intPopulation) = vecDistroHit;
		cellSaveMultiDimDistanceDistros(3,:,intPopulation) = vecDistroMiss;
		cellSaveMultiDimDistanceDistros(4,:,intPopulation) = vecMirrorDistroAll;
		cellSaveMultiDimDistanceDistros(5,:,intPopulation) = vecMirrorDistroHit;
		cellSaveMultiDimDistanceDistros(6,:,intPopulation) = vecMirrorDistroMiss;
		
		
		%save figure
		if sParams.boolSavePlots
			drawnow;
			jFig = get(handle(gcf), 'JavaFrame');
			jFig.setMaximized(true);
			figure(gcf);
			drawnow;
			strFig = sprintf('%sdistance_multidimensional_neural_response_space_pop%d',strSes,intPopulation);
			export_fig([strFig '.tif']);
			export_fig([strFig '.pdf']);
		end
		
		%% remove mean, heterogeneity, none or both
		matMeanRemovedTestResps = bsxfun(@minus,matTestResps,mean(matTestResps,1));
		matHeteroRemovedTestResps = bsxfun(@rdivide,matTestResps,sqrt(sum(matTestResps.^2)));%matMeanRemovedTestResps
		%matHeteroMeanRemovedTestResps = bsxfun(@rdivide,matMeanRemovedTestResps,sqrt(sum(matMeanRemovedTestResps.^2)));
		matHeteroMeanRemovedTestResps = bsxfun(@minus,matHeteroRemovedTestResps,mean(matHeteroRemovedTestResps,1));
		
		figure
		colormap(hot)
		subplot(2,2,1)
		imagesc(matTestResps)
		title('Original data')
		xlabel('Trial')
		ylabel('Neuron')
		
		subplot(2,2,2)
		imagesc(matMeanRemovedTestResps)
		title('Mean removed')
		xlabel('Trial')
		ylabel('Neuron')
		
		subplot(2,2,3)
		imagesc(matHeteroRemovedTestResps)
		title('Heterogeneity removed')
		xlabel('Trial')
		ylabel('Neuron')
		
		subplot(2,2,4)
		imagesc(matHeteroMeanRemovedTestResps)
		title('Mean+heterogeneity removed')
		xlabel('Trial')
		ylabel('Neuron')
		
		%save figure
		if sParams.boolSavePlots
			drawnow;
			jFig = get(handle(gcf), 'JavaFrame');
			jFig.setMaximized(true);
			figure(gcf);
			drawnow;
			strFig = sprintf('%sneural_responses_het_mean_removed_pop%d',strSes,intPopulation);
			export_fig([strFig '.tif']);
			export_fig([strFig '.pdf']);
		end
		
		%% decode hit/miss
		%raw data (mean+heterogeneity+multidim)
		[dblPerformanceRawML,vecDecodedIndexCV,matPosteriorProbabilityCV] = doCrossValidatedDecodingML(matTestResps',indTestRespTrials,1);
		[dummy,vecCI_Raw] = binofit(dblPerformanceRawML*intTestTrials,intTestTrials);
		
		%heterogeneity removed (mean+multidim)
		[dblPerformanceNoHetML,vecDecodedIndexCV,matPosteriorProbabilityCV] = doCrossValidatedDecodingML(matHeteroRemovedTestResps',indTestRespTrials,1);
		[dummy,vecCI_NoHet] = binofit(dblPerformanceNoHetML*intTestTrials,intTestTrials);
		
		%mean removed (heterogeneity+multidim)
		[dblPerformanceNoMeanML,vecDecodedIndexCV,matPosteriorProbabilityCV] = doCrossValidatedDecodingML(matMeanRemovedTestResps',indTestRespTrials,1);
		[dummy,vecCI_NoMean] = binofit(dblPerformanceNoMeanML*intTestTrials,intTestTrials);
		
		%mean+heterogeneity removed (only multidim)
		[dblPerformanceNoHetMeanML,vecDecodedIndexCV,matPosteriorProbabilityCV] = doCrossValidatedDecodingML(matHeteroMeanRemovedTestResps',indTestRespTrials,1);
		[dummy,vecCI_NoHetMean] = binofit(dblPerformanceNoHetMeanML*intTestTrials,intTestTrials);
		
		vecMean = [dblPerformanceRawML dblPerformanceNoMeanML dblPerformanceNoHetML dblPerformanceNoHetMeanML];
		figure
		errorbar(1:4,vecMean,([vecCI_Raw(1) vecCI_NoMean(1) vecCI_NoHet(1) vecCI_NoHetMean(1)]-vecMean)/sqrt(intTestTrials),([vecCI_Raw(2) vecCI_NoMean(2) vecCI_NoHet(2) vecCI_NoHetMean(2)]-vecMean)/sqrt(intTestTrials),'xb');
		ylim([0.5 0.7])
		set(gca,'xtick',1:4,'xticklabel',{'Original Data','Mean removed','Het removed','Both removed'});
		ylabel('Hit/miss decoding accuracy');
		
		%save figure
		if sParams.boolSavePlots
			drawnow;
			strFig = sprintf('%sdec_perf_neural_responses_het_mean_removed_pop%d',strSes,intPopulation);
			export_fig([strFig '.tif']);
			export_fig([strFig '.pdf']);
		end
		
		%save data
		cellSaveHitMissDecoding(:,intPopulation) = vecMean;
		
		%% =================================================================================> compare stimulus presence decoding with hit/miss behavior
		
		
		%% analysis
		%pre-allocate
		vecContrasts = unique(cellMultiSes{1}.structStim.Contrast)*100;
		vecContrasts(1) = 0.2;
		intContrasts = length(vecContrasts);
		intRepetitions = max(vecTrialRepetition);
		vecBehavRepD_p = nan(1,intRepetitions);
		matCohenRepHet = nan(intContrasts,intRepetitions);
		matCohenRepAct = nan(intContrasts,intRepetitions);
		
		for intRep=1:intRepetitions
			indRepetition = vecTrialRepetition == intRep;
			vecC0 = structStim.vecTrialResponse(indRepetition & cellSelectContrasts{1});
			vecC100 = structStim.vecTrialResponse(indRepetition & cellSelectContrasts{end});
			vecBehavRepD_p(intRep) = sum(vecC100)/numel(vecC100) - sum(vecC0)/numel(vecC0);
			
			for intContrastIndex=1:intContrasts
				%get resp het
				vecHeteroHit = vecHeterogeneity(cellSelectContrasts{intContrastIndex} & indSelectRespTrials & indRepetition);
				vecHeteroMiss = vecHeterogeneity(cellSelectContrasts{intContrastIndex} & ~indSelectRespTrials & indRepetition);
				
				%get resp act
				vecActHit = vecActPrefPop(cellSelectContrasts{intContrastIndex} & indSelectRespTrials & indRepetition);
				vecActMiss = vecActPrefPop(cellSelectContrasts{intContrastIndex} & ~indSelectRespTrials & indRepetition);
				
				%cohen's d
				matCohenRepAct(intContrastIndex,intRep) = getCohensD(vecActHit',vecActMiss');
				matCohenRepHet(intContrastIndex,intRep) = getCohensD(vecHeteroHit',vecHeteroMiss');
			end
		end
		vecCohenRepAct = nanmean(matCohenRepAct(2:5,:),1);
		vecCohenRepHet = nanmean(matCohenRepHet(2:5,:),1);
		
		cellSaveCohenBlocks{intPopulation,1} = [vecBehavRepD_p;vecCohenRepAct;vecCohenRepHet];
		cellSaveCohenBlocks{intPopulation,2} = matCohenRepAct;
		cellSaveCohenBlocks{intPopulation,3} = matCohenRepHet;
		
		%{
		scatter(vecBehavRepD_p,vecCohenRepHet,'kx');
		hold on
		scatter(vecBehavRepD_p,vecCohenRepAct,'rx');
		hold off;
		drawnow
		%pause
		%}
		
		%% no reps, across all
		%pre-allocate
		vecCohenRepHetAll = nan(intContrasts,intRepetitions);
		vecCohenRepActAll = nan(intContrasts,intRepetitions);
		
		
		vecC0 = structStim.vecTrialResponse(cellSelectContrasts{1});
		vecC100 = structStim.vecTrialResponse(cellSelectContrasts{end});
		dblBehavRepD_p = sum(vecC100)/numel(vecC100) - sum(vecC0)/numel(vecC0);
		
		for intContrastIndex=1:intContrasts
			%get resp het
			vecHeteroHit = vecHeterogeneity(cellSelectContrasts{intContrastIndex} & indSelectRespTrials);
			vecHeteroMiss = vecHeterogeneity(cellSelectContrasts{intContrastIndex} & ~indSelectRespTrials);
			
			%get resp act
			vecActHit = vecActPrefPop(cellSelectContrasts{intContrastIndex} & indSelectRespTrials);
			vecActMiss = vecActPrefPop(cellSelectContrasts{intContrastIndex} & ~indSelectRespTrials);
			
			%cohen's d
			vecCohenRepActAll(intContrastIndex) = getCohensD(vecActHit',vecActMiss');
			vecCohenRepHetAll(intContrastIndex) = getCohensD(vecHeteroHit',vecHeteroMiss');
		end
		
		dblCohenAct = nanmean(vecCohenRepActAll(2:5));
		dblCohenHet = nanmean(vecCohenRepHetAll(2:5));
		
		cellSaveCohenBlocks{intPopulation,1} = [vecBehavRepD_p;vecCohenRepAct;vecCohenRepHet];
		cellSaveCohenBlocks{intPopulation,2} = [dblBehavRepD_p;dblCohenAct;dblCohenHet];
		cellSaveCohenBlocks{intPopulation,3} = matCohenRepAct;
		cellSaveCohenBlocks{intPopulation,4} = matCohenRepHet;
		continue
		%scatter(vecBehavRepD_p,vecCohenRepHet,'kx');
		%hold on
		%scatter(vecBehavRepD_p,vecCohenRepAct,'rx');
		%hold off;
		%drawnow
		%pause
		
		%% z-drift analysis
		%[3 4 5 7 8]
		%if is file
		if ~exist(['D:\Data\Results\stimdetection\' strSes(end-7:end) '_zdrift_aggregate.mat'],'file')
			continue;
		end
		
		%load file
		intNumFrames = length(cellMultiSes{intPopulation}.neuron(1).dFoF);
		sLoad = load(['D:\Data\Results\stimdetection\' strSes(end-7:end) '_zdrift_aggregate.mat']);
		vecDriftZ = sLoad.vecDriftZ;
		%split if necessary
		if length(vecBlockTypes) > 1
			if intPopulation == 1
				vecDriftZ = vecDriftZ(1:intNumFrames);
			else
				vecDriftZ = vecDriftZ((end-intNumFrames+1):end);
			end
		end
		vecShiftZ = abs(diff(vecDriftZ));
		vecOffsetZ = abs(vecDriftZ - mean(vecDriftZ(:)));
		
		
		
		%get hit-response times
		vecHitTimes = cellMultiSes{intPopulation}.structStim.FrameOff(indSelectRespTrials);
		sEvents.vecOn = vecHitTimes;
		sEvents.dblFrameRate = 25.4;
		sEvents.vecWindowSecs = [-5 3];
		sEvents.handleFig = 1;
		
		%plot
		figure
		subplot(2,2,1);
		[cellHandles,sOutOffsetZ] = doPEP(sEvents,vecOffsetZ);
		ylabel('Offset from mean z-plane (micron)')
		ylim([0 2]);
		xlabel('Time after hit response (s)');
		
		subplot(2,2,2);
		[cellHandles,sOutShiftZ] = doPEP(sEvents,vecShiftZ);
		ylabel('Inter-frame z-shift (d(micron))')
		ylim([0 1]);
		xlabel('Time after hit response (s)');
		
		
		%make PSTH of z-shifts
		matAct = zeros(intNeurons,intNumFrames);
		for intNeuron=1:intNeurons
			matAct(intNeuron,:) = cellMultiSes{intPopulation}.neuron(intNeuron).dFoF;
		end
		[vecHeterogeneity,vecActivity] = calcMatRespHeteroGen(matAct);
		
		subplot(2,2,3);
		[cellHandles,sOutdFoF] = doPEP(sEvents,vecActivity);
		ylabel('Mean pop. dF/F0')
		%ylim([0 1]);
		xlabel('Time after hit response (s)');
		
		subplot(2,2,4);
		[cellHandles,sOutHet] = doPEP(sEvents,vecHeterogeneity);
		ylabel('Heterogeneity')
		%ylim([0 1]);
		xlabel('Time after hit response (s)');
		
		%save data
		cellSaveLickZDrift{intPopulation,1} = sOutOffsetZ.vecLineY;
		cellSaveLickZDrift{intPopulation,2} = sOutShiftZ.vecLineY;
		cellSaveLickZDrift{intPopulation,3} = sOutdFoF.vecLineY;
		cellSaveLickZDrift{intPopulation,4} = sOutHet.vecLineY;
		
		%save plot
		if sParams.boolSavePlots
			drawnow;
			jFig = get(handle(gcf), 'JavaFrame');
			jFig.setMaximized(true);
			figure(gcf);
			drawnow;
			strFig = sprintf('%spop%d_lick_induced_zdrift',strSes,intPopulation);
			export_fig([strFig '.tif']);
			export_fig([strFig '.pdf']);
		end
	end
	
	%% save data structures
	if sParams.boolSaveData
		vecClock = fix(clock);
		strDate = num2str(vecClock(1:3));
		strFile = ['data_aggregate' strSes '_' strrep(strDate,'    ','_')];
		%save(['D:\Data\Results\stimdetection\' strFile],'-v7.3');
		%if exist(['D:\Data\Results\stimdetection\' strFile '.mat'],'file')
		%	load(['D:\Data\Results\stimdetection\' strFile]);
		%end
		cellVariables=who('cellSave*');
		if isempty(cellVariables)
			warning([mfilename ':NoVariablesToSave'],'No variables found to be saved to disk');
		else
			save(['D:\Data\Results\stimdetection\' strFile],'cellSave*','-v7.3');
		end
		cd(strOldDir);
	end
end
%{
figures to save and take mean over animals
- cellSavePreStim: dF/F detect/no-detect preceding stimulus frames [5 (c) x 2 (d/nd) cell array] x n (animals) with every cell = vector of values
- cellSaveBehavDetect: behavioral detection [5 (c) x n (animals)] cell array with every cell = vector of detect/no detect values (1/0)
- cellSaveBehavRT: behavioral RT [5 (c) x n (animals)] cell array with every cell = vector of RTs
- cellSaveBaseSubtrStimAct: dF/F detect/no-detect baseline-subtracted activity [5 (c) x 2 (d/nd) cell array] x n (animals) with every cell = vector of values
- cellSaveDuringStim: dF/F detect/no-detect during stimulus frames [5 (c) x 2 (d/nd) cell array] x n (animals) with every cell = vector of values
- cellSaveTraceAct: dF/F trace detect/no-detect pref/non-pref [2 (pref/no-pref) x 2 (detect/no-detect) x 6 (c)] x n (animals) with every cell = [Nreps x Nframes] matrix of trace dF/F values
-  AD: cellSaveMatrices: blob matrices in [1 x 4 (detect/no-detect/diff/significance)] cell array with [5 (c) x 31 (bins)] matrix
-  AD: cellSaveHCAENeuronConsistency: [1 x i (neuronal population)] cell array with vector of shuffled correlation values, where first value is actual data
-  AD: cellSaveSignalCorrs: signal correlations [3 (within-HCAE, between, within-non-HCAE) x n (animals/blocks)] with every cell = vector of pairwise correlation values
-  AD: cellSaveNoiseCorrs: noise correlations [3 (within-HCAE, between, within-non-HCAE) x n (animals/blocks)] with every cell = vector of pairwise correlation values
-  AD: cellSaveHCAE: normalized increase in dF/F for detection trials [5 (c) x n (animals/blocks)] cell array with each cell = vector of neurons with normalized dF/F increase for detection trials
-  AD: cellSaveSignalCorrsBiDir: signal correlations [3 (within-HCAE, between, within-non-HCAE) x n (animals/blocks)] with every cell = vector of pairwise correlation values
-  AD: cellSaveNoiseCorrsBiDir: noise correlations [3 (within-HCAE, between, within-non-HCAE) x n (animals/blocks)] with every cell = vector of pairwise correlation values
-  AD: cellSaveNormActDissim: within-group z-scored activation dissimilarity [2 (high/low HCAE) x 2 (detect/no-detect)] x n (animals/blocks) with every cell = vector of dissimilarity values
-  AD: cellSaveDissimCorrITD: inter-trial-distance dependence of assembly consistency [2 (high/low HCAE) x 2 (hits/misses) x n (animals/blocks)] with every cell = [p x 2] matrix with intertrial distance and correlation values for each trial pair
	- extra graph of above with means comparison of hits/misses high/low HCAE

%}