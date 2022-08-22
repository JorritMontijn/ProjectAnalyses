%function vecPredictorAssemblyZscores = doPlotActivityDistributions(cellContrastDataAD,sMetaDataAD,sPC,sParams)

%if ~exist('sParams','var'),sParams = struct;end
%if isfield(sParams,'intStartSes'),intStartSes=sParams.intStartSes;else intStartSes=1;end
%if isfield(sParams,'intStopSes'),intStopSes=sParams.intStopSes;else intStopSes=length(cellContrastDataAD);end
for intMouse=1
	close all
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
		load(['D:\Data\Results\stimdetection\dataRawPre_aggregate' strSes '.mat']);
	end
	
	%load separate ses files
	for intFile=1:numel(cellAggregate)
		%load data
		sLoad = load([cellAggregate{intFile} '.mat']);
		ses = sLoad.ses;

		%transform orientation to 
		ses.structStim.Orientation = mod(ses.structStim.Orientation,180);
		
		%assign to multi-cell
		cellSes{intFile} = ses;
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
	sParams.boolSaveData = true;
	strOldDir = cd(sParams.strFigDir);
	
	for intPopulation = vecBlockTypes
		%% GET RESPONSES AFTER COMBINING BLOCKS OF SAME POPULATION INTO ONE
		%take opposite directions as the same
		cellMultiSes{intPopulation}.structStim.Orientation = mod(cellMultiSes{intPopulation}.structStim.Orientation,180);
		vecOrientations = unique(cellMultiSes{intPopulation}.structStim.Orientation);
		intNeurons = numel(cellMultiSes{intPopulation}.neuron);
		
		%get neuronal responses per trial
		[matTrialResponse,cellSelectContrasts] = getTrialResponseData(cellMultiSes{intPopulation},cellMultiSes{intPopulation}.structStim);
		intContrasts = length(cellSelectContrasts);
		
		%get orientation-based trial selection vectors
		sTypesOri = getStimulusTypes(cellMultiSes{intPopulation},{'Orientation'});
		cellSelectOri = getSelectionVectors(cellMultiSes{intPopulation}.structStim,sTypesOri);
		
		%make normalized activity matrix per contrast
		matTrialNormResponse = zeros(size(matTrialResponse));
		for intContrast=1:intContrasts
			matTrialNormResponse(:,cellSelectContrasts{intContrast}) = imnorm(matTrialResponse(:,cellSelectContrasts{intContrast}));
		end
		
		%calculate preferred stimulus orientation for all neurons
		vecOriPrefContrasts = 4:6; % 8% - 100%
		matPref = nan(length(vecOriPrefContrasts),intNeurons);
		intCounter = 0;
		for intContrast = vecOriPrefContrasts
			intCounter = intCounter + 1;
			matResp = matTrialNormResponse(:,cellSelectContrasts{intContrast});
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
		vecNeuronPrefStim(isnan(vecNeuronPrefStim)) = [];
		intNeurons = numel(cellMultiSes{intPopulation}.neuron);
		intTrials = length(cellMultiSes{intPopulation}.structStim.Orientation);
		
		%get updated response matrices
		%get neuronal responses per trial
		[matTrialResponse,cellSelectContrasts] = getTrialResponseData(cellMultiSes{intPopulation},cellMultiSes{intPopulation}.structStim);
		
		%make normalized activity matrix per contrast
		matTrialNormResponse = zeros(size(matTrialResponse));
		for intContrast=1:intContrasts
			matTrialNormResponse(:,cellSelectContrasts{intContrast}) = imnorm(matTrialResponse(:,cellSelectContrasts{intContrast}));
		end
		
		
		%% calculate hit-correlated activity enhancement
		%{
	loop through contrasts, select preferred stimulus trials per neuron,
	calculate df/f activity for pref stim response for hits and for misses,
	calculate hit/miss difference per neuron per contrast per trial and
	plot binned matrices
	
	calculate correlation between 20% most hit-correlated activity enhanced
	neurons for different contrasts, make plot
	
	calculate overall HCAE for all neurons; vecStimDetectActInc
	select neurons with overall 20% highest HCAE as mean over different
	contrasts; indSelectHitCorrelatedNeurons and 20% lowest HCAE;
	indSelectHitAnticorrelatedNeurons
		%}
		
		%pre-allocate
		matHCAE = nan(intNeurons,length(cellSelectContrasts));
		cellHitResp = cell(6,1);
		cellMissResp = cell(6,1);
		indHitTrials = logical(cellMultiSes{intPopulation}.structStim.vecTrialResponse);
		for intContrast=1:intContrasts
			indSelectContrast = cellSelectContrasts{intContrast};
			intCounter=intContrast;
			intHitTrialCounter = 1;
			intMissTrialCounter = 1;
			cellHitResp{intContrast} = nan(1,intNeurons*intTrials);
			cellMissResp{intContrast} = nan(1,intNeurons*intTrials);
			%select pref stim trials per neuron
			for intNeuron=1:intNeurons
				intPrefStim = vecNeuronPrefStim(intNeuron);
				indPrefStimTrials = cellSelectOri{intPrefStim};
				indPrefHitTrials = indPrefStimTrials & indHitTrials & indSelectContrast;
				indPrefMissTrials = indPrefStimTrials & ~indHitTrials & indSelectContrast;
				
				%get data
				vecHitResp = matTrialResponse(intNeuron,indPrefHitTrials);
				vecMissResp = matTrialResponse(intNeuron,indPrefMissTrials);
				cellHitResp{intContrast}(intHitTrialCounter:(intHitTrialCounter+sum(indPrefHitTrials)-1)) = vecHitResp;
				cellMissResp{intContrast}(intMissTrialCounter:(intMissTrialCounter+sum(indPrefMissTrials)-1)) = vecMissResp;
				intHitTrialCounter = intHitTrialCounter + sum(indPrefHitTrials);
				intMissTrialCounter = intMissTrialCounter + sum(indPrefMissTrials);
				
				%get normalized response to calculate hit correlated activity enhancement (range [-1 1])
				vecHitNormResp = matTrialNormResponse(intNeuron,indPrefHitTrials);
				vecMissNormResp = matTrialNormResponse(intNeuron,indPrefMissTrials);
				dblMeanHitResp = mean(vecHitNormResp);if dblMeanHitResp < 0,dblMeanHitResp = 0;end
				dblMeanMissResp = mean(vecMissNormResp);if dblMeanMissResp < 0,dblMeanMissResp = 0;end
				matHCAE(intNeuron,intContrast) = (dblMeanHitResp-dblMeanMissResp)/(dblMeanHitResp+dblMeanMissResp);
			end
			%remove nans
			cellHitResp{intContrast}(isnan(cellHitResp{intContrast})) = [];
			cellMissResp{intContrast}(isnan(cellMissResp{intContrast})) = [];
		end
		%set nans to 0
		matHCAE(isnan(matHCAE)) = 0;
		vecStimDetectActInc = mean(matHCAE(:,2:5),2);
		
		%define correlated/uncorrelated neurons
		dblFrac = 1/3;
		intNumNeurons = round(intNeurons * dblFrac);
		[vecHCAE,vecNeuronsHCAE] = findmax(vecStimDetectActInc,intNumNeurons);
		[vecHUAE,vecNeuronsHAAE] = findmin(vecStimDetectActInc,intNumNeurons);
		indNeuronsHCAE = false(size(vecStimDetectActInc));
		indNeuronsHCAE(vecNeuronsHCAE) = true;
		indNeuronsHAAE = false(size(vecStimDetectActInc));
		indNeuronsHAAE(vecNeuronsHAAE) = true;
		
		indSelectHitCorrelatedNeurons = false(size(vecStimDetectActInc));
		indSelectHitCorrelatedNeurons(vecNeuronsHCAE) = true;
		indSelectHitUncorrelatedNeurons = false(size(vecStimDetectActInc));
		indSelectHitUncorrelatedNeurons(~indNeuronsHCAE & ~indNeuronsHAAE) = true;
		
		%get which neurons are most correlated for each contrast
		matMostCorrelated = false(size(matHCAE));
		for intContrast=1:intContrasts
			indMembers = false(1,intNeurons);
			[vecHCAE_Temp,vecNeuronsHCAE_Temp] = findmax(matHCAE(:,intContrast),intNumNeurons);
			indMembers(vecNeuronsHCAE_Temp) = true;
			matMostCorrelated(:,intContrast) = indMembers;
		end
		matHCN_Consistency = corr(matMostCorrelated);
		%set auto correlation to 0
		matHCN_Consistency(diag(diag(true(size(matHCN_Consistency))))) = 0;
		
		
		%% create hit-correlated activity region plot
		%{
	loop through contrasts, select preferred stimulus trials per neuron,
	calculate df/f activity for pref stim response for hits and for misses,
	calculate hit/miss difference per neuron per contrast per trial and
	plot binned matrices
	
	calculate correlation between 20% most hit-correlated activity enhanced
	neurons for different contrasts, make plot
	
	calculate overall HCAE for all neurons; vecStimDetectActInc
	select neurons with overall 20% highest HCAE as mean over different
	contrasts; indSelectHitCorrelatedNeurons and 20% lowest HCAE;
	indSelectHitAnticorrelatedNeurons
		%}
		
		%retrieve separate ses-based activity
		intRepRecs = length(cellSes);
		cellHitPrefResp = cell(intContrasts,intRepRecs);
		cellMissPrefResp = cell(intContrasts,intRepRecs);
		for intRepRec=1:intRepRecs
			
			%get data
			ses = cellSes{intRepRec};
			
			vecOrientations_sba = unique(ses.structStim.Orientation);
			intNeurons_sba = numel(ses.neuron);

			%get neuronal responses per trial
			[matTrialResponse_sba,cellSelectContrasts_sba] = getTrialResponseData(ses,ses.structStim);

			%get contrast & orientation-based trial selection vectors
			sTypesOri_sba = getStimulusTypes(ses,{'Orientation'});
			cellSelectOri_sba = getSelectionVectors(ses.structStim,sTypesOri_sba);

			sTypesC_sba = getStimulusTypes(ses,{'Contrast'});
			cellSelectC_sba = getSelectionVectors(ses.structStim,sTypesC_sba);
			
			%get preferred orientation for all neurons
			indHitTrials_sba = logical(ses.structStim.vecTrialResponse);
			
			%for which contrasts
			cellSelectTuning_sba = cell(size(cellSelectOri_sba));
			indSelectContrasts = cellSelectC_sba{4} | cellSelectC_sba{5} | cellSelectC_sba{6};
			for intOri=1:length(cellSelectOri_sba)
				cellSelectTuning_sba{intOri} = cellSelectOri_sba{intOri} & indSelectContrasts;
			end
			sTuning_sba = calcTuningRespMat(matTrialResponse_sba,cellSelectTuning_sba,vecOrientations);
			vecNeuronPrefStim_sba = sTuning_sba.vecPrefIndex;
			
			for intContrast=1:intContrasts
				indSelectContrast_sba = cellSelectC_sba{intContrast};
				intHitTrialCounter_sba = 1;
				intMissTrialCounter_sba = 1;
				cellHitPrefResp{intContrast,intRepRec} = nan(1,intNeurons*length(ses.structStim.Orientation));
				cellMissPrefResp{intContrast,intRepRec} = nan(1,intNeurons*length(ses.structStim.Orientation));
				for intNeuron=1:intNeurons_sba
					intPrefStim = vecNeuronPrefStim_sba(intNeuron);
					indPrefStimTrials = cellSelectOri_sba{intPrefStim};
					indPrefHitTrials = indPrefStimTrials & indHitTrials_sba & indSelectContrast_sba;
					indPrefMissTrials = indPrefStimTrials & ~indHitTrials_sba & indSelectContrast_sba;

					%get data
					vecHitResp = matTrialResponse_sba(intNeuron,indPrefHitTrials);
					vecMissResp = matTrialResponse_sba(intNeuron,indPrefMissTrials);
					cellHitPrefResp{intContrast,intRepRec}(intHitTrialCounter_sba:(intHitTrialCounter_sba+sum(indPrefHitTrials)-1)) = vecHitResp;
					cellMissPrefResp{intContrast,intRepRec}(intMissTrialCounter_sba:(intMissTrialCounter_sba+sum(indPrefMissTrials)-1)) = vecMissResp;
					intHitTrialCounter_sba = intHitTrialCounter_sba + sum(indPrefHitTrials);
					intMissTrialCounter_sba = intMissTrialCounter_sba + sum(indPrefMissTrials);
				end
				%remove trailing nans
				cellHitPrefResp{intContrast,intRepRec}(isnan(cellHitPrefResp{intContrast,intRepRec})) = [];
				cellMissPrefResp{intContrast,intRepRec}(isnan(cellMissPrefResp{intContrast,intRepRec})) = [];
			end
		end
		
		%data
		intSwitchZ = 1;
		if intSwitchZ
			vecBins = -1:0.1:2;
			vecTickY = [1 11 21 31];
			matHit = zeros(length(vecBins),intContrasts);
			matMiss = zeros(length(vecBins),intContrasts);
			for intRepRec = 1:intRepRecs
				for intContrast=1:intContrasts
					matHit(:,intContrast) = matHit(:,intContrast) + hist(zscore(cellHitPrefResp{intContrast,intRepRec}),vecBins)';
					matMiss(:,intContrast) = matMiss(:,intContrast) + hist(zscore(cellMissPrefResp{intContrast,intRepRec}),vecBins)';
				end
			end
			strLabelY = sprintf('Z-scored dF/F0 distribution');
		else
			vecBins = -0.05:0.01:0.25;
			vecTickY = [6 16 26];
			matHit = zeros(length(vecBins),intContrasts);
			matMiss = zeros(length(vecBins),intContrasts);
			for intRepRec = 1:intRepRecs
				for intContrast=1:intContrasts
					matHit(:,intContrast) = matHit(:,intContrast) + hist(zscore(cellHitPrefResp{intContrast,intRepRec}),vecBins)';
					matMiss(:,intContrast) = matMiss(:,intContrast) + hist(zscore(cellMissPrefResp{intContrast,intRepRec}),vecBins)';
				end
			end
			strLabelY = sprintf('Normalized dF/F0 distribution');
		end
		%if intMouse == 8 %only 2 hit resps at 100%, so remove from analysis
		%	matMiss(:,end) = mean(matHit(:,end));
		%end
		
		%normalize
		matHit = conv2(matHit,(1/3)*ones(3,1), 'same');
		matMiss = conv2(matMiss,(1/3)*ones(3,1), 'same');
		matHit = matHit./repmat(sum(matHit,1),[size(matHit,1) 1]);
		matMiss = matMiss./repmat(sum(matMiss,1),[size(matMiss,1) 1]);
		matDiff = matHit ./ (matHit + matMiss);
		matDiff = conv2(matDiff,(1/3)*ones(3,1), 'same');
		cellSaveMatrices{1,intPopulation} = matHit;
		cellSaveMatrices{2,intPopulation} = matMiss;
		cellSaveMatrices{3,intPopulation} = matDiff;
		
		%vars
		vecContrasts = [0 0.5 2 8 32 100];
		vecTickLabelY = vecBins(vecTickY);
		vecTickLabelX = vecContrasts;
		
		%plot over contrasts
		hActMat = figure;
		set(hActMat,'Color',[1 1 1]);
		figure(hActMat);
		
		%hits
		subplot(2,2,1)
		imagesc(matHit);
		axis xy;
		title('Pref / Hit; z=normalized count')
		set(gca,'YTick',vecTickY,'YTickLabel',vecTickLabelY)
		set(gca,'XTick',1:length(vecTickLabelX),'XTickLabel',vecTickLabelX)
		ylabel(strLabelY)
		xlabel('Contrast (%)')
		colormap(hot(128));colorbar;drawnow;freezeColors;cbfreeze;
		
		%misses
		subplot(2,2,2)
		imagesc(matMiss);
		axis xy;
		title('Pref / Miss; z=normalized count')
		set(gca,'YTick',vecTickY,'YTickLabel',vecTickLabelY)
		set(gca,'XTick',1:length(vecTickLabelX),'XTickLabel',vecTickLabelX)
		ylabel(strLabelY)
		xlabel('Contrast (%)')
		colormap(hot(128));colorbar;drawnow;freezeColors;cbfreeze;
		
		%difference
		subplot(2,2,3)
		imagesc(matDiff,[0 1]);
		axis xy;
		title('Diff')
		set(gca,'YTick',vecTickY,'YTickLabel',vecTickLabelY)
		set(gca,'XTick',1:length(vecTickLabelX),'XTickLabel',vecTickLabelX)
		ylabel(strLabelY)
		xlabel('Contrast (%)')
		colormap(redblue(128));colorbar;drawnow;freezeColors;cbfreeze;
		
		%select which bins to plot
		intTestType = 3;
		intBins = length(matDiff(:,1));
		matH = false(size(matDiff));
		matP = nan(size(matDiff));
		if intTestType < 3
			for intContrastIndex=1:size(matDiff,2)
				if intTestType == 1
					%bonferroni corrected ttest
					alpha=0.05/intBins;
					for intBin=1:length(matDiff(:,intContrastIndex))
						[h,p]=ttest(matDiff(:,intContrastIndex),matDiff(intBin,intContrastIndex),alpha);
						matP(intBin,intContrastIndex) = p;
						matH(intBin,intContrastIndex) = h;
					end
				elseif intTestType == 2
					%z-score bin count percentage, then take  >2sd's
					matP(:,intContrastIndex) = zscore(matDiff(:,intContrastIndex));
					matH(:,intContrastIndex) = matP(:,intContrastIndex)>1.5;
				end
			end
		elseif intTestType == 3
			matDiffT = matDiff;
			matDiffT(isnan(matDiffT)) = nanmean(matDiff(:));
			matP = reshape(zscore(matDiffT(:)),size(matDiff));
			matH = matP>2;
		end
		
		%select for calculation all bins that are significant on any contrast
		matPlotP = ones(size(matDiff));
		matPlotP(matH) = matP(matH);
		subplot(2,2,4)
		imagesc(matH)
		cellSaveMatrices{4,intPopulation} = matH;
		axis xy
		title('Significant outlier bins')
		set(gca,'YTick',vecTickY,'YTickLabel',vecTickLabelY)
		set(gca,'XTick',1:length(vecTickLabelX),'XTickLabel',vecTickLabelX)
		ylabel(strLabelY)
		xlabel('Contrast (%)')
		colormap(hot(2));colorbar;drawnow;freezeColors;cbfreeze;
		
		if sParams.boolSavePlots
			drawnow;
			jFig = get(handle(gcf), 'JavaFrame');
			jFig.setMaximized(true);
			figure(gcf);
			drawnow;
			strFig = sprintf('%sagg_detectcorrelated_heatmaps_pop%d_raw',strSes,intPopulation);
			export_fig([strFig '.tif']);
			export_fig([strFig '.pdf']);
		end
		
		%% plot highest hit-correlated neuron consistency
		figure
		imagesc(matHCN_Consistency,[-0.4 0.4]);
		colormap(redblue(128));colorbar;drawnow;freezeColors;cbfreeze;
		set(gca,'XTick',1:length(vecTickLabelX),'XTickLabel',vecTickLabelX)
		set(gca,'YTick',1:length(vecTickLabelX),'YTickLabel',vecTickLabelX)
		ylabel('Contrast (%)')
		xlabel('Contrast (%)')
		title(sprintf('Correlation of highest %.0f%% of hit-correlated neurons',dblFrac*100))
		
		%save data
		cellSaveNeuronConsistency{intPopulation} = matHCN_Consistency;
		
		if sParams.boolSavePlots
			drawnow;
			strFig = sprintf('%sagg_HCN_consistency_pop%d_raw',strSes,intPopulation);
			export_fig([strFig '.tif']);
			export_fig([strFig '.pdf']);
		end
		
		%% calculate noise + signal correlations
		[dummy,vecNeuronSortingIndex] = sort(vecStimDetectActInc,'descend');
		
		%get selection vectors&data
		sCorr = calcStimCorrs(cellMultiSes{intPopulation});
		intBorderSize = ceil(log10(intNeurons)+1);
		
		%signal correlations
		matReorderedSC = sCorr.matSignalCorrs(vecNeuronSortingIndex,vecNeuronSortingIndex);
		matPlotSC = [matReorderedSC; nan(intBorderSize,size(matReorderedSC,intPopulation))];
		matPlotSC = [matPlotSC nan(size(matPlotSC,1),intBorderSize)];
		hCorr = figure;
		subplot(2,2,1)
		imagesc(matPlotSC);
		nancolorbar(matPlotSC);
		
		%make borders
		vecReorderedNeuronType = 2*ones(size(indSelectHitCorrelatedNeurons));
		vecReorderedNeuronType(1:sum(indSelectHitCorrelatedNeurons)) = 0;
		vecReorderedNeuronType((sum(indSelectHitCorrelatedNeurons)+1):(sum(indSelectHitCorrelatedNeurons)+sum(indSelectHitUncorrelatedNeurons))) = 1;
		%vecReorderedNeuronType = [true(1,intNumNeurons) false(1,intNeurons-intNumNeurons)];
		hold on;
		cmap = colormap(jet(3));
		colormap(redblue(128));
		intStartNeuron = 1;
		intStopNeuron = intNumNeurons;
		fill([intStartNeuron intStartNeuron intStopNeuron intStopNeuron],[intNeurons+1 intNeurons+intBorderSize intNeurons+intBorderSize intNeurons+1],cmap(1,:),'EdgeColor',cmap(1,:));
		fill([intNeurons+1 intNeurons+intBorderSize intNeurons+intBorderSize intNeurons+1],[intStartNeuron intStartNeuron intStopNeuron intStopNeuron],cmap(1,:),'EdgeColor',cmap(1,:));
		intStartNeuron = intNumNeurons+1;
		intStopNeuron = intStartNeuron+intNumNeurons;
		fill([intStartNeuron intStartNeuron intStopNeuron intStopNeuron],[intNeurons+1 intNeurons+intBorderSize intNeurons+intBorderSize intNeurons+1],cmap(2,:),'EdgeColor',cmap(2,:));
		fill([intNeurons+1 intNeurons+intBorderSize intNeurons+intBorderSize intNeurons+1],[intStartNeuron intStartNeuron intStopNeuron intStopNeuron],cmap(2,:),'EdgeColor',cmap(2,:));
		intStartNeuron = intStopNeuron;
		intStopNeuron = intNeurons;
		fill([intStartNeuron intStartNeuron intStopNeuron intStopNeuron],[intNeurons+1 intNeurons+intBorderSize intNeurons+intBorderSize intNeurons+1],cmap(3,:),'EdgeColor',cmap(3,:));
		fill([intNeurons+1 intNeurons+intBorderSize intNeurons+intBorderSize intNeurons+1],[intStartNeuron intStartNeuron intStopNeuron intStopNeuron],cmap(3,:),'EdgeColor',cmap(3,:));
		hold off;
		title(sprintf('Signal correlation for block %d; blue=response correlated',intPopulation));
		
		xlabel('Neuron # [grouped by assembly membership]');
		ylabel('Neuron # [grouped by assembly membership]');
		
		%plot mean signal correlation within assembly, between assembly/non-assembly and within non-assembly
		subplot(2,2,2)
		[matMeanSC,matStdSC,cellValsSC] = getMatrixBlockMeans(matReorderedSC,vecReorderedNeuronType);
		vecY = [matMeanSC(1,1) matMeanSC(2,1) matMeanSC(2,2)];%within assembly/between/within non-assembly
		vecE = [matStdSC(1,1)/sqrt(length(cellValsSC{1,1})) matStdSC(2,1)/sqrt(length(cellValsSC{2,1})) matStdSC(2,2)/sqrt(length(cellValsSC{2,2}))];
		[h,dblP_AN] = ttest2(cellValsSC{1,1},cellValsSC{2,2}); %assembly/non-assembly
		[h,dblP_AB] = ttest2(cellValsSC{1,1},cellValsSC{2,1});%assembly/between
		[h,dblP_BN] = ttest2(cellValsSC{2,1},cellValsSC{2,2});%between/non-assembly
		errorbar(vecY,vecE,'ob','LineStyle','none');
		set(gca,'XTick',[1 2 3],'XTickLabel',{'Assembly','Between','Non-assembly'})
		ylabel('Mean Signal Correlation')
		title(sprintf('Mean +/- st err of signal correlation; A-N, p=%.3f; A-B, p=%.3f; B-N, p=%.3f',dblP_AN,dblP_AB,dblP_BN))
		
		%put in output
		cellSaveSignalCorrs{1,intPopulation} = cellValsSC{1,1};%signal correlations [3 (within-HCAR, between, within-non-HCAR) x n (animals)] with every cell = vector of pairwise correlation values
		cellSaveSignalCorrs{2,intPopulation} = cellValsSC{2,1};
		cellSaveSignalCorrs{3,intPopulation} = cellValsSC{2,2};
		
		%noise correlations
		subplot(2,2,3)
		matReorderedNC = sCorr.matNoiseCorrs(vecNeuronSortingIndex,vecNeuronSortingIndex);
		matPlotNC = [matReorderedNC; nan(intBorderSize,size(matReorderedNC,intPopulation))];
		matPlotNC = [matPlotNC nan(size(matPlotNC,1),intBorderSize)];
		imagesc(matPlotNC);
		nancolorbar(matPlotNC);
		
		%make borders
		hold on;
		cmap = colormap(jet(3));
		colormap(redblue(128));
		intStartNeuron = 1;
		intStopNeuron = intNumNeurons;
		fill([intStartNeuron intStartNeuron intStopNeuron intStopNeuron],[intNeurons+1 intNeurons+intBorderSize intNeurons+intBorderSize intNeurons+1],cmap(1,:),'EdgeColor',cmap(1,:));
		fill([intNeurons+1 intNeurons+intBorderSize intNeurons+intBorderSize intNeurons+1],[intStartNeuron intStartNeuron intStopNeuron intStopNeuron],cmap(1,:),'EdgeColor',cmap(1,:));
		intStartNeuron = intNumNeurons+1;
		intStopNeuron = intStartNeuron+intNumNeurons;
		fill([intStartNeuron intStartNeuron intStopNeuron intStopNeuron],[intNeurons+1 intNeurons+intBorderSize intNeurons+intBorderSize intNeurons+1],cmap(2,:),'EdgeColor',cmap(2,:));
		fill([intNeurons+1 intNeurons+intBorderSize intNeurons+intBorderSize intNeurons+1],[intStartNeuron intStartNeuron intStopNeuron intStopNeuron],cmap(2,:),'EdgeColor',cmap(2,:));
		intStartNeuron = intStopNeuron;
		intStopNeuron = intNeurons;
		fill([intStartNeuron intStartNeuron intStopNeuron intStopNeuron],[intNeurons+1 intNeurons+intBorderSize intNeurons+intBorderSize intNeurons+1],cmap(3,:),'EdgeColor',cmap(3,:));
		fill([intNeurons+1 intNeurons+intBorderSize intNeurons+intBorderSize intNeurons+1],[intStartNeuron intStartNeuron intStopNeuron intStopNeuron],cmap(3,:),'EdgeColor',cmap(3,:));
		hold off;
		title(sprintf('Noise correlation for block %d; blue=response correlated',intPopulation));
		
		xlabel('Neuron # [grouped by assembly membership]');
		ylabel('Neuron # [grouped by assembly membership]');
		
		%plot mean noise correlation within assembly, between assembly/non-assembly and within non-assembly
		subplot(2,2,4)
		[matMeanNC,matStdNC,cellValsNC] = getMatrixBlockMeans(matReorderedNC,vecReorderedNeuronType);
		vecY = [matMeanNC(1,1) matMeanNC(2,1) matMeanNC(2,2)];%within assembly/between/within non-assembly
		vecE = [matStdNC(1,1)/sqrt(length(cellValsNC{1,1})) matStdNC(2,1)/sqrt(length(cellValsNC{2,1})) matStdNC(2,2)/sqrt(length(cellValsNC{2,2}))];
		[h,dblP_AN] = ttest2(cellValsNC{1,1},cellValsNC{2,2}); %assembly/non-assembly
		[h,dblP_AB] = ttest2(cellValsNC{1,1},cellValsNC{2,1});%assembly/between
		[h,dblP_BN] = ttest2(cellValsNC{2,1},cellValsNC{2,2});%between/non-assembly
		errorbar(vecY,vecE,'ob','LineStyle','none');
		set(gca,'XTick',[1 2 3],'XTickLabel',{'Assembly','Between','Non-assembly'})
		ylabel('Mean Noise Correlation')
		title(sprintf('Mean +/- st err of noise correlation; A-N, p=%.3f; A-B, p=%.3f; B-N, p=%.3f',dblP_AN,dblP_AB,dblP_BN))
		
		%put in output
		cellSaveNoiseCorrs{1,intPopulation} = cellValsNC{1,1};%signal correlations [3 (within-HCAR, between, within-non-HCAR) x n (animals)] with every cell = vector of pairwise correlation values
		cellSaveNoiseCorrs{2,intPopulation} = cellValsNC{2,1};
		cellSaveNoiseCorrs{3,intPopulation} = cellValsNC{2,2};
		
		if sParams.boolSavePlots
			drawnow;
			jFig = get(handle(gcf), 'JavaFrame');
			jFig.setMaximized(true);
			figure(gcf);
			drawnow;
			strFig = sprintf('%sagg_HCAEneurons_correlations_pop%d_raw',strSes,intPopulation);
			export_fig([strFig '.tif']);
			export_fig([strFig '.pdf']);
		end
		
		
		%% calculation based on HCAE neurons [hit-correlated activity enhanced]
		%group per pref stim
		vecReorder = [];
		vecReorderedPrefStim = [];
		for intStim=unique(vecNeuronPrefStim)
			%neurons
			vecSubPop = find(vecNeuronPrefStim==intStim);
			[vecVals,vecSubSort] = sort(vecStimDetectActInc(vecSubPop),'descend');
			vecSortOrig = vecSubPop(vecSubSort);
			vecReorder = [vecReorder vecSortOrig];
			vecReorderedPrefStim = [vecReorderedPrefStim ones(1,length(vecSubPop))*intStim];
		end
		
		%loop through contrasts & stim types for trial reordering
		vecTrialOriVector = [];
		vecTrialContrastVector = [];
		
		for intStim=unique(vecNeuronPrefStim)
			for intContrastIndex=1:length(cellSelectContrasts)
				%trials
				indSelectOri = cellSelectOri{intStim};
				vecTrialOriVector = [vecTrialOriVector find(indSelectOri)];
				
				vecSelectContrastTrials = cellSelectContrasts{intContrastIndex} & indSelectOri;
				vecTrialContrastVector = [vecTrialContrastVector find(vecSelectContrastTrials)];
			end
		end
		
		%normalize per contrast
		matRespNormPerContrast = nan(size(matTrialResponse));
		for intContrastIndex=1:length(cellSelectContrasts)
			vecSelectContrastTrials = cellSelectContrasts{intContrastIndex};
			matRespNormPerContrast(:,vecSelectContrastTrials) = zscore(matTrialResponse(:,vecSelectContrastTrials),[],2);
		end
		
		%sort neurons by detection correlated activity enhancement
		[dummy,vecSortNeurons] = sort(vecStimDetectActInc,'descend');
		
		%plot
		figure
		colormap(hot(256))
		%imagesc(matRespNormPerContrast)
		imagesc(matRespNormPerContrast(vecSortNeurons,vecTrialContrastVector))
		set(gca,'XTick',1:sum(vecSelectContrastTrials):length(vecTrialContrastVector))
		colorbar
		title('Z-scored activation level per neuron per trial')
		ylabel('Neuron number sorted by HCAE')
		xlabel('Trial number, grouped by contrast, sorted by stimulus ori')
		%save plot
		if sParams.boolSavePlots
			drawnow;
			strFig = sprintf('%sagg_detectcorrelated_actdissimilarity_pop%d_raw',strSes,intPopulation);
			export_fig([strFig '.tif']);
			export_fig([strFig '.pdf']);
		end
		
		%high/low DCAE
		%indSelectHitCorrelatedNeurons = vecStimDetectActInc>(mean(vecStimDetectActInc)+std(vecStimDetectActInc));
		%indSelectHitUncorrelatedNeurons = vecStimDetectActInc<(mean(vecStimDetectActInc)-std(vecStimDetectActInc));
		%indSelectHitCorrelatedNeurons = vecStimDetectActInc>(mean(vecStimDetectActInc));
		%indSelectHitUncorrelatedNeurons = vecStimDetectActInc<(mean(vecStimDetectActInc));
		intNeuronsHigh = sum(indSelectHitCorrelatedNeurons);
		intNeuronsLow = sum(indSelectHitUncorrelatedNeurons);
		intValuesHigh = (intNeuronsHigh*intNeuronsHigh-intNeuronsHigh)/2;
		intValuesLow = (intNeuronsLow*intNeuronsLow-intNeuronsLow)/2;
		intValuesTotal = (intNeurons*intNeurons-intNeurons)/2;
		matHeterogeneity = nan(intTrials,intValuesTotal);
		matRespDistNoRemHitCorrelated = nan(intTrials,intValuesHigh);
		matRespDistNoRemHitUncorrelated = nan(intTrials,intValuesLow);
		
		matRespDistLoZRemHitCorrelated = nan(intTrials,intValuesHigh);
		matRespDistLoZRemHitUncorrelated = nan(intTrials,intValuesLow);
		
		matRespDistHiZRemHitCorrelated = nan(intTrials,intValuesHigh);
		matRespDistHiZRemHitUncorrelated = nan(intTrials,intValuesLow);
		
		matRespDistLoARemHitCorrelated = nan(intTrials,intValuesHigh);
		matRespDistLoARemHitUncorrelated = nan(intTrials,intValuesLow);
		
		matRespDistHiARemHitCorrelated = nan(intTrials,intValuesHigh);
		matRespDistHiARemHitUncorrelated = nan(intTrials,intValuesLow);
		
		matTrialIndexHigh = nan(intTrials,intValuesHigh);
		matTrialIndexLow = nan(intTrials,intValuesLow);
		%calculate inverse correlation matrix per trial
		for intTrial=1:intTrials
			%general
			matTrialIndexHigh(intTrial,:) = intTrial;
			matTrialIndexLow(intTrial,:) = intTrial;
			
			%% STANDARD
			%perform calculation for all neurons
			vecActivity = matRespNormPerContrast(:,intTrial);
			matZ1 = repmat(vecActivity,[1 intNeurons]);
			matZ2 = repmat(vecActivity',[intNeurons 1]);
			matDistAll = abs(matZ1 - matZ2);
			
			%do calculation for high response correlated neurons
			vecActivity = matRespNormPerContrast(indSelectHitCorrelatedNeurons,intTrial);
			matZ1 = repmat(vecActivity,[1 intNeuronsHigh]);
			matZ2 = repmat(vecActivity',[intNeuronsHigh 1]);
			matDistHigh = abs(matZ1 - matZ2);
			
			%do calculation for low response correlated neurons
			vecActivity = matRespNormPerContrast(indSelectHitUncorrelatedNeurons,intTrial);
			matZ1 = repmat(vecActivity,[1 intNeuronsLow]);
			matZ2 = repmat(vecActivity',[intNeuronsLow 1]);
			matDistLow = abs(matZ1 - matZ2);
			
			%save data as vectors in matrix
			matSelectAll = tril(true(size(matDistAll)),-1);
			matSelectHigh = tril(true(size(matDistHigh)),-1);
			matSelectLow = tril(true(size(matDistLow)),-1);
			vecRespDistHighHCAE = matDistHigh(matSelectHigh);
			vecRespDistLowHCAE = matDistLow(matSelectLow);
			vecRespDistAll = matDistAll(matSelectAll);
			
			matRespDistNoRemHitCorrelated(intTrial,:) = vecRespDistHighHCAE;
			matRespDistNoRemHitUncorrelated(intTrial,:) = vecRespDistLowHCAE;
			matHeterogeneity(intTrial,:) = vecRespDistAll;
			
			%% REMOVE MOST ACTIVE NEURONS (z-score)
			%do calculation for high response correlated neurons
			vecActivity = matRespNormPerContrast(indSelectHitCorrelatedNeurons,intTrial);
			vecActivity = findmin(vecActivity,round(length(vecActivity)*0.8))'; %remove 20% most active cells
			matZ1 = repmat(vecActivity,[1 length(vecActivity)]);
			matZ2 = repmat(vecActivity',[length(vecActivity) 1]);
			matDistHigh = abs(matZ1 - matZ2);
			
			%do calculation for low response correlated neurons
			vecActivity = matRespNormPerContrast(indSelectHitUncorrelatedNeurons,intTrial);
			vecActivity = findmin(vecActivity,round(length(vecActivity)*0.8))'; %remove 20% most active cells
			matZ1 = repmat(vecActivity,[1 length(vecActivity)]);
			matZ2 = repmat(vecActivity',[length(vecActivity) 1]);
			matDistLow = abs(matZ1 - matZ2);
			
			%save data as vectors in matrix
			matSelectHigh = tril(true(size(matDistHigh)),-1);
			matSelectLow = tril(true(size(matDistLow)),-1);
			vecRespDistHighHCAE = matDistHigh(matSelectHigh);
			vecRespDistLowHCAE = matDistLow(matSelectLow);
			
			matRespDistHiZRemHitCorrelated(intTrial,:) = [vecRespDistHighHCAE;nan(intValuesHigh-length(vecRespDistHighHCAE),1)];
			matRespDistHiZRemHitUncorrelated(intTrial,:) = [vecRespDistLowHCAE;nan(intValuesLow-length(vecRespDistLowHCAE),1)];
			
			%% REMOVE LEAST ACTIVE NEURONS (z-score)
			%do calculation for high response correlated neurons
			vecActivity = matRespNormPerContrast(indSelectHitCorrelatedNeurons,intTrial);
			vecActivity = findmax(vecActivity,round(length(vecActivity)*0.8))'; %remove 20% least active cells
			matZ1 = repmat(vecActivity,[1 length(vecActivity)]);
			matZ2 = repmat(vecActivity',[length(vecActivity) 1]);
			matDistHigh = abs(matZ1 - matZ2);
			
			%do calculation for low response correlated neurons
			vecActivity = matRespNormPerContrast(indSelectHitUncorrelatedNeurons,intTrial);
			vecActivity = findmax(vecActivity,round(length(vecActivity)*0.8))'; %remove 20% least active cells
			matZ1 = repmat(vecActivity,[1 length(vecActivity)]);
			matZ2 = repmat(vecActivity',[length(vecActivity) 1]);
			matDistLow = abs(matZ1 - matZ2);
			
			%save data as vectors in matrix
			matSelectHigh = tril(true(size(matDistHigh)),-1);
			matSelectLow = tril(true(size(matDistLow)),-1);
			vecRespDistHighHCAE = matDistHigh(matSelectHigh);
			vecRespDistLowHCAE = matDistLow(matSelectLow);
			
			matRespDistLoZRemHitCorrelated(intTrial,:) = [vecRespDistHighHCAE;nan(intValuesHigh-length(vecRespDistHighHCAE),1)];
			matRespDistLoZRemHitUncorrelated(intTrial,:) = [vecRespDistLowHCAE;nan(intValuesLow-length(vecRespDistLowHCAE),1)];
			
			%% REMOVE MOST ACTIVE NEURONS (df/f)
			%do calculation for high response correlated neurons
			vecRawActivity = matTrialResponse(indSelectHitCorrelatedNeurons,intTrial);
			vecActivity = matRespNormPerContrast(indSelectHitCorrelatedNeurons,intTrial);
			[dummy,vecIndex] = findmin(vecRawActivity,round(length(vecActivity)*0.8)); %remove 20% most active cells
			vecActivity = vecActivity(vecIndex);
			matZ1 = repmat(vecActivity,[1 length(vecActivity)]);
			matZ2 = repmat(vecActivity',[length(vecActivity) 1]);
			matDistHigh = abs(matZ1 - matZ2);
			
			%do calculation for low response correlated neurons
			vecRawActivity = matTrialResponse(indSelectHitUncorrelatedNeurons,intTrial);
			vecActivity = matRespNormPerContrast(indSelectHitUncorrelatedNeurons,intTrial);
			[dummy,vecIndex] = findmin(vecRawActivity,round(length(vecActivity)*0.8)); %remove 20% most active cells
			vecActivity = vecActivity(vecIndex);
			matZ1 = repmat(vecActivity,[1 length(vecActivity)]);
			matZ2 = repmat(vecActivity',[length(vecActivity) 1]);
			matDistLow = abs(matZ1 - matZ2);
			
			%save data as vectors in matrix
			matSelectHigh = tril(true(size(matDistHigh)),-1);
			matSelectLow = tril(true(size(matDistLow)),-1);
			vecRespDistHighHCAE = matDistHigh(matSelectHigh);
			vecRespDistLowHCAE = matDistLow(matSelectLow);
			
			matRespDistHiARemHitCorrelated(intTrial,:) = [vecRespDistHighHCAE;nan(intValuesHigh-length(vecRespDistHighHCAE),1)];
			matRespDistHiARemHitUncorrelated(intTrial,:) = [vecRespDistLowHCAE;nan(intValuesLow-length(vecRespDistLowHCAE),1)];
			
			%% REMOVE LEAST ACTIVE NEURONS (df/f)
			%do calculation for high response correlated neurons
			vecRawActivity = matTrialResponse(indSelectHitCorrelatedNeurons,intTrial);
			vecActivity = matRespNormPerContrast(indSelectHitCorrelatedNeurons,intTrial);
			[dummy,vecIndex] = findmax(vecRawActivity,round(length(vecActivity)*0.8)); %remove 20% most active cells
			vecActivity = vecActivity(vecIndex);
			matZ1 = repmat(vecActivity,[1 length(vecActivity)]);
			matZ2 = repmat(vecActivity',[length(vecActivity) 1]);
			matDistHigh = abs(matZ1 - matZ2);
			
			%do calculation for low response correlated neurons
			vecRawActivity = matTrialResponse(indSelectHitUncorrelatedNeurons,intTrial);
			vecActivity = matRespNormPerContrast(indSelectHitUncorrelatedNeurons,intTrial);
			[dummy,vecIndex] = findmax(vecRawActivity,round(length(vecActivity)*0.8)); %remove 20% most active cells
			vecActivity = vecActivity(vecIndex);
			matZ1 = repmat(vecActivity,[1 length(vecActivity)]);
			matZ2 = repmat(vecActivity',[length(vecActivity) 1]);
			matDistLow = abs(matZ1 - matZ2);
			
			%save data as vectors in matrix
			matSelectHigh = tril(true(size(matDistHigh)),-1);
			matSelectLow = tril(true(size(matDistLow)),-1);
			vecRespDistHighHCAE = matDistHigh(matSelectHigh);
			vecRespDistLowHCAE = matDistLow(matSelectLow);
			
			matRespDistLoARemHitCorrelated(intTrial,:) = [vecRespDistHighHCAE;nan(intValuesHigh-length(vecRespDistHighHCAE),1)];
			matRespDistLoARemHitUncorrelated(intTrial,:) = [vecRespDistLowHCAE;nan(intValuesLow-length(vecRespDistLowHCAE),1)];
			
		end
		
		for intRemType=1:5
			if intRemType == 1
				%no rem
				matRespDistHitCorrelated = matRespDistNoRemHitCorrelated;
				matRespDistHitUncorrelated = matRespDistNoRemHitUncorrelated;
				strRemType = 'No';
			elseif intRemType == 2
				%rem hi
				indKeep = ~isnan(matRespDistHiZRemHitCorrelated(1,:));
				matRespDistHiZRemHitCorrelated = matRespDistHiZRemHitCorrelated(:,indKeep);
				matRespDistHiZRemHitUncorrelated = matRespDistHiZRemHitUncorrelated(:,indKeep);
				matRespDistHitCorrelated = matRespDistHiZRemHitCorrelated;
				matRespDistHitUncorrelated = matRespDistHiZRemHitUncorrelated;
				strRemType = 'High Z';
			elseif intRemType == 3
				%rem lo
				indKeep = ~isnan(matRespDistLoZRemHitCorrelated(1,:));
				matRespDistLoZRemHitCorrelated = matRespDistLoZRemHitCorrelated(:,indKeep);
				matRespDistLoZRemHitUncorrelated = matRespDistLoZRemHitUncorrelated(:,indKeep);
				matRespDistHitCorrelated = matRespDistLoZRemHitCorrelated;
				matRespDistHitUncorrelated = matRespDistLoZRemHitUncorrelated;
				strRemType = 'Low Z';
			elseif intRemType == 4
				%rem hi
				indKeep = ~isnan(matRespDistHiARemHitCorrelated(1,:));
				matRespDistHiARemHitCorrelated = matRespDistHiARemHitCorrelated(:,indKeep);
				matRespDistHiARemHitUncorrelated = matRespDistHiARemHitUncorrelated(:,indKeep);
				matRespDistHitCorrelated = matRespDistHiARemHitCorrelated;
				matRespDistHitUncorrelated = matRespDistHiARemHitUncorrelated;
				strRemType = 'High A';
			elseif intRemType == 5
				%rem lo
				indKeep = ~isnan(matRespDistLoARemHitCorrelated(1,:));
				matRespDistLoARemHitCorrelated = matRespDistLoARemHitCorrelated(:,indKeep);
				matRespDistLoARemHitUncorrelated = matRespDistLoARemHitUncorrelated(:,indKeep);
				matRespDistHitCorrelated = matRespDistLoARemHitCorrelated;
				matRespDistHitUncorrelated = matRespDistLoARemHitUncorrelated;
				strRemType = 'Low A';
			end
			
			%{
		for intTrial=1:intTrials
			%do calculation for high response correlated neurons
			vecActivity = matRespNormPerContrast(:,intTrial);
			matZ1 = repmat(vecActivity,[1 intNeurons]);
			matZ2 = repmat(vecActivity',[intNeurons 1]);
			matDist = abs(matZ1 - matZ2);

			figure
			imagesc(matDist)
			colormap(flipud(jet))
			colorbar
			xlabel('Neuron ID')
			ylabel('Neuron ID')
			title(sprintf('Neuronal pairwise distance in z-scored activation level\nExample trial (# %d)',intTrial))
			pause
			close
		end
			%}
			
			%plot activity dissimilarity over trials
			figure
			subplot(2,2,1)
			%hit correlated neurons
			vecX = 1:intTrials;
			vecInvX = intTrials:-1:1;
			vecMeanH = nanmean(matRespDistHitCorrelated,2)';
			vecErrH = nanstd(matRespDistHitCorrelated,[],2)'/sqrt(intValuesHigh);
			vecMinTrace = vecMeanH-vecErrH;
			vecMaxTrace = vecMeanH+vecErrH;
			
			fill([vecX vecInvX],[vecMaxTrace vecMinTrace(vecInvX)],[1.0 0.7 0.7],'EdgeColor','none');
			hold on
			plot(vecX,vecMeanH,'-','LineWidth',2,'Color',[1 0 0]);
			
			%hit uncorrelated neurons
			vecMeanL = nanmean(matRespDistHitUncorrelated,2)';
			vecErrL = nanstd(matRespDistHitUncorrelated,[],2)'/sqrt(intValuesLow);
			vecMinTrace = vecMeanL-vecErrL;
			vecMaxTrace = vecMeanL+vecErrL;
			
			fill([vecX vecInvX],[vecMinTrace vecMaxTrace(vecInvX)],[0.7 0.7 1],'EdgeColor','none');
			plot(vecX,vecMeanL,'-','LineWidth',2,'Color',[0 0 1]);
			hold off
			[h,dblP] = ttest2(matRespDistHitCorrelated(:),matRespDistHitUncorrelated(:));
			title(sprintf('%s rem; High HCAE (red): mean=%.3f; Low HCAE (blue): mean=%.3f; ttest p=%.3f',strRemType,mean(matRespDistHitCorrelated(:)),mean(matRespDistHitUncorrelated(:)),dblP))
			ylabel(sprintf('Mean population activity dissimilarity'))
			xlabel('Trial number')
			ylim([0 2])
			
			%split data set
			indSelectContrasts = cellMultiSes{intPopulation}.structStim.Contrast == 0.005 | cellMultiSes{intPopulation}.structStim.Contrast == 0.02;
			indSelectResp = cellMultiSes{intPopulation}.structStim.vecTrialResponse;
			
			%hit trials only
			subplot(2,2,3)
			indSelectRespTrials = indSelectContrasts & indSelectResp;
			%hit correlated neurons
			intRespTrials = sum(indSelectRespTrials);
			matRespDistHitcorrelatedSub = matRespDistHitCorrelated(indSelectRespTrials,:);
			vecX = 1:intRespTrials;
			vecInvX = intRespTrials:-1:1;
			vecMeanRespH = nanmean(matRespDistHitcorrelatedSub,2)';
			vecErrRespH = nanstd(matRespDistHitcorrelatedSub,[],2)'/sqrt(intValuesHigh);
			vecMinTraceRespH = vecMeanRespH-vecErrRespH;
			vecMaxTraceRespH = vecMeanRespH+vecErrRespH;
			
			fill([vecX vecInvX],[vecMaxTraceRespH vecMinTraceRespH(vecInvX)],[1.0 0.7 0.7],'EdgeColor','none');
			hold on
			plot(vecX,vecMeanRespH,'-','LineWidth',2,'Color',[1 0 0]);
			
			%hit uncorrelated neurons
			matRespDistHituncorrelatedSub = matRespDistHitUncorrelated(indSelectRespTrials,:);
			vecX = 1:intRespTrials;
			vecInvX = intRespTrials:-1:1;
			vecMeanRespL = nanmean(matRespDistHituncorrelatedSub,2)';
			vecErrRespL = nanstd(matRespDistHituncorrelatedSub,[],2)'/sqrt(intValuesLow);
			vecMinTraceRespL = vecMeanRespL-vecErrRespL;
			vecMaxTraceRespL = vecMeanRespL+vecErrRespL;
			
			fill([vecX vecInvX],[vecMaxTraceRespL vecMinTraceRespL(vecInvX)],[0.7 0.7 1],'EdgeColor','none');
			plot(vecX,vecMeanRespL,'-','LineWidth',2,'Color',[0 0 1]);
			hold off
			[h,dblP] = ttest2(matRespDistHituncorrelatedSub(:),matRespDistHitcorrelatedSub(:));
			dblMeanRH = nanmean(matRespDistHitcorrelatedSub(:));
			dblSdRH = nanstd(matRespDistHitcorrelatedSub(:));
			intValsRH = length(matRespDistHitcorrelatedSub(:));
			dblMeanRL = nanmean(matRespDistHituncorrelatedSub(:));
			dblSdRL = nanstd(matRespDistHituncorrelatedSub(:));
			intValsRL = length(matRespDistHituncorrelatedSub(:));
			title(sprintf('%s rem; Resp; High HCAE (red): mean=%.3f; Low HCAE (blue): mean=%.3f; ttest p=%.3f',strRemType,dblMeanRH,dblMeanRL,dblP))
			ylabel(sprintf('Mean population activity dissimilarity'))
			xlabel('Trial number')
			ylim([0 2])
			
			%miss trials only
			indSelectNoRespTrials = indSelectContrasts & ~indSelectResp;
			subplot(2,2,4)
			%hit correlated neurons
			intNoRespTrials = sum(indSelectNoRespTrials);
			matNoRespHigh = matRespDistHitCorrelated(indSelectNoRespTrials,:);
			vecX = 1:intNoRespTrials;
			vecInvX = intNoRespTrials:-1:1;
			vecMeanNoRespH = nanmean(matNoRespHigh,2)';
			vecErrNoRespH = nanstd(matNoRespHigh,[],2)'/sqrt(intValuesHigh);
			vecMinTraceNoRespH = vecMeanNoRespH-vecErrNoRespH;
			vecMaxTraceNoRespH = vecMeanNoRespH+vecErrNoRespH;
			
			fill([vecX vecInvX],[vecMaxTraceNoRespH vecMinTraceNoRespH(vecInvX)],[1.0 0.7 0.7],'EdgeColor','none');
			hold on
			plot(vecX,vecMeanNoRespH,'-','LineWidth',2,'Color',[1 0 0]);
			
			%hit uncorrelated neurons
			matNoRespLow = matRespDistHitUncorrelated(indSelectNoRespTrials,:);
			vecX = 1:intNoRespTrials;
			vecInvX = intNoRespTrials:-1:1;
			vecMeanNoRespL = nanmean(matNoRespLow,2)';
			vecErrNoRespL = nanstd(matNoRespLow,[],2)'/sqrt(intValuesLow);
			vecMinTraceNoRespL = vecMeanNoRespL-vecErrNoRespL;
			vecMaxTraceNoRespL = vecMeanNoRespL+vecErrNoRespL;
			
			fill([vecX vecInvX],[vecMaxTraceNoRespL vecMinTraceNoRespL(vecInvX)],[0.7 0.7 1],'EdgeColor','none');
			plot(vecX,vecMeanNoRespL,'-','LineWidth',2,'Color',[0 0 1]);
			hold off
			[h,dblP] = ttest2(matNoRespLow(:),matNoRespHigh(:));
			dblMeanNRH = nanmean(matNoRespHigh(:));
			dblSdNRH = nanstd(matNoRespHigh(:));
			intValsNRH = length(matNoRespHigh(:));
			dblMeanNRL = nanmean(matNoRespLow(:));
			dblSdNRL = nanstd(matNoRespLow(:));
			intValsNRL = length(matNoRespLow(:));
			
			title(sprintf('%s rem; No Resp; High HCAE (red): mean=%.3f; Low HCAE (blue): mean=%.3f; ttest p=%.3f',strRemType,dblMeanNRH,dblMeanNRL,dblP))
			ylabel(sprintf('Mean population activity dissimilarity'))
			xlabel('Trial number')
			ylim([0 2])
			
			%comparison of hit/miss and hit corr/uncorr
			subplot(2,2,2)
			dblErrNRH = dblSdNRH/sqrt(intValsNRH);
			dblErrNRL = dblSdNRL/sqrt(intValsNRL);
			dblErrRH = dblSdRH/sqrt(intValsRH);
			dblErrRL = dblSdRL/sqrt(intValsRL);
			errorbar(1:4,[dblMeanNRH dblMeanNRL dblMeanRH dblMeanRL],[dblErrNRH dblErrNRL dblErrRH dblErrRL],'Linestyle','none','Marker','x');
			set(gca,'XTick',1:4,'XTickLabel',{'No Resp High','No Resp Low','Resp High','Resp Low'})
			xlim([0.5 4.5])
			ylabel('Mean within-group z-scored activation dissimilarity')
			title(sprintf('%s rem; Block %d',strRemType,intPopulation));
			
			%put in output
			cellSaveNormActDissim{1,1,intPopulation,intRemType} = matRespDistHitcorrelatedSub(:); %within-group z-scored activation dissimilarity [2 (high/low HCAE) x 2 (detect/no-detect)] x n (animals) with every cell = vector of dissimilarity values
			cellSaveNormActDissim{1,2,intPopulation,intRemType} = matNoRespHigh(:);
			cellSaveNormActDissim{2,1,intPopulation,intRemType} = matRespDistHituncorrelatedSub(:);
			cellSaveNormActDissim{2,2,intPopulation,intRemType} = matNoRespLow(:);
			
			if sParams.boolSavePlots
				drawnow;
				jFig = get(handle(gcf), 'JavaFrame');
				jFig.setMaximized(true);
				figure(gcf);
				drawnow;
				strFig = sprintf('%sagg_detectcorrelated_actdissimilarity_pop%d_%srem_raw',strSes,intPopulation,strRemType);
				export_fig([strFig '.tif']);
				export_fig([strFig '.pdf']);
			end
		end
		%no rem
		matRespDistHitCorrelated = matRespDistNoRemHitCorrelated;
		matRespDistHitUncorrelated = matRespDistNoRemHitUncorrelated;
		strRemType = 'No';
		
		%% calculate the face of god
		%group per pref stim
		[vecPrefStimC,vecReorderC] = sort(vecNeuronPrefStim(indSelectHitCorrelatedNeurons),'ascend');
		[vecPrefStimU,vecReorderU] = sort(vecNeuronPrefStim(indSelectHitUncorrelatedNeurons),'ascend');
		
		%loop through contrasts & stim types for trial reordering
		vecTrialOriVector = [];
		vecTrialContrastVector = [];
		
		intStimType = 0;
		vecStimType = zeros(1,length(cellSelectOri{1}));
		for intStim=unique(vecNeuronPrefStim)
			for intContrastIndex=2:length(cellSelectContrasts)
				intStimType = intStimType + 1;
				%trials
				indSelectTrials = cellSelectOri{intStim} & cellSelectContrasts{intContrastIndex};
				vecStimType(indSelectTrials) = intStimType;
			end
		end
		
		intCorrType = 2;
		if intCorrType == 1 %all neurons
			%calculate trial-by-trial correlations of population activity
			matCorrC=corr(matRespNormPerContrast(indSelectHitCorrelatedNeurons,:));
			matCorrU=corr(matRespNormPerContrast(indSelectHitUncorrelatedNeurons,:));
		else %only pref-pop
			for intPrefType=1:4
				indPrefNeurons = vecNeuronPrefStim' == intPrefType;
				indSelectC = indSelectHitCorrelatedNeurons & indPrefNeurons;
				indSelectU = indSelectHitUncorrelatedNeurons & indPrefNeurons;
				
				%calculate trial-by-trial correlations of population activity
				if sum(indSelectC) > 0,cellCorrC{intPrefType}=corr(matRespNormPerContrast(indSelectC,:));else cellCorrC{intPrefType} = nan(intTrials,intTrials);end
				if sum(indSelectU) > 0,cellCorrU{intPrefType}=corr(matRespNormPerContrast(indSelectU,:));else cellCorrU{intPrefType} = nan(intTrials,intTrials);end
			end
			matCorrU = nan(size(cellCorrU{1}));
			matCorrC = nan(size(cellCorrC{1}));
			
			for intStimType=vecStimType
				intPrefType = mod(intStimType,4);
				if intPrefType == 0,intPrefType = 4;end
				
				%hit
				indSelectHitTrials = vecStimType == intStimType & indSelectResp;
				
				%get trials
				matCorrU(indSelectHitTrials,indSelectHitTrials) = cellCorrU{intPrefType}(indSelectHitTrials,indSelectHitTrials);
				matCorrC(indSelectHitTrials,indSelectHitTrials) = cellCorrC{intPrefType}(indSelectHitTrials,indSelectHitTrials);
				
				%miss
				indSelectMissTrials = vecStimType == intStimType & ~indSelectResp;
				
				%get trials
				matCorrU(indSelectMissTrials,indSelectMissTrials) = cellCorrU{intPrefType}(indSelectMissTrials,indSelectMissTrials);
				matCorrC(indSelectMissTrials,indSelectMissTrials) = cellCorrC{intPrefType}(indSelectMissTrials,indSelectMissTrials);
				
				%hit-miss
				%get trials
				matCorrU(indSelectMissTrials,indSelectHitTrials) = cellCorrU{intPrefType}(indSelectMissTrials,indSelectHitTrials);
				matCorrU(indSelectHitTrials,indSelectMissTrials) = cellCorrU{intPrefType}(indSelectHitTrials,indSelectMissTrials);
				matCorrC(indSelectMissTrials,indSelectHitTrials) = cellCorrC{intPrefType}(indSelectMissTrials,indSelectHitTrials);
				matCorrC(indSelectHitTrials,indSelectMissTrials) = cellCorrC{intPrefType}(indSelectHitTrials,indSelectMissTrials);
			end
		end
		
		
		%compute absolute trial distance
		matTrialDist = abs(repmat(1:intTrials,[intTrials 1])-repmat((1:intTrials)',[1 intTrials]));
		%matSelect = tril(true(size(matTrialDist)),-1);
		
		%separate across-trial correlation correlations for different trial types (1-4)
		%also split for hit/miss trials
		%select only 0.5 & 2%
		intHitTrials = sum(indSelectResp);
		intMissTrials = sum(~indSelectResp);
		vecDistHit = [];
		vecDistMiss = [];
		vecCorrHitC = [];
		vecCorrHitU = [];
		vecCorrMissC = [];
		vecCorrMissU = [];
		vecDistBetween = [];
		vecCorrBetweenC = [];
		vecCorrBetweenU = [];
		for intStimType=vecStimType
			%% hit; select trials
			indSelectHitTrials = vecStimType == intStimType & indSelectResp;
			if sum(indSelectHitTrials) > 2
				%get trials
				matDistSubHit = matTrialDist(indSelectHitTrials,indSelectHitTrials);
				matCorrSubHitLow = matCorrU(indSelectHitTrials,indSelectHitTrials);
				matCorrSubHitHigh = matCorrC(indSelectHitTrials,indSelectHitTrials);
				
				%subselection matrix hit
				matSubSelectHit = tril(true(size(matDistSubHit)),-1);
				
				%subselect
				vecDistHit = [vecDistHit matDistSubHit(matSubSelectHit)'];
				vecCorrHitU = [vecCorrHitU matCorrSubHitLow(matSubSelectHit)'];
				vecCorrHitC = [vecCorrHitC matCorrSubHitHigh(matSubSelectHit)'];
			end
			
			%% miss; select trials
			indSelectMissTrials = vecStimType == intStimType & ~indSelectResp;
			if sum(indSelectMissTrials) > 2
				%get trials
				matDistSubMiss = matTrialDist(indSelectMissTrials,indSelectMissTrials);
				matCorrSubMissLow = matCorrU(indSelectMissTrials,indSelectMissTrials);
				matCorrSubMissHigh = matCorrC(indSelectMissTrials,indSelectMissTrials);
				
				%subselection matrix hit
				matSubSelectMiss = tril(true(size(matDistSubMiss)),-1);
				
				%subselect
				vecDistMiss = [vecDistMiss matDistSubMiss(matSubSelectMiss)'];
				vecCorrMissU = [vecCorrMissU matCorrSubMissLow(matSubSelectMiss)'];
				vecCorrMissC = [vecCorrMissC matCorrSubMissHigh(matSubSelectMiss)'];
			end
			
			%% hit-miss; select trials
			if sum(indSelectMissTrials | indSelectHitTrials) > 3
				%get trials
				matDistSubBetween = matTrialDist(indSelectMissTrials,indSelectHitTrials);
				matCorrSubBetweenLow = matCorrU(indSelectMissTrials,indSelectHitTrials);
				matCorrSubBetweenHigh = matCorrC(indSelectMissTrials,indSelectHitTrials);
				
				%subselection matrix hit
				matSubSelectBetween = tril(true(size(matDistSubBetween)),-1);
				
				%subselect
				vecDistBetween = [vecDistBetween matDistSubBetween(matSubSelectBetween)'];
				vecCorrBetweenU = [vecCorrBetweenU matCorrSubBetweenLow(matSubSelectBetween)'];
				vecCorrBetweenC = [vecCorrBetweenC matCorrSubBetweenHigh(matSubSelectBetween)'];
			end
		end
		
		%plot
		figure
		subplot(2,2,1)
		imagesc(matCorrU);
		colorbar;
		title(sprintf('Correlation between trials of z-scored activation distance (Low HCAE; block %d)',intPopulation))
		xlabel('Trial')
		ylabel('Trial')
		
		subplot(2,2,2)
		imagesc(matCorrC);
		colorbar;
		title(sprintf('Correlation between trials of z-scored activation distance (High HCAE; block %d)',intPopulation))
		xlabel('Trial')
		ylabel('Trial')
		
		subplot(2,2,3)
		%high HCAE
		intStep = 25;
		intMax = ceil(intTrials/intStep)*intStep;
		vecBinX = 0:intStep:intMax;
		vecPlotX = (intStep/2):intStep:(intMax-intStep/2);
		[nVec,meanVec,stdVec,cellVals,cellIDs] = makeBins(vecDistHit,vecCorrHitC,vecBinX);
		errorbar(vecPlotX,meanVec,stdVec./nVec, 'Color', 'r')
		sStatsH=regstats(vecCorrHitC,vecDistHit,'linear',{'beta','rsquare','tstat'});
		vecX = get(gca,'XLim');
		vecY = polyval(sStatsH.beta([2 1]),vecX);
		hold on
		plot(vecX,vecY,'r')
		
		%low HCAE
		vecBinX = 0:intStep:intMax;
		vecPlotX = (intStep/2):intStep:(intMax-intStep/2);
		[nVec,meanVec,stdVec,cellVals,cellIDs] = makeBins(vecDistHit,vecCorrHitU,vecBinX);
		errorbar(vecPlotX,meanVec,stdVec./nVec, 'Color', 'b')
		sStatsL=regstats(vecCorrHitU,vecDistHit,'linear',{'beta','rsquare','tstat'});
		vecX = get(gca,'XLim');
		vecY = polyval(sStatsL.beta([2 1]),vecX);
		plot(vecX,vecY,'b')
		hold off
		title(sprintf('Hits; Lin reg high HCAE (red): slope=%.3f; p=%.3f; intercept=%.3f; p=%.3f; R^2=%.3f\nLin reg low HCAE (blue): slope=%.3f; p=%.3f; intercept=%.3f; p=%.3f; R^2=%.3f',...
			sStatsH.beta(2),sStatsH.tstat.pval(2),sStatsH.beta(1),sStatsH.tstat.pval(1),sStatsH.rsquare,sStatsL.beta(2),sStatsL.tstat.pval(2),sStatsL.beta(1),sStatsL.tstat.pval(1),sStatsL.rsquare))
		ylabel(sprintf('Assembly consistency over time\n(Inter-trial correlation of mean population activity)'))
		xlabel('Inter-trial distance')
		
		subplot(2,2,4)
		%high HCAE
		vecBinX = 0:intStep:intMax;
		vecPlotX = (intStep/2):intStep:(intMax-intStep/2);
		[nVec,meanVec,stdVec,cellVals,cellIDs] = makeBins(vecDistMiss,vecCorrMissC,vecBinX);
		errorbar(vecPlotX,meanVec,stdVec./nVec, 'Color', 'm')
		sStatsH=regstats(vecCorrMissC,vecDistMiss,'linear',{'beta','rsquare','tstat'});
		vecX = get(gca,'XLim');
		vecY = polyval(sStatsH.beta([2 1]),vecX);
		hold on
		plot(vecX,vecY,'m')
		
		%low HCAE
		vecBinX = 0:intStep:intMax;
		vecPlotX = (intStep/2):intStep:(intMax-intStep/2);
		[nVec,meanVec,stdVec,cellVals,cellIDs] = makeBins(vecDistMiss,vecCorrMissU,vecBinX);
		errorbar(vecPlotX,meanVec,stdVec./nVec, 'Color', 'c')
		sStatsL=regstats(vecCorrMissU,vecDistMiss,'linear',{'beta','rsquare','tstat'});
		vecX = get(gca,'XLim');
		vecY = polyval(sStatsL.beta([2 1]),vecX);
		plot(vecX,vecY,'c')
		hold off
		title(sprintf('Misses; Lin reg high HCAE (magenta): slope=%.3f; p=%.3f; intercept=%.3f; p=%.3f; R^2=%.3f\nLin reg low HCAE (cyan): slope=%.3f; p=%.3f; intercept=%.3f; p=%.3f; R^2=%.3f',...
			sStatsH.beta(2),sStatsH.tstat.pval(2),sStatsH.beta(1),sStatsH.tstat.pval(1),sStatsH.rsquare,sStatsL.beta(2),sStatsL.tstat.pval(2),sStatsL.beta(1),sStatsL.tstat.pval(1),sStatsL.rsquare))
		ylabel(sprintf('Assembly consistency over time\n(Inter-trial correlation of mean population activity)'))
		xlabel('Inter-trial distance')
		
		%%%% OUTPUT
		cellSaveCorrITD{1,1,intPopulation} = vecDistHit;
		cellSaveCorrITD{1,2,intPopulation} = vecCorrHitC;
		cellSaveCorrITD{1,3,intPopulation} = vecCorrHitU;
		cellSaveCorrITD{2,1,intPopulation} = vecDistMiss;
		cellSaveCorrITD{2,2,intPopulation} = vecCorrMissC;
		cellSaveCorrITD{2,3,intPopulation} = vecCorrMissU;
		cellSaveCorrITD{3,1,intPopulation} = vecDistBetween;
		cellSaveCorrITD{3,2,intPopulation} = vecCorrBetweenC;
		cellSaveCorrITD{3,3,intPopulation} = vecCorrBetweenU;
		
		%plot
		if sParams.boolSavePlots
			drawnow;
			jFig = get(handle(gcf), 'JavaFrame');
			jFig.setMaximized(true);
			figure(gcf);
			drawnow;
			strFig = sprintf('%sagg_HCAEneurons_dissimilaritycorrelations_pop%d_%srem_raw',strSes,intPopulation,strRemType);
			export_fig([strFig '.tif']);
			export_fig([strFig '.pdf']);
		end
		
		
		%comparison of means
		figure
		dblErrCCHH = nanstd(vecCorrHitC)/sqrt(length(vecCorrHitC));
		dblErrCCHL = nanstd(vecCorrHitU)/sqrt(length(vecCorrHitU));
		dblErrCCMH = nanstd(vecCorrMissC)/sqrt(length(vecCorrMissC));
		dblErrCCML = nanstd(vecCorrMissU)/sqrt(length(vecCorrMissU));
		dblErrCCBH = nanstd(vecCorrBetweenC)/sqrt(length(vecCorrBetweenC));
		dblErrCCBL = nanstd(vecCorrBetweenU)/sqrt(length(vecCorrBetweenU));
		errorbar(1:6,[nanmean(vecCorrMissC) nanmean(vecCorrMissU) nanmean(vecCorrHitC) nanmean(vecCorrHitU) nanmean(vecCorrBetweenC) nanmean(vecCorrBetweenU)],[dblErrCCMH dblErrCCML dblErrCCHH dblErrCCHL dblErrCCBH dblErrCCBL],'Linestyle','none','Marker','x');
		set(gca,'XTick',1:6,'XTickLabel',{'Miss HCN','Miss HUN','Hit HCN','Hit HUN','Between HCN','Between HUN'})
		xlim([0.5 6.5])
		ylabel('Mean inter-trial correlation of population activity')
		title(sprintf('Block %d',intPopulation));
		
		if sParams.boolSavePlots
			drawnow;
			strFig = sprintf('%sagg_HCAEneurons_intertrialcorr_pop%d_%srem_raw',strSes,intPopulation,strRemType);
			export_fig([strFig '.tif']);
			export_fig([strFig '.pdf']);
		end
		
		%plot debiased activity distance increases
		figure
		subplot(1,2,1)
		dblMeanNRHd = dblMeanNRH / dblMeanNRH; %debiased
		dblErrNRHd = dblErrNRH / dblMeanNRH;
		dblMeanRHd = dblMeanRH / dblMeanNRH;
		dblErrRHd = dblErrRH/ dblMeanNRH;
		dblMeanNRLd = dblMeanNRL / dblMeanNRL;
		dblErrNRLd = dblErrNRL / dblMeanNRL;
		dblMeanRLd = dblMeanRL / dblMeanNRL;
		dblErrRLd = dblErrRL/ dblMeanNRL;
		
		errorbar(1:4,[dblMeanNRHd dblMeanNRLd dblMeanRHd dblMeanRLd],[dblErrNRHd dblErrNRLd dblErrRHd dblErrRLd],'Linestyle','none','Marker','x');
		set(gca,'XTick',1:4,'XTickLabel',{'No Resp High','No Resp Low','Resp High','Resp Low'})
		xlim([0.5 4.5])
		ylabel('Mean debiased within-group z-scored activation dissimilarity')
		title(sprintf('%s rem; Block %d',strRemType,intPopulation))
		%relative changes
		
		subplot(1,2,2)
		vecHighResp = matRespDistHitcorrelatedSub(:);
		vecLowResp = matRespDistHituncorrelatedSub(:);
		vecHighRespInv = ((vecHighResp - dblMeanNRH) / dblMeanNRH)*100;
		vecLowRespInv = ((vecLowResp - dblMeanNRL) / dblMeanNRL)*100;
		dblMeanHRI = mean(vecHighRespInv);
		dblErrHRI = std(vecHighRespInv) / sqrt(length(vecHighRespInv));
		dblMeanLRI = mean(vecLowRespInv);
		dblErrLRI = std(vecLowRespInv) / sqrt(length(vecLowRespInv));
		
		errorbar(1:2,[dblMeanHRI dblMeanLRI],[dblErrHRI dblErrLRI],'Linestyle','none','Marker','x');
		set(gca,'XTick',1:2,'XTickLabel',{'High HCAE neurons','Low HCAE neurons'})
		xlim([0.5 2.5])
		ylabel('% increase of within-group z-scored distances during detection trials')
		if sParams.boolSavePlots
			drawnow;
			jFig = get(handle(gcf), 'JavaFrame');
			jFig.setMaximized(true);
			figure(gcf);
			drawnow;
			strFig = sprintf('%sagg_debiased_dissimilarity_pop%d_raw',strSes,intPopulation);
			export_fig([strFig '.tif']);
			export_fig([strFig '.pdf']);
		end
		
		%plot dependency on reaction times
		vecRTs = cellMultiSes{intPopulation}.structStim.vecTrialRespSecs(indSelectRespTrials) - cellMultiSes{intPopulation}.structStim.SecsOn(indSelectRespTrials); %get RTs for hit trials and selected contrasts
		[vecRTsSorted,vecRTSortIndex] = sort(vecRTs,'ascend');
		matHitcorrRTSorted = matRespDistHitcorrelatedSub(vecRTSortIndex,:);
		matHituncorrRTSorted = matRespDistHituncorrelatedSub(vecRTSortIndex,:);
		vecHitcorr = nanmean(matHitcorrRTSorted,2);
		vecHituncorr = nanmean(matHituncorrRTSorted,2);
		
		%Z activity
		matActZHCAE = matRespNormPerContrast(indSelectHitCorrelatedNeurons,indSelectRespTrials)';
		matActZNonHCAE = matRespNormPerContrast(indSelectHitUncorrelatedNeurons,indSelectRespTrials)';
		
		matZHCAERTSorted = matActZHCAE(vecRTSortIndex,:);
		matZNonHCAERTSorted = matActZNonHCAE(vecRTSortIndex,:);
		vecZHCAE = mean(matZHCAERTSorted,2);
		vecZNonHCAE = mean(matZNonHCAERTSorted,2);
		
		%dF/F
		matActHCAE = matTrialResponse(indSelectHitCorrelatedNeurons,indSelectRespTrials)';
		matActNonHCAE = matTrialResponse(indSelectHitUncorrelatedNeurons,indSelectRespTrials)';
		
		matAHCAERTSorted = matActHCAE(vecRTSortIndex,:);
		matANonHCAERTSorted = matActNonHCAE(vecRTSortIndex,:);
		vecAHCAE = mean(matAHCAERTSorted,2);
		vecANonHCAE = mean(matANonHCAERTSorted,2);
		
		%save
		cellSaveRTDependency{1,intPopulation} = vecRTsSorted;
		cellSaveRTDependency{2,intPopulation} = vecHitcorr;
		cellSaveRTDependency{3,intPopulation} = vecHituncorr;
		cellSaveRTDependency{4,intPopulation} = vecZHCAE;
		cellSaveRTDependency{5,intPopulation} = vecZNonHCAE;
		cellSaveRTDependency{6,intPopulation} = vecAHCAE;
		cellSaveRTDependency{7,intPopulation} = vecANonHCAE;
		
		%plot
		figure
		subplot(2,2,1);
		scatter(vecRTsSorted,vecHitcorr,'r')
		hold on
		scatter(vecRTsSorted,vecHituncorr,'b')
		drawnow
		
		%perform regressions
		sStatsC=regstats(vecHitcorr,vecRTsSorted,'linear',{'beta','rsquare','tstat'});
		vecX = get(gca,'XLim');
		vecY = polyval(sStatsC.beta([2 1]),vecX);
		plot(vecX,vecY,'r')
		
		sStatsU=regstats(vecHituncorr,vecRTsSorted,'linear',{'beta','rsquare','tstat'});
		vecY = polyval(sStatsU.beta([2 1]),vecX);
		
		%plot reaction time vs activity dissimilarity correlation
		plot(vecX,vecY,'b')
		hold off
		
		title(sprintf('Pop heterogeneity; Lin reg high HCAE (blue): slope=%.3f; p=%.3f; intercept=%.3f; p=%.3f; R^2=%.3f\nLin reg low HCAE (red): slope=%.3f; p=%.3f; intercept=%.3f; p=%.3f; R^2=%.3f',...
			sStatsC.beta(2),sStatsC.tstat.pval(2),sStatsC.beta(1),sStatsC.tstat.pval(1),sStatsC.rsquare,sStatsU.beta(2),sStatsU.tstat.pval(2),sStatsU.beta(1),sStatsU.tstat.pval(1),sStatsU.rsquare))
		
		xlabel('Reaction Time (s)')
		ylabel('Mean population activity dissimilarity')
		
		%z-scored activity
		subplot(2,2,2);
		scatter(vecRTsSorted,vecZHCAE,'r')
		hold on
		scatter(vecRTsSorted,vecZNonHCAE,'b')
		drawnow
		
		%perform regressions
		sStatsC=regstats(vecZHCAE,vecRTsSorted,'linear',{'beta','rsquare','tstat'});
		vecX = get(gca,'XLim');
		vecY = polyval(sStatsC.beta([2 1]),vecX);
		plot(vecX,vecY,'r')
		
		sStatsU=regstats(vecZNonHCAE,vecRTsSorted,'linear',{'beta','rsquare','tstat'});
		vecY = polyval(sStatsU.beta([2 1]),vecX);
		
		%plot reaction time vs activity dissimilarity correlation
		plot(vecX,vecY,'b')
		hold off
		
		title(sprintf('Z-scored act; Lin reg high HCAE (blue): slope=%.3f; p=%.3f; intercept=%.3f; p=%.3f; R^2=%.3f\nLin reg low HCAE (red): slope=%.3f; p=%.3f; intercept=%.3f; p=%.3f; R^2=%.3f',...
			sStatsC.beta(2),sStatsC.tstat.pval(2),sStatsC.beta(1),sStatsC.tstat.pval(1),sStatsC.rsquare,sStatsU.beta(2),sStatsU.tstat.pval(2),sStatsU.beta(1),sStatsU.tstat.pval(1),sStatsU.rsquare))
		
		xlabel('Reaction Time (s)')
		ylabel('Mean z-scored population activity')
		
		%dF/F activity
		subplot(2,2,3);
		scatter(vecRTsSorted,vecAHCAE,'r')
		hold on
		scatter(vecRTsSorted,vecANonHCAE,'b')
		drawnow
		
		%perform regressions
		sStatsC=regstats(vecAHCAE,vecRTsSorted,'linear',{'beta','rsquare','tstat'});
		vecX = get(gca,'XLim');
		vecY = polyval(sStatsC.beta([2 1]),vecX);
		plot(vecX,vecY,'r')
		
		sStatsU=regstats(vecANonHCAE,vecRTsSorted,'linear',{'beta','rsquare','tstat'});
		vecY = polyval(sStatsU.beta([2 1]),vecX);
		
		%plot reaction time vs activity dissimilarity correlation
		plot(vecX,vecY,'b')
		hold off
		
		title(sprintf('Z-scored act; Lin reg HCAE (red): slope=%.3f; p=%.3f; intercept=%.3f; p=%.3f; R^2=%.3f\nLin reg non-HCAE (blue): slope=%.3f; p=%.3f; intercept=%.3f; p=%.3f; R^2=%.3f',...
			sStatsC.beta(2),sStatsC.tstat.pval(2),sStatsC.beta(1),sStatsC.tstat.pval(1),sStatsC.rsquare,sStatsU.beta(2),sStatsU.tstat.pval(2),sStatsU.beta(1),sStatsU.tstat.pval(1),sStatsU.rsquare))
		
		xlabel('Reaction Time (s)')
		ylabel('Mean dF/F population activity')
		
		if sParams.boolSavePlots
			drawnow;
			jFig = get(handle(gcf), 'JavaFrame');
			jFig.setMaximized(true);
			figure(gcf);
			drawnow;
			strFig = sprintf('%sagg_behaviorallycorrelated_act_z_dissimilarity_pop%d_raw',strSes,intPopulation);
			export_fig([strFig '.tif']);
			export_fig([strFig '.pdf']);
		end
		
		%% heterogeneity for different contrasts, split for hit/miss (like dF/F)
		%figure
		hHetCon = figure;
		set(hHetCon,'Color',[1 1 1]);
		figure(hHetCon);
		
		%pre-allocate
		vecContrasts = sMetaDataAD{1}.vecContrasts(2:end)*100;
		intContrasts = length(vecContrasts);
		vecP = nan(1,intContrasts);
		vecMetaHitY = nan(1,intContrasts);
		vecMetaHitE = nan(1,intContrasts);
		vecMetaMissY = nan(1,intContrasts);
		vecMetaMissE = nan(1,intContrasts);
		intContrastCounter = 0;
		for intContrastIndex=1:intContrasts
			%get contrast
			dblContrast = vecContrasts(intContrastIndex);
			if dblContrast == 0,continue;end
			intContrastCounter = intContrastCounter + 1;
			
			%get resp
			vecHeteroHit = matHeterogeneity(cellSelectContrasts{intContrastCounter+1} & indSelectResp,:);
			vecHeteroHit = vecHeteroHit(:);
			vecHeteroMiss = matHeterogeneity(cellSelectContrasts{intContrastCounter+1} & ~indSelectResp,:);
			vecHeteroMiss = vecHeteroMiss(:);
			
			%perform analyses per contrast
			[h,p,ci] = ttest2(vecHeteroHit,vecHeteroMiss);
			vecY = [mean(vecHeteroHit) mean(vecHeteroMiss)];
			
			%put in output
			cellSaveMatContHetero{intPopulation}(intContrastCounter,1) = vecY(1);
			cellSaveMatContHetero{intPopulation}(intContrastCounter,2) = vecY(2);
			
			%put in meta vector
			vecP(intContrastIndex) = p;
			vecMetaHitY(intContrastIndex) = vecY(1);
			vecMetaHitE(intContrastIndex) = std(vecHeteroHit)/sqrt(length(vecHeteroHit));
			vecMetaMissY(intContrastIndex) = vecY(2);
			vecMetaMissE(intContrastIndex) = std(vecHeteroMiss)/sqrt(length(vecHeteroMiss));
		end
		
		%pre-compute variables
		vecWindow = [1 5];
		vecWindowSelect = vecWindow(1):vecWindow(end);
		intWL = length(vecWindowSelect);
		vecWindowInv = intWL:-1:1;
		vecWindowPlotInv = vecWindow(end):-1:vecWindow(1);
		vecX = [vecContrasts(vecWindowSelect) vecContrasts(vecWindowPlotInv)];
		vecLineX = vecContrasts(vecWindowSelect);
		
		%get data
		for intResp=[0 1]
			if intResp == 1
				vecMeanTrace = vecMetaHitY(vecWindowSelect);
				vecSE = vecMetaHitE(vecWindowSelect);
				vecColorFill = [0.7 1 0.7];
				vecColorLine = [0 1 0];
			else
				vecColorLine = [1 0 0];
				vecColorFill = [1 0.7 0.7];
				vecMeanTrace = vecMetaMissY(vecWindowSelect);
				vecSE = vecMetaMissE(vecWindowSelect);
			end
			
			vecMinTrace = vecMeanTrace-vecSE;
			vecMaxTrace = vecMeanTrace+vecSE;
			vecY = [vecMinTrace vecMaxTrace(vecWindowInv)];
			
			%plot
			hold on
			fill(vecX,vecY,vecColorFill,'EdgeColor','none');
			plot(vecLineX,vecMeanTrace,'-','LineWidth',2,'Color',vecColorLine);
			hold off
		end
		set(gca,'XScale','log','YScale','linear')
		set(gca,'XTick',vecLineX,'XTickLabel',vecLineX)
		title(sprintf('Pop het [%d]; p-vals: 0.5: %.3f; 2: %.3f; 8: %.3f; 32: %.3f; 100: %.3f',intPopulation,vecP))
		grid on
		xlabel('Contrast')
		ylabel('Population response heterogeneity')
		xlim(vecContrasts(vecWindow))
		%ylim([-0.01 0.06])
		legend({'SEM','Miss','SEM','Hit'},'Location','Best')
		
		drawnow;
		
		if sParams.boolSavePlots
			drawnow;
			strFig = sprintf('%s_heterogeneity_over_contrasts_pop%d_raw',strSes,intPopulation);
			export_fig([strFig '.tif']);
			export_fig([strFig '.pdf']);
		end
		
		
		
		%% dF/F over time + heterogeneity over time; also split for fast/slow trials
		sParamsHet.intWindowLength = 1;
		[vecHeterogeneity,vecActivity_dFoF] = calcSlidingHeteroGen(cellMultiSes{intPopulation},sParamsHet);
		
		%get stim data
		cellFieldsC = {'Contrast'};
		sTypesC = getStimulusTypes(cellMultiSes{intPopulation},cellFieldsC);
		cellSelectC = getSelectionVectors(cellMultiSes{intPopulation}.structStim,sTypesC);
		vecC = sTypesC.matTypes;
		
		%get resp data
		indMiss = cellMultiSes{intPopulation}.structStim.vecTrialResponse == 0;
		indFast = cellMultiSes{intPopulation}.structStim.vecTrialRespSecs < nanmean(cellMultiSes{intPopulation}.structStim.vecTrialRespSecs);
		indSlow = cellMultiSes{intPopulation}.structStim.vecTrialRespSecs >= nanmean(cellMultiSes{intPopulation}.structStim.vecTrialRespSecs);
		
		%general
		vecWindowSecs = [-3 5];
		vecWindow = round(vecWindowSecs*cellMultiSes{intPopulation}.samplingFreq);
		vecLineX = (vecWindow(1):vecWindow(end))/cellMultiSes{intPopulation}.samplingFreq;
		vecWindowInv = length(vecLineX):-1:1;
		vecX = [vecLineX vecLineX(vecWindowInv)];
		
		%plot
		cellSaveActTime = {};
		cellSaveHetTime = {};
		hHetTime = figure;
		hActTime = figure;
		for intC=1:length(vecC)
			indTrials = cellSelectC{intC};
			for intAct=[0 1]
				if intAct == 0,figure(hHetTime);else figure(hActTime);end
				subplot(2,3,intC);
				for intType=1:3
					if intType == 1 %miss
						vecTrials = find(indTrials&indMiss);
						strType = 'Miss';
						vecColorLine = [1 0 0];
						vecColorFill = [1 0.7 0.7];
					elseif intType == 2 %slow
						vecTrials = find(indTrials&indSlow);
						strType = 'Slow';
						vecColorLine = [1 1 0];
						vecColorFill = [1 1 0.7];
					else %fast
						vecTrials = find(indTrials&indFast);
						strType = 'Fast';
						vecColorLine = [0 1 0];
						vecColorFill = [0.7 1 0.7];
					end
					matRawData = zeros(length(vecTrials),length(vecWindow(1):vecWindow(end)));
					
					%get data
					vecStimOn = cellMultiSes{intPopulation}.structStim.FrameOn(vecTrials);
					for intTrial=1:length(vecTrials)
						intStart = vecStimOn(intTrial)+vecWindow(1);
						intStop = vecStimOn(intTrial)+vecWindow(end);
						if intAct == 0,
							matRawData(intTrial,:) = vecHeterogeneity(intStart:intStop);
						else
							matRawData(intTrial,:) = vecActivity_dFoF(intStart:intStop);
						end
						
					end
					cellType{intType} = matRawData;
					
					vecMeanTrace = mean(matRawData,1);
					vecSE = std(matRawData,[],1);
					vecMinTrace = vecMeanTrace-vecSE;
					vecMaxTrace = vecMeanTrace+vecSE;
					vecY = [vecMinTrace vecMaxTrace(vecWindowInv)];
					
					%plot
					hold on
					%fill(vecX,vecY,vecColorFill,'EdgeColor','none');
					plot(vecLineX,vecMeanTrace,'-','LineWidth',2,'Color',vecColorLine);
					hold off
				end
				
				%save data&labels
				if intAct == 0,
					cellSaveHetTime{intC,intPopulation} = cellType;
					title(sprintf('Pop Het [%d]; Contrast %.1f',intPopulation,vecC(intC)*100))
					ylabel('Population response heterogeneity')
				else
					cellSaveActTime{intC,intPopulation} = cellType;
					title(sprintf('Pop Act [%d]; Contrast %.1f',intPopulation,vecC(intC)*100))
					ylabel('Population dF/F')
				end
				grid on
				xlabel('Time after stim onset (s)')
				xlim(vecWindowSecs)
				%ylim([-0.01 0.06])
				drawnow;
			end
		end
		
		if sParams.boolSavePlots
			figure(hHetTime)
			drawnow;
			strFig = sprintf('%s_heterogeneity_over_time_pop%d_raw',strSes,intPopulation);
			export_fig([strFig '.tif']);
			export_fig([strFig '.pdf']);
			
			figure(hActTime)
			drawnow;
			strFig = sprintf('%s_activity_over_time_pop%d_raw',strSes,intPopulation);
			export_fig([strFig '.tif']);
			export_fig([strFig '.pdf']);
		end
		
		%% decode if trial will be fast, slow or miss based on dF/F and heterogeneity
		%{
	per animal, sum miss, slow and fast trials; trake min(m,s,f)-1 trials for
	likelihood construction and take one other random trial per condition to
	decode [make procedure work over all time points]
	=> then bootstrap procedure
		%}
		
		%get data
		if ~exist('vecHeterogeneity','var') || isempty(vecHeterogeneity)
			sParamsHet.intWindowLength = 1;%round(cellMultiSes{1}.samplingFreq);
			[vecHeterogeneity,vecActivity_dFoF] = calcSlidingHeteroGen(cellMultiSes{intPopulation},sParamsHet);
		end
		
		%get resp data
		indMiss = isnan(cellMultiSes{intPopulation}.structStim.vecTrialRespSecs);
		indFast = cellMultiSes{intPopulation}.structStim.vecTrialRespSecs < nanmean(cellMultiSes{intPopulation}.structStim.vecTrialRespSecs);
		indSlow = cellMultiSes{intPopulation}.structStim.vecTrialRespSecs >= nanmean(cellMultiSes{intPopulation}.structStim.vecTrialRespSecs);
		
		%general
		vecWindowSecs = [-3 5];
		vecWindow = round(vecWindowSecs*cellMultiSes{intPopulation}.samplingFreq);
		vecLineX = (vecWindow(1):vecWindow(end))/cellMultiSes{intPopulation}.samplingFreq;
		vecSelectSecs = [-3 0];
		vecSelect = round(vecSelectSecs*cellMultiSes{intPopulation}.samplingFreq);
		intStart = find(vecLineX>=vecSelectSecs(1),1);
		intStop = find(vecLineX>=vecSelectSecs(2),1);
		intFramesPerStim = length(vecWindow(1):vecWindow(end));
		vecSelectFrames = intStart:intStop;
		
		%save vars
		cellSaveActDecodeTime = {};
		cellSaveHetDecodeTime = {};
		hDecRespTypeTime = figure;
		
		%get nr of trials
		intMiss = sum(indMiss);
		vecMiss = find(indMiss);
		intFast = sum(indFast);
		vecFast = find(indFast);
		intSlow = sum(indSlow);
		vecSlow = find(indSlow);
		intNrTrials = min([intMiss intFast intSlow]);
		
		%transform data to easy-access format
		%pre-allocate
		matRawDataAct = zeros(length(indMiss),intFramesPerStim);
		matRawDataHet = zeros(length(indMiss),intFramesPerStim);
		
		%get data
		vecStimOn = cellMultiSes{intPopulation}.structStim.FrameOn;
		for intTrial=1:length(indMiss)
			intStart = vecStimOn(intTrial)+vecWindow(1);
			intStop = vecStimOn(intTrial)+vecWindow(end);
			
			matRawDataHet(intTrial,:) = vecHeterogeneity(intStart:intStop);
			matRawDataAct(intTrial,:) = vecActivity_dFoF(intStart:intStop);
		end
		
		
		%bootstrap decoding
		intIters = 1000;
		
		%pre-allocate output
		matPostAct = nan(3,3,intIters); %[actual stim] x [decoded stim] x [time point] x [trials]; value = probability; order= miss / slow / fast
		matPostHet = nan(3,3,intIters); %[actual stim] x [decoded stim] x [time point] x [trials]; value = probability; order= miss / slow / fast
		
		%loop
		for intIter=1:intIters
			%get random trials
			vecMissRand = randperm(intMiss,intNrTrials);
			vecMissLikelihood = vecMiss(vecMissRand(1:(end-1)));
			intDecodeTrialMiss = vecMiss(vecMissRand(end));
			
			vecFastRand = randperm(intFast,intNrTrials);
			vecFastLikelihood = vecFast(vecFastRand(1:(end-1)));
			intDecodeTrialFast = vecFast(vecFastRand(end));
			
			vecSlowRand = randperm(intSlow,intNrTrials);
			vecSlowLikelihood = vecSlow(vecSlowRand(1:(end-1)));
			intDecodeTrialSlow = vecSlow(vecSlowRand(end));
			
			%get likelihoods df/f
			dblLikeActMissMean = mean(mean(matRawDataAct(vecMissLikelihood,vecSelectFrames),2));
			dblLikeActMissSD = std(mean(matRawDataAct(vecMissLikelihood,vecSelectFrames),2));
			
			dblLikeActFastMean = mean(mean(matRawDataAct(vecFastLikelihood,vecSelectFrames),2));
			dblLikeActFastSD = std(mean(matRawDataAct(vecFastLikelihood,vecSelectFrames),2));
			
			dblLikeActSlowMean = mean(mean(matRawDataAct(vecSlowLikelihood,vecSelectFrames),2));
			dblLikeActSlowSD = std(mean(matRawDataAct(vecSlowLikelihood,vecSelectFrames),2));
			
			%get likelihoods heterogeneity
			dblLikeHetMissMean = mean(mean(matRawDataHet(vecMissLikelihood,vecSelectFrames),2));
			dblLikeHetMissSD = std(mean(matRawDataHet(vecMissLikelihood,vecSelectFrames),2));
			
			dblLikeHetFastMean = mean(mean(matRawDataHet(vecFastLikelihood,vecSelectFrames),2));
			dblLikeHetFastSD = std(mean(matRawDataHet(vecFastLikelihood,vecSelectFrames),2));
			
			dblLikeHetSlowMean = mean(mean(matRawDataHet(vecSlowLikelihood,vecSelectFrames),2));
			dblLikeHetSlowSD = std(mean(matRawDataHet(vecSlowLikelihood,vecSelectFrames),2));
			
			
			%get act probabilities for miss trial
			dblActMissTrial = mean(matRawDataAct(intDecodeTrialMiss,vecSelectFrames));
			matPostAct(1,1,intIter) = normpdf(dblActMissTrial,dblLikeActMissMean,dblLikeActMissSD);
			matPostAct(1,2,intIter) = normpdf(dblActMissTrial,dblLikeActSlowMean,dblLikeActSlowSD);
			matPostAct(1,3,intIter) = normpdf(dblActMissTrial,dblLikeActFastMean,dblLikeActFastSD);
			
			%get act probabilities for slow trial
			dblActSlowTrial = mean(matRawDataAct(intDecodeTrialSlow,vecSelectFrames));
			matPostAct(2,1,intIter) = normpdf(dblActSlowTrial,dblLikeActMissMean,dblLikeActMissSD);
			matPostAct(2,2,intIter) = normpdf(dblActSlowTrial,dblLikeActSlowMean,dblLikeActSlowSD);
			matPostAct(2,3,intIter) = normpdf(dblActSlowTrial,dblLikeActFastMean,dblLikeActFastSD);
			
			%get act probabilities for fast trial
			dblActFastTrial = mean(matRawDataAct(intDecodeTrialFast,vecSelectFrames));
			matPostAct(3,1,intIter) = normpdf(dblActFastTrial,dblLikeActMissMean,dblLikeActMissSD);
			matPostAct(3,2,intIter) = normpdf(dblActFastTrial,dblLikeActSlowMean,dblLikeActSlowSD);
			matPostAct(3,3,intIter) = normpdf(dblActFastTrial,dblLikeActFastMean,dblLikeActFastSD);
			
			
			%get het probabilities for miss trial
			vecHetMissTrial = mean(matRawDataHet(intDecodeTrialMiss,vecSelectFrames));
			matPostHet(1,1,intIter) = normpdf(vecHetMissTrial,dblLikeHetMissMean,dblLikeHetMissSD);
			matPostHet(1,2,intIter) = normpdf(vecHetMissTrial,dblLikeHetSlowMean,dblLikeHetSlowSD);
			matPostHet(1,3,intIter) = normpdf(vecHetMissTrial,dblLikeHetFastMean,dblLikeHetFastSD);
			
			%get het probabilities for slow trial
			vecHetSlowTrial = mean(matRawDataHet(intDecodeTrialSlow,vecSelectFrames));
			matPostHet(2,1,intIter) = normpdf(vecHetSlowTrial,dblLikeHetMissMean,dblLikeHetMissSD);
			matPostHet(2,2,intIter) = normpdf(vecHetSlowTrial,dblLikeHetSlowMean,dblLikeHetSlowSD);
			matPostHet(2,3,intIter) = normpdf(vecHetSlowTrial,dblLikeHetFastMean,dblLikeHetFastSD);
			
			%get het probabilities for fast trial
			vecHetFastTrial = mean(matRawDataHet(intDecodeTrialFast,vecSelectFrames));
			matPostHet(3,1,intIter) = normpdf(vecHetFastTrial,dblLikeHetMissMean,dblLikeHetMissSD);
			matPostHet(3,2,intIter) = normpdf(vecHetFastTrial,dblLikeHetSlowMean,dblLikeHetSlowSD);
			matPostHet(3,3,intIter) = normpdf(vecHetFastTrial,dblLikeHetFastMean,dblLikeHetFastSD);
		end
		
		%normalize probabilities
		matPostHet(:,1,:) = matPostHet(:,1,:)/mean(reshape(matPostHet(:,1,:),[1 numel(matPostHet(:,1,:))]));
		matPostHet(:,2,:) = matPostHet(:,2,:)/mean(reshape(matPostHet(:,2,:),[1 numel(matPostHet(:,2,:))]));
		matPostHet(:,3,:) = matPostHet(:,3,:)/mean(reshape(matPostHet(:,3,:),[1 numel(matPostHet(:,3,:))]));
		matPostAct(:,1,:) = matPostAct(:,1,:)/mean(reshape(matPostAct(:,1,:),[1 numel(matPostAct(:,1,:))]));
		matPostAct(:,2,:) = matPostAct(:,2,:)/mean(reshape(matPostAct(:,2,:),[1 numel(matPostAct(:,2,:))]));
		matPostAct(:,3,:) = matPostAct(:,3,:)/mean(reshape(matPostAct(:,3,:),[1 numel(matPostAct(:,3,:))]));
		
		%get fraction correct
		%act
		[dummy,matDecodedActMiss]=max(matPostAct(1,:,:),[],2);
		vecMissActCorrect = sum(squeeze(matDecodedActMiss == 1))/intIters;
		[dummy,matDecodedActSlow]=max(matPostAct(2,:,:),[],2);
		vecSlowActCorrect = sum(squeeze(matDecodedActSlow == 2))/intIters;
		[dummy,matDecodedActFast]=max(matPostAct(3,:,:),[],2);
		vecFastActCorrect = sum(squeeze(matDecodedActFast == 3))/intIters;
		
		%het
		[dummy,matDecodedHetMiss]=max(matPostHet(1,:,:),[],2);
		vecMissHetCorrect = sum(squeeze(matDecodedHetMiss == 1))/intIters;
		[dummy,matDecodedHetSlow]=max(matPostHet(2,:,:),[],2);
		vecSlowHetCorrect = sum(squeeze(matDecodedHetSlow == 2))/intIters;
		[dummy,matDecodedHetFast]=max(matPostHet(3,:,:),[],2);
		vecFastHetCorrect = sum(squeeze(matDecodedHetFast == 3))/intIters;
		
		
		%% plot decoding performance
		%{
	make angle for m/s/f separated by 2/3pi. Calc resultant vector angle of
	triplet decoding, graph aspolar plot
		%}
		cellSaveRespDecode{intPopulation,1} = zeros(2,2,3);%vectors; [act/het] x [theta/rho] x [miss/slow/fast]
		cellSaveRespDecode{intPopulation,2} = zeros(2,3);%distance; [act/het] x [miss/slow/fast]
		for intRespType=1:3
			if intRespType == 1
				strTitle = 'Miss';
				matDecodedAct = matDecodedActMiss;
				matDecodedHet = matDecodedHetMiss;
				thetaCorr = (4/3)*pi;
				cellColor = {'r','k','k'};
			elseif intRespType == 2
				strTitle = 'Slow';
				matDecodedAct = matDecodedActSlow;
				matDecodedHet = matDecodedHetSlow;
				thetaCorr = (0/3)*pi;
				cellColor = {'k','y','k'};
			elseif intRespType == 3
				strTitle = 'Fast';
				matDecodedAct = matDecodedActFast;
				matDecodedHet = matDecodedHetFast;
				thetaCorr = (2/3)*pi;
				cellColor = {'k','k','g'};
			end
			subplot(2,2,intRespType)
			%calc corr vector
			[xCorr,yCorr] = pol2cart(thetaCorr,1);
			
			%calc resultant vector act
			[x1,y1] = pol2cart((4/3)*pi,sum(squeeze(matDecodedAct == 1))/intIters);%miss
			[x2,y2] = pol2cart((0/3)*pi,sum(squeeze(matDecodedAct == 2))/intIters);%slow
			[x3,y3] = pol2cart((2/3)*pi,sum(squeeze(matDecodedAct == 3))/intIters);%fast
			x=x1+x2+x3;
			y=y1+y2+y3;
			[thetaAct,rhoAct]=cart2pol(x,y);
			[thetaDiffAct,rhoDiffAct] = cart2pol(xCorr-x,yCorr-y);
			
			%resultant vector het
			[x1,y1] = pol2cart((4/3)*pi,sum(squeeze(matDecodedHet == 1))/intIters);%miss
			[x2,y2] = pol2cart((0/3)*pi,sum(squeeze(matDecodedHet == 2))/intIters);%slow
			[x3,y3] = pol2cart((2/3)*pi,sum(squeeze(matDecodedHet == 3))/intIters);%fast
			x=x1+x2+x3;
			y=y1+y2+y3;
			[thetaHet,rhoHet]=cart2pol(x,y);
			[thetaDiffHet,rhoDiffHet] = cart2pol(xCorr-x,yCorr-y);
			
			%save data
			cellSaveRespDecode{intPopulation,1}(1,1,intRespType) = thetaAct;%vectors;[act/het] x [theta/rho] x [miss/slow/fast]
			cellSaveRespDecode{intPopulation,1}(1,2,intRespType) = rhoAct;
			cellSaveRespDecode{intPopulation,1}(2,1,intRespType) = thetaHet;
			cellSaveRespDecode{intPopulation,1}(2,2,intRespType) = rhoHet;
			cellSaveRespDecode{intPopulation,2}(1,intRespType) = rhoDiffAct;%distance; [act/het] x [miss/slow/fast]
			cellSaveRespDecode{intPopulation,2}(2,intRespType) = rhoDiffHet;%distance; [act/het] x [miss/slow/fast]
			
			%make polar background
			[x,y] = pol2cart((4/3)*pi,0.9);%miss
			text(x,y,'Miss','Color',cellColor{1});
			hold on
			[x,y] = pol2cart((0/3)*pi,0.7);%slow
			text(x,y,'Slow','Color',cellColor{2});
			[x,y] = pol2cart((2/3)*pi,0.9);%fast
			text(x,y,'Fast','Color',cellColor{3});
			[x,y] = pol2cart((1/3)*pi,1);%slow/fast
			plot([0 x],[0 y],'k');
			[x,y] = pol2cart((3/3)*pi,1);%miss/fast
			plot([0 x],[0 y],'k');
			[x,y] = pol2cart((5/3)*pi,1);%slow/miss
			plot([0 x],[0 y],'k');
			ang=0:0.001:2*pi;
			x=cos(ang);
			y=sin(ang);
			plot(x,y, 'k');
			
			%plot data
			[x,y]=pol2cart(thetaAct,rhoAct);
			scatter(x,y,'bx');
			[x,y]=pol2cart(thetaHet,rhoHet);
			scatter(x,y,'kx');
			hold off;
			xlim([-1 1]);
			ylim([-1 1]);
			title(sprintf('%s trials',strTitle))
		end
		
		%plot summary
		subplot(2,2,4);
		matPerformance = imnorm(1-abs(cellSaveRespDecode{intPopulation,2}));
		vecActPerfomance = matPerformance(1,:);
		vecHetPerfomance = matPerformance(2,:);
		dblOffset=0.1;
		scatter((1-dblOffset)*ones(size(vecActPerfomance)),vecActPerfomance,'bx')
		hold on
		dblMeanAct=mean(vecActPerfomance);
		dblSDAct = std(vecActPerfomance);
		errorbar((1+dblOffset),dblMeanAct,dblSDAct/sqrt(length(vecActPerfomance)),'xb')
		dblMeanHet=mean(vecHetPerfomance);
		dblSDHet = std(vecHetPerfomance);
		scatter((2-dblOffset)*ones(size(vecHetPerfomance)),vecHetPerfomance,'kx')
		errorbar((2+dblOffset),dblMeanHet,dblSDHet/sqrt(length(vecHetPerfomance)),'xk')
		hold off
		title('Normalized resp type decoding')
		ylabel('Norm decod perf')
		set(gca,'XTick',[1 2],'XTickLabel',{'dF/F0', 'Heterogen'})
		
		%save
		if sParams.boolSavePlots
			drawnow;
			strFig = sprintf('%s_resptypedecoding_over_time_pop%d_raw',strSes,intPopulation);
			export_fig([strFig '.tif']);
			export_fig([strFig '.pdf']);
		end
	end
	
	%% save data structures
	if sParams.boolSaveData
		vecClock = fix(clock);
		strDate = num2str(vecClock(1:3));
		strRecs = num2str(vecRecordings);
		strFile = ['data_aggregateAD' strSes '_' strrep(strDate,'    ','_') '_' strRecs(~isspace(strRecs))];
		%save(['D:\Data\Results\stimdetection\' strFile],'-v7.3');
		%if exist(['D:\Data\Results\stimdetection\' strFile '.mat'],'file')
		%	load(['D:\Data\Results\stimdetection\' strFile]);
		%end
		save(['D:\Data\Results\stimdetection\' strFile],'cellSave*','indSelectHitCorrelatedNeurons','indSelectHitUncorrelatedNeurons','-v7.3');
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