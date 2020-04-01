%{
figures to save and take mean over animals
- cellSavePreStim: dF/F detect/no-detect preceding stimulus frames [5 (c) x 2 (d/nd) cell array] x n (animals) with every cell = vector of values
- cellSaveBehavDetect: behavioral detection [5 (c) x n (animals)] cell array with every cell = vector of detect/no detect values (1/0)
- cellSaveBehavRT: behavioral RT [5 (c) x n (animals)] cell array with every cell = vector of RTs
- cellSaveBaseSubtrStimAct: dF/F detect/no-detect baseline-subtracted activity [5 (c) x 2 (d/nd) cell array] x n (animals) with every cell = vector of values
- cellSaveDuringStim: dF/F detect/no-detect during stimulus frames [5 (c) x 2 (d/nd) cell array] x n (animals) with every cell = vector of values
- cellSaveTraceAct: dF/F trace detect/no-detect pref/non-pref [2 (pref/no-pref) x 2 (detect/no-detect) x 6 (c)] x n (animals) with every cell = [Nreps x Nframes] matrix of trace dF/F values
-  AD: cellSaveMatrices: blob matrices in [1 x 4 (detect/no-detect/diff/significance)] cell array with [5 (c) x 31 (bins)] matrix
-  AD: cellSaveHCARNeuronConsistency: [1 x i (neuronal population)] cell array with vector of shuffled correlation values, where first value is actual data
-  AD: cellSaveSignalCorrs: signal correlations [3 (within-HCAR, between, within-non-HCAR) x n (animals/blocks)] with every cell = vector of pairwise correlation values
-  AD: cellSaveNoiseCorrs: noise correlations [3 (within-HCAR, between, within-non-HCAR) x n (animals/blocks)] with every cell = vector of pairwise correlation values
-  AD: cellSaveHCAR: normalized increase in dF/F for detection trials [5 (c) x n (animals/blocks)] cell array with each cell = vector of neurons with normalized dF/F increase for detection trials
-  AD: cellSaveSignalCorrsBiDir: signal correlations [3 (within-HCAR, between, within-non-HCAR) x n (animals/blocks)] with every cell = vector of pairwise correlation values
-  AD: cellSaveNoiseCorrsBiDir: noise correlations [3 (within-HCAR, between, within-non-HCAR) x n (animals/blocks)] with every cell = vector of pairwise correlation values
-  AD: cellSaveNormActDissim: within-group z-scored activation dissimilarity [2 (high/low HCAR) x 2 (detect/no-detect)] x n (animals/blocks) with every cell = vector of dissimilarity values
-  AD: cellSaveDissimCorrITD: inter-trial-distance dependence of assembly consistency [2 (high/low HCAR) x 2 (hits/misses) x n (animals/blocks)] with every cell = [p x 2] matrix with intertrial distance and correlation values for each trial pair
	- extra graph of above with means comparison of hits/misses high/low HCAR

%}
close all
clear all
%% parameters
strFigDir = 'D:\Data\ResultsAstroAnalysis\_meta';
strDataDir = 'D:\Data\ResultsAstroAnalysis';
cellInclude = {'20140207','20140314','20140425','20140507','20140530','20140604','20140711','20140715'};%exclude:20140430
%cellInclude = {'20140207','20140507','20140711','20140715'};%exclude:20140430
boolSavePlots = true;
strAnalyzeType = 'astrocyte';
boolNPS = 0; 
boolEvents = false; 
intFigCounter = 0;
strTag = '';
if boolNPS
	cellInclude = cellfun(@strcat,cellInclude,cellfill('NPS',size(cellInclude)),'UniformOutput',false);
	strFigDir = [strFigDir 'NPS'];
	strTag = ['NPS' strTag];
end
if boolEvents
	cellInclude = cellfun(@strcat,cellfill('Ev',size(cellInclude)),cellInclude,'UniformOutput',false);
	strTag = ['Ev' strTag];
end
cellInclude = cellfun(@strcat,cellfill(strAnalyzeType,size(cellInclude)),cellInclude,'UniformOutput',false);
strTag = [strAnalyzeType strTag];

%% pre-allocate
intAnimal = 0;
vecAnimalSource = [];
vecSubPopulationSource = [];

%% create aggregate data structures
sDir = dir(strDataDir);
strOldDir = cd(strDataDir);
intCounter = 0;
cellMeanOTCR = [];
matMeanOTCR = [];
for intFile=1:numel(sDir)
	%check if data part 2
	strFile = sDir(intFile).name;
	strRec = getFlankedBy(strFile,['data_aggregate_'],'_');
	intRec = find(strcmp(strRec,cellInclude),1);
	if ~isempty(intRec)
		%load file
		load(strFile);
		fprintf('Loaded %s [%s]\n',strFile,getTime);
		
		%update population source
		intNrPops = length(cellSaveMeanHit);
		intAnimal = intAnimal + 1;
		vecAnimalSource(end+1) = intAnimal;
		vecSubPopulationSource(end+1) = intAnimal;
		if intNrPops > 1 %group populations
			
			%concatenate cellSaveMean
			cellSaveMeanHit = {cat(1,cellSaveMeanHit{1},cellSaveMeanHit{2})};
			cellSaveMeanMiss = {cat(1,cellSaveMeanMiss{1},cellSaveMeanMiss{2})};
		end
			
		%assign data
		matOTCR = cat(4,cellSaveMeanHit{1},cellSaveMeanMiss{1});%[object x time x contrast x response]
		cellMeanOTCR{intAnimal} = matOTCR; 
		matMeanOTCR = cat(1,matMeanOTCR,matOTCR); %[object x time x contrast x response]
		
		
		intCounter = intCounter + 1;
	end
end

%% meta analyses part 1
cd(strFigDir);
%close all;
intAnimals = intCounter;
vecContrasts = [0 0.5 2 8 32 100];
vecLimY = [-0.01 0.04];
vecDiffLimY = [-0.02 0.02];
cellLegend(1) = {'r=FA,b=CR'};
cellLegend(2:6) = {'r=hit,b=miss'};
cellLegendDiff(1) = {'diff FA/CR'};
cellLegendDiff(2:6) = {'diff hit/miss'};

%% plot over cells
matMeanHitOTC = matMeanOTCR(:,:,:,1);
matMeanMissOTC = matMeanOTCR(:,:,:,2);	
intT = size(matMeanHitOTC,2);
vecX = ((1:intT)-round(intT/2))/25.4;

figure
for intContrast=1:size(matMeanHitOTC,3)
	dblC = vecContrasts(intContrast);
	
	subplot(2,3,intContrast);
	intN = size(matMeanHitOTC,1);
	
	vecMeanHit = nanmean(matMeanHitOTC(:,:,intContrast),1);
	vecSEMHit = nanstd(matMeanHitOTC(:,:,intContrast),[],1)/sqrt(intN);
	
	vecMeanMiss = nanmean(matMeanMissOTC(:,:,intContrast),1);
	vecSEMMiss = nanstd(matMeanMissOTC(:,:,intContrast),[],1)/sqrt(intN);
	
	errorfill(vecX,vecMeanHit,vecSEMHit,[1 0 0],[1 0.7 0.7]);
	hold on
	errorfill(vecX,vecMeanMiss,vecSEMMiss,[0 0 1],[0.7 0.7 1]);
	hold off
	ylim(vecLimY);
	xlabel('Time after stimulus onset (s)');
	ylabel(sprintf('Mean dF/F0 over %ss',strAnalyzeType));
	title(sprintf('Stim contrast %.1f%%; %s',vecContrasts(intContrast),cellLegend{intContrast}))
end
%save fig
strFigTitle = ['ActAcrossCells' strTag];
if boolSavePlots
	intFigCounter = intFigCounter + 1;
	drawnow;
	jFig = get(handle(gcf), 'JavaFrame');
	jFig.setMaximized(true);
	figure(gcf);
	drawnow;
	strFig = sprintf('Meta%d_%s_raw',intFigCounter,strFigTitle);
	export_fig([strFig '.tif']);
	export_fig([strFig '.pdf']);
end

%% plot difference
matMeanHitOTC = matMeanOTCR(:,:,:,1);
matMeanMissOTC = matMeanOTCR(:,:,:,2);	
matMeaDiffOTC = matMeanHitOTC - matMeanMissOTC;
vecX = ((1:intT)-round(intT/2))/25.4;

figure
for intContrast=1:size(matMeanHitOTC,3)
	dblC = vecContrasts(intContrast);
	
	subplot(2,3,intContrast);
	intN = size(matMeanHitOTC,1);
	
	vecMeanDiff = nanmean(matMeaDiffOTC(:,:,intContrast),1);
	vecSEMDiff = nanstd(matMeaDiffOTC(:,:,intContrast),[],1)/sqrt(intN);
	
	errorfill(vecX,vecMeanDiff,vecSEMDiff,[0.5 0.5 0.5],[0.7 0.7 0.7]);
	ylim(vecDiffLimY);
	xlabel('Time after stimulus onset (s)');
	ylabel(sprintf('Mean %s(dF/F0) over %ss',getGreek(4),strAnalyzeType));
	title(sprintf('Stim contrast %.1f%%; %s',vecContrasts(intContrast),cellLegendDiff{intContrast}))
end
%save fig
strFigTitle = ['DiffAcrossCells' strTag];
if boolSavePlots
	intFigCounter = intFigCounter + 1;
	drawnow;
	jFig = get(handle(gcf), 'JavaFrame');
	jFig.setMaximized(true);
	figure(gcf);
	drawnow;
	strFig = sprintf('Meta%d_%s_raw',intFigCounter,strFigTitle);
	export_fig([strFig '.tif']);
	export_fig([strFig '.pdf']);
end

%% plot over animals
matMeanATCR = cell2mat(cellfun(@mean,cellMeanOTCR,cellfill(1,size(cellMeanOTCR)),'UniformOutput',false)');
matMeanHitATC = matMeanATCR(:,:,:,1);
matMeanMissATC = matMeanATCR(:,:,:,2);	
intT = size(matMeanHitATC,2);
vecX = ((1:intT)-round(intT/2))/25.4;

figure
for intContrast=1:size(matMeanHitATC,3)
	dblC = vecContrasts(intContrast);
	
	subplot(2,3,intContrast);
	intN = size(matMeanHitATC,1);
	
	vecMeanHit = nanmean(matMeanHitATC(:,:,intContrast),1);
	vecSEMHit = nanstd(matMeanHitATC(:,:,intContrast),[],1)/sqrt(intN);
	
	vecMeanMiss = nanmean(matMeanMissATC(:,:,intContrast),1);
	vecSEMMiss = nanstd(matMeanMissATC(:,:,intContrast),[],1)/sqrt(intN);
	
	errorfill(vecX,vecMeanHit,vecSEMHit,[1 0 0],[1 0.7 0.7]);
	hold on
	errorfill(vecX,vecMeanMiss,vecSEMMiss,[0 0 1],[0.7 0.7 1]);
	hold off
	ylim(vecLimY);
	xlabel('Time after stimulus onset (s)');
	ylabel(sprintf('Mean %s dF/F0 over animals ',strAnalyzeType));
	title(sprintf('Stim contrast %.1f%%; %s',vecContrasts(intContrast),cellLegend{intContrast}))
end
%save fig
strFigTitle = ['ActAcrossAnim' strTag];
if boolSavePlots
	intFigCounter = intFigCounter + 1;
	drawnow;
	jFig = get(handle(gcf), 'JavaFrame');
	jFig.setMaximized(true);
	figure(gcf);
	drawnow;
	strFig = sprintf('Meta%d_%s_raw',intFigCounter,strFigTitle);
	export_fig([strFig '.tif']);
	export_fig([strFig '.pdf']);
end

%% plot difference over animals
matMeanDiffATC = matMeanHitATC - matMeanMissATC;
vecX = ((1:intT)-round(intT/2))/25.4;

figure
for intContrast=1:size(matMeanDiffATC,3)
	dblC = vecContrasts(intContrast);
	
	subplot(2,3,intContrast);
	intN = size(matMeanDiffATC,1);
	
	vecMeanDiff = nanmean(matMeanDiffATC(:,:,intContrast),1);
	vecSEMDiff = nanstd(matMeanDiffATC(:,:,intContrast),[],1)/sqrt(intN);
	
	errorfill(vecX,vecMeanDiff,vecSEMDiff,[0.5 0.5 0.5],[0.7 0.7 0.7]);
	ylim(vecDiffLimY);
	xlabel('Time after stimulus onset (s)');
	ylabel(sprintf('Mean %s %s(dF/F0) over animals',strAnalyzeType,getGreek(4)));
	title(sprintf('Stim contrast %.1f%%; %s',vecContrasts(intContrast),cellLegendDiff{intContrast}))
end

%save fig
strFigTitle = ['DiffAcrossAnim' strTag];
if boolSavePlots
	intFigCounter = intFigCounter + 1;
	drawnow;
	jFig = get(handle(gcf), 'JavaFrame');
	jFig.setMaximized(true);
	figure(gcf);
	drawnow;
	strFig = sprintf('Meta%d_%s_raw',intFigCounter,strFigTitle);
	export_fig([strFig '.tif']);
	export_fig([strFig '.pdf']);
end

%% clean up
cd(strOldDir);