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
strFigDir = 'D:\Data\Results\stimdetection\metaNS_Supp';
strDataDir = 'D:\Data\Results\stimdetection';
cellInclude = {'20140207','20140314','20140425','20140507','20140530','20140604','20140711','20140715'};%exclude:20140430
boolSavePlots = true;
intFigCounter = 0;
boolUseNeuropilSubtraction = true;
if boolUseNeuropilSubtraction
	strFigDir = [strFigDir 'NPS'];
end

%% pre-allocate
intAnimal = 0;
vecAnimalSource = [];
vecSubPopulationSource = [];

cellHitOverTime = {};
matContAct = [];
matDetect = [];
matBehavDetect = [];
cellHetTime = {};
cellActTime = {};
matQuintileRemoval = []; %contrast x quintile x hit/miss x animal
%matNeuronConsistency = [];
%vecRecordingStability = [];
cellRTDependency = {};
matDetectLowHighActHetRNR = [];
matHetChange = []; %[contrast x behav x quintile x pre/post]
		
%% create aggregate data structures
sDir = dir(strDataDir);
strOldDir = cd(strDataDir);
intCounterF1 = 0;
intCounterF2 = 0;

for intFile=1:numel(sDir)
	%check if data part 1
	strFile = sDir(intFile).name;
	strRec = getFlankedBy(strFile,'data_aggregateNS_Supp','_');
	if length(strRec) > 3 && boolUseNeuropilSubtraction && strcmp(strRec(1:3),'NPS')
		intRec = find(strcmp(strRec(4:end),cellInclude),1);
	elseif ~boolUseNeuropilSubtraction
		intRec = find(strcmp(strRec,cellInclude),1);
	else
		intRec = [];
	end
	
	if ~isempty(intRec)
		strRec
		%load file
		load(strFile);
		
		%update population source
		intNrPops = size(cellSaveMatContAct,1);
		intAnimal = intAnimal + 1;
		vecAnimalSource(end+1) = intAnimal;
		vecSubPopulationSource(end+1) = intAnimal;
		if intNrPops > 1 %group populations
			%take mean over matrices
			cellSaveMatContAct = {nanmean(cat(3,cellSaveMatContAct{1},cellSaveMatContAct{2}),3)};
			cellSaveMatDetect = {nanmean(cat(3,cellSaveMatDetect{1},cellSaveMatDetect{2}),3)};
			cellSaveMatDetectLowHighActHet = {mean(cat(5,cellSaveMatDetectLowHighActHet{1},cellSaveMatDetectLowHighActHet{2}),5)};
			
			cellSaveHetQs = {cat(2,cellSaveHetQs{1},cellSaveHetQs{2})};
			cellSaveHetAll = {cat(2,cellSaveHetAll{1},cellSaveHetAll{2})};
			cellSaveRespTimeTrials = {cat(2,cellSaveRespTimeTrials{1},cellSaveRespTimeTrials{2})};
			cellSaveContrastTrials = {cat(2,cellSaveContrastTrials{1},cellSaveContrastTrials{2})};
			%cellSaveNeuronConsistency = {cat(2,cellSaveNeuronConsistency{1},cellSaveNeuronConsistency{2})};
			%cellSaveMaxEpochs = {cellSaveMaxEpochs{1} cellSaveMaxEpochs{2}};
			
			for intEl=1:size(cellSaveRTDependency,1)
				cellSaveRTDependency{intEl,1} = [cellSaveRTDependency{intEl,1}; cellSaveRTDependency{intEl,2}];
			end
			cellSaveRTDependency(:,2) = [];
			
			cellSaveChange = {nanmean(cat(5,cellSaveChange{1},cellSaveChange{1}),5)}; %[contrast x behav x quintile x pre/post]
		
		end
		
		%quintile pre-pro
		matHet = [cellSaveHetQs{1};cellSaveHetAll{1}];
		matQ = nan(6,6,2);%contrast x quintile x hit/miss
		vecC = unique(cellSaveContrastTrials{1});
		for intContrast=1:6
			indC = cellSaveContrastTrials{1} == vecC(intContrast);
			for intQuintile=1:6
				for intHitMiss=1:2
					if intHitMiss == 1,indR = ~isnan(cellSaveRespTimeTrials{1});
					else indR = isnan(cellSaveRespTimeTrials{1});end
					if sum(indC&indR) == 0
						vecR=find(indC);
						[dummy,vecI]=findmin(cellSaveRespTimeTrials{1}(vecR),2);
						cellSaveRespTimeTrials{1}(vecR(vecI)) = nan;
						if intHitMiss == 1,indR = ~isnan(cellSaveRespTimeTrials{1});
						else indR = isnan(cellSaveRespTimeTrials{1});end
					end
					matQ(intContrast,intQuintile,intHitMiss) = mean(matHet(intQuintile,indC&indR));
				end
			end
		end
		
		%assign data
		cellHitOverTime = cat(1,cellHitOverTime,cellSaveHitOverTime);
		matDetect = cat(3,matDetect,cellSaveMatDetect{1});
		matDetectLowHighActHetRNR = cat(5,matDetectLowHighActHetRNR,cellSaveMatDetectLowHighActHet{1});
		matContAct = cat(3,matContAct,cellSaveMatContAct{1});
		matBehavDetect = cat(2,matBehavDetect,cellfun(@mean,cellSaveBehavDetect));
		matQuintileRemoval = cat(4,matQuintileRemoval,matQ);
		%matNeuronConsistency = cat(3,matNeuronConsistency,cellSaveNeuronConsistency{1});
		%vecRecordingStability = cat(2,vecRecordingStability,cellSaveMaxEpochs{1});
		cellRTDependency = cat(2,cellRTDependency,cellSaveRTDependency);
		matHetChange = cat(5,matHetChange,cellSaveChange{1}); %[contrast x behav x quintile x pre/post]
		intCounterF1 = intCounterF1 + 1;
	end
	
	%check if data part 2
	strFile = sDir(intFile).name;
	strRec = getFlankedBy(strFile,'data_aggregateNS','_');
	
	if length(strRec) > 3 && boolUseNeuropilSubtraction && strcmp(strRec(1:3),'NPS')
		intRec = find(strcmp(strRec(4:end),cellInclude),1);
	else
		intRec = find(strcmp(strRec,cellInclude),1);
	end
	if ~isempty(intRec)
		%load file
		load(strFile);
		
		%update population source
		intNrPops = size(cellSaveHetTime,2);
		if intNrPops > 1 %group populations
			%concatenate heterogen-time
			for intC=1:size(cellSaveHetTime,1)
				cellSaveHetTime{intC,1} = cat(1,cellSaveHetTime{intC,1},cellSaveHetTime{intC,2});
				cellSaveActTime{intC,1} = cat(1,cellSaveActTime{intC,1},cellSaveActTime{intC,2});
			end
			cellSaveHetTime = cellSaveHetTime(:,1);
			cellSaveActTime = cellSaveActTime(:,1);
		end
		
		%assign data
		cellHetTime = cat(2,cellHetTime,cellSaveHetTime);
		cellActTime = cat(2,cellActTime,cellSaveActTime);
		intCounterF2 = intCounterF2 + 1;
	end
end

%check for inconsistencies
%if ~((intCounterF1 == intCounterF2) && (intCounterF2 == intCounterF3) && (intCounterF3 == intCounterF4) && (intCounterF4 == intCounterF5))
%	warning([mfilename ':InconsistentFileNumbers'],'File counters inconsistent; F1=%d; F2=%d; F3=%d; F4=%d; F5=%d',intCounterF1,intCounterF2,intCounterF3,intCounterF4,intCounterF5)
%end

%% meta analyses part 1
cd(strFigDir);
%close all;
intAnimals = intCounterF1;
vecContrasts = [0 0.5 2 8 32 100];

%% heterogeneity change
matHetChangeMean = nanmean(matHetChange,5);%[contrast x behav x quintile x pre/post]
matHetChangeSD = nanstd(matHetChange,[],5);
matNrAnimals = sum(~isnan(matHetChange),5);

intSegments = size(matHetChangeMean,3);
hHetChange = figure;
for intC=1:size(matHetChangeMean,1)
	%plot increase
	subplot(2,3,intC)
	hold on;
	
	%miss
	meanVecPre = squeeze(matHetChangeMean(intC,1,:,1));
	meanVecStim = squeeze(matHetChangeMean(intC,1,:,2));
	plot(repmat([1 2],[intSegments 1])',[meanVecPre meanVecStim]','rx-')
	
	%slow
	meanVecPre = squeeze(matHetChangeMean(intC,2,:,1));
	meanVecStim = squeeze(matHetChangeMean(intC,2,:,2));
	plot(repmat([3 4],[intSegments 1])',[meanVecPre meanVecStim]','bx-')
	
	%fast
	meanVecPre = squeeze(matHetChangeMean(intC,3,:,1));
	meanVecStim = squeeze(matHetChangeMean(intC,3,:,2));
	plot(repmat([5 6],[intSegments 1])',[meanVecPre meanVecStim]','gx-')
	
	
	
	hold off;
	ylabel('Heterogeneity');
	ylim([0.5 1.5]);
	title(sprintf('Contrast %.1f%%; n-trials, M,T1=%d,T2=%d,T3=%d; \nS,T1=%d,T2=%d,T3=%d; F,T1=%d,T2=%d,T3=%d',vecContrasts(intC),matNrAnimals(intC,1,:,1),matNrAnimals(intC,2,:,1),matNrAnimals(intC,3,:,1)));
	
end

if boolSavePlots
	figure(hHetChange);
	intFigCounter = intFigCounter + 1;
	drawnow;
	jFig = get(handle(gcf), 'JavaFrame');
	jFig.setMaximized(true);
	figure(gcf);
	drawnow;
	strFig = sprintf('Meta%d_heterogeneity_change_raw',intFigCounter);
	export_fig([strFig '.tif']);
	export_fig([strFig '.pdf']);
end

%% reaction-time dependency
vecX = [0 3];
hRTDependenceFig = figure;
hMeasuresR2 = figure;
intActTypes = size(cellRTDependency,1)-1;
matR2 = nan(intActTypes,intAnimals);
vecPs = nan(1,size(cellRTDependency,1)-1);
for intActType=1:intActTypes
	figure(hRTDependenceFig);
	if intActType == 1
		%heterogeneity
		intHCAR = 2;
		intNonHCAR = 3;
		strLabelY = 'Heterogeneity (mean sigma)';
		strC = 'k';
	elseif intActType == 2
		%z-scored act
		strLabelY = 'Z-scored dF/F0 (sigma)';
		strC = 'r';
	elseif intActType == 3
		%dF/F act
		strLabelY = 'dF/F0';
		strC = 'b';
	elseif intActType == 4
		%variance
		strLabelY = 'Variance (sigma^2)';
		strC = 'g';
	elseif intActType == 5
		%sparseness
		strLabelY = 'Sparseness (kurtosis)';
		strC = 'y';
	elseif intActType == 6
		%sparseness
		strLabelY = 'Pearson-like (r)';
		strC = 'c';
	elseif intActType == 7
		%sparseness
		strLabelY = 'SD Pearson-like (sigma r)';
		strC = 'm';
	elseif intActType == 8
		%sparseness
		strLabelY = 'Pref-pop dF/F0';
		strC = 'b';
	elseif intActType == 9
		%sparseness
		strLabelY = 'Pref-pop Z-scored dF/F0 (sigma)';
		strC = 'r';
	elseif intActType == 10
		%sparseness
		strLabelY = 'Mean sliding correlation (semblance)';
		strC = 'c';
	elseif intActType == 11
		%sparseness
		strLabelY = 'SD sliding correlation (semblance)';
		strC = 'm';
	end

	matReg = zeros(intAnimals,length(vecX));
	vecSlopes = zeros(intAnimals,1);
	vecR2 = zeros(intAnimals,1);
	subplot(2,intActTypes,intActType);
	hold on;
	vecAggRTs = [];
	vecAggAct = [];
	for intPopulation=1:intAnimals
		vecAggRTs = [vecAggRTs cellRTDependency{1,intPopulation}'];
		vecAggAct = [vecAggAct cellRTDependency{intActType+1,intPopulation}'];
		
		%do regressions per population
		sStatsC=regstats(cellRTDependency{intActType+1,intPopulation}',cellRTDependency{1,intPopulation},'linear');
		vecY = polyval(sStatsC.beta([2 1]),vecX);
		matReg(intPopulation,:) = vecY;
		vecSlopes(intPopulation) = sStatsC.beta(2);
		vecR2(intPopulation) = sStatsC.rsquare;
		plot(vecX,vecY,'Color',[0.5 0.5 1.0]);
	end
	matR2(intActType,:) = vecR2;
	
	%plot means
	vecMeanY = mean(matReg,1);
	plot(vecX,vecMeanY,'Color',strC,'LineWidth',2)
	
	%ttest
	[h,dblPSlope] = ttest(vecSlopes);
	[corrected_p, h]=bonf_holm([ones([1 intActTypes*2-1]) dblPSlope],0.05);
	dblPSlopeCorr = corrected_p(end);
	vecPs(intActType) = dblPSlope;
	title(sprintf('Mean of linear regressions over animals; slope-p=%.3f',dblPSlopeCorr))
	xlabel('Reaction Time (s)')
	ylabel(strLabelY)
	xlim([0 3]);
	
	%plot aggregate
	subplot(2,intActTypes,intActType+intActTypes);
	scatter(vecAggRTs,vecAggAct,'bx')
	hold on
	drawnow
	
	%perform regressions
	sStatsC=regstats(vecAggAct,vecAggRTs,'linear');
	vecX = [0 3];
	vecY = polyval(sStatsC.beta([2 1]),vecX);
	plot(vecX,vecY,'r')
	hold off
	xlim([0 3]);
	
	title(sprintf('Aggregate data set; slope=%.3f; p=%.3f; R^2=%.3f',...
		sStatsC.beta(2),sStatsC.tstat.pval(2),sStatsC.rsquare))
	%vecPs((intActType-1)*2+2) = sStatsC.tstat.pval(2);
	xlabel('Reaction Time (s)')
	ylabel(strLabelY)
	
	%plot R2
	figure(hMeasuresR2);
	hold on;
	errorbar(intActType,mean(vecR2),std(vecR2)/sqrt(intAnimals),'x')
	cellLabels{intActType} = strLabelY;
	hold off;
end
set(gca,'XTick',1:intActTypes,'XTickLabel',cellLabels);
ylabel('Explained variance (R^2)');
%ylim([0 0.2])

[corrected_p, h]=bonf_holm(vecPs,0.05);
fprintf('\nCorrected p-values RT dependence: %.3f\n',corrected_p)

matP = nan(intActTypes,intActTypes);
for intActType1=1:(intActTypes-1)
	for intActType2=(intActType1+1):intActTypes
		[h,matP(intActType2,intActType1)] = ttest(matR2(intActType1,:),matR2(intActType2,:));
	end
end
vecP = matP(~isnan(matP));
[h,crit,vecFDR]=fdr_bh(matP(~isnan(matP)),0.05);

matV = ~isnan(matP);
matP_FDR = nan(size(matP));
intCounter=0;
for intActType1=1:(intActTypes-1)
	for intActType2=(intActType1+1):intActTypes
		if matV(intActType2,intActType1) == 1
			intCounter = intCounter + 1;
			matP_FDR(intActType2,intActType1) = vecFDR(intCounter);
		end
	end
end
matP_FDR=matP_FDR';
matP = matP';
%return

if boolSavePlots
	figure(hRTDependenceFig);
	intFigCounter = intFigCounter + 1;
	drawnow;
	jFig = get(handle(gcf), 'JavaFrame');
	jFig.setMaximized(true);
	figure(gcf);
	drawnow;
	strFig = sprintf('Meta%d_activationdissimilarity_RTDependency_raw',intFigCounter);
	export_fig([strFig '.tif']);
	export_fig([strFig '.pdf']);
	
	figure(hMeasuresR2);
	intFigCounter = intFigCounter + 1;
	drawnow;
	strFig = sprintf('Meta%d_RTDependency_R2_raw',intFigCounter);
	export_fig([strFig '.tif']);
	export_fig([strFig '.pdf']);
end

%% recording stability
%{
h=figure;
vecBins = 0:0.04:1.1;
vecPlotBins = 0.02:0.04:1.12;
vecV = hist(vecRecordingStability,vecBins);
bar(vecPlotBins,vecV,1)
xlabel('Maximum below threshold epoch duration (s)');
ylabel('Number of neurons');
xlim([0 1]);
if boolSavePlots
	intFigCounter = intFigCounter + 1;
	drawnow;
	strFig = sprintf('Meta%d_recording_stability_raw',intFigCounter);
	export_fig([strFig '.tif']);
	export_fig([strFig '.pdf']);
end
%}
		
%% HCN consistency
%{
hConsistencyFigure = figure;
vecConsistencyC = [0.5 2 8 32 100];
vecConsistencySelect = 2:6;
matNeuronConsistency = matNeuronConsistency(vecConsistencySelect,vecConsistencySelect,:);
matMeanNeuronConsistency = mean(matNeuronConsistency,3);
imagesc(matMeanNeuronConsistency,[-0.15 0.15]);
colormap(redblue);
colorbar;
set(gca,'XTick',1:length(vecConsistencyC),'XTickLabel',vecConsistencyC)
set(gca,'YTick',1:length(vecConsistencyC),'YTickLabel',vecConsistencyC)
ylabel('Contrast (%)')
xlabel('Contrast (%)')
matMeanNeuronConsistency(diag(diag(true(size(matMeanNeuronConsistency))))) = nan;

%statistics
matSelect = tril(true(size(matMeanNeuronConsistency)),-1);
matSelectAll = repmat(matSelect,[1 1 size(matNeuronConsistency,3)]);
vecVals = matMeanNeuronConsistency(matSelect);
vecValsAll = matNeuronConsistency(matSelectAll);
[boolAllH,dblAllP] = ttest(vecValsAll);
[boolH,dblP] = ttest(vecVals);
title(sprintf('Mean corr of highest 33%% of hit-corr neurons; ttest rho~=0, p=%.3f (all); p=%.3f (means)',dblAllP,dblP))

if boolSavePlots
	intFigCounter = intFigCounter + 1;
	drawnow;
	strFig = sprintf('Meta%d_HCN_consistency_raw',intFigCounter);
	export_fig([strFig '.tif']);
	export_fig([strFig '.pdf']);
end
%}

%% hit over time
%regression per animal
vecSlopes = nan(1,intAnimals);
vecIntercepts = nan(1,intAnimals);
h=figure;

%resp over rep (C)
subplot(2,2,1)
hold on
for intAnimal=1:intAnimals
	%resp C
	matBehav = cellHitOverTime{intAnimal,1};
	matRep = meshgrid(1:size(cellHitOverTime{intAnimal,1},2),1:size(cellHitOverTime{intAnimal,1},1));
	
	sStats=regstats(matBehav(:),matRep(:),'linear');
	vecX = [1 max(matRep(:))];
	vecY = polyval(sStats.beta([2 1]),vecX);
	plot(vecX,vecY,'Color',[0.5 0.5 0.5])
	
	vecIntercepts(intAnimal) = sStats.beta(1);
	vecSlopes(intAnimal) = sStats.beta(2);
end
vecY = polyval([mean(vecSlopes) mean(vecIntercepts)],vecX);
plot(vecX,vecY,'Color','b');

%mean over animals + raw data
hold off
ylim([0 1])

[h,p]=ttest(vecSlopes);
[corrected_p, h]=bonf_holm([p 1],0.05/2);
xlabel('Stimulus repetition')
ylabel('Response proportion')
title(sprintf('T-test slopes vs 0; p=%.3f',corrected_p(1)))

%RT over rep (C)
vecSlopes = nan(1,intAnimals);
vecIntercepts = nan(1,intAnimals);
subplot(2,2,3)
hold on
for intAnimal=1:intAnimals
	%resp C
	matBehav = cellHitOverTime{intAnimal,2};
	matRep = meshgrid(1:size(cellHitOverTime{intAnimal,2},2),1:size(cellHitOverTime{intAnimal,2},1));
	
	sStats=regstats(matBehav(:),matRep(:),'linear');
	vecX = [1 max(matRep(:))];
	vecY = polyval(sStats.beta([2 1]),vecX);
	plot(vecX,vecY,'Color',[0.5 0.5 0.5])
	
	vecIntercepts(intAnimal) = sStats.beta(1);
	vecSlopes(intAnimal) = sStats.beta(2);
end
vecY = polyval([mean(vecSlopes) mean(vecIntercepts)],vecX);
plot(vecX,vecY,'Color','b');
ylim([0 2])
%mean over animals + raw data
hold off

[h,p]=ttest(vecSlopes);
[corrected_p, h]=bonf_holm([p 1],0.05/2);
xlabel('Stimulus repetition')
ylabel('Reaction Time (s)')
title(sprintf('T-test slopes vs 0; p=%.3f',corrected_p(1)))


if boolSavePlots
	intFigCounter = intFigCounter + 1;
	drawnow;
	jFig = get(handle(gcf), 'JavaFrame');
	jFig.setMaximized(true);
	figure(gcf);
	drawnow;
	strFig = sprintf('Meta%d_behavior_repetition_dependency_raw',intFigCounter);
	export_fig([strFig '.tif']);
	export_fig([strFig '.pdf']);
end


%% dF/F over contrast for hit/miss and still/move
matContActNorm = zeros(size(matContAct));
for intAnimal=1:size(matContAct,3)
	matTempAct = matContAct(:,:,intAnimal);
	matContActNorm(:,:,intAnimal) = matTempAct;% ./ mean(matTempAct(:));
end
hActCon=figure;
for intAct=[0 1]
	if intAct == 0
		intMatOffset = 0;
		vecLimY = [-0.02 0.06];
	else
		intMatOffset = 8;
		vecLimY = [0.7 1.1];
	end
	for intMove=[0 1]
		if intMove == 0
			matHitAct = matContActNorm(:,1+intMatOffset,:);
			matMissAct = matContActNorm(:,2+intMatOffset,:);
			strMove = 'Still';
		else
			matHitAct = matContActNorm(:,3+intMatOffset,:);
			matMissAct = matContActNorm(:,4+intMatOffset,:);
			strMove = 'Move';
		end
		
		vecHitY = nanmean(matHitAct,3)';
		vecHitE = nanstd(matHitAct,[],3)'/sqrt(size(matContActNorm,3));
		vecMissY = nanmean(matMissAct,3)';
		vecMissE = nanstd(matMissAct,[],3)'/sqrt(size(matContActNorm,3));
		
		%get sig
		vecP = zeros(1,size(matContActNorm,1));
		for intC=1:size(matContActNorm,1)
			[h,p,ci] = ttest(matHitAct(intC,1,:),matMissAct(intC,1,:));
			
			%put in vector
			vecP(intC) = p;
		end
		
		%overall t-test
		matActHits = matHitAct(2:5,1,:);
		matActMisses = matMissAct(2:5,1,:);
		[h,dblP_All,ci] = ttest(matActHits(:),matActMisses(:));
		
		%pre-compute variables
		vecContrasts = [0.2 0.5 2 8 32 100];
		vecWindow = [1 length(vecContrasts)];
		vecWindowSelect = vecWindow(1):vecWindow(end);
		intWL = length(vecWindowSelect);
		vecWindowInv = intWL:-1:1;
		vecWindowPlotInv = vecWindow(end):-1:vecWindow(1);
		vecX = [vecContrasts(vecWindowSelect) vecContrasts(vecWindowPlotInv)];
		vecLineX = vecContrasts(vecWindowSelect);
		
		%get data
		subplot(2,2,intMove+1+intAct*2)
		for intResp=[0 1]
			if intResp == 1
				vecMeanTrace = vecHitY(vecWindowSelect);
				vecSE = vecHitE(vecWindowSelect);
				vecColorFill = [0.7 1 0.7];
				vecColorLine = [0 1 0];
			else
				vecColorLine = [1 0 0];
				vecColorFill = [1 0.7 0.7];
				vecMeanTrace = vecMissY(vecWindowSelect);
				vecSE = vecMissE(vecWindowSelect);
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
		set(gca,'XTick',vecLineX,'XTickLabel',[0 vecLineX(2:end)])
		dblFracMoved = mean(mean(matContAct(:,6,:)));
		title(sprintf('%s; Pop act during stim [frac moved: %.3f]; p-vals: 0: %.3f; 0.5: %.3f; 2: %.3f; 8: %.3f; 32: %.3f; 100: %.3f; Overall,p=%.3f',strMove,dblFracMoved,[vecP dblP_All]))
		grid on
		xlabel('Stimulus contrast (%)')
		ylabel('dF/F0 during stimulus')
		xlim(vecContrasts(vecWindow))
		%ylim(vecLimY)
		legend({'SEM','Miss','SEM','Hit'},'Location','Best')
		
	end
	vecUseC = 2:6;
	vecHitS = reshape(matContActNorm(vecUseC,1+intMatOffset,:),[1 numel(matContActNorm(vecUseC,1+intMatOffset,:))]);
	vecMissS = reshape(matContActNorm(vecUseC,2+intMatOffset,:),[1 numel(matContActNorm(vecUseC,2+intMatOffset,:))]);
	vecStill = [vecHitS vecMissS];
	
	vecHitM = reshape(matContActNorm(vecUseC,3+intMatOffset,:),[1 numel(matContActNorm(vecUseC,3+intMatOffset,:))]);
	vecMissM = reshape(matContActNorm(vecUseC,4+intMatOffset,:),[1 numel(matContActNorm(vecUseC,4+intMatOffset,:))]);
	vecMove = [vecHitM vecMissM];
	[h,p]=ttest(vecStill,vecMove);
	fprintf('Act %d; T-test still/move, p=%.3f\n',intAct,p)
end

%save
drawnow;
strFigTitle = [strMove '_act_over_contrasts_split_by_movement'];
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

%% normal plots of dF/F and heterogeneity
h=figure;
%plot dF/F over contrasts
for intActType=1:2
	subplot(2,2,intActType)
	if intActType==1
		strActType = 'dF/F';
		matHitAct = matContActNorm(:,7,:);
		matMissAct = matContActNorm(:,8,:);
		%vecLimY = [-0.02 0.06];
	else
		strActType = 'Heterogeneity';
		matHitAct = matContActNorm(:,13,:);
		matMissAct = matContActNorm(:,14,:);
		%for intAnimal=1:size(matHitAct,3)
		%	matTempActH = matHitAct(:,:,intAnimal);
		%	matTempActM = matMissAct(:,:,intAnimal);
		%	matHitAct(:,:,intAnimal) = matTempActH ./ mean([matTempActH(:);matTempActM(:)]);
		%	matMissAct(:,:,intAnimal) = matTempActM ./ mean([matTempActH(:);matTempActM(:)]);
		%end
		%vecLimY = [0.75 1];
	end
	
	vecHitY = nanmean(matHitAct,3)';
	vecHitE = nanstd(matHitAct,[],3)'/sqrt(size(matContActNorm,3));
	vecMissY = nanmean(matMissAct,3)';
	vecMissE = nanstd(matMissAct,[],3)'/sqrt(size(matContActNorm,3));
	
	%get sig
	vecP = zeros(1,size(matContActNorm,1));
	for intC=1:size(matContActNorm,1)
		[h,p,ci] = ttest(matHitAct(intC,1,:),matMissAct(intC,1,:));
		
		%put in vector
		vecP(intC) = p;
	end
	
	%overall t-test
	matActHits = matHitAct(2:5,1,:);
	matActMisses = matMissAct(2:5,1,:);
	[h,dblP_All,ci] = ttest(matActHits(:),matActMisses(:));
	
	%get data
	for intResp=[0 1]
		if intResp == 1
			vecMeanTrace = vecHitY(vecWindowSelect);
			vecSE = vecHitE(vecWindowSelect);
			vecColorFill = [0.7 1 0.7];
			vecColorLine = [0 1 0];
		else
			vecColorLine = [1 0 0];
			vecColorFill = [1 0.7 0.7];
			vecMeanTrace = vecMissY(vecWindowSelect);
			vecSE = vecMissE(vecWindowSelect);
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
	set(gca,'XTick',vecLineX,'XTickLabel',[0 vecLineX(2:end)])
	dblFracMoved = mean(mean(matContAct(:,6,:)));
	title(sprintf('Pop act during stim; p-vals: 0: %.3f; 0.5: %.3f; 2: %.3f; 8: %.3f; 32: %.3f; 100: %.3f; Overall,p=%.3f',[vecP dblP_All]))
	grid on
	xlabel('Stimulus contrast (%)');
	ylabel(sprintf('%s during stimulus',strActType));
	xlim(vecContrasts(vecWindow))
	%ylim(vecLimY);
	legend({'SEM','Miss','SEM','Hit'},'Location','Best')
end

drawnow;
strFigTitle = [strMove '_act_over_contrasts'];
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

%% quintile removal
%matQuintileRemoval = [6,6,2,8]; %contrast x quintile x hit/miss x animal
h=figure;
matQ_H = squeeze(nanmean(matQuintileRemoval(2:end,:,1,:),1));
matQ_M = squeeze(nanmean(matQuintileRemoval(2:end,:,2,:),1));

%t-tests
[h,p]=ttest(matQ_M',matQ_H')

%plot
errorbar(1:6,mean(matQ_M,2),std(matQ_M,[],2)/sqrt(size(matQ_M,2)),'r')
hold on
errorbar(1:6,mean(matQ_H,2),std(matQ_M,[],2)/sqrt(size(matQ_H,2)),'g')
hold off
set(gca,'XTick',1:6,'XTickLabel',{'0-20%','20-40%','40-60%','60-80%','80-100%','None'})
xlabel('Activity percentile of removed neurons')
ylabel('Heterogeneity')

drawnow;
strFigTitle = [strMove '_quintiles_removed'];
if boolSavePlots
	intFigCounter = intFigCounter + 1;
	drawnow;
	strFig = sprintf('Meta%d_%s_raw',intFigCounter,strFigTitle);
	export_fig([strFig '.tif']);
	export_fig([strFig '.pdf']);
end

%% decoding of stimulus presence
%decoded stim presence probability + 95% CI
hDecodeStimPresence = figure;
matMean = mean(matDetect,3); %[contrasts] x [hit/miss] x [animal]
matErr = std(matDetect,[],3)/sqrt(size(matDetect,3));


%plot
subplot(2,2,1)
vecWindow = 1:length(vecContrasts);
vecWindowSelect = vecWindow(1):vecWindow(end);
intWL = length(vecWindowSelect);
vecWindowInv = intWL:-1:1;
vecWindowPlotInv = vecWindow(end):-1:vecWindow(1);
vecX = [vecContrasts(vecWindowSelect) vecContrasts(vecWindowPlotInv)];
vecLineX = vecContrasts(vecWindowSelect);

for intDetect=[0 1]
	if intDetect == 1
		vecMeanTrace = matMean(:,1);
		vecMinTrace = matMean(:,1)-matErr(:,1);
		vecMaxTrace = matMean(:,1)+matErr(:,1);
		vecColorFill = [0.7 1.0 0.7];
		vecColorLine = [0 1 0];
	else
		vecMeanTrace = matMean(:,2);
		vecMinTrace = matMean(:,2)-matErr(:,2);
		vecMaxTrace = matMean(:,2)+matErr(:,2);
		vecColorLine = [1 0 0];
		vecColorFill = [1 0.7 0.7];
	end
	vecY = [vecMinTrace; vecMaxTrace(vecWindowInv)];
	
	%plot
	hold on
	fill(vecX,vecY,vecColorFill,'EdgeColor','none');
	plot(vecLineX,vecMeanTrace,'-','LineWidth',2,'Color',vecColorLine);
	hold off
end

%t-tests
vecP = nan(1,length(vecContrasts));
for intContrast=1:length(vecContrasts)
	[h,vecP(intContrast)] = ttest(matDetect(intContrast,1,:),matDetect(intContrast,2,:));
end
%overall t-test
matHits = matDetect(2:5,1,:);
matMisses = matDetect(2:5,2,:);
[h,dblP_All,ci] = ttest(matHits(:),matMisses(:));

strTitle = 'Stim Presence decoding; T-test p: ';
for intC=1:length(vecP)
	strTitle = [strTitle sprintf('C%.1f: ',vecLineX(intC)) sprintf('%.3f; ',vecP(intC))];
end
strTitle = [strTitle '; overall, p=' sprintf('%.3f',dblP_All)];

title(strTitle)


set(gca,'XScale','log','YScale','linear')
set(gca,'XTick',vecLineX,'XTickLabel',[0 vecLineX(2:end)])
grid on
xlabel('Contrast')
ylabel('Decoded stimulus presence')
xlim([min(vecContrasts(vecWindow)) max(vecContrasts(vecWindow))])
ylim([0 1])
legend({'StErr','Miss','StErr','Hit'},'Location','Best')


%% CALCULATE CORRELATION DECODER PERFORMANCE WITH BEHAVIORAL PERFORMANCE
% shuffle intermediate contrasts to get confidence interval; show higher
% decoder-behavior correlation for intermediate contrasts than chance
intIters = 1000;
matDecDetect = squeeze(mean(matDetect,2));
vecRealCorrs = nan(1,intAnimals);
matShuffCorrs = nan(intIters,intAnimals);
matDec = matDecDetect;%(2:(end-1),:);
matBeh = matBehavDetect;%(2:(end-1),:);
vecShuffMean = nan(1,intAnimals);
vecShuffUpper = nan(1,intAnimals);
vecShuffLower = nan(1,intAnimals);
for intAnimal=1:intAnimals
	%get real correlations
	vecRealCorrs(intAnimal) = corr(matDec(:,intAnimal),matBeh(:,intAnimal));
	
	%shuffle
	for intIter=1:intIters
		vecDec = matDec(randperm(size(matDec,1)),intAnimal);
		vecBeh = matBeh(randperm(size(matBeh,1)),intAnimal);
		matShuffCorrs(intIter,intAnimal) = corr(vecDec,vecBeh);
	end
	vecShuffMean(intAnimal) = mean(matShuffCorrs(:,intAnimal));
end
[h,p] = ttest(vecShuffMean,vecRealCorrs);

%plot
subplot(2,2,2)
errorbar(1.2,mean(vecShuffMean),-std(vecShuffMean)/sqrt(intAnimals),+std(vecShuffMean)/sqrt(intAnimals),'xb')
hold on
scatter(ones(1,intAnimals)*0.8,vecShuffMean,[],[0.5 0.5 0.5],'x');
errorbar(2.2,mean(vecRealCorrs),-std(vecRealCorrs)/sqrt(intAnimals),+std(vecRealCorrs)/sqrt(intAnimals),'xb')
scatter(ones(1,intAnimals)*1.8,vecRealCorrs,[],[0.5 0.5 0.5],'x');
plot([0.5 2.5],[0 0],'k--')
hold off
xlim([0 3])
ylim([-1 1])
ylabel('Similarity behavioral & decoding performance')
set(gca,'XTick',[1 2],'XTickLabel',{'Shuffled','Real'})
title(sprintf('Behavioral similarity decoder; p=%.3f',p))

%matDetectLowHighActHetRNR(intContrasts,intHighLow,intActTypes,intResp_NoResp);
%high/low act/act
for intActType = [1 2]
	if intActType == 1,strType='dF/F0';else strType='Het';end
	subplot(2,2,2+intActType)
	matData = squeeze(matDetectLowHighActHetRNR(:,:,intActType,:,:));
	
	for intLowHigh=[1 2]
		dblBlue = intLowHigh-1;
		matMean = squeeze(nanmean(matData(:,intLowHigh,:,:),4));
		matErr = squeeze(nanstd(matData(:,intLowHigh,:,:),[],4))/sqrt(intAnimals);
		for intResp=[1 2]
			if intResp == 2
				vecMeanTrace = matMean(:,1);
				vecMinTrace = matMean(:,1)-matErr(:,1);
				vecMaxTrace = matMean(:,1)+matErr(:,1);
				vecColorFill = [0.7 1.0 dblBlue];
				vecColorLine = [0 1 dblBlue];
			else
				vecMeanTrace = matMean(:,2);
				vecMinTrace = matMean(:,2)-matErr(:,2);
				vecMaxTrace = matMean(:,2)+matErr(:,2);
				vecColorLine = [1 0 dblBlue];
				vecColorFill = [1 0.7 dblBlue];
			end
			vecY = [vecMinTrace; vecMaxTrace(vecWindowInv)];
			
			%plot
			hold on
			fill(vecX,vecY,vecColorFill,'EdgeColor','none');
			plot(vecLineX,vecMeanTrace,'-','LineWidth',2,'Color',vecColorLine);
			hold off
		end
	end
	%{
	%t-tests
	vecP = nan(1,length(vecContrasts));
	for intContrast=1:length(vecContrasts)
		[h,vecP(intContrast)] = ttest(matData(intContrast,1,:),matData(intContrast,2,:));
	end
	%overall t-test
	matHits = matData(2:5,1,:);
	matMisses = matData(2:5,2,:);
	[h,dblP_All,ci] = ttest(matHits(:),matMisses(:));
	
	strTitle = sprintf('Presence decoding %s; T-test p: ',strType);
	for intC=1:length(vecP)
		strTitle = [strTitle sprintf('C%.1f: ',vecLineX(intC)) sprintf('%.3f; ',vecP(intC))];
	end
	strTitle = [strTitle '; overall, p=' sprintf('%.3f',dblP_All)];
	
	title(strTitle)
	
	%}
	set(gca,'XScale','log','YScale','linear')
	set(gca,'XTick',vecLineX,'XTickLabel',[0 vecLineX(2:end)])
	grid on
	xlabel('Contrast')
	ylabel('Decoded stimulus presence')
	xlim([min(vecContrasts(vecWindow)) max(vecContrasts(vecWindow))])
	ylim([0 1])
	legend({'StErr','Low NoResp','StErr','Low Resp','StErr','High NoResp','StErr','High Resp'},'Location','Best')
end


%matDetectLowHighActHetRNR(intContrasts,intHighLow,intActTypes,intResp_NoResp);
cellData{1} = squeeze(mean(matDetectLowHighActHetRNR(:,2,1,:) - matDetectLowHighActHetRNR(:,1,1,:),1));
cellData{2} = squeeze(mean(matDetectLowHighActHetRNR(:,2,2,:) - matDetectLowHighActHetRNR(:,1,2,:),1));
for intActType=[1 2]
	if intActType == 1,strType='dF/F0';else strType='Het';end
	subplot(2,2,2+intActType)

	scatter(ones(size(cellData{intActType})),cellData{intActType})
	hold on
	errorbar(2,mean(cellData{intActType}),std(cellData{intActType})/sqrt(intAnimals),'x')
	hold off
	ylim([0 0.6+eps])
	xlim([0 3])
end
[h,p]=ttest(cellData{1},cellData{2})


%save fig
strFigTitle = [strMove '_decodeSplitStimPresence'];
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


%% plot resp type predictability
%{
vecHetMiss = [];
vecHetHit = [];
vecActMiss = [];
vecActHit = [];

%which data points?
vecWindowSecs = [-3 5];
vecWindow = round(vecWindowSecs*25.4);
vecLineX = (vecWindow(1):vecWindow(end))/25.4;
vecSelectSecs = [-3 0];
vecSelect = round(vecSelectSecs*25.4);
intStart = find(vecLineX>=vecSelectSecs(1),1);
intStop = find(vecLineX>=vecSelectSecs(2),1);

%get data
for intAnimal=1:size(cellHetTime,2)
	matTemp = [];
	for intC=2:size(cellHetTime,1)
		matTemp = [matTemp;cellHetTime{intC,intAnimal}{1};cellHetTime{intC,intAnimal}{2};cellHetTime{intC,intAnimal}{3}];
	end
	dblMean = mean(matTemp(:)); %mean over all values of animal
	
	%apply normalization
	for intC=2:size(cellHetTime,1)
		cellHetTime{intC,intAnimal}{1} = cellHetTime{intC,intAnimal}{1} ./ dblMean;
		cellHetTime{intC,intAnimal}{2} = cellHetTime{intC,intAnimal}{2} ./ dblMean;
		cellHetTime{intC,intAnimal}{3} = cellHetTime{intC,intAnimal}{3} ./ dblMean;
	end
	
	%prep data
	matTempHetM = [];
	matTempHetH = [];
	matTempActM = [];
	matTempActH = [];
	for intC=2:size(cellHetTime,1)
		matTempHetM = [matTempHetM;cellHetTime{intC,intAnimal}{1}];
		matTempHetH = [matTempHetH;cellHetTime{intC,intAnimal}{2}];
		matTempActM = [matTempActM;cellActTime{intC,intAnimal}{1}];
		matTempActH = [matTempActH;cellActTime{intC,intAnimal}{2}];
	end
	vecHetMiss(intAnimal) = mean(mean(matTempHetM(:,intStart:intStop),1));
	vecHetHit(intAnimal) = mean(mean(matTempHetH(:,intStart:intStop),1));
	vecActMiss(intAnimal) = mean(mean(matTempActM(:,intStart:intStop),1));
	vecActHit(intAnimal) = mean(mean(matTempActH(:,intStart:intStop),1));
end
% check separability
%{
calculate per point intra-cluster vs inter-cluster difference (ratio) for
dF/F and for heterogeneity dimension; will give distribution of
separability for each cluster; then perform ANOVA/ttests for
between-cluster distribution comparison
%}

%get cluster means
dblMeanActMiss = mean(vecActMiss);
dblMeanActHit = mean(vecActHit);
dblMeanHetMiss = mean(vecHetMiss);
dblMeanHetHit = mean(vecHetHit);
intAnimals = length(vecActMiss);

%separability; inter vs intra cluster distance
vecSepRatioActMiss_w_Hit = abs(vecActMiss - dblMeanActHit) ./ (abs(vecActMiss - dblMeanActMiss) + abs(vecActMiss - dblMeanActHit));%inter / (inter+intra) distance
vecSepRatioActHit_w_Miss = abs(vecActHit - dblMeanActMiss) ./ (abs(vecActHit - dblMeanActHit) + abs(vecActHit - dblMeanActMiss));

vecSepRatioHetMiss_w_Hit = abs(vecHetMiss - dblMeanHetHit) ./ (abs(vecHetMiss - dblMeanHetMiss) + abs(vecHetMiss - dblMeanHetHit));
vecSepRatioHetHit_w_Miss = abs(vecHetHit - dblMeanHetMiss) ./ (abs(vecHetHit - dblMeanHetHit) + abs(vecHetHit - dblMeanHetMiss));

%act agg
vecActHitMiss = [vecSepRatioActMiss_w_Hit vecSepRatioActHit_w_Miss];

%het agg
vecHetHitMiss = [vecSepRatioHetMiss_w_Hit vecSepRatioHetHit_w_Miss];

%fig
hFigRTP = figure;

%plot
subplot(2,2,1)
scatter(vecHetMiss,vecActMiss,'rx')
hold on
scatter(vecHetHit,vecActHit,'gx')
hold off
title('Stim; Response type predictability')
xlabel('Heterogeneity preceding stimulus')
ylabel('dF/F preceding stimulus')

subplot(2,2,2)

dblOffset = 0.1;
plot([0.5 1.5],[0.5 0.5],'k--')
%plot act
hold on
errorbar(1-dblOffset,mean(vecActHitMiss),std(vecActHitMiss)/sqrt(length(vecActHitMiss)),'bx')
%plot het
errorbar(1+dblOffset,mean(vecHetHitMiss),std(vecHetHitMiss)/sqrt(length(vecHetHitMiss)),'kx')
hold off
set(gca,'XTick',1,'XTickLabel',{'H-M'})
ylabel(sprintf('Separability \n(inter/(intra+inter) cluster distance)'))
[h,pActHM] = ttest(vecActHitMiss,0.5);
[h,pHetHM] = ttest(vecHetHitMiss,0.5);
[h,pHM] = ttest(vecActHitMiss,vecHetHitMiss);
title(sprintf(['Stim; Sep act; HM,p=%.3f;\n'...
	'Sep het; HM,p=%.3f;\n',...
	'Diff act-het; HM,p=%.3f'],...
	pActHM,...
	pHetHM,...
	pHM));
ylim([0 1])
xlim([0.5 3.5])
legend('dF/F','Heterogeneity')

if boolSavePlots
	strFigTitle = 'predict_hitmiss';
	intFigCounter = intFigCounter + 1;
	drawnow;
	strFig = sprintf('Meta%d_%s_raw',intFigCounter,strFigTitle);
	export_fig([strFig '.tif']);
	export_fig([strFig '.pdf']);
end
%}
%% clean up
cd(strOldDir);