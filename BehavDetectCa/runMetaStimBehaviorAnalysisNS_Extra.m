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
boolUseNeuropilSubtraction = false;
strFigDir = 'D:\Data\Results\stimdetection\metaNS_Supp2';
strDataDir = 'D:\Data\Results\stimdetection';
cellInclude = {'20140207','20140314','20140425','20140507','20140530','20140604','20140711','20140715'};%exclude:20140430
boolSavePlots = true;
intFigCounter = 0;

%% pre-allocate
intAnimal = 0;
vecAnimalSource = [];
vecSubPopulationSource = [];
matCohen = [];
matCohenAll = [];
matDriftZLicking = [];
matHitMissROC = [];
cellRTDependency = [];
matHitMissCohensD = [];
matDetect = [];
cellStimDetect = [];
matITC = [];
matITC_Shuffled = [];
cellRTDependency2 = [];
matHitMissCohensD2 = [];
		
%% create aggregate data structures
sDir = dir(strDataDir);
strOldDir = cd(strDataDir);
intCounterF1 = 0;
intCounterF2 = 0;
intAnimZ = 0;

if boolUseNeuropilSubtraction
	strNPS = 'NPS';
else
	strNPS = '';
end

for intFile=1:numel(sDir)
	%check if data part 1
	strFile = sDir(intFile).name;
	strRec = getFlankedBy(strFile,['data_aggregateNS_Extra_' strNPS],'_');
	intRec = find(strcmp(strRec,cellInclude),1);
	if ~isempty(intRec)
		%load file
		load(strFile);
		intCounterF1 = intCounterF1 + 1;
		
		if size(cellSaveSeparability,1) > 1
			%cellSaveCohenBlocks{1,1} = [cellSaveCohenBlocks{1,1} cellSaveCohenBlocks{2,1}];
			%cellSaveCohenBlocks{1,2} = mean(cat(2,cellSaveCohenBlocks{1,2},cellSaveCohenBlocks{2,2}),2);
			
			if exist('cellSaveLickZDrift','var')
				for intType=1:4
					cellSaveLickZDrift{1,intType} = mean(cat(3,cellSaveLickZDrift{1,intType},cellSaveLickZDrift{2,intType}),3);
				end
				cellSaveLickZDrift(2,:) = [];
			end
			
			cellSaveSeparability = mean(cellSaveSeparability,1);
			
			
			%distros
			cellSaveMultiDimDistanceDistros = mean(cellSaveMultiDimDistanceDistros,3);
			
			%hit/miss decoding
			cellSaveHitMissDecoding = mean(cellSaveHitMissDecoding,2);
			%{
			%RT dep
			for intEl=1:size(cellSaveRTDependency,1)
				cellSaveRTDependency{intEl,1} = [cellSaveRTDependency{intEl,1}; cellSaveRTDependency{intEl,2}];
			end
			cellSaveRTDependency(:,2) = [];
			
			%hit-miss cohen's D
			cellSaveCohensD_HitMiss{1} = mean(cat(3,cellSaveCohensD_HitMiss{1},cellSaveCohensD_HitMiss{2}),3);
			cellSaveCohensD_HitMiss(2) = [];
			%}
			
			%RT dep 2
			for intEl=1:size(cellSaveRTDependency2,1)
				cellSaveRTDependency2{intEl,1} = [cellSaveRTDependency2{intEl,1}; cellSaveRTDependency2{intEl,2}];
			end
			cellSaveRTDependency2(:,2) = [];
			
			%hit-miss cohen's D
			cellSaveCohensD_HitMiss2{1} = mean(cat(3,cellSaveCohensD_HitMiss2{1},cellSaveCohensD_HitMiss2{2}),3);
			cellSaveCohensD_HitMiss2(2) = [];
			
			
			%stim presence decoding
			cellSaveMatDetect{1,1} = mean(cat(4,cellSaveMatDetect{1,1},cellSaveMatDetect{2,1}),4); %[contrasts] x [hit/miss]
			for intEl=2:size(cellSaveMatDetect,2)
				cellSaveMatDetect{1,intEl} = [cellSaveMatDetect{1,intEl} cellSaveMatDetect{2,intEl}];
			end
			cellSaveMatDetect(2,:) = [];
			
			%inter-trial correlations
			%for intEl=1:size(cellSaveITC,2)
			%	cellSaveITC{1,1} = mean(cat(3,cellSaveITC{1,1},cellSaveITC{2,1}),3);
			%end
			%cellSaveITC(2,:) = [];
		end
		
		%assign data
		%matCohen = [matCohen cellSaveCohenBlocks{1,1}];
		%matCohenAll = [matCohenAll cellSaveCohenBlocks{1,2}];
		if exist('cellSaveLickZDrift','var')
			intAnimZ = intAnimZ + 1;
			for intType=1:4
				matDriftZLicking(intAnimZ,intType,:) = cellSaveLickZDrift{1,1};
			end
			clear cellSaveLickZDrift;
		end
		matHitMissROC(intCounterF1,:) = cellSaveSeparability;
		
		%matrix multidimensional distance distributions
		matMDDD(:,:,intCounterF1) = cellSaveMultiDimDistanceDistros;
		
		%hit/miss decoding accuracy
		matHMD(:,intCounterF1) = cellSaveHitMissDecoding; %raw no-mean no-het no-het-mean
		
		%RT dependency
		%cellRTDependency = cat(2,cellRTDependency,cellSaveRTDependency);
		
		%hit/miss cohen's d
		%matHitMissCohensD = cat(3,matHitMissCohensD,cellSaveCohensD_HitMiss{1});
		
		%RT dependency
		cellRTDependency2 = cat(2,cellRTDependency2,cellSaveRTDependency2);
		
		%hit/miss cohen's d
		matHitMissCohensD2 = cat(3,matHitMissCohensD2,cellSaveCohensD_HitMiss2{1});
		
		%stim presence decoding
		matDetect = cat(4,matDetect,cellSaveMatDetect{1,1});
		cellStimDetect{1,intCounterF1} = cellSaveMatDetect{1,2};
		cellStimDetect{2,intCounterF1} = cellSaveMatDetect{1,3};
		cellStimDetect{3,intCounterF1} = cellSaveMatDetect{1,4};
		
		%inter-trial correlations
		%matITC = cat(3,matITC,cellSaveITC{1,1});
		%matITC_Shuffled = cat(3,matITC_Shuffled,cellSaveITC{1,2});
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

%% het-type comparison
% reaction-time dependency
vecX = [0 3];
hRTDependenceFig = figure;
hMeasuresR2 = figure;
intActTypes = size(cellRTDependency2,1)-1;
matR2 = nan(intActTypes,intAnimals);
vecPs = nan(1,size(cellRTDependency2,1)-1);
for intActType=1:intActTypes
	figure(hRTDependenceFig);
	if intActType == 1
		%heterogeneity
		strLabelY = 'dF/F0';
		strC = 'k';
	elseif intActType == 2
		%z-scored act
		strLabelY = 'Heterogeneity';
		strC = 'r';
	elseif intActType == 3
		%dF/F act
		strLabelY = 'Heterogeneity MD';
		strC = 'm';
	elseif intActType == 4
		%variance
		strLabelY = 'Heterogeneity MD non-z';
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
		vecAggRTs = [vecAggRTs cellRTDependency2{1,intPopulation}'];
		vecAggAct = [vecAggAct cellRTDependency2{intActType+1,intPopulation}'];
		
		%do regressions per population
		sStatsC=regstats(cellRTDependency2{intActType+1,intPopulation}',cellRTDependency2{1,intPopulation},'linear');
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
	strFig = sprintf(['Meta%d_HetTypeComp_RTDependency_raw' strNPS],intFigCounter);
	export_fig([strFig '.tif']);
	export_fig([strFig '.pdf']);
	
	figure(hMeasuresR2);
	intFigCounter = intFigCounter + 1;
	drawnow;
	strFig = sprintf(['Meta%d_HetTypeComp_RTDependency_R2_raw' strNPS],intFigCounter);
	export_fig([strFig '.tif']);
	export_fig([strFig '.pdf']);
end


%% hit/miss effect sizes for same metrics as RT
matCohenTest2 = squeeze(mean(matHitMissCohensD2(:,2:5,:),2));
matCohenProbe2 = squeeze(mean(matHitMissCohensD2(:,[1 6],:),2));
%matCohenTest(:,end) = [];
vecCohenMean2_HM = nanmean(matCohenTest2,2);
vecCohenSEM2_HM = nanstd(matCohenTest2,[],2)./sqrt(intAnimals);


figure
hold on
intMeasures = length(vecCohenMean2_HM);
vecX = [1:intMeasures]/intMeasures;
errorbar(vecX,vecCohenMean2_HM,vecCohenSEM2_HM,'bx')
hold off
cellMeasures = {'dF/F0','Het','Het MD','Het MD non-z'};
xlim([-0.1 1.1])
ylabel('Cohen''s D hit-miss difference')
set(gca,'xtick',vecX,'xticklabel',cellMeasures)

%t-tests
[h,pAllvs0]=ttest(matCohenTest2,[],[],[],2);
[h,pHvsAll]=ttest(matCohenTest2,repmat(matCohenTest2(1,:),[intMeasures 1]),[],[],2);
pAllvs0
pHvsAll

[d,d,pAllvs0_adj] = fdr_bh(pAllvs0);
[d,d,pHvsAll_adj] = fdr_bh(pHvsAll);

drawnow;
strFigTitle = ['hit_miss_HetComp_cohens_d' strNPS];
if boolSavePlots
	intFigCounter = intFigCounter + 1;
	drawnow;
	strFig = sprintf('Meta%d_%s_raw',intFigCounter,strFigTitle);
	export_fig([strFig '.tif']);
	export_fig([strFig '.pdf']);
end

%{
%% reaction-time dependency
%
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
	strFig = sprintf(['Meta%d_activationdissimilarity_RTDependency_raw' strNPS],intFigCounter);
	export_fig([strFig '.tif']);
	export_fig([strFig '.pdf']);
	
	figure(hMeasuresR2);
	intFigCounter = intFigCounter + 1;
	drawnow;
	strFig = sprintf(['Meta%d_RTDependency_R2_raw' strNPS],intFigCounter);
	export_fig([strFig '.tif']);
	export_fig([strFig '.pdf']);
end
%}
%{
%% hit/miss effect sizes for same metrics as RT
matCohenTest = squeeze(mean(matHitMissCohensD(:,2:5,:),2));
matCohenProbe = squeeze(mean(matHitMissCohensD(:,[1 6],:),2));
%matCohenTest(:,end) = [];
matCohenTest([5 end],:) = -matCohenTest([5 end],:,:);
matCohenMean_HM = nanmean(matCohenTest,2);
matCohenSEM_HM = nanstd(matCohenTest,[],2)./sqrt(intAnimals);


figure
hold on
vecX = [1:11]/11;
errorbar(vecX,matCohenMean_HM,matCohenSEM_HM,'bx')
hold off
cellMeasures = {'Het','Z all','Act all','var','sparse','P-like','SD P-like','PP act','PP Z','Mean SWC','SD SWC'};
xlim([-0.1 1.1])
ylabel('Cohen''s D hit-miss difference')
set(gca,'xtick',vecX,'xticklabel',cellMeasures)

%t-tests
[h,pAllvs0]=ttest(matCohenTest,[],[],[],2);
[h,pHvsAll]=ttest(matCohenTest,repmat(matCohenTest(1,:),[11 1]),[],[],2);
pAllvs0
pHvsAll

%[d,d,pAllvs0_adj] = fdr_bh(pAllvs0)
%[d,d,pHvsAll_adj] = fdr_bh(pHvsAll)

drawnow;
strFigTitle = ['hit_miss_cohens_d' strNPS];
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
%}
%% stim presence decoding vs behavior
%cellStimDetect %[decoded/behavior/contrast] x [animals]
%matDetect; %[up/mean/low] x [contrasts] x [hit/miss/all/behavior] x [animals]

vecDecoded = cell2mat(cellStimDetect(1,:));
vecBehavior = cell2mat(cellStimDetect(2,:));

[table,chi2,p,labels] = crosstab(vecDecoded,vecBehavior);


%% inter-trial correlations + shuffle control
%{
%for whole pop or only preferred pop
matMean = median(matITC,3);
matErr = std(matITC,[],3)/sqrt(intAnimals);
matMeanS = median(matITC_Shuffled,3);
matErrS = std(matITC_Shuffled,[],3)/sqrt(intAnimals);


hPopCorrPlot = figure;

matErrLow = (matMean-quantile(matITC,normcdf(-1,0,1),3))/sqrt(intAnimals);
matErrHigh = (quantile(matITC,normcdf(1,0,1),3)-matMean)/sqrt(intAnimals);

matErrLowS = (matMeanS-quantile(matITC_Shuffled,normcdf(-1,0,1),3))/sqrt(intAnimals);
matErrHighS = (quantile(matITC_Shuffled,normcdf(1,0,1),3)-matMeanS)/sqrt(intAnimals);

subplot(2,2,1)
cellColor={'r','m','g'};
hold on;
for intResp=1:3
	errorbar(intResp,matMean(1,intResp),matErr(1,intResp),['x' cellColor{intResp}]);
	errorbar(intResp,matMeanS(1,intResp),matErrS(1,intResp),['xk']);
end
hold off
title('Non-preferred population');
cellLabels = {'Miss','Slow','Fast'};
set(gca,'xtick',1:3,'xticklabel',cellLabels)
ylabel('Inter-trial correlation');
ylim([0 0.12]);

subplot(2,2,2)
hold on;
for intResp=1:3
	errorbar(intResp,matMean(2,intResp),matErr(2,intResp),['x' cellColor{intResp}]);
	errorbar(intResp,matMeanS(2,intResp),matErrS(2,intResp),['xk']);
end
hold off
title('Preferred population');
set(gca,'xtick',1:3,'xticklabel',cellLabels)
ylabel('Inter-trial correlation');
ylim([0 0.12]);

[h,dblP_NP_MS] = ttest(matITC(1,1,:),matITC(1,2,:));
[h,dblP_NP_MF] =  ttest(matITC(1,1,:),matITC(1,3,:));
[h,dblP_NP_SF] =  ttest(matITC(1,2,:),matITC(1,3,:));

[h,dblP_P_MS] = ttest(matITC(2,1,:),matITC(2,2,:));
[h,dblP_P_MF] =  ttest(matITC(2,1,:),matITC(2,3,:));
[h,dblP_P_SF] =  ttest(matITC(2,2,:),matITC(2,3,:));

[h,pNPM]=ttest(matITC(1,1,:));
[h,pNPS]=ttest(matITC(1,2,:));
[h,pNPF]=ttest(matITC(1,3,:));
[h,pPM]=ttest(matITC(2,1,:));
[h,pPS]=ttest(matITC(2,2,:));
[h,pPF]=ttest(matITC(2,3,:));

[h,pNPM_S]=ttest(matITC(1,1,:),matITC_Shuffled(1,1,:));
[h,pNPS_S]=ttest(matITC(1,2,:),matITC_Shuffled(1,2,:));
[h,pNPF_S]=ttest(matITC(1,3,:),matITC_Shuffled(1,3,:));
[h,pPM_S]=ttest(matITC(2,1,:),matITC_Shuffled(2,1,:));
[h,pPS_S]=ttest(matITC(2,2,:),matITC_Shuffled(2,2,:));
[h,pPF_S]=ttest(matITC(2,3,:),matITC_Shuffled(2,3,:));

vecP_adjP = bonf_holm([pPM pPS pPF]);
vecP_adjNP = bonf_holm([pNPM pNPS pNPF]);

vecP_adjShuffled = bonf_holm([pNPM_S pNPS_S pNPF_S pPM_S pPS_S pPF_S]);


if boolSavePlots
	intFigCounter = intFigCounter + 1;
	figure(hPopCorrPlot);
	drawnow;
	jFig = get(handle(hPopCorrPlot), 'JavaFrame');
	jFig.setMaximized(true);
	figure(hPopCorrPlot);
	drawnow;
	strFig = sprintf('Meta%d_intertrialcorr_%s_pop%d_raw',intFigCounter);
	export_fig([strFig '.tif']);
	export_fig([strFig '.pdf']);
end
%}
%% distributions multdimensional inter-trial distance for all/hit/miss trials
%cellSaveMultiDimDistanceDistros(1,:,intPopulation) = vecDistroAll;
%cellSaveMultiDimDistanceDistros(2,:,intPopulation) = vecDistroHit;
%cellSaveMultiDimDistanceDistros(3,:,intPopulation) = vecDistroMiss;
%cellSaveMultiDimDistanceDistros(4,:,intPopulation) = vecMirrorDistroAll;
%cellSaveMultiDimDistanceDistros(5,:,intPopulation) = vecMirrorDistroHit;
%cellSaveMultiDimDistanceDistros(6,:,intPopulation) = vecMirrorDistroMiss;
%matMDDD

matMeanDistros = mean(matMDDD,3);
matSEMDistros = std(matMDDD,[],3)/sqrt(intAnimals);
intType = 2;
if intType == 1
	%binning vector
	dblStep = 0.02; %dF/F: 0.02
	vecBins = 0:dblStep:1.4;
	strLabelX = 'Distance (dF/F0)';
else
	%binning vector
	dblStep = 0.4; %dF/F: 0.02
	vecBins = 0:dblStep:30;
	strLabelX = 'Distance (z-score units)';
end

vecCoM_All = nan(1,8);
vecCoM_Hit = nan(1,8);
vecCoM_Miss = nan(1,8);
for intAnimal=1:intAnimals
	vecCoM_All(intAnimal) = calcCenterOfMass(matMDDD(1,:,intAnimal),2)*dblStep;
	vecCoM_Hit(intAnimal) = calcCenterOfMass(matMDDD(2,:,intAnimal),2)*dblStep;
	vecCoM_Miss(intAnimal) = calcCenterOfMass(matMDDD(3,:,intAnimal),2)*dblStep;
end
vecAH = vecCoM_All - vecCoM_Hit;
vecAM = vecCoM_All - vecCoM_Miss;
vecHM = vecCoM_Hit - vecCoM_Miss;

%hit-miss difference original distribution
[h,pAH]=ttest(vecAH);
[h,pAM]=ttest(vecAM);
[h,pHM]=ttest(vecHM);
[d,d,vecP_adj] = fdr_bh([pAH pAM pHM]);

figure
subplot(2,2,1);
stairs(vecBins,matMeanDistros(1,:),'-','Color',[0 0 0]);
hold on
stairs(vecBins,matMeanDistros(2,:),'-','Color',[0 1 0]);
stairs(vecBins,matMeanDistros(3,:),'-','Color',[1 0 0]);
hold off

subplot(2,2,2);
%effect of mirroring
matHitMirrored = matMDDD(2,:,:) + matMDDD(5,:,:);
matMissMirrored = matMDDD(3,:,:) + matMDDD(6,:,:);
vecCoM_Hit_Mirrored = nan(1,8);
vecCoM_Miss_Mirrored = nan(1,8);
for intAnimal=1:intAnimals
	vecCoM_Hit_Mirrored(intAnimal) = calcCenterOfMass(matHitMirrored(1,:,intAnimal),2)*dblStep;
	vecCoM_Miss_Mirrored(intAnimal) = calcCenterOfMass(matMissMirrored(1,:,intAnimal),2)*dblStep;
end
vecHiMiDi = vecCoM_Hit_Mirrored - vecCoM_Hit;
vecMiMiDi = vecCoM_Miss_Mirrored - vecCoM_Miss;

[h,pHiMi]=ttest(vecHiMiDi);
[h,pMiMi]=ttest(vecMiMiDi);
[vecP_adj_mi] = bonf_holm([pHiMi pMiMi]);

[h,pHM_HiMi]=ttest(vecHM,vecHiMiDi);
[h,pHM_MiMi]=ttest(vecHM,vecMiMiDi);
[vecP_adj_mi] = bonf_holm([pHM_HiMi pHM_MiMi]);

errorbar(1:3,[mean(vecHM) mean(vecHiMiDi) mean(vecMiMiDi)],[std(vecHM) std(vecHiMiDi) std(vecMiMiDi)]/sqrt(intAnimals),'x')
ylim([-0.5 0.5])
ylabel('d(CoM) of pairwise interpoint distance (dF/F0)')
set(gca,'xtick',1:3,'xticklabel',{'Hit vs. Miss','Hit vs. Mirrored','Miss vs. Mirrored'})

subplot(2,2,3)
errorbar(1:3,[mean(vecHM) mean(vecHiMiDi) mean(vecMiMiDi)],[std(vecHM) std(vecHiMiDi) std(vecMiMiDi)]/sqrt(intAnimals),'x')
ylim([-0.001 0.001])
ylabel('d(CoM) of pairwise interpoint distance (dF/F0)')
set(gca,'xtick',1:3,'xticklabel',{'Hit vs. Miss','Hit vs. Mirrored','Miss vs. Mirrored'})

%{
%errorbar(vecBins,matMeanDistros(1,:),matSEMDistros(1,:),'-','Color',[0 0 0]);
%hold on
%errorbar(vecBins,matMeanDistros(2,:),matSEMDistros(2,:),'-','Color',[0 1 0]);
%errorbar(vecBins,matMeanDistros(3,:),matSEMDistros(3,:),'-','Color',[1 0 0]);
%hold off
xlabel(strLabelX);
ylabel('Mean normalized count of trial-pairs');
%title(sprintf('Mean distances, all=%.3f, hit=%.3f, miss=%.3f',mean(matPairwiseDistances(matSelectTrials)),mean(matPairwiseDistances(matSelectHitTrials)),mean(matPairwiseDistances(matSelectMissTrials))));
	
subplot(2,2,2);
errorfill(vecBins,matMeanDistros(4,:),matSEMDistros(4,:),[0 0 0],[0.5 0.5 0.5]);
xlabel(strLabelX);
ylabel('Normalized count of trial-pairs');
title(sprintf('Mean distance All trials; Mirrored-Raw'));%, data=%.3f, mirrored=%.3f',mean(matPairwiseDistances(matSelectTrials)),mean(matPairwiseDistances(matSelectHitTrials)),mean(matPairwiseDistances(matSelectMissTrials))));

subplot(2,2,3);
errorfill(vecBins,matMeanDistros(5,:),matSEMDistros(5,:),[0 1 0],[0.5 1 0.5]);
xlabel(strLabelX);
ylabel('Normalized count of trial-pairs');
title(sprintf('Mean distance Hit trials; Mirrored-Raw'));%, data=%.3f, mirrored=%.3f',mean(matPairwiseDistances(matSelectTrials)),mean(matPairwiseDistances(matSelectHitTrials)),mean(matPairwiseDistances(matSelectMissTrials))));

subplot(2,2,4);
errorfill(vecBins,matMeanDistros(6,:),matSEMDistros(6,:),[1 0 0],[1 0.5 0.5]);
xlabel(strLabelX);
ylabel('Normalized count of trial-pairs');
title(sprintf('Mean distance Miss trials; Mirrored-Raw'));%, data=%.3f, mirrored=%.3f',mean(matPairwiseDistances(matSelectTrials)),mean(matPairwiseDistances(matSelectHitTrials)),mean(matPairwiseDistances(matSelectMissTrials))));
%}
%save
drawnow;
strFigTitle = ['multidim_distance_all-hit-miss_trials' strNPS];
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
%}
%% hit/miss decoding
%cellSaveHitMissDecoding(:,intPopulation) = [dblPerformanceRawML dblPerformanceNoMeanML dblPerformanceNoHetML dblPerformanceNoHetMeanML];
%matHMD
vecMeanHMD = mean(matHMD,2);
vecSEMHMD = std(matHMD,[],2)/sqrt(intAnimals);

%t-tests
[h,pOM]=ttest(matHMD(1,:),matHMD(2,:));
[h,pOH]=ttest(matHMD(1,:),matHMD(3,:));
[h,pOB]=ttest(matHMD(1,:),matHMD(4,:));
[h,pMH]=ttest(matHMD(2,:),matHMD(3,:));
[h,pMB]=ttest(matHMD(2,:),matHMD(4,:));
[h,pHB]=ttest(matHMD(3,:),matHMD(4,:));

figure
subplot(2,2,1)
errorbar(1:4,vecMeanHMD,vecSEMHMD,'xb');
ylim([0.5 0.7])
set(gca,'xtick',1:4,'xticklabel',{'Original Data','Mean removed','Het removed','Both removed'});
ylabel('Hit/miss decoding accuracy');
title(sprintf('T-tests; OM,p=%.3f; OH,p=%.3f; OB,p=%.3f; MH,p=%.3f; MB,p=%.3f; HB,p=%.3f',[pOM pOH pOB pMH pMB pHB]))


%% single-trial hit/miss linear separability ROC
subplot(2,2,2)
%matHitMissROC(8,:) = [];
[h,pH]=ttest(matHitMissROC(:,1),0.5);
[h,pA]=ttest(matHitMissROC(:,2),0.5);
[h,pHA]=ttest(matHitMissROC(:,1),matHitMissROC(:,2));

errorbar([1 2],[mean(matHitMissROC(:,1)) mean(matHitMissROC(:,2))],[std(matHitMissROC(:,1)) std(matHitMissROC(:,2))]/sqrt(intAnimals),'x');
ylim([0.5 0.65])
set(gca,'xtick',[1 2],'xticklabel',{'Heterogeneity','mean dF/F0'});
title(sprintf('Single-trial hit/miss separability; het vs 0.5; p=%.3f; dF/F vs 0.5; p=%.3f; het vs dF/F; p=%.3f',pH,pA,pHA));
ylabel('Area under ROC curve')

%save
drawnow;
strFigTitle = ['hit-miss_decoding' strNPS];
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


%% cohens d
%{
figure
subplot(2,2,1)
vecBehavRepD_p = matCohen(1,:);
vecCohenRepAct = matCohen(2,:);
vecCohenRepHet = matCohen(3,:);

scatter(vecBehavRepD_p,vecCohenRepHet,'kx');
hold on
scatter(vecBehavRepD_p,vecCohenRepAct,'rx');
hold off;
drawnow
xlabel('Behavioral performance (P(hit) - P(FA))')
ylabel('Cohen''s d hit-miss')
legend('Heterogeneity','dF/F0')
fixfig

subplot(2,2,2)
vecBehavRepD_p = matCohenAll(1,:);
vecCohenRepAct = matCohenAll(2,:);
vecCohenRepHet = matCohenAll(3,:);

scatter(vecBehavRepD_p,vecCohenRepHet,'kx');
hold on
scatter(vecBehavRepD_p,vecCohenRepAct,'rx');
hold off;
drawnow
xlabel('Behavioral performance (P(hit) - P(FA)')
ylabel('Cohen''s d hit-miss')
legend('Heterogeneity','dF/F0')
fixfig

%save
drawnow;
strFigTitle = 'performance';
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
%}
%% z-drift mean
%{
%plot
figure
subplot(2,2,1);
vecPlotX = linspace(-5,3,size(matDriftZLicking,3));
errorfill(vecPlotX,squeeze(mean(matDriftZLicking(:,1,:),1)),squeeze(std(matDriftZLicking(:,1,:),[],1)/sqrt(intAnimZ)),[0 0 1],[0.5 0.5 1]);
ylabel('Offset from mean z-plane (micron)')
ylim([0 2]);
xlabel('Time after hit response (s)');
title('Mean +/- st. err. across animals (n=5)');

subplot(2,2,2);
errorfill(vecPlotX,squeeze(mean(matDriftZLicking(:,2,:),1)),squeeze(std(matDriftZLicking(:,2,:),[],1)/sqrt(intAnimZ)),[0 0 1],[0.5 0.5 1]);
ylabel('Inter-frame z-shift (d(micron))')
ylim([0 1]);
xlabel('Time after hit response (s)');
title('Mean +/- st. err. across animals (n=5)');


%make PSTH of z-shifts
subplot(2,2,3);
errorfill(vecPlotX,squeeze(mean(matDriftZLicking(:,3,:),1)),squeeze(std(matDriftZLicking(:,3,:),[],1)/sqrt(intAnimZ)),[0 0 1],[0.5 0.5 1]);
ylabel('Mean pop. dF/F0')
%ylim([0 1]);
xlabel('Time after hit response (s)');
title('Mean +/- st. err. across animals (n=5)');

subplot(2,2,4);
errorfill(vecPlotX,squeeze(mean(matDriftZLicking(:,4,:),1)),squeeze(std(matDriftZLicking(:,4,:),[],1)/sqrt(intAnimZ)),[0 0 1],[0.5 0.5 1]);
ylabel('Heterogeneity')
%ylim([0 1]);
xlabel('Time after hit response (s)');
title('Mean +/- st. err. across animals (n=5)');


%save
drawnow;
strFigTitle = 'zdrifts';
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
%}
%% clean up
cd(strOldDir);