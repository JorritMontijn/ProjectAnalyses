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
strFigDir = 'D:\Data\Results\stimdetection\metaNS';
strDataDir = 'D:\Data\Results\stimdetection';
cellInclude = {'20140207','20140314','20140425','20140507','20140530','20140604','20140711','20140715'};%exclude:20140430
boolSavePlots = true;
intFigCounter = 0;
boolUseNeuropilSubtraction = false;
if boolUseNeuropilSubtraction
	strFigDir = [strFigDir 'NPS'];
end

%% pre-allocate
intAnimal = 0;
vecAnimalSource = [];
vecSubPopulationSource = [];

matSlopes = [];
matRiseTimes = [];

%% create aggregate data structures
sDir = dir(strDataDir);
strOldDir = cd(strDataDir);
intCounterF1 = 0;
intCounterF2 = 0;

for intFile=1:numel(sDir)
	%check if data part 1
	strFile = sDir(intFile).name;
	strRec = getFlankedBy(strFile,'data_aggregateX','_');
	if length(strRec) > 3 && boolUseNeuropilSubtraction && strcmp(strRec(1:3),'NPS')
		intRec = find(strcmp(strRec(4:end),cellInclude),1);
	elseif ~boolUseNeuropilSubtraction
		intRec = find(strcmp(strRec,cellInclude),1);
	else
		intRec = [];
	end
	
	if ~isempty(intRec)
		%load file
		load(strFile);
		
		%update population source
		intNrPops = size(cellSaveNeuralRespSpeed,1);
		intAnimal = intAnimal + 1;
		vecAnimalSource(end+1) = intAnimal;
		vecSubPopulationSource(end+1) = intAnimal;
		if intNrPops > 1 %group populations
			%take mean over matrices
			cellSaveNeuralRespSpeed{1,1} = nanmean(cat(4,cellSaveNeuralRespSpeed{1,1},cellSaveNeuralRespSpeed{2,1}),4);
			cellSaveNeuralRespSpeed{1,2} = nanmean(cat(4,cellSaveNeuralRespSpeed{1,2},cellSaveNeuralRespSpeed{2,2}),4);
		end
		
		
		%assign data
		matSlopes = cat(4,matSlopes,cellSaveNeuralRespSpeed{1,1});
		matRiseTimes = cat(4,matRiseTimes,cellSaveNeuralRespSpeed{1,2});
		intCounterF1 = intCounterF1 + 1;
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

%matSlopes(intC,intType,intAct,intAnimal);
matActSlopes = squeeze(matSlopes(:,:,2,:));
matHetSlopes = squeeze(matSlopes(:,:,1,:));
matActSlopesT = matActSlopes(2:(end-1),:,:);
matHetSlopesT = matHetSlopes(2:(end-1),:,:);

matActSS = squeeze(nanmean(matActSlopesT,1));
matHetSS = squeeze(nanmean(matHetSlopesT,1));

[h,pAMS]=ttest(matActSS(1,:),matActSS(2,:));
[h,pAMF]=ttest(matActSS(1,:),matActSS(3,:));
[h,pASF]=ttest(matActSS(2,:),matActSS(3,:));
[a,b,vecActP_corr] = fdr_bh([pAMS pAMF pASF])

[h,pHMS]=ttest(matHetSS(1,:),matHetSS(2,:));
[h,pHMF]=ttest(matHetSS(1,:),matHetSS(3,:));
[h,pHSF]=ttest(matHetSS(2,:),matHetSS(3,:));
[a,b,vecHetP_corr] = fdr_bh([pHMS pHMF pHSF])

%matRiseTimes(intC,intType,intAct,intAnimal);
matActRiseTimes = squeeze(matRiseTimes(:,:,2,:));
matHetRiseTimes = squeeze(matRiseTimes(:,:,1,:));
matActRiseTimesT = matActRiseTimes(2:(end-1),:,:);
matHetRiseTimesT = matHetRiseTimes(2:(end-1),:,:);

matActRS = squeeze(nanmean(matActRiseTimesT,1));
matHetRS = squeeze(nanmean(matHetRiseTimesT,1));

[h,pAMS]=ttest(matActRS(1,:),matActRS(2,:));
[h,pAMF]=ttest(matActRS(1,:),matActRS(3,:));
[h,pASF]=ttest(matActRS(2,:),matActRS(3,:));

[h,pHMS]=ttest(matHetRS(1,:),matHetRS(2,:));
[h,pHMF]=ttest(matHetRS(1,:),matHetRS(3,:));
[h,pHSF]=ttest(matHetRS(2,:),matHetRS(3,:));

h=figure;
errorbar(0.3,mean(matHetRS(1,:)),std(matHetRS(1,:))/sqrt(intAnimals),'LineStyle','none','Marker','x','Color',[0.7 0 0])
hold on
errorbar(0.5,mean(matHetRS(2,:)),std(matHetRS(2,:))/sqrt(intAnimals),'LineStyle','none','Marker','x','Color',[0.5 0 0.5])
errorbar(0.7,mean(matHetRS(3,:)),std(matHetRS(3,:))/sqrt(intAnimals),'LineStyle','none','Marker','x','Color',[0 0.7 0])
hold off
ylim([0 max(get(gca,'ylim'))]);
xlim([0 1])
ylabel('Rise time heterogeneity (s)')
set(gca,'xtick',[0.3 0.5 0.7],'xticklabel',{'Miss','Slow','Fast'});

if boolSavePlots
	strFigTitle = 'risetime_heterogeneity';
	intFigCounter = intFigCounter + 1;
	drawnow;
	strFig = sprintf('Meta%d_%s_raw',intFigCounter,strFigTitle);
	export_fig([strFig '.tif']);
	export_fig([strFig '.pdf']);
end
%% clean up
cd(strOldDir);