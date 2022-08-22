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
strFigDir = 'D:\Data\Results\stimdetection\meta';
strDataDir = 'D:\Data\Results\stimdetection';
cellInclude = {'20140207','20140314','20140425','20140507','20140530','20140604','20140711','20140715'};%exclude:20140430
boolSavePlots = true;
intFigCounter = 0;

%% pre-allocate
intAnimal = 0;
vecAnimalSource = [];
vecSubPopulationSource = [];
cellPreStim = {};%: dF/F detect/no-detect preceding stimulus frames [5 (c) x 2 (d/nd) cell array] x n (animals) with every cell = vector of values
cellBehavDetect = {};%: behavioral detection [5 (c) x n (animals)] cell array with every cell = vector of detect/no detect values (1/0)
cellBehavRT = {};%: behavioral RT [5 (c) x n (animals)] cell array with every cell = vector of RTs
cellBaseSubtrStimAct = {};%: dF/F detect/no-detect baseline-subtracted activity [5 (c) x 2 (d/nd) cell array] x n (animals) with every cell = vector of values
cellDuringStim = {};%: dF/F detect/no-detect during stimulus frames [5 (c) x 2 (d/nd) cell array] x n (animals) with every cell = vector of values
cellTraceAct = {};%: dF/F trace detect/no-detect pref/non-pref [2 (pref/no-pref) x 2 (detect/no-detect) x 6 (c)] x n (animals) with every cell = [Nreps x Nframes] matrix of trace dF/F values
cellMatrices = {};%: blob matrices in [1 x 4 (detect/no-detect/diff/significance)] cell array with [5 (c) x 31 (bins)] matrix
cellSignalCorrs = {};%: signal correlations [3 (within-HCAR, between, within-non-HCAR) x n (animals/blocks)] with every cell = vector of pairwise correlation values
cellNoiseCorrs = {};%: noise correlations [3 (within-HCAR, between, within-non-HCAR) x n (animals/blocks)] with every cell = vector of pairwise correlation values
cellDCAE = {};%: normalized increase in dF/F for detection trials [5 (c) x n (animals/blocks)] cell array with each cell = vector of neurons with normalized dF/F increase for detection trials
cellSignalCorrsBiDir = {};%: signal correlations [3 (within-HCAR, between, within-non-HCAR) x n (animals/blocks)] with every cell = vector of pairwise correlation values
cellNoiseCorrsBiDir = {};%: noise correlations [3 (within-HCAR, between, within-non-HCAR) x n (animals/blocks)] with every cell = vector of pairwise correlation values
cellNormActDissim = {};%: within-group z-scored activation dissimilarity [2 (high/low HCAR) x 2 (detect/no-detect)] x n (animals/blocks) with every cell = vector of dissimilarity values
cellCorrITD = {};%: inter-trial-distance dependence of assembly consistency [2 (hit/miss) x 3 (ITD/HCAR/non-HCAR) x n (animals/blocks)] with every cell = vector with correlation values
cellRTDependency = {};%: RT dependence of activity dissimilarity [3 (RT/HCAR/non-HCAR) x n (animals/blocks)] with every cell = vector of trials
matMixBootstrappedDecodingOutput1 = [];
matMixBootstrappedDecodingOutput2 = [];
matContHetero = [];
cellHetTime = {};
cellActTime = {};
matRespDecode = [];
matRespDecodeDist = [];
matNeuronConsistency = [];%[corr matrix (X x Y)] x n (animals)

%% create aggregate data structures
sDir = dir(strDataDir);
strOldDir = cd(strDataDir);
intCounterF1 = 0;
intCounterF2 = 0;
intCounterF3 = 0;
intCounterF4 = 0;
for intFile=1:numel(sDir)
	%check if data part 1
	strFile = sDir(intFile).name;
	strRec = getFlankedBy(strFile,'data_aggregate','_');
	intRec = find(strcmp(strRec,cellInclude),1);
	if ~isempty(intRec)
		%load file
		load(strFile);
		
		%assign data
		cellBehavDetect = cat(2,cellBehavDetect,cellSaveBehavDetect);%: behavioral detection [6 (c) x n (animals)] cell array with every cell = vector of detect/no detect values (1/0)
		cellBehavRT = cat(2,cellBehavRT,cellSaveBehavRT);%: behavioral RT [5 (c) x n (animals)] cell array with every cell = vector of RTs
		cellPreStim = cat(3,cellPreStim,cellSavePreStim);%: dF/F detect/no-detect preceding stimulus frames [5 (c) x 2 (d/nd) cell array] x n (animals) with every cell = vector of values
		cellBaseSubtrStimAct = cat(3,cellBaseSubtrStimAct,cellSaveBaseSubtrStimAct);%: dF/F detect/no-detect baseline-subtracted activity [5 (c) x 2 (d/nd) cell array] x n (animals) with every cell = vector of values
		cellDuringStim = cat(3,cellDuringStim,cellSaveDuringStim);%: dF/F detect/no-detect during stimulus frames [5 (c) x 2 (d/nd) cell array] x n (animals) with every cell = vector of values
		cellTraceAct = cat(4,cellTraceAct,cellSaveTraceAct);%: dF/F trace detect/no-detect pref/non-pref [2 (pref/no-pref) x 2 (detect/no-detect) x 6 (c)] x n (animals) with every cell = [Nreps x Nframes] matrix of trace dF/F values
		intCounterF1 = intCounterF1 + 1;
	end
	
	%check if data part 2
	strFile = sDir(intFile).name;
	strRec = getFlankedBy(strFile,'data_aggregateAD','_');
	intRec = find(strcmp(strRec,cellInclude),1);
	if ~isempty(intRec)
		%load file
		load(strFile);
		
		%update population source
		intNrPops = size(cellSaveMatrices,2);
		intAnimal = intAnimal + 1;
		vecAnimalSource(end+1) = intAnimal;
		vecSubPopulationSource(end+1) = intAnimal;
		if intNrPops > 1 %group populations
			%take mean over matrices
			cellSaveMatrices{1,1} = nanmean(cat(3,cellSaveMatrices{1,1},cellSaveMatrices{1,2}),3);
			cellSaveMatrices{2,1} = nanmean(cat(3,cellSaveMatrices{2,1},cellSaveMatrices{2,2}),3);
			cellSaveMatrices{3,1} = nanmean(cat(3,cellSaveMatrices{3,1},cellSaveMatrices{3,2}),3);
			cellSaveMatrices{4,1} = cellSaveMatrices{4,1} | cellSaveMatrices{4,2};
			cellSaveMatrices(:,2) = [];
			
			%mean of consistencies
			cellSaveNeuronConsistency{1}(:,:,2) = cellSaveNeuronConsistency{2};
			cellSaveNeuronConsistency{1} = mean(cellSaveNeuronConsistency{1},3);
			cellSaveNeuronConsistency(2) = [];
			
			%concatenate signal correlations
			cSC{1,1} = [cellSaveSignalCorrs{1,1};cellSaveSignalCorrs{1,2}];
			cSC{2,1} = [cellSaveSignalCorrs{2,1};cellSaveSignalCorrs{2,2}];
			cSC{3,1} = [cellSaveSignalCorrs{3,1};cellSaveSignalCorrs{3,2}];
			cellSaveSignalCorrs = cSC;
			%concatenate noise correlations
			cNC{1,1} = [cellSaveNoiseCorrs{1,1};cellSaveNoiseCorrs{1,2}];
			cNC{2,1} = [cellSaveNoiseCorrs{2,1};cellSaveNoiseCorrs{2,2}];
			cNC{3,1} = [cellSaveNoiseCorrs{3,1};cellSaveNoiseCorrs{3,2}];
			cellSaveNoiseCorrs = cNC;
			
			%concatenate heterogeneity measures
			for intDim1=1:size(cellSaveNormActDissim,1)
				for intDim2=1:size(cellSaveNormActDissim,2)
					for intDim4=1:size(cellSaveNormActDissim,4)
						cellSaveNormActDissim{intDim1,intDim2,1,intDim4} = [cellSaveNormActDissim{intDim1,intDim2,1,intDim4};cellSaveNormActDissim{intDim1,intDim2,2,intDim4}];
					end
				end
			end
			cellSaveNormActDissim(:,:,2,:) = [];
			
			%concatenate ITD correlations
			for intDim1=1:size(cellSaveCorrITD,1)
				for intDim2=1:size(cellSaveCorrITD,2)
					cellSaveCorrITD{intDim1,intDim2,1} = [cellSaveCorrITD{intDim1,intDim2,1} cellSaveCorrITD{intDim1,intDim2,2}];
				end
			end
			cellSaveCorrITD(:,:,2) = [];
			
			%concatenate RT dependency
			cellSaveRTDependency{1,1} = [cellSaveRTDependency{1,1} cellSaveRTDependency{1,2}];
			for intDim1=2:size(cellSaveRTDependency,1)
				cellSaveRTDependency{intDim1,1} = [cellSaveRTDependency{intDim1,1}; cellSaveRTDependency{intDim1,2}];
			end
			cellSaveRTDependency(:,2) = [];
			
			%take mean over heterogen
			metTemp = cellSaveMatContHetero{1};
			matTempHet(:,:,2) = cellSaveMatContHetero{2};
			cellSaveMatContHetero = {mean(matTempHet,3)};
			
			%concatenate heterogen-time
			for intC=1:size(cellSaveHetTime,1)
				cellSaveHetTime{intC,1} = cat(1,cellSaveHetTime{intC,1},cellSaveHetTime{intC,2});
				cellSaveActTime{intC,1} = cat(1,cellSaveActTime{intC,1},cellSaveActTime{intC,2});
			end
			cellSaveHetTime = cellSaveHetTime(:,1);
			cellSaveActTime = cellSaveActTime(:,1);
			
			%concatenate resp type decoding
			cellSaveRespDecode{1,1} = mean(cat(4,cellSaveRespDecode{1,1},cellSaveRespDecode{2,1}),4);
			cellSaveRespDecode{1,2} = mean(cat(4,cellSaveRespDecode{1,2},cellSaveRespDecode{2,2}),4);
			cellSaveRespDecode(2,:) = [];
		end
		
		%assign data
		cellMatrices = cat(2,cellMatrices,cellSaveMatrices);%: blob matrices in [1 x 4 (detect/no-detect/diff/significance)] cell array with [5 (c) x 31 (bins)] matrix
		matNeuronConsistency = cat(3,matNeuronConsistency,cellSaveNeuronConsistency{1});%: [1 x i (neuronal population)] cell array with vector of shuffled correlation values, where first value is actual data
		cellSignalCorrs = cat(2,cellSignalCorrs,cellSaveSignalCorrs);%: signal correlations [3 (within-HCAR, between, within-non-HCAR) x n (animals/blocks)] with every cell = vector of pairwise correlation values
		cellNoiseCorrs = cat(2,cellNoiseCorrs,cellSaveNoiseCorrs);%: noise correlations [3 (within-HCAR, between, within-non-HCAR) x n (animals/blocks)] with every cell = vector of pairwise correlation values
		%cellDCAE = cat(2,cellDCAE,cellSaveDCAE');%: normalized increase in dF/F for detection trials [5 (c) x n (animals/blocks)] cell array with each cell = vector of neurons with normalized dF/F increase for detection trials
		%cellSignalCorrsBiDir = cat(2,cellSignalCorrsBiDir,cellSaveSignalCorrsBiDir);%: signal correlations [3 (within-HCAR, between, within-non-HCAR) x n (animals/blocks)] with every cell = vector of pairwise correlation values
		%cellNoiseCorrsBiDir = cat(2,cellNoiseCorrsBiDir,cellSaveNoiseCorrsBiDir);%: noise correlations [3 (within-HCAR, between, within-non-HCAR) x n (animals/blocks)] with every cell = vector of pairwise correlation values
		cellNormActDissim = cat(3,cellNormActDissim,cellSaveNormActDissim);%: within-group z-scored activation dissimilarity [2 (high/low HCAR) x 2 (detect/no-detect)] x n (animals/blocks) with every cell = vector of dissimilarity values
		cellCorrITD = cat(3,cellCorrITD,cellSaveCorrITD);%: inter-trial-distance dependence of assembly consistency [2 (high/low HCAR) x 2 (hits/misses) x n (animals/blocks)] with every cell = [p x 2] matrix with intertrial distance and correlation values for each trial pair
		cellRTDependency = cat(2,cellRTDependency,cellSaveRTDependency);%: RT dependence of activity dissimilarity [3 (RT/HCAR/non-HCAR) x n (animals/blocks)] with every cell = vector of trials
		matContHetero = cat(3,matContHetero,cellSaveMatContHetero{1});
		cellHetTime = cat(2,cellHetTime,cellSaveHetTime);
		cellActTime = cat(2,cellActTime,cellSaveActTime);
		matRespDecode = cat(4,matRespDecode,cellSaveRespDecode{1,1});
		matRespDecodeDist = cat(3,matRespDecodeDist,cellSaveRespDecode{1,2});
		intCounterF2 = intCounterF2 + 1;
	end
	
	%check for decoding
	strFile = sDir(intFile).name;
	strRec = getFlankedBy(strFile,'data_aggregateSD','_');
	intRec = find(strcmp(strRec,cellInclude),1);
	if ~isempty(intRec)
		%load file
		load(strFile);
		
		%update population source
		intNrPops = numel(cellSaveDecoding)/2;
		if intNrPops > 1 %group populations
			%take mean over matrices
			cellSaveDecoding(1,1).matMixBootstrappedDecodingOutput = nanmean(cat(4,cellSaveDecoding(1,1).matMixBootstrappedDecodingOutput,cellSaveDecoding(1,2).matMixBootstrappedDecodingOutput),4);
			cellSaveDecoding(2,1).matMixBootstrappedDecodingOutput = nanmean(cat(4,cellSaveDecoding(2,1).matMixBootstrappedDecodingOutput,cellSaveDecoding(2,2).matMixBootstrappedDecodingOutput),4);
			cellSaveDecoding(:,2) = [];
		end
		
		%assign data
		matMixBootstrappedDecodingOutput1 = cat(1,matMixBootstrappedDecodingOutput1,mean(cellSaveDecoding(1).matMixBootstrappedDecodingOutput,1));
		matMixBootstrappedDecodingOutput2 = cat(1,matMixBootstrappedDecodingOutput2,mean(cellSaveDecoding(2).matMixBootstrappedDecodingOutput,1));
		intCounterF3 = intCounterF3 + 1;
	end
	
	%check for decoding 2
	strFile = sDir(intFile).name;
	strRec = getFlankedBy(strFile,'data_aggregateSD2_','_');
	intRec = find(strcmp(strRec,cellInclude),1);
	if ~isempty(intRec)
		%load file
		load(strFile);
		
		%update population source
		intNrPops = numel(cellSaveDecoding2);
		if intNrPops > 1 %group populations
			%take mean over matrices
			cellSaveDecoding2(1).matHD_Orientation = cat(1,cellSaveDecoding2(1).matHD_Orientation,cellSaveDecoding2(2).matHD_Orientation);
			cellSaveDecoding2(1).matHD_Contrast = cat(1,cellSaveDecoding2(1).matHD_Contrast,cellSaveDecoding2(2).matHD_Contrast);
			cellSaveDecoding2(1).matHD_ReactionTime = cat(1,cellSaveDecoding2(1).matHD_ReactionTime,cellSaveDecoding2(2).matHD_ReactionTime);
			cellSaveDecoding2(1).matHD_DecodingAccuracy = cat(1,cellSaveDecoding2(1).matHD_DecodingAccuracy,cellSaveDecoding2(2).matHD_DecodingAccuracy);
			cellSaveDecoding2(1).matHD_HeterogeneityPre = cat(1,cellSaveDecoding2(1).matHD_HeterogeneityPre,cellSaveDecoding2(2).matHD_HeterogeneityPre);
			cellSaveDecoding2(1).matHD_HeterogeneityDuring = cat(1,cellSaveDecoding2(1).matHD_HeterogeneityDuring,cellSaveDecoding2(2).matHD_HeterogeneityDuring);
			cellSaveDecoding2(1).matHD_ActivityPre = cat(1,cellSaveDecoding2(1).matHD_ActivityPre,cellSaveDecoding2(2).matHD_ActivityPre);
			cellSaveDecoding2(1).matHD_ActivityDuring = cat(1,cellSaveDecoding2(1).matHD_ActivityDuring,cellSaveDecoding2(2).matHD_ActivityDuring);
			cellSaveDecoding2(2) = [];
		end
		
		%assign data
		cellHD_Orientation{intCounterF4+1} = cellSaveDecoding2(1).matHD_Orientation;
		cellHD_Contrast{intCounterF4+1} = cellSaveDecoding2(1).matHD_Contrast;
		cellHD_ReactionTime{intCounterF4+1} = cellSaveDecoding2(1).matHD_ReactionTime;
		cellHD_DecodingAccuracy{intCounterF4+1} = cellSaveDecoding2(1).matHD_DecodingAccuracy;
		cellHD_HeterogeneityPre{intCounterF4+1} = cellSaveDecoding2(1).matHD_HeterogeneityPre;
		cellHD_HeterogeneityDuring{intCounterF4+1} = cellSaveDecoding2(1).matHD_HeterogeneityDuring;
		cellHD_ActivityPre{intCounterF4+1} = cellSaveDecoding2(1).matHD_ActivityPre;
		cellHD_ActivityDuring{intCounterF4+1} = cellSaveDecoding2(1).matHD_ActivityDuring;
		
		intCounterF4 = intCounterF4 + 1;
	end
end

%check for inconsistencies
if ~(intCounterF1 == intCounterF2 && intCounterF2 == intCounterF3 && intCounterF3 == intCounterF4)
	warning([mfilename ':InconsistentFileNumbers'],'File counters inconsistent; F1=%d; F2=%d; F3=%d; F4=%d',intCounterF1,intCounterF2,intCounterF3,intCounterF4)
end

%% meta analyses part 1
cd(strFigDir);
close all;
intAnimals = size(cellBehavDetect,2);
vecContrasts = [0 0.5 2 8 32 100];

%% behavioral analyses
%% overview graph 0%/100% detection rates with significance calculation
matRTAgg = ones(5,intAnimals);
matDetectAgg = ones(6,intAnimals);
matDetect_C0 = ones(3,intAnimals); %upper, mean, lower
matDetect_C100 = ones(3,intAnimals);
for intAnimal=1:intAnimals
	%get data
	vecC0 = cellBehavDetect{1,intAnimal};
	vecC100 = cellBehavDetect{6,intAnimal};
	for intC=1:6
		matDetectAgg(intC,intAnimal) = mean(cellBehavDetect{intC,intAnimal});
		if intC>1,matRTAgg(intC-1,intAnimal) = mean(cellBehavRT{intC-1,intAnimal});end
	end
	
	%calculate confidence intervals using Clopper-Pearson method
	[dblP_C0,dblCI_C0] = binofit(sum(vecC0),length(vecC0));
	matDetect_C0(1,intAnimal) = dblCI_C0(2); %upper
	matDetect_C0(2,intAnimal) = dblP_C0; %mean
	matDetect_C0(3,intAnimal) = dblCI_C0(1); %lower
	[dblP_C100,dblCI_C100] = binofit(sum(vecC100),length(vecC100));
	matDetect_C100(1,intAnimal) = dblCI_C100(2); %upper
	matDetect_C100(2,intAnimal) = dblP_C100; %mean
	matDetect_C100(3,intAnimal) = dblCI_C100(1); %lower
end
hDetect_C100_C0 = figure;
errorbar(0.4:1:intAnimals,matDetect_C0(2,:),matDetect_C0(3,:)-matDetect_C0(2,:),matDetect_C0(1,:)-matDetect_C0(2,:),'LineStyle','none','Marker','x','Color','b')
hold on
errorbar(0.6:1:intAnimals,matDetect_C100(2,:),matDetect_C100(3,:)-matDetect_C100(2,:),matDetect_C100(1,:)-matDetect_C100(2,:),'LineStyle','none','Marker','x','Color','r')
hold off
legend('0% stimulus contrast','100% stimulus contrast')
ylim([0 1])
ylabel('Behavioral response probability')
xlabel('Animal')
set(gca,'XTick',0.5:1:intAnimals,'XTickLabel',1:intAnimals)
title('Response proportion with 95% confidence intervals using Clopper-Pearson method')
if boolSavePlots
	intFigCounter = intFigCounter + 1;
	drawnow;
	strFig = sprintf('Meta%d_behavPerformance_raw',intFigCounter);
	export_fig([strFig '.tif']);
	export_fig([strFig '.pdf']);
end

%% mean detection rate over contrasts (mean = blue, single animals = grey lines)
hDetect_AllContrasts = figure;
set(gca,'XScale','log','YScale','linear')
vecX = vecContrasts;
vecLabelsX = vecX;
vecX(1) = 0.25;
vecLabelsX(1) = 0;
semilogx(vecX(2:end),matDetectAgg(2:end,:),'-x','Color',[0.5 0.5 0.5])
hold on
semilogx(vecX(1),matDetectAgg(1,:),'-x','Color',[0.5 0.5 0.5])
semilogx(vecX(1),mean(matDetectAgg(1,:),2),'-x','Color',[0 0 1])
semilogx(vecX(2:end),mean(matDetectAgg(2:end,:),2),'-x','Color',[0 0 1])
hold off
set(gca,'XScale','log','YScale','linear')
set(gca,'XTick',vecX,'XTickLabel',vecLabelsX)
ylabel('Fraction correct')
ylim([0 1])
xlabel('Stimulus contrast (%)')
grid on


%regression analysis
%regression
matFitY = matDetectAgg(2:end,:)';
vecFitX = log(vecX(2:end));
matFitX = repmat(vecFitX,[size(matFitY,1) 1]);
stats = regstats(matFitY(:),matFitX(:),'linear',{'tstat'});

%nlinfit
[beta,resid,J,COVB,mse] = nlinfit(matFitX(:),matFitY(:),@linfit,[stats.tstat.beta(2) stats.tstat.beta(1)]);
[ypred,delta] = nlpredci(@linfit,vecFitX,beta,resid,'covar',COVB);

hold on
vecFillX = [vecX(2:end) vecX(end:-1:2)];
vecFillY = [ypred+delta' ypred(end:-1:1)-delta(end:-1:1)'];
fill(vecFillX,vecFillY,[0.5 0.5 0.5]);
semilogx(vecX(2:end),vecFitX * stats.tstat.beta(2) + stats.tstat.beta(1),'k');
hold off
title(sprintf('Behavioral stimulus detection per animal (grey) and mean (blue);regress p=%.9f',stats.tstat.pval(2)))
if boolSavePlots
	intFigCounter = intFigCounter + 1;
	drawnow;
	strFig = sprintf('Meta%d_behavDetect_raw',intFigCounter);
	export_fig([strFig '.tif']);
	export_fig([strFig '.pdf']);
end

%%  mean RT over contrasts (mean = blue, single animals = grey lines)
matRTs = matRTAgg*1000;
hRT_AllContrasts = figure;
set(gca,'XScale','log','YScale','linear')
vecX = vecContrasts(2:end);
vecLabelsX = vecX;
semilogx(vecX,matRTs,'-x','Color',[0.5 0.5 0.5])
hold on
semilogx(vecX,mean(matRTs,2),'-x','Color',[0 0 1])
hold off
set(gca,'XScale','log','YScale','linear')
set(gca,'XTick',vecX,'XTickLabel',vecLabelsX)
ylabel('Mean Reaction Time (ms)')
xlabel('Stimulus contrast (%)')
grid on

%regression analysis
%regression
matFitY = matRTs';
vecFitX = log(vecX);
matFitX = repmat(vecFitX,[size(matFitY,1) 1]);
stats = regstats(matFitY(:),matFitX(:),'linear',{'tstat'});

%nlinfit
[beta,resid,J,COVB,mse] = nlinfit(matFitX(:),matFitY(:),@linfit,[stats.tstat.beta(2) stats.tstat.beta(1)]);
[ypred,delta] = nlpredci(@linfit,vecFitX,beta,resid,'covar',COVB);

hold on
vecFillX = [vecX vecX(end:-1:1)];
vecFillY = [ypred+delta' ypred(end:-1:1)-delta(end:-1:1)'];
fill(vecFillX,vecFillY,[0.5 0.5 0.5]);
semilogx(vecX,vecFitX * stats.tstat.beta(2) + stats.tstat.beta(1),'k');
hold off
title(sprintf('Behavioral RT for hit trials per animal (grey) and mean (blue);regress p=%.9f',stats.tstat.pval(2)))

if boolSavePlots
	intFigCounter = intFigCounter + 1;
	drawnow;
	strFig = sprintf('Meta%d_behavRT_raw',intFigCounter);
	export_fig([strFig '.tif']);
	export_fig([strFig '.pdf']);
end

%% dF/F over contrasts for detect/no detect
%for pre-stim, during-stim, pre-subtracted-stim
%mean dF/F +/- st err over animals for detect/nodetect trials
for intType =1:3
	%get type
	if intType==1
		hActPreStim = figure;
		strFigTitle = 'ActPreStim';
		strLabelY = 'dF/F preceding stimulus';
		cellData = cellPreStim;
	elseif intType==2
		hActDurStim = figure;
		strFigTitle = 'ActDurStim';
		strLabelY = 'dF/F during stimulus';
		cellData = cellDuringStim;
	elseif intType==3
		hActPBSStim = figure;
		strFigTitle = 'ActBaSuStim';
		strLabelY = 'Baseline-subtracted dF/F during stimulus';
		cellData = cellBaseSubtrStimAct;
	end
	
	%get data
	matActHit = ones(5,intAnimals);
	matActMiss = ones(5,intAnimals);
	vecP = ones(5,1);
	
	for intC=1:size(cellPreStim,1)
		for intAnimal=1:intAnimals
			matActHit(intC,intAnimal) = mean(cellData{intC,1,intAnimal});
			matActMiss(intC,intAnimal) = mean(cellData{intC,2,intAnimal});
		end
		
		%ttest
		[h,vecP(intC),ci] = ttest(matActHit(intC,:),matActMiss(intC,:));
	end
	
	%overall t-test
	matHits = matActHit(2:4,:);
	matMisses = matActMiss(2:4,:);
	[h,dblP_All,ci] = ttest(matHits(:),matMisses(:));
	
	%pre-compute variables
	vecWindow = [2 6];
	vecWindowSelect = vecWindow(1):vecWindow(end);
	intWL = length(vecWindowSelect);
	vecWindowInv = intWL:-1:1;
	vecWindowPlotInv = vecWindow(end):-1:vecWindow(1);
	vecX = [vecContrasts(vecWindowSelect) vecContrasts(vecWindowPlotInv)];
	vecLineX = vecContrasts(vecWindowSelect);
	
	%get data
	for intDetect=[0 1]
		if intDetect == 1
			vecMeanTrace = mean(matActHit,2)';
			vecSE = std(matActHit,[],2)'/sqrt(size(matActHit,2));
			vecColorFill = [0.7 1.0 0.7];
			vecColorLine = [0 1 0];
		else
			vecColorLine = [1 0 0];
			vecColorFill = [1 0.7 0.7];
			vecMeanTrace = mean(matActMiss,2)';
			vecSE = std(matActMiss,[],2)'/sqrt(size(matActMiss,2));
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
	strTitle = 'T-test p: ';
	for intC=1:length(vecP)
		strTitle = [strTitle sprintf('C%.1f: ',vecLineX(intC)) sprintf('%.3f; ',vecP(intC))];
	end
	strTitle = [strTitle '; overall, p=' sprintf('%.3f',dblP_All)];
	
	title(strTitle)
	grid on
	xlabel('Contrast')
	ylabel(strLabelY)
	xlim(vecContrasts(vecWindow))
	ylim([-0.02 0.08])
	legend({'StErr','Miss','StErr','Hit'},'Location','Best')
	
	if boolSavePlots
		intFigCounter = intFigCounter + 1;
		drawnow;
		strFig = sprintf('Meta%d_%s_raw',intFigCounter,strFigTitle);
		export_fig([strFig '.tif']);
		export_fig([strFig '.pdf']);
	end
end

%% trace activity plots
%2 (detect/nodetect) x 2 (pref neurons/non-pref neurons) plot for each contrast
%mean (blue line) and single animals (grey lines)
hTraceAct = zeros(1,6);
cellStrDetect = {'Hit','Miss'};
cellStrPopType = {'Pref','Non-pref'};
for intC=1:6
	hTraceAct(intC) = figure;
	intPlotCounter = 0;
	for intDetect=1:2
		strDetect = cellStrDetect{intDetect};
		for intPrefPop=1:2
			strPopType = cellStrPopType{intPrefPop};
			%get data
			intFrames = size(cellTraceAct{intPrefPop,intDetect,intC,1},2);
			matAct = zeros(intAnimals,intFrames);
			for intAnimal = 1:intAnimals
				matAct(intAnimal,:) = nanmean(cellTraceAct{intPrefPop,intDetect,intC,intAnimal},1);
			end
			
			%plot
			intPlotCounter = intPlotCounter + 1;
			subplot(2,2,intPlotCounter);
			
			%pre-compute variables
			dblSamplingFreq = 25.36;
			vecWindow = round([-3*dblSamplingFreq 5*dblSamplingFreq]);
			vecWindowSelect = vecWindow(1):vecWindow(end);
			intWL = length(vecWindowSelect);
			vecWindowInv = intWL:-1:1;
			vecWindowPlotInv = length(vecWindowSelect):-1:1;
			vecX = [vecWindowSelect vecWindowSelect(vecWindowPlotInv)]/dblSamplingFreq;
			vecLineX = vecWindowSelect/dblSamplingFreq;
			
			
			vecColorLine = [0 0 1];
			vecColorFill = [0.7 0.7 1];
			vecMeanTrace = mean(matAct,1);
			vecSE = std(matAct,[],1)/sqrt(size(matAct,1));
			
			
			vecMinTrace = vecMeanTrace-vecSE;
			vecMaxTrace = vecMeanTrace+vecSE;
			vecY = [vecMinTrace vecMaxTrace(vecWindowInv)];
			
			%plot
			hold on
			
			fill(vecX,vecY,vecColorFill,'EdgeColor','none');
			plot(vecLineX,vecMeanTrace,'-','LineWidth',2,'Color',vecColorLine);
			hold off
			xlim([-3 5])
			ylim([-0.02 0.08])
			title(sprintf('Contrast %.1f%%; %s %s',vecContrasts(intC),strDetect,strPopType))
			grid on
			xlabel('Time after stimulus onset (s)')
			ylabel('Mean pop dF/F over animals')
		end
	end
	if boolSavePlots
		intFigCounter = intFigCounter + 1;
		drawnow;
		strFig = sprintf('Meta%d_traceAct_C%.1f_raw',intFigCounter,vecContrasts(intC));
		export_fig([strFig '.tif']);
		export_fig([strFig '.pdf']);
	end
end


%% meta analyses part 2
%% blob matrices
%2 x 2 heat map plot; detect/no-detect/diff/significance
%mean over animals
%cellMatrices = cat(2,cellMatrices,cellSaveMatrices');
%: blob matrices in [1 x 4 (detect/no-detect/diff/significance)] cell array with [5 (c) x 31 (bins)] matrix
intSwitchZ = 1;
if intSwitchZ == 1
	vecBins = -1:0.1:2;
	vecContrasts = [0 0.5 2 8 32 100];
	vecTickY = [1 11 21 31];
	vecTickLabelY = vecBins(vecTickY);
	vecTickLabelX = vecContrasts;
	strLabelY = sprintf('Z-scored dF/F0 distribution');
else
	vecBins = -0.05:0.01:0.25;
	vecContrasts = [0 0.5 2 8 32 100];
	vecTickY = [6 16 26];%[6 16 26 36];%[1 6 11 16 20];
	vecTickLabelY = vecBins(vecTickY);
	vecTickLabelX = vecContrasts;
	strLabelY = sprintf('Normalized dF/F0 distribution');
end

hFigHeatMaps = figure;
for intType = 1:4
	if intType == 1
		strTitle = sprintf('Hit; z=normalized count; mean over animals (n=%d)',size(cellMatrices,2));
	elseif intType == 2
		strTitle = sprintf('Miss; z=normalized count; mean over animals (n=%d)',size(cellMatrices,2));
	elseif intType == 3
		strTitle = sprintf('Difference; z=normalized count; mean over animals (n=%d)',size(cellMatrices,2));
	elseif intType == 4
		strTitle = sprintf('Significance; z=normalized count; mean over animals (n=%d)',size(cellMatrices,2));
	end
	
	subplot(2,2,intType)
	matAgg = zeros(size(cellMatrices{intType,1},1),size(cellMatrices{intType,1},2),size(cellMatrices,2));
	for intAnimal=1:size(cellMatrices,2)
		matAgg(:,:,intAnimal) = double(cellMatrices{intType,intAnimal});
	end
	cellMatMean{intType} = nanmean(matAgg,3);
	if intType == 4
		cellMatMean{intType} = cellMatMean{intType}*size(cellMatrices,2);
		matDiffT = matDiff;
		matDiffT(isnan(matDiffT)) = nanmean(matDiff(:));
		matP = reshape(zscore(matDiffT(:)),size(matDiff));
		cellMatMean{intType} = matP>1.5;
	end
	if intType < 3,imagesc(cellMatMean{intType},[0 0.075]);colormap(hot);colorbar;drawnow;freezeColors;cbfreeze;
	elseif intType == 3
		matDiff = conv2(cellMatMean{1} ./ (cellMatMean{1} + cellMatMean{2}),(1/3)*ones(3,1), 'same');
		imagesc(matDiff,[0.4-eps 0.6+eps]);
		colormap(redblue);colorbar;drawnow;freezeColors;cbfreeze;
	else imagesc(cellMatMean{intType});colormap(hot);colorbar;drawnow;freezeColors;cbfreeze;
	end
	
	%plot
	axis xy
	title(strTitle)
	set(gca,'YTick',vecTickY,'YTickLabel',vecTickLabelY)
	set(gca,'XTick',1:length(vecTickLabelX),'XTickLabel',vecTickLabelX)
	ylabel(strLabelY)
	xlabel('Contrast (%)')
end
if boolSavePlots
	intFigCounter = intFigCounter + 1;
	drawnow;
	jFig = get(handle(gcf), 'JavaFrame');
	jFig.setMaximized(true);
	figure(gcf);
	drawnow;
	strFig = sprintf('Meta%d_heatmaps_raw',intFigCounter);
	export_fig([strFig '.tif']);
	export_fig([strFig '.pdf']);
end

%% hit-correlated neuron consistency
%matNeuronConsistency

hConsistencyFigure = figure;

matMeanNeuronConsistency = mean(matNeuronConsistency,3);
imagesc(matMeanNeuronConsistency,[-0.2 0.2]);
colormap(redblue);
colorbar;
set(gca,'XTick',1:length(vecTickLabelX),'XTickLabel',vecTickLabelX)
set(gca,'YTick',1:length(vecTickLabelX),'YTickLabel',vecTickLabelX)
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

%% signal correlations
%1 x 3 plot; within-HCAR, between, within-non-HCAR neurons
%mean + standard error over animals
%cellSignalCorrs = cat(2,cellSignalCorrs,cellSaveSignalCorrs);
%: signal correlations [3 (within-HCAR, between, within-non-HCAR) x n (populations)] with every cell = vector of pairwise correlation values
hSignalCorrs = figure;
intPopulations = size(cellSignalCorrs,2);
matGroup = nan(3,intPopulations);
matDataSC = nan(3,intPopulations);
for intPopulation = 1:intPopulations
	matGroup(:,intPopulation) = 1:3;
	matDataSC(1,intPopulation) = mean(cellSignalCorrs{1,intPopulation});
	matDataSC(2,intPopulation) = mean(cellSignalCorrs{2,intPopulation});
	matDataSC(3,intPopulation) = mean(cellSignalCorrs{3,intPopulation});
end
vecMean = mean(matDataSC,2);
vecErr = std(matDataSC,[],2)/sqrt(intPopulations);

[h,dblP_AN] = ttest2(matDataSC(1,:),matDataSC(3,:)); %assembly/non-assembly
[h,dblP_AB] = ttest2(matDataSC(1,:),matDataSC(2,:));%assembly/between
[h,dblP_BN] = ttest2(matDataSC(2,:),matDataSC(3,:));%between/non-assembly

errorbar(1:3+0.1,vecMean,vecErr,'ob','LineStyle','none');
hold on
scatter(matGroup(:)-0.1,matDataSC(:),'rx')
hold off
set(gca,'XTick',[1 2 3],'XTickLabel',{'HCAR neurons','Between','Non-HCAR neurons'})
ylabel('Mean Signal Correlation')
title(sprintf('Mean +/- st err of noise correlation; C-N, p=%.3f; C-B, p=%.3f; B-N, p=%.3f',dblP_AN,dblP_AB,dblP_BN))
if boolSavePlots
	intFigCounter = intFigCounter + 1;
	drawnow;
	strFig = sprintf('Meta%d_signalcorrelations_raw',intFigCounter);
	export_fig([strFig '.tif']);
	export_fig([strFig '.pdf']);
end
vecIncreaseHCAR_SC = matDataSC(1,:)-matDataSC(3,:);

%% noise correlations
%1 x 3 plot; within-HCAR, between, within-non-HCAR neurons
%mean + standard error over animals
%cellNoiseCorrs = cat(2,cellNoiseCorrs,cellSaveNoiseCorrs);
%: noise correlations [3 (within-HCAR, between, within-non-HCAR) x n (populations)] with every cell = vector of pairwise correlation values

intPopulations = size(cellNoiseCorrs,2);
if boolSavePlots,intFigCounter = intFigCounter + 1;end
%for boolNormalize = [false true]
boolNormalize = false;
hNoiseCorrs = figure;
matGroup = nan(3,intPopulations);
matDataNC = nan(3,intPopulations);
for intPopulation = 1:intPopulations
	matGroup(:,intPopulation) = 1:3;
	matDataNC(1,intPopulation) = mean(cellNoiseCorrs{1,intPopulation});
	matDataNC(2,intPopulation) = mean(cellNoiseCorrs{2,intPopulation});
	matDataNC(3,intPopulation) = mean(cellNoiseCorrs{3,intPopulation});
end
if boolNormalize,matDataNC = matDataNC./repmat(matDataNC(2,:), [3 1]);end
vecMean = mean(matDataNC,2);
vecErr = std(matDataNC,[],2)/sqrt(intPopulations);

[h,dblP_AN] = ttest2(matDataNC(1,:),matDataNC(3,:)); %assembly/non-assembly
[h,dblP_AB] = ttest2(matDataNC(1,:),matDataNC(2,:));%assembly/between
[h,dblP_BN] = ttest2(matDataNC(2,:),matDataNC(3,:));%between/non-assembly

errorbar(1:3+0.1,vecMean,vecErr,'ob','LineStyle','none');
hold on
scatter(matGroup(:)-0.1,matDataNC(:),'rx')
hold off
set(gca,'XTick',[1 2 3],'XTickLabel',{'HCAR neurons','Between','Non-HCAR neurons'})
if ~boolNormalize,ylabel('Mean Noise Correlation'),else ylabel('Increase in noise correlation relative to between');end
title(sprintf('Mean +/- st err of noise correlation; C-N, p=%.3f; C-B, p=%.3f; B-N, p=%.3f',dblP_AN,dblP_AB,dblP_BN))
if boolSavePlots
	drawnow;
	strFig = sprintf('Meta%d_noisecorrelations%d_raw',intFigCounter,boolNormalize);
	export_fig([strFig '.tif']);
	export_fig([strFig '.pdf']);
end
vecIncreaseHCAR_NC = matDataNC(1,:)-matDataNC(3,:);
%end

%% relative change in correlations
hDiffCorrs = figure;
scatter(vecIncreaseHCAR_NC,vecIncreaseHCAR_SC,'rx')
hold on
dblLimX = max(abs(get(gca,'xlim')));
plot([-dblLimX dblLimX],[0 0],'k')
set(gca,'xlim',[-dblLimX dblLimX]);
dblLimY = max(abs(get(gca,'ylim')));
plot([0 0],[-dblLimY dblLimY],'k')
set(gca,'ylim',[-dblLimY dblLimY]);

%ttest
[hSC,pSC,ciSC] = ttest(vecIncreaseHCAR_SC);
[hNC,pNC,ciNC] = ttest(vecIncreaseHCAR_NC);

plot([mean(vecIncreaseHCAR_NC) mean(vecIncreaseHCAR_NC)],ciSC,'r')
plot(ciNC,[mean(vecIncreaseHCAR_SC) mean(vecIncreaseHCAR_SC)],'r')
hold off

xlabel('Increase in Noise Correlation for HCAR neurons vs Non-HCAR neurons')
ylabel('Increase in Signal Correlation for HCAR neurons vs Non-HCAR neurons')
title(sprintf('Difference in corrs; SC, p=%.3f; NC, p=%.3f; cross is 95%% CI',pSC,pNC))
if boolSavePlots
	intFigCounter = intFigCounter + 1;
	drawnow;
	strFig = sprintf('Meta%d_diffcorrelations%d_raw',intFigCounter,boolNormalize);
	export_fig([strFig '.tif']);
	export_fig([strFig '.pdf']);
end

%% neuron-based DCAE distribution
%sum over animals, then normalize peak to 1, over contrasts
%cellDCAE = cat(2,cellDCAE,cellSaveDCAE');
%: normalized increase in dF/F for detection trials [5 (c) x n (populations)] cell array with each cell = vector of neurons with normalized dF/F increase for detection trials


%% within-group z-scored activation dissimilarity
%% reaction-time dependency

vecX = [0 3];
hRTDependenceFig = figure;

for intActType=1:3
	if intActType == 1
		%dissimilarity
		intHCAR = 2;
		intNonHCAR = 3;
		strLabelY = 'Mean population activity heterogeneity';
	elseif intActType == 2
		%z-scored act
		intHCAR = 4;
		intNonHCAR = 5;
		strLabelY = 'Mean population z-scored activity';
	elseif intActType == 3
		%dF/F act
		intHCAR = 6;
		intNonHCAR = 7;
		strLabelY = 'Mean population dF/F activity';
	end
	
	matRegHCAR = zeros(intPopulations,length(vecX));
	matRegNonHCAR = zeros(intPopulations,length(vecX));
	vecSlopesHCAR = zeros(intPopulations,1);
	vecSlopesNonHCAR = zeros(intPopulations,1);
	
	subplot(2,3,intActType);
	hold on;
	vecAggRTs = [];
	vecAggActHCAR = [];
	vecAggActNonHCAR = [];
	for intPopulation=1:intPopulations
		vecAggRTs = [vecAggRTs cellRTDependency{1,intPopulation}];
		vecAggActHCAR = [vecAggActHCAR cellRTDependency{intHCAR,intPopulation}'];
		vecAggActNonHCAR = [vecAggActNonHCAR cellRTDependency{intNonHCAR,intPopulation}'];
		
		
		%do regressions per population
		sStatsC=regstats(cellRTDependency{intHCAR,intPopulation}',cellRTDependency{1,intPopulation},'linear');
		vecY = polyval(sStatsC.beta([2 1]),vecX);
		matRegHCAR(intPopulation,:) = vecY;
		vecSlopesHCAR(intPopulation) = sStatsC.beta(2);
		plot(vecX,vecY,'Color',[1.0 0.5 0.5]);
		
		sStatsU=regstats(cellRTDependency{intNonHCAR,intPopulation}',cellRTDependency{1,intPopulation},'linear');
		vecY = polyval(sStatsU.beta([2 1]),vecX);
		matRegNonHCAR(intPopulation,:) = vecY;
		vecSlopesNonHCAR(intPopulation) = sStatsU.beta(2);
		plot(vecX,vecY,'Color',[0.5 0.5 1.0]);
	end
	
	%plot means
	vecMeanY = mean(matRegHCAR,1);
	plot(vecX,vecMeanY,'Color',[1 0 0],'LineWidth',2)
	
	vecMeanY = mean(matRegNonHCAR,1);
	plot(vecX,vecMeanY,'Color',[0 0 1],'LineWidth',2)
	
	%ttest
	[h,dblPSlopeHCAR] = ttest(vecSlopesHCAR);
	[h,dblPSlopeNonHCAR] = ttest(vecSlopesNonHCAR);
	
	title(sprintf('Mean of linear regressions over animals; HCAR (red): slope-p=%.3f; Non-HCAR (blue): slope-p=%.3f',...
		dblPSlopeHCAR,dblPSlopeNonHCAR))
	
	xlabel('Reaction Time (s)')
	ylabel(strLabelY)
	
	%plot aggregate
	subplot(2,3,intActType+3);
	scatter(vecAggRTs,vecAggActHCAR,'r')
	hold on
	scatter(vecAggRTs,vecAggActNonHCAR,'b')
	drawnow
	
	%perform regressions
	sStatsC=regstats(vecAggActHCAR,vecAggRTs,'linear');
	vecX = [0 3];
	vecY = polyval(sStatsC.beta([2 1]),vecX);
	plot(vecX,vecY,'r')
	
	sStatsU=regstats(vecAggActNonHCAR,vecAggRTs,'linear');
	vecY = polyval(sStatsU.beta([2 1]),vecX);
	plot(vecX,vecY,'b')
	hold off
	
	title(sprintf('Aggregate data set; Lin reg HCAR (red): slope=%.3f; p=%.3f; intercept=%.3f; p=%.3f; R^2=%.3f\nLin reg non-HCAR (blue): slope=%.3f; p=%.3f; intercept=%.3f; p=%.3f; R^2=%.3f',...
		sStatsC.beta(2),sStatsC.tstat.pval(2),sStatsC.beta(1),sStatsC.tstat.pval(1),sStatsC.rsquare,sStatsU.beta(2),sStatsU.tstat.pval(2),sStatsU.beta(1),sStatsU.tstat.pval(1),sStatsU.rsquare))
	
	xlabel('Reaction Time (s)')
	ylabel(strLabelY)
end
if boolSavePlots
	intFigCounter = intFigCounter + 1;
	drawnow;
	jFig = get(handle(gcf), 'JavaFrame');
	jFig.setMaximized(true);
	figure(gcf);
	drawnow;
	strFig = sprintf('Meta%d_activationdissimilarity_RTDependency_raw',intFigCounter);
	export_fig([strFig '.tif']);
	export_fig([strFig '.pdf']);
end

%% population activity heterogeneity (similarity distance) with removal most/least active active
hZscoredActivityDissimilarity = figure;
intPopulations = size(cellNormActDissim,3);
for intRemType=1:5
	if intRemType == 1
		%no rem
		strRemType = 'No';
		intPlotNr = 1;
	elseif intRemType == 2
		%rem hi
		strRemType = 'High Z';
		intPlotNr = 2;
	elseif intRemType == 3
		%rem lo
		strRemType = 'Low Z';
		intPlotNr = 3;
	elseif intRemType == 4
		%rem hi
		strRemType = 'High A';
		intPlotNr = 5;
	elseif intRemType == 5
		%rem lo
		strRemType = 'Low A';
		intPlotNr = 6;
	end
	matGroup = nan(4,intPopulations);
	matData = nan(4,intPopulations);
	for intPopulation = 1:intPopulations
		matGroup(:,intPopulation) = 1:4;
		matData(4,intPopulation) = mean(cellNormActDissim{1,1,intPopulation,intRemType}); %high/hit
		matData(3,intPopulation) = mean(cellNormActDissim{1,2,intPopulation,intRemType});%high/miss
		matData(2,intPopulation) = mean(cellNormActDissim{2,1,intPopulation,intRemType});%low/hit
		matData(1,intPopulation) = mean(cellNormActDissim{2,2,intPopulation,intRemType});%low/miss
	end
	
	% perform ttests
	matP = nan(4,4);
	[h,matP(1,2)] = ttest(matData(1,:),matData(2,:)); %assembly/non-assembly
	[h,matP(1,3)] = ttest(matData(1,:),matData(3,:)); %assembly/non-assembly
	[h,matP(1,4)] = ttest(matData(1,:),matData(4,:)); %assembly/non-assembly
	[h,matP(2,3)] = ttest(matData(2,:),matData(3,:)); %assembly/non-assembly
	[h,matP(2,4)] = ttest(matData(2,:),matData(4,:)); %assembly/non-assembly
	[h,matP(3,4)] = ttest(matData(3,:),matData(4,:)); %assembly/non-assembly
	
	%normalize
	matData = matData./repmat(mean(matData,1),[4 1]);
	vecMean = nanmean(matData,2);
	vecErr = nanstd(matData,[],2)/sqrt(intPopulations);
	
	%plot
	subplot(2,3,intPlotNr)
	errorbar((1:4)+0.1,vecMean,vecErr,'ob','LineStyle','none');
	hold on
	scatter(matGroup(:)-0.1,matData(:),'rx')
	hold off
	set(gca,'XTick',1:4,'XTickLabel',{'Miss Non-HCAR','Hit Non-HCAR','Miss HCAR','Hit HCAR'})
	ylabel('Normalized mean within-group activation dissimilarity')
	title(sprintf('%s rem; T-test p-values, 1-2=%.3f; 1-3=%.3f; 1-4=%.3f; 2-3=%.3f; 2-4=%.3f; 3-4=%.3f;',strRemType,matP(1,2),matP(1,3),matP(1,4),matP(2,3),matP(2,4),matP(3,4)))
end

if boolSavePlots
	intFigCounter = intFigCounter + 1;
	drawnow;
	jFig = get(handle(gcf), 'JavaFrame');
	jFig.setMaximized(true);
	figure(gcf);
	drawnow;
	strFig = sprintf('Meta%d_activationdissimilarity_raw',intFigCounter);
	export_fig([strFig '.tif']);
	export_fig([strFig '.pdf']);
end

%% inter-trial-distance dependence of assembly consistency
%2 (high/low HCAR) x 2 (hits/misses) plot
%	cellSaveCorrITD{1,1,intPopulation} = vecDistHit;
%	cellSaveCorrITD{1,2,intPopulation} = vecCorrHitC;
%	cellSaveCorrITD{1,3,intPopulation} = vecCorrHitU;
%	cellSaveCorrITD{2,1,intPopulation} = vecDistMiss;
%	cellSaveCorrITD{2,2,intPopulation} = vecCorrMissC;
%	cellSaveCorrITD{2,3,intPopulation} = vecCorrMissU;

%get data
vecDistHit = cell2mat(shiftdim(squeeze(cellCorrITD(1,1,:)),1));
vecCorrHitC = cell2mat(shiftdim(squeeze(cellCorrITD(1,2,:)),1));
vecCorrHitU = cell2mat(shiftdim(squeeze(cellCorrITD(1,3,:)),1));
vecDistMiss = cell2mat(shiftdim(squeeze(cellCorrITD(2,1,:)),1));
vecCorrMissC = cell2mat(shiftdim(squeeze(cellCorrITD(2,2,:)),1));
vecCorrMissU = cell2mat(shiftdim(squeeze(cellCorrITD(2,3,:)),1));
vecDistBetween = cell2mat(shiftdim(squeeze(cellCorrITD(3,1,:)),1));
vecCorrBetweenC = cell2mat(shiftdim(squeeze(cellCorrITD(3,2,:)),1));
vecCorrBetweenU = cell2mat(shiftdim(squeeze(cellCorrITD(3,3,:)),1));

%plot
hCorrITD = figure;

subplot(2,2,1)
%high HCAR
intStep = 25;
intMax = ceil(max(max(vecDistHit),max(vecDistMiss))/intStep)*intStep;
vecBinX = 0:intStep:intMax;
vecPlotX = (intStep/2):intStep:(intMax-intStep/2);
[nVec,meanVec,stdVec,cellVals,cellIDs] = makeBins(vecDistHit,vecCorrHitC,vecBinX);
errorbar(vecPlotX,meanVec,stdVec./nVec, 'Color', 'r')
sStatsH=regstats(vecCorrHitC,vecDistHit,'linear',{'beta','rsquare','tstat'});
vecX = get(gca,'XLim');
vecY = polyval(sStatsH.beta([2 1]),vecX);
hold on
plot(vecX,vecY,'r')

%low HCAR
vecBinX = 0:intStep:intMax;
vecPlotX = (intStep/2):intStep:(intMax-intStep/2);
[nVec,meanVec,stdVec,cellVals,cellIDs] = makeBins(vecDistHit,vecCorrHitU,vecBinX);
errorbar(vecPlotX,meanVec,stdVec./nVec, 'Color', 'b')
sStatsL=regstats(vecCorrHitU,vecDistHit,'linear',{'beta','rsquare','tstat'});
vecX = get(gca,'XLim');
vecY = polyval(sStatsL.beta([2 1]),vecX);
plot(vecX,vecY,'b')
hold off
title(sprintf('Hits; Lin reg high HCAR (red): slope=%.3f; p=%.3f; intercept=%.3f; p=%.3f; R^2=%.3f\nLin reg low HCAR (blue): slope=%.3f; p=%.3f; intercept=%.3f; p=%.3f; R^2=%.3f',...
	sStatsH.beta(2),sStatsH.tstat.pval(2),sStatsH.beta(1),sStatsH.tstat.pval(1),sStatsH.rsquare,sStatsL.beta(2),sStatsL.tstat.pval(2),sStatsL.beta(1),sStatsL.tstat.pval(1),sStatsL.rsquare))
ylabel(sprintf('Assembly consistency over time\n(Inter-trial correlation of mean population activity)'))
xlabel('Inter-trial distance')

subplot(2,2,3)
%high HCAR
vecBinX = 0:intStep:intMax;
vecPlotX = (intStep/2):intStep:(intMax-intStep/2);
[nVec,meanVec,stdVec,cellVals,cellIDs] = makeBins(vecDistMiss,vecCorrMissC,vecBinX);
errorbar(vecPlotX,meanVec,stdVec./nVec, 'Color', 'm')
sStatsH=regstats(vecCorrMissC,vecDistMiss,'linear',{'beta','rsquare','tstat'});
vecX = get(gca,'XLim');
vecY = polyval(sStatsH.beta([2 1]),vecX);
hold on
plot(vecX,vecY,'m')

%low HCAR
vecBinX = 0:intStep:intMax;
vecPlotX = (intStep/2):intStep:(intMax-intStep/2);
[nVec,meanVec,stdVec,cellVals,cellIDs] = makeBins(vecDistMiss,vecCorrMissU,vecBinX);
errorbar(vecPlotX,meanVec,stdVec./nVec, 'Color', 'c')
sStatsL=regstats(vecCorrMissU,vecDistMiss,'linear',{'beta','rsquare','tstat'});
vecX = get(gca,'XLim');
vecY = polyval(sStatsL.beta([2 1]),vecX);
plot(vecX,vecY,'c')
hold off
title(sprintf('Misses; Lin reg high HCAR (magenta): slope=%.3f; p=%.3f; intercept=%.3f; p=%.3f; R^2=%.3f\nLin reg low HCAR (cyan): slope=%.3f; p=%.3f; intercept=%.3f; p=%.3f; R^2=%.3f',...
	sStatsH.beta(2),sStatsH.tstat.pval(2),sStatsH.beta(1),sStatsH.tstat.pval(1),sStatsH.rsquare,sStatsL.beta(2),sStatsL.tstat.pval(2),sStatsL.beta(1),sStatsL.tstat.pval(1),sStatsL.rsquare))
ylabel(sprintf('Assembly consistency over time\n(Inter-trial correlation of mean population activity)'))
xlabel('Inter-trial distance')


if boolSavePlots
	intFigCounter = intFigCounter + 1;
	drawnow;
	strFig = sprintf('Meta%d_activation_ITD_correlation_raw',intFigCounter);
	export_fig([strFig '.tif']);
	export_fig([strFig '.pdf']);
end


%% overall correlation of population activation
%1 x 4 plot; 1 data point = mean correlation per population for
%high/low/hits/misses combination; plot is mean +/- st err of data points
hAssemblyConsistency = figure;
intPopulations = size(cellCorrITD,3);
matGroup = nan(6,intPopulations);
matData = nan(6,intPopulations);
for intPopulation = 1:intPopulations
	matGroup(:,intPopulation) = 1:6;
	matData(1,intPopulation) = nanmean(cellCorrITD{2,3,intPopulation}); %miss HUN
	matData(2,intPopulation) = nanmean(cellCorrITD{1,3,intPopulation}); %hit HUN
	matData(3,intPopulation) = nanmean(cellCorrITD{2,2,intPopulation}); %miss HCN
	matData(4,intPopulation) = nanmean(cellCorrITD{1,2,intPopulation}); %hit HCN
	matData(5,intPopulation) = nanmean(cellCorrITD{3,2,intPopulation}); %between HCN
	matData(6,intPopulation) = nanmean(cellCorrITD{3,3,intPopulation}); %between HUN
end


dblErr1 = std(matData(1,:)-matData(6,:))/sqrt(length(matData(1,:)));
dblErr2 = std(matData(2,:)-matData(6,:))/sqrt(length(matData(2,:)));
dblErr3 = std(matData(3,:)-matData(5,:))/sqrt(length(matData(3,:)));
dblErr4 = std(matData(4,:)-matData(5,:))/sqrt(length(matData(4,:)));
errorbar(1:4,[mean(matData(1,:)-matData(6,:)) mean(matData(2,:)-matData(6,:)) mean(matData(3,:)-matData(5,:)) mean(matData(4,:)-matData(5,:))],[dblErr3 dblErr4 dblErr1 dblErr2],'Linestyle','none','Marker','x');
set(gca,'XTick',1:4,'XTickLabel',{'Miss HUN','Hit HUN','Miss HCN','Hit HCN'})
xlim([0.5 4.5])
ylabel('Mean inter-trial correlation of population activity')

% perform ttests
intGroups=4;
matP = nan(intGroups,intGroups);
strTitle = 't-tests';
for intGroup1=1:(intGroups-1)
	for intGroup2=(intGroup1+1):intGroups
		[h,matP(intGroup1,intGroup2)] = ttest(matData(intGroup1,:),matData(intGroup2,:)); %assembly/non-assembly
		strTitle = [strTitle sprintf(';%d-%d=%.3f',intGroup1,intGroup2,matP(intGroup1,intGroup2))];
	end
end

%hit-miss
[h,dblRespP] = ttest2([matData(2,:)-matData(6,:) matData(4,:)-matData(5,:)],[matData(1,:)-matData(6,:) matData(3,:)-matData(5,:)]);
strTitle = [strTitle sprintf(';H-M=%.3f',dblRespP)];
[h,dblGroupP] = ttest2([matData(1,:)-matData(6,:) matData(2,:)-matData(6,:)],[matData(3,:)-matData(5,:) matData(4,:)-matData(5,:)]);
strTitle = [strTitle sprintf(';C-U=%.3f',dblGroupP)];
title(strTitle);

if boolSavePlots
	intFigCounter = intFigCounter + 1;
	drawnow;
	strFig = sprintf('Meta%d_activation_correlation_raw',intFigCounter);
	export_fig([strFig '.tif']);
	export_fig([strFig '.pdf']);
end

%% decoding plot (hit/miss HCAR/random)
matMeanCorrectDecoded = nan(4,6); %rand/miss rand/hit HCAR/miss HCAR/hit
matSdCorrectDecoded = nan(4,6);

%calc mean
matMeanCorrectDecoded(4,:) = nanmean(matMixBootstrappedDecodingOutput1(:,:,1),1); %HCAR hit
matMeanCorrectDecoded(3,:) = nanmean(matMixBootstrappedDecodingOutput1(:,:,2),1); %HCAR miss
matMeanCorrectDecoded(2,:) = nanmean(matMixBootstrappedDecodingOutput2(:,:,1),1); %rand hit
matMeanCorrectDecoded(1,:) = nanmean(matMixBootstrappedDecodingOutput2(:,:,2),1); %rand miss

%calc sd
matSdCorrectDecoded(4,:) = nanstd(matMixBootstrappedDecodingOutput1(:,:,1),[],1); %HCAR hit
matSdCorrectDecoded(3,:) = nanstd(matMixBootstrappedDecodingOutput1(:,:,2),[],1); %HCAR miss
matSdCorrectDecoded(2,:) = nanstd(matMixBootstrappedDecodingOutput2(:,:,1),[],1); %rand hit
matSdCorrectDecoded(1,:) = nanstd(matMixBootstrappedDecodingOutput2(:,:,2),[],1); %rand miss

%calc over contrasts
matMeanOverContrasts(4,:) = nanmean(matMixBootstrappedDecodingOutput1(:,:,1),2);
matMeanOverContrasts(3,:) = nanmean(matMixBootstrappedDecodingOutput1(:,:,2),2);
matMeanOverContrasts(2,:) = nanmean(matMixBootstrappedDecodingOutput2(:,:,1),2);
matMeanOverContrasts(1,:) = nanmean(matMixBootstrappedDecodingOutput2(:,:,2),2);

%set vars
vecLineX = 1:6;
vecWindowInv = 6:-1:1;
vecX = [vecLineX vecWindowInv];

%plot
intN = size(matMixBootstrappedDecodingOutput1,1);
h=figure;
for intPlotType=1:4
	
	vecMeanTrace = matMeanCorrectDecoded(intPlotType,:);
	vecSE = matSdCorrectDecoded(intPlotType,:)./sqrt(intN);
	
	if intPlotType == 1
		vecColorFill = [1.0 0.7 0.7];
		vecColorLine = [1 0 0];
	elseif intPlotType == 2
		vecColorFill = [0.7 1.0 0.7];
		vecColorLine = [0 1 0];
	elseif intPlotType == 3
		vecColorFill = [1.0 0.7 1.0];
		vecColorLine = [1 0 1];
	elseif intPlotType == 4
		vecColorFill = [0.7 1.0 1.0];
		vecColorLine = [0 1 1];
	end
	
	vecMinTrace = vecMeanTrace-vecSE;
	vecMaxTrace = vecMeanTrace+vecSE;
	vecY = [vecMinTrace vecMaxTrace(vecWindowInv)];
	
	%plot
	hold on
	fill(vecX,vecY,vecColorFill,'EdgeColor','none');
	plot(vecLineX,vecMeanTrace,'-','LineWidth',2,'Color',vecColorLine);
	
	
	%scatter mean
	dblX = max(vecX)+0.5+0.2*intPlotType;
	dblY = mean(matMeanOverContrasts(intPlotType,:),2);
	scatter(dblX,dblY,[],vecColorLine);
	dblErr = std(matMeanOverContrasts(intPlotType,:),[],2)./sqrt(intN);
	errorbar(dblX,dblY,dblErr,'Color',vecColorLine);
	
	hold off
	%end
end
legend(gca,'Miss Random','Hit Random','Miss HCAR','Hit HCAR','Location', 'Best')
set(gca,'XTick',[1 2 3 4 5 6],'XTickLabel',[0 0.5 2 8 32 100]);
xlabel('Stimulus Contrast (%)')
ylabel('Decoding performance (proportion correct)')
ylim([0 1])
xlim([0.5 7.5])

hold on
plot(get(gca,'XLim'),[1/21 1/21],'k--')
hold off

[h,p12]=ttest(matMeanOverContrasts(1,:),matMeanOverContrasts(2,:));
[h,p13]=ttest(matMeanOverContrasts(1,:),matMeanOverContrasts(3,:));
[h,p14]=ttest(matMeanOverContrasts(1,:),matMeanOverContrasts(4,:));
[h,p23]=ttest(matMeanOverContrasts(2,:),matMeanOverContrasts(3,:));
[h,p24]=ttest(matMeanOverContrasts(2,:),matMeanOverContrasts(4,:));
[h,p34]=ttest(matMeanOverContrasts(3,:),matMeanOverContrasts(4,:));
[h,dblGroup]=ttest([matMeanOverContrasts(1,:) matMeanOverContrasts(2,:)],[matMeanOverContrasts(3,:) matMeanOverContrasts(4,:)]);
[h,dblHitMiss]=ttest([matMeanOverContrasts(1,:) matMeanOverContrasts(3,:)],[matMeanOverContrasts(2,:) matMeanOverContrasts(4,:)]);

%random
[h,pR0]=ttest(matMixBootstrappedDecodingOutput2(:,1,1),matMixBootstrappedDecodingOutput2(:,1,2));
[h,pR05]=ttest(matMixBootstrappedDecodingOutput2(:,2,1),matMixBootstrappedDecodingOutput2(:,2,2));
[h,pR2]=ttest(matMixBootstrappedDecodingOutput2(:,3,1),matMixBootstrappedDecodingOutput2(:,3,2));
[h,pR8]=ttest(matMixBootstrappedDecodingOutput2(:,4,1),matMixBootstrappedDecodingOutput2(:,4,2));
[h,pR32]=ttest(matMixBootstrappedDecodingOutput2(:,5,1),matMixBootstrappedDecodingOutput2(:,5,2));
[h,pR100]=ttest(matMixBootstrappedDecodingOutput2(:,6,1),matMixBootstrappedDecodingOutput2(:,6,2));

%HCAR
[h,pH0]=ttest(matMixBootstrappedDecodingOutput1(:,1,1),matMixBootstrappedDecodingOutput1(:,1,2));
[h,pH05]=ttest(matMixBootstrappedDecodingOutput1(:,2,1),matMixBootstrappedDecodingOutput1(:,2,2));
[h,pH2]=ttest(matMixBootstrappedDecodingOutput1(:,3,1),matMixBootstrappedDecodingOutput1(:,3,2));
[h,pH8]=ttest(matMixBootstrappedDecodingOutput1(:,4,1),matMixBootstrappedDecodingOutput1(:,4,2));
[h,pH32]=ttest(matMixBootstrappedDecodingOutput1(:,5,1),matMixBootstrappedDecodingOutput1(:,5,2));
[h,pH100]=ttest(matMixBootstrappedDecodingOutput1(:,6,1),matMixBootstrappedDecodingOutput1(:,6,2));

title(sprintf('(1)Red=Miss Random;(2)Green=Hit Random;(3)Purple=Miss HCAR;(4)Cyan=Hit HCAR\nttest p-vals;1-2: %.3f; 1-3: %.3f; 1-4: %.3f; 2-3: %.3f; 2-4: %.3f; 3-4: %.3f; group: %.3f;hit-miss: %.3f\nHit/miss p-vals per contrast; Random, 0: %.3f; 0.5: %.3f; 2: %.3f; 8: %.3f; 32: %.3f; 100: %.3f; HCAR, 0: %.3f; 0.5: %.3f; 2: %.3f; 8: %.3f; 32: %.3f; 100: %.3f;'...
	,p12,p13,p14,p23,p24,p34,dblGroup,dblHitMiss,pR0,pR05,pR2,pR8,pR32,pR100,pH0,pH05,pH2,pH8,pH32,pH100))
strFigTitle = 'decoding_performance';
if boolSavePlots
	intFigCounter = intFigCounter + 1;
	drawnow;
	strFig = sprintf('Meta%d_%s_raw',intFigCounter,strFigTitle);
	export_fig([strFig '.tif']);
	export_fig([strFig '.pdf']);
end

%% heterogeneity over contrasts
%normalize heterogeneity per animal
matContHeteroNorm = zeros(size(matContHetero));
for intAnimal=1:size(matContHetero,3)
	matTempHet = matContHetero(:,:,intAnimal);
	matContHeteroNorm(:,:,intAnimal) = matTempHet ./ mean(matTempHet(:));
end
matHitHet = matContHeteroNorm(:,1,:);
matMissHet = matContHeteroNorm(:,2,:);
vecHitY = mean(matHitHet,3)';
vecHitE = std(matHitHet,[],3)'/sqrt(size(matContHeteroNorm,3));
vecMissY = mean(matMissHet,3)';
vecMissE = std(matMissHet,[],3)'/sqrt(size(matContHeteroNorm,3));

%get sig
vecP = zeros(1,size(matContHeteroNorm,1));
for intC=1:size(matContHeteroNorm,1)
	[h,p,ci] = ttest(matHitHet(intC,1,:),matMissHet(intC,1,:));
	
	%put in vector
	vecP(intC) = p;
end


%pre-compute variables
vecContrasts = [0.5 2 8 32 100];
vecWindow = [1 5];
vecWindowSelect = vecWindow(1):vecWindow(end);
intWL = length(vecWindowSelect);
vecWindowInv = intWL:-1:1;
vecWindowPlotInv = vecWindow(end):-1:vecWindow(1);
vecX = [vecContrasts(vecWindowSelect) vecContrasts(vecWindowPlotInv)];
vecLineX = vecContrasts(vecWindowSelect);

%get data
hHetCon=figure;
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
set(gca,'XTick',vecLineX,'XTickLabel',vecLineX)
title(sprintf('Pop het during stim; p-vals: 0.5: %.3f; 2: %.3f; 8: %.3f; 32: %.3f; 100: %.3f',vecP))
grid on
xlabel('Contrast')
ylabel('Normalized population response heterogeneity')
xlim(vecContrasts(vecWindow))
%ylim([-0.01 0.06])
legend({'SEM','Miss','SEM','Hit'},'Location','Best')
drawnow;
strFigTitle = 'hetero_over_contrasts';
if boolSavePlots
	intFigCounter = intFigCounter + 1;
	drawnow;
	strFig = sprintf('Meta%d_%s_raw',intFigCounter,strFigTitle);
	export_fig([strFig '.tif']);
	export_fig([strFig '.pdf']);
end

%% heterogeneity+dF/F over time
%cellHetTime{intC,intAnimal} = {trial x time}[1=miss] {trial x time}[2=slow] {trial x time}[3=fast]
%normalize heterogeneity per animal

%pre-allocate for later usage
matDiffFS = nan(2,size(cellHetTime,1),size(cellHetTime,2),size(cellHetTime{1,1}{1},2));
for intAct = [0 1];
	if intAct == 0
		cellDataTime = cellHetTime;
	else
		cellDataTime = cellActTime;
	end
	
	cellRawDataMiss = {};
	cellRawDataSlow = {};
	cellRawDataFast = {};
	
	
	for intAnimal=1:size(cellDataTime,2)
		if intAct == 0
			matTempHet = [];
			for intC=1:size(cellDataTime,1)
				matTempHet = [matTempHet;cellDataTime{intC,intAnimal}{1};cellDataTime{intC,intAnimal}{2};cellDataTime{intC,intAnimal}{3}];
			end
			dblMean = mean(matTempHet(:)); %mean over all values of animal
		end
		%apply normalization
		for intC=1:size(cellDataTime,1)
			if intAct == 0
				cellDataTime{intC,intAnimal}{1} = cellDataTime{intC,intAnimal}{1} ./ dblMean;
				cellDataTime{intC,intAnimal}{2} = cellDataTime{intC,intAnimal}{2} ./ dblMean;
				cellDataTime{intC,intAnimal}{3} = cellDataTime{intC,intAnimal}{3} ./ dblMean;
			end
			cellRawDataMiss{intC}(intAnimal,:) = mean(cellDataTime{intC,intAnimal}{1},1);
			cellRawDataSlow{intC}(intAnimal,:) = mean(cellDataTime{intC,intAnimal}{2},1);
			cellRawDataFast{intC}(intAnimal,:) = mean(cellDataTime{intC,intAnimal}{3},1);
			
			%put in output
			matDiffFS(intAct+1,intC,intAnimal,:) =  cellRawDataFast{intC}(intAnimal,:) - cellRawDataSlow{intC}(intAnimal,:);
		end
	end
	
	%general
	vecC = [0 0.5 2 8 32 100];
	vecWindowSecs = [-3 5];
	vecWindow = round(vecWindowSecs*25.4);
	vecLineX = (vecWindow(1):vecWindow(end))/25.4;
	vecWindowInv = length(vecLineX):-1:1;
	vecX = [vecLineX vecLineX(vecWindowInv)];
	
	%plot
	hHetTime = figure;
	for intC=1:length(vecC)
		subplot(2,3,intC);
		for intType=1:3
			if intType == 1 %miss
				matRawData = cellRawDataMiss{intC};
				strType = 'Miss';
				vecColorLine = [1 0 0];
				vecColorFill = [1 0.7 0.7];
			elseif intType == 2 %slow
				matRawData = cellRawDataSlow{intC};
				strType = 'Slow';
				vecColorLine = [1 1 0];
				vecColorFill = [1 1 0.7];
			else %fast
				matRawData = cellRawDataFast{intC};
				strType = 'Fast';
				vecColorLine = [0 1 0];
				vecColorFill = [0.7 1 0.7];
			end
			
			
			%get data
			vecMeanTrace = mean(matRawData,1);
			vecSE = std(matRawData,[],1)./sqrt(size(cellDataTime,2));
			vecMinTrace = vecMeanTrace-vecSE;
			vecMaxTrace = vecMeanTrace+vecSE;
			vecY = [vecMinTrace vecMaxTrace(vecWindowInv)];
			
			%plot
			hold on
			fill(vecX,vecY,vecColorFill,'EdgeColor','none');
			plot(vecLineX,vecMeanTrace,'-','LineWidth',2,'Color',vecColorLine);
			hold off
		end
		
		%labels
		if intAct == 0
			title(sprintf('Pop heterogeneity over time; Contrast %.1f',vecC(intC)))
			ylabel('Normalized population response heterogeneity')
			ylim([0.8 1.6])
			strFigTitle = 'hetero_over_time';
		else
			title(sprintf('Pop dF/F0 over time; Contrast %.1f',vecC(intC)))
			ylabel('Normalized population dF/F0')
			ylim([-0.02 0.06])
			strFigTitle = 'activity_over_time';
		end
		grid on
		xlabel('Time after stim onset (s)')
		xlim(vecWindowSecs)
		drawnow;
	end
	
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
end

%% plot resp type predictability
vecHetMiss = [];
vecHetSlow = [];
vecHetFast = [];
vecActMiss = [];
vecActSlow = [];
vecActFast = [];

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
	for intC=1:size(cellHetTime,1)
		matTemp = [matTemp;cellHetTime{intC,intAnimal}{1};cellHetTime{intC,intAnimal}{2};cellHetTime{intC,intAnimal}{3}];
	end
	dblMean = mean(matTemp(:)); %mean over all values of animal
	
	%apply normalization
	for intC=1:size(cellHetTime,1)
		cellHetTime{intC,intAnimal}{1} = cellHetTime{intC,intAnimal}{1} ./ dblMean;
		cellHetTime{intC,intAnimal}{2} = cellHetTime{intC,intAnimal}{2} ./ dblMean;
		cellHetTime{intC,intAnimal}{3} = cellHetTime{intC,intAnimal}{3} ./ dblMean;
	end
	
	%prep data
	matTempHetM = [];
	matTempHetS = [];
	matTempHetF = [];
	matTempActM = [];
	matTempActS = [];
	matTempActF = [];
	for intC=1:size(cellHetTime,1)
		matTempHetM = [matTempHetM;cellHetTime{intC,intAnimal}{1}];
		matTempHetS = [matTempHetS;cellHetTime{intC,intAnimal}{2}];
		matTempHetF = [matTempHetF;cellHetTime{intC,intAnimal}{3}];
		matTempActM = [matTempActM;cellActTime{intC,intAnimal}{1}];
		matTempActS = [matTempActS;cellActTime{intC,intAnimal}{2}];
		matTempActF = [matTempActF;cellActTime{intC,intAnimal}{3}];
	end
	vecHetMiss(intAnimal) = mean(mean(matTempHetM(:,intStart:intStop),1));
	vecHetSlow(intAnimal) = mean(mean(matTempHetS(:,intStart:intStop),1));
	vecHetFast(intAnimal) = mean(mean(matTempHetF(:,intStart:intStop),1));
	vecActMiss(intAnimal) = mean(mean(matTempActM(:,intStart:intStop),1));
	vecActSlow(intAnimal) = mean(mean(matTempActS(:,intStart:intStop),1));
	vecActFast(intAnimal) = mean(mean(matTempActF(:,intStart:intStop),1));
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
dblMeanActSlow = mean(vecActSlow);
dblMeanActFast = mean(vecActFast);
dblMeanHetMiss = mean(vecHetMiss);
dblMeanHetSlow = mean(vecHetSlow);
dblMeanHetFast = mean(vecHetFast);
intAnimals = length(vecActMiss);

%separability; inter vs intra cluster distance
vecSepRatioActMiss_w_Slow = abs(vecActMiss - dblMeanActSlow) ./ (abs(vecActMiss - dblMeanActMiss) + abs(vecActMiss - dblMeanActSlow));%inter / (inter+intra) distance
vecSepRatioActMiss_w_Fast = abs(vecActMiss - dblMeanActFast) ./ (abs(vecActMiss - dblMeanActMiss) + abs(vecActMiss - dblMeanActFast));
vecSepRatioActSlow_w_Miss = abs(vecActSlow - dblMeanActMiss) ./ (abs(vecActSlow - dblMeanActSlow) + abs(vecActSlow - dblMeanActMiss));
vecSepRatioActSlow_w_Fast = abs(vecActSlow - dblMeanActFast) ./ (abs(vecActSlow - dblMeanActSlow) + abs(vecActSlow - dblMeanActFast));
vecSepRatioActFast_w_Slow = abs(vecActFast - dblMeanActSlow) ./ (abs(vecActFast - dblMeanActFast) + abs(vecActFast - dblMeanActSlow));
vecSepRatioActFast_w_Miss = abs(vecActFast - dblMeanActMiss) ./ (abs(vecActFast - dblMeanActFast) + abs(vecActFast - dblMeanActMiss));
vecSepRatioHetMiss_w_Slow = abs(vecHetMiss - dblMeanHetSlow) ./ (abs(vecHetMiss - dblMeanHetMiss) + abs(vecHetMiss - dblMeanHetSlow));
vecSepRatioHetMiss_w_Fast = abs(vecHetMiss - dblMeanHetFast) ./ (abs(vecHetMiss - dblMeanHetMiss) + abs(vecHetMiss - dblMeanHetFast));
vecSepRatioHetSlow_w_Miss = abs(vecHetSlow - dblMeanHetMiss) ./ (abs(vecHetSlow - dblMeanHetSlow) + abs(vecHetSlow - dblMeanHetMiss));
vecSepRatioHetSlow_w_Fast = abs(vecHetSlow - dblMeanHetFast) ./ (abs(vecHetSlow - dblMeanHetSlow) + abs(vecHetSlow - dblMeanHetFast));
vecSepRatioHetFast_w_Slow = abs(vecHetFast - dblMeanHetSlow) ./ (abs(vecHetFast - dblMeanHetFast) + abs(vecHetFast - dblMeanHetSlow));
vecSepRatioHetFast_w_Miss = abs(vecHetFast - dblMeanHetMiss) ./ (abs(vecHetFast - dblMeanHetFast) + abs(vecHetFast - dblMeanHetMiss));

%act agg
vecActSlowMiss = [vecSepRatioActMiss_w_Slow vecSepRatioActSlow_w_Miss];
vecActSlowFast = [vecSepRatioActFast_w_Slow vecSepRatioActSlow_w_Fast];
vecActFastMiss = [vecSepRatioActMiss_w_Fast vecSepRatioActFast_w_Miss];

%het agg
vecHetSlowMiss = [vecSepRatioHetMiss_w_Slow vecSepRatioHetSlow_w_Miss];
vecHetSlowFast = [vecSepRatioHetFast_w_Slow vecSepRatioHetSlow_w_Fast];
vecHetFastMiss = [vecSepRatioHetMiss_w_Fast vecSepRatioHetFast_w_Miss];

%fig
hFigRTP = figure;

%plot
subplot(1,2,1)
scatter(vecHetMiss,vecActMiss,'rx')
hold on
scatter(vecHetSlow,vecActSlow,'yx')
scatter(vecHetFast,vecActFast,'gx')
hold off
title('Response type predictability')
xlabel('Heterogeneity during 1s preceding stimulus')
ylabel('dF/F during 1s preceding stimulus')

subplot(1,2,2)

dblOffset = 0.1;
plot([0.5 3.5],[0.5 0.5],'k--')
%plot act
hold on
errorbar(1-dblOffset,mean(vecActSlowMiss),std(vecActSlowMiss)/sqrt(length(vecActSlowMiss)),'bx')
errorbar(2-dblOffset,mean(vecActSlowFast),std(vecActSlowFast)/sqrt(length(vecActSlowFast)),'bx')
errorbar(3-dblOffset,mean(vecActFastMiss),std(vecActFastMiss)/sqrt(length(vecActFastMiss)),'bx')
%plot het
errorbar(1+dblOffset,mean(vecHetSlowMiss),std(vecHetSlowMiss)/sqrt(length(vecHetSlowMiss)),'kx')
errorbar(2+dblOffset,mean(vecHetSlowFast),std(vecHetSlowFast)/sqrt(length(vecHetSlowFast)),'kx')
errorbar(3+dblOffset,mean(vecHetFastMiss),std(vecHetFastMiss)/sqrt(length(vecHetFastMiss)),'kx')
hold off
set(gca,'XTick',1:3,'XTickLabel',{'S-M','S-F','F-M'})
ylabel(sprintf('Separability \n(inter/(intra+inter) cluster distance)'))
[h,pActSM] = ttest(vecActSlowMiss,0.5);
[h,pActSF] = ttest(vecActSlowFast,0.5);
[h,pActFM] = ttest(vecActFastMiss,0.5);
[h,pHetSM] = ttest(vecHetSlowMiss,0.5);
[h,pHetSF] = ttest(vecHetSlowFast,0.5);
[h,pHetFM] = ttest(vecHetFastMiss,0.5);
[h,pSM] = ttest(vecActSlowMiss,vecHetSlowMiss);
[h,pSF] = ttest(vecActSlowFast,vecHetSlowFast);
[h,pFM] = ttest(vecActFastMiss,vecHetFastMiss);
title(sprintf(['Sep act; SM,p=%.3f; SF,p=%.3f; FM,p=%.3f\n'...
	'Sep het; SM,p=%.3f; SF,p=%.3f; FM,p=%.3f\n',...
	'Diff act-het; SM,p=%.3f; SF,p=%.3f; FM,p=%.3f'],...
	pActSM,pActSF,pActFM,...
	pHetSM,pHetSF,pHetFM,...
	pSM,pSF,pFM));
ylim([0 1])
xlim([0.5 3.5])
legend('dF/F','Heterogeneity')

if boolSavePlots
	strFigTitle = 'RespTypePredictability';
	intFigCounter = intFigCounter + 1;
	drawnow;
	strFig = sprintf('Meta%d_%s_raw',intFigCounter,strFigTitle);
	export_fig([strFig '.tif']);
	export_fig([strFig '.pdf']);
end


%% resp type decoding
hFigRTD = figure;
for intRespType=1:3
	if intRespType == 1
		strTitle = 'Miss';
		cellColor = {'r','k','k'};
	elseif intRespType == 2
		strTitle = 'Slow';
		cellColor = {'k','y','k'};
	elseif intRespType == 3
		strTitle = 'Fast';
		cellColor = {'k','k','g'};
	end
	subplot(2,2,intRespType)
	
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
	for intAnimal=1:size(matRespDecode,4)
		thetaAct = matRespDecode(1,1,intRespType,intAnimal);
		rhoAct = matRespDecode(1,2,intRespType,intAnimal);
		thetaHet = matRespDecode(2,1,intRespType,intAnimal);
		rhoHet = matRespDecode(2,2,intRespType,intAnimal);
		
		[x,y]=pol2cart(thetaAct,rhoAct);
		scatter(x,y,'bx');
		[x,y]=pol2cart(thetaHet,rhoHet);
		scatter(x,y,'kx');
	end
	hold off;
	xlim([-1 1]);
	ylim([-1 1]);
	title(sprintf('%s trials',strTitle))
end

%plot summary
subplot(2,2,4);
matPerformance = 1-abs(matRespDecodeDist);
vecActPerfomance = matPerformance(1,:);
vecHetPerfomance = matPerformance(2,:);

%statistics
[hH,pH]=ttest(vecHetPerfomance);
[hA,pA]=ttest(vecActPerfomance);
[hHA,pHA]=ttest(vecActPerfomance,vecHetPerfomance);
vecSubHetP = zeros(1,3);
vecSubActP = zeros(1,3);
vecSubDiffP = zeros(1,3);

%plot
dblOffset=0.1;
dblMeanAct=mean(vecActPerfomance);
dblSDAct = std(vecActPerfomance);
errorbar((1+dblOffset),dblMeanAct,dblSDAct/sqrt(length(vecActPerfomance)),'xb')
hold on
dblMeanHet=mean(vecHetPerfomance);
dblSDHet = std(vecHetPerfomance);
errorbar((2+dblOffset),dblMeanHet,dblSDHet/sqrt(length(vecHetPerfomance)),'xk')
plot([0.5 2.5],[0 0],'k--')
for intRespType=1:3
	if intRespType == 1
		strColor = 'r';
	elseif intRespType == 2
		strColor = 'y';
	elseif intRespType == 3
		strColor = 'g';
	end
	vecThisAct = vecActPerfomance(intRespType:3:end);
	vecThisHet = vecHetPerfomance(intRespType:3:end);
	scatter((1-dblOffset)*ones(size(vecThisAct)),vecThisAct,[strColor 'x']);
	scatter((2-dblOffset)*ones(size(vecThisHet)),vecThisHet,[strColor 'x']);
	
	[dummy,vecSubHetP(intRespType)] = ttest(vecThisHet);
	[dummy,vecSubActP(intRespType)] = ttest(vecThisAct);
	[dummy,vecSubDiffP(intRespType)] = ttest(vecThisAct,vecThisHet);
end
hold off
title(sprintf('Normalized resp type decoding; pA=%.3f;pH=%.3f;pHA=%.3f;\nM:pA=%.3f;pH=%.3f;pHA=%.3f S:pA=%.3f;pH=%.3f;pHA=%.3f F:pA=%.3f;pH=%.3f;pHA=%.3f',...
	pA,pH,pHA,vecSubActP(1),vecSubHetP(1),vecSubDiffP(1),vecSubActP(2),vecSubHetP(2),vecSubDiffP(2),vecSubActP(3),vecSubHetP(3),vecSubDiffP(3)))
ylabel('Norm decod perf')
set(gca,'XTick',[1 2],'XTickLabel',{'dF/F0', 'Heterogen'})
ylim([-1 1])
xlim([0.5 2.5])
if boolSavePlots
	strFigTitle = 'RespTypeDecoding';
	intFigCounter = intFigCounter + 1;
	drawnow;
	strFig = sprintf('Meta%d_%s_raw',intFigCounter,strFigTitle);
	export_fig([strFig '.tif']);
	export_fig([strFig '.pdf']);
end


%% difference fast/slow response trials for dF/F + heterogeneity over time
%intAct=1 = het; intAct=2 = dF/F
%matDiffFS(intAct,intC,intAnimal,intT)


%vars
vecC = [0 0.5 2 8 32 100];
vecWindowSecs = [-3 5];
vecWindow = round(vecWindowSecs*25.4);
vecLineX = (vecWindow(1):vecWindow(end))/25.4;
vecWindowInv = length(vecLineX):-1:1;
vecX = [vecLineX vecLineX(vecWindowInv)];

%loop
for intAct=1:2
	if intAct == 1 %miss
		strType = 'Heterogeneity';
		vecColorLine = [0 0 0];
		vecColorFill = [0.7 0.7 0.7];
	else %fast
		strType = 'dF/F0';
		vecColorLine = [0 0 1];
		vecColorFill = [0.7 0.7 1];
	end
	
	%plot
	hFigDFST = figure;
	for intC=1:length(vecC)
		subplot(2,3,intC);
		
		%pre-allocate
		intMaxT = size(matDiffFS,4);
		vecSignificant = zeros(1,intMaxT);
		vecP = zeros(1,intMaxT);
		matCI = zeros(2,intMaxT);
		for intT=1:intMaxT
			%get data
			vecData = matDiffFS(intAct,intC,:,intT);
			
			%test
			[boolH,dblP,vecCI] = ttest(vecData);
			matCI(:,intT) = vecCI;
			vecP(intT) = dblP;
			vecSignificant(intT) = boolH;
		end
		vecMeanTrace = squeeze(mean(matDiffFS(intAct,intC,:,:),3));
		
		%plot curves
		vecMinTrace = vecMeanTrace - squeeze(std(matDiffFS(intAct,intC,:,:),[],3))/sqrt(size(matDiffFS,3));
		vecMaxTrace = vecMeanTrace + squeeze(std(matDiffFS(intAct,intC,:,:),[],3))/sqrt(size(matDiffFS,3));
		vecY = [vecMinTrace' vecMaxTrace(vecWindowInv)'];
		
		%plot
		hold on
		fill(vecX,vecY,vecColorFill,'EdgeColor','none');
		plot(vecLineX,vecMeanTrace,'-','LineWidth',2,'Color',vecColorLine);
		
		%plot significance
		vecLimY = get(gca,'YLim');
		dblPlotY = vecLimY(2) - (vecLimY(2) - vecLimY(1))*(1/20);
		boolLastPlotted = false;
		intCounterT = 1;
		while ~boolLastPlotted
			if vecSignificant(1)
				intNextT = find(vecSignificant==0,1,'first');
				if isempty(intNextT),boolLastPlotted=true;break;end
				plot([vecLineX(intCounterT) vecLineX(intNextT+intCounterT-1)],[dblPlotY dblPlotY]);
			else
				intNextT = find(vecSignificant==1,1,'first');
				if isempty(intNextT),boolLastPlotted=true;break;end
			end
			
			%update counters
			intCounterT = intCounterT + intNextT -1;
			vecSignificant = vecSignificant(intNextT:end);
		end
		hold off
		
		%labels
		if intAct == 1
			title(sprintf('Diff Fast-slow heterogeneity over time; Contrast %.1f',vecC(intC)))
			ylabel('Difference in heterogeneity')
			ylim([-0.1 0.4])
			strFigTitle = 'diffFS_hetero_over_time';
		else
			title(sprintf('Diff Fast-slow dF/F0 over time; Contrast %.1f',vecC(intC)))
			ylabel('Difference in dF/F0')
			ylim([-0.01 0.04])
			strFigTitle = 'diffFS_activity_over_time';
		end
		grid on
		xlabel('Time after stim onset (s)')
		xlim(vecWindowSecs)
		drawnow;
	end
	
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
end

%% decoding vs heterogeneity
%decoding vs contrast; decoding vs RT; decoding vs hetero-pre; decoding vs
%hetero-dur; decoding vs act-pre; decoding vs act-dur
hDecode2 = figure;

%contrast
subplot(2,3,1);
vecContrasts = unique(cellHD_Contrast{1});
vecContrastsPlot = [0.0025; vecContrasts(2:end)];
intContrasts = length(vecContrasts);
matDecCon = nan(intAnimals,intContrasts,2);
matDecConP = nan(intContrasts,2);
strTitle = sprintf('\n');
for intContrast=1:intContrasts
	for intAnimal=1:intAnimals
		indContrast = cellHD_Contrast{intAnimal}(:,1) == vecContrasts(intContrast);
		matDecCon(intAnimal,intContrast,1) = mean(cellHD_DecodingAccuracy{intAnimal}(indContrast,1));
		matDecCon(intAnimal,intContrast,2) = mean(cellHD_DecodingAccuracy{intAnimal}(indContrast,2));
	end
	[h,matDecConP(intContrast,1)] = ttest(matDecCon(:,intContrast,1),0.25);
	[h,matDecConP(intContrast,2)] = ttest(matDecCon(:,intContrast,2),0.25);
	strTitle = [strTitle sprintf(';C%.3f: HCN-P=%.3f;HUN-P=%.3f',vecContrasts(intContrast),matDecConP(intContrast,1),matDecConP(intContrast,2))];
end
vecDecConMeanHCN = mean(matDecCon(:,:,1),1);
vecDecConSDHCN = std(matDecCon(:,:,1),[],1);
vecDecConMeanHUN = mean(matDecCon(:,:,2),1);
vecDecConSDHUN = std(matDecCon(:,:,2),[],1);

errorbar(vecContrastsPlot,vecDecConMeanHCN,vecDecConSDHCN/sqrt(intAnimals),'r');
hold on
errorbar(vecContrastsPlot,vecDecConMeanHUN,vecDecConSDHUN/sqrt(intAnimals),'b');
plot([vecContrastsPlot(1) vecContrastsPlot(end)],[0.25 0.25],'--k')
hold off
title(['Cross-validated orientation decoding; red=HCN;blue=HUN' strTitle])
xlabel('Stimulus Contrast')
ylabel('Orientation decoding performance')
set(gca,'XScale','log','YScale','linear')

%RT; fast/slow/miss?
subplot(2,3,4);


%remove catch trials & normalize activity
for intAnimal=1:intAnimals
	%remove catch trials	
	cellHD_DecodingAccuracy2{intAnimal} = reshape(cellHD_DecodingAccuracy{intAnimal}(cellHD_Contrast{intAnimal} ~= 0),sum((cellHD_Contrast{intAnimal}(:,1) ~= 0)),2);
	cellHD_Orientation2{intAnimal} = reshape(cellHD_Orientation{intAnimal}(cellHD_Contrast{intAnimal} ~= 0),sum((cellHD_Contrast{intAnimal}(:,1) ~= 0)),2);
	cellHD_Contrast2{intAnimal} = reshape(cellHD_Contrast{intAnimal}(cellHD_Contrast{intAnimal} ~= 0),sum((cellHD_Contrast{intAnimal}(:,1) ~= 0)),2);
	cellHD_ReactionTime2{intAnimal} = reshape(cellHD_ReactionTime{intAnimal}(cellHD_Contrast{intAnimal} ~= 0),sum((cellHD_Contrast{intAnimal}(:,1) ~= 0)),2);
	cellHD_HeterogeneityPre2{intAnimal} = reshape(cellHD_HeterogeneityPre{intAnimal}(cellHD_Contrast{intAnimal} ~= 0),sum((cellHD_Contrast{intAnimal}(:,1) ~= 0)),2);
	cellHD_HeterogeneityDuring2{intAnimal} = reshape(cellHD_HeterogeneityDuring{intAnimal}(cellHD_Contrast{intAnimal} ~= 0),sum((cellHD_Contrast{intAnimal}(:,1) ~= 0)),2);
	cellHD_ActivityPre2{intAnimal} = reshape(cellHD_ActivityPre{intAnimal}(cellHD_Contrast{intAnimal} ~= 0),sum((cellHD_Contrast{intAnimal}(:,1) ~= 0)),2);
	cellHD_ActivityDuring2{intAnimal} = reshape(cellHD_ActivityDuring{intAnimal}(cellHD_Contrast{intAnimal} ~= 0),sum((cellHD_Contrast{intAnimal}(:,1) ~= 0)),2);
	
	%norm act
	cellHD_HeterogeneityPre2{intAnimal} = cellHD_HeterogeneityPre2{intAnimal}/mean(cellHD_HeterogeneityPre2{intAnimal}(:));
	cellHD_HeterogeneityDuring2{intAnimal} = cellHD_HeterogeneityDuring2{intAnimal}/mean(cellHD_HeterogeneityDuring2{intAnimal}(:));
	%cellHD_ActivityPre2{intAnimal} = cellHD_ActivityPre2{intAnimal}/mean(cellHD_ActivityPre2{intAnimal}(:));
	%cellHD_ActivityDuring2{intAnimal} = cellHD_ActivityDuring2{intAnimal}/mean(cellHD_ActivityDuring2{intAnimal}(:));
	
end

%dur-hetero
vecHetDurCorrHCN = nan(1,intAnimals);
vecHetDurIncorrHCN = nan(1,intAnimals);
vecHetDurCorrHUN = nan(1,intAnimals);
vecHetDurIncorrHUN = nan(1,intAnimals);
for intAnimal=1:intAnimals
	vecHetDurCorrHCN(intAnimal) = mean(cellHD_HeterogeneityDuring2{intAnimal}(cellHD_DecodingAccuracy2{intAnimal}(:,1)==1,1));
	vecHetDurIncorrHCN(intAnimal) = mean(cellHD_HeterogeneityDuring2{intAnimal}(cellHD_DecodingAccuracy2{intAnimal}(:,1)==0,1));
	vecHetDurCorrHUN(intAnimal) = mean(cellHD_HeterogeneityDuring2{intAnimal}(cellHD_DecodingAccuracy2{intAnimal}(:,2)==1,2));
	vecHetDurIncorrHUN(intAnimal) = mean(cellHD_HeterogeneityDuring2{intAnimal}(cellHD_DecodingAccuracy2{intAnimal}(:,2)==0,2));
end
%t-tests
[h,p12] = ttest(vecHetDurCorrHUN,vecHetDurIncorrHUN);
[h,p13] = ttest(vecHetDurCorrHUN,vecHetDurCorrHCN);
[h,p24] = ttest(vecHetDurIncorrHUN,vecHetDurIncorrHCN);
[h,p34] = ttest(vecHetDurCorrHCN,vecHetDurIncorrHCN);
[h,p12_34] = ttest(vecHetDurCorrHUN-vecHetDurIncorrHUN,vecHetDurCorrHCN-vecHetDurIncorrHCN);
[h,pCI] = ttest([vecHetDurCorrHUN vecHetDurCorrHCN],[vecHetDurIncorrHUN vecHetDurIncorrHCN]);
subplot(2,3,2);
errorbar([3 4],[mean(vecHetDurCorrHCN) mean(vecHetDurIncorrHCN)],[std(vecHetDurCorrHCN) std(vecHetDurIncorrHCN)]/sqrt(intAnimals),'Linestyle','none','Marker','x','Color','r');
hold on
errorbar([1 2],[mean(vecHetDurCorrHUN) mean(vecHetDurIncorrHUN)],[std(vecHetDurCorrHUN) std(vecHetDurIncorrHUN)]/sqrt(intAnimals),'Linestyle','none','Marker','x','Color','b');
hold off
set(gca,'XTick',1:4,'XTickLabel',{'Corr HUN','Incorr HUN','Corr HCN','Incorr HCN'})
title(sprintf('Heterogeneity during stim correlated with decoding performance;n=8 animals\nT-tests; 1-2,p=%.3f; 1-3,p=%.3f; 2-4,p=%.3f; 3-4,p=%.3f; d12-d34,p=%.3f; C-I,p=%.3f',p12,p13,p24,p34,p12_34,pCI))
ylabel('Normalized heterogeneity during stimulus')

%pre-hetero
vecHetPreCorrHCN = nan(1,intAnimals);
vecHetPreIncorrHCN = nan(1,intAnimals);
vecHetPreCorrHUN = nan(1,intAnimals);
vecHetPreIncorrHUN = nan(1,intAnimals);
for intAnimal=1:intAnimals
	vecHetPreCorrHCN(intAnimal) = mean(cellHD_HeterogeneityPre2{intAnimal}(cellHD_DecodingAccuracy2{intAnimal}(:,1)==1,1));
	vecHetPreIncorrHCN(intAnimal) = mean(cellHD_HeterogeneityPre2{intAnimal}(cellHD_DecodingAccuracy2{intAnimal}(:,1)==0,1));
	vecHetPreCorrHUN(intAnimal) = mean(cellHD_HeterogeneityPre2{intAnimal}(cellHD_DecodingAccuracy2{intAnimal}(:,2)==1,2));
	vecHetPreIncorrHUN(intAnimal) = mean(cellHD_HeterogeneityPre2{intAnimal}(cellHD_DecodingAccuracy2{intAnimal}(:,2)==0,2));
end
%t-tests
[h,p12] = ttest(vecHetPreCorrHUN,vecHetPreIncorrHUN);
[h,p13] = ttest(vecHetPreCorrHUN,vecHetPreCorrHCN);
[h,p24] = ttest(vecHetPreIncorrHUN,vecHetPreIncorrHCN);
[h,p34] = ttest(vecHetPreCorrHCN,vecHetPreIncorrHCN);
[h,p12_34] = ttest(vecHetPreCorrHUN-vecHetPreIncorrHUN,vecHetPreCorrHCN-vecHetPreIncorrHCN);
[h,pCI] = ttest([vecHetPreCorrHUN vecHetPreCorrHCN],[vecHetPreIncorrHUN vecHetPreIncorrHCN]);
subplot(2,3,3);
errorbar([3 4],[mean(vecHetPreCorrHCN) mean(vecHetPreIncorrHCN)],[std(vecHetPreCorrHCN) std(vecHetPreIncorrHCN)]/sqrt(intAnimals),'Linestyle','none','Marker','x','Color','r');
hold on
errorbar([1 2],[mean(vecHetPreCorrHUN) mean(vecHetPreIncorrHUN)],[std(vecHetPreCorrHUN) std(vecHetPreIncorrHUN)]/sqrt(intAnimals),'Linestyle','none','Marker','x','Color','b');
hold off
set(gca,'XTick',1:4,'XTickLabel',{'Corr HUN','Incorr HUN','Corr HCN','Incorr HCN'})
title(sprintf('Heterogeneity preceding stim correlated with decoding performance;n=8 animals\nT-tests; 1-2,p=%.3f; 1-3,p=%.3f; 2-4,p=%.3f; 3-4,p=%.3f; d12-d34,p=%.3f; C-I,p=%.3f',p12,p13,p24,p34,p12_34,pCI))
ylabel('Normalized heterogeneity preceding stimulus')

%dur-act
vecActDurCorrHCN = nan(1,intAnimals);
vecActDurIncorrHCN = nan(1,intAnimals);
vecActDurCorrHUN = nan(1,intAnimals);
vecActDurIncorrHUN = nan(1,intAnimals);
for intAnimal=1:intAnimals
	vecActDurCorrHCN(intAnimal) = mean(cellHD_ActivityDuring2{intAnimal}(cellHD_DecodingAccuracy2{intAnimal}(:,1)==1,1));
	vecActDurIncorrHCN(intAnimal) = mean(cellHD_ActivityDuring2{intAnimal}(cellHD_DecodingAccuracy2{intAnimal}(:,1)==0,1));
	vecActDurCorrHUN(intAnimal) = mean(cellHD_ActivityDuring2{intAnimal}(cellHD_DecodingAccuracy2{intAnimal}(:,2)==1,2));
	vecActDurIncorrHUN(intAnimal) = mean(cellHD_ActivityDuring2{intAnimal}(cellHD_DecodingAccuracy2{intAnimal}(:,2)==0,2));
end
%t-tests
[h,p12] = ttest(vecActDurCorrHUN,vecActDurIncorrHUN);
[h,p13] = ttest(vecActDurCorrHUN,vecActDurCorrHCN);
[h,p24] = ttest(vecActDurIncorrHUN,vecActDurIncorrHCN);
[h,p34] = ttest(vecActDurCorrHCN,vecActDurIncorrHCN);
[h,p12_34] = ttest(vecActDurCorrHUN-vecActDurIncorrHUN,vecActDurCorrHCN-vecActDurIncorrHCN);
[h,pCI] = ttest([vecActDurCorrHUN vecActDurCorrHCN],[vecActDurIncorrHUN vecActDurIncorrHCN]);
subplot(2,3,5);
errorbar([3 4],[mean(vecActDurCorrHCN) mean(vecActDurIncorrHCN)],[std(vecActDurCorrHCN) std(vecActDurIncorrHCN)]/sqrt(intAnimals),'Linestyle','none','Marker','x','Color','r');
hold on
errorbar([1 2],[mean(vecActDurCorrHUN) mean(vecActDurIncorrHUN)],[std(vecActDurCorrHUN) std(vecActDurIncorrHUN)]/sqrt(intAnimals),'Linestyle','none','Marker','x','Color','b');
hold off
set(gca,'XTick',1:4,'XTickLabel',{'Corr HUN','Incorr HUN','Corr HCN','Incorr HCN'})
title(sprintf('Activity during stim correlated with decoding performance;n=8 animals\nT-tests; 1-2,p=%.3f; 1-3,p=%.3f; 2-4,p=%.3f; 3-4,p=%.3f; d12-d34,p=%.3f; C-I,p=%.3f',p12,p13,p24,p34,p12_34,pCI))
ylabel('Normalized activity during stimulus')

%pre-act
vecActPreCorrHCN = nan(1,intAnimals);
vecActPreIncorrHCN = nan(1,intAnimals);
vecActPreCorrHUN = nan(1,intAnimals);
vecActPreIncorrHUN = nan(1,intAnimals);
for intAnimal=1:intAnimals
	vecActPreCorrHCN(intAnimal) = mean(cellHD_ActivityPre2{intAnimal}(cellHD_DecodingAccuracy2{intAnimal}(:,1)==1,1));
	vecActPreIncorrHCN(intAnimal) = mean(cellHD_ActivityPre2{intAnimal}(cellHD_DecodingAccuracy2{intAnimal}(:,1)==0,1));
	vecActPreCorrHUN(intAnimal) = mean(cellHD_ActivityPre2{intAnimal}(cellHD_DecodingAccuracy2{intAnimal}(:,2)==1,2));
	vecActPreIncorrHUN(intAnimal) = mean(cellHD_ActivityPre2{intAnimal}(cellHD_DecodingAccuracy2{intAnimal}(:,2)==0,2));
end
%t-tests
[h,p12] = ttest(vecActPreCorrHUN,vecActPreIncorrHUN);
[h,p13] = ttest(vecActPreCorrHUN,vecActPreCorrHCN);
[h,p24] = ttest(vecActPreIncorrHUN,vecActPreIncorrHCN);
[h,p34] = ttest(vecActPreCorrHCN,vecActPreIncorrHCN);
[h,p12_34] = ttest(vecActPreCorrHUN-vecActPreIncorrHUN,vecActPreCorrHCN-vecActPreIncorrHCN);
[h,pCI] = ttest([vecActPreCorrHUN vecActPreCorrHCN],[vecActPreIncorrHUN vecActPreIncorrHCN]);
subplot(2,3,6);
errorbar([3 4],[mean(vecActPreCorrHCN) mean(vecActPreIncorrHCN)],[std(vecActPreCorrHCN) std(vecActPreIncorrHCN)]/sqrt(intAnimals),'Linestyle','none','Marker','x','Color','r');
hold on
errorbar([1 2],[mean(vecActPreCorrHUN) mean(vecActPreIncorrHUN)],[std(vecActPreCorrHUN) std(vecActPreIncorrHUN)]/sqrt(intAnimals),'Linestyle','none','Marker','x','Color','b');
hold off
set(gca,'XTick',1:4,'XTickLabel',{'Corr HUN','Incorr HUN','Corr HCN','Incorr HCN'})
title(sprintf('Activity preceding stim correlated with decoding performance;n=8 animals\nT-tests; 1-2,p=%.3f; 1-3,p=%.3f; 2-4,p=%.3f; 3-4,p=%.3f; d12-d34,p=%.3f; C-I,p=%.3f',p12,p13,p24,p34,p12_34,pCI))
ylabel('Normalized activity preceding stimulus')

strFigTitle = 'decoding_dependence';
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