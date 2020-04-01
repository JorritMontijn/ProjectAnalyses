%{
figures to save and take mean over animals
- cellSavePreStim: dF/F detect/no-detect preceding stimulus frames [5 (c) x 2 (d/nd) cell array] x n (animals) with every cell = vector of values
- cellSaveBehavDetect: behavioral detection [5 (c) x n (animals)] cell array with every cell = vector of detect/no detect values (1/0)
- cellSaveBehavRT: behavioral RT [5 (c) x n (animals)] cell array with every cell = vector of RTs
- cellSaveBaseSubtrStimAct: dF/F detect/no-detect baseline-subtracted activity [5 (c) x 2 (d/nd) cell array] x n (animals) with every cell = vector of values
- cellSaveDuringStim: dF/F detect/no-detect during stimulus frames [5 (c) x 2 (d/nd) cell array] x n (animals) with every cell = vector of values
- cellSaveTraceAct: dF/F trace detect/no-detect pref/non-pref [2 (pref/no-pref) x 2 (detect/no-detect) x 6 (c)] x n (animals) with every cell = [Nreps x Nframes] matrix of trace dF/F values
-  AD: cellSaveMatrices: blob matrices in [1 x 4 (detect/no-detect/diff/significance)] cell array with [5 (c) x 31 (bins)] matrix
-  AD: cellSaveDCAENeuronConsistency: [1 x i (neuronal population)] cell array with vector of shuffled correlation values, where first value is actual data
-  AD: cellSaveSignalCorrs: signal correlations [3 (within-DCAE, between, within-non-DCAE) x n (animals/blocks)] with every cell = vector of pairwise correlation values
-  AD: cellSaveNoiseCorrs: noise correlations [3 (within-DCAE, between, within-non-DCAE) x n (animals/blocks)] with every cell = vector of pairwise correlation values
-  AD: cellSaveDCAE: normalized increase in dF/F for detection trials [5 (c) x n (animals/blocks)] cell array with each cell = vector of neurons with normalized dF/F increase for detection trials
-  AD: cellSaveSignalCorrsBiDir: signal correlations [3 (within-DCAE, between, within-non-DCAE) x n (animals/blocks)] with every cell = vector of pairwise correlation values
-  AD: cellSaveNoiseCorrsBiDir: noise correlations [3 (within-DCAE, between, within-non-DCAE) x n (animals/blocks)] with every cell = vector of pairwise correlation values
-  AD: cellSaveNormActDissim: within-group z-scored activation dissimilarity [2 (high/low DCAE) x 2 (detect/no-detect)] x n (animals/blocks) with every cell = vector of dissimilarity values
-  AD: cellSaveDissimCorrITD: inter-trial-distance dependence of assembly consistency [2 (high/low DCAE) x 2 (hits/misses) x n (animals/blocks)] with every cell = [p x 2] matrix with intertrial distance and correlation values for each trial pair
	- extra graph of above with means comparison of hits/misses high/low DCAE

%}
clear all
%% parameters
strFigDir = 'D:\Data\Results\stimdetection\meta';
strDataDir = 'D:\Data\Results\stimdetection';
cellInclude = {'20140207','20140314','20140425','20140430','20140507'};
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
cellDCAENeuronConsistency = {};%: [1 x i (neuronal population)] cell array with vector of shuffled correlation values, where first value is actual data
cellSignalCorrs = {};%: signal correlations [3 (within-DCAE, between, within-non-DCAE) x n (animals/blocks)] with every cell = vector of pairwise correlation values
cellNoiseCorrs = {};%: noise correlations [3 (within-DCAE, between, within-non-DCAE) x n (animals/blocks)] with every cell = vector of pairwise correlation values
cellDCAE = {};%: normalized increase in dF/F for detection trials [5 (c) x n (animals/blocks)] cell array with each cell = vector of neurons with normalized dF/F increase for detection trials
cellSignalCorrsBiDir = {};%: signal correlations [3 (within-DCAE, between, within-non-DCAE) x n (animals/blocks)] with every cell = vector of pairwise correlation values
cellNoiseCorrsBiDir = {};%: noise correlations [3 (within-DCAE, between, within-non-DCAE) x n (animals/blocks)] with every cell = vector of pairwise correlation values
cellNormActDissim = {};%: within-group z-scored activation dissimilarity [2 (high/low DCAE) x 2 (detect/no-detect)] x n (animals/blocks) with every cell = vector of dissimilarity values
cellDissimCorrITD = {};%: inter-trial-distance dependence of assembly consistency [2 (high/low DCAE) x 2 (hits/misses) x n (animals/blocks)] with every cell = [p x 2] matrix with intertrial distance and correlation values for each trial pair



%% create aggregate data structures
sDir = dir(strDataDir);
strOldDir = cd(strDataDir);
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
	end
	
	%check if data part 2
	strFile = sDir(intFile).name;
	strRec = getFlankedBy(strFile,'data_aggregateAD','_');
	intRec = find(strcmp(strRec,cellInclude),1);
	if ~isempty(intRec)
		%load file
		load(strFile);
		
		%update population source
		intNrPops = length(cellSaveDCAENeuronConsistency);
		intAnimal = intAnimal + 1;
		for intPop=1:intNrPops
			vecAnimalSource(end+1) = intAnimal;
			vecSubPopulationSource(end+1) = intPop;
		end
		
		%assign data
		cellMatrices = cat(2,cellMatrices,cellSaveMatrices');%: blob matrices in [1 x 4 (detect/no-detect/diff/significance)] cell array with [5 (c) x 31 (bins)] matrix
		cellDCAENeuronConsistency = cat(1,cellDCAENeuronConsistency,cellSaveDCAENeuronConsistency');%: [1 x i (neuronal population)] cell array with vector of shuffled correlation values, where first value is actual data
		cellSignalCorrs = cat(2,cellSignalCorrs,cellSaveSignalCorrs);%: signal correlations [3 (within-DCAE, between, within-non-DCAE) x n (animals/blocks)] with every cell = vector of pairwise correlation values
		cellNoiseCorrs = cat(2,cellNoiseCorrs,cellSaveNoiseCorrs);%: noise correlations [3 (within-DCAE, between, within-non-DCAE) x n (animals/blocks)] with every cell = vector of pairwise correlation values
		cellDCAE = cat(2,cellDCAE,cellSaveDCAE');%: normalized increase in dF/F for detection trials [5 (c) x n (animals/blocks)] cell array with each cell = vector of neurons with normalized dF/F increase for detection trials
		cellSignalCorrsBiDir = cat(2,cellSignalCorrsBiDir,cellSaveSignalCorrsBiDir);%: signal correlations [3 (within-DCAE, between, within-non-DCAE) x n (animals/blocks)] with every cell = vector of pairwise correlation values
		cellNoiseCorrsBiDir = cat(2,cellNoiseCorrsBiDir,cellSaveNoiseCorrsBiDir);%: noise correlations [3 (within-DCAE, between, within-non-DCAE) x n (animals/blocks)] with every cell = vector of pairwise correlation values
		cellNormActDissim = cat(3,cellNormActDissim,cellSaveNormActDissim);%: within-group z-scored activation dissimilarity [2 (high/low DCAE) x 2 (detect/no-detect)] x n (animals/blocks) with every cell = vector of dissimilarity values
		cellDissimCorrITD = cat(3,cellDissimCorrITD,cellSaveDissimCorrITD);%: inter-trial-distance dependence of assembly consistency [2 (high/low DCAE) x 2 (hits/misses) x n (animals/blocks)] with every cell = [p x 2] matrix with intertrial distance and correlation values for each trial pair
	end
end



%% meta analyses part 1
cd(strFigDir);
close all;
intAnimals = size(cellBehavDetect,2);
vecContrasts = [0 0.5 2 8 32 100];

%% behavioral analyses
%overview graph 0%/100% detection rates with significance calculation
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

%mean detection rate over contrasts (mean = blue, single animals = grey lines)
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
title('Behavioral stimulus detection per animal (grey) and mean (blue)')
if boolSavePlots
	intFigCounter = intFigCounter + 1;
	drawnow;
	strFig = sprintf('Meta%d_behavDetect_raw',intFigCounter);
	export_fig([strFig '.tif']);
	export_fig([strFig '.pdf']);
end

%mean RT over contrasts (mean = blue, single animals = grey lines)
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
title('Behavioral RT for hit trials per animal (grey) and mean (blue)')
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
	
	for intC=1:5
		for intAnimal=1:intAnimals
			matActHit(intC,intAnimal) = mean(cellData{intC,1,intAnimal});
			matActMiss(intC,intAnimal) = mean(cellData{intC,2,intAnimal});
		end
		
		%ttest
		[h,vecP(intC),ci] = ttest2(matActHit(intC,:),matActMiss(intC,:));
	end
	
	%plot
		
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
	title(strTitle)
	grid on
	xlabel('Contrast')
	ylabel(strLabelY)
	xlim(vecContrasts(vecWindow))
	%ylim([-0.01 0.06])
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
				matAct(intAnimal,:) = mean(cellTraceAct{intPrefPop,intDetect,intC,intAnimal},1);
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
			ylim([-0.02 0.06])
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
intBinStepSize = 0.1;
intBinStart = -1; %-1
intBinStop = 2; %2
vecBins = (intBinStart-intBinStepSize):intBinStepSize:(intBinStop+intBinStepSize);
vecTickY = [1 11 21 31];%[6 16 26 36];%[1 6 11 16 20];
vecTickLabelY = vecBins(vecTickY+1);
vecTickLabelX = vecContrasts(2:end);
strLabelY = sprintf('Z-scored dF/F (sigma)');

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
	matMean = mean(matAgg,3);
	if intType == 4, matMean = matMean*size(cellMatrices,2);end
	if intType < 3,imagesc(matMean);else imagesc(matMean);end
	axis xy
	title(strTitle)
	set(gca,'YTick',vecTickY,'YTickLabel',vecTickLabelY)
	set(gca,'XTick',1:length(vecTickLabelX),'XTickLabel',vecTickLabelX)
	ylabel(strLabelY)
	xlabel('Contrast (%)')
	colorbar
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

%% blob consistency
%shuffled confidence interval of DCAE-bins neuron membership consistency + actual data plot per neuronal population
%cellDCAENeuronConsistency = cat(1,cellDCAENeuronConsistency,cellSaveDCAENeuronConsistency');
%: [1 x i (neuronal population)] cell array with vector of shuffled correlation values, where first value is actual data
hConsistencyFigure = figure;
intPopulations = size(cellDCAENeuronConsistency,1);
intIterations = length(cellDCAENeuronConsistency{1})-1;
dblAlpha = 0.05;
intLow = round(intIterations*dblAlpha);
intHigh = round(intIterations*(1-dblAlpha));
vecData = nan(1,intPopulations);
vecX = nan(1,intPopulations);
vecY = nan(1,intPopulations);
vecL = nan(1,intPopulations);
vecH = nan(1,intPopulations);
cellLabelsX = cell(1,intPopulations);
for intPopulation = 1:intPopulations
	intAnimal = vecAnimalSource(intPopulation);
	intSubPop = vecSubPopulationSource(intPopulation);
	
	dblMean = mean(cellDCAENeuronConsistency{intPopulation}(2:end));
	vecShuffSorted = sort(cellDCAENeuronConsistency{intPopulation}(2:end),'ascend');
	dblLowVal = vecShuffSorted(intLow);
	dblHighVal =  vecShuffSorted(intHigh);
	
	vecData(intPopulation) = cellDCAENeuronConsistency{intPopulation}(1);
	vecX(intPopulation) = intPopulation;
	vecY(intPopulation) = dblMean;
	vecL(intPopulation) = dblMean-dblLowVal;
	vecH(intPopulation) = dblHighVal-dblMean;
	cellLabelsX{intPopulation} = [num2str(intAnimal) ' - ' num2str(intSubPop)];
end
errorbar(1.4,mean(vecData),std(vecData)/sqrt(intPopulations),'LineStyle','none','Marker','o')
hold on
scatter(1.6*ones(size(vecData)),vecData, 'rx', 'SizeData',60)
xlim([1 2])
plot(get(gca,'XLim'),[0 0],'k--')
hold off
set(gca,'XTick',[])
set(gca,'XTickLabel',[])
ylabel('Pearson correlation')
title('Neuronal membership consistency in Hits-specific activity bin ')
ylim([-1 1])
if boolSavePlots
	intFigCounter = intFigCounter + 1;
	drawnow;
	strFig = sprintf('Meta%d_DCAEmembership_consistency_raw',intFigCounter);
	export_fig([strFig '.tif']);
	export_fig([strFig '.pdf']);
end

%% signal correlations
%1 x 3 plot; within-DCAE, between, within-non-DCAE neurons
%mean + standard error over animals
%cellSignalCorrs = cat(2,cellSignalCorrs,cellSaveSignalCorrs);
%: signal correlations [3 (within-DCAE, between, within-non-DCAE) x n (populations)] with every cell = vector of pairwise correlation values
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
	strFig = sprintf('Meta%d_noisecorrelations_raw',intFigCounter);
	export_fig([strFig '.tif']);
	export_fig([strFig '.pdf']);
end
vecIncreaseHCAR_SC = matDataSC(1,:)-matDataSC(3,:);

%% noise correlations
%1 x 3 plot; within-DCAE, between, within-non-DCAE neurons
%mean + standard error over animals
%cellNoiseCorrs = cat(2,cellNoiseCorrs,cellSaveNoiseCorrs);
%: noise correlations [3 (within-DCAE, between, within-non-DCAE) x n (populations)] with every cell = vector of pairwise correlation values

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
	drawnow;
	strFig = sprintf('Meta%d_diffcorrelations%d_raw',intFigCounter,boolNormalize);
	export_fig([strFig '.tif']);
	export_fig([strFig '.pdf']);
end

%% neuron-based DCAE distribution
%sum over animals, then normalize peak to 1, over contrasts
%cellDCAE = cat(2,cellDCAE,cellSaveDCAE');
%: normalized increase in dF/F for detection trials [5 (c) x n (populations)] cell array with each cell = vector of neurons with normalized dF/F increase for detection trials
		
%% bi-dir signal correlations
%1 x 3 plot; within-DCAE, between, within-non-DCAE neurons
%mean + standard error over animals
%cellSignalCorrsBiDir = cat(2,cellSignalCorrsBiDir,cellSaveSignalCorrsBiDir);
%: signal correlations [3 (within-DCAE, between, within-non-DCAE) x n (populations)] with every cell = vector of pairwise correlation values
		
%% bi-dir noise correlations
%1 x 3 plot; within-DCAE, between, within-non-DCAE neurons
%mean + standard error over animals
%cellNoiseCorrsBiDir = cat(2,cellNoiseCorrsBiDir,cellSaveNoiseCorrsBiDir);
%: noise correlations [3 (within-DCAE, between, within-non-DCAE) x n (populations)] with every cell = vector of pairwise correlation values
		
%% within-group z-scored activation dissimilarity
%2 (high/low DCAE) x 2 (detect/no-detect)] plot
%mean + standard error over animals
%cellNormActDissim = cat(3,cellNormActDissim,cellSaveNormActDissim);
%: within-group z-scored activation dissimilarity [2 (high/low DCAE) x 2 (hits/misses)] x n (populations) with every cell = vector of dissimilarity values

%1 x 4 plot; 1 data point = mean correlation per population for
%high/low/hits/misses combination; plot is mean +/- st err of data points
hZscoredActivityDissimilarity = figure;
intPopulations = size(cellNormActDissim,3);
matGroup = nan(4,intPopulations);
matData = nan(4,intPopulations);
for intPopulation = 1:intPopulations
	matGroup(:,intPopulation) = 1:4;
	matData(4,intPopulation) = mean(cellNormActDissim{1,1,intPopulation}); %high/hit
	matData(3,intPopulation) = mean(cellNormActDissim{1,2,intPopulation});%high/miss
	matData(2,intPopulation) = mean(cellNormActDissim{2,1,intPopulation});%low/hit
	matData(1,intPopulation) = mean(cellNormActDissim{2,2,intPopulation});%low/miss
end
vecMean = mean(matData,2);
vecErr = std(matData,[],2)/sqrt(intPopulations);

% perform ttests
matP = nan(4,4);
[h,matP(1,2)] = ttest2(matData(1,:),matData(2,:)); %assembly/non-assembly
[h,matP(1,3)] = ttest2(matData(1,:),matData(3,:)); %assembly/non-assembly
[h,matP(1,4)] = ttest2(matData(1,:),matData(4,:)); %assembly/non-assembly
[h,matP(2,3)] = ttest2(matData(2,:),matData(3,:)); %assembly/non-assembly
[h,matP(2,4)] = ttest2(matData(2,:),matData(4,:)); %assembly/non-assembly
[h,matP(3,4)] = ttest2(matData(3,:),matData(4,:)); %assembly/non-assembly

%plot
errorbar((1:4)+0.1,vecMean,vecErr,'ob','LineStyle','none');
hold on
scatter(matGroup(:)-0.1,matData(:),'rx')
hold off
set(gca,'XTick',1:4,'XTickLabel',{'Miss Non-DCAE','Hit Non-DCAE','Miss DCAE','Hit DCAE'})
ylabel('Mean within-group activation dissimilarity')
title(sprintf('T-test p-values, 1-2=%.3f; 1-3=%.3f; 1-4=%.3f; 2-3=%.3f; 2-4=%.3f; 3-4=%.3f;',matP(1,2),matP(1,3),matP(1,4),matP(2,3),matP(2,4),matP(3,4)))
if boolSavePlots
	intFigCounter = intFigCounter + 1;
	drawnow;
	strFig = sprintf('Meta%d_activationdissimilarity_raw',intFigCounter);
	export_fig([strFig '.tif']);
	export_fig([strFig '.pdf']);
end

%% inter-trial-distance dependence of assembly consistency 
%2 (high/low DCAE) x 2 (hits/misses) plot
%aggregate of all data points from all populations
%cellDissimCorrITD = cat(3,cellDissimCorrITD,cellSaveDissimCorrITD);
%: inter-trial-distance dependence of assembly consistency [2 (high/low DCAE) x 2 (hits/misses) x n (populations)] with every cell = [p x 2] matrix with intertrial distance and correlation values for each trial pair


%1 x 4 plot; 1 data point = mean correlation per population for
%high/low/hits/misses combination; plot is mean +/- st err of data points
hAssemblyConsistency = figure;
intPopulations = size(cellDissimCorrITD,3);
matGroup = nan(4,intPopulations);
matData = nan(4,intPopulations);
for intPopulation = 1:intPopulations
	matGroup(:,intPopulation) = 1:4;
	matData(4,intPopulation) = mean(cellDissimCorrITD{1,1,intPopulation}(:,end)); %high/hit
	matData(3,intPopulation) = mean(cellDissimCorrITD{1,2,intPopulation}(:,end));%high/miss
	matData(2,intPopulation) = mean(cellDissimCorrITD{2,1,intPopulation}(:,end));%low/hit
	matData(1,intPopulation) = mean(cellDissimCorrITD{2,2,intPopulation}(:,end));%low/miss
end
vecMean = mean(matData,2);
vecErr = std(matData,[],2)/sqrt(intPopulations);

% perform ttests
matP = nan(4,4);
[h,matP(1,2)] = ttest2(matData(1,:),matData(2,:)); %assembly/non-assembly
[h,matP(1,3)] = ttest2(matData(1,:),matData(3,:)); %assembly/non-assembly
[h,matP(1,4)] = ttest2(matData(1,:),matData(4,:)); %assembly/non-assembly
[h,matP(2,3)] = ttest2(matData(2,:),matData(3,:)); %assembly/non-assembly
[h,matP(2,4)] = ttest2(matData(2,:),matData(4,:)); %assembly/non-assembly
[h,matP(3,4)] = ttest2(matData(3,:),matData(4,:)); %assembly/non-assembly

%plot
errorbar((1:4)+0.1,vecMean,vecErr,'ob','LineStyle','none');
hold on
scatter(matGroup(:)-0.1,matData(:),'rx')
hold off
set(gca,'XTick',1:4,'XTickLabel',{'Miss Non-DCAE','Hit Non-DCAE','Miss DCAE','Hit DCAE'})
ylabel('Mean inter-trial correlation of activation dissimilarity')
title(sprintf('T-test p-values, 1-2=%.3f; 1-3=%.3f; 1-4=%.3f; 2-3=%.3f; 2-4=%.3f; 3-4=%.3f;',matP(1,2),matP(1,3),matP(1,4),matP(2,3),matP(2,4),matP(3,4)))
if boolSavePlots
	intFigCounter = intFigCounter + 1;
	drawnow;
	strFig = sprintf('Meta%d_assemblyconsistency_raw',intFigCounter);
	export_fig([strFig '.tif']);
	export_fig([strFig '.pdf']);
end

%% clean up
cd(strOldDir);