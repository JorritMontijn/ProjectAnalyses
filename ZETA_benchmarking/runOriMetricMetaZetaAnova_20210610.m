clear all;
%close all;
strPath = 'F:\Data\Processed\ZETA\Inclusion\';
strFigPath = 'F:\Data\Results\ZETA\Inclusion\';
cellUniqueAreas = {...
	'V1',...Area 1
	'SC',...Area 2
	...'Poisson',...Area 3
	'Retina',...Area 4
	...%'CaNM',...Area 5
	'CaDG',...Area 6
	'lateral geniculate',...Area 7
	'Primary visual',...Area 8
	'Lateral posterior nucleus',...Area 9
	'Anterior pretectal nucleus',...Area 10
	'Nucleus of the optic tract',...Area 11
	'Superior colliculus',...Area 12
	'Anteromedial visual',...Area 13
	'posteromedial visual',...%,...Area 14
	'Anterolateral visual',...Area 15
	'Lateral visual',...Area 16
	'Rostrolateral area',...Area 17
	'Anterior area'};%,...Area 18
cellRunStim = {...
	'',...Stim 1
	'RunDriftingGratings',...Stim 2
	%'RunNaturalMovie'...Stim 3
	};
cellRunRand = {...
	'',...Rand 1
	'-Rand',...Rand 2
	};
cellRepStr = {...
	'RunDriftingGratings','';...
	'RunNaturalMovie','-NM';...
	'lateral geniculate','LGN';...
	'Primary visual','V1';...
	'Lateral posterior nucleus','LP';...
	'Anterior pretectal nucleus','APN';...
	'Nucleus of the optic tract','NOT';...
	'Superior colliculus','SC';...
	'Anteromedial visual','AM';...
	'posteromedial visual','PM';...
	'Anterolateral visual','AL';...
	'Lateral visual','L';...
	'Rostrolateral area','RL';...
	'Anterior area','A';...
	'Subiculum','H-Sub';...
	'Field CA1','H-CA1';...
	'Field CA2','H-CA2';...
	'Field CA3','H-CA3';...
	'Dentate gyrus','H-DG';...
	'Retrosplenial','RetSpl';...
	};

%% prep
cellDatasetNames = {};
matZeta =[];
matNumCells = [];
matSignifZ = [];
matSignifTtest = [];
matSignifAnova = [];

vecAUC_Zeta = [];
vecAUC_Ttest = [];
vecAUC_Anova = [];


cellComputTimeAnova = [];
cellComputTimeZETA = [];
cellTtestP = [];
cellAnovaP = [];
cellNumSpikes = [];
cellZetaP = [];
intIdxNpx = 0;
intIdx = 0;
for intArea=1:numel(cellUniqueAreas)
	strArea = cellUniqueAreas{intArea}; %V1, SC, Retina, Poisson, GCaMP
	if intArea < 5%7
		vecRunStims = 1;
	else
		vecRunStims = 2:numel(cellRunStim);
	end
	for intStimType=vecRunStims
		intIdx = intIdx + 1;
		strStim = cellRunStim{intStimType};
		strName = replace([strArea strStim],cellRepStr(:,1),cellRepStr(:,2));
		cellDatasetNames{intIdx} = strName;
		
		for intRandType=1:2
			%set var
			strRand = cellRunRand{intRandType};
			
			
			%% load data
			strRunType = [strArea strRand strStim];
			sDir=dir([strPath 'ZetaDataAnova' strRunType '*']);
			intFiles=numel(sDir);
			for intFile=1:intFiles
				strFile = sDir(intFile).name;
				intResampNum = str2double(getFlankedBy(strFile,'Resamp','.mat'));
				sLoad=load([strPath strFile]);
				
				%remove nans
				matSignifZ(intIdx,intRandType) = sum(sLoad.vecZetaP<0.05)/numel(sLoad.vecZetaP);
				matSignifTtest(intIdx,intRandType) = sum(sLoad.vecTtestP<0.05)/numel(sLoad.vecTtestP);
				matSignifAnova(intIdx,intRandType) = sum(sLoad.vecAnovaP<0.05)/numel(sLoad.vecAnovaP);
				
				cellNumSpikes{intIdx,intRandType} = sLoad.vecNumSpikes;
				cellZetaP{intIdx,intRandType} = sLoad.vecZetaP;
				cellAnovaP{intIdx,intRandType} = sLoad.vecAnovaP;
				cellTtestP{intIdx,intRandType} = sLoad.vecTtestP;
				cellZetaP{intIdx,intRandType}(isnan(cellZetaP{intIdx,intRandType}))=1;
				cellAnovaP{intIdx,intRandType}(isnan(cellZetaP{intIdx,intRandType}))=1;
				cellTtestP{intIdx,intRandType}(isnan(cellZetaP{intIdx,intRandType}))=1;
				
				sLoad.vecZetaTime;
			end
			
		end
		%% plot ROC
		if size(cellZetaP,1) >= intIdx && ~isempty(cellZetaP{intIdx,1})% && intStimType == 2
			intIdxNpx = intIdxNpx + 1;
			figure;
			%vecH(intIdxNpx) = subplot(4,3,intIdxNpx);
			subplot(2,3,1)
			maxfig;
			hold on;
			vecPlot = 1:3;
			vecAUC = nan(size(vecPlot));
			cellColorPlot = {'k','b','r',[0.7 0 0.7]};
			cellTests = {'Mean-rate t-test','ZETA','ANOVA'};
			for intPlotType=vecPlot
				if intPlotType == 1
					cellData = cellTtestP;
				elseif intPlotType == 2
					cellData = cellZetaP;
				elseif intPlotType == 3
					cellData = cellAnovaP;
				end
				
				vecBothData = cat(2,cellData{intIdx,1},cellData{intIdx,2});
				vecBothLabels = cat(2,zeros(size(cellData{intIdx,1})),ones(size(cellData{intIdx,2})));
				%remove nans
				indRem = isnan(vecBothData);
				vecBothData(indRem) = [];
				vecBothLabels(indRem) = [];
				vecThresholds = sort(vecBothData);
				vecTP = sum(cellData{intIdx,1}<=vecThresholds',2)/sum(~isnan(cellData{intIdx,1}));
				vecFP = sum(cellData{intIdx,2}<=vecThresholds',2)/sum(~isnan(cellData{intIdx,2}));
				
				plot(vecFP,vecTP,'Color',cellColorPlot{intPlotType});
				
				[dblAUC,Aci] = auc(cat(1,vecBothLabels,vecBothData)');
				vecAUC(intPlotType) = dblAUC;
			end
			hold off;
			xlabel('False positive fraction');
			ylabel('Inclusion fraction');
			fixfig;
			legend(cellTests(vecPlot),'location','best')
			title(sprintf('%s, N=%d, AUCs: t=%.3f; Z=%.3f; A=%.3f',strArea,size(cellData{intIdx,1},2),vecAUC))
			
			%save AUCs
			vecAUC_Zeta(intArea) = vecAUC(2);
			vecAUC_Ttest(intArea) = vecAUC(1);
			vecAUC_Anova(intArea) = vecAUC(3);
			
			%shuffled p MIMI
			subplot(2,3,2)
			dblStep = 0.05;
			vecBinsP = 0:dblStep:1;
			vecBinsPlot = (dblStep/2):dblStep:(1-dblStep/2);
			vecCountsHzP = histcounts(cellTtestP{intIdx,2},vecBinsP);
			vecCountsZETAP = histcounts(cellZetaP{intIdx,2},vecBinsP);
			vecCountsAnovaP = histcounts(cellAnovaP{intIdx,2},vecBinsP);
			hold on
			plot(vecBinsPlot,vecCountsHzP,'k');
			plot(vecBinsPlot,vecCountsZETAP,'b');
			plot(vecBinsPlot,vecCountsAnovaP,'r');
			hold off
			title(sprintf('False positive p-value distribution'));
			xlabel('p-value of shuffled control');
			ylabel('Number of cells (count)');
			fixfig;
			legend(cellTests(vecPlot),'location','best')
			
			%normal p
			subplot(2,3,3)
			dblStep = 0.05;
			vecBinsP = 0:dblStep:1;vecBinsP(1)=vecBinsP(1)-eps;vecBinsP(end)=vecBinsP(end)+eps;
			vecBinsPlot = (dblStep/2):dblStep:(1-dblStep/2);
			vecCountsHzP = histcounts(cellTtestP{intIdx,1},vecBinsP);
			vecCountsZETAP = histcounts(cellZetaP{intIdx,1},vecBinsP);
			vecCountsAnovaP = histcounts(cellAnovaP{intIdx,1},vecBinsP);
			hold on
			plot(vecBinsPlot,vecCountsHzP,'k');
			plot(vecBinsPlot,vecCountsZETAP,'b');
			plot(vecBinsPlot,vecCountsAnovaP,'r');
			hold off
			title(sprintf('Real data p-value distribution'));
			xlabel('p-value of real data');
			ylabel('Number of cells (count)');
			fixfig;
			legend(cellTests(vecPlot),'location','best')
			
			% save figure
			drawnow;
			strFigName = sprintf('DG_ZA_ROC_%s',strArea);
			export_fig([strFigPath strFigName '.tif']);
			export_fig([strFigPath strFigName '.pdf']);
			
		end
	end
	
end

%% overall test
%means
dblMeanIncT = mean(matSignifTtest(vecUseAreas,1));
dblMeanIncA = mean(matSignifAnova(vecUseAreas,1));
dblMeanIncZ = mean(matSignifZ(vecUseAreas,1));

dblMeanFaT = mean(matSignifTtest(vecUseAreas,2));
dblMeanFaA = mean(matSignifAnova(vecUseAreas,2));
dblMeanFaZ = mean(matSignifZ(vecUseAreas,2));

dblMeanAucT = mean(vecAUC_Ttest);
dblMeanAucA = mean(vecAUC_Anova);
dblMeanAucZ = mean(vecAUC_Zeta);

%at alpha=0.05
vecUseAreas = find(vecAUC_Zeta>0);
[h,pFivePercAZ]=ttest(matSignifAnova(vecUseAreas,1),matSignifZ(vecUseAreas,1));
[h,pFivePercAT]=ttest(matSignifAnova(vecUseAreas,1),matSignifTtest(vecUseAreas,1));
[h,pFivePercZT]=ttest(matSignifZ(vecUseAreas,1),matSignifTtest(vecUseAreas,1));

%ROC-AUC
[h,pAUC_AZ]=ttest(vecAUC_Anova(vecUseAreas),vecAUC_Zeta(vecUseAreas));
[h,pAUC_AT]=ttest(vecAUC_Anova(vecUseAreas),vecAUC_Ttest(vecUseAreas));
[h,pAUC_ZT]=ttest(vecAUC_Zeta(vecUseAreas),vecAUC_Ttest(vecUseAreas));
		
figure
plot(1:3,[matSignifTtest(:,1) matSignifAnova(:,1) matSignifZ(:,1)],'g')
hold on
plot(4:6,[matSignifTtest(:,2) matSignifAnova(:,2) matSignifZ(:,2)],'r')
hold off

set(gca,'xtick',1:6,'xticklabel',{'T-test','ANOVA','ZETA','T-test','ANOVA','ZETA'})

%FA
[h,pFivePercAZ_FA_alpha]=ttest(matSignifAnova(vecUseAreas,2),0.05);
[h,pFivePercAT_FA_alpha]=ttest(matSignifTtest(vecUseAreas,2),0.05);
[h,pFivePercZT_FA_alpha]=ttest(matSignifZ(vecUseAreas,2),0.05);
	
[h,pFivePercAZ_FA]=ttest(matSignifAnova(vecUseAreas,2),matSignifZ(vecUseAreas,2));
[h,pFivePercAT_FA]=ttest(matSignifAnova(vecUseAreas,2),matSignifTtest(vecUseAreas,2));
[h,pFivePercZT_FA]=ttest(matSignifZ(vecUseAreas,2),matSignifTtest(vecUseAreas,2));

drawnow;
strFigName = sprintf('Fig1G_5perc');
export_fig([strFigPath strFigName '.tif']);
export_fig([strFigPath strFigName '.pdf']);
