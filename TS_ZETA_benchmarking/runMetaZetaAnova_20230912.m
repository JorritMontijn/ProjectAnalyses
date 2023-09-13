clear all;
%close all;
strPath = 'F:\Data\Processed\ZETAv2\Inclusion\';
strFigPath = 'F:\Data\Results\ZETAv2\Inclusion\';
cellUniqueAreas = {...
	'HeteroPoissonPeak',...Area 1
	'TriPhasic',...Area 2
	'',...Area 3
	'',...Area 4
	'',...Area 5
	'',...Area 6
	'lateral geniculate',...Area 7
	'Primary visual',...Area 8
	};

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
matNumCells = [];
matSignifZetaOld = [];
matSignifZetaNew = [];
matSignifTtest = [];
matSignifAnova = [];

cellTtestP = [];
cellAnovaP = [];
cellNumSpikes = [];
cellZetaNewP = [];
cellZetaOldP = [];
intIdxNpx = 0;
intIdx = 0;
for intArea=2%1:numel(cellUniqueAreas)
	strArea = cellUniqueAreas{intArea}; %V1, SC, Retina, Poisson, GCaMP
	if intArea < 5%7
		vecRunStims = 1;
	else
		vecRunStims = 2:numel(cellRunStim);
	end
	for intStimType=vecRunStims
		intIdx = intIdx + 1;
		strStim = cellRunStim{intStimType};
		cellDatasetNames{intIdx} = [strArea strStim];
		
	for intRandType=1:2
		%set var
		strRand = cellRunRand{intRandType};
		
		
		%% load data
		strRunType = [strArea strRand strStim];
		sDir=dir([strPath 'ZetaDataAnova' strRunType 'Resamp*.mat']);
		intFiles=numel(sDir);
		for intFile=1:intFiles
			strFile = sDir(intFile).name;
			
            sLoad=load([strPath strFile]);
			
			%remove cells with too few spikes
            indRemCells = sLoad.vecNumSpikes < 100;
            cellFields = fieldnames(sLoad);
            for intField=1:numel(cellFields)
            varData = sLoad.(cellFields{intField});
            sLoad.(cellFields{intField}) = varData(~indRemCells);
            end

			matSignifZetaOld(intIdx,intRandType) = sum(sLoad.vecZetaP_old<0.05)/numel(sLoad.vecZetaP_old);
			matSignifZetaNew(intIdx,intRandType) = sum(sLoad.vecZetaP_new<0.05)/numel(sLoad.vecZetaP_new);
			matSignifTtest(intIdx,intRandType) = sum(sLoad.vecTtestP<0.05)/numel(sLoad.vecTtestP);
			matSignifAnova(intIdx,intRandType) = sum(sLoad.vecAnovaP<0.05)/numel(sLoad.vecAnovaP);
			
			cellTtestP{intIdx,intRandType} = sLoad.vecTtestP;
			cellNumSpikes{intIdx,intRandType} = sLoad.vecNumSpikes;
			cellZetaOldP{intIdx,intRandType} = sLoad.vecZetaP_old;
			cellZetaNewP{intIdx,intRandType} = sLoad.vecZetaP_new;
			cellAnovaP{intIdx,intRandType} = sLoad.vecAnovaP;

		end
		
	end
	%plot ROC
    matAUCp = nan(4,4);
	if size(cellTtestP,1) >= intIdx && ~isempty(cellTtestP{intIdx,1})
		intIdxNpx = intIdxNpx + 1;
		figure;
		%vecH(intIdxNpx) = subplot(4,3,intIdxNpx);
		subplot(2,3,1)
		maxfig;
		hold on;
		vecAUC = nan(1,4);
		cellColorPlot = {'k',[0.7 0 0.7],'b','r'};
        cellLegend = {'Mean-rate t-test','ANOVA','Old ZETA','New zeta-test'};
        vecPlotOrder = 1:4;
        for intPlotType=vecPlotOrder
            if intPlotType == 1
      		     cellData = cellTtestP;
            elseif intPlotType == 2
                cellData = cellAnovaP;
            elseif intPlotType == 3
                cellData = cellZetaOldP;
            elseif intPlotType == 4
                cellData = cellZetaNewP;
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

			[dblAUC,Aci,Ase,pAuc] = getAuc(vecTP,vecFP);

			vecAUC(intPlotType) = dblAUC;
            vecAUC_se(intPlotType) = Ase;

            for intCompAuc=1:(intPlotType-1)
                % Observed data
                m0 = vecAUC(intCompAuc) - vecAUC(intPlotType);
                s0 = (vecAUC_se(intCompAuc) + vecAUC_se(intPlotType))/2;
                z = m0/s0;
                matAUCp(intCompAuc,intPlotType) = 1 - abs(normcdf(z)-normcdf(-z));
            end
        end

		hold off;
		xlabel('False positive fraction');
		ylabel('Inclusion fraction');
		fixfig;
		legend(cellLegend(vecPlotOrder),'location','best')
		title(sprintf('%s, N=%d, AUCs: t=%.3f; A=%.3f; Zo=%.3f; Zn=%.3f',strArea,size(cellData{intIdx,1},2),vecAUC))
		
		% save figure
		drawnow;
		strFigName = sprintf('MetaComp_ROC_%s',strArea);
		export_fig([strFigPath strFigName '.tif']);
		export_fig([strFigPath strFigName '.pdf']);

	end
	end

end

