clear all;
%close all;
if isfolder('F:\Drive\MontijnHeimel_TimeseriesZeta')
	strPath = 'F:\Drive\MontijnHeimel_TimeseriesZeta';
else
	strPath = 'C:\Drive\MontijnHeimel_TimeseriesZeta';
end
strDataPath = fullfile(strPath,'\Data\');
strFigPath = fullfile(strPath,'\Figs\');

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

matSignifZetaUniStitch = [];
matSignifZetaUniNoStitch = [];
matSignifZetaLinStitch = [];
matSignifZetaLinNoStitch = [];
			
matSignifTtest = [];
matSignifAnova = [];

cellTtestP = [];
cellAnovaP = [];
cellNumSpikes = [];
cellUniStitchP = [];
cellUniNoStitchP = [];
cellLinStitchP = [];
cellLinNoStitchP = [];

cellZetaOldP = [];
intIdxNpx = 0;
intIdx = 0;
for intArea=[1 2 8]%1:numel(cellUniqueAreas)
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
		sDir=dir([strDataPath 'ZetaDataAnova' strRunType 'Resamp*.mat']);
		intFiles=numel(sDir);
		for intFile=1:intFiles
			strFile = sDir(intFile).name;
			
            sLoad=load([strDataPath strFile]);
			
			%remove cells with too few spikes
            indRemCells = sLoad.vecNumSpikes < 3;
            cellFields = fieldnames(sLoad);
            for intField=1:numel(cellFields)
            varData = sLoad.(cellFields{intField});
            sLoad.(cellFields{intField}) = varData(~indRemCells);
            end

			matSignifZetaOld(intIdx,intRandType) = sum(sLoad.vecZetaP_old<0.05)/numel(sLoad.vecZetaP_old);
			
			matSignifZetaUniStitch(intIdx,intRandType) = sum(sLoad.vecZetaP_UniStitch<0.05)/numel(sLoad.vecZetaP_old);
			matSignifZetaUniNoStitch(intIdx,intRandType) = sum(sLoad.vecZetaP_UniNoStitch<0.05)/numel(sLoad.vecZetaP_old);
			matSignifZetaLinStitch(intIdx,intRandType) = sum(sLoad.vecZetaP_LinStitch<0.05)/numel(sLoad.vecZetaP_old);
			matSignifZetaLinNoStitch(intIdx,intRandType) = sum(sLoad.vecZetaP_LinNoStitch<0.05)/numel(sLoad.vecZetaP_old);
			
			matSignifTtest(intIdx,intRandType) = sum(sLoad.vecTtestP<0.05)/numel(sLoad.vecTtestP);
			matSignifAnova(intIdx,intRandType) = sum(sLoad.vecAnovaP<0.05)/numel(sLoad.vecAnovaP);
			
			cellTtestP{intIdx,intRandType} = sLoad.vecTtestP;
			cellNumSpikes{intIdx,intRandType} = sLoad.vecNumSpikes;
			cellZetaOldP{intIdx,intRandType} = sLoad.vecZetaP_old;
			
			cellUniStitchP{intIdx,intRandType} = sLoad.vecZetaP_UniStitch;
			cellUniNoStitchP{intIdx,intRandType} = sLoad.vecZetaP_UniNoStitch;
			cellLinStitchP{intIdx,intRandType} = sLoad.vecZetaP_LinStitch;
			cellLinNoStitchP{intIdx,intRandType} = sLoad.vecZetaP_LinNoStitch;
			
			cellAnovaP{intIdx,intRandType} = sLoad.vecAnovaP;
		end
	end
	
	%plot ROC
    matAUCp = [];
	matAUC_dprime = [];
	if size(cellTtestP,1) >= intIdx && ~isempty(cellTtestP{intIdx,1})
		intIdxNpx = intIdxNpx + 1;
		figure;
		%vecH(intIdxNpx) = subplot(4,3,intIdxNpx);
		subplot(2,3,1)
		maxfig;
		hold on;
		cellColorPlot = {'k',[0.7 0 0.7],'b','r',[0.1 0 0.9],[0.8 0 0.2],[0.3 0 0.7]};
        cellNames = {'Ttest','Anova','ZetaOld','UniStitch','UniNoStitch','LinStitch','LinNoStitch'};
        cellLegend = cellNames;
        vecPlotOrder = [1 2 4 5];
        vecAUC = nan(1,numel(cellLegend));
		for intPlotType=vecPlotOrder
            eval(['cellData = cell' cellNames{intPlotType} 'P;']);
			
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
                matAUCp(intCompAuc,intPlotType) = 2*normcdf(abs(z),'upper');%1 - abs(normcdf(z)-normcdf(-z));
				
				v1 = (vecAUC_se(intCompAuc)*sqrt(2*numel(cellData{1}))).^2;
				v2 = (vecAUC_se(intPlotType)*sqrt(2*numel(cellData{1}))).^2;
				dPrime = (vecAUC(intCompAuc) - vecAUC(intPlotType))/sqrt(0.5*(v1+v2));
				matAUC_dprime(intCompAuc,intPlotType) = dPrime;
			end
			cellLegend{intPlotType} = [cellLegend{intPlotType} sprintf('=%.3f',dblAUC)];
        end

		hold off;
		xlabel('False positive fraction');
		ylabel('Inclusion fraction');
		fixfig;
		legend(cellLegend(vecPlotOrder),'location','best')
		title(sprintf('%s, N=%d'...
			,strArea,size(cellData{intIdx,1},2)))
		
		% save figure
		drawnow;
		strFigName = sprintf('MetaComp_ROC_%s',strArea);
		export_fig([strFigPath strFigName '.tif']);
		export_fig([strFigPath strFigName '.pdf']);

	end
	end

end

