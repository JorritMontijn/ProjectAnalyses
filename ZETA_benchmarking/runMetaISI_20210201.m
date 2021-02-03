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
boolUseGumbel = true;
intUseISItype = 2;%0=ISI,1=GISI,2=BISI
cellDatasetNames = {};
matZeta =[];
matNumCells = [];
matSignifZ = [];
matSignifHz = [];
matSignifISI = [];

if intUseISItype == 1
	strTest = 'GISI';
elseif intUseISItype == 2
	strTest = 'BISI';
else
	strTest = 'ISI';
end

cellComputTimeISI = [];
cellComputTimeZETA = [];
cellHzP = [];
cellISIP = [];
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
		sDir=dir([strPath 'Zeta' strTest strRunType 'Resamp100*']);
		intFiles=numel(sDir);
		for intFile=1:intFiles
			strFile = sDir(intFile).name;
			intResampNum = str2double(getFlankedBy(strFile,'Resamp','.mat'));
			sLoad=load([strPath strFile]);
			
			
			matSignifZ(intIdx,intRandType) = sum(sLoad.vecZetaP<0.05)/numel(sLoad.vecZetaP);
			matSignifHz(intIdx,intRandType) = sum(sLoad.vecHzP<0.05)/numel(sLoad.vecHzP);
			matSignifISI(intIdx,intRandType) = sum(sLoad.vecISIP<0.05)/numel(sLoad.vecISIP);
			
			cellComputTimeISI{intIdx,intRandType} = sLoad.vecComputTimeISI;
			cellComputTimeZETA{intIdx,intRandType} = sLoad.vecComputTimeZETA;
			cellHzP{intIdx,intRandType} = sLoad.vecHzP;
			cellISIP{intIdx,intRandType} = sLoad.vecISIP;
			cellNumSpikes{intIdx,intRandType} = sLoad.vecNumSpikes;
			cellZetaP{intIdx,intRandType} = sLoad.vecZetaP;
			
			
		end
		
	end
	%plot ROC
	if size(cellHzP,1) >= intIdx && ~isempty(cellHzP{intIdx,1}) && intStimType == 2
		intIdxNpx = intIdxNpx + 1;
		figure;
		%vecH(intIdxNpx) = subplot(4,3,intIdxNpx);
		subplot(2,2,1)
		maxfig;
		hold on;
		vecAUC = nan(1,3);
		cellColorPlot = {'k','r','b'};
		for intPlotType=1:3
			if intPlotType == 1
				cellData = cellHzP;
			elseif intPlotType == 2
				cellData = cellISIP;
			else
				cellData = cellZetaP;
			end
		intCells = numel(cellData{intIdx,1});
		vecBothData = cat(2,cellData{intIdx,1},cellData{intIdx,2});
		vecBothLabels = cat(2,zeros(size(cellData{intIdx,1})),ones(size(cellData{intIdx,2})));
		vecThresholds = sort(vecBothData);
		vecTP = sum(cellData{intIdx,1}<vecThresholds',2)/intCells;
		vecFP = sum(cellData{intIdx,2}<vecThresholds',2)/intCells;
		
		plot(vecFP,vecTP,cellColorPlot{intPlotType});
		
		 [dblAUC,Aci] = auc(cat(1,vecBothLabels,vecBothData)');
		 vecAUC(intPlotType) = dblAUC;
		end
		hold off;
		xlabel('False positive fraction');
		ylabel('Inclusion fraction');
		fixfig;
		legend({'Mean-rate t-test',[strTest ' shuffle'],'ZETA'},'location','best')
		title(sprintf('%s, AUCs: t=%.3f; ISI=%.3f, Z=%.3f',strArea,vecAUC))
		
		%shuffled p MIMI
		subplot(2,2,2)
		dblStep = 0.05;
		vecBinsP = 0:dblStep:1;
		vecBinsPlot = (dblStep/2):dblStep:(1-dblStep/2);
		vecCountsHzP = histcounts(cellHzP{intIdx,2},vecBinsP);
		vecCountsISIP = histcounts(cellISIP{intIdx,2},vecBinsP);
		vecCountsZETAP = histcounts(cellZetaP{intIdx,2},vecBinsP);
		hold on
		plot(vecBinsPlot,vecCountsHzP,'k');
		plot(vecBinsPlot,vecCountsISIP,'r');
		plot(vecBinsPlot,vecCountsZETAP,'b');
		hold off
		title(sprintf('False positive p-value distribution'));
		xlabel('p-value of shuffled control');
		ylabel('Number of cells (count)');
		fixfig;
		legend({'Mean-rate t-test',[strTest ' shuffle'],'ZETA'},'location','best')
		
		%comput time
		subplot(2,2,3)
		
		vecBinsC = logspace(log10(min(cellComputTimeZETA{intIdx,1}))-1,log10(max(cellComputTimeISI{intIdx,1}))+1,31);
		vecBinsPlotC = vecBinsC;
		vecCountsISIC = histcounts(cellComputTimeISI{intIdx,1},vecBinsC);
		vecCountsZETAC = histcounts(cellComputTimeZETA{intIdx,1},vecBinsC);
		hold on
		plot(vecBinsPlotC(2:end),vecCountsISIC,'r');
		plot(vecBinsPlotC(2:end),vecCountsZETAC,'b');
		hold off
		title(sprintf('Computation time comparison'));
		set(gca,'xscale','log');
		xlabel('Computation time per cell (s)');
		ylabel('Number of cells (count)');
		fixfig;
		legend({[strTest ' shuffle'],'ZETA'},'location','best')
		
		% save figure
		drawnow;
		strFigName = sprintf('%s_ROC_%s',strTest,strArea);
		export_fig([strFigPath strFigName '.tif']);
		export_fig([strFigPath strFigName '.pdf']);

	end
	end

end

