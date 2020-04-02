clear all;
%close all;
strPath = 'D:\Data\Results\OriMetric\Data\';
strFigPath = 'D:\Data\Results\OriMetric\';
cellUniqueAreas = {...
	'V1',...Area 1
	'SC',...Area 2
	'Poisson',...Area 3
	'Retina',...Area 4
	'CaNM',...Area 5
	'CaDG',...Area 6
	'lateral geniculate',...Area 7
	'Primary visual',...Area 8
	'Lateral posterior nucleus',...Area 9
	...'Anterior pretectal nucleus',...Area 10
	'Nucleus of the optic tract',...Area 11
	...'Superior colliculus',...Area 12
	'Anteromedial visual',...Area 13
	'posteromedial visual',...Area 14
	'Anterolateral visual',...Area 15
	'Lateral visual',...Area 16
	'Rostrolateral area',...Area 17
	'Anterior area',...Area 18
	'Subiculum',...Area 19
	'Field CA1',...Area 20
	'Field CA2',...Area 21
	'Field CA3',...Area 22
	'Dentate gyrus',...Area 23
	'Retrosplenial'...Area 24
	};

cellRunStim = {...
	'',...Stim 1 
	'RunDriftingGratings',...Stim 2 
	'RunNaturalMovie'...Stim 3
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
matSignifHz = [];
intIdx = 0;
for intArea=1:numel(cellUniqueAreas)
	if intArea==3 || intArea==4 || intArea==5 || intArea==6,continue;end
	strArea = cellUniqueAreas{intArea}; %V1, SC, Retina, Poisson, GCaMP
	if intArea < 7
		continue;
		vecRunStims = 1;
	elseif intArea > 9
		%continue;
	else
		vecRunStims = [2];
	end
	for intStimType=vecRunStims
		intIdx = intIdx + 1;
		strStim = cellRunStim{intStimType};
		strName = replace([strArea strStim],cellRepStr(:,1),cellRepStr(:,2));
		cellDatasetNames{intIdx} = strName;
		
		%% load data
		strRunType = [strArea strStim];
		sDir=dir([strPath 'ZetaDataMSD' strRunType 'Resamp100*']);
		if isempty(sDir),continue;end
		strFile = sDir(1).name;
		intResampNum = str2double(getFlankedBy(strFile,'Resamp','.mat'));
		sLoad=load([strPath strFile]);
		
		
		intNumN = numel(sLoad.cellDeriv);
		
		vecOnset = nan(1,intNumN);
		vecP = nan(1,intNumN);
		vecOnZ = nan(1,intNumN);
		for intN=1:intNumN
			vecT = sLoad.cellInterpT{intN};
			vecD = sLoad.cellDeriv{intN};
			vecP(intN) = sLoad.vecP(intN);
			if isempty(vecT) || vecP(intN) > 0.01,continue;end
			[dblPeakRate,dblPeakTime,dblPeakWidth,vecPeakStartStop,intPeakLoc,vecPeakStartStopIdx] = ...
				getPeak(vecD,vecT,[0 inf]);
			[dblOnset,dblValue] = getOnset(vecD,vecT,dblPeakTime,[0 inf]);
			if dblOnset > 0.5 || dblOnset < 0.015
				vecD = mean(vecD) - vecD;
				vecD(vecD<0)=0;
				[dblPeakRate,dblPeakTime,dblPeakWidth,vecPeakStartStop,intPeakLoc,vecPeakStartStopIdx] = ...
					getPeak(vecD,vecT,[0.015 1]);
				[dblOnset,dblValue] = getOnset(vecD,vecT,dblPeakTime,[0.015 1]);
			end
			vecOnset(intN) = dblOnset;
			%z
			vecZ = zscore(vecD);
			intOnsetIdx = find(vecT>dblOnset,1);
			dblOnZ = vecZ(intOnsetIdx);
			if isempty(dblOnZ),dblOnZ = nan;end
			vecOnZ(intN) = dblOnZ;

			%{
			%plot
			clf;
			plot(vecT,vecD)
			hold on
			scatter(dblOnset,dblValue)
			hold off
			title(sprintf('%s, N%d; p=%.3f, Onset=%d, z=%.1f',strArea,intN,sLoad.vecP(intN),round(dblOnset*1000),dblOnZ))
			pause
%}
		end
		vecOnset(vecP>0.001) = 0.5;
		vecOnset(vecOnset>0.2) = 0.5;
		vecOnset(vecOnset<0.035) = 0.5;
		
		
		cellOnZ{intIdx} = vecOnZ;
		cellP{intIdx} = vecP;
		cellOnset{intIdx} = vecOnset;
			
	end
end

%% prep data
cellOnset2 = cellOnset;
cellOnZ2 = cellOnZ;
cellDatasetNames2 = cellDatasetNames;
indRem = cellfun(@isempty,cellOnset2);
indRem(10:end) = true;
cellOnset2(indRem) = [];
cellOnZ2(indRem) = [];
cellDatasetNames2(indRem) = [];
cellOnset2 = cellfun(@(x) min(x,0.5),cellOnset2,'uniformoutput',false);
cellOnset3 = cellfun(@(x) x(x<0.5),cellOnset2,'uniformoutput',false);

vecSdOnset = cellfun(@std,cellOnset3);
vecOrderOnset = cellfun(@mean,cellOnset3)-vecSdOnset;
[d,vecReorder]=sort(vecOrderOnset,'ascend');
cellOnset2 = cellOnset2(vecReorder);
cellOnset3 = cellOnset3(vecReorder);
cellDatasetNames2 = cellDatasetNames2(vecReorder);

%scatter(abs(cell2vec(cellOnZ)),cell2vec(cellOnset))

%% plot
figure
subplot(2,2,1);
hold on;
for intRec=1:numel(cellOnset2)
	stairs(1000*sort(cellOnset2{intRec}),linspace(0,1,numel(cellOnset2{intRec})));
end
hold off;
ylabel('Fraction of cells');
xlabel('Onset latency (ms)');
fixfig;

subplot(2,2,2);
hold on;
for intRec=1:numel(cellOnset3)
	vecUseOnsets = cellOnset3{intRec};
	vecUseOnsets(vecUseOnsets>=0.5) = [];
	stairs(sort(vecUseOnsets)*1000,linspace(0,1,numel(vecUseOnsets)));
end
hold off;
legend(cellDatasetNames2,'location','best');
ylabel('Fraction of cells');
xlabel('Onset latency (ms)');
fixfig;


subplot(2,2,3);
hold on
matC=lines(numel(cellDatasetNames2));
vecAllOnsets=[];
cellAreas = {};
for intArea=1:numel(cellDatasetNames2)
vecUseOnsets = cellOnset3{intArea};
	vecUseOnsets(vecUseOnsets>=0.5) = [];
	
vecAllOnsets((end+1):(end+numel(vecUseOnsets))) = vecUseOnsets;
cellAreas((end+1):(end+numel(vecUseOnsets))) = cellfill(cellDatasetNames2{intArea},[1 numel(vecUseOnsets)]);
errorbar(intArea,1000*cellfun(@mean,cellOnset3(intArea)),1000*cellfun(@std,cellOnset3(intArea)),'x','Color',matC(intArea,:))
end

set(gca,'xtick',1:numel(cellDatasetNames2),'xticklabel',cellDatasetNames2);
xlim([0 numel(cellDatasetNames2)+1]);
ylim([0 120])
ylabel('Onset latency (ms)');
fixfig
maxfig;

return
%% save
drawnow;
export_fig(sprintf('%sMetaOnsetsFig.tif',strFigPath));
export_fig(sprintf('%sMetaOnsetsFig.pdf',strFigPath));
print(gcf,'-dpdf', sprintf('%sMetaOnsetsFig.pdf',strFigPath));