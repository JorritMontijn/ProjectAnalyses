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
	'Superior colliculus',...Area 12
	'Anteromedial visual',...Area 13
	'posteromedial visual',...Area 14
	'Anterolateral visual',...Area 15
	'Lateral visual',...Area 16
	'Rostrolateral area',...Area 17
	...'Anterior area',...Area 18
	...'Subiculum',...Area 19
	'Field CA1',...Area 20
	'Field CA2',...Area 21
	...'Field CA3',...Area 22
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
intRandIdx=1; %1=normal,2=rand
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
		strRand = cellRunRand{intRandIdx};
		strStim = cellRunStim{intStimType};
		strName = replace([strArea strRand strStim],cellRepStr(:,1),cellRepStr(:,2));
		cellDatasetNames{intIdx} = strName;
		
		%% load data
		strRunType = [strArea strRand strStim];
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
			if isempty(vecT) || (vecP(intN) > 0.01 && intRandIdx == 1),
				continue;
			end
			%%
			[dblOnset,dblValueOn,dblBaseVal,dblPeakT] = getOnset(smooth(vecD,5),vecT,[],[0 1],1);
			[dblOffset,dblValueOff,dblBaseVal,dblOffPeakT] = getOnset(smooth(vecD,5),vecT,[],[1 1.5],1);
			%z
			vecZ = zscore(vecD);
			intOnsetIdx = find(vecT>dblOnset,1);
			dblOnZ = vecZ(intOnsetIdx);
			if isempty(dblOnZ),dblOnZ = nan;end
			vecOnZ(intN) = dblOnZ;

			intOffsetIdx = find(vecT>dblOffset,1);
			dblOffZ = vecZ(intOffsetIdx);
			if isempty(dblOffZ),dblOffZ = nan;end
			vecOffZ(intN) = dblOffZ;

			%{
			%plot
			clf;
			plot(vecT,vecD)
			hold on
			scatter(dblOnset,dblValueOn,'o');
			scatter(dblOffset,dblValueOff,'x');
			hold off
			title(sprintf('%s, N%d; p=%.4f, Onset=%d, z=%.1f, Offset=%d, z=%.1f',strArea,intN,sLoad.vecP(intN),round(dblOnset*1000),dblOnZ,round(dblOffset*1000),dblOffZ))
			pause
			%}
			if abs(dblOnZ) > abs(dblOffZ)
				vecOnset(intN) = dblOnset;
			else
				vecOnset(intN) = dblOffset-1;
			end
		end
		%vecOnset(vecP>0.01) = 0.5;
		
		cellOnZ{intIdx} = vecOnZ;
		cellP{intIdx} = vecP;
		cellOnset{intIdx} = vecOnset;
			
	end
end

%% prep data
cellP2 = cellP;
cellOnset2 = cellOnset;
cellOnZ2 = cellOnZ;
cellDatasetNames2 = cellDatasetNames;
indRem = cellfun(@isempty,cellOnset2);
%indRem(7:end) = true;
cellOnset2(indRem) = [];
cellOnZ2(indRem) = [];
cellP2(indRem) = [];
cellDatasetNames2(indRem) = [];
cellOnset2 = cellfun(@(x) min(x,0.5),cellOnset2,'uniformoutput',false);
%cellOnset2 = cellfun(@(x,y) x(y<0.001),cellOnset2,cellP2,'uniformoutput',false);
%cellOnset3 = cellfun(@(x) min(x,0.1),cellOnset2,'uniformoutput',false);
vecRunSubset = [1:7];

vecSdOnset = cellfun(@std,cellOnset2);
vecOrderOnset = cellfun(@mean,cellOnset2(vecRunSubset))-vecSdOnset(vecRunSubset);
[d,vecReorder]=sort(vecOrderOnset,'ascend');
%cellOnset2 = cellOnset2(vecReorder);
%cellDatasetNames2 = cellDatasetNames2(vecReorder);

%% plot
matC=lines(numel(cellDatasetNames2));
figure
subplot(2,2,1);
hold on;
for intRec=1:numel(cellOnset2)
	stairs(1000*sort(cellOnset2{intRec}),linspace(0,1,numel(cellOnset2{intRec})),'Color',matC(intRec,:));
end
hold off;
legend(cellDatasetNames2,'location','best');
ylabel('Fraction of cells');
xlabel('Onset latency (ms)');
fixfig;

subplot(2,2,2);
hold on;
for intRec=vecRunSubset
	vecUseOnsets = cellOnset2{intRec};
	vecUseOnsets(vecUseOnsets>=0.5) = [];
	stairs(sort(vecUseOnsets)*1000,linspace(0,1,numel(vecUseOnsets)));
end
hold off;
ylabel('Fraction of cells');
xlabel('Onset latency (ms)');
fixfig;


subplot(2,2,3);
hold on
vecAllOnsets=[];
cellAreas = {};
intC=0;
cellAreaPlot={};
for intArea=vecRunSubset
	intC=intC+1;
	cellAreaPlot{intC} = cellDatasetNames2{intArea};
	
	vecUseOnsets = cellOnset2{intArea};
	vecUseOnsets(vecUseOnsets>=0.5) = [];
	
	vecAllOnsets((end+1):(end+numel(vecUseOnsets))) = vecUseOnsets;
	cellAreas((end+1):(end+numel(vecUseOnsets))) = cellfill(cellDatasetNames2{intArea},[1 numel(vecUseOnsets)]);
	errorbar(intC,1000*mean(vecUseOnsets),1000*std(vecUseOnsets),'x','Color',matC(intArea,:))
	
end

set(gca,'xtick',1:numel(vecRunSubset),'xticklabel',cellAreaPlot);
xlim([0 numel(vecRunSubset)+1]);
%ylim([0 120])
ylabel('Onset latency (ms)');
fixfig


subplot(2,2,4);
hold on
vecAllOnsets=[];
cellAreas = {};
intC=0;
cellAreaPlot={};
for intArea=vecRunSubset
	intC=intC+1;
	cellAreaPlot{intC} = cellDatasetNames2{intArea};
	
	vecUseOnsets = cellOnset2{intArea};
	vecUseOnsets(vecUseOnsets>=0.5) = [];
	
	vecAllOnsets((end+1):(end+numel(vecUseOnsets))) = vecUseOnsets;
	cellAreas((end+1):(end+numel(vecUseOnsets))) = cellfill(cellDatasetNames2{intArea},[1 numel(vecUseOnsets)]);
	bplot(vecUseOnsets*1000,intC);
	
end

set(gca,'xtick',1:numel(vecRunSubset),'xticklabel',cellAreaPlot);
xlim([0 numel(vecRunSubset)+1]);
%ylim([0 120])
ylabel('Onset latency (ms)');
fixfig
maxfig;
return
%% save
drawnow;
export_fig(sprintf('%sMetaOnsetsFig.tif',strFigPath));
export_fig(sprintf('%sMetaOnsetsFig.pdf',strFigPath));
print(gcf,'-dpdf', sprintf('%sMetaOnsetsFig.pdf',strFigPath));