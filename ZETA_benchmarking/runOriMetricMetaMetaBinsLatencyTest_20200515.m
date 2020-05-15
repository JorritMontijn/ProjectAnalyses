clear all;
%close all;
strPath = 'F:\Data\Processed\ZETA\Latencies\';
strFigPath = 'F:\Data\Results\ZETA\Latencies\';
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
	'Anterior pretectal nucleus',...Area 10
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
	...'Field CA2',...Area 21
	...'Field CA3',...Area 22
	...'Dentate gyrus',...Area 23
	...'Retrosplenial'...Area 24
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
intRandIdx=1; %1=normal,2=rand
intIdx = 0;
for intArea=8:numel(cellUniqueAreas)
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
	intStimType=vecRunStims;
	intIdx = intIdx + 1;
	strRand = cellRunRand{intRandIdx};
	strName = replace(strArea,cellRepStr(:,1),cellRepStr(:,2));
	cellDatasetNames{intIdx} = strName;
	
	%% load data
	strRunType = [strArea strRand];
	sDir=dir([strPath 'ZetaDataBinsLatencies2' strRunType '*']);
	if isempty(sDir),continue;end
	strFile = sDir(1).name;
	sLoad=load([strPath strFile]);
	
	%%
	for intNeuron=1:numel(sLoad.vecZetaLatencies)
		clf;
		hold on
		scatter(sLoad.vecBinDurs,sLoad.matBinLatencies(:,intNeuron),'k');
		plot(sLoad.vecBinDurs,ones(size(sLoad.vecBinDurs))*sLoad.vecZetaLatencies(intNeuron),'b');
		hold off
		title(sprintf('%s, N%d, Onset=%.3f',strArea,intNeuron,sLoad.vecZetaLatencies(intNeuron)))
		set(gca,'xscale','log')
		pause
	end
	return
end
%% list
%[6 9]
%% prep data
vecRunSubset = [1:8];
matAllAct2 = matAllAct;
cellDatasetNames2 = cellDatasetNames;
indRem = sum(matAllAct,2)==0;
matAllAct2(indRem,:) = [];
matAllAct2 = matAllAct2(vecRunSubset,:);
vecOnset2 =vecOnset;

matNormCumSum = cumsum(matAllAct-matAllAct(:,1),2);
matNormCumSum2 = matNormCumSum./sum(matNormCumSum,2);

matNormCumSumB = cumsum(matAllAct(:,1:100)-matAllAct(:,1),2);
matNormCumSumB2 = matNormCumSumB./sum(matNormCumSumB,2);

matUseData = matNormCumSumB2;
vecUseBins = vecBins(1:100);
%% plot
matC=jet(numel(vecRunSubset));
figure
subplot(2,2,1);
hold on;
for intRec=vecRunSubset
	[dblOnset,dblValueOn,dblBaseVal,dblPeakT] = getOnset(matUseData(intRec,:),vecUseBins,[],[0 1],0);
	vecOnsets(intRec) = dblOnset;
	stairs(vecUseBins,matUseData(intRec,:),'Color',matC(intRec,:));
end
hold off;
legend(cellDatasetNames2,'location','best');
ylabel('Normalized activity');
xlabel('Time (ms)');
fixfig;

%reorder
[vecOnset3,vecReorder] = sort(vecOnsets(vecRunSubset));
subplot(2,2,2);
scatter(1:numel(vecOnset3),vecOnset3);
set(gca,'xtick',1:numel(vecOnset3),'xticklabel',cellDatasetNames2(vecReorder));

return
%% save
drawnow;
export_fig(sprintf('%sMetaOnsetsFig%.3f-%.3f.tif',strFigPath,dblFloorThresh,dblCeilThresh));
export_fig(sprintf('%sMetaOnsetsFig%.3f-%.3f.pdf',strFigPath,dblFloorThresh,dblCeilThresh));
print(gcf,'-dpdf', sprintf('%sMetaOnsetsFig%.3f-%.3f.pdf',strFigPath,dblFloorThresh,dblCeilThresh));
%}