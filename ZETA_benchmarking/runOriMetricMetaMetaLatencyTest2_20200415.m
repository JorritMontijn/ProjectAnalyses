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
dblBinSize = 0.001;
vecBinEdges = 0:dblBinSize:0.5;
vecBins = vecBinEdges(2:end) - dblBinSize/2;
vecVals = zeros([1 numel(vecBinEdges)-1]);
matAllAct = nan([numel(cellUniqueAreas) numel(vecVals)]);
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
	intStimType=vecRunStims;
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
	intC = 0;
	for intN=1:intNumN
		vecT = sLoad.cellInterpT{intN};
		vecD = sLoad.cellDeriv{intN};
		vecP(intN) = sLoad.vecP(intN);
		if isempty(vecT) || (vecP(intN) > 0.01 && intRandIdx == 1),
			continue;
		end
		intC = intC + 1;
		[a,b]=unique(vecT);
		vecNewVals = interp1(vecT(b),vecD(b),vecBins);
		vecVals = vecVals + vecNewVals(:)';
	end
	vecVals = vecVals / intC;
	matAllAct(intIdx,:) = vecVals;
	
	%% onset
	[dblOnset,dblValueOn,dblBaseVal,dblPeakT] = getOnset(vecVals,vecBins,[],[0 1],0);
	%z
	vecZ = zscore(vecVals);
	intOnsetIdx = find(vecBins>dblOnset,1);
	dblOnZ = vecZ(intOnsetIdx);
	if isempty(dblOnZ),dblOnZ = nan;end
	vecOnZ(intIdx) = dblOnZ;
	
	%plot
	clf;
	plot(vecBins,vecVals)
	hold on
	scatter(dblOnset,dblValueOn,'x');
	hold off
	title(sprintf('%s, N=%d; Onset=%d, z=%.1f',strArea,intC,round(dblOnset*1000),dblOnZ))
	%pause
	
	%save
	vecOnset(intIdx) = dblOnset;
end

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