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
	'Primary visual'...Area 8
	'Lateral posterior nucleus'...Area 9
	'Field CA3'...Area 10
	'Field CA1'...Area 11
	'Dentate Gyrus'...Area 12
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
	'lateral geniculate','LGN';...
	'Primary visual','V1';...
	'RunDriftingGratings','-DG';...
	'RunNaturalMovie','-NM';...
	'Lateral posterior nucleus','LP';...
	'Field CA3','Hip-CA3';...
	'Field CA1','Hip-CA1';...
	'Dentate Gyrus','DentGyr';...
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
		strFile = sDir(1).name;
		intResampNum = str2double(getFlankedBy(strFile,'Resamp','.mat'));
		sLoad=load([strPath strFile]);
		
		
		intNumN = numel(sLoad.cellDeriv);
		
		vecOnset = nan(1,intNumN);
		vecP = nan(1,intNumN);
		
		for intN=1:intNumN
			vecT = sLoad.cellInterpT{intN};
			vecD = sLoad.cellDeriv{intN};
			
			[dblPeakRate,dblPeakTime,dblPeakWidth,vecPeakStartStop,intPeakLoc,vecPeakStartStopIdx] = ...
				getPeak(vecD,vecT,[0 inf]);
			[dblOnset,dblValue] = getOnset(vecD,vecT,dblPeakTime,[0 inf]);
			
			vecOnset(intN) = dblOnset;
			vecP(intN) = sLoad.vecP(intN);
		end
		
		cellP{intIdx} = vecP;
		vecOnset(vecP>0.01) = 0.5;
		cellOnset{intIdx} = vecOnset;
			
	end
end

%% prep data
cellOnset2 = cellfun(@(x) min(x,0.5),cellOnset,'uniformoutput',false);

vecAllOnsets = cell2vec(cellOnset2);

%% plot
figure
hold on;
for intRec=1:numel(cellOnset)
	stairs(sort(cellOnset2{intRec}),linspace(0,1,numel(cellOnset2{intRec})));
end
hold off;
legend(cellDatasetNames,'location','best');
fixfig;
return
%% save
drawnow;
export_fig(sprintf('%sMetaOnsetsFig.tif',strFigPath));
export_fig(sprintf('%sMetaOnsetsFig.pdf',strFigPath));
print(gcf,'-dpdf', sprintf('%sMetaOnsetsFig.pdf',strFigPath));