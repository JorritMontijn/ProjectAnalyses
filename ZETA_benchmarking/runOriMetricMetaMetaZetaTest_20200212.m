clear all;
%close all;
strPath = 'D:\Data\Results\OriMetric\Data\';
strFigPath = 'D:\Data\Results\OriMetric\';
cellUniqueAreas = {...
	'V1',...Area 1
	'SC',...Area 2
	'Poisson',...Area 3
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
	'posteromedial visual'};%,...Area 14
%{
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
%}
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
matSignifHz = [];
intIdx = 0;
for intArea=1:numel(cellUniqueAreas)
	strArea = cellUniqueAreas{intArea}; %V1, SC, Retina, Poisson, GCaMP
	if intArea < 7
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
		sDir=dir([strPath 'ZetaDataMSD' strRunType 'Resamp100*']);
		intFiles=numel(sDir);
		for intFile=1:intFiles
			strFile = sDir(intFile).name;
			intResampNum = str2double(getFlankedBy(strFile,'Resamp','.mat'));
			sLoad=load([strPath strFile]);
			vecZeta = abs(sLoad.vecZeta);
			vecZP=1-(normcdf(abs(vecZeta))-normcdf(-abs(vecZeta)));
			matNumCells(intIdx,intRandType) = numel(vecZP);
			matSignifZ(intIdx,intRandType) = sum(vecZP<0.05);
			matSignifHz(intIdx,intRandType) = sum(sLoad.vecHzP<0.05);
		end
	end
	end
end
%%
indRem = any(matNumCells < 20,2);
matSignifZ(indRem,:) = [];
matNumCells(indRem,:) = [];
matSignifHz(indRem,:) = [];
cellDatasetNames(indRem) = [];

[vecP_ZI,matCI_ZI] = binofit(matSignifZ(:,1),matNumCells(:,1),0.25);
[vecP_ZFA,matCI_ZFA] = binofit(matSignifZ(:,2),matNumCells(:,2),0.25);

[vecP_HzI,matCI_HzI] = binofit(matSignifHz(:,1),matNumCells(:,1),0.25);
[vecP_HzFA,matCI_HzFA] = binofit(matSignifHz(:,2),matNumCells(:,2),0.25);
figure
for intDataset=1:size(matSignifZ,1)
	matTable = zeros(2,2);
	matTable(1,1) = matSignifZ(intDataset,1);
	matTable(1,2) = matNumCells(intDataset,1)-matSignifZ(intDataset,1);
	matTable(2,1) = matSignifHz(intDataset,1);
	matTable(2,2) = matNumCells(intDataset,1)-matSignifHz(intDataset,1);
	[hI,pI] = fishertest(matTable);
	
	matTable = zeros(2,2);
	matTable(1,1) = matSignifZ(intDataset,2);
	matTable(1,2) = matNumCells(intDataset,2)-matSignifZ(intDataset,2);
	matTable(2,1) = matSignifHz(intDataset,2);
	matTable(2,2) = matNumCells(intDataset,2)-matSignifHz(intDataset,2);
	[hF,pF] = fishertest(matTable);
	
	subplot(4,6,intDataset)
	errorbar(2,vecP_ZI(intDataset),matCI_ZI(intDataset,1)-vecP_ZI(intDataset),matCI_ZI(intDataset,2)-vecP_ZI(intDataset),'bx')
	hold on
	errorbar(4,vecP_ZFA(intDataset),matCI_ZFA(intDataset,1)-vecP_ZFA(intDataset),matCI_ZFA(intDataset,2)-vecP_ZFA(intDataset),'bx')
	errorbar(1,vecP_HzI(intDataset),matCI_HzI(intDataset,1)-vecP_HzI(intDataset),matCI_HzI(intDataset,2)-vecP_HzI(intDataset),'kx')
	errorbar(3,vecP_HzFA(intDataset),matCI_HzFA(intDataset,1)-vecP_HzFA(intDataset),matCI_HzFA(intDataset,2)-vecP_HzFA(intDataset),'kx')
	hold off
	set(gca,'xtick',1:4,'xticklabel',{getGreek('Mu'),getGreek('Zeta'),['FA ' getGreek('Mu')],['FA ' getGreek('Zeta')]})
	%title(sprintf('Incl p=%.3f, FA p=%.3f',pI,pF));
	title(sprintf('%s; N=%d',cellDatasetNames{intDataset},matNumCells(intDataset,1)));
	ylabel('Fraction of significant cells')
	xlim([0.5 4.5])
	ylim([0 1])
	fixfig;
end

[h,pI]=ttest(vecP_ZI,vecP_HzI);
[h,pF]=ttest(vecP_ZFA,vecP_HzFA);
subplot(4,6,24)
plot(repmat([1 2],[intDataset 1])',[vecP_HzI vecP_ZI]','g')
hold on
plot(repmat([3 4],[intDataset 1])',[vecP_HzFA vecP_ZFA]','r')
title(sprintf('Incl p=%.3f, FA p=%.3f',pI,pF));
hold off	
xlim([0.5 4.5])
ylim([0 1])
set(gca,'xtick',1:4,'xticklabel',{getGreek('Mu'),getGreek('Zeta'),['FA ' getGreek('Mu')],['FA' getGreek('Zeta')]})
ylabel('Fraction of significant cells')
fixfig;
maxfig()
return
drawnow;
export_fig(sprintf('%sMetaSummaryFig.tif',strFigPath));
export_fig(sprintf('%sMetaSummaryFigEF.pdf',strFigPath));
print(gcf,'-dpdf', sprintf('%sMetaSummaryFig.pdf',strFigPath));