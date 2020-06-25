clear all;
%close all;
strDisk = 'F:';
strDataSource = [strDisk '\Data\Processed\ZETA\TrialNum\'];
strFigPath = [strDisk '\Data\Results\ZETA\TrialNum\'];
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
for intArea=7:numel(cellUniqueAreas)
	strArea = cellUniqueAreas{intArea}; %V1, SC, Retina, Poisson, GCaMP
	strName = replace(strArea,cellRepStr(:,1),cellRepStr(:,2));
	intIdx = intIdx + 1;
	for intRandType=1:2
		%set var
		%% load data
		strRand = cellRunRand{intRandType};
		strRunType = [strArea strRand];
		sDir=dir([strDataSource 'ZetaDataTrialNum*' strRunType '*']);
		if ~isempty(sDir) && isempty(strRand)
			cellDatasetNames{intIdx} = strName;
			sDir = sDir(~contains({sDir.name},'Rand'));
		end
		%load data
		intFiles=numel(sDir);
		for intFile=1:intFiles
			strFile = sDir(intFile).name;
			strPath = sDir(intFile).folder;
			intResampNum = str2double(getFlankedBy(strFile,'Resamp','.mat'));
			sLoad=load([strPath filesep strFile]);
			
			matZeta=sLoad.matZeta;
			matZetaP=sLoad.matZetaP;
			matHzP=sLoad.matHzP;
			matHzD=sLoad.matHzD;
			vecSubsampleTrials = sLoad.vecSubsampleTrials;
			intN = size(matZetaP,1);
			matNumCells(intIdx,intRandType) = intN;
			cellTrialNums{intIdx,intRandType} = vecSubsampleTrials;
			cellIncludeZ{intIdx,intRandType} = sum(matZetaP<0.05,1)/intN;
			cellIncludeM{intIdx,intRandType} = sum(matHzP<0.05,1)/intN;
		end
		
	end
end

%%
indRem = any(matNumCells < 20,2);
cellTrialNums(indRem,:) = [];
matNumCells(indRem,:) = [];
cellIncludeZ(indRem,:) = [];
cellIncludeM(indRem,:) = [];
cellDatasetNames(indRem) = [];


figure
for intDataset=1:size(matSignifZ,1)
	[vecP_ZI,matCI_ZI] = binofit(cellIncludeZ{intDataset,1}*matNumCells(intDataset,1),matNumCells(intDataset,1),0.25);
	[vecP_ZFA,matCI_ZFA] = binofit(cellIncludeZ{intDataset,2}*matNumCells(intDataset,2),matNumCells(intDataset,2),0.25);
	
	[vecP_HzI,matCI_HzI] = binofit(cellIncludeM{intDataset,1}*matNumCells(intDataset,1),matNumCells(intDataset,1),0.25);
	[vecP_HzFA,matCI_HzFA] = binofit(cellIncludeM{intDataset,2}*matNumCells(intDataset,2),matNumCells(intDataset,2),0.25);
	


	errorbar(cellTrialNums{intDataset,1}/24,vecP_ZI,matCI_ZI(:,1)'-vecP_ZI,matCI_ZI(:,2)'-vecP_ZI,'-bx')
	hold on
	errorbar(cellTrialNums{intDataset,2}/24,vecP_ZFA,matCI_ZFA(:,1)'-vecP_ZFA,matCI_ZFA(:,2)'-vecP_ZFA,'-bo')
	errorbar(cellTrialNums{intDataset,1}/24-1,vecP_HzI,matCI_HzI(:,1)'-vecP_HzI,matCI_HzI(:,2)'-vecP_HzI,'-kx')
	errorbar(cellTrialNums{intDataset,2}/24+1,vecP_HzFA,matCI_HzFA(:,1)'-vecP_HzFA,matCI_HzFA(:,2)'-vecP_HzFA,'-ko')
	hold off
	legend({'Z','Z FA','M','M FA'},'location','best');
	xlabel('# of stimulus repetitions');
	ylabel('Fraction of significant cells')
	%title(sprintf('Incl p=%.3f, FA p=%.3f',pI,pF));
	title(sprintf('%s; N=%d',cellDatasetNames{intDataset},matNumCells(intDataset,1)));
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