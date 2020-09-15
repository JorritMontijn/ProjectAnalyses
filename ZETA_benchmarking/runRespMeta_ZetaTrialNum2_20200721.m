clear all;
%close all;
strDisk = 'F:';
strAnalysisType = '1A'; %A=500ms, B=300ms; 1=pooled,2=split,3=split,uncorr
strZetaType = 'v2'; %{'','v2'}
strDataSource = [strDisk '\Data\Processed\ZETA\TrialNum\'];
strFigPath = [strDisk '\Data\Results\ZETA\Inclusion\'];
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
cellNumCells = [];
cellNumCells2 = [];
cellTrialNums2 = [];
			
intIdx = 0;
for intArea=8%7:numel(cellUniqueAreas)
	strArea = cellUniqueAreas{intArea}; %V1, SC, Retina, Poisson, GCaMP
	strName = replace(strArea,cellRepStr(:,1),cellRepStr(:,2));
	intIdx = intIdx + 1;
	for intRandType=1:2
		%set var
		%% load data 1
		strRand = cellRunRand{intRandType};
		strRunType = [strArea strRand];
		sDir=dir([strDataSource 'ZetaData' strZetaType 'TrialNum' strRunType '*']);
		if ~isempty(sDir) && isempty(strRand)
			cellDatasetNames{intIdx} = strName;
			sDir = sDir(~contains({sDir.name},'Rand'));
		end
		%load data
		intFiles=numel(sDir);
		for intFile=1:intFiles
			strFile = sDir(intFile).name;
			strPath = sDir(intFile).folder;
			sLoad=load([strPath filesep strFile]);
			
			matHzP=sLoad.matHzP;
			matZeta=sLoad.matZeta;
			matZetaP=sLoad.matZetaP;
			vecSubsampleTrials = sLoad.vecSubsampleTrials;
			intN = size(matZetaP,1);
			cellNumCells{intIdx,intRandType} = sum(~isnan(matZeta),1);
			cellTrialNums{intIdx,intRandType} = vecSubsampleTrials;
			cellIncludeZ{intIdx,intRandType} = sum(matZetaP<0.05,1)/intN;
			cellIncludeM{intIdx,intRandType} = sum(matHzP<0.05,1)/intN;
			
		end
		%% load data 2
		strRand = cellRunRand{intRandType};
		strRunType = [strrep(strAnalysisType,'3','2') '_' strArea strRand];
		sDir=dir([strDataSource 'ZetaData' strZetaType 'TrialNum*' strRunType '*']);
		if ~isempty(sDir) && isempty(strRand)
			cellDatasetNames{intIdx} = strName;
			sDir = sDir(~contains({sDir.name},'Rand'));
		end
		%load data
		intFiles=numel(sDir);
		for intFile=1:intFiles
			strFile = sDir(intFile).name;
			strPath = sDir(intFile).folder;
			sLoad=load([strPath filesep strFile]);
			
			matCorrP=sLoad.matCorrP;
			matAllP = sLoad.matAllP;
			matMinP = min(matAllP,[],3);
			vecSubsampleTrials = sLoad.vecSubsampleTrials;
			intN = size(matCorrP,1);
			cellNumCells2{intIdx,intRandType} = sum(~isnan(matZeta),1);
			cellTrialNums2{intIdx,intRandType} = vecSubsampleTrials;
			%cellIncludeM{intIdx,intRandType} = sum(matMinP<0.05,1)/intN;
            if contains(strAnalysisType,'1') || contains(strAnalysisType,'2') %corr
                cellIncludeM{intIdx,intRandType} = sum(matCorrP<0.05,1)/intN;
            elseif contains(strAnalysisType,'3') %remove corr
                cellIncludeM{intIdx,intRandType} = sum((matCorrP/24)<0.05,1)/intN;
            else
               error 
            end
		end
	end
end

%%
if isempty(cellNumCells2),cellNumCells2=cellNumCells;end
indRem = any(cellfun(@max,cellNumCells) < 20,2);
indRem2 = any(cellfun(@max,cellNumCells2) < 20,2);
indRem2(numel(indRem)+1:end) = true;
indRem(numel(indRem2)+1:end) = true;

cellTrialNums(indRem,:) = [];
cellNumCells(indRem,:) = [];
cellIncludeZ(indRem,:) = [];
cellIncludeM(indRem2,:) = [];
cellDatasetNames(indRem2) = [];



for intDataset=1:size(cellIncludeZ,1)
	intMaxN = cellNumCells{intDataset,1}(1);
	vecUseData = cellNumCells{intDataset,1} == intMaxN;
	[vecP_ZI,matCI_ZI] = binofit(cellIncludeZ{intDataset,1}(vecUseData)*intMaxN,cellNumCells{intDataset,1}(vecUseData),0.25);
	[vecP_ZFA,matCI_ZFA] = binofit(cellIncludeZ{intDataset,2}(vecUseData)*intMaxN,cellNumCells{intDataset,2}(vecUseData),0.25);
	
	[vecP_HzI,matCI_HzI] = binofit(cellIncludeM{intDataset,1}(vecUseData)*intMaxN,cellNumCells{intDataset,1}(vecUseData),0.25);
	[vecP_HzFA,matCI_HzFA] = binofit(cellIncludeM{intDataset,2}(vecUseData)*intMaxN,cellNumCells{intDataset,2}(vecUseData),0.25);
	
[vecP_ZvsM,vecZ]=bino2test(cellIncludeZ{intDataset,1}(vecUseData).*intMaxN,cellNumCells{intDataset,1}(vecUseData),cellIncludeM{intDataset,1}(vecUseData).*intMaxN,cellNumCells{intDataset,1}(vecUseData));
[dummy,dummy,vecP_ZvsM_corr] = fdr_bh(vecP_ZvsM);

figure
	plot([0 cellTrialNums{intDataset,1}(vecUseData)],[0 0 0; vecP_ZI' matCI_ZI],'-b')
	hold on
	plot([0 cellTrialNums{intDataset,2}(vecUseData)],[0 0 0; vecP_ZFA' matCI_ZFA],'-m')
	plot([0 cellTrialNums{intDataset,1}(vecUseData)],[0 0 0; vecP_HzI' matCI_HzI],'-k')
	plot([0 cellTrialNums{intDataset,2}(vecUseData)],[0 0 0; vecP_HzFA' matCI_HzFA],'-r')
	hold off
	legend({'Z','Z FA','M','M FA'},'location','best');
	xlabel('# of stimulus repetitions');
	ylabel('Fraction of significant cells')
	%title(sprintf('Incl p=%.3f, FA p=%.3f',pI,pF));
	title(sprintf('%s; %s, N=%d-%d',cellDatasetNames{intDataset},strAnalysisType,min(cellNumCells{intDataset,1}),max(cellNumCells{intDataset,1})));
	fixfig;
end
%%
strTag = [strAnalysisType '_' cellDatasetNames{intDataset}];
drawnow;
export_fig(sprintf('%sMeta%s.tif',strFigPath,strTag));
export_fig(sprintf('%sMeta%sEF.pdf',strFigPath,strTag));
print(gcf,'-dpdf', sprintf('%sMeta%s.pdf',strFigPath,strTag));
return
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
