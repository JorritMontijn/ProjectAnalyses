%% load data
if exist('E:\DataPreProcessed','dir')
	strDataPath = 'E:\DataPreProcessed';
	%strDataPath = 'C:\Drive\Neuropixels\Old';
else
	strDataPath = 'F:\Data\Processed\Neuropixels';
	%strDataPath = 'F:\Drive\Neuropixels\Old';
end
fprintf('Loading %s [%s]\n',strDataPath,getTime);
sFiles = dir(fullpath(strDataPath,'*_AP.mat'));
sFiles(contains({sFiles.name},'MP')) = [];
if ~exist('sExp','var') || isempty(sExp)
	sExp = [];
	for intFile=1:numel(sFiles)
		fprintf('Loading %d/%d: %s [%s]\n',intFile,numel(sFiles),sFiles(intFile).name,getTime);
		sLoad = load(fullpath(sFiles(intFile).folder,sFiles(intFile).name));
		sLoad.sAP.Name = sFiles(intFile).name;
		if ~isfield(sLoad.sAP,'sPupil')
			sLoad.sAP.sPupil = [];
		end
		if isempty(sExp)
			sExp = sLoad.sAP;
		else
			sExp(end+1) = sLoad.sAP;
		end
	end
end

%MP_20200115 eye tracking remove last stimulus (gunk in eye)
cellUseForEyeTrackingMP = {'Topo'}; %don't forget to set high vid lum as blinks
cellUseForEyeTrackingMA = {'20210212','20210215','20210218','20210220','20210225','20210301'};
cellUseForEyeTracking = cat(2,cellUseForEyeTrackingMA,cellUseForEyeTrackingMP);
strTargetPath = 'D:\Data\Results\AlbinoProject';

%best rec BL6: 20191216B5 (rec 17)
%best rec DBA: 20210212B2 (rec 5)

%% define area categories
%atlas location
if exist('E:\AllenCCF','dir')
	strAllenCCFPath = 'E:\AllenCCF';
else
	strAllenCCFPath = 'F:\Data\AllenCCF';
end

%cortex
%cellUseAreas{1} = {'Primary visual','Posteromedial visual','anteromedial visual'};
cellUseAreas{1} = {'Primary visual','Posteromedial visual'};
%NOT
cellUseAreas{2} = {'nucleus of the optic tract'};
%hippocampus
cellUseAreas{3} = {'Field CA1','Field CA2','Field CA3'};
%cellUseAreas{3} = {'Hippocampal','Field CA1','Field CA2','Field CA3','subiculum','Postsubiculum','Prosubiculum','dentate gyrus'};
%cellUseAreas{3} = {'superior colliculus'};
cellAreaGroups = {'Vis. ctx','NOT','Hippocampus'};
cellAreaGroupsAbbr = {'Ctx','NOT','Hip'};
cellSubjectGroups = {'BL6','DBA'};