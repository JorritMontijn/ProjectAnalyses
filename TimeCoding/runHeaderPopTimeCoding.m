%clear data
clearvars -except cellDataTypes intRunDataType strRunType cellTypes strRunStim boolFixSpikeGroupSize sAggABI dblRemOnset;

%set paths
if isfolder('F:\Drive\PopTimeCoding') && isfolder('F:\Data\Processed\Neuropixels\')
	strDataPathABI = '';
	strDataPathSim = 'F:\Data\Processed\PopTimeCoding\';
	strDataPathSimT0 = 'F:\Drive\PopTimeCoding\data\';
	strDataPath = 'F:\Data\Processed\Neuropixels\';
	strDataPathNora = 'F:\Data\Processed\PopTimeCoding\Noradata';
	strFigurePathSR = 'F:\Drive\PopTimeCoding\single_recs';
	strFigurePath = 'F:\Drive\PopTimeCoding\figures\';
	strTargetDataPath = 'F:\Drive\PopTimeCoding\data\';
else
	strDataPathABI = 'E:\AllenBrainVisualEphys';
	strDataPathSim = 'D:\Data\Processed\PopTimeCoding\';
	strDataPathSimT0 = 'C:\Drive\PopTimeCoding\data\';
	strDataPath = 'E:\DataPreProcessed\';
	strDataPathNora = 'D:\Data\Processed\PopTimeCoding\Noradata';
	strFigurePathSR = 'C:\Drive\PopTimeCoding\single_recs';
	strFigurePath = 'C:\Drive\PopTimeCoding\figures\';
	strTargetDataPath = 'C:\Drive\PopTimeCoding\data\';
end
%% define data
if ~exist('cellDataTypes','var')
	cellDataTypes = {'Npx','Sim','ABI','SWN'};%topo, model, allen, nora
	intRunDataType = 1;
end
if exist('intRunDataType','var') && intRunDataType == 4
	%whisking nora s1
	cellUseAreas = {...
		'Primary somatosensory area',...
		...'posteromedial visual area',...
		};
	strArea = cellUseAreas{1};
	strRunType = 'SWN'; %somatosensory whisking nora
	strRunStim = 'WS';%WS=whisker stimulation
else
	cellUseAreas = {...
		'Primary visual area',...
		...'posteromedial visual area',...
		};
	strArea = cellUseAreas{1};
end
strRunType = cellDataTypes{intRunDataType};

%% pre-allocate matrices
if exist('cellTypes','var')
	intNumTypes = numel(cellTypes);
end
if ~exist('dblRemOnset','var')
	dblRemOnset = 0;
end
if strcmp(strRunType,'ABI')
	runLoadABI;
elseif strcmp(strRunType,'Sim')
	%get prepped sim files
	sSimRecs = dir(fullpath(strDataPathSim,'SimDG18_MatchedTo*.mat'));
	intRecNum = numel(sSimRecs);
elseif strcmp(strRunType,'SWN')
	runLoadNpxNora;
elseif strcmp(strRunType,'Npx')
	runLoadNpx;
else
	error unknown
end

%% onset string
if dblRemOnset == 0
	strOnset = '';
else
	strOnset = sprintf('%.2f',dblRemOnset);
end