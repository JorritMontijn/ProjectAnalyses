%% set variables
cellUseAreas = {...
	'Primary visual area',...
	...'posteromedial visual area',...
	};
strArea = cellUseAreas{1};
strAreaAbbr = 'V1';
if strcmp(strRunStim,'DG')
	strLoadStim = 'driftinggrating';
else
	strLoadStim = 'naturalmovie';
end
%% load data
if ~exist('sAggStim','var') || isempty(sAggStim)
	[sAggStim,sAggNeuron,sAggSources]=loadDataNpx(strArea,strLoadStim,strDataPath);
end
sAggNeuron(strcmpi({sAggNeuron.SubjectType},'DBA')) = [];
intRecNum = numel(sAggStim);