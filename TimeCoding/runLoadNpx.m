%% set variables
cellUseAreas = {...
	'Primary visual area',...
	...'posteromedial visual area',...
	};
strArea = cellUseAreas{1};
strAreaAbbr = 'V1';	

%% load data
if ~exist('sAggStim','var') || isempty(sAggStim)
	[sAggStim,sAggNeuron,sAggSources]=loadDataNpx('','driftinggrating',strDataPath);
end
sAggNeuron(strcmpi({sAggNeuron.SubjectType},'DBA')) = [];
intRecNum = numel(sAggStim);