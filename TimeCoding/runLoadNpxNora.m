%% set variables
cellUseAreas = {...
	'Primary somatosensory area',...
	...'posteromedial visual area',...
	};
strArea = cellUseAreas{1};
strAreaAbbr = 'S1';
if strcmp(strRunStim,'WS')
	strLoadStim = 'RunOptoWhiskerStim';
end

%% load data
if ~exist('sAggStim','var') || isempty(sAggStim)
	[sAggStim,sAggNeuron,sAggSources]=loadDataNpx(strArea,strLoadStim,strDataPathNora);
end
%indKeepCells = [sAggNeuron.KilosortGood] | [sAggNeuron.Contamination] < 0.1;
%sAggNeuron = sAggNeuron(indKeepCells);
sAggNeuron(strcmpi({sAggNeuron.SubjectType},'DBA')) = [];
intRecNum = numel(sAggStim);