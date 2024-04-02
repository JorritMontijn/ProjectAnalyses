%% set variables
if strcmp(strRunStim,'DG')
	strLoadStim = 'driftinggrating';
	strLoadFile= 'Simulation_xAreaDistributed_SG18_2019-07-04.mat';
else
	strLoadFile='Simulation_xAreaDistributed_IndRetNoise0_0_2019-08-28.mat';
	strLoadStim = '';
end

%get prepped sim files
cellUseAreas = {...
	'SimV1',...
	...'posteromedial visual area',...
	};
strArea = cellUseAreas{1};
strAreaAbbr = 'V1';

sSimRecs = dir(fullpath(strDataPathSim,'SimDG18_MatchedTo*.mat'));
intRecNum = numel(sSimRecs);