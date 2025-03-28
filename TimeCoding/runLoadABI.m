%{'VISal'}    {'VISam'}    {'VISl'}    {'VISp'}    {'VISpm'}    {'VISrl'}
strRunArea = 'VISp';%'posteromedial visual area' 'Primary visual area'
cellUseAreas = {strRunArea};
strArea = strRunArea;
strAreaAbbr = 'VISp';

%% check if necessary
if exist('sAggABI','var') && ~isempty(sAggABI)
	fprintf('ABI data are already loaded\n');
	return;
end
%% load ABI v1
strFileV1 = 'AggSes2022-03-28.mat';
if exist(fullpath(strDataPathABI,strFileV1),'file')
	fprintf('Loading Allen Brain EcEphys data... [%s]\n',getTime);
		sLoad = load(fullpath(strDataPathABI,'AggSes2022-03-28.mat'));
		sAggABI = sLoad.sSes;
		clear sLoad;
	vecUseRec = 1:numel(sAggABI);%find(contains({sAggABA.Exp},'MP'));
	intRecNum = numel(sAggABI);
else
	%% v2, load all sessions separately
	fprintf('Loading Allen Brain Visual Behaviour data... [%s]\n',getTime);
	sFiles = dir(fullpath(strDataPathABI,'AggSes_ecephys_session*.mat'));
	clear sAggABI;
	for i=1:numel(sFiles)
		strFile = fullpath(sFiles(i).folder,sFiles(i).name);
		sLoad = load(strFile);
		sAggABI(i) = sLoad.sSes;
		fprintf('Loaded %d/%d (%s) [%s]\n',i,numel(sFiles),sFiles(i).name,getTime);
	end
end