%{'VISal'}    {'VISam'}    {'VISl'}    {'VISp'}    {'VISpm'}    {'VISrl'}
strRunArea = 'VISp';%'posteromedial visual area' 'Primary visual area'
cellUseAreas = {strRunArea};
strArea = strRunArea;
strAreaAbbr = 'VISp';	

%% load ABI
if ~exist('sAggABI','var') || isempty(sAggABI)
	fprintf('Loading Allen Brain EcEphys data... [%s]\n',getTime);
	sLoad = load(fullpath(strDataPathABI,'AggSes2022-03-28.mat'));
	sAggABI = sLoad.sSes;
	clear sLoad;
end
vecUseRec = 1:numel(sAggABI);%find(contains({sAggABA.Exp},'MP'));
intRecNum = numel(sAggABI);