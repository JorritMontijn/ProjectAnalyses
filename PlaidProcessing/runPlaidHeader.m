
%% load data
if (~exist('sData','var') && ~exist('sSimRun','var')) || ~exist('strLastSim','var') || ~strcmp(strLastSim,strSimulation)
	%msg
	fprintf('Loading plaid data set <%s> (G:xyt%02d/P:xyt%02d) [%s]\n',strDate,intBlockGratings,intBlockPlaids,getTime);
	
	%load data
	sLoadGratings = load(strGratingFile);
	sesG = sLoadGratings.ses;
	clear sLoadGratings;
	sLoadPlaids = load(strPlaidFile);
	sesP = sLoadPlaids.ses;
	clear sLoadPlaids;
end

%% concatenate
%add contrast field to gratings
sesG.structStim.Contrast = 100*ones(size(sesG.structStim.Orientation));
sSesAggregate = buildMultiSesAggregate(sesG,sesP);

%remove neurons that are absent
for intCellType=1:numel(sSesAggregate.cellType)
	strCellType = sSesAggregate.cellType{intCellType};
	if isfield(sesG,strCellType)
		indKeepG = ~cellfun(@strcmpi,cellfill('absent',[1 numel(sesG.(strCellType))]),{sesG.(strCellType).strPresence});
		indKeepP = ~cellfun(@strcmpi,cellfill('absent',[1 numel(sesP.(strCellType))]),{sesP.(strCellType).strPresence});
		
		%check which to keep
		vecKeepAgg = zeros([1 numel(sSesAggregate.(strCellType))]);
		vecKeepAgg(indKeepG) = 1;
		vecKeepAgg(indKeepP) = vecKeepAgg(indKeepP) + 1;
		indKeepAgg = vecKeepAgg==2;
		%remove entries
		sSesAggregate.(strCellType) = sSesAggregate.(strCellType)(indKeepAgg);
	end
end

%% done, set switch
boolLoad=false;