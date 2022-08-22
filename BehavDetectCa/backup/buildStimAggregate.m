function sAggregate = buildStimAggregate(sIn)
	%UNTITLED2 Summary of this function goes here
	%   Detailed explanation goes here
	
	if ~exist('sIn','var'), sIn = struct;end
	if isfield(sIn,'strSes'), strSes = sIn.strSes;else strSes = '20140207';end
	if isfield(sIn,'strMasterPath'), strMasterPath = sIn.strMasterPath;else strMasterPath = 'D:\Data\Processed\imagingdata';end
	if isfield(sIn,'vecRecordings'), vecRecordings = sIn.vecRecordings;else vecRecordings = 1:8;end
	
	%create links to files
	intCounter = 0;
	cellAggregate = cell(1,length(vecRecordings));
	for intRec=vecRecordings
		intCounter = intCounter + 1;
		cellAggregate{intCounter} = sprintf('%s%s%s%sxyt%02d%s%sxyt%02d_ses.mat',strMasterPath,filesep,strSes,filesep,intRec,filesep,strSes,intRec);
	end
	
	%pre-allocate
	sAggregate = struct;
	
	%aggregate all behavioral sessions
	for intFile=1:numel(cellAggregate)
		%load data
		sLoad = load(cellAggregate{intFile});
		
		%get saved structure
		cellFields = fieldnames(sLoad);
		sTemp = sLoad.(cellFields{1});
		
		%get data location
		if isfield(sTemp,'structStim') && isfield(sTemp.structStim,'Contrast')
			sData = sTemp.structStim;
		elseif isfield(sTemp,'Contrast')
			sData = sTemp;
		end
		
		%add structure to aggregate
		sAggregate = doAggregateStructures(sAggregate,sData);
	end
end

