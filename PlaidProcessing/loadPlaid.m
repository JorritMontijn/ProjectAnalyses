
%% define master source path
strDataSource = 'D:\Data\Processed\imagingdata\';

%% set vars for different data sets
if intLoadSet == 1 && boolLoad
	strDate = '20130612';
	intBlockGratings = 2;
	intBlockPlaids = 1;
	
elseif intLoadSet == 2 && boolLoad
	strDate = '20130625';
	intBlockGratings = 2;
	intBlockPlaids = 1;
	
elseif intLoadSet == 3 && boolLoad
	strDate = '20131016';
	intBlockGratings = 1;
	intBlockPlaids = 2;
	
elseif intLoadSet == 4 && boolLoad
	strDate = '20131022';
	intBlockGratings = 1;
	intBlockPlaids = 2;
	
elseif intLoadSet == 5 && boolLoad
	strDate = '20140129';
	intBlockGratings = 1;
	intBlockPlaids = 2;
	
elseif intLoadSet == 6 && boolLoad
	strDate = '20140314';
	intBlockGratings = 8;
	intBlockPlaids = 9;
	
elseif intLoadSet == 7 && boolLoad
	strDate = '20140423';
	intBlockGratings = 2;
	intBlockPlaids = 1;
	
elseif intLoadSet == 8 && boolLoad
	strDate = '20140425';
	intBlockGratings = 9;
	intBlockPlaids = 10;
end

%% construct target file locations
strBlockGratings = sprintf('xyt%02d',intBlockGratings);
strBlockPlaids = sprintf('xyt%02d',intBlockPlaids);
strGratingFile = strcat(strDataSource,strDate,filesep,strBlockGratings,filesep,strDate,strBlockGratings,'_ses.mat');
strPlaidFile = strcat(strDataSource,strDate,filesep,strBlockPlaids,filesep,strDate,strBlockPlaids,'_ses.mat');
