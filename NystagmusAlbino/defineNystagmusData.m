%% define data
sSession = [];
%% 20190315
%query tags
sSession(end+1).MouseType = 'WT';
sSession(end).Area = 'SC';
%identifier tags
sSession(end).Mouse = 'MB2';
sSession(end).Date = '20190315';
sSession(end).DepthPerCh = (32:-1:1)*25;
%pop 1
sSession(end).Recording(1).Location = 1;
sSession(end).Recording(end).Depth = 2000;
sSession(end).Recording(end).Stim = 'RF';
sSession(end).Recording(end).Block = 3;
sSession(end).Recording(end+1).Location = sSession(end).Recording(end).Location;
sSession(end).Recording(end).Depth = sSession(end).Recording(end-1).Depth;
sSession(end).Recording(end).Stim = 'DG';
sSession(end).Recording(end).Block = 4;
%pop 2
sSession(end).Recording(end+1).Location = 1;
sSession(end).Recording(end).Depth = 2600;
sSession(end).Recording(end).Stim = 'RF';
sSession(end).Recording(end).Block = 6;
sSession(end).Recording(end+1).Location = sSession(end).Recording(end).Location;
sSession(end).Recording(end).Depth = sSession(end).Recording(end-1).Depth;
sSession(end).Recording(end).Stim = 'DG';
sSession(end).Recording(end).Block = 7;
%pop 3
sSession(end).Recording(end+1).Location = 2;
sSession(end).Recording(end).Depth = 2000;
sSession(end).Recording(end).Stim = 'RF';
sSession(end).Recording(end).Block = 9;
sSession(end).Recording(end+1).Location = sSession(end).Recording(end).Location;
sSession(end).Recording(end).Depth = sSession(end).Recording(end-1).Depth;
sSession(end).Recording(end).Stim = 'DG';
sSession(end).Recording(end).Block = 10;

%% 20190314
%query tags
sSession(end+1).MouseType = 'WT';
sSession(end).Area = 'V1';
%identifier tags
sSession(end).Mouse = 'MB2';
sSession(end).Date = '20190314';
sSession(end).DepthPerCh = (32:-1:1)*25;
%pop 1
sSession(end).Recording(1).Location = 1;
sSession(end).Recording(end).Depth = 700;
sSession(end).Recording(end).Stim = 'RF';
sSession(end).Recording(end).Block = 1;
sSession(end).Recording(end+1).Location = sSession(end).Recording(end).Location;
sSession(end).Recording(end).Depth = sSession(end).Recording(end-1).Depth;
sSession(end).Recording(end).Stim = 'DG';
sSession(end).Recording(end).Block = 2;

%% 20190508
%query tags
sSession(end+1).MouseType = 'WT';
sSession(end).Area = 'SC';
%identifier tags
sSession(end).Mouse = 'MB5';
sSession(end).Date = '20190508';
sSession(end).DepthPerCh = (32:-1:1)*25;
%pop 1; STILL TO CLUSTER CUT!!!
sSession(end).Recording(1).Location = 1;
sSession(end).Recording(end).Depth = 2200;
sSession(end).Recording(end).Stim = 'RF';
sSession(end).Recording(end).Block = '9_still_to_cut';
%pop 2
sSession(end).Recording(1).Location = 1;
sSession(end).Recording(end).Depth = 2200;
sSession(end).Recording(end).Stim = 'DG';
sSession(end).Recording(end).Block = '11-12';

%% 20190510
%query tags
sSession(end+1).MouseType = 'WT';
sSession(end).Area = 'SC';
%identifier tags
sSession(end).Mouse = 'MB5';
sSession(end).Date = '20190510';
sSession(end).DepthPerCh = (32:-1:1)*25;
%pop 1;
sSession(end).Recording(1).Location = 1;
sSession(end).Recording(end).Depth = 2500;
sSession(end).Recording(end).Stim = 'DG-RF';
sSession(end).Recording(end).Block = '15-16-17';

%% 20190515
%query tags
sSession(end+1).MouseType = 'WT';
sSession(end).Area = 'V1';
%identifier tags
sSession(end).Mouse = 'MB4';
sSession(end).Date = '20190515';
sSession(end).DepthPerCh = (32:-1:1)*25;
%pop 1, ori tuned cells!
sSession(end).Recording(1).Location = 1;
sSession(end).Recording(end).Depth = 800;
sSession(end).Recording(end).Stim = 'DG';
sSession(end).Recording(end).Block = '2';

%query tags
sSession(end+1).MouseType = 'WT';
sSession(end).Area = 'NOT';
%identifier tags
sSession(end).Mouse = 'MB4';
sSession(end).Date = '20190515';
sSession(end).DepthPerCh = (32:-1:1)*25;
%pop 1, unstable RF
sSession(end).Recording(1).Location = 1;
sSession(end).Recording(end).Depth = 2300;
sSession(end).Recording(end).Stim = 'RF';
sSession(end).Recording(end).Block = '8';
%pop 2, a couple of cells
sSession(end).Recording(end+1).Location = 1;
sSession(end).Recording(end).Depth = 2300;
sSession(end).Recording(end).Stim = 'DG';
sSession(end).Recording(end).Block = '9-11';

%% 20190516
%query tags
sSession(end+1).MouseType = 'WT';
sSession(end).Area = 'V1';
%identifier tags
sSession(end).Mouse = 'MB4';
sSession(end).Date = '20190516';
sSession(end).DepthPerCh = (32:-1:1)*25;
%pop 1
sSession(end).Recording(1).Location = 1;
sSession(end).Recording(end).Depth = 800;
sSession(end).Recording(end).Stim = 'Rff';
sSession(end).Recording(end).Block = '1';
%pop 2
sSession(end).Recording(end+1).Location = 1;
sSession(end).Recording(end).Depth = 800;
sSession(end).Recording(end).Stim = 'DG';
sSession(end).Recording(end).Block = '2';
%pop 2
sSession(end).Recording(end+1).Location = 1;
sSession(end).Recording(end).Depth = 800;
sSession(end).Recording(end).Stim = 'RF';
sSession(end).Recording(end).Block = '3';