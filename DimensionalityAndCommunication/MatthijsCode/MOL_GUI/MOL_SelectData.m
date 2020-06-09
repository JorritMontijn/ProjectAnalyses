%Get data script - Should be supergeneral
function [Data] = MOL_SelectData() %#ok<STOUT>
RootDataDir = 'E:\Data\';

%% Create Dataselector Figure
global dataselector dataselector_handles
close all
dataselector = figure('Name','',...
    'units','normalized',...
    'Position',[0 0 0.4 0.91],...
    'Color',[1 1 1],...
    'Resize','off',...
    'NumberTitle','off'...
    ,'tag','_main');
dataselector_handles = struct();
dataselector_handles.header_string = uicontrol('Parent',dataselector,'Style','text','units','normalized',...
    'Position',[0.1 0.9 0.8 0.1],'String','SELECT DATA',...
    'BackgroundColor',[1 1 1],'fontsize',20);

%% Which-Projects-Button
Projects = {'PMAL' 'CHDET' 'MODDISCR' 'SIFI' 'UPDOWN'};

x_pos   = 0.1;
y_pos   = 0.675;
x_width = 0.25;
y_width = 0.1;

for proj = 1:length(Projects)
    callbackstr = strcat('global dataselector_handles; dataselector_handles.project = ''',Projects{proj},'''; uiresume(gcf);');
    dataselector_handles.(Projects{proj})    = uicontrol('Parent',dataselector,'Style', 'push','units','normalized',...
        'Position',[x_pos y_pos x_width y_width],'String',Projects{proj},'fontsize',14,...
        'backgroundcolor',[0.2,0.6,0.9],'Callback',callbackstr); x_pos = x_pos+x_width*1.1;
    if x_pos+x_width>1
        x_pos   = 0.1;
        y_pos   = y_pos - y_width*1.1;
    end
end

%Wait for data button to be pressed: callback is continue with script
uiwait(dataselector)

for proj = 1:length(Projects)
    delete(dataselector_handles.(Projects{proj}))
end; pause(0.05)

%% Which-Experiment-Button
RootProjectDir  = fullfile(RootDataDir,dataselector_handles.project);
files           = dir(RootProjectDir);
files(1:2)      = [];
dirFlags        = [files.isdir];
Experiments     = {files(dirFlags).name};
Experiments     = Experiments(~strcmp(Experiments,'RawData'));

x_pos   = 0.1;
y_pos   = 0.675;
x_width = 0.25;
y_width = 0.1;

for exp = 1:length(Experiments)
    callbackstr = strcat('global dataselector_handles; dataselector_handles.experiment = ''',Experiments{exp},'''; uiresume(gcf);');
    dataselector_handles.(Experiments{exp})    = uicontrol('Parent',dataselector,'Style', 'push','units','normalized',...
        'Position',[x_pos y_pos x_width y_width],'String',Experiments{exp},'fontsize',14,...
        'backgroundcolor',[0.2,0.9,0.4],'Callback',callbackstr); x_pos = x_pos+x_width*1.1;
    if x_pos+x_width>1
        x_pos   = 0.1;
        y_pos   = y_pos - y_width*1.1;
    end
end

%Wait for data button to be pressed: callback is continue with script
uiwait(dataselector)

for exp = 1:length(Experiments)
    delete(dataselector_handles.(Experiments{exp}))
end; pause(0.05)

% Get a list of all files and folders in this folder.
RootExpDir = fullfile(RootDataDir,dataselector_handles.project,dataselector_handles.experiment);
files       = dir(RootExpDir);
files(1:2)  = [];
dirFlags    = [files.isdir];
dataselector_handles.AllMice     = {files(dirFlags).name};

%% Selection boxes for spikeData and lfpData
x_pos   = 0.1;
y_pos   = 0.9;
x_width = 0.3;
y_width = 0.05;

dataselector_handles.tickboxspikedata = uicontrol('Parent',dataselector,'Style','checkbox','units','normalized','enable','on',...
    'BackgroundColor',[1 1 1],'Position',[x_pos y_pos x_width y_width],'String','Load Spikes','fontsize',10,'value',0);
x_pos = x_pos + x_width;
dataselector_handles.tickboxlfpdata = uicontrol('Parent',dataselector,'Style','checkbox','units','normalized','enable','on',...
    'BackgroundColor',[1 1 1],'Position',[x_pos y_pos x_width y_width],'String','Load LFP','fontsize',10,'value',0);
x_pos = x_pos + x_width;
dataselector_handles.tickboxpupildata = uicontrol('Parent',dataselector,'Style','checkbox','units','normalized','enable','on',...
    'BackgroundColor',[1 1 1],'Position',[x_pos y_pos x_width y_width],'String','Load Pupil','fontsize',10,'value',0);

%% Show mouse data overview in a table:

for mouse = 1:length(dataselector_handles.AllMice);
    if exist(fullfile(RootExpDir,dataselector_handles.AllMice{mouse}),'dir')
        % Get a list of all files and folders in this folder.
        files       = dir(fullfile(RootExpDir,dataselector_handles.AllMice{mouse}));
        files(1:2)  = [];
        dirFlags    = [files.isdir];
        dataselector_handles.Sessions(1:sum(dirFlags),mouse)     = {files(dirFlags).name};
    end
end


SesTable = uitable('Parent', dataselector, 'Data', dataselector_handles.Sessions,'units','normalized',...
    'Position', [0.05 0.05 0.9 0.85],'ColumnName',dataselector_handles.AllMice,'RowName','numbered',...
    'CellSelectionCallback',@(src,evnt)set(src,'UserData',evnt.Indices));
    
%% Load Button + Wait for press

dataselector_handles.load_data = uicontrol('Parent',dataselector,'Style', 'push','units','normalized',...
    'Position',[0.7 0.075 0.25 0.1],'String','Load this data','fontsize',14,...
    'backgroundcolor',[0.2,0.5,0.9],'Callback','uiresume(gcf);');
uiwait(dataselector) %Wait for load button to be pressed: callback is continue with script

%% Continue with getting the data
% Initialize types of data: (DataTypes   = {'sessionData' 'trialData' 'spikeData' 'lfpData' 'pupilData'};
sessionData     = struct();
trialData       = struct();
spikeData       = struct();
lfpData         = struct();
pupilData       = struct();

%% Get the selected data from the selected table:
SesSelec = get(SesTable,'UserData');

for ses = 1:size(SesSelec,1)
    curdir = fullfile(RootExpDir,dataselector_handles.AllMice{SesSelec(ses,2)},dataselector_handles.Sessions{SesSelec(ses,1),SesSelec(ses,2)},dataselector_handles.experiment);
    if ~exist(curdir,'dir')
        curdir = fullfile(RootExpDir,dataselector_handles.AllMice{SesSelec(ses,2)},dataselector_handles.Sessions{SesSelec(ses,1),SesSelec(ses,2)});
    end
    
    if exist(curdir,'dir'); %It's selected and available
        
        %sessionData
        if exist(fullfile(curdir,'sessionData.mat'),'file')
            loadstruct              = load(fullfile(curdir,'sessionData.mat'));
            tempsessionData         = loadstruct.sessionData;
            sessionData             = AppendStruct(sessionData,tempsessionData);
        end
        
        %trialData
        if exist(fullfile(curdir,'trialData.mat'),'file')
            loadstruct              = load(fullfile(curdir,'trialData.mat'));
            temptrialData           = loadstruct.trialData;
            trialData               = AppendStruct(trialData,temptrialData);
        end
        
        %spikeData
        if get(dataselector_handles.tickboxspikedata,'value') && exist(fullfile(curdir,'spikeData.mat'),'file')
            loadstruct              = load(fullfile(curdir,'spikeData.mat'));
            tempspikeData         = loadstruct.spikeData;
            spikeData             = AppendStruct(spikeData,tempspikeData);
        end
        
        %lfpData
        if get(dataselector_handles.tickboxlfpdata,'value') && exist(fullfile(curdir,'lfpData.mat'),'file')
            loadstruct          = load(fullfile(curdir,'lfpData.mat'));
            templfpData         = loadstruct.lfpData; 
            clear loadstruct;
            
            %to append struct of different signal length
            if ~isa(templfpData.signal,'cell')
               templfpData.signal= mat2cell(templfpData.signal, ones(1,length(templfpData.signal(:,1))), length(templfpData.signal(1,:)));
            end
            
            lfpData             = AppendStruct(lfpData,templfpData);
        end
        
        %pupilData
        if get(dataselector_handles.tickboxpupildata,'value') && exist(fullfile(curdir,'pupilData.mat'),'file')
            loadstruct          = load(fullfile(curdir,'pupilData.mat'));
            temppupilData       = loadstruct.pupilData;
            pupilData           = AppendStruct(pupilData,temppupilData);
        end
        
    end
end

%% Print output:
nMice           = numel(unique(sessionData.mousename));
nSessions       = size(sessionData.session_ID,1);
if isfield(trialData,'session_ID')
    nTrials         =  size(trialData.session_ID,1);
else nTrials = 0;
end
if isfield(spikeData,'session_ID')
    nCells = size(spikeData.session_ID,1);
else nCells = 0;
end
fprintf('\n\nLoaded %s - %s:\n %d Sessions for %d Mice\n A total of %d trials and %d neurons\n\n',dataselector_handles.project,dataselector_handles.experiment,nSessions,nMice, nTrials, nCells)

%% Assign output argument Data
Datatypes = {'sessionData' 'trialData' 'spikeData' 'lfpData' 'pupilData'}; 
for dt = 1:length(Datatypes)
    if exist(Datatypes{dt},'var') && eval(strcat('~all(structfun(@isempty,',Datatypes{dt},'))'))
        eval(strcat('Data.',Datatypes{dt},' = ',Datatypes{dt},';'))
    end
end

%% Close the data selection figure
close(dataselector)
clear global

end
