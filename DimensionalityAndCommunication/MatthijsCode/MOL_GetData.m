%Get data script - Should be supergeneral
function [Data] = MOL_GetData(RootDir,Projects,Protocols,Mice,Sessions,Datatypes) %#ok<STOUT>
% Input:
% RootDir           = Harddrive, e.g. 'E:'
% Projects          = all projects which to load from (string or cell array) 
% Protocols         = all protocols which to load from (string or cell array) 
% Mice              = all mice which to load from (string or cell array) 
% Sessions          = subselection of sessions which to load (string or cell array) 
% Datatypes         = cell array with the requested dataypes {'sessionData' 'trialData' 'spikeData' 'lfpData' 'pupilData'}
% Output:
% Data              = struct with all concatenated struct for each datatype requested
%
% Example: [Data] = MOL_GetData('E:','CHDET',{'ChangeDetectionConflict' 'ChangeDetectionConflictDecor'},[],[],{'sessionData' 'trialData'});
% Matthijs oude Lohuis 2019

%% Initialize types of data: (DataTypes   = {'sessionData' 'trialData' 'spikeData' 'lfpData' 'pupilData'};
for iDt = 1:length(Datatypes)
    eval(sprintf('%s = struct();',Datatypes{iDt}))
end

%% Loop over all projects and sessions and get data:
RootDataDir                     = fullfile(RootDir,'Data');
ProjectFolders                  = subselectfolders(RootDataDir,Projects);
for iPr = 1:length(ProjectFolders)
    ProjectDir                  = fullfile(RootDataDir,ProjectFolders{iPr});
    ProtocolFolders           = subselectfolders(ProjectDir,Protocols);
    for iExp = 1:length(ProtocolFolders)
        ProtocolDir           = fullfile(ProjectDir,ProtocolFolders{iExp});
        MouseFolders            = subselectfolders(ProtocolDir,Mice);
        for iMouse = 1:length(MouseFolders)
            MouseDir            = fullfile(ProtocolDir,MouseFolders{iMouse});
            SessionFolders      = subselectfolders(MouseDir,Sessions);
            for iSes = 1:length(SessionFolders)
                SessionDir      = fullfile(MouseDir,SessionFolders{iSes}); %#ok<NASGU>
                for iDt = 1:length(Datatypes)
                    string1     = sprintf('if exist(fullfile(SessionDir,''%s.mat''),''file'');',Datatypes{iDt});
                    string2     = sprintf('loadstruct = load(fullfile(SessionDir,''%s.mat''));',Datatypes{iDt});
                    string3     = sprintf('tempData = loadstruct.%s;',Datatypes{iDt});
                    string4     = sprintf('%s = AppendStruct(%s,tempData); end',Datatypes{iDt},Datatypes{iDt});
					eval([string1 string2 string3 string4])
					Datatypes{iDt};
                end
            end
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
if exist('spikeData','var')
    nCells = size(spikeData.session_ID,1);
else nCells = 0;
end
fprintf('\n\nLoaded:\n n=%d sessions\n n=%d mice\n n=%d trials\n n=%d neurons\n\n',nSessions,nMice,nTrials,nCells)

%% Assign output argument Data
for dt = 1:length(Datatypes)
    eval(strcat('Data.',Datatypes{dt},' = ',Datatypes{dt},';'))
end

end


function outfolders = subselectfolders(rootdir,selection)
files               = dir(rootdir);
if isempty(selection) %Take all the available subfolders:
    files(1:2)      = [];
    dirFlags        = [files.isdir];
    outfolders      = {files(dirFlags).name};
else %Take the available subfolders that match the selection:
    dirFlags        = [files.isdir];
    outfolders      = intersect({files(dirFlags).name}, selection);
end

end
