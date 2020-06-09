function [sessionData,trialData,spikeData,lfpData,pupilData] = MOL_getDataInputs(params)

% Get a list of all files and folders in this folder.
RootExpDir  = fullfile(params.RootDataDir,params.Project,params.Protocol);

% Initialize types of data: (DataTypes   = {'sessionData' 'trialData' 'spikeData' 'lfpData' 'pupilData'};
sessionData     = struct();
trialData       = struct();
spikeData       = struct();
lfpData         = struct();
pupilData       = struct();

for iMouse = 1:length(params.Animals)
    if exist(fullfile(RootExpDir,params.Animals{iMouse}),'dir')
        
        if ~isfield(params,'Sessions')
            files       = dir(fullfile(RootExpDir,params.Animals{iMouse}));
            files(1:2)  = [];
            dirFlags    = [files.isdir];
            params.Sessions    = {files(dirFlags).name};
        end
        
        for iSes = 1:length(params.Sessions)
            curdir = fullfile(RootExpDir,params.Animals{iMouse},params.Sessions{iSes});
            
            if exist(curdir,'dir'); %It's selected and available
                
                %sessionData
                if exist(fullfile(curdir,'sessionData.mat'),'file') && params.loadsessionData
                    loadstruct              = load(fullfile(curdir,'sessionData.mat'));
                    tempsessionData         = loadstruct.sessionData;
                    sessionData             = AppendStruct(sessionData,tempsessionData);
                end
                
                %trialData
                if exist(fullfile(curdir,'trialData.mat'),'file') && params.loadtrialData
                    loadstruct              = load(fullfile(curdir,'trialData.mat'));
                    temptrialData           = loadstruct.trialData;
                    trialData               = AppendStruct(trialData,temptrialData);
                end
                
                %spikeData
                if exist(fullfile(curdir,'spikeData.mat'),'file') && params.loadspikeData
                    loadstruct              = load(fullfile(curdir,'spikeData.mat'));
                    tempspikeData         = loadstruct.spikeData;
                    spikeData             = AppendStruct(spikeData,tempspikeData);
                end
                
                %lfpData
                if exist(fullfile(curdir,'lfpData.mat'),'file') && params.loadlfpData
                    loadstruct          = load(fullfile(curdir,'lfpData.mat'));
                    templfpData         = loadstruct.lfpData;
                    lfpData             = AppendStruct(lfpData,templfpData);
                end
                
                %pupilData
                if exist(fullfile(curdir,'pupilData.mat'),'file') && params.loadpupilData
                    loadstruct          = load(fullfile(curdir,'pupilData.mat'));
                    temppupilData       = loadstruct.pupilData;
                    pupilData           = AppendStruct(pupilData,temppupilData);
                end
                
            end
        end
    end
end

end