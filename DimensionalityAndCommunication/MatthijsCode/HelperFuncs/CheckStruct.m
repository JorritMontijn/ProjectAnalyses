function [Flag,FlagMessage] = CheckStruct(struct)
Flag = false;

%% Get struct type
if isfield(struct,'mousename')
    structtype = 'session';
elseif isfield(struct,'trial_ID')
    structtype = 'trial';
elseif isfield(struct,'session_ID')
    structtype = 'lfp';
elseif isfield(struct,'cell_ID')
    structtype = 'spike';
else
    Flag = 1;
    FlagMessage = 'Unknown struct type. Neither session, trial, lfp nor spike';
end

%% Check general struct demands:

%Struct should contain a reference to the session it belongs to:
if ~isfield(struct,'session_ID')
    Flag = 1;
    FlagMessage = 'struct does not contain session_ID field';
end

%% Check struct specific demands:
switch structtype
    case 'session'
        ObligFields = {'mousename'};
        
        
    case 'trial'
        ObligFields = {'trialnum' 'trialStart' 'trialEnd'};
        
        %All fields should be of equal length:
        trialfields = fieldnames(struct);
        ntrials = length(struct.trialnum);
        for trialfield = 1:length(trialfields)
            if ~length(struct.(trialfields{trialfield})) == ntrials
                Flag = 1;
                FlagMessage = sprintf('%s field is not of equal length as the rest.',trialfields{trialfield});
            end
        end
        
    case 'lfp'
        ObligFields = {'fs' 'signal'};
        
        
    case 'spike'
        ObligFields = {'ts' 'wf' 'ch' 'cell_ID'};
        
end

for field = 1:length(ObligFields)
    if ~isfield(struct,ObligFields{field})
        Flag = 1;
        FlagMessage = sprintf('%s field is missing.',ObligFields{field});
    end
end


%% Output:
if ~Flag
    FlagMessage = 'Struct is Okay';
else
    disp(FlagMessage)
end

end