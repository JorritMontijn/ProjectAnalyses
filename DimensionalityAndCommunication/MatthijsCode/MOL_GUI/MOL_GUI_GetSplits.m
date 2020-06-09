function [splits,colors] = MOL_GUI_GetSplits(temptrialData)
global plot_fig_handles_left

%% Set standard colors:

for i = 1:5
    Colorset.visualInt{i}       = [1-i*0.2 1-i*0.2 1];
    Colorset.audioInt{i}        = [1 1-i*0.2 1-i*0.2];
    Colorset.visualIntNorm{i}   = [1-i*0.2 1-i*0.2 1];
    Colorset.audioIntNorm{i}    = [1 1-i*0.2 1-i*0.2];
end

Colorset.Standard{1} = [0 0 0];
for i = 2:30
    Colorset.Standard{i} = mod([i*0.25 i*0.7 i*0.4],1);
end

Labelset.trialType      = {'A'          'Y'         'V'         'X'         'C'         'Q'         'R'             'P'};
Colorset.trialType      = {[1 0.2 0.2] [1 0.2 0.2] [0.2 0.2 1] [0.2 0.2 1] [0.9 0 0.9] [0.5 0.5 0.5] [0.5 0.5 0.5] [0.5 0.5 0.5]};

Labelset.area           = {'V1'             'PM'            'AL'            'PPC'       'CG1'};
Colorset.area           = {[0.5 0.3 0.3] [0 0.55 0.22] [0.5 0.22 0.55] [0.9 0.6 0] [0.8 0 0.7]};

Colorset.correctResponse= {[0.5 0.5 0.5] [0.2 0.9 0.2]};
Colorset.hasphotostim   = {[0 0 0] [0 0 1]};
Colorset.visualSpeed    = {[0 0.55 0.22] [0.5 0.3 0.3] [0.5 0.22 0.55]};
Colorset.visualOri      = {[0.58 0 0.83] [0.3 0 0.51] [0 0 1] [0 1 0] [1 1 0] [1 0.5 0] [1 0  0]};

Datatype                = 'trialData';

ToggleFieldNames        = fieldnames(plot_fig_handles_left.Toggle.(Datatype));
SplitIdx.(Datatype)     = {};
ColorIdx.(Datatype)     = {};
LabelIdx.(Datatype)     = {};
%         normalizex = {};
splitrow                = 1;
splitcolumn             = 1;
for splitfield = 1:length(ToggleFieldNames)
    if length(plot_fig_handles_left.Checkboxes.(Datatype).(ToggleFieldNames{splitfield}))>1
        if get(plot_fig_handles_left.Toggle.(Datatype).(ToggleFieldNames{splitfield}),'Value') ...
                && sum(cell2mat(get(plot_fig_handles_left.Checkboxes.(Datatype).(ToggleFieldNames{splitfield})(1:end),'Value')))>=2
            for checkbox = 1:length(plot_fig_handles_left.Checkboxes.(Datatype).(ToggleFieldNames{splitfield}))
                if get(plot_fig_handles_left.Checkboxes.(Datatype).(ToggleFieldNames{splitfield})(checkbox),'Value')
                    
                    val = get(plot_fig_handles_left.Checkboxes.(Datatype).(ToggleFieldNames{splitfield})(checkbox),'String');
                    if iscell(temptrialData.(ToggleFieldNames{splitfield}))
                        SplitIdx.(Datatype){splitrow,splitcolumn} = strcmp(temptrialData.(ToggleFieldNames{splitfield}),val);
                        LabelIdx.(Datatype){splitrow,splitcolumn} = sprintf('%s%s',ToggleFieldNames{splitfield}(1:3),val);
                        
                        try
                            if isfield(Labelset,ToggleFieldNames{splitfield})
                            idx = strcmp(val,Labelset.(ToggleFieldNames{splitfield}));
                            ColorIdx.(Datatype){splitrow,splitcolumn} = Colorset.(ToggleFieldNames{splitfield}){idx};
                            else
                                ColorIdx.(Datatype){splitrow,splitcolumn} = Colorset.(ToggleFieldNames{splitfield}){checkbox};
                            end
                        catch
                            ColorIdx.(Datatype){splitrow,splitcolumn} = Colorset.Standard{checkbox}; %Use standard color if no color is preset
                        end
                    elseif isnumeric(temptrialData.(ToggleFieldNames{splitfield}))
                        
                        SplitIdx.(Datatype){splitrow,splitcolumn} = ismember(temptrialData.(ToggleFieldNames{splitfield}),str2double(val));
                        LabelIdx.(Datatype){splitrow,splitcolumn} = sprintf('%s%d',ToggleFieldNames{splitfield}(1:3),str2double(val));
                        try ColorIdx.(Datatype){splitrow,splitcolumn} = Colorset.(ToggleFieldNames{splitfield}){checkbox};
                        catch
                            ColorIdx.(Datatype){splitrow,splitcolumn} = Colorset.Standard{checkbox}; %Use standard color if no color is preset
                        end
                    end
                    splitcolumn = splitcolumn+1;
                end
            end
            splitcolumn = 1;
            splitrow = splitrow+1;
        end
    end
end

if ~isempty(SplitIdx.(Datatype))
    splits = SplitIdx.(Datatype);
    colors = ColorIdx.(Datatype);
else
    splits{1} = true(length(temptrialData.(ToggleFieldNames{1,1})),1);
    colors{1} = Colorset.Standard{1};
end
end