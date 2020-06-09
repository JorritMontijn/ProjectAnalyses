function [plotydata,plotxaxisdata,colors,labels,errorflag] = MOL_GUI_PrepareDataPlot()

global Data plot_fig_handles_left
errorflag = 0;
colors = {};
labels = {};

%% Set standard colors:
for i = 1:5
    figsettings.Colorset.visualInt{i}   = [1-i*0.2 1-i*0.2 1];
    figsettings.Colorset.audioInt{i}    = [1 1-i*0.2 1-i*0.2];
    figsettings.Colorset.visualIntNorm{i} = [1-i*0.2 1-i*0.2 1];
    figsettings.Colorset.audioIntNorm{i} = [1 1-i*0.2 1-i*0.2];
end

for i = 1:30
    figsettings.Colorset.Standard{i} = mod([i*0.25 i*0.3 i*0.4],1);
end

figsettings.Colorset.trialType      = {[1 0.2 0.2] [0.9 0 0.9] [0.5 0.5 0.5] [0.2 0.2 1]};
figsettings.Colorset.trialType      = {[1 0.2 0.2] [0.5 0.5 0.5] [0.2 0.2 1]};
figsettings.Colorset.trialType      = {[0.5 0.5 0.5] [0.2 0.2 1] [1 0.2 0.2] };

% still to append:
% -visualOri
% -audioFreq
% -correctResponse

% trialtypes      = {'A'          'V'          'P'            'C'         'Y'              'X'};
% trialcolors     = {[1 0.2 0.2], [0.2 0.2 1], [0.5 0.5 0.5], [0.9 0 0.9] ,[0.95 0.5 0.5] ,[0.5 0.5 0.95]};

%% Get the variables of interest:
figsettings.Data_type   = plot_fig_handles_left.data_type.String{get(plot_fig_handles_left.data_type,'Value')};
figsettings.Ydatavar    = plot_fig_handles_left.YDATA.String{get(plot_fig_handles_left.YDATA,'Value')};
figsettings.Xaxisvar    = plot_fig_handles_left.XDATA.String{get(plot_fig_handles_left.XDATA,'Value')};
if ~isfield(Data.(figsettings.Data_type),figsettings.Ydatavar) || ~isfield(Data.(figsettings.Data_type),figsettings.Xaxisvar)
    errorflag = 1;
    errortext = 'Update X Y variable selection';
end

% if ~(get(plot_fig_handles_left.Normalization,'Value')==1);
%     Normalize = 1;
% else Normalize = 0;
% end

%% Go through all session selections (Checkboxes)
Datatypes = fieldnames(Data);
if ~errorflag
    for i = 1:length(Datatypes)
        Datatype            = Datatypes{i};
        CheckBoxFieldNames  = fieldnames(plot_fig_handles_left.Checkboxes.(Datatype));
        selectedindex       = {};
        selectedindex{1}    = true(length((Data.(Datatype).(CheckBoxFieldNames{1}))),1);
        j = 2;
        for checkbox = 1:length(CheckBoxFieldNames)
            try checkthischeckbox = sum(cell2mat(get(plot_fig_handles_left.Checkboxes.(Datatype).(CheckBoxFieldNames{checkbox})(1:end),'Value')))<length(plot_fig_handles_left.Checkboxes.(Datatype).(CheckBoxFieldNames{checkbox}))...
                    || get(plot_fig_handles_left.Toggle.(Datatype).(CheckBoxFieldNames{checkbox}),'Value');
            catch ME %#ok<NASGU>
                checkthischeckbox = 0;
            end
            if  checkthischeckbox %Get only if toggle or not all checkboxes
                checkboxindex = false;
                for checkboxcond = 1:length(plot_fig_handles_left.Checkboxes.(Datatype).(CheckBoxFieldNames{checkbox}))
                    if get(plot_fig_handles_left.Checkboxes.(Datatype).(CheckBoxFieldNames{checkbox})(checkboxcond),'Value')
                        
                        if isfield(Data.(Datatype),CheckBoxFieldNames{checkbox})
                            if iscell(Data.(Datatype).(CheckBoxFieldNames{checkbox}))
                                val = get(plot_fig_handles_left.Checkboxes.(Datatype).(CheckBoxFieldNames{checkbox})(checkboxcond),'String');
                                tempindex = strcmp(Data.(Datatype).(CheckBoxFieldNames{checkbox}),val);
                            elseif isnumeric(Data.(Datatype).(CheckBoxFieldNames{checkbox}))
                                val = str2double(get(plot_fig_handles_left.Checkboxes.(Datatype).(CheckBoxFieldNames{checkbox})(checkboxcond),'String'));
                                tempindex = ismember(Data.(Datatype).(CheckBoxFieldNames{checkbox}),val);
                            end
                        end
                        checkboxindex = checkboxindex | tempindex;
                    end
                end
                selectedindex{j} = checkboxindex; %#ok<AGROW>
                j = j+1;
            end
        end
        SelecIdx.(Datatype) = all(cell2mat(selectedindex),2);
    end
end

%% Trial Data selection:
% Select only trials 1:end-20 or something
lastn = 20;
if ~errorflag && isfield(Data,'trialData') && get(plot_fig_handles_left.lastn_button,'value')
    sessionsel = unique(Data.trialData.session_ID(SelecIdx.trialData));
    for ses = sessionsel'
        totaltrials = max(Data.trialData.trialNum(strcmp(Data.trialData.session_ID,ses)));
        SelecIdx.trialData(strcmp(Data.trialData.session_ID,ses) & ismember(Data.trialData.trialNum,(totaltrials-lastn-1):totaltrials)) = 0;
    end
end

%% Combine the selection from sessionData to other datatypes:
if ~errorflag
    if strcmp(figsettings.Data_type,'sessionData')
        SelecIdx = rmfield(SelecIdx,Datatypes(~ismember(Datatypes,{'sessionData'})));
    else
        if ~isempty(SelecIdx.sessionData)
            SelecIdx.(figsettings.Data_type) = SelecIdx.(figsettings.Data_type) & ...
                ismember(Data.(figsettings.Data_type).session_ID,...
                Data.sessionData.session_ID(SelecIdx.sessionData));
        end
        SelecIdx = rmfield(SelecIdx,Datatypes(~ismember(Datatypes,figsettings.Data_type)));
    end
end
%% Combine all positive indices and make index per session
if ~errorflag
    nSession = 1;
    sessionsel = unique(Data.(figsettings.Data_type).session_ID(SelecIdx.(figsettings.Data_type)));
    totalselectedindexpersession = cell(1,length(sessionsel));
    for ses = sessionsel'
        totalselectedindexpersession{nSession} = SelecIdx.(figsettings.Data_type) & strcmp(Data.(figsettings.Data_type).session_ID,ses);
        nSession = nSession + 1;
    end
end


%% Do check whether anything is selected at all:
if ~errorflag
    try errorflag = ~any(any(cell2mat(totalselectedindexpersession)));
    catch ME %#ok<NASGU>
        errorflag = 1;
    end
    if errorflag
        errortext = sprintf('Probably no data selected');
    end
end


%% Go through all splittings (Toggles)
%Make an index per toggle combination which data to select
if ~errorflag
    for i = 1:length(Datatypes)
        Datatype            = Datatypes{i};
        
        ToggleFieldNames = fieldnames(plot_fig_handles_left.Toggle.(Datatype));
        SplitIdx.(Datatype) = {};
        ColorIdx.(Datatype) = {};
        LabelIdx.(Datatype) = {};
        %         normalizex = {};
        splitrow = 1;
        splitcolumn = 1;
        for splitfield = 1:length(ToggleFieldNames)
            if length(plot_fig_handles_left.Checkboxes.(Datatype).(ToggleFieldNames{splitfield}))>1
                if get(plot_fig_handles_left.Toggle.(Datatype).(ToggleFieldNames{splitfield}),'Value') ...
                        && sum(cell2mat(get(plot_fig_handles_left.Checkboxes.(Datatype).(ToggleFieldNames{splitfield})(1:end),'Value')))>=2
                    for checkbox = 1:length(plot_fig_handles_left.Checkboxes.(Datatype).(ToggleFieldNames{splitfield}))
                        if get(plot_fig_handles_left.Checkboxes.(Datatype).(ToggleFieldNames{splitfield})(checkbox),'Value')
                            
                            val = get(plot_fig_handles_left.Checkboxes.(Datatype).(ToggleFieldNames{splitfield})(checkbox),'String');
                            if iscell(Data.(Datatype).(ToggleFieldNames{splitfield}))
                                SplitIdx.(Datatype){splitrow,splitcolumn} = strcmp(Data.(Datatype).(ToggleFieldNames{splitfield}),val);
                                LabelIdx.(Datatype){splitrow,splitcolumn} = sprintf('%s%s',ToggleFieldNames{splitfield}(1:3),val);
                                try ColorIdx.(Datatype){splitrow,splitcolumn} = figsettings.Colorset.(ToggleFieldNames{splitfield}){checkbox};
                                catch
                                    ColorIdx.(Datatype){splitrow,splitcolumn} = figsettings.Colorset.Standard{checkbox}; %Use standard color if no color is preset
                                end
                            elseif isnumeric(Data.(Datatype).(ToggleFieldNames{splitfield}))
                                
                                SplitIdx.(Datatype){splitrow,splitcolumn} = ismember(Data.(Datatype).(ToggleFieldNames{splitfield}),str2double(val));
                                LabelIdx.(Datatype){splitrow,splitcolumn} = sprintf('%s%d',ToggleFieldNames{splitfield}(1:3),str2double(val));
                                try ColorIdx.(Datatype){splitrow,splitcolumn} = figsettings.Colorset.(ToggleFieldNames{splitfield}){checkbox};
                                catch
                                    ColorIdx.(Datatype){splitrow,splitcolumn} = figsettings.Colorset.Standard{checkbox}; %Use standard color if no color is preset
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
    end
end

%% Combine the splits for sessionData splits
if ~errorflag
    if strcmp(figsettings.Data_type,'sessionData')
        SplitIdx = rmfield(SplitIdx,Datatypes(~ismember(Datatypes,{'sessionData'})));
    else
        if ~isempty(SplitIdx.sessionData)
            nsplitrows = size(SplitIdx.(figsettings.Data_type),1);
            for i = 1:length(SplitIdx.sessionData)
                SplitIdx.(figsettings.Data_type){nsplitrows+1,i} = ismember(Data.(figsettings.Data_type).session_ID,...
                    Data.sessionData.session_ID(SplitIdx.sessionData{i}));
                LabelIdx.(figsettings.Data_type){nsplitrows+1,i} = LabelIdx.sessionData{i};
            end
        end
        SplitIdx = rmfield(SplitIdx,Datatypes(~ismember(Datatypes,figsettings.Data_type)));
    end
end

%% Now get data according to both selected and splitted preferences
if ~errorflag
    [splitrow,splitcol] = size(SplitIdx.(figsettings.Data_type));
%     MouseIDs = {};
    if splitrow == 0                            %If no splitting is selected
        colors{1} = [0.3 0.3 0.5];
        labels = {figsettings.Ydatavar};
        for nSession = 1:length(totalselectedindexpersession)
            totalindex = totalselectedindexpersession{nSession};
            if any(totalindex)
                currYdata = Data.(figsettings.Data_type).(figsettings.Ydatavar)(totalindex)';
                currXdata = Data.(figsettings.Data_type).(figsettings.Xaxisvar)(totalindex)';
                
                try
                    plotydata{1} = ([plotydata{1}; currYdata NaN(1,10000-length(currYdata))]);
                catch
                    plotydata{1} = [currYdata NaN(1,10000-length(currYdata))];
                end
                
                try
                    plotxaxisdata{1} = ([plotxaxisdata{1}; currXdata NaN(1,10000-length(currXdata))]);
                catch
                    plotxaxisdata{1} = [currXdata NaN(1,10000-length(currXdata))];
                end

            end
        end

    elseif splitrow == 1    %If 1 splitting is selected, loop through all the condition indices
        for splitc = 1:splitcol
            colors{1,splitc} = ColorIdx.(figsettings.Data_type){1,splitc};
            labels{1,splitc} = LabelIdx.(figsettings.Data_type){1,splitc};
            for nSession = 1:length(totalselectedindexpersession)
                totalindex = totalselectedindexpersession{nSession} & SplitIdx.(figsettings.Data_type){1,splitc};
                if any(totalindex)
                    currYdata = Data.(figsettings.Data_type).(figsettings.Ydatavar)(totalindex)';
                    currXdata = Data.(figsettings.Data_type).(figsettings.Xaxisvar)(totalindex)';

                    try
                        plotydata{1,splitc} = ([plotydata{1,splitc}; currYdata NaN(1,10000-length(currYdata))]);
                    catch
                        plotydata{1,splitc} = ([currYdata NaN(1,10000-length(currYdata))]);
                    end
                    
                    try
                        plotxaxisdata{1,splitc} = ([plotxaxisdata{1,splitc}; currXdata NaN(1,10000-length(currXdata))]);
                    catch
                        plotxaxisdata{1,splitc} = ([currXdata NaN(1,10000-length(currXdata))]);
                    end
                end
            end
        end

    elseif splitrow == 2   %If two splittings are selected, make combinatorial indices

        %Take longest split first:
        for splitr = 1:splitrow
            nConditions(splitr) = sum(~cellfun(@isempty, SplitIdx.(figsettings.Data_type)(splitr,:)));
        end
        [nConditions,sorted] = sort(nConditions,'descend');
        SplitIdx.(figsettings.Data_type) = SplitIdx.(figsettings.Data_type)(sorted,:);
        ColorIdx.(figsettings.Data_type) = ColorIdx.(figsettings.Data_type)(sorted,:);
        LabelIdx.(figsettings.Data_type) = LabelIdx.(figsettings.Data_type)(sorted,:);

        for nCondition1 = 1:nConditions(1)
            for nCondition2 = 1:nConditions(2)
                colors{nCondition2,nCondition1} = mean([ColorIdx{1,nCondition1}; ColorIdx{2,nCondition2}]); %Average the 2 colors
                labels{nCondition2,nCondition1} = strcat(LabelIdx.(figsettings.Data_type){1,nCondition1},LabelIdx.(figsettings.Data_type){2,nCondition2}); %concatenate labels
                for nSession = 1:length(totalselectedindexpersession)
                    totalindex = totalselectedindexpersession{nSession} & SplitIdx.(figsettings.Data_type){1,nCondition1} & SplitIdx.(figsettings.Data_type){2,nCondition2};
                    if any(totalindex)
                        currYdata = Data.(figsettings.Data_type).(figsettings.Ydatavar)(totalindex)';
                        currXdata = Data.(figsettings.Data_type).(figsettings.Xaxisvar)(totalindex)';
                        
                        try
                            plotydata{nCondition2,nCondition1} = ([plotydata{nCondition2,nCondition1}; currYdata NaN(1,10000-length(currYdata))]);
                        catch
                            plotxaxisdata{nCondition2,nCondition1} = ([plotxaxisdata{nCondition2,nCondition1}; currXdata NaN(1,10000-length(currXdata))]);
                        end
                        try
                            plotydata{nCondition2,nCondition1} = ([currYdata NaN(1,10000-length(currYdata))]);
                        catch
                            plotxaxisdata{nCondition2,nCondition1} = ([currXdata NaN(1,10000-length(currXdata))]);
                        end
                    end
                end
            end
        end
    else errorflag = 1;
        errortext = 'More than two splits selected';
    end
end

%% Section to test whether the selected data is valid
if ~errorflag %If not already a flag
    [splitrow,~] = size(plotydata);
    if splitrow == 0
        errorflag = all(all(isnan(plotydata{1})));
        if errorflag
            errortext = 'No data selected';
        end
    end
end

%% If errorflag was set, show errortext in stead of header text
if errorflag
    set(plot_fig_handles_left.header,'String',errortext,'foregroundcolor',[1 0 0],'fontsize',14);
    pause(2)
    set(plot_fig_handles_left.header,'String','Figure Settings','foregroundcolor',[0 0 0],'fontsize',20);
    plotydata = {};
    plotxaxisdata = {};
    colors = {};
    labels = {};
end

end


