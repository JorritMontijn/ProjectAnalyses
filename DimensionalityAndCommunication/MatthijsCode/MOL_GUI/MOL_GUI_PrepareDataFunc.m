function [sessionData,trialData,spikeData,errorflag] = MOL_GUI_PrepareDataFunc()

lastn = 20;

global Data plot_fig_handles_left
errorflag = 0;
% colors = {};
% labels = {};

% figsettings.Data_type   = plot_fig_handles_left.data_type.String{get(plot_fig_handles_left.data_type,'Value')};
% figsettings.Ydatavar    = plot_fig_handles_left.YDATAVARS{get(plot_fig_handles_left.YDATA,'Value')};
% figsettings.Xaxisvar    = plot_fig_handles_left.XDATAVARS{get(plot_fig_handles_left.XDATA,'Value')};

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
        Indexes.(Datatype) = all(cell2mat(selectedindex),2);
    end
end

%% Trial Data selection:
%Select only trials 1:end-20 or something
if ~errorflag && isfield(Data,'trialData') && get(plot_fig_handles_left.lastn_button,'value')
    sessionsel = unique(Data.trialData.session_ID(Indexes.trialData));
    for ses = sessionsel'
        totaltrials = max(Data.trialData.trialNum(strcmp(Data.trialData.session_ID,ses)));
        Indexes.trialData(strcmp(Data.trialData.session_ID,ses) & ismember(Data.trialData.trialNum,(totaltrials-lastn+1):totaltrials)) = 0;
    end
end

%% Convert sessionselection to trialData spikeData and lfpData selection:
for i = 1:length(Datatypes)
    Datatype            = Datatypes{i};
    %Convert sessionselection to trialData spikeData and lfpData selection:
    Indexes.(Datatype)  = Indexes.(Datatype) & ismember(Data.(Datatype).session_ID,Data.sessionData.session_ID(Indexes.sessionData));
    
    fields = fieldnames(Data.(Datatype));
    for fld = 1:length(fields)
        SelecData.(Datatype).(fields{fld}) = Data.(Datatype).(fields{fld})(Indexes.(Datatype));
    end
end

%% Get output variables:
sessionData     = SelecData.sessionData; %Always
trialData       = SelecData.trialData;%Always
if isfield(SelecData,'spikeData') %If present
    spikeData       = SelecData.spikeData;
else
    spikeData = struct();
end
if isfield(SelecData,'lfpData') %If present
    lfpData       = SelecData.lfpData;
else
    lfpData = struct();
end

% if ~errorflag
%     if isfield(Data.(figsettings.Data_type),'trialNum')
%         selectedindex{end+1} = ismember(Data.(figsettings.Data_type).trialNum,figsettings.trial_range);
%         if ~any(selectedindex{end}) %Check whether trials are not out of range
%             errorflag = 1;
%             errortext = sprintf('Trials out of range for all mice');
%         end
%     end
% end
%
%
%
% selectedindex{end+1} = ismember(Data.(figsettings.Data_type).Mouse,find(selectedmice));
%
% %% Combine all positive indices and make index per session
% if ~errorflag
%     totalselectedindex = all(cell2mat(selectedindex),2);
%     nSession = 1;
%     for nDay = figsettings.day_range
%         for nMouse = figsettings.mouse_range
%             if any(ismember(Data.(figsettings.Data_type).Day,nDay)) && any(ismember(Data.(figsettings.Data_type).Mouse,nMouse))
%                 totalselectedindexpersession{nSession} = totalselectedindex & ismember(Data.(figsettings.Data_type).Day,nDay) & ismember(Data.(figsettings.Data_type).Mouse,nMouse);
%             else
%                 totalselectedindexpersession{nSession} = zeros(length(totalselectedindex),1);
%             end
%             MouseIDpersession(nSession) = nMouse;
%             nSession = nSession + 1;
%         end
%     end
% end
%
% if ~errorflag
%     try errorflag = ~any(any(cell2mat(totalselectedindexpersession)));
%     catch ME %#ok<NASGU>
%         errorflag = 1;
%     end
%     if errorflag
%         errortext = sprintf('Probably no data selected');
%     end
% end
%
% %% Go through all splittings (Toggles)
% %Make an index per toggle combination which data to select
% if ~errorflag
%     ToggleFieldNames = fieldnames(plot_fig_handles_left.Toggle);
%     splitindex = {};
%     colorindex = {};
%     labelindex = {};
%     normalizex = {};
%     splitrow = 1;
%     splitcolumn = 1;
%     for splitfield = 1:length(ToggleFieldNames)
%         if length(plot_fig_handles_left.Checkboxes.(CheckBoxFieldNames{splitfield}))>1
%             if get(plot_fig_handles_left.Toggle.(ToggleFieldNames{splitfield}),'Value') && sum(cell2mat(get(plot_fig_handles_left.Checkboxes.(CheckBoxFieldNames{splitfield})(1:end),'Value')))>=2
%                 for checkbox = 1:length(plot_fig_handles_left.Checkboxes.(ToggleFieldNames{splitfield}))
%                     if get(plot_fig_handles_left.Checkboxes.(CheckBoxFieldNames{splitfield})(checkbox),'Value')
%                         val = str2double(get(plot_fig_handles_left.Checkboxes.(CheckBoxFieldNames{splitfield})(checkbox),'String'));
%                         splitindex{splitrow,splitcolumn} = ismember(Data.(figsettings.Data_type).(CheckBoxFieldNames{splitfield}),val);
%                         labelindex{splitrow,splitcolumn} = sprintf('%s%d',CheckBoxFieldNames{splitfield}(1:3),val);
%                         try colorindex{splitrow,splitcolumn} = figsettings.Colorset.(CheckBoxFieldNames{splitfield}){checkbox};
%                         catch
%                             colorindex{splitrow,splitcolumn} = figsettings.Colorset.Standard{checkbox}; %Use standard color if no color is preset
%                         end
%                         splitcolumn = splitcolumn+1;
%                     end
%                 end
%                 splitcolumn = 1;
%                 splitrow = splitrow+1;
%             end
%         end
%     end
% end
%
% %% For normalization
%
% if ~errorflag
%     [splitrow,splitcol] = size(splitindex);
%     normalizex = zeros(splitrow,splitcol);
%     normalizex(1,1) = 1; %Standard normalization, to first conditions
%     if get(plot_fig_handles_left.Toggle.RewardProb,'Value') && get(plot_fig_handles_left.Toggle.FlippingGamma,'Value')
%         normalizex = [0 1; 0 0]; %90/30 as baseline
%     end
% end
%
% %% Now get data according to both selected and splitted preferences
% if ~errorflag
%     [splitrow,splitcol] = size(splitindex);
%     MouseIDs = {};
%     if splitrow == 0                            %If no splitting is selected
%         colors{1} = [0.3 0.3 0.5];
%         labels = {figsettings.Ydatavar};
%         for nSession = 1:length(totalselectedindexpersession)
%             totalindex = totalselectedindexpersession{nSession};
%             if any(totalindex)
%                 currYdata = Data.(figsettings.Data_type).(figsettings.Ydatavar)(totalindex)';
%                 currXdata = Data.(figsettings.Data_type).(figsettings.Xaxisvar)(totalindex)';
%                 try
%                     plotydata{1} = ([plotydata{1}; currYdata NaN(1,10000-length(currYdata))]);
%                     plotxaxisdata{1} = ([plotxaxisdata{1}; currXdata NaN(1,10000-length(currXdata))]);
%                     MouseIDs{1}(end+1,1) = MouseIDpersession(nSession);
%                 catch
%                     plotydata{1} = [currYdata NaN(1,10000-length(currYdata))];
%                     plotxaxisdata{1} = [currXdata NaN(1,10000-length(currXdata))];
%                     MouseIDs{1}(1,1) = MouseIDpersession(nSession);
%                 end
%
%             end
%         end
%
%     elseif splitrow == 1    %If 1 splitting is selected, loop through all the condition indices
%         for splitc = 1:splitcol
%             colors{1,splitc} = colorindex{1,splitc};
%             labels{1,splitc} = labelindex{1,splitc};
%             for nSession = 1:length(totalselectedindexpersession)
%                 totalindex = totalselectedindexpersession{nSession} & splitindex{1,splitc};
%                 if any(totalindex)
%                     currYdata = Data.(figsettings.Data_type).(figsettings.Ydatavar)(totalindex)';
%                     currXdata = Data.(figsettings.Data_type).(figsettings.Xaxisvar)(totalindex)';
%
%                     try
%                         plotydata{1,splitc} = ([plotydata{1,splitc}; currYdata NaN(1,10000-length(currYdata))]);
%                         plotxaxisdata{1,splitc} = ([plotxaxisdata{1,splitc}; currXdata NaN(1,10000-length(currXdata))]);
%                         MouseIDs{1,splitc}(end+1,1) = MouseIDpersession(nSession);
%                     catch
%                         plotydata{1,splitc} = ([currYdata NaN(1,10000-length(currYdata))]);
%                         plotxaxisdata{1,splitc} = ([currXdata NaN(1,10000-length(currXdata))]);
%                         MouseIDs{1,splitc}(1,1) = MouseIDpersession(nSession);
%                     end
%                 end
%             end
%         end
%
%     elseif splitrow == 2   %If two splittings are selected, make combinatorial indices
%
%         %Take longest split first:
%         for splitr = 1:splitrow
%             nConditions(splitr) = sum(~cellfun(@isempty,splitindex(splitr,:)));
%         end
%         [nConditions,sorted] = sort(nConditions,'descend');
%         splitindex = splitindex(sorted,:);
%         colorindex = colorindex(sorted,:);
%         labelindex = labelindex(sorted,:);
%
%         for nCondition1 = 1:nConditions(1)
%             for nCondition2 = 1:nConditions(2)
%                 colors{nCondition2,nCondition1} = mean([colorindex{1,nCondition1}; colorindex{2,nCondition2}]); %Average the 2 colors
%                 labels{nCondition2,nCondition1} = strcat(labelindex{1,nCondition1},labelindex{2,nCondition2}); %concatenate labels
%                 for nSession = 1:length(totalselectedindexpersession)
%                     totalindex = totalselectedindexpersession{nSession} & splitindex{1,nCondition1} & splitindex{2,nCondition2};
%                     if any(totalindex)
%                         currYdata = Data.(figsettings.Data_type).(figsettings.Ydatavar)(totalindex)';
%                         currXdata = Data.(figsettings.Data_type).(figsettings.Xaxisvar)(totalindex)';
%
%                         try
%                             plotydata{nCondition2,nCondition1} = ([plotydata{nCondition2,nCondition1}; currYdata NaN(1,10000-length(currYdata))]);
%                             plotxaxisdata{nCondition2,nCondition1} = ([plotxaxisdata{nCondition2,nCondition1}; currXdata NaN(1,10000-length(currXdata))]);
%                             MouseIDs{nCondition2,nCondition1}(end+1,1) = MouseIDpersession(nSession);
%                         catch
%                             plotydata{nCondition2,nCondition1} = ([currYdata NaN(1,10000-length(currYdata))]);
%                             plotxaxisdata{nCondition2,nCondition1} = ([currXdata NaN(1,10000-length(currXdata))]);
%                             MouseIDs{nCondition2,nCondition1}(1,1) = MouseIDpersession(nSession);
%                         end
%                     end
%                 end
%             end
%         end
%     else errorflag = 1;
%         errortext = 'More than two splits selected';
%     end
% end
%
% %% Section to test whether the selected data is valid
% if ~errorflag %If not already a flag
%     [splitrow,~] = size(plotydata);
%     if splitrow == 0
%         errorflag = all(all(isnan(plotydata{1})));
%         if errorflag
%             errortext = 'No data selected';
%         end
%     end
% end
%
% %% Normalization
% if Normalize
%     [rows,columns] = size(plotydata);
%     [NormRow, NormCol] = find(normalizex);
%
%     for nUniqueMice = unique(MouseIDs{NormRow,NormCol})';
%         NormAvg(nUniqueMice) = nanmean(nanmean(plotydata{NormRow,NormCol}(MouseIDs{NormRow,NormCol}==nUniqueMice,:),2));
%     end
%
%     for nRow = 1:rows
%         for nCol = 1:columns
%             for nUniqueMice = unique(MouseIDs{nRow,nCol})';
%                 plotydata{nRow,nCol}(MouseIDs{nRow,nCol}==nUniqueMice,:) = plotydata{nRow,nCol}(MouseIDs{nRow,nCol}==nUniqueMice,:)-NormAvg(nUniqueMice);
%             end
%         end
%     end
% end
%
% %% If errorflag was set, show errortext in stead of header text
% if errorflag
%     set(plot_fig_handles_left.header,'String',errortext,'foregroundcolor',[1 0 0],'fontsize',14);
%     pause(2)
%     set(plot_fig_handles_left.header,'String','Figure Settings','foregroundcolor',[0 0 0],'fontsize',20);
%     plotydata = {};
%     plotxaxisdata = {};
%     colors = {};
%     labels = {};
% end

end


