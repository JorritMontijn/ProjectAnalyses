function MOL_GUI_Build

global Data plot_fig plot_fig_handles_left

plot_fig_handles_left= [];
%Define the functions that can be called. The first one is a general
%variable plotting function using either bars, histogram, or xyhistogram
plot_fig_handles_left.plot_types = {'MOL_PlotVar' 'MOL_2ADC_Full' 'MOL_2ADC_Psy' 'MOL_Psy_2Sided'...
    'MOL_Psy_2Sided_SplitHistory' 'MOL_OptoV1_Behavior' 'MOL_PlotTuning' 'MOL_CorrEarlyLate' ...
    'MOL_LickHistogram' 'MOL_SessionRast' 'MOL_PlotPSTH' 'MOL_SnakePlot' 'MOL_ResponseSnake'...
    'MOL_Sortquality' 'MOL_NeuroMetric' 'MOL_PopulationOriTuning' 'MOL_ConflictDominance'...
    'MOL_PCA_trialtypes' 'MOL_dPCA' 'MOL_OptoOnly' 'MOL_OptoMeanResponse'...
    'MOL_PopulationRecField' 'PMAL_PopulationOriTuning' 'MOL_GLM_History'};

% s = what('MOL_Analysis');
% dirinfo = dir(s.path);
% dirinfo(~[dirinfo.isdir]) = [];  %remove non-directories
% dirinfo = dirinfo(3:end);  %remove non-directories
% 
% subdirinfo = cell(length(dirinfo),1);
% for K = 1 : length(dirinfo)
%   thisdir = fullfile(s.path,dirinfo(K).name);
%   
%   subdirinfo{K} = dir(fullfile(thisdir, 'MOL_*'));
%   
%   subdirinfo{K} = dir(fullfile(thisdir, 'MOL_*.m'));
% end
% 
% for K = 1 : length(dirinfo)
%     plot_fig_handles_left.plot_types(end+1) = subdirinfo{K}.name
%     
% end
% subdirinfo{K}

%The variables from either sessionData, trialData, spikeData or lfpData
%that can be selected or split upon to plot or analyze a subset of
%sessions/trials/cells/channels etc.
ConditionalVars = {'mousename' 'Date' 'Rec_datetime' 'State' 'PostChangeOptoStart'...
    'trialType' 'vecResponse' 'correctResponse' 'hasvisualchange' 'hasaudiochange'...
    'visualInt' 'audioInt' 'visualOri' 'audioFreq'...
    'visualIntNorm' 'audioIntNorm' 'visualOriChangeNorm' 'audioFreqChangeNorm' 'visualSpeed' 'hasphotostim' 'photostimPow'...
    'visualOriPostNorm' 'audioFreqPostNorm'... %'visualOriChange' 'audioFreqChange' 
    'cell_ID' 'area' 'celltype'};
maxtickboxes    = 8; %Maximum number of tickboxes displayed

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting figure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plot_fig=figure('Name','Analysis GUI',...
    'units','normalized',...
    'Position',[0 0 0.25 1],...
    'Color',[1 1 1],...
    'Resize','off',...
    'NumberTitle','off',...
    'tag','_main');

%% Left-hand section to select settings
plot_fig_handles_left.panel = uipanel('Parent',plot_fig,'units','normalized',...
    'pos',[0.015 0.015 0.975 0.975],'BackgroundColor',[0.8 0.3 0.9]);

plot_fig_handles_left.header = uicontrol('Parent',plot_fig_handles_left.panel,'Style','text','units','normalized',...
    'Position',[0.025 0.9 0.95 0.1],'String','Figure Settings:','fontsize',20,'BackgroundColor',[0.8 0.3 0.9]);

%% General Settings for all figures:

x_pos=0.025;
y_pos=0.9;
x_width=0.4;
y_width=0.03;

%These colors are used for the toggles, fields and boxes:
plot_fig_handles_left.Backgroundcolor1 = [0 1 0];
plot_fig_handles_left.Backgroundcolor2 = [1 1 1];
plot_fig_handles_left.Backgroundcolor3 = [0.8 1 0.7];

%Select Plot type
plot_fig_handles_left.plot_type_header = uicontrol('Parent',plot_fig_handles_left.panel,'Style', 'text','units','normalized',...
    'Position',[x_pos y_pos x_width y_width],'String','ANALYSIS TYPE:',...
    'backgroundcolor',plot_fig_handles_left.Backgroundcolor1);
plot_fig_handles_left.plot_type = uicontrol('Parent',plot_fig_handles_left.panel,'Style','popupmenu','units','normalized',...
    'Position',[x_pos+x_width y_pos x_width y_width],'String',plot_fig_handles_left.plot_types,...
    'backgroundcolor',plot_fig_handles_left.Backgroundcolor2,...
    'Callback','global plot_fig_handles_left; if get(plot_fig_handles_left.plot_type,''Value'')<2; set(plot_fig_handles_left.plotvarpanel,''visible'',''on''); else; set(plot_fig_handles_left.plotvarpanel,''visible'',''off'');end;');
%     'backgroundcolor',plot_fig_handles_left.Backgroundcolor2);

%% Plot Variables panel (displayed when plot type is MOL_PlotVar)
plot_fig_handles_left.plotvarpanel = uipanel('Parent',plot_fig_handles_left.panel,'units','normalized',...
    'pos',[0.025 0.775 0.95 0.125],'BackgroundColor',[0.3 0.6 0.9]);

x_pos=0.0;
y_pos=0.725;
x_width=0.425;
y_width=0.25;

%Select Plot style:
plot_fig_handles_left.plot_styles = {'BAR' 'HIST' 'CUM' 'HAZARD' 'XYDATA'};
plot_fig_handles_left.plot_style_header = uicontrol('Parent',plot_fig_handles_left.plotvarpanel,'Style', 'text','units','normalized',...
    'Position',[x_pos y_pos x_width y_width],'String','PLOT STYLE',...
    'backgroundcolor',plot_fig_handles_left.Backgroundcolor1);
plot_fig_handles_left.plot_style = uicontrol('Parent',plot_fig_handles_left.plotvarpanel,'Style','popupmenu','units','normalized',...
    'Position',[x_pos+x_width y_pos x_width y_width],'String',plot_fig_handles_left.plot_styles,...
    'backgroundcolor',plot_fig_handles_left.Backgroundcolor2); y_pos=y_pos-y_width;

%Select Data type:
plot_fig_handles_left.data_types = fieldnames(Data);
plot_fig_handles_left.data_type_header = uicontrol('Parent',plot_fig_handles_left.plotvarpanel,'Style', 'text','units','normalized',...
    'Position',[x_pos y_pos x_width y_width],'String','DATA TYPE',...
    'backgroundcolor',plot_fig_handles_left.Backgroundcolor1);
plot_fig_handles_left.data_type = uicontrol('Parent',plot_fig_handles_left.plotvarpanel,'Style','popupmenu','units','normalized',...
    'Position',[x_pos+x_width y_pos x_width y_width],'String',plot_fig_handles_left.data_types,...
    'backgroundcolor',plot_fig_handles_left.Backgroundcolor2,'Callback','global plot_fig_handles_left; set(plot_fig_handles_left.YDATA,''Value'',1); set(plot_fig_handles_left.YDATA,''String'',fieldnames(Data.(plot_fig_handles_left.data_types{get(plot_fig_handles_left.data_type,''value'')}))); set(plot_fig_handles_left.XDATA,''Value'',1); set(plot_fig_handles_left.XDATA,''String'',fieldnames(Data.(plot_fig_handles_left.data_types{get(plot_fig_handles_left.data_type,''value'')})))'); y_pos=y_pos-y_width;

% Choose from available variables (Which vars can you plot):
plot_fig_handles_left.YDATAVARS = fieldnames(Data.(plot_fig_handles_left.data_types{get(plot_fig_handles_left.data_type,'value')}));
plot_fig_handles_left.YDATAHEADER = uicontrol('Parent',plot_fig_handles_left.plotvarpanel,'Style', 'text','units','normalized',...
    'Position',[x_pos y_pos x_width y_width],'String','SELECT Y DATA','backgroundcolor',plot_fig_handles_left.Backgroundcolor1);
plot_fig_handles_left.YDATA = uicontrol('Parent',plot_fig_handles_left.plotvarpanel,'Style','popupmenu','units','normalized',...
    'Position',[x_pos+x_width y_pos x_width y_width],'String',plot_fig_handles_left.YDATAVARS,...
    'backgroundcolor',plot_fig_handles_left.Backgroundcolor2);y_pos=y_pos-1*y_width;

plot_fig_handles_left.XDATAVARS = fieldnames(Data.(plot_fig_handles_left.data_types{get(plot_fig_handles_left.data_type,'value')}));
plot_fig_handles_left.XDATAHEADER = uicontrol('Parent',plot_fig_handles_left.plotvarpanel,'Style','text','units','normalized',...
    'Position',[x_pos y_pos x_width y_width],'String','SELECT X AXIS','backgroundcolor',plot_fig_handles_left.Backgroundcolor1);
plot_fig_handles_left.XDATA = uicontrol('Parent',plot_fig_handles_left.plotvarpanel,'Style','popupmenu','units','normalized',...
    'Position',[x_pos+x_width y_pos x_width y_width],'String',plot_fig_handles_left.XDATAVARS,...
    'backgroundcolor',plot_fig_handles_left.Backgroundcolor2); %y_pos=y_pos-1.5*y_width;

%% Extra button to remove last 20 trials for behavioral sessions:
plot_fig_handles_left.lastn_button = uicontrol('Parent',plot_fig_handles_left.panel,'Style','checkbox','units','normalized',...
    'Position',[x_pos+x_width 0.75 0.45 0.025],'String','rem. last 20 trials',...
    'backgroundcolor',plot_fig_handles_left.Backgroundcolor2);

%% Selection panel
%The variables from either sessionData, trialData, spikeData or lfpData
%that can be selected or split upon to plot or analyze a subset of
%sessions/trials/cells/channels etc.

x_pos       = 0.025;
y_pos       = 0.925;
% x_width_box = 0.17;
y_width      = 0.05;
x_width_header = 0.3;

plot_fig_handles_left.Conditional=uipanel('Parent',plot_fig_handles_left.panel,'units','normalized',...
    'pos',[0.025 0.25 0.95 0.5],'BackgroundColor',[0.7 0.9 0.1]);

Datatypes = fieldnames(Data);
for i = 1:length(Datatypes)
    Datatype = Datatypes{i};
    AllVars = fieldnames(Data.(Datatype))';
    
    plot_fig_handles_left.data_type_header = uicontrol('Parent',plot_fig_handles_left.Conditional,'Style', 'text','units','normalized',...
        'Position',[x_pos y_pos 0.95 y_width],'String',sprintf('Select %s variables:',Datatype),...
        'backgroundcolor',plot_fig_handles_left.Backgroundcolor3); y_pos=y_pos-y_width;
    
    if strcmp(Datatype,'spikeData')
        plot_fig_handles_left.filterneuronbutton = uicontrol('Parent',plot_fig_handles_left.Conditional,'Style', 'push','units','normalized',...
        'Position',[0.75 y_pos+y_width 0.2 y_width],'String','Filter neurons',...
        'backgroundcolor',plot_fig_handles_left.Backgroundcolor2,'Callback','global Data; Data.spikeData = MOL_filterNeurons(Data.sessionData,Data.trialData,Data.spikeData);');
    end
    
    for var = 1:length(AllVars)
        if any(strcmp(AllVars{var},ConditionalVars))
            if iscell(Data.(Datatype).(AllVars{var}))
                nonemptyind = ~cellfun(@isempty,(Data.(Datatype).(AllVars{var})));
                uniqueTickBoxes = unique(Data.(Datatype).(AllVars{var})(nonemptyind));
            elseif isnumeric(Data.(Datatype).(AllVars{var}))
                nonnandataind = ~isnan(Data.(Datatype).(AllVars{var}));
                uniqueTickBoxes = unique(Data.(Datatype).(AllVars{var})(nonnandataind));
            else error('Unknown field format')
            end
            
            if numel(uniqueTickBoxes)>0 %(if any unique conditions are present)
                plot_fig_handles_left.Toggle.(Datatype).(AllVars{var}) = uicontrol('Parent',plot_fig_handles_left.Conditional,'Style','toggle','units','normalized',...
                    'Position',[x_pos y_pos x_width_header y_width],'String',AllVars{var},'backgroundcolor',plot_fig_handles_left.Backgroundcolor1);
                
                if strcmp(AllVars{var},'cell_ID')
                    plot_fig_handles_left.Toggle.(Datatype).(AllVars{var}).Callback = 'MOL_GUI_CellSelector()';
                end
                
                for box = 1:length(uniqueTickBoxes) %(make tickbox for all unique values)
                    if iscell(uniqueTickBoxes(box))
                        checkboxstring = uniqueTickBoxes{box};
                    elseif isnumeric(uniqueTickBoxes(box))
                        checkboxstring = num2str(uniqueTickBoxes(box));
                    end
                    
                    x_width_box = 0.07 + 0.017*numel(checkboxstring);
                    if x_width_box >0.15
                        x_width_box = 0.15;
                    end
                    if x_pos+x_width_header+x_width_box > 1 && length(uniqueTickBoxes)<=maxtickboxes
                        x_pos = 0.025; y_pos=y_pos-y_width;
                    end
                    plot_fig_handles_left.Checkboxes.(Datatype).(AllVars{var})(box) = uicontrol('Parent',plot_fig_handles_left.Conditional,'Style','checkbox','units','normalized',...
                        'Position',[x_pos+x_width_header y_pos x_width_box y_width],'String',checkboxstring,...
                        'backgroundcolor',plot_fig_handles_left.Backgroundcolor2,'value',1);x_pos=x_pos+x_width_box;
                    if length(uniqueTickBoxes)>maxtickboxes
                        set(plot_fig_handles_left.Checkboxes.(Datatype).(AllVars{var})(box),'Visible','off')
                    end
                end
                y_pos=y_pos-y_width;
                x_pos = 0.025;
            end
            %         YVars = YVars(~strcmp(YVars,AllVars{var}));
            %         XVars = XVars(~strcmp(XVars,AllVars{var}));
        end
    end
end

%% Lower Buttons:

% Button to plot with current settings:
plot_fig_handles_left.plotbutton = uicontrol('Parent',plot_fig_handles_left.panel,'Style','push','units','normalized',...
    'pos',[0.025 0.125 0.45 0.1],'String','PLOT','fontsize',14,...
    'backgroundcolor',[1,0.8,0.6],'Callback','MOL_GUI_PlotFigure');

% Button to save current figure:
plot_fig_handles_left.savebutton = uicontrol('Parent',plot_fig_handles_left.panel,'Style','push','units','normalized',...
    'pos',[0.525 0.125 0.45 0.1],'String','SAVE','fontsize',14,...
    'backgroundcolor',[0.1,0.2,0.9],'Callback','MOL_SaveFig');

% Button to close all figures:
plot_fig_handles_left.closebutton = uicontrol('Parent',plot_fig_handles_left.panel,'Style','push','units','normalized',...
    'pos',[0.025 0.025 0.45 0.1],'String','CLOSE ALL','fontsize',14,...
    'backgroundcolor',[1,0.1,0.9],'Callback','global plot_fig; set(plot_fig, ''HandleVisibility'', ''off''); close all; set(plot_fig, ''HandleVisibility'', ''on'');');

% Option to reload new data (discard old)
plot_fig_handles_left.loadbutton = uicontrol('Parent',plot_fig_handles_left.panel,'Style', 'push','units','normalized',...
    'Position',[0.525 0.025 0.45 0.1],'String','RELOAD DATA','fontsize',14,...
    'backgroundcolor',[0.7,0.5,0.2],'Callback','Data = MOL_SelectData(); global Data; MOL_GUI_Build;');

end



