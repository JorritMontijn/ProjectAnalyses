function MOL_GUI_PlotFigure

global plot_fig_handles_left

script = plot_fig_handles_left.plot_types{get(plot_fig_handles_left.plot_type,'value')};

%% Get the data and call the right plotting function:
if strcmp(script,'MOL_PlotVar')
    [plotydata,plotxaxisdata,colors,labels,errorflag] = MOL_GUI_PrepareDataPlot();

    if ~errorflag %If errorflag was set, no plotting will occur.. too bad
        Plot_Style = plot_fig_handles_left.plot_styles{get(plot_fig_handles_left.plot_style,'value')};
        switch Plot_Style
            case 'BAR'
                MOL_PlotBar(plotydata,colors,labels);
            otherwise
                MOL_PlotHist(plotydata,plotxaxisdata,colors,labels,Plot_Style);
        end
    end
else
    [sessionData,trialData,spikeData,errorflag] = MOL_GUI_PrepareDataFunc(); %#ok<ASGLU>
    % Execute the right script function:
    if ~errorflag
        eval(sprintf('%s(sessionData,trialData,spikeData);',script));
    end
end


end