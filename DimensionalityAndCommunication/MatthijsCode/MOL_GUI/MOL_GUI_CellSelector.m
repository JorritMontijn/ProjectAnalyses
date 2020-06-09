function MOL_GUI_CellSelector()
global plot_fig_handles_left

if get(plot_fig_handles_left.Toggle.spikeData.cell_ID,'value')
    cell_ids = get(plot_fig_handles_left.Checkboxes.spikeData.cell_ID,'string');
    Cell_ID_selector = figure;
    set(Cell_ID_selector,'units','normalized','Position',[0.2 0.05 0.2 0.85],'color','w')
    plot_fig_handles_left.Cell_ID_Table = uitable('Parent', Cell_ID_selector, 'Data', cell_ids,'units','normalized',...
        'Position', [0.05 0.05 0.9 0.85],'ColumnName','Cell_ID','RowName','numbered',...
        'CellSelectionCallback',@(src,evnt)set(src,'UserData',evnt.Indices));
    set(plot_fig_handles_left.Cell_ID_Table,'ColumnWidth', {135,52,83,83})
    
    Cell_ID_selector_all_cells = uicontrol('Parent',Cell_ID_selector,'Style', 'push','units','normalized',...
        'Position',[0.1 0.9 0.4 0.1],'String','Select ALL Cells','fontsize',14,...
        'backgroundcolor',[0.8,0.3,0.6],'Callback','global plot_fig_handles_left; set(plot_fig_handles_left.Cell_ID_Table,''UserData'',true(length(get(plot_fig_handles_left.Cell_ID_Table,''Data'')),1))');
    
    Cell_ID_selector_load_data = uicontrol('Parent',Cell_ID_selector,'Style', 'push','units','normalized',...
        'Position',[0.5 0.9 0.4 0.1],'String','Select Cells','fontsize',14,...
        'backgroundcolor',[0.2,0.5,0.9],'Callback','uiresume(gcf);');
    uiwait(Cell_ID_selector) %Wait for load button to be pressed: callback is continue with script
    
    idx = get(plot_fig_handles_left.Cell_ID_Table,'UserData'); idx = idx(:,1); %Get selected cells
    for iChBx = 1:length(plot_fig_handles_left.Checkboxes.spikeData.cell_ID)
        set(plot_fig_handles_left.Checkboxes.spikeData.cell_ID(iChBx),'value',ismember(iChBx,idx))
    end
    
    set(plot_fig_handles_left.Toggle.spikeData.cell_ID,'value',0)
    close(Cell_ID_selector);
    
end

end