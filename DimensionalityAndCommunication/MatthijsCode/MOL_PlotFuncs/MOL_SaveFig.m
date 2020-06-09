function MOL_SaveFig()

Par.SaveAs = {'eps' 'bmp'}; %'bmp' 'svg' 'pdf' 'jpg' 'eps'

%File folder as given by user input
% [filename,pathname] = uiputfile('E:\Documents\PhD\Figures',{'*.m';'*.mdl';'*.mat';'*.*'},'Give file name:');
[filename,pathname] = uiputfile('E:\Documents\PhD\Figures\*.*','Give file name:');
[~,filename,~] = fileparts(filename); %Remove extension

%% Save figure as:
if filename
    figHandles = get(groot,'Children'); %Get all figures
    figure(figHandles(2))

    if any(strcmp(Par.SaveAs,'pdf'))
        print(figure(figHandles(2)),fullfile(pathname,strcat(filename,'.pdf')),'-dpdf','-bestfit');
        export_fig(fullfile(pathname,filename),'-dpdf');
    end
    
    if any(strcmp(Par.SaveAs,'eps'))
        export_fig(fullfile(pathname,filename),'-eps');
    end
    
    if any(strcmp(Par.SaveAs,'jpg'))
        print(figure(figHandles(2)),fullfile(pathname,strcat(filename,'.jpg')),'-djpg');
    end
    
    if any(strcmp(Par.SaveAs,'svg'))
      saveas(figHandles(2),fullfile(pathname,strcat(filename,'.svg')), 'svg')
    end
    
    if any(strcmp(Par.SaveAs,'bmp'))
      saveas(figHandles(2),fullfile(pathname,strcat(filename,'.bmp')), 'bmp')
    end
    
    if any(strcmp(Par.SaveAs,'png'))
        export_fig(fullfile(pathname,filename),'-png');
    end
end

%% Print info about the picture onto the pdf?


end