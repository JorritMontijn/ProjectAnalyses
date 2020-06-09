function [] = MOL_Plot_Histmat(edges,hist_mat,splitidxs,colors,params,figtitle)
%% MOL 2018
% plots several psths split based on variable splitidxs
% splitidxs is a cell array (1xN) with each cell a logical column vector of the number
% of trials (same length hist_mat) and N is the number of different splits

%% Create figure PSTH:
figure; hold all;
set(gcf,'units','normalized','Position',[0.2 0.1 0.5 0.7],'color','w')
% title(cell_ID,'FontSize',25)
title(sprintf('%s',figtitle))
for sp = 1:size(splitidxs,2)
    %% Get the mean, std and sem:
    hist_mat_mean   = mean(hist_mat(splitidxs{sp},:),1);
    hist_mat_std    = std(hist_mat(splitidxs{sp},:),1);
    hist_mat_sem    = hist_mat_std/sqrt(size(hist_mat(splitidxs{sp},:),1));
    
    switch params.plottype
        case 'bins'
            bar(edges/1e6,hist_mat_mean);
        case 'errorline'   %hist as line with error margins:
            CurrLineHandle = shadedErrorBar(edges,hist_mat_mean,hist_mat_sem,{'-','LineWidth',1.5,'Color',colors{sp}},0);
            %                     shadedErrorBar(edges,hist_mat_mean,hist_mat_sem,{'-','LineWidth',1.5,'Color',colors{sp}},1);
            LineHandles(sp) = CurrLineHandle.mainLine;
    end
end

%Figure makeup:
xlabel('Time (sec)','FontSize',15)    %Label x-axis
ylabel('Hz','FontSize',15)            %Label y-axis
set(gca, 'XTick', edges(1:500:length(edges)), 'XTickLabels', edges(1:500:length(edges))/1e6,'FontSize', 20)
set(gca, 'FontSize', 20)
xlim([edges(1) edges(end)]);
% ylim([0 18]);
% legend(LineHandles,{'Imperceptible' 'Subthreshold' 'Threshold' 'Suprathreshold' 'Full'},'FontSize',15);
legend(LineHandles,{'Conflict' 'Visual' 'Auditory'},'FontSize',15);
% legend(LineHandles,{'Miss' 'Hit'},'FontSize',15);

legend(gca,'boxoff');

end