function MOL_ResponseSnake(varargin)
%% Get input arguments:
sessionData     = varargin{1};
trialData       = varargin{2};
spikeData       = varargin{3};

%% SnakePlot Parameters:
snkprm.SortBy           = 'peakLatency';
% snkprm.SortBy           = 'maxResponse';
% snkprm.SortBy           = 'respLatency';
% snkprm.AlignOn          = 'stimStart';      %On which timestamp to align as t=0
snkprm.AlignOn          = 'stimChange';      %On which timestamp to align as t=0
% snkprm.AlignOn          = 'lick_firstofresponse_R';
% snkprm.AlignOn          = 'lick_firstofresponse_L';
% snkprm.AlignOn          = 'lick_firstofresponse';
% snkprm.AlignOn          = 'lick_firstofitibout_L';
% snkprm.AlignOn          = 'lick_firstofitibout_R';
% snkprm.AlignOn          = 'lick_firstofitibout';

snkprm.AlignRespTime    = 0;

snkprm.Normalize        = 0;                %Whether firing rate are normalized to maximum for that neuron (each scaled between 0-1)
% snkprm.cscale           = [0 1]; 
snkprm.cscale           = [-1 5]; 

%% Parameter settings for PSTH
params                      = params_histresponse(); % All time is in microseconds
params.histmethod           = 'individual'; %Whether to bin (and smooth) individual trials 'individual' or or all trials together ('total')
params.zscore               = 1; 

if strfind(snkprm.AlignOn,'lick')
    params.t_pre                = params.t_pre-1e6;
    params.twin_baseline_start  = params.twin_baseline_start- 1e6;         %Start window for calculating baseline
    params.twin_baseline_stop   = params.twin_baseline_stop- 1e6;       %End window baseline
end

%% Main loop to get psth matrix:
neuroncounter               = 1;
for sesid = unique(sessionData.session_ID)'
    %% Get the relevant data for each session individually:
    [temptrialData,tempspikeData] = MOL_getTempPerSes(sesid,trialData,spikeData);
    
    if strfind(snkprm.AlignOn,'lick')
        events_ts                   = MOL_getLicks(temptrialData,snkprm.AlignOn);
    else
        events_ts = temptrialData.(snkprm.AlignOn);
    end
    
    %% construct normalized psth matrix:
    for iNeuron = 1:length(tempspikeData.ts)
        %% Get the spikes for this neuron:
        spikes_ts = tempspikeData.ts{iNeuron};

        %% Get histogram per cell
        [edges,hist_mat]  = calc_psth(events_ts,spikes_ts,params);
        
        if snkprm.AlignRespTime
            hist_mat = MOL_WarpPsth(hist_mat,edges,0,2e6,temptrialData.responseLatency,nanmedian(trialData.responseLatency));
        end
        
        %% Normalize:
        if snkprm.Normalize
            hist_mat                  = hist_mat/max(mean(hist_mat,1));
        end
        
        %% Sort snakemat
        if isfield(snkprm,'SortBy')
            [maxresp,maxidx]          = max(hist_mat(:,edges>params.twin_resp_start & edges<=params.twin_resp_stop),[],2);
            if strcmp(snkprm.SortBy,'peakLatency')
                [~,sortidx]         = sort(maxidx,1,'descend');
                hist_mat            = hist_mat(sortidx,:);
                resptime            = temptrialData.responseLatency(sortidx,:);
            elseif strcmp(snkprm.SortBy,'respLatency')
                [~,sortidx]         = sort(temptrialData.responseLatency,1,'descend');
                hist_mat            = hist_mat(sortidx,:);
                resptime            = temptrialData.responseLatency(sortidx,:);
            elseif strcmp(snkprm.SortBy,'maxResponse')
                [~,sortidx]         = sort(maxresp,1,'descend');
                hist_mat            = hist_mat(sortidx,:);
                resptime            = temptrialData.responseLatency(sortidx,:);
            end
        else
            resptime            = temptrialData.responseLatency;
        end
        
        %% Make figure:
        figure; set(gcf,'units','normalized','Position',[0.1 0.1 0.5 0.7],'color','w')
        imagesc(hist_mat,snkprm.cscale); hold on;
        suptitle(sprintf('Response Snake - Neuron %16.0f',tempspikeData.cell_ID(iNeuron)))
        plot([find(edges == 0) find(edges == 0)], [0 size(hist_mat,1)+0.5],'k','LineWidth',5);
        set(gca, 'XTick', 1:500:length(edges), 'XTickLabels', edges(1:500:length(edges))/1e6,'FontSize', 20)
        set(gca, 'YTick', [1 size(hist_mat,1)], 'YTickLabels', [1 size(hist_mat,1)],'FontSize', 20)
        xlabel('Time (s)','FontSize', 20)
        ylabel('Trials','FontSize', 20)
        c = colorbar;
        c.Label.String = 'Z-scored firing rate';
        for iTr=1:size(hist_mat,1)
            if ~isnan(resptime(iTr))
                respidx = find(edges>resptime(iTr),1,'first');
                plot([respidx respidx],[iTr-0.5 iTr+0.5],'k','LineWidth',5);
            end            
        end
%         snakemat_train(neuroncounter,:)     = hist_mean; %Store in matrix overview
%         snakemat_test(neuroncounter,:)      = hist_mean_test; %Store in matrix overview
%         neuroncounter                       = neuroncounter+1; %add 1 to neuroncounter
        
    end
end




% %% Make figure of the mean:
% figure; set(gcf,'units','normalized','Position',[0.4 0.3 0.4 0.3],'color','w')
% shadedErrorBar(edges,mean(snakemat_test,1),std(snakemat_test,1)/sqrt(size(snakemat_test,1)),{'r','markerfacecolor','r'});    
% set(gca, 'XTick', edges(1:500:length(edges)), 'XTickLabels', edges(1:500:length(edges))/1e6,'FontSize', 20)
% xlim([edges(1) edges(end)]);
% ylabel('Z-scored firing rate','FontSize', 20)
% xlabel('Time (s)','FontSize', 20)


end