function MOL_SnakePlot(varargin)
%% Get input arguments:
sessionData     = varargin{1};
trialData       = varargin{2};
spikeData       = varargin{3};

%% SnakePlot Parameters:
% snkprm.SortBy           = 'peakLatency';
snkprm.SortBy           = 'maxResponse';
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
snkprm.TrainFraction    = 0;              %Fraction of trials used for training order, rest for testing
% snkprm.cscale           = [0 1]; 
snkprm.cscale           = [-1 2]; 

%% Parameter settings for PSTH
params                      = params_histresponse(); % All time is in microseconds
params.histmethod           = 'individual'; %Whether to bin (and smooth) individual trials 'individual' or or all trials together ('total')
params.zscore               = 1; 

if strfind(snkprm.AlignOn,'lick')
    params.t_pre                = params.t_pre-1e6;
    params.twin_baseline_start  = params.twin_baseline_start- 1e6;         %Start window for calculating baseline
    params.twin_baseline_stop   = params.twin_baseline_stop- 1e6;       %End window baseline
end

% [splits,colors]     = MOL_GUI_GetSplits(trialData);             %Get splits on the basis of the GUI

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

        %% Take random subset of events_ts for sorting and for display:
        idx_train   = rand(length(events_ts),1)<snkprm.TrainFraction;
        idx_test    = ~idx_train;
        
        %% Get histogram per cell
        [edges,hist_mat_train]  = calc_psth(events_ts(idx_train),spikes_ts,params);
        [edges,hist_mat_test]   = calc_psth(events_ts(idx_test),spikes_ts,params);
        
        if snkprm.AlignRespTime
            hist_mat_test = MOL_WarpPsth(hist_mat_test,edges,0,2e6,temptrialData.responseLatency,nanmedian(trialData.responseLatency));
            hist_mat_train = MOL_WarpPsth(hist_mat_train,edges,0,2e6,temptrialData.responseLatency,nanmedian(trialData.responseLatency));
        end
        
        %% Take the mean of all responses:
        hist_mean_train                 = mean(hist_mat_train,1);
        hist_mean_test                  = mean(hist_mat_test,1);
        
        %% Normalize:
        if snkprm.Normalize
            hist_mean_train                  = hist_mean_train/max(hist_mean_train);
            hist_mean_test                   = hist_mean_test/max(hist_mean_test);
        end
        
        snakemat_train(neuroncounter,:)     = hist_mean_train; %Store in matrix overview
        snakemat_test(neuroncounter,:)      = hist_mean_test; %Store in matrix overview
        neuroncounter                       = neuroncounter+1; %add 1 to neuroncounter
        if snkprm.TrainFraction && any(isnan(snakemat_train(neuroncounter-1,:)))
            neuroncounter = neuroncounter-1;
        elseif any(isnan(snakemat_test(neuroncounter-1,:)))
            neuroncounter = neuroncounter-1;
        end
    end
end

%% Sort snakemat
if isfield(snkprm,'SortBy')
    if snkprm.TrainFraction
        [maxresp,maxidx]          = max(snakemat_train(:,edges>params.twin_resp_start & edges<=params.twin_resp_stop),[],2);
    else
        [maxresp,maxidx]          = max(snakemat_test(:,edges>params.twin_resp_start & edges<=params.twin_resp_stop),[],2);
    end
    %         [maxresp,maxidx]          = max(snakemat_train,[],2);
    if strcmp(snkprm.SortBy,'peakLatency')
        [~,sortidx]             = sort(maxidx,1,'descend');
        snakemat_test            = snakemat_test(sortidx,:);
    elseif strcmp(snkprm.SortBy,'maxResponse')
        [~,sortidx]         = sort(maxresp,1,'descend');
        snakemat_test            = snakemat_test(sortidx,:);
    end
end

snakemat_test = snakemat_test(~any(snakemat_test>10,2),:);

%% Make figure:
figure; set(gcf,'units','normalized','Position',[0.1 0.1 0.5 0.7],'color','w')
imagesc(snakemat_test,snkprm.cscale); hold on;
plot([find(edges == 0) find(edges == 0)], [0 size(snakemat_test,1)+0.5],'k','LineWidth',5);
set(gca, 'XTick', 1:500:length(edges), 'XTickLabels', edges(1:500:length(edges))/1e6,'FontSize', 20)
set(gca, 'YTick', [1 size(snakemat_test,1)], 'YTickLabels', [1 size(snakemat_test,1)],'FontSize', 20)
xlabel('Time (s)','FontSize', 20)
ylabel('Neuron','FontSize', 20)
c = colorbar;
colormap(jet)
c.Label.String = 'Z-scored firing rate';

%% Make figure of the mean:
figure; set(gcf,'units','normalized','Position',[0.4 0.3 0.4 0.3],'color','w')
shadedErrorBar(edges,nanmean(snakemat_test,1),nanstd(snakemat_test,1)/sqrt(size(snakemat_test,1)),{'r','markerfacecolor','r'});    
set(gca, 'XTick', edges(1:500:length(edges)), 'XTickLabels', edges(1:500:length(edges))/1e6,'FontSize', 20)
xlim([edges(1) edges(end)]);
ylabel('Z-scored firing rate','FontSize', 20)
xlabel('Time (s)','FontSize', 20)


end