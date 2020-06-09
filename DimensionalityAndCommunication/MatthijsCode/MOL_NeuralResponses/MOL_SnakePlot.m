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
snkprm.clipzscore       = 15;
snkprm.minTrialCond     = 5;
snkprm.Normalize        = 0;                %Whether firing rate are normalized to maximum for that neuron (each scaled between 0-1)
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

[splits,colors]     = MOL_GUI_GetSplits(trialData);             %Get splits on the basis of the GUI

%% Remove sessions with too few trials in a give condition:
uniqueses = unique(sessionData.session_ID);
for iSes = 1:length(uniqueses)
    for iSplit = 1:length(splits)
        nTrials(iSes,iSplit) = sum(strcmp(trialData.session_ID,uniqueses(iSes)) & splits{iSplit});
    end
end
if any(any(nTrials<snkprm.minTrialCond,2))
    fprintf('Removed sessions due to insufficient trials in one split conditions\n');
    for iSplit = 1:length(splits)
        splits{iSplit} = splits{iSplit}(ismember(trialData.session_ID,uniqueses(~any(nTrials<snkprm.minTrialCond,2))));
    end
    [sessionData,trialData,spikeData] = MOL_getTempPerSes(uniqueses(~any(nTrials<snkprm.minTrialCond,2)),sessionData,trialData,spikeData);
end

%% Main loop to get psth matrix:
for iSplit = 1:length(splits)
    
    datafields = fieldnames(trialData);
    splittrialData = struct();
    for field = 1:length(datafields)
        splittrialData.(datafields{field}) = trialData.(datafields{field})(splits{iSplit});
    end
    
    neuroncounter               = 1;
    snakemat                    = [];
    for sesid = unique(sessionData.session_ID)'
        %% Get the relevant data for each session individually:
        [temptrialData,tempspikeData] = MOL_getTempPerSes(sesid,splittrialData,spikeData);
        
        if strfind(snkprm.AlignOn,'lick')
            events_ts                   = MOL_getLicks(temptrialData,snkprm.AlignOn);
        else
            events_ts = temptrialData.(snkprm.AlignOn);
            if isempty(events_ts)
                error('no events found to align to')
            end
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
            
            %% Take the mean of all responses:
            hist_mean                       = nanmean(hist_mat,1);
            
            %% Normalize:
            if snkprm.Normalize
                hist_mean                   = hist_mean/max(hist_mean);
            end
            
            snakemat(neuroncounter,:)       = hist_mean; %Store in matrix overview
            neuroncounter                   = neuroncounter+1; %add 1 to neuroncounter
            %         if snkprm.TrainFraction && any(isnan(snakemat(neuroncounter-1,:)))
            %             neuroncounter = neuroncounter-1;
            %         elseif any(isnan(snakemat_test(neuroncounter-1,:)))
            %             neuroncounter = neuroncounter-1;
            %         end
        end
    end
    snakemat_splits{iSplit} = snakemat;
end

%% Sort snakemat
if isfield(snkprm,'SortBy')
%     [maxresp,maxidx]          = max(snakemat_splits{1}(:,edges>params.twin_resp_start & edges<=params.twin_resp_stop),[],2);
    [maxresp,maxidx]          = max(snakemat_splits{1}(:,edges>params.twin_resp_start & edges<=params.twin_resp_stop),[],2);
    %         [maxresp,maxidx]          = max(snakemat_train,[],2);
    if strcmp(snkprm.SortBy,'peakLatency')
        [~,sortidx]             = sort(maxidx,1,'descend');
        for iSplit = 1:length(splits)
            snakemat_splits{iSplit}            = snakemat_splits{iSplit}(sortidx,:);
        end
    elseif strcmp(snkprm.SortBy,'maxResponse')
        [~,sortidx]             = sort(maxresp,1,'descend');
        for iSplit = 1:length(splits)
            snakemat_splits{iSplit}            = snakemat_splits{iSplit}(sortidx,:);
        end
    end
end

if length(splits)==1
    snakemat_splits{1} = snakemat_splits{1}(~any(snakemat_splits{1}>snkprm.clipzscore,2),:);
else
    for iSplit = 1:length(splits)
        snakemat_splits{1} = snakemat_splits{1}(~any(snakemat_splits{1}>snkprm.clipzscore,2),:);
        snakemat_splits{2} = snakemat_splits{2}(~any(snakemat_splits{1}>snkprm.clipzscore,2),:);
        
        snakemat_splits{1} = snakemat_splits{1}(~any(snakemat_splits{2}>snkprm.clipzscore,2),:);
        snakemat_splits{2} = snakemat_splits{2}(~any(snakemat_splits{2}>snkprm.clipzscore,2),:);
    end
end

%% Make figure:
figure; set(gcf,'units','normalized','Position',[0.1 0.1 0.6 0.5],'color','w')
for iSplit = 1:length(splits)
    subplot(1,length(splits),iSplit)
    imagesc(snakemat_splits{iSplit},snkprm.cscale); hold on;
    plot([find(edges == 0) find(edges == 0)], [0 size(snakemat_splits{iSplit},1)+0.5],'k','LineWidth',5);
    set(gca, 'XTick', 1:1000:length(edges), 'XTickLabels', edges(1:1000:length(edges))/1e6,'FontSize', 20)
    set(gca, 'YTick', [1 size(snakemat_splits{iSplit},1)], 'YTickLabels', [1 size(snakemat_splits{iSplit},1)],'FontSize', 20)
    xlabel('Time (s)','FontSize', 20)
    ylabel('Neuron','FontSize', 20)
    c = colorbar;
    colormap(parula);
    c.Label.String = 'Z-scored firing rate';
end

%% Make figure of the mean:
figure; set(gcf,'units','normalized','Position',[0.4 0.3 0.4 0.3],'color','w')
for iSplit = 1:length(splits)
    color_char = colornames('MATLAB',colors{iSplit});
%     shadedErrorBar(edges,nanmean(snakemat_splits{iSplit},1),nanstd(snakemat_splits{iSplit},1)/sqrt(size(snakemat_splits{iSplit},1)),{color_char{1}(1),'markerfacecolor',colors{iSplit}});
    shadedErrorBar(edges,nanmean(snakemat_splits{iSplit},1),nanstd(snakemat_splits{iSplit},1)/sqrt(size(snakemat_splits{iSplit},1)),{color_char{1}(1),'markerfacecolor',colors{iSplit}},0);
%     shadedErrorBar(edges,nanmean(snakemat_splits{iSplit},1),nanstd(snakemat_splits{iSplit},1)/sqrt(size(snakemat_splits{iSplit},1)));
    hold all;
end
set(gca, 'XTick', edges(1:500:length(edges)), 'XTickLabels', edges(1:500:length(edges))/1e6,'FontSize', 20)
xlim([edges(1) edges(end)]);
ylabel('Z-scored firing rate','FontSize', 20)
xlabel('Time (s)','FontSize', 20)


end