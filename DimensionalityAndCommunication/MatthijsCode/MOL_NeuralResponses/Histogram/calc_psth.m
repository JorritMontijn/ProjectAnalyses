function [edges,hist_mat] = calc_psth(events_ts,spikes_ts,params)
% input:
% - M timestamps for specific events
% - timestamps of all spikes
% - parameters (Parameter settings are to be defined elsewhere)
% output:
% psth matrix of size (M by N) where M is number of events, and N is
% number of edges
% MOL 2017

%% Get spiking times:
trial_ts            = cell(1,length(events_ts));
for ev=1:length(events_ts) % Get spikes within window around event ev:
    trial_ts{ev}     = spikes_ts(spikes_ts>events_ts(ev)+params.t_pre & spikes_ts<=events_ts(ev)+params.t_post)-events_ts(ev);
end

%% Histogram:
edges               = [params.t_pre:params.binsize:params.t_post] - params.binsize/2; %Define edges
% switch params.histmethod
%     case 'total' %Option 1: Take the mean response:
%         hist_mat                = histc([trial_ts{:}],edges)/length(events_ts) * 1e6/params.binsize;
%         hist_mat                = hist_mat(1:end-1); %Remove last value (=exact end)
%     case 'individual' %Option 2: bin individual trials and later the take the mean e.g.:
hist_mat = NaN(length(events_ts),length(edges)-1);
for ev=1:length(events_ts)
    temp                = histc(trial_ts{ev},edges) * 1e6/params.binsize;
    if ~isempty(temp) %No spikes whatsover in the trial
        hist_mat(ev,:)      = temp(1:end-1); %Store + Remove last value (=exactly values at endbin)
    end
end
% end
edges = edges(1:end-1) + params.binsize/2; %Correct edges

%% Apply smoothing if requested:
if params.smoothing
    if strcmp(params.conv_win,'flat') %Construct smoothing window:
        win                 = ones(1,params.conv_twin/params.binsize)/(params.conv_twin/params.binsize);
    elseif strcmp(params.conv_win,'gaussian')
        N                   = params.conv_twin/params.binsize;
        alpha               = ((N-1)/(params.conv_sigma/params.binsize))/2; %? = (N – 1)/(2?)
        win                 = gausswin(N,alpha); %convolution with gaussian
        win                 = win/sum(win); %normalized
    end
    %Smooth either the total or the individual trials:
    %     switch params.histmethod
    %         case 'total'
    %             hist_mat          = padarray(hist_mat,[0 round(length(win)/2)],'symmetric','both'); %pad the array on both sides for convolution
    %             hist_mat          = conv(hist_mat,win,'valid'); %Take only the valid overlapping center of convolution
    %             hist_mat          = hist_mat(1:length(edges)); %slight correction to get same size (edge vs convolution)
    %         case 'individual'
    hist = cell(1,length(events_ts));
    for ev = 1:length(events_ts)
        hist{ev}                = padarray(hist_mat(ev,:),[0 round(length(win)/2)],'symmetric','both'); %pad the array on both sides for convolution
        hist{ev}                = conv(hist{ev},win,'valid'); %Take only the valid overlapping center of convolution
        hist_mat(ev,:)          = hist{ev}(1:length(edges)); %slight correction to get same size (edge vs convolution)
    end
    %     end
end

%% Z-SCORE
if params.zscore
    %Get baseline period:
    bsl_idx = edges>params.twin_baseline_start & edges<=params.twin_baseline_stop;
    %Get std per trial baseline period:
    std_psth_bsl_win = std(hist_mat(:,bsl_idx),0,2);
    std_psth_bsl_win = std_psth_bsl_win(~isnan(std_psth_bsl_win));
    %Each trial minus its mean baseline activity divided by mean standard
    %deviation during baseline all trials:
    if mean(std_psth_bsl_win)==0
        warning('zero spikes found during baseline, no z-score computed, z set to zero')
        hist_mat(:,:) = 0;
    else
        for ev = 1:size(hist_mat,1)
            hist_mat(ev,:) = (hist_mat(ev,:) - mean(hist_mat(ev,bsl_idx)))/mean(std_psth_bsl_win);
        end
    end
end

end
