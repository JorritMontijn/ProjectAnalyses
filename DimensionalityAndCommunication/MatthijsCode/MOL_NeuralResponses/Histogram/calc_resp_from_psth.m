function [resp,tmax] = calc_resp_from_psth(edges,hist_mat,params)
% input:
% - edges of size N (time series corresponding to the bins of firing rates)
% - psth matrix of size (M by N) where M is number of events, and N is
% number of edges
% - parameters (Parameter settings are to be defined elsewhere)
% output:
% - responses (of size M)
% - latency of response (of size M)
% MOL 2017

%% Switch between two methods: If the response is taken from the mean of all input trials: average the psth, otherwise calculater per trial
switch params.trialrespmethod
    case 'total'
        hist_mat = repmat(mean(hist_mat,1),size(hist_mat,1),1);
    case 'individual'
        %do nothing
end

%% Calculate response:
resp            = NaN(size(hist_mat,1),1); %Preallocate
tmax            = NaN(size(hist_mat,1),1); %Preallocate
baseline        = NaN(size(hist_mat,1),1); %Preallocate
maxresp_stim    = NaN(size(hist_mat,1),1); %Preallocate

for ev = 1:size(hist_mat,1)
    if isfield(params,'twin_resp_stop_trial') %different response window per trial (variable stimulus length)
        params.twin_resp_stop = params.twin_resp_stop_trial(ev);
    end
    
    %Calculate baseline:
    baseline(ev)                = mean(hist_mat(ev,edges>params.twin_baseline_start & edges<=params.twin_baseline_stop));
    %Get maximum response:
    maxresp_stim(ev)            = max(hist_mat(ev,edges>params.twin_resp_start & edges<=params.twin_resp_stop));
    %Get timing of maximum response:
    tmax(ev)                    = mean(edges(hist_mat(ev,:)==maxresp_stim(ev) & edges>params.twin_resp_start & edges<=params.twin_resp_stop));
    
    switch params.respcalcmethod
        case 'max' %Calculate response as max minus baseline:
            resp(ev)    = maxresp_stim(ev);
%             resp(ev)    = maxresp_stim(ev) - baseline(ev);
        case 'mean' %Calculate response as the mean during response window:
            resp(ev)    = mean(hist_mat(ev,edges>params.twin_resp_start & edges<=params.twin_resp_stop));
        case 'meanpeak' %Calculate response as the mean around response peak:
            maxpeak     = max(mean(hist_mat(:,edges>params.twin_resp_start & edges<=params.twin_resp_stop),1));
            t_peak      = mean(edges(mean(hist_mat,1)==maxpeak & edges>params.twin_resp_start & edges<=params.twin_resp_stop));
            resp(ev)    = mean(hist_mat(ev,edges>t_peak-params.peakwin/2 & edges<=t_peak+params.peakwin/2));
        case 'div' %Calculate response as max over baseline:
            resp(ev)    = maxresp_stim(ev) / baseline(ev);
        case 'AUC' %Calculate response as area under curve:
            resp(ev)    = sum(hist_mat(ev,edges>params.twin_resp_start & edges<=params.twin_resp_stop)-baseline(ev)) / params.binsize;
    end
    
    if params.subtr_baseline
        resp(ev)    = resp(ev) - baseline(ev);
    end
    
%     if resp(ev)<0
%         resp(ev) = 0;        
%     end
    
    if isnan(resp(ev))
        resp(ev)    = 0;
        tmax(ev)    = NaN;
        %       %Get negative response:
        %     resp   = min(conv_hist_1ms(edges_1ms>params.twin_resp_start & edges_1ms<=params.twin_resp_stop)) - baseline;
    end
    
end

