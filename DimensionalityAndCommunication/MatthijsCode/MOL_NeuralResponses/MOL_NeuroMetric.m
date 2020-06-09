function MOL_NeuroMetric(varargin)
%% Get input arguments:
sessionData     = varargin{1};
trialData       = varargin{2};
spikeData       = varargin{3};

%% Parameters for neuronal responses:
params = params_histresponse(); % All time is in microseconds

%% Parameters for neurometric performance:
params.eventofinterest      = 'stimStart';
params.tuningofinterest     = 'visualInt';
% params.tuningofinterest   = 'audioInt';
params.allowsplitdata       = 1;
params.normresponses        = 1;
params.mintrials            = 4; %minimum amount of trials per condition to include

%% Get responses per trial for each neuron
neuroncounter = 1;

for sesid = unique(sessionData.session_ID)'
    %% Get the relevant data for each session individually:
    [temptrialData,tempspikeData] = MOL_getTempPerSes(sesid,trialData,spikeData);
    
    if params.allowsplitdata
        [splitidxs,colors]     = MOL_GUI_GetSplits(temptrialData);                 %Get splits on the basis of the GUI
    else  colors{1} = 'k';
    end
    
    allTuningConds          = unique(temptrialData.(params.tuningofinterest)(~isnan(temptrialData.(params.tuningofinterest))));
    for sp = 1:size(splitidxs,2) %Preallocate
        resp_perTuningCond{sp}      = NaN(length(spikeData.ts),length(allTuningConds));
%         sem_perTuningCond{sp}       = NaN(length(spikeData.ts),length(allTuningConds));
    end

    %% Loop over neurons
    for iNeuron = 1:length(tempspikeData.ts)
        events_ts = temptrialData.(params.eventofinterest)(temptrialData.session_ID == tempspikeData.session_ID(iNeuron));
        spikes_ts = tempspikeData.ts{iNeuron};
        
        %Get histogram per neuron
        [edges,hist_mat] = calc_psth(events_ts,spikes_ts,params);
        %Get responses per neuron
        [resp,~]            = calc_resp_from_psth(edges,hist_mat,params);

        %Get tuning per neuron:
        for sp = 1:size(splitidxs,2)
            for iTuningCond = 1:length(allTuningConds)
                idx_trials = temptrialData.(params.tuningofinterest)==allTuningConds(iTuningCond) & splitidxs{sp};
                if sum(idx_trials)>params.mintrials
                    resp_perTuningCond{sp}(neuroncounter,iTuningCond) = mean(resp(idx_trials));
                    sem_perTuningCond{sp}(neuroncounter,iTuningCond) = std(resp(idx_trials))/sqrt(sum(idx_trials));
                end
            end
        end
        neuroncounter                       = neuroncounter+1; %add 1 to neuroncounter
    end
    
end

%% Normalize if requested
if params.normresponses
   for sp = 1:size(splitidxs,2)
       for iNeuron = 1:size(resp_perTuningCond{sp},1)
           resp_perTuningCond{sp}(iNeuron,:) = resp_perTuningCond{sp}(iNeuron,:)/max(resp_perTuningCond{sp}(iNeuron,:));
       end
   end
end

%% Make the figure
figure; set(gcf,'units','normalized','Position',[0.2 0.1 0.5 0.7],'color','w')
for sp = 1:size(splitidxs,2)
    meanresp{sp} = nanmean(resp_perTuningCond{sp},1);
    sem_resp{sp} = nanstd(resp_perTuningCond{sp},1)/sqrt(size(resp_perTuningCond{sp},1));
    errorbar(allTuningConds,meanresp{sp},sem_resp{sp},'Color',colors{sp},'LineWidth',3); hold on;
end
%pimping figure:
title(sprintf('%s tuning - Average Normalized response',params.tuningofinterest))
set(gca,'XTick',allTuningConds,'XTickLabels',allTuningConds*100,'FontSize', 15)
set(gca, 'XScale', 'log')
xlim([min(allTuningConds) max(allTuningConds)])
ylim([0.1 1]);
ylabel('Response','Fontsize',20);
xlabel(params.tuningofinterest,'Fontsize',20);


        
end
