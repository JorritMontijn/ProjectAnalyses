function MOL_PlotTuning(varargin)
%% Get input arguments:
sessionData     = varargin{1};
trialData       = varargin{2};
spikeData       = varargin{3};

%% Get parameters for hist and responses:
if strfind(sessionData.Experiment{1},'Bars')
    params                  = params_histresponse_pmal(); % All time is in microseconds
elseif strfind(sessionData.Experiment{1},'Gratings')
    params                  = params_histresponse_pmal_gratings(); % All time is in microseconds
else params                 = params_histresponse(); % All time is in microseconds
end

%% Extra parameters:
if strfind(sessionData.Experiment{1},'ChangeDetection')
    params.eventofinterest = 'stimChange';
else
    params.eventofinterest = 'stimStart';
end
% params.eventofinterest      = 'stimStart';
% params.tuningofinterest     = 'visualOri';
params.tuningofinterest   = 'visualOriPostNorm';
% params.tuningofinterest   = 'visualOriChangeNorm';

params.tuningofinterest   = 'audioFreqPostNorm';
% params.tuningofinterest   = 'audioInt';
params.allowsplitdata       = 1;
params.plottype             = 'linear';
% params.plottype             = 'circular';

for sesid = unique(sessionData.session_ID)'
    %% Get the relevant data for each session individually:
    [temptrialData,tempspikeData] = MOL_getTempPerSes(sesid,trialData,spikeData);
    
    if params.allowsplitdata
        [splitidxs,colors]     = MOL_GUI_GetSplits(temptrialData);                 %Get splits on the basis of the GUI
    else  colors{1} = 'k';
    end
    
    %% Loop over neurons
    for iNeuron = 1:length(tempspikeData.ts)
        events_ts = temptrialData.(params.eventofinterest)(strcmp(temptrialData.session_ID,tempspikeData.session_ID(iNeuron)));
        spikes_ts = tempspikeData.ts{iNeuron};
        
        %% Get spiking times for raster plot:
        trial_ts            = cell(1,length(events_ts));
        for ev=1:length(events_ts) % Get spikes within window around event ev:
            trial_ts{ev}     = spikes_ts(spikes_ts>events_ts(ev)+params.t_pre & spikes_ts<=events_ts(ev)+params.t_post)-events_ts(ev);
        end
        
        %Get histogram per neuron
        [edges,hist_mat]    = calc_psth(events_ts,spikes_ts,params);
        
        switch params.plottype
            case 'linear'
                figure; set(gcf,'units','normalized','Position',[0.2 0.1 0.5 0.7],'color','w')
                
                %Compute response per tuning condition per neuron:
                allTuningConds          = unique(temptrialData.(params.tuningofinterest)(~isnan(temptrialData.(params.tuningofinterest))));
                resp_perTuningCond      = NaN(length(allTuningConds),1); sem_perTuningCond = NaN(length(allTuningConds),1);
                for sp = 1:size(splitidxs,2)
                    for iTuningCond = 1:length(allTuningConds)
                        %Get selection of trials for this condition:
                        idx_trials = temptrialData.(params.tuningofinterest)==allTuningConds(iTuningCond) & splitidxs{sp};
                        if strcmp(sessionData.Experiment(strcmp(sessionData.session_ID,sesid)),'Bars') %Get correct time windows:
                            params.twin_resp_stop_trial = temptrialData.stimEnd(idx_trials) - temptrialData.stimStart(idx_trials);  %different response window per trial (variable stimulus length)
                        end
                        %Get responses:
                        [resp,~]                        = calc_resp_from_psth(edges,hist_mat(idx_trials,:),params);
                        resp_perTuningCond(iTuningCond) = mean(resp);
                        sem_perTuningCond(iTuningCond)  = std(resp)/sqrt(sum(idx_trials));
                    end
                    errorbar(allTuningConds,resp_perTuningCond,sem_perTuningCond,'Color',colors{sp},'LineWidth',3); hold on;
                end
                
                suptitle(sprintf('%s tuning - Neuron %d',params.tuningofinterest,tempspikeData.cell_ID{iNeuron}))
                %         set(gca,'XTick',allTuningConds,'XTickLabels',allTuningConds*100,'FontSize', 15)
                %         set(gca, 'XScale', 'log')
                
                set(gca,'XTick',allTuningConds,'XTickLabels',allTuningConds,'FontSize', 15)
%                 xlim([min(allTuningConds) max(allTuningConds)])
                xlim([min(allTuningConds)-10 max(allTuningConds)+10])
                ylim([0 max(resp_perTuningCond)*1.2+0.01]);
                ylabel('Response (Hz)','Fontsize',20);
                xlabel(params.tuningofinterest,'Fontsize',20);
                
            case 'circular'
                circfig = figure; set(gcf,'units','normalized','Position',[0.2 0.1 0.5 0.7],'color','w')
                plrwidth = 0.25;
                plrplot = polaraxes('Position', [0.52-plrwidth/2 0.52-plrwidth/2 plrwidth plrwidth]);
                subplotorder = [6 3 2 1 4 7 8 9];

                %Compute response per tuning condition per neuron:
                allTuningConds          = unique(temptrialData.(params.tuningofinterest)(~isnan(temptrialData.(params.tuningofinterest))));
                resp_perTuningCond      = NaN(length(allTuningConds),1); sem_perTuningCond = NaN(length(allTuningConds),1);
                for sp = 1:size(splitidxs,2)
                    for iTuningCond = 1:length(allTuningConds)
                        splots(iTuningCond) = subplot(3,3,subplotorder(iTuningCond));
                        idx_trials = temptrialData.(params.tuningofinterest)==allTuningConds(iTuningCond) & splitidxs{sp};
                        plot(edges,mean(hist_mat(idx_trials,:),1),'Color',colors{sp},'LineWidth',2); hold on;
                        
                        if strcmp(sessionData.Experiment(sessionData.session_ID==sesid),'Bars') %Get correct time windows:
                            params.twin_resp_stop_trial = temptrialData.stimEnd(idx_trials) - temptrialData.stimStart(idx_trials);  %different response window per trial (variable stimulus length)
                        end
                        %Get responses:
                        [resp,~]                        = calc_resp_from_psth(edges,hist_mat(idx_trials,:),params);
                        resp_perTuningCond(iTuningCond) = mean(resp);
                        sem_perTuningCond(iTuningCond)  = std(resp)/sqrt(sum(idx_trials));
                        
                        title(allTuningConds(iTuningCond))
                        xlim([min(edges) max(edges)])
                        
                        set(gca,'XTick',[],'XTickLabels',[],'FontSize', 15)
                        if ismember(subplotorder(iTuningCond),[1 4 7])
                            ylabel('Firing rate (Hz)','Fontsize',15);
                        else
                            set(gca,'YTick',[],'YTickLabels',[],'FontSize', 15)
                        end
                        set(gca,'XTick',edges(1:1000:end),'XTickLabels',edges(1:1000:end)*1e-6,'FontSize', 15)
                        if ismember(subplotorder(iTuningCond),[7])
                            xlabel('Time (s)','Fontsize',15);
                        end
                    end
                    axes(plrplot); hold on;
                    polarplot(plrplot,deg2rad([allTuningConds; allTuningConds(1)]),[resp_perTuningCond; resp_perTuningCond(1)],'Color',colors{sp},'LineWidth',3); hold on;
                end
                
                for iTuningCond = 1:length(allTuningConds)
                    subplot(3,3,subplotorder(iTuningCond));
                    ylim([0 max([splots.YLim])/1.2+0.01]);
                end
                suptitle(sprintf('%s tuning - Neuron %d',params.tuningofinterest,tempspikeData.cell_ID(iNeuron)))

%                 suptitle(sprintf('%s tuning - Neuron %d',params.tuningofinterest,tempspikeData.cell_ID(iNeuron)))
                %         set(gca,'XTick',allTuningConds,'XTickLabels',allTuningConds*100,'FontSize', 15)
                %         set(gca, 'XScale', 'log')
                
        end
        
    end
    
end
distFig('Not',1) %distribute figures over the screen
end
