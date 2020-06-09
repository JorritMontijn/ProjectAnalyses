function MOL_PopulationOriTuning(varargin)
%% Get input arguments:
sessionData     = varargin{1};
trialData       = varargin{2};
spikeData       = varargin{3};

%% Get parameters for hist and responses:
if strfind(sessionData.Experiment{1},'Bars')
    params          = params_histresponse_pmal(); % All time is in microseconds
else    params          = params_histresponse(); % All time is in microseconds
end

%% Extra parameters:
params.eventofinterest      = 'stimStart';
params.tuningofinterest     = 'visualOri';
% params.tuningofinterest   = 'audioInt';
params.allowsplitdata       = 1;
params.mintrials            = 0; %minimum amount of trials per condition to include
params.ShowFig              = 0;
params.fitmethod            = 'Gaussian'; %Gaussian or VonMises
params.FitTuningCurve       = 0; %Gaussian or VonMises

neuroncounter = 1;

for sesid = unique(sessionData.session_ID)'
    %% Get the relevant data for each session individually:
    [tempsessionData,temptrialData,tempspikeData] = MOL_getTempPerSes(sesid,sessionData,trialData,spikeData);
        
    %% Selection of neurons:
    [tempspikeData] = MOL_filterNeurons(tempsessionData,temptrialData,tempspikeData);
    
    %%
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
        
        %Get histogram per neuron:
        events_ts = temptrialData.(params.eventofinterest)(temptrialData.session_ID == tempspikeData.session_ID(iNeuron));
        spikes_ts = tempspikeData.ts{iNeuron};
        
        [edges,hist_mat]    = calc_psth(events_ts,spikes_ts,params);
        
        %Get tuning per neuron:
        for sp = 1:size(splitidxs,2)
            for iTuningCond = 1:length(allTuningConds)
                %Get selection of trials for this condition:
                idx_trials = temptrialData.(params.tuningofinterest)==allTuningConds(iTuningCond) & splitidxs{sp};
                if strcmp(sessionData.Experiment(sessionData.session_ID==sesid),'Bars') %Get correct time windows:
                    params.twin_resp_stop_trial = temptrialData.stimEnd(idx_trials) - temptrialData.stimStart(idx_trials);  %different response window per trial (variable stimulus length)
                end
                %Get responses:
                [resp,~]                        = calc_resp_from_psth(edges,hist_mat(idx_trials,:),params);
                if sum(idx_trials)>params.mintrials
                    resp_perTuningCond{sp}(neuroncounter,iTuningCond) = mean(resp);
                    sem_perTuningCond{sp}(neuroncounter,iTuningCond) = std(resp)/sqrt(sum(idx_trials));
                end
            end
            
            [OSI(sp,neuroncounter),DSI(sp,neuroncounter),gOSI(sp,neuroncounter)] = calc_OSIDSI(allTuningConds',resp_perTuningCond{sp}(neuroncounter,:),0);
            
            if params.FitTuningCurve
                switch params.fitmethod                    %get fit:
                    case 'Gaussian'
                        [Fit(sp,neuroncounter)]      = calc_gaussfit(allTuningConds,resp_perTuningCond{sp}(neuroncounter,:),sem_perTuningCond{sp}(neuroncounter,:),params.ShowFig);
                    case 'VonMises'
                        try
                            [Fit(sp,neuroncounter)]   = calc_VonMisesfit_mean(allTuningConds,resp_perTuningCond{sp}(neuroncounter,:),sem_perTuningCond{sp}(neuroncounter,:),params.ShowFig);
                        catch
                        end
                end
            end
            
        end
        neuroncounter = neuroncounter +1;
    end
end

%% Make the OSI figure
figure; set(gcf,'units','normalized','Position',[0.2 0.1 0.5 0.7],'color','w')
scatter(OSI(1,:),OSI(2,:),200,'b.'); hold on;
% Add textCell
for ii = 1:size(OSI,2)
    text(OSI(1,ii)+.02, OSI(2,ii)+.02,num2str(ii),'FontSize',10);
end
figuremax = 1;
xlim([0 figuremax])
ylim([0 figuremax])
plot([0 figuremax],[0 figuremax],'k','LineWidth',2)
title('Orientation tuning','Fontsize',20)
ylabel('OSI Condition 1','Fontsize',20);
xlabel('OSI Condition 0','Fontsize',20);

%% Make the HWHM figure
if params.FitTuningCurve
    figure; set(gcf,'units','normalized','Position',[0.2 0.1 0.5 0.7],'color','w')
    scatter([Fit(1,:).HWHM],[Fit(2,:).HWHM],200,'r.'); hold on;
    figuremax = 100;
    xlim([0 figuremax])
    ylim([0 figuremax])
    plot([0 figuremax],[0 figuremax],'k','LineWidth',2)
    title('Orientation tuning - Detected vs Non Detected','Fontsize',20)
    ylabel('Half Width at Half Maximum (Condition 1)','Fontsize',20);
    xlabel('Half Width at Half Maximum (Condition 0)','Fontsize',20);
end

end
