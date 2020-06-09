function MOL_OptoMeanResponse(varargin)
%% Get input arguments:
sessionData     = varargin{1};
trialData       = varargin{2};
spikeData       = varargin{3};

%% Analysis of the effects of optical stimulation only:
% ToDo:
% -

%%
% spikeData.cellType  = calc_neurontype(spikeData); %Get Interneurons and Pyramidal Neurons:

%% Parameter settings for PSTH
params = params_histresponse_pmal(); 

% params.histmethod           = 'individual'; %Whether to bin (and smooth) individual trials 'individual' or or all trials together ('total')
params.eventofinterest      = 'stimStart';

meanbaseline            = [];
meanopto                = [];

for sesid = unique(sessionData.session_ID)'
    %% Get the relevant data for each session individually:
    [temptrialData,tempspikeData] = MOL_getTempPerSes(sesid,trialData,spikeData);
    if ~isempty(tempspikeData.ts)
        if strcmp(sessionData.Experiment(sessionData.session_ID==sesid),'Bars')
            params.twin_resp_stop_trial = temptrialData.stimEnd - temptrialData.stimStart;  %different response window per trial (variable stimulus length)
        end
        
        hist_mat = [];
        resp = [];
        for iNeuron = 1:length(tempspikeData.ts)
            events_ts           = temptrialData.(params.eventofinterest);          %Get events
            spikes_ts           = tempspikeData.ts{iNeuron};                    %Get spikes for this neuron
            
            %Get histogram per neuron
            [edges,hist_mat(:,:,iNeuron)]       = calc_psth(events_ts,spikes_ts,params);
            %Get responses per neuron
            [resp(:,iNeuron),~]                 = calc_resp_from_psth(edges,hist_mat(:,:,iNeuron),params);
        end
    end
    meanbaseline            = [meanbaseline nanmean(resp(temptrialData.hasphotostim==0,:))];
    meanopto                = [meanopto nanmean(resp(temptrialData.hasphotostim==1,:))];
end

%% Cap values at 0.01 Hz:
meanopto(meanopto<0.01) = 0.012;
meanbaseline(meanbaseline<0.01) = 0.012;

%% Convert to log scale:
meanbaseline            = log10(meanbaseline);
meanopto                = log10(meanopto);

%% Make scatter figure
figure; set(gcf,'units','normalized','Position',[0.2 0.1 0.5 0.7],'color','w'); hold on;

scatter(meanbaseline,meanopto,40,'b','filled'); %Plot dots
% errorbar(meanbaseline, meanopto, semopto, 'k.');
% herrorbar(meanbaseline, meanopto, sembaseline,'k.');
plot([-5 5],[-5 5],'k','LineWidth',3); %Plot equity line

%Make up of figure:
grid on
ticks = [0.01 0.03 0.1 0.3 1 3 10 30 100];
set(gca,'XTick',log10(ticks),'XTickLabels',ticks,'FontSize', 20)
set(gca,'YTick',log10(ticks),'YTickLabels',ticks,'FontSize', 20)
set(gca,'linewidth',3)
xlim([min(log10(ticks)) max(log10(ticks))*1.1])
ylim([min(log10(ticks)) max(log10(ticks))*1.1])
xlabel('Response OPTO OFF(Hz)','FontSize', 20)
ylabel('Response OPTO ON (Hz)','FontSize', 20)

fprintf('n = %d neurons\n',numel(meanbaseline));

end

