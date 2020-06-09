function MOL_PopulationRecField(varargin)
%% Get input arguments:
sessionData     = varargin{1};
trialData       = varargin{2};
spikeData       = varargin{3};

%% Get parameters for hist and responses:
params          = params_histresponse_pmal(); % All time is in microseconds

%% Extra parameters:
params.eventofinterest      = 'stimStart';
params.tuningofinterest     = 'visualOri';
params.allowsplitdata       = 1;
params.ShowFig              = 0;
params.latencyCorrection    = 150e3;

%%

for sesid = unique(sessionData.session_ID)'
    %% Get the relevant data for each session individually:
    [temptrialData,tempspikeData] = MOL_getTempPerSes(sesid,trialData,spikeData);
    
    %% Loop over neurons
    for iNeuron = 1:length(tempspikeData.ts)
        
        %Get histogram per neuron:
        events_ts = temptrialData.(params.eventofinterest)(temptrialData.session_ID == tempspikeData.session_ID(iNeuron));
        spikes_ts = tempspikeData.ts{iNeuron};
        
        if params.allowsplitdata
            [splitidxs,colors]     = MOL_GUI_GetSplits(temptrialData);                 %Get splits on the basis of the GUI
        else  colors{1} = 'k';
        end
        
        [edges,hist_mat]        = calc_psth(events_ts,spikes_ts,params);
        %Apply latency correction:
        edges                   = edges-params.latencyCorrection;
        
        allTuningConds          = unique(temptrialData.(params.tuningofinterest)(~isnan(temptrialData.(params.tuningofinterest))));
        
        %Get tuning per neuron:
        for sp = 1:size(splitidxs,2)
            
            %align all psths:
            maxlen                      = ceil(max(temptrialData.stimEnd-temptrialData.stimStart)/1000)*1000;
            hist_array                  = NaN(size(hist_mat,1),maxlen/1000);
            
            for iTuningCond = 1:length(allTuningConds)
                %Get selection of trials for this condition:
                idx_trials = temptrialData.(params.tuningofinterest)==allTuningConds(iTuningCond) & splitidxs{sp};
                params.twin_resp_stop_trial = temptrialData.stimEnd(idx_trials) - temptrialData.stimStart(idx_trials);  %different response window per trial (variable stimulus length)
                
                hist_cells{iTuningCond} = mean(hist_mat(idx_trials,edges>0 & edges<mean(params.twin_resp_stop_trial)));
            end
            
            %%
            maxlen                      = max(cellfun(@length,hist_cells));
            hist_array                  = NaN(length(allTuningConds),maxlen);
            
            for iTuningCond = 1:length(allTuningConds)
                padding = (maxlen-length(hist_cells{iTuningCond}))/2;
                hist_array(iTuningCond,:) = [zeros(1,floor(padding)) hist_cells{iTuningCond} zeros(1,ceil(padding))];
            end
            
            map = back_project(hist_array, mod(allTuningConds+180,360));
            figure;
            imagesc(map)
            
            Fit2dGaussian(X, 
            A0 = [1,0,50,0,50,0];   % Inital (guess) parameters

            
        end
    end
    
end

end