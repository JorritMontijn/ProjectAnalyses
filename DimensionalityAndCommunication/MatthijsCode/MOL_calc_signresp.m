function [sessionData,trialData,spikeData] = MOL_calc_signresp(varargin)

%% e.g. 
% [sessionData,trialData,spikeData] = MOL_calc_signresp(sessionData,trialData,spikeData);

%% Get input arguments:
if nargin==3
    sessionData     = varargin{1};
    trialData       = varargin{2};
    spikeData       = varargin{3};
else
    [Data] = MOL_GetData('E:','CHDET',{'ChangeDetectionConflict'},{'2003' '2009' '2010' '2011' '2012' '2013'},[],{'sessionData' 'trialData' 'spikeData'});
    sessionData     = Data.sessionData;
    trialData       = Data.trialData;
    spikeData       = Data.spikeData;
    %% Remove last 20 trials:
    trialData = MOL_RemoveLastnTrials(trialData,20);
end

%% General parameters:
params                      = params_histresponse(); % All time is in microseconds
params.histmethod           = 'individual'; %Whether to bin (and smooth) individual trials 'individual' or or all trials together ('total')
params.zscore               = 0;
params.smoothing            = 0;

params.AlignOn                = 'stimChange';      %On which timestamp to align as t=0

% parameters for window of interest:
params.t_start_0              = -0.4e6;
params.t_stop_0               = -0.2e6;
params.t_start_1              = 0e6;
params.t_stop_1               = 0.2e6;
params.t_start_2              = 0.2e6;
params.t_stop_2               = 0.4e6;
% params.t_window_dur         = params.t_post-params.t_pre;
params.ttestalpha             = 0.05;

%% Get index of trialtypes to compute significance of:
idx     = true(size(trialData.trialType));
idx     = idx       & strcmp(trialData.trialType,'X');
idx     = idx       & trialData.correctResponse==1;
% idx     = idx       & trialData.visualOriChangeNorm==2;

%% Init output fields:
spikeData.sign_firstbump    = NaN(size(spikeData.session_ID));
spikeData.sign_secondbump   = NaN(size(spikeData.session_ID));

%% For each session for each neuron compute sign resp to idx:
for sesid = unique(sessionData.session_ID)'
    %% Get the relevant data for each session individually:
    [temptrialData,tempspikeData] = MOL_getTempPerSes(sesid,trialData,spikeData);
    
    %% Construct tensor:
    events_ts       = temptrialData.(params.AlignOn);
    nNeurons        = length(tempspikeData.ts);
    nTimebins       = length(params.t_pre:params.binsize:params.t_post-params.binsize);
    nTrials         = length(events_ts);
    tensor          = NaN(nNeurons,nTimebins,nTrials);
    for iNeuron = 1:nNeurons
        %Get histogram per neuron
        [edges,hist_mat]        = calc_psth(events_ts,tempspikeData.ts{iNeuron},params);
        if params.smoothing 
            tensor(iNeuron,:,:)     = hist_mat';
        else
            tensor(iNeuron,:,:)     = hist_mat' / (1e6/params.binsize); %Keep actual spike counts
        end
    end
    
    %% Compute response and baseline during specified timewindow:
    baseline    = NaN(nNeurons,nTrials);
    resp_1      = NaN(nNeurons,nTrials);
    resp_2      = NaN(nNeurons,nTrials);
    
    for iNeuron = 1:nNeurons
        baseline(iNeuron,:)       = squeeze(sum(tensor(iNeuron,edges>params.t_start_0 & edges<params.t_stop_0,:)));
        resp_1(iNeuron,:)         = squeeze(sum(tensor(iNeuron,edges>params.t_start_1 & edges<params.t_stop_1,:)));
        resp_2(iNeuron,:)         = squeeze(sum(tensor(iNeuron,edges>params.t_start_2 & edges<params.t_stop_2,:)));
    end
    
    %% Compute significant response:
    sestrialidx                 = idx(strcmp(trialData.session_ID,sesid));
    sign_firstbump              = NaN(nNeurons,1);
    sign_secondbump              = NaN(nNeurons,1);
    
    for iNeuron = 1:nNeurons
        sign_firstbump(iNeuron) = ttest(baseline(iNeuron,sestrialidx),resp_1(iNeuron,sestrialidx),'alpha',params.ttestalpha);
        sign_secondbump(iNeuron) = ttest(baseline(iNeuron,sestrialidx),resp_2(iNeuron,sestrialidx),'alpha',params.ttestalpha);
    end
    
    spikeData.sign_firstbump(strcmp(spikeData.session_ID,sesid)) = sign_firstbump;
    spikeData.sign_secondbump(strcmp(spikeData.session_ID,sesid)) = sign_secondbump;
    
end

%Convert to logical:
spikeData.sign_firstbump = spikeData.sign_firstbump==1;
spikeData.sign_secondbump = spikeData.sign_secondbump==1;

fprintf('Calculated significant 1st or 2nd bump response for %d neurons\n\n', length(spikeData.sign_firstbump))

end