function [edges,AUC,MI] = MOL_FeatureDiscr_timebins(varargin)

%% e.g.
% [outputvar] = MOL_AUC_orientation(sessionData,trialData,spikeData);

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
params.binsize              = 0.1e6;          %Size of the bins
params.t_pre                = -1e6;
params.t_post               = 2e6;
params.smoothing            = 0;
params.zscore               = 0;            %Whether to convert psth matrix to zscore (x-meanxbaseline/sdbaseline)

params.AlignOn              = 'stimChange';      %On which timestamp to align as t=0

params.nshuffle             = 100;

params.min_ntrial = 5;

%% Get index of trialtypes to compute significance of:
% params.feature              = 'visualOriPostChangeNorm';
% idx     = true(size(trialData.trialType));
% idx     = idx       & strcmp(trialData.trialType,'X');
% % idx     = idx       & trialData.correctResponse==1;
% idx     = idx       & trialData.visualOriChangeNorm==3;

%%
params.feature              = 'audioFreqPostChangeNorm';
idx     = true(size(trialData.trialType));
idx     = idx       & strcmp(trialData.trialType,'Y');
idx     = idx       & trialData.audioFreqChangeNorm==3;

%% Init output fields:
edges               = [params.t_pre:params.binsize:params.t_post]; % - params.binsize/2; %Define edges

MI                  = NaN(length(spikeData.session_ID),length(edges));
AUC                 = NaN(length(spikeData.session_ID),length(edges));

%%
neuroncounter = 1;
% AUC = NaN(nNeurons,length(edges),3);

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
    
    %% Get the right indices:
    sestrialidx                 = idx(strcmp(trialData.session_ID,sesid));
    
    %Feature to discriminate:
    feature         = temptrialData.(params.feature)(sestrialidx);
    
    if numel(unique(feature))==2
        feat_2levels    = feature;
%     elseif numel(unique(feature))==3
%         keyboard();
    elseif numel(unique(feature))==4 || numel(unique(feature))==3
        feat_2levels    = feature;
        feat_2levels(ismember(feat_2levels,[1 2])) = 1;
        feat_2levels(ismember(feat_2levels,[3 4])) = 2;
    else
        feat_2levels    = feature;
        feat_2levels = feat_2levels(ismember(feat_2levels,[1 2]));
        %         if unique(feature(sestrialidx))>4
        %         keyboard();
        %         feat_2levels    = feature;
        %         feat_2levels(ismember(feat_2levels,[3 4])) = 2;
        %         feat_idx         = ismember(feat_2levels,[1 2]) & sestrialidx;
    end
    
     %% MI and AUC over time:
     MI_shuffle = NaN(params.nshuffle,1);
     AUC_shuffle = NaN(params.nshuffle,1);
     
    if sum(sestrialidx)>params.min_ntrial && sum(feat_2levels==1)>params.min_ntrial && sum(feat_2levels==2)>params.min_ntrial
        for iNeuron = 1:nNeurons
            for iTime = 1:length(edges)
                if params.nshuffle
                    for it = 1:params.nshuffle
                        [MI_shuffle(it)]           = MutualInformation(feature(randperm(length(feature))),squeeze(tensor(iNeuron,iTime,sestrialidx)));
                        [~,~,~,AUC_shuffle(it)]    = perfcurve(feat_2levels(randperm(length(feat_2levels))),squeeze(tensor(iNeuron,iTime,sestrialidx)),1);
                    end
                    [MI(neuroncounter,iTime)]       = MutualInformation(feature,squeeze(tensor(iNeuron,iTime,sestrialidx))) - mean(MI_shuffle);
                    
                    [~,~,~,AUC_temp]                = perfcurve(feat_2levels,squeeze(tensor(iNeuron,iTime,sestrialidx)),1);
                    AUC(neuroncounter,iTime)        = abs(AUC_temp-0.5) - mean(abs(AUC_shuffle-0.5)); %absolutize AUC, values <0.5 discriminate the other way around
                    
                else
                    [MI(neuroncounter,iTime)]           = MutualInformation(feature,squeeze(tensor(iNeuron,iTime,sestrialidx)));
                    [~,~,~,AUC(neuroncounter,iTime)]    = perfcurve(feat_2levels,squeeze(tensor(iNeuron,iTime,sestrialidx)),1);
                end
            end
            neuroncounter = neuroncounter + 1;
        end
    else
        neuroncounter = neuroncounter + length(tempspikeData.session_ID);
    end
    
end

%% Time correction:
AUC = AUC(:,1:size(AUC,2)-1);
MI = MI(:,1:size(MI,2)-1);
edges = edges + params.binsize/2;

end