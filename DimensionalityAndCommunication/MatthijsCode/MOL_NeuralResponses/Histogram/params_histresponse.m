function params = params_histresponse() 
% Parameter settings calculating responses and showing histogram
% All time is in microseconds
params                      = struct();

%% parameters for window of interest:
params.t_pre                = -1e6;
params.t_post               = 2.5e6;

%% parameters for binning:
%Method can be either large squared bins  or smoothing with 1 ms bins with window ('smoothing') 

params.histmethod           = 'individual'; %Whether to bin (and smooth) individual trials 'individual' or or all trials together ('total')
params.binsize              = 1e3;          %Size of the bins
params.smoothing            = 1;    
params.conv_win             = 'gaussian';   %Type of window used for smoothing {flat, gaussian)
params.conv_twin            = 0.3e6;        %Window size for smoothing
params.conv_sigma           = 0.025e6;        %sd of gaussian window for smoothing

%% Z-score?
params.zscore               = 1;            %Whether to convert psth matrix to zscore (x-meanxbaseline/sdbaseline)

%% parameters for calculating response
params.respcalcmethod       = 'mean';        %Which method to calculate response {'max','mean','div','AUC'}
params.trialrespmethod      = 'individual';
params.twin_baseline_start  = -1e6;         %Start window for calculating baseline
params.twin_baseline_stop   = -0.2e6;       %End window baseline
params.twin_resp_start      = 0e6;          %Start window for calculating response
params.twin_resp_stop       = 1e6;          %End window response %note that if stimulus is of variable length give input as 'twin_resp_stop_trial'
params.subtr_baseline       = 1;            %Subtract baseline or not

%% parameters for PSTH plot
params.plottype             = 'errorline';  %Either plot 'errorline' with mean and SEM, or plot 'bins' to plot classical bins 
% params.plottype             = 'bins';  %Either plot 'errorline' with mean and SEM, or plot 'bins' to plot classical bins 

end