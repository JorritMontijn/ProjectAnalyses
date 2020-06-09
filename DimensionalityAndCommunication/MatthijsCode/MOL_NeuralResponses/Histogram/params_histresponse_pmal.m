function params = params_histresponse_pmal() 
% Parameter settings calculating responses and showing histogram
% All time is in microseconds
params                      = struct();

%% parameters for window of interest:
params.t_pre                = -2e6;
params.t_post               = 7e6;

%% parameters for binning:
%Method can be either large squared bins  or smoothing with 1 ms bins with window ('smoothing') 

params.binsize              = 1e3;          %Size of the bins
params.smoothing            = 1;    
params.conv_win             = 'gaussian';   %Type of window used for smoothing {flat, gaussian)
params.conv_twin            = 0.3e6;        %Window size for smoothing
params.conv_sigma           = 0.05e6;       %sd of gaussian window for smoothing

%% Z-score?
params.zscore               = 0;            %Whether to convert psth matrix to zscore (x-meanxbaseline/sdbaseline)

%% parameters for calculating response
% params.trialrespmethod      = 'total';    %Whether to bin (and smooth) individual trials 'individual' or or all trials together ('total')
% params.respcalcmethod       = 'max';      %Which method to calculate response {'max','mean','div','AUC'}
params.trialrespmethod      = 'individual'; %Whether to bin (and smooth) individual trials 'individual' or or all trials together ('total')
params.respcalcmethod       = 'meanpeak';   %Which method to calculate response {'max','mean','div','AUC'}
params.twin_baseline_start  = -2e6;         %Start window for calculating baseline
params.twin_baseline_stop   = -0.2e6;       %End window baseline
params.twin_resp_start      = 0.1e6;        %Start window for calculating response
params.twin_resp_stop       = 2e6;          %End window response %note that if stimulus is of variable length give input as 'twin_resp_stop_trial'
params.subtr_baseline       = 1;            %Subtract baseline or not
params.peakwin              = 0.3e6;        %time window around peak that firing rate is averaged

%% parameters for PSTH plot
params.plottype             = 'errorline';  %Either plot 'errorline' with mean and SEM, or plot 'bins' to plot classical bins 
% params.plottype             = 'bins';  %Either plot 'errorline' with mean and SEM, or plot 'bins' to plot classical bins 

end