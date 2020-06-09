function params = params_timefreq() 
% Parameter settings for time-frequency analyses:
% All time is in microseconds
params                      = struct();

%% parameters for window of interest:
params.t_pre                = -1e6;
params.t_post               = 1.5e6;

%% parameters 

% params.blabla

%% Z-score?
params.zscore               = 0;            %Whether to convert to z-scored power

%% parameters for calculating response
params.respcalcmethod       = 'mean';        %Which method to calculate response {'max','mean','div','AUC'}
params.twin_baseline_start  = -1e6;         %Start window for calculating baseline
params.twin_baseline_stop   = -0.2e6;       %End window baseline
params.twin_resp_start      = 0e6;          %Start window for calculating response
params.twin_resp_stop       = 0.3e6;        %End window response
params.subtr_baseline       = 1;            %Subtract baseline or not

%% parameters for plot
params.plottype             = 'errorline';  %Either plot 'errorline' with mean and SEM, or plot 'bins' to plot classical bins 
% params.plottype             = 'bins';  %Either plot 'errorline' with mean and SEM, or plot 'bins' to plot classical bins 

end