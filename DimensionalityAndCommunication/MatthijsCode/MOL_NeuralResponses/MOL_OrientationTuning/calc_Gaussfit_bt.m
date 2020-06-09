function [MU, SD, R, BASELINE, OSI, DSI, EDGES, FIT,succesflag] = calc_Gaussfit_bt(ori,resp,showFig)
% function [MU, SD, R, BASELINE, OSI, DSI, EDGES, FIT,succesflag] = calc_Gaussfit_bt(ori,resp,showFig)

%% Define parameters:
params.nBootstrap       = 20;   %n iterations of bootstrapping
params.fBtTrials        = 0.90; %Fraction of trials taken randomly each bootstrap iteration
params.degTolerance     = 22.5; %Degrees from original fit that bootstrap fits are regarded as the same
params.fracFitTolerance = 0.90; %fraction of bootstrap iterations that needs to fit same preferred orientation for success

%% Get fit of mean of ALL trials:
all_oris                = unique(ori); %Get orientations
resp_perOri_all         = NaN(size(all_oris)); %init temp vector
std_perOri_all          = NaN(size(all_oris)); %init temp vector

for iOri = 1:length(all_oris) %Get mean and std of the response to each orientation:
    resp_perOri_all(iOri)   = mean(resp(ori==all_oris(iOri)));
    std_perOri_all(iOri)    = std(resp(ori==all_oris(iOri)));
end
resp_baseline           = min(resp_perOri_all);
resp_perOri_all         = resp_perOri_all - resp_baseline;

% [Gaussfit_ALL]           = calc_Gaussfit_mean(all_oris,resp_perOri_all,std_perOri_all,1);

%% Bootstrapping: Get fit n times for subselection of trials:
btMu            = NaN(params.nBootstrap,2); %init output vars
btSD            = NaN(params.nBootstrap,2);
btR             = NaN(params.nBootstrap,2);
btBASELINE      = NaN(params.nBootstrap,1);

resp_perOri             = NaN(size(all_oris)); %init temp vector

for iBt = 1:params.nBootstrap %loop over bootstrap iterations:
    randselec       = rand(length(resp),1)<params.fBtTrials;    %Get random fraction of trials
    for iOri = 1:length(all_oris)                               %Get mean response for each orientation
        resp_perOri(iOri) = mean(resp(ori==all_oris(iOri) & randselec));
    end
    resp_perOri         = resp_perOri - min(resp_perOri_all); %subtract baseline response)

    [Gaussfit]          = calc_Gaussfit_mean(all_oris,resp_perOri,zeros(size(all_oris)),0);

    btMu(iBt,:)         = [Gaussfit.OriPref Gaussfit.OriNull]; %store output vars
    btSD(iBt,:)         = [Gaussfit.TuningWidth Gaussfit.TuningWidth];
    btR(iBt,:)          = [Gaussfit.Apref Gaussfit.Anull];
    btBASELINE(iBt,:)   = Gaussfit.Baseline;
end

%% Check consistency of the fit:

fracFitAccepted     = sum(btMu(:,1) > mean(btMu(:,1))-params.degTolerance & btMu(:,1) <  mean(btMu(:,1))+params.degTolerance) / params.nBootstrap;
if fracFitAccepted>params.fracFitTolerance
    succesflag = 1;
else 
    succesflag = 0;
end

%% Final fit output parameters:
MU                  = mean(btMu);
SD                  = mean(btSD);
R                   = mean(btR);
BASELINE            = mean(btBASELINE);

%% create line for the plot:
nBinsFit = 100;
EDGES = linspace(0,360,nBinsFit);
FIT = Gaussfit.Baseline+Gaussfit.Apref*exp(-((EDGES-Gaussfit.OriPref)/Gaussfit.TuningWidth).^2)...
    +Gaussfit.Anull*exp(-((EDGES-Gaussfit.OriPref-180)/Gaussfit.TuningWidth).^2);

%Compute OSI and DSI from fit:
[OSI,DSI,~]             = calc_OSIDSI(EDGES,FIT,0);

% Plot figure if requested
if showFig
    %% Make figure:
    figure; 
    set(gcf,'color','w'); 
    %Plot original responses:
    errorbar(all_oris,resp_perOri_all,std_perOri_all,'LineWidth',2); hold on;
    %Plot fit with increased x-axis resolution:
    plot(EDGES,FIT,'LineWidth',3)
    %Figure make-up:
    title('Orientation tuning fit')
    set(gca,'XTick',all_oris,'XTickLabels',all_oris,'FontSize', 15)
    xlim([min(all_oris)-10 max(all_oris)+10])
    ylim([0 max(resp_perOri_all+std_perOri_all)*1.1+0.01]);
    ylabel('Response (Hz)','Fontsize',20);
    xlabel('DEGREE OF ORIENTATION'); 
end

%OUtput is [MU, SD, R, BASELINE, OSI, DSI, EDGES, FIT]
end