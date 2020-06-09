function [MU, KAPPA, R, BASELINE, OSI, DSI, EDGES, FIT, GOF, succesflag] = calc_VonMisesfit_bt(ori,resp,showFig)

%% Define parameters:
params.nBootstrap       = 50;   %n iterations of bootstrapping
params.fBtTrials        = 0.8;  %Fraction of trials taken randomly each bootstrap iteration
params.degTolerance     = 22.5;   %Degrees from original fit that bootstrap fits are regarded as the same
params.fracFitTolerance = 0.90; %fraction of bootstrap iterations that needs to fit same preferred orientation for success

params.VMfitMethod      = 2;

%% Get fit of mean of ALL trials:
all_oris                = unique(ori); %Get orientations
resp_perOri_all         = NaN(size(all_oris)); %init temp vector
std_perOri_all          = NaN(size(all_oris)); %init temp vector

for iOri = 1:length(all_oris) %Get mean and std of the response to each orientation:
    resp_perOri_all(iOri)   = mean(resp(ori==all_oris(iOri)));
    std_perOri_all(iOri)    = std(resp(ori==all_oris(iOri)));
end
resp_baseline               = min(resp_perOri_all);

switch params.VMfitMethod
    case 1
        [VMfit_ALL]           = calc_VonMisesfit_mean(all_oris,resp_perOri_all-resp_baseline,std_perOri_all,0);
    case 2
        [VMfit_ALL]           = calc_VonMisesfit_mean_method2(all_oris,resp_perOri_all,std_perOri_all,0);
end

% [VMfit_ALL]           = calc_VonMisesfit_mean(all_oris,resp_perOri_all,std_perOri_all,1);

%% Bootstrapping: Get fit n times for subselection of trials:
btMu            = NaN(params.nBootstrap,2); %init output vars
btKappa         = NaN(params.nBootstrap,2);
btR             = NaN(params.nBootstrap,2);
btBaseline      = NaN(params.nBootstrap,1);


for iBt = 1:params.nBootstrap %loop over bootstrap iterations:
    resp_perOri             = NaN(size(all_oris)); %init temp vector
    
    while any(isnan(resp_perOri))
        randselec               = rand(length(resp),1)<params.fBtTrials;    %Get random fraction of trials
        for iOri = 1:length(all_oris)                               %Get mean response for each orientation
            resp_perOri(iOri) = mean(resp(ori==all_oris(iOri) & randselec));
        end
    end

    switch params.VMfitMethod
        case 1
            [VMfit]           = calc_VonMisesfit_mean(all_oris,resp_perOri,zeros(size(all_oris)),0);
        case 2
            [VMfit]           = calc_VonMisesfit_mean_method2(all_oris,resp_perOri,zeros(size(all_oris)),0);
    end
    
    btMu(iBt,:)       = VMfit.mu; %store output vars
    btKappa(iBt,:)    = VMfit.kappa;
    btR(iBt,:)        = VMfit.componentProportion;
    btBaseline(iBt,:) = VMfit.Baseline;
end

%% Sort the components by size:
[~,maxidx_bt]       = max(btR,[],2);                        %Get the largest of the two components
idx                 = bsxfun(@eq, cumsum(ones(size(btR)),2), maxidx_bt); %convert to logical index

btR                 = [sum(btR.*idx, 2) sum(btR.*~idx, 2)];
btKappa             = [sum(btKappa.*idx, 2) sum(btKappa.*~idx, 2)];
btMu                = [sum(btMu.*idx, 2) sum(btMu.*~idx, 2)];

%% Sort the overall fit on the basis of component size:
[~,maxidx_all]      = max(VMfit_ALL.componentProportion); %Get the largest of the two components
MeanFitDegs         = VMfit_ALL.mu(maxidx_all);            %Get the preferred orientation (mu of largest component in radians)

%% Check consistency of the fit:

fracFitAccepted     = sum(mod(btMu(:,1),180) > mod(MeanFitDegs,180)-params.degTolerance &...
                          mod(btMu(:,1),180) < mod(MeanFitDegs,180)+params.degTolerance) / params.nBootstrap;

if fracFitAccepted>params.fracFitTolerance
    succesflag = 1;
else 
    succesflag = 0;
end

%% Compute output fit:
nBinsFit            = 1000;
switch params.VMfitMethod
    case 1
        edges               = linspace(-pi,pi,nBinsFit)'; % The pdf() function expects a column-vector as input
        likelihoods         = VMfit_ALL.pdf(edges);
        fit_resp            = likelihoods/max(likelihoods) * max(resp_perOri_all); % Go back to firing rates
        FIT                 = fit_resp + resp_baseline; %Add baseline response again
        EDGES               = rad2deg(edges+pi);
    case 2
        edges               = linspace(0,360,nBinsFit)'; % The pdf() function expects a column-vector as input
        ft                  = fittype('r+a1.*exp(kappa*cos((x-b1)*pi/180)) + a2.*exp(kappa*cos((x-b1+pi)*pi/180))');
        FIT                 = feval(ft,VMfit.componentProportion(1),VMfit.componentProportion(2),VMfit.mu(1),VMfit.kappa(1),VMfit.Baseline,edges);
        
%         fit_resp            = likelihoods/max(likelihoods) * max(resp_perOri_all); % Go back to firing rates
%         FIT                 = fit_resp + resp_baseline; %Add baseline response again
        EDGES               = edges;

end
    

%% Final fit output parameters:
MU                      = median(btMu);
KAPPA                   = median(btKappa);
BASELINE                = median(btBaseline);
R                       = median(btR);

ft                      = fittype('r+a1.*exp(kappa*cos((x-b1)*pi/180)) + a2.*exp(kappa*cos((x-b1+180)*pi/180))');
FIT                     = feval(ft,R(1),R(2),MU(1),KAPPA(1),BASELINE,edges);
EDGES                   = edges;

%Compute OSI and DSI from fit:
[OSI,DSI,~]             = calc_OSIDSI(EDGES,FIT,0);

if showFig
    figure; set(gcf,'color','w');
    hold all;

    %Plot original responses:
    errorbar([all_oris; 360],[resp_perOri_all; resp_perOri_all(1)],[std_perOri_all; std_perOri_all(1)],'Color',[0 0 0.4],'LineWidth',2); hold on;
    
    for iBt = 1:params.nBootstrap %loop over bootstrap iterations:
        btfit = feval(ft,btR(iBt,1),btR(iBt,2),btMu(iBt,1),btKappa(iBt,1),btBaseline(iBt,1),EDGES);
        plot(EDGES,btfit,'Color',[0.7 0.7 0.7],'LineWidth',0.5)
    end
    
    %Plot mean fit:
    plot(EDGES,FIT,'Color',[0 0 0],'LineWidth',2)

    %Figure makeup
    title('Orientation tuning fit')
    set(gca,'XTick',all_oris,'XTickLabels',all_oris,'FontSize', 15)
    xlim([min(all_oris)-10 max(all_oris)+10])
    ylim([0 max(resp_perOri_all+std_perOri_all)*1.1+0.01]);
    ylabel('Response (Hz)','Fontsize',20);
    xlabel('DEGREE OF ORIENTATION');
    
end

%% Calculate Goodness of fit:
resp_fit                        = feval(ft,R(1),R(2),MU(1),KAPPA(1),BASELINE,all_oris);
GOF                             = sum((resp_perOri_all - resp_fit).^2); %Take goodness of fit as Sum of Squared Errors

%% OUT
%Output is: [MU, KAPPA, R, BASELINE, OSI, DSI, EDGES, GOF, FIT]

end