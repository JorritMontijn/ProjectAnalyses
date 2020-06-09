function [VonMisesfit] = calc_VonMisesfit_mean_method2(all_oris,ori_resp,ori_var,showFig)
%% Make everything column vectors:
all_oris    = all_oris(:); ori_resp    = ori_resp(:); ori_var     = ori_var(:);
all_ori_rad = deg2rad(all_oris);

%% Get preferred and opposite direction:
nOrientations   = numel(all_oris);
ori_pref        = find(ori_resp==max(ori_resp)); ori_pref = ori_pref(1); %control if exactly similar response in two conditions
ori_oppo        = mod(ori_pref+nOrientations/2,nOrientations); if ori_oppo==0;    ori_oppo = nOrientations;   end

%% Parameters
params.nComponents = 2;

%% Fit

% ft              = fittype('r+a1*exp(-((x-b1)/c)^2)+a2*exp(-((x-b1-180)/c)^2)');
% ft              = fittype('r+a1*exp(-((x-b1)^2/(2*c)^2))+a2*exp(-((x-b2)^2/(2*c)^2))');
% ft              = fittype('r+a1.*exp(kappa*cos((x-b1)*pi/180)-1) + a2.*exp(kappa*cos(((x-b1)+pi)*pi/180)-1)');
% ft              = fittype('r+a1.*exp(kappa*cos((x-b1)*pi/180)-1) + a2.*exp(kappa*cos((x-b1+180)*pi/180)-1)');
ft              = fittype('r+a1.*exp(kappa*cos((x-b1)*pi/180)) + a2.*exp(kappa*cos((x-b1+180)*pi/180))');

% ft              = fittype('r+a1.*exp(kappa*cos((x-b1))-1) + a2.*exp(kappa*cos((x-b1+pi))-1)');

% where x is the stim direction, b1 is preferred direction, b2 is opposite direction,
% r is the baseline response, kappa is the dispersion parameter and a1 and a2 are the peak amplitude

%Parameters     = ['a1'         'a2'    'b1'                        'kappa'     'r'];
StartPoints     = [0.5,         0.5,    all_oris(ori_pref),         2,     min(ori_resp)];
% Lower           = [0,           0,      all_oris(ori_pref)-45,      0.1,   0];
% Upper           = [rmax*1.5,    rmax*1.5, all_oris(ori_pref)+45,    8,     80];
Lower           = [0,           0,      0,                          0,      0];
Upper           = [20,          20,     360,                        8,     80];

% %Parameters     = ['a1'         'a2'    'b1'                        'c'     'r'];
% StartPoints     = [0.5,         0.5,    all_oris(ori_pref),         20,     min(ori_resp)];
% Lower           = [0,           0,      all_oris(ori_pref)-60,      5,      min(ori_resp)];
% Upper           = [rmax*1.5,    rmax*1.5, all_oris(ori_pref)+60,    180,    min(ori_resp)];

[CurvObj,GOF]                   = fit(all_oris,ori_resp,ft,'StartPoint',StartPoints,'Lower',Lower,'Upper',Upper, 'Robust', 'LAR','Display','off');
fitparams                       = coeffvalues(CurvObj);

% VonMisesfit.OriPref                 = rad2deg(fitparams(3));
% VonMisesfit.OriNull                = mod(VonMisesfit.OriPref+180,360);
% VonMisesfit.Apref                  = fitparams(1);
% VonMisesfit.Anull                  = fitparams(2);
% VonMisesfit.TuningWidth            = fitparams(4);
% VonMisesfit.Baseline               = fitparams(5);
% VonMisesfit.GOF                    = GOF;
% VonMisesfit.HWHM                   = VonMisesfit.TuningWidth * sqrt(2*log(2));         %% Calculate tuning width, Half Width at half maximum

VonMisesfit.mu                     = [fitparams(3) mod(fitparams(3)+180,360)];
VonMisesfit.kappa                  = [fitparams(4) fitparams(4)];
VonMisesfit.componentProportion    = [fitparams(1) fitparams(2)];
VonMisesfit.Baseline               = fitparams(5);
VonMisesfit.GOF                    = GOF;

% VonMisesfit.Apref                  = fitparams(1);
% VonMisesfit.Anull                  = fitparams(2);
% VonMisesfit.TuningWidth            = fitparams(4);
% VonMisesfit.Baseline               = fitparams(5);
% VonMisesfit.GOF                    = GOF;
% VonMisesfit.HWHM                   = VonMisesfit.TuningWidth * sqrt(2*log(2));         %% Calculate tuning width, Half Width at half maximum


%% Make the figure:
if showFig
    figure;
    set(gcf,'color','w');
    all_oris    = [all_oris; 360];
    ori_resp    = [ori_resp; ori_resp(1)];
    ori_var     = [ori_var; ori_var(1)];
    
    %Plot original responses:
    errorbar(all_oris,ori_resp,ori_var,'LineWidth',2); hold on;
    
    %% Plot von mises fit:
    
    h = plot(CurvObj);
    set(h,'LineWidth',2)

%     edges               = linspace(-pi,pi,1000)'; % The pdf() function expects a column-vector as input
%     likelihoods         = VonMisesfit.pdf(edges);
%     fit_resp            = likelihoods/max(likelihoods) * max(ori_resp); % Go back to firing rates
%     edges               = edges + pi;
%     edges               = edges*180/pi; % Back to angles from radians
%     plot(edges,fit_resp,'LineWidth',2);
    
    %Figure makeup
    title('Orientation tuning fit')
    set(gca,'XTick',all_oris,'XTickLabels',all_oris,'FontSize', 15)
    xlim([min(all_oris)-10 max(all_oris)+10])
    ylim([0 max(ori_resp+ori_var)*1.1+0.01]);
    ylabel('Response (Hz)','Fontsize',20);
    xlabel('DEGREE OF ORIENTATION');
    
    legend(['nComponents ' num2str(params.nComponents)])
end


end