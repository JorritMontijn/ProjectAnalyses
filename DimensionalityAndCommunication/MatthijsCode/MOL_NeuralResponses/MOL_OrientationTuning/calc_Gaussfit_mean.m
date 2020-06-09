function [Gaussfit] = calc_Gaussfit_mean(all_oris,ori_resp,ori_var,showFig)
params.shiftoris = 1;
params.targetidx = 2;

%% Make everything column vectors:
all_oris    = all_oris(:); ori_resp    = ori_resp(:); ori_var     = ori_var(:);

%% Get preferred and opposite direction:
nOrientations   = numel(all_oris);
ori_pref        = find(ori_resp==max(ori_resp)); ori_pref = ori_pref(1); %control if exactly similar response in two conditions
ori_oppo        = mod(ori_pref+nOrientations/2,nOrientations); if ori_oppo==0;    ori_oppo = nOrientations;   end

%% Shift
if params.shiftoris %center ori_pref at 2:
    shift = params.targetidx-ori_pref;
    shift_idx = mod([1:nOrientations]+shift-1,nOrientations)+1;
    % all_oris(shift_idx)     = all_oris;
    ori_resp(shift_idx)     = ori_resp;
    ori_var(shift_idx)      = ori_var;
    ori_pref = params.targetidx;
    ori_oppo = params.targetidx+4;
end

gaussfitmethod = 4;
Gaussfit = struct();

%% Fit double Gaussian to orientation tuning curve.
switch gaussfitmethod
    case 1        % Basic fit double gaussian, no contraints
        [CurvObj,GOF] = fit(all_oris,ori_resp,'gauss2'); % Y = a1*exp(-((x-b1)/c1)^2)+a2*... exp(-((x-b2)/c2)^2)
        
        fitparams                        = coeffvalues(CurvObj);
        Gaussfit.Apref           = fitparams(1);
        Gaussfit.Anull           = fitparams(4);
        Gaussfit.TuningWidth                 = mean([fitparams(3) fitparams(6)]);
        Gaussfit.Baseline          = [];
        Gaussfit.GOF                     = GOF;
        
    case 2        % With setting lower and upper bound constraints, i.e. fix pref and opposite directions:

        %       [a1,    b1,                 c1,     a2,     b2,                 c2];
        Lower = [-inf,  all_oris(ori_pref), -inf,   -inf,   all_oris(ori_oppo), -inf];
        Upper = [inf,   all_oris(ori_pref), inf,    inf,    all_oris(ori_oppo), inf];
        
        [CurvObj,GOF] = fit(all_oris,ori_resp,'gauss2','Lower',Lower,'Upper',Upper); % Y = a1*exp(-((x-b1)/c1)^2)+a2*... exp(-((x-b2)/c2)^2)
        Gaussfit.GOF                     = GOF;

    case 3        % Gaussian fitting including a baseline response and equal variance of both gaussians:
        
        ft              = fittype('r+a1*exp(-((x-b1)/c)^2)+a2*exp(-((x-b2)/c)^2)');
        % ft              = fittype('r+a1*exp(-((x-b1)^2/(2*c)^2))+a2*exp(-((x-b2)^2/(2*c)^2))');
        
        % where x is the stim direction, b1 is preferred direction, b2 is opposite direction,
        % r is the baseline response, c is the SD of the gaussians, and a1 and a2 are the peak amplitude
        
        %Parameters     = ['a1'     'a2'    'b1'                'b2'                'c'     'r'];
        StartPoints     = [0.5,     0.5,    all_oris(ori_pref), all_oris(ori_oppo), 20,     0];
        Lower           = [0,       0,      all_oris(ori_pref), all_oris(ori_oppo), 0,      0];
        Upper           = [inf,     inf,    all_oris(ori_pref), all_oris(ori_oppo), inf,    inf];
        
        [CurvObj,GOF] = fit(all_oris,ori_resp,ft,'StartPoint',StartPoints,'Lower',Lower,'Upper',Upper);
        
        fitparams                        = coeffvalues(CurvObj);
        Gaussfit.Apref           = fitparams(1);
        Gaussfit.Anull           = fitparams(2);
        Gaussfit.TuningWidth                 = fitparams(5);
        Gaussfit.Baseline          = fitparams(6);
        Gaussfit.GOF                     = GOF;
        
    case 4        % Gaussian fitting including a baseline response and equal variance of both gaussians:
        
        ft              = fittype('r+a1*exp(-((x-b1)/c)^2)+a2*exp(-((x-b1-180)/c)^2)');
        % ft              = fittype('r+a1*exp(-((x-b1)^2/(2*c)^2))+a2*exp(-((x-b2)^2/(2*c)^2))');
        
        % where x is the stim direction, b1 is preferred direction, b2 is opposite direction,
        % r is the baseline response, c is the SD of the gaussians, and a1 and a2 are the peak amplitude
        rmax = max(ori_resp)-min(ori_resp);
        
        %Parameters     = ['a1'         'a2'    'b1'                        'c'     'r'];
        StartPoints     = [0.5,         0.5,    all_oris(ori_pref),         20,     3];
        Lower           = [0,           0,      all_oris(ori_pref)-45,      5,      0];
        Upper           = [rmax*1.5,    rmax*1.5, all_oris(ori_pref)+45,    180,    80];
        
        %Parameters     = ['a1'         'a2'    'b1'                        'c'     'r'];
        StartPoints     = [0.5,         0.5,    all_oris(ori_pref),         20,     min(ori_resp)];
        Lower           = [0,           0,      all_oris(ori_pref)-60,      5,      min(ori_resp)];
        Upper           = [rmax*1.5,    rmax*1.5, all_oris(ori_pref)+60,    180,    min(ori_resp)];
        
        [CurvObj,GOF] = fit(all_oris,ori_resp,ft,'StartPoint',StartPoints,'Lower',Lower,'Upper',Upper, 'Robust', 'LAR');
        fitparams                       = coeffvalues(CurvObj);
        if params.shiftoris
            Gaussfit.OriPref = mod(fitparams(3)-shift*45,360);
        else                Gaussfit.OriPref                = fitparams(3);
        end
        Gaussfit.OriNull                = mod(Gaussfit.OriPref+180,360);
        Gaussfit.Apref                  = fitparams(1);
        Gaussfit.Anull                  = fitparams(2);
        Gaussfit.TuningWidth            = fitparams(4);
        Gaussfit.Baseline               = fitparams(5);
        Gaussfit.GOF                    = GOF;
        Gaussfit.HWHM                   = Gaussfit.TuningWidth * sqrt(2*log(2));         %% Calculate tuning width, Half Width at half maximum
end

if params.shiftoris %shift responses back:
    shift_idx = mod([1:nOrientations]-shift-1,nOrientations)+1;
    % all_oris(shift_idx)     = all_oris;
    ori_resp(shift_idx)     = ori_resp;
    ori_var(shift_idx)      = ori_var;   
end

%% Plot figure if requested
if showFig
    %% create line for the plot:
    fit_xspacing = 100;
    fit_oris = linspace(0,360,fit_xspacing);
    fit_resp = Gaussfit.Baseline+Gaussfit.Apref*exp(-((fit_oris-Gaussfit.OriPref)/Gaussfit.TuningWidth).^2)...
        +Gaussfit.Anull*exp(-((fit_oris-Gaussfit.OriPref-180)/Gaussfit.TuningWidth).^2);
    %% Make figure:
    figure; 
    set(gcf,'color','w'); 
    %Plot original responses:
    errorbar(all_oris,ori_resp,ori_var,'LineWidth',2); hold on;
    %Plot fit with increased x-axis resolution:
    plot(fit_oris,fit_resp,'LineWidth',3)
    %Figure make-up:
    title('Orientation tuning fit')
    set(gca,'XTick',all_oris,'XTickLabels',all_oris,'FontSize', 15)
    xlim([min(all_oris)-10 max(all_oris)+10])
    ylim([0 max(ori_resp+ori_var)*1.1+0.01]);
    ylabel('Response (Hz)','Fontsize',20);
    xlabel('DEGREE OF ORIENTATION'); 
end


end




% From Atallah et al 2012:

% Orientation selectivity index (OSI) was calculated as 1 - circular variance 
% (Ringach et al., 1997). Responses to the 12 grating directions were fit with 
% orientation tuning curves i.e., a sum-of-Gaussians (Figure 1, Figure 3 and 
% Figure 4). The Gaussians are forced to peak 180 degrees apart, and to have 
% the same tuning sharpness (?) but can have unequal height (Apref and Anull, 
% to account for direction selectivity), and a constant baseline B. The tuning 
% sharpness was measured as ? (2 ln(2))1/2, i.e., the half-width at half height 
% (HWHH). Direction selectivity index (DSI) was calculated as (Rpref – Rnull) 
% / (Rpref + Rnull), where Rpref is the response at the preferred direction 
% and Rnull is the response 180 degrees away from the preferred direction. 
% Contrast-response curves were fit with the hyperbolic ratio equation (Albrecht 
% and Hamilton, 1982): R(C) = Rmax cn / (C50n + cn) + Roffset, where c is contrast, 
% C50 is the semisaturation contrast, and n is a fitting exponent that describes 
% the shape of the curve, Rmax determines the gain, and Roffset is the baseline 
% response.

