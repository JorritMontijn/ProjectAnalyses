function [VonMisesfit] = calc_VonMisesfit_mean(all_oris,ori_resp,ori_var,showFig)
%% Make everything column vectors:
all_oris    = all_oris(:); ori_resp    = ori_resp(:); ori_var     = ori_var(:);
all_ori_rad = deg2rad(all_oris);

%% Parameters
params.nComponents = 2;

%% von Mises Fit:
norm_test               = ori_resp/sum(ori_resp); % Normalize to unity
norm_round_test         = round(norm_test*1e4); % Get integers values so ones() function
% can accept it as an argument

samples                 = []; %Init var
for i=1:length(all_ori_rad)
    samples=[samples,ones(1,norm_round_test(i))*all_ori_rad(i)]; % In every iteration
    %     creates a row vector of length = norm_round_test(i) filled with
    %     all_ori_rad(i), and finally add it to the previous row vector
end

samples             = (samples - pi)'; % fitmvmdist expects a column vector between -pi and pi
% figure; plot(samples);
VonMisesfit         = fitmvmdist(samples,params.nComponents); % Construct von Mises mixture model with 2 components

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
    edges               = linspace(-pi,pi,1000)'; % The pdf() function expects a column-vector as input
    likelihoods         = VonMisesfit.pdf(edges);
    fit_resp            = likelihoods/max(likelihoods) * max(ori_resp); % Go back to firing rates
    edges               = edges + pi;
    edges               = edges*180/pi; % Back to angles from radians
    plot(edges,fit_resp,'LineWidth',2);
    
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