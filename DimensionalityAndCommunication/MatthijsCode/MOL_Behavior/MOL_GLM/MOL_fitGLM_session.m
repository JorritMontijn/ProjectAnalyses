function [B,CV_Perf] = MOL_fitGLM_session(allX,Y)


params.crossval     = 4;

K                   = numel(unique(Y));             %Number of outcome levels (categorical response options)
[N,P]               = size(allX);                   %N = num observations, P = num predictors
trialchunks         = round(linspace(1,N,params.crossval+1)); %Divide trials into chunks for training and testing
B                   = NaN(P+1,K-1,params.crossval); %Init coeff matrix
CV_Perf             = NaN(1,params.crossval);       %Init crossvalidated performance output

for iCV = 1:params.crossval %loop over test/training set alternations
    idx_test            = trialchunks(iCV):trialchunks(iCV+1);
    idx_train           = ~ismember(1:N,idx_test);
    
    %Fit model:
    [B(:,:,iCV),dev,stats]       = mnrfit(allX(idx_train,:),Y(idx_train),'model','nominal','Interactions','on');
    %Evaluate model for test predictors:
    [Y_prob_est]        = mnrval(B(:,:,iCV),allX(idx_test,:),stats,'model','nominal');
    %Compute percentage correct predictions of held out test set responses:
    [~,Y_est]           = max(Y_prob_est,[],2);
    CV_Perf(iCV)        = sum(Y_est == Y(idx_test))/numel(idx_test);
end
%Average over crossvalidated chunks:
B           = mean(B,3);
CV_Perf     = mean(CV_Perf);
    
% [B,dev,stats] = mnrfit(X,Y,'model','nominal','Interactions','on');
% 
% LL = stats.beta - 1.96.*stats.se;
% UL = stats.beta + 1.96.*stats.se;
% 
% figure;
% bar(B','grouped');
% 
% %% FIGURE:
% 
% figure; set(gcf,'color','w','units','normalized','Position', [0.1 0.1 .8 .8]); hold all;
% visprobepos = 0.3;
% auprobepos = 0.5;
% 
% subplot(1,2,1); hold all;
% % Labelset.trialType      = {'Y'         'X'         'P'          'C'         };
% Colorset.trialType      = {[1 0.2 0.2] [0.2 0.2 1] [0.5 0.5 0.5] [0.9 0 0.9]};
% xvec                    = [exp(linspace(log(visprobepos),log(visprobepos),100))' exp(linspace(log(auprobepos),log(4000),100))'];
% [pp_au,pp_ci_au,~]      = mnrval(B,xvec,stats,'model','nominal');
% 
% for iResp=1:3
%     CurrLineHandle = shadedErrorBar(xvec(:,2),pp_au(:,iResp),pp_ci_au(:,iResp),{'-','LineWidth',1,'Color',Colorset.trialType{iResp}},0.3);
% end
% 
% %Figure Make up
% auprobepos = 0.5;
% auticks = [10 50 100 1000 4000];
% xlim([auprobepos 4000*1.02])
% ylim([0 1])
% xlabel('Change in frequency (Hz)','FontSize', 20)
% ylabel('% Response Auditory','FontSize', 20)
% set(gca,'FontSize',15,'YAxisLocation','left','linewidth',2)
% set(gca,'XScale','log');
% XTickLabels = ['Probe' num2cell(auticks)];
% set(gca,'Xdir','reverse','XTick',[auprobepos auticks],'XTickLabels',XTickLabels);
% box on
% 
% %% 
% subplot(1,2,2); hold all;
% 
% xvec                    = [exp(linspace(log(visprobepos),log(90),100))' exp(linspace(log(auprobepos),log(auprobepos),100))'];
% [pp_vis,pp_ci_vis,~]    = mnrval(B,xvec,stats,'model','nominal');
% 
% 
% for iResp=1:3
%     CurrLineHandle = shadedErrorBar(xvec(:,1),pp_vis(:,iResp),pp_ci_vis(:,iResp),{'-','LineWidth',1,'Color',Colorset.trialType{iResp}},0.3);
% end
% 
% % Make up
% visprobepos = 0.3;
% visticks = [2 5 15 30 90];
% xlim([visprobepos 90*1.02])
% ylim([0 1])
% xlabel('Change in orientation (degrees)','FontSize', 20)
% ylabel('% Response Visual','FontSize', 20)
% set(gca,'FontSize',15,'YAxisLocation','right')
% set(gca,'linewidth',2)
% set(gca,'XScale','log');
% XTickLabels = ['Probe' num2cell(visticks)];
% set(gca,'Xdir','normal','XTick',[visprobepos visticks],'XTickLabels',XTickLabels);
% box on
% 
% %%
% xvec = exp(linspace(log(auprobepos),log(4000),100))';
% 
% ln(audio/probe) = B(1,1) + B(1,2)*xvec;
% h = B(1,1) + B(2,1)*xvec;
% exp(B(1,1) + B(2,1)*xvec);

end



%% Separate Logistic Regression
% % Construct Predictor matrix:
% idx = Y_Response~=1;
% X = X_visualOriChange(idx);
% Y = Y_Response(idx);
% Y = Y==2;
% [logitCoef_vis,dev,STATS] = glmfit(X,Y,'binomial','link','logit');
% 
% % Construct Predictor matrix:
% idx = Y_Response~=2;
% % X = X_visualOriChange(idx);
% X = X_audioFreqChange(idx);
% 
% Y = Y_Response(idx);
% % Y = Y==2;
% [logitCoef_au,dev,STATS] = glmfit(X,Y,'binomial','link',linkfunction);
% % [logitCoef,dev] = glmfit(x,[y' n'],'binomial','link',linkfunction);
% 
% %% Figure:
% close all; figure;
% newx = linspace(min(X),max(X),100);
% logitFit = glmval(logitCoef,newx,linkfunction);
% plot(X,Y,'bs',newx,logitFit,'r-');
% xlabel('Visual Intensity'); ylabel('Percentage left choice');
% ylim([0 1]);