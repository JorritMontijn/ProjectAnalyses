function MOL_Psy_2Sided_SplitHistory(varargin)
%% Get input arguments:
sessionData     = varargin{1};
trialData       = varargin{2};

%% General settings:
% showIndFig          = 0;
% showResFig          = 1;
% sessioncounter      = 0;

historyshift            = 1;

%% Loop over mice and then sessions:
mouseids = unique(sessionData.mousename)';
for mou = 1:length(mouseids)
    
    %Give title based on sessiondata information:
    mouseid         = mouseids{mou};
    sesselec        = unique(sessionData.session_ID(strcmp(sessionData.mousename,mouseid)))';
    
    [tempsessionData,temptrialData] = MOL_getTempPerSes(sesselec,sessionData,trialData);%Get the sessionData for each mouse individually:
    
    if tempsessionData.VisualLeftCorrectSide
        X_visualCorrect     = temptrialData.leftCorrect & strcmp(temptrialData.responseSide,'L');
        X_visualMiss        = strcmp(temptrialData.trialType,'X') & ~strcmp(temptrialData.responseSide,'L');
        X_audioCorrect      = temptrialData.rightCorrect & strcmp(temptrialData.responseSide,'R');
        X_audioMiss         = strcmp(temptrialData.trialType,'Y') & ~strcmp(temptrialData.responseSide,'R');
    else
        X_visualCorrect     = temptrialData.leftCorrect & strcmp(temptrialData.responseSide,'R');
        X_visualMiss        = strcmp(temptrialData.trialType,'X') & ~strcmp(temptrialData.responseSide,'R');
        X_audioCorrect      = temptrialData.rightCorrect & strcmp(temptrialData.responseSide,'L');
        X_audioMiss         = strcmp(temptrialData.trialType,'Y') & ~strcmp(temptrialData.responseSide,'L');
    end
    
    idx             = [false(historyshift,1); X_visualCorrect(1:end-historyshift)];
    datafields      = fieldnames(temptrialData);
    for field = 1:length(datafields)
        selectrialData.(datafields{field}) = temptrialData.(datafields{field})(idx);
    end
    [x_vis_viscorr(mou,:),y_vis_viscorr(mou,:),vispar_viscorr(mou),x_au_viscorr(mou,:),y_au_viscorr(mou,:),aupar_viscorr(mou)] = MOL_Psy_Trials(tempsessionData,selectrialData);
    
    idx             = [false(historyshift,1); X_visualMiss(1:end-historyshift)];
    for field = 1:length(datafields)
        selectrialData.(datafields{field}) = temptrialData.(datafields{field})(idx);
    end
    [x_vis_vismiss(mou,:),y_vis_vismiss(mou,:),vispar_vismiss(mou),x_au_vismiss(mou,:),y_au_vismiss(mou,:),aupar_vismiss(mou)] = MOL_Psy_Trials(tempsessionData,selectrialData); %#ok<*AGROW,ASGLU>
    
    idx             = [false(historyshift,1); X_audioCorrect(1:end-historyshift)];
    for field = 1:length(datafields)
        selectrialData.(datafields{field}) = temptrialData.(datafields{field})(idx);
    end
    [x_vis_aucorr(mou,:),y_vis_aucorr(mou,:),vispar_aucorr(mou),x_au_aucorr(mou,:),y_au_aucorr(mou,:),aupar_aucorr(mou)] = MOL_Psy_Trials(tempsessionData,selectrialData);
    
    idx             = [false(historyshift,1); X_audioMiss(1:end-historyshift)];
    for field = 1:length(datafields)
        selectrialData.(datafields{field}) = temptrialData.(datafields{field})(idx);
    end
    [x_vis_aumiss(mou,:),y_vis_aumiss(mou,:),vispar_aumiss(mou),x_au_aumiss(mou,:),y_au_aumiss(mou,:),aupar_aumiss(mou)] = MOL_Psy_Trials(tempsessionData,selectrialData);
    
end

% for iM = 1:size(y_vis,1)
%     plot(x_vis(iM,:),y_vis(iM,:),'Color',[0.75,0.6,0.6],'LineWidth',3); hold on;
% end
% plot(mean(x_vis,1),mean(y_vis,1),'r','LineWidth',3); hold on;

%% Figure (Previous trial VISUAL correct or miss)
figure; set(gcf,'color','w','units','normalized','Position', [0.1 0.1 .8 .8]);
h = suptitle('Psychometric curve depending on previous VISUAL trial Correct (G) or Miss (R)');
subplot(1,2,1); hold all;

for iM = 1:size(x_au_viscorr,1)
    plot(aupar_viscorr(iM).curve(:,1),aupar_viscorr(iM).curve(:,2),'Color',[0.6,0.75,0.6],'LineWidth',3);
end
curveconc = [aupar_viscorr(:).curve];
plot(aupar_viscorr(1).curve(:,1),mean(curveconc(:,2:2:end),2),'g','LineWidth',3);

for iM = 1:size(x_au_vismiss,1)
    plot(aupar_vismiss(iM).curve(:,1),aupar_vismiss(iM).curve(:,2),'Color',[0.75,0.6,0.6],'LineWidth',3); 
end
curveconc = [aupar_vismiss(:).curve];
plot(aupar_vismiss(1).curve(:,1),mean(curveconc(:,2:2:end),2),'r','LineWidth',3);

%Figure Make up
auprobepos = 0.5;
auticks = [10 50 100 1000 4000];
xlim([auprobepos 4000*1.02])
ylim([0 1])
xlabel('Change in frequency (Hz)','FontSize', 20)
ylabel('% Response Auditory','FontSize', 20)
set(gca,'FontSize',15,'YAxisLocation','left','linewidth',2)
set(gca,'XScale','log');
XTickLabels = ['Probe' num2cell(auticks)];
set(gca,'Xdir','reverse','XTick',[auprobepos auticks],'XTickLabels',XTickLabels);
box on

%%
subplot(1,2,2); hold all;
for iM = 1:size(x_vis_viscorr,1)
    plot(vispar_viscorr(iM).curve(:,1),vispar_viscorr(iM).curve(:,2),'Color',[0.6,0.75,0.6],'LineWidth',3);
end
curveconc = [vispar_viscorr(:).curve];
plot(vispar_viscorr(1).curve(:,1),mean(curveconc(:,2:2:end),2),'g','LineWidth',3);

for iM = 1:size(x_au_vismiss,1)
    plot(vispar_vismiss(iM).curve(:,1),vispar_vismiss(iM).curve(:,2),'Color',[0.75,0.6,0.6],'LineWidth',3); 
end
curveconc = [vispar_vismiss(:).curve];
plot(vispar_vismiss(1).curve(:,1),mean(curveconc(:,2:2:end),2),'r','LineWidth',3);

% Make up
visprobepos = 0.3;
visticks = [2 5 15 30 90];
xlim([visprobepos 90*1.02])
ylim([0 1])
xlabel('Change in orientation (degrees)','FontSize', 20)
ylabel('% Response Visual','FontSize', 20)
set(gca,'FontSize',15,'YAxisLocation','right')
set(gca,'linewidth',2)
set(gca,'XScale','log');
XTickLabels = ['Probe' num2cell(visticks)];
set(gca,'Xdir','normal','XTick',[visprobepos visticks],'XTickLabels',XTickLabels);
box on

%% Figure (Previous trial AUDITORY correct or miss)
figure; set(gcf,'color','w','units','normalized','Position', [0.1 0.1 .8 .8]);
h = suptitle('Psychometric curve depending on previous AUDITORY trial Correct (G) or Miss (R)');
subplot(1,2,1); hold all;

for iM = 1:size(x_au_viscorr,1)
    plot(aupar_aucorr(iM).curve(:,1),aupar_aucorr(iM).curve(:,2),'Color',[0.6,0.75,0.6],'LineWidth',3);
end
curveconc = [aupar_aucorr(:).curve];
plot(aupar_aucorr(1).curve(:,1),mean(curveconc(:,2:2:end),2),'g','LineWidth',3);

for iM = 1:size(x_au_vismiss,1)
    plot(aupar_aumiss(iM).curve(:,1),aupar_aumiss(iM).curve(:,2),'Color',[0.75,0.6,0.6],'LineWidth',3); 
end
curveconc = [aupar_aumiss(:).curve];
plot(aupar_aumiss(1).curve(:,1),mean(curveconc(:,2:2:end),2),'r','LineWidth',3);

%Figure Make up
auprobepos = 0.5;
auticks = [10 50 100 1000 4000];
xlim([auprobepos 4000*1.02])
ylim([0 1])
xlabel('Change in frequency (Hz)','FontSize', 20)
ylabel('% Response Auditory','FontSize', 20)
set(gca,'FontSize',15,'YAxisLocation','left','linewidth',2)
set(gca,'XScale','log');
XTickLabels = ['Probe' num2cell(auticks)];
set(gca,'Xdir','reverse','XTick',[auprobepos auticks],'XTickLabels',XTickLabels);
box on

%%
subplot(1,2,2); hold all;
for iM = 1:size(x_vis_viscorr,1)
    plot(vispar_aucorr(iM).curve(:,1),vispar_aucorr(iM).curve(:,2),'Color',[0.6,0.75,0.6],'LineWidth',3);
end
curveconc = [vispar_aucorr(:).curve];
plot(vispar_aucorr(1).curve(:,1),mean(curveconc(:,2:2:end),2),'g','LineWidth',3);

for iM = 1:size(x_au_vismiss,1)
    plot(vispar_aumiss(iM).curve(:,1),vispar_aumiss(iM).curve(:,2),'Color',[0.75,0.6,0.6],'LineWidth',3); 
end
curveconc = [vispar_aumiss(:).curve];
plot(vispar_aumiss(1).curve(:,1),mean(curveconc(:,2:2:end),2),'r','LineWidth',3);

% Make up
visprobepos = 0.3;
visticks = [2 5 15 30 90];
xlim([visprobepos 90*1.02])
ylim([0 1])
xlabel('Change in orientation (degrees)','FontSize', 20)
ylabel('% Response Visual','FontSize', 20)
set(gca,'FontSize',15,'YAxisLocation','right')
set(gca,'linewidth',2)
set(gca,'XScale','log');
XTickLabels = ['Probe' num2cell(visticks)];
set(gca,'Xdir','normal','XTick',[visprobepos visticks],'XTickLabels',XTickLabels);
box on


end