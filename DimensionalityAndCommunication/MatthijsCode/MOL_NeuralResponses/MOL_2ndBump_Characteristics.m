[Data] = MOL_GetData('E:','CHDET',{'ChangeDetectionConflict'},{'2003' '2009' '2010' '2011' '2012' '2013'},[],{'sessionData' 'trialData' 'spikeData'});
sessionData     = Data.sessionData;
trialData       = Data.trialData;
spikeData       = Data.spikeData;

%% Remove last 20 trials:
trialData = MOL_RemoveLastnTrials(trialData,20);

%% Filter good neurons:
spikeData = MOL_filterNeurons(sessionData,trialData,spikeData);

%% Filter only V1 neurons:
spikeFields = fieldnames(spikeData);
for fld = 1:length(spikeFields)
    spikeData.(spikeFields{fld}) = spikeData.(spikeFields{fld})(strcmp(spikeData.area,'V1'));
end

%% Get significant first and second bump responses:
[sessionData,trialData,spikeData] = MOL_calc_signresp(sessionData,trialData,spikeData);

%% Venn diagram of population overlap:
nTotalNeurons = length(spikeData.session_ID);
Bump_0 = sum(~spikeData.sign_firstbump & ~spikeData.sign_secondbump) / nTotalNeurons;
Bump_1 = sum(spikeData.sign_firstbump) / nTotalNeurons;
Bump_2 = sum(spikeData.sign_secondbump) / nTotalNeurons;
Bump_12 = sum(spikeData.sign_firstbump & spikeData.sign_secondbump) / nTotalNeurons;

%For a 3-circle venn, Z is a 7 element vector [z1 z2 z3 z12 z13 z23 z123]
Z = [Bump_0 Bump_1 Bump_2 0 0 Bump_12 0];
figure;
[H] = venn(Z,'FaceAlpha', 1);
fprintf('No significant response: %3.0f neurons, %2.1f %%\n',Bump_0*nTotalNeurons,Bump_0*100)
fprintf('1st Bump: %3.0f neurons, %2.1f %%\n',Bump_1*nTotalNeurons,Bump_1*100)
fprintf('2nd Bump: %3.0f neurons, %2.1f %%\n',Bump_2*nTotalNeurons,Bump_2*100)
fprintf('Overlap: %3.0f neurons, %2.1f %%\n',Bump_12*nTotalNeurons,Bump_12*100)

%% Make figure of laminar distribution of first and second bump:

%Histogram with binning on depth:
edges           = 0:100:1200;
histY_all       = histc(spikeData.ChannelY,edges);
histY_bump_1    = histc(spikeData.ChannelY(spikeData.sign_firstbump),edges);
histY_bump_2    = histc(spikeData.ChannelY(spikeData.sign_secondbump),edges);

laminardepthfig = figure; set(gcf,'units','normalized','Position',[0.2 0.2 0.4 0.6],'color','w');
set(laminardepthfig, 'DefaultLineLineWidth', 2);

hold all;
plot(histY_all,edges,'-k')
plot(histY_bump_1,edges,'-b')
plot(histY_bump_2,edges,'-r')

ylabel('Depth from dura (um)','FontSize', 15)
xlabel('Number of sign units','FontSize',15)
set(gca,'YDir','reverse','linewidth',2)
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'fontsize',15,'FontWeight','bold')

legend({'Overall' '1st' '2nd'})
legend boxoff

% normalize to the amount of units at that depth:
histY_bump_1_norm    = histY_bump_1 ./ histY_all;
histY_bump_2_norm    = histY_bump_2 ./ histY_all;

laminardepthfignorm = figure; set(gcf,'units','normalized','Position',[0.2 0.2 0.4 0.6],'color','w');
set(laminardepthfignorm, 'DefaultLineLineWidth', 2);

hold all;
plot(histY_bump_1_norm,edges,'-b')
plot(histY_bump_2_norm,edges,'-r')

ylabel('Depth from dura (um)','FontSize', 15)
xlabel('Percentage of sign units','FontSize',15)
set(gca,'YDir','reverse','linewidth',2)
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'fontsize',15,'FontWeight','bold')

legend({'1st' '2nd'})
legend boxoff

%% AUC and Mutual Information discrimination score for 1st vs 2nd bump neurons:

[edges,AUC,MI] = MOL_FeatureDiscr_timebins(sessionData,trialData,spikeData);
  
%% 
AUC = AUC_visual;
MI = MI_visual;
%%
AUC = AUC_auditory;
MI = MI_auditory;

%% 
Idx_0           = ~spikeData.sign_firstbump & ~spikeData.sign_secondbump;
Idx_1_only      = spikeData.sign_firstbump & ~spikeData.sign_secondbump;
Idx_2_only      = ~spikeData.sign_firstbump & spikeData.sign_secondbump;
Idx_1           = spikeData.sign_firstbump;
Idx_2           = spikeData.sign_secondbump;
Idx_12          = spikeData.sign_firstbump & spikeData.sign_secondbump;

figure; hold all; set(gcf,'units','normalized','Position',[0.1 0.1 0.5 0.7],'color','w'); LineHandles = [];
meanforplot = nanmean(AUC(Idx_1_only,:));
semforplot = nanstd(AUC(Idx_1_only,:)) / sqrt(sum(Idx_1_only));
CurrLineHandle = shadedErrorBar(edges,meanforplot,semforplot,{'-','LineWidth',1.5,'Color',[0 0 1]},1);
LineHandles(1) = CurrLineHandle.mainLine;
meanforplot = nanmean(AUC(Idx_2_only,:));
semforplot = nanstd(AUC(Idx_2_only,:)) / sqrt(sum(Idx_2_only));
CurrLineHandle = shadedErrorBar(edges,meanforplot,semforplot,{'-','LineWidth',1.5,'Color',[1 0 0]},1);
LineHandles(2) = CurrLineHandle.mainLine;
meanforplot = nanmean(AUC(Idx_12,:));
semforplot = nanstd(AUC(Idx_12,:)) / sqrt(sum(Idx_12));
CurrLineHandle = shadedErrorBar(edges,meanforplot,semforplot,{'-','LineWidth',1.5,'Color',[0 1 0]},1);
LineHandles(3) = CurrLineHandle.mainLine;
meanforplot = nanmean(AUC(Idx_0,:));
semforplot = nanstd(AUC(Idx_0,:)) / sqrt(sum(Idx_0));
CurrLineHandle = shadedErrorBar(edges,meanforplot,semforplot,{'-','LineWidth',1.5,'Color',[0.5 0.5 0.5]},1);
LineHandles(4) = CurrLineHandle.mainLine;

set(gca, 'XTick', -1e6:0.2e6:2e6, 'XTickLabels', [-1e6:0.2e6:2e6]/1e6,'FontSize', 20)
xlabel('Time (s)','FontSize', 20)
ylabel('AUC (improvement over shuffle)','FontSize', 18)
legend(LineHandles,{'1st bump excl' '2nd bump excl' 'Both' 'None'},'FontSize',15);
legend boxoff;

%% 
figure; hold all; set(gcf,'units','normalized','Position',[0.1 0.1 0.5 0.7],'color','w'); LineHandles = [];
meanforplot = nanmean(MI(spikeData.sign_firstbump & ~spikeData.sign_secondbump,:));
semforplot = nanstd(MI(spikeData.sign_firstbump & ~spikeData.sign_secondbump,:)) / sqrt(sum(spikeData.sign_firstbump & ~spikeData.sign_secondbump));
CurrLineHandle = shadedErrorBar(edges,meanforplot,semforplot,{'-','LineWidth',1.5,'Color',[0 0 1]},1);
LineHandles(1) = CurrLineHandle.mainLine;
meanforplot = nanmean(MI(~spikeData.sign_firstbump & spikeData.sign_secondbump,:));
semforplot = nanstd(MI(~spikeData.sign_firstbump & spikeData.sign_secondbump,:)) / sqrt(sum(~spikeData.sign_firstbump & spikeData.sign_secondbump));
CurrLineHandle = shadedErrorBar(edges,meanforplot,semforplot,{'-','LineWidth',1.5,'Color',[1 0 0]},1);
LineHandles(2) = CurrLineHandle.mainLine;

set(gca, 'XTick', -1e6:0.2e6:2e6, 'XTickLabels', [-1e6:0.2e6:2e6]/1e6,'FontSize', 20)
xlabel('Time (s)','FontSize', 20)
ylabel('Mutual information (bits) (Rate x Feature)','FontSize', 20)
legend(LineHandles,{'1st bump neurons' '2nd bump neurons'},'FontSize',15);
legend boxoff;

%% 
figure; hold all; set(gcf,'units','normalized','Position',[0.1 0.1 0.5 0.7],'color','w'); LineHandles = [];
meanforplot = nanmean(AUC(:,:));
semforplot = nanstd(AUC(:,:)) / sqrt(sum(~isnan(AUC(:,1))));
CurrLineHandle = shadedErrorBar(edges,meanforplot,semforplot,{'-','LineWidth',1.5,'Color',[0 0 1]},1);
LineHandles(1) = CurrLineHandle.mainLine;

set(gca, 'XTick', -1e6:0.2e6:2e6, 'XTickLabels', [-1e6:0.2e6:2e6]/1e6,'FontSize', 20)
xlabel('Time (s)','FontSize', 20)
ylabel('AUC (improvement over shuffle)','FontSize', 15)
legend(LineHandles,{'All neurons'},'FontSize',15);
legend boxoff;


