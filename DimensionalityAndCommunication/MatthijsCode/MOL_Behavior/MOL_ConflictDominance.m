function MOL_ConflictDominance(varargin)
%% Get input arguments:
sessionData     = varargin{1};
trialData       = varargin{2};

%% General settings:
showIndFig          = 0;
showMeanFig         = 1;

set(0,'defaultAxesFontSize',20);

sessioncounter = 0;

%% Loop over mice and then sessions:
mouseids = unique(sessionData.mousename)';
for mou = 1:length(mouseids)
    
    %Give title based on sessiondata information:
    mouseid         = mouseids{mou};
    sesselec        = unique(sessionData.session_ID(strcmp(sessionData.mousename,mouseid)))';
    
    SesDominanceMat = [];
    SesnTrialsMat   = [];
    
    [tempsessionData,temptrialData] = MOL_getTempPerSes(sesselec,sessionData,trialData);%Get the sessionData for each session individually:
    
    [x_vis,x_au,FullRespMat,FullnTrialsMat] = MOL_Psy_GetTrialMat(tempsessionData,temptrialData);
    
    DominanceIndexMat = FullRespMat(:,:,1)./(FullRespMat(:,:,1) + FullRespMat(:,:,2));
    
    SesDominanceMat(:,:,mou) = DominanceIndexMat;
    SesnTrialsMat(:,:,mou)   = FullnTrialsMat;
    
    %     for ses = 1:length(sesselec)
    %         sessioncounter  = sessioncounter+1;         %Add one to the counter
    %         sesid           = sesselec(ses);            %Get sessionID for this session
    %         [tempsessionData,temptrialData] = MOL_getTempPerSes(sesid,sessionData,trialData);%Get the sessionData for each session individually:
    %
    %         [x_vis,x_au,FullRespMat] = MOL_Psy_GetTrialMat(tempsessionData,temptrialData)
    %
    %
    %         DominanceIndexMat = FullRespMat(:,:,1)./(FullRespMat(:,:,1) + FullRespMat(:,:,2));
    %
    %         %%
    %         DominanceIndexMat   = flipud(DominanceIndexMat);
    %         FullnTrialsMat      = flipud(FullnTrialsMat);
    %
    %         SesDominanceMat(:,:,ses) = DominanceIndexMat;
    %         SesnTrialsMat(:,:,ses)   = FullnTrialsMat;
    %
    %     end
    
    %% Make the figure:
    if showMeanFig
        plotDominanceMap(SesnTrialsMat(:,:,mou),SesDominanceMat(:,:,mou));
        suptitle(mouseid)
    end
end


end

function plotDominanceMap(FullnTrialsMat,DominanceIndexMat)

DominanceIndexMat   = flipud(DominanceIndexMat);
FullnTrialsMat      = flipud(FullnTrialsMat);

nCond               = size(DominanceIndexMat,1)-0.5;
params.colormap     = 'redblue';

f = figure; set(gcf,'color','w','units','normalized','Position', [0.1 0.4 .8 .4]);

switch size(DominanceIndexMat,1)
    case 2
        axislabels = {'Probe' 'Threshold'};
    case 5
        axislabels = {'Probe' 'Sub' 'Threshold' 'Supra' 'Full'};
end

%% Show table with number of trials per condition:
t = uitable(f,'Data',FullnTrialsMat,... %Create the uitable
    'ColumnName',axislabels,...
    'RowName',fliplr(axislabels),...
    'ColumnWidth',{80}, 'ColumnWidth',{80},'FontS',25);
subplot(1,2,1) %,plot(3)
pos = get(subplot(1,2,1),'position');
title('Number of trials')
set(subplot(1,2,1),'yTick',[])
set(subplot(1,2,1),'xTick',[])
set(t,'units','normalized')
set(t,'position',pos)

%% Show Heatmap with dominance indices:
subplot(1,2,2)
imagesc(DominanceIndexMat,[0 1]); hold on;
probex      = [0.5 1.5 1.5 0.5 0.5];
probey      = [nCond nCond nCond+1 nCond+1 nCond];
plot(probex,probey,'k-', 'LineWidth', 3);
audiox      = [0.5 1.5 1.5 0.5 0.5];
audioy      = [0.5 0.5 nCond nCond 0.5];
plot(audiox,audioy,'r-', 'LineWidth', 3);
visualx     = [1.5 1.5 nCond+1 nCond+1 1.5];
visualy     = [nCond nCond+1 nCond+1 nCond nCond];
plot(visualx,visualy,'b-', 'LineWidth', 3);
confx       = [1.5 nCond+1 nCond+1 1.5 1.5];
confy      = [0.5 0.5 nCond nCond 0.5];
plot(confx,confy,'m-', 'LineWidth', 3);

switch params.colormap
    case 'redblue'
        h       = coolwarm + repmat(0.1-abs(linspace(-0.1,.1,64))',1,3);
        h(h>1)  = 1;
        colormap(h);
    case 'parula'
        colormap(parula);
end
set(gca,'XTickLabels',axislabels,'YTickLabels',fliplr(axislabels),'FontSize',15)
xlabel('Visual')
ylabel('Auditory')

c = colorbar;
c.Label.String = 'Auditory Dominance Index';
end