function [aupar,vispar] = MOL_Opto_V1PPC_Behavior(varargin)
%% Get input arguments:
if nargin==2
    sessionData     = varargin{1};
    trialData       = varargin{2};
else
    [Data] = MOL_GetData('E:','CHDET',{'ChangeDetectionConflict' 'BehaviorConflict'},{'2009' '2010' '2011' '2012' '2013' '2019' '2020' '2021' '2022' '2023'},[],{'sessionData' 'trialData'});
    sessionData     = Data.sessionData;
    trialData       = Data.trialData;
    %% Remove last 20 trials:
    trialData = MOL_RemoveLastnTrials(trialData,20);
end

%% General settings:
Par.showMeanFig         = 0;
Par.showDprimeFig       = 1;
Par.normSessions        = 0;

Par.sumTrialsCondition  = 30;
Par.minPerf             = 0.4;

set(0,'defaultAxesFontSize',20)
sessioncounter      = 0;

%% Initialize structure for saving output fit parameters:
vispar  = struct();
aupar   = struct();

% fulldata_viscorr        = NaN(10,10,10,10); %Init matrix for storing all data
% fulldata_visincorr      = NaN(10,10,10,10); %Init matrix for storing all data
% fulldata_aucorr         = NaN(10,10,10,10); %Init matrix for storing all data
% fulldata_auincorr       = NaN(10,10,10,10); %Init matrix for storing all data

dVis_V1                    = NaN(10,10,10,10); %Init matrix for storing all data
dAud_V1                    = NaN(10,10,10,10); %Init matrix for storing all data
dVis_PPC                    = NaN(10,10,10,10); %Init matrix for storing all data
dAud_PPC                    = NaN(10,10,10,10); %Init matrix for storing all data

%% Trim post change delayed inhibition:
trialfields = fieldnames(trialData);
for iF = 1:length(trialfields)
    trialData.(trialfields{iF})    = trialData.(trialfields{iF})(~(trialData.PostChangeOptoStart > 0));
end

%% Loop over mice and then sessions for V1:
mouseids = unique(sessionData.mousename)';
for mou = 1:length(mouseids)
    
    mouseid                 = mouseids{mou};
    sesidx                  = strcmp(sessionData.mousename,mouseid) & sessionData.UseOpto==1 ...
        & strcmp(sessionData.PhotostimArea,'V1')...
        & sessionData.Photostimpower > 2;
    sesselec                = unique(sessionData.session_ID(sesidx));
    
    for ses = 1:length(sesselec)
        sesid               = sesselec(ses);
        sessioncounter      = sessioncounter+1;         %Add one to the overall counter
        
        [tempsessionData,temptrialData] = MOL_getTempPerSes(sesid,sessionData,trialData);%Get the sessionData for each session individually:
        
        [visconditions,auconditions,FullRespMat,FullnTrialsMat] = MOL_Psy_GetTrialMatOpto(tempsessionData,temptrialData);
        
        if numel(visconditions) > 2
            FullRespMat = FullRespMat(:,[1 3 end],:,:); visconditions  = visconditions([1 3 end]);
            FullnTrialsMat = FullnTrialsMat(:,[1 3 end],:); 
        end
        if numel(auconditions) > 2
            FullRespMat = FullRespMat([1 3 end],:,:,:); auconditions  = auconditions([1 3 end]);
            FullnTrialsMat = FullnTrialsMat([1 3 end],:,:); auconditions  = auconditions([1 3 end]);
        end
        
        %         fulldata_viscorr(mou,ses,1:size(selectionopto,2),1:size(selectionvis,2))        = percsidecorrectvisual; %Store in matrix
        %         fulldata_visincorr(mou,ses,1:size(selectionopto,2),1:size(selectionvis,2))      = percsideincorrectvisual; %Store in matrix
        
        %if session meets criteria for performance:
        if FullRespMat(end,1,1,1) > Par.minPerf && FullRespMat(1,end,1,2) > Par.minPerf
            %% Compute d-prime for each condition:
            FullRespMat_trialn = FullRespMat .* repmat(FullnTrialsMat,1,1,1,3);
            
            for iTrial = 1:2
                for iOpto = 1:size(FullRespMat_trialn,3)
                    outcome         = NaN(3,3);
                    outcome(1,:)    = squeeze(FullRespMat_trialn(1,1+iTrial,iOpto,:)); %visual trials
                    outcome(2,:)    = squeeze(FullRespMat_trialn(1+iTrial,1,iOpto,:)); %audio trials
                    outcome(3,:)    = squeeze(FullRespMat_trialn(1,1,iOpto,:)); %probe trials
                    outcome         = outcome(:,[2 1 3]); %Swap response coding, visual first.
                    
                    if sum(outcome(:))>Par.sumTrialsCondition
                        [dVis_V1(mou,ses,iTrial,iOpto),dAud_V1(mou,ses,iTrial,iOpto),~,~] = ...
                            MOL_Fit_2ADC_Full_Session(tempsessionData,trialData,0,outcome);
                    end
                end
            end
        end
    end
end

%% Loop over mice and then sessions for PPC:
mouseids = unique(sessionData.mousename)';
for mou = 1:length(mouseids)
    
    mouseid                 = mouseids{mou};
    sesidx                  = strcmp(sessionData.mousename,mouseid) & sessionData.UseOpto==1 ...
        & strcmp(sessionData.PhotostimArea,'PPC')...
        & sessionData.Photostimpower > 2;
    sesselec                = unique(sessionData.session_ID(sesidx));
    
    for ses = 1:length(sesselec)
        sesid               = sesselec(ses);
        sessioncounter      = sessioncounter+1;         %Add one to the overall counter
        
        [tempsessionData,temptrialData] = MOL_getTempPerSes(sesid,sessionData,trialData);%Get the sessionData for each session individually:
        
        [visconditions,auconditions,FullRespMat,FullnTrialsMat] = MOL_Psy_GetTrialMatOpto(tempsessionData,temptrialData);
        
        if numel(visconditions) > 2
            FullRespMat = FullRespMat(:,[1 3 end],:,:); visconditions  = visconditions([1 3 end]);
            FullnTrialsMat = FullnTrialsMat(:,[1 3 end],:); 
        end
        if numel(auconditions) > 2
            FullRespMat = FullRespMat([1 3 end],:,:,:); auconditions  = auconditions([1 3 end]);
            FullnTrialsMat = FullnTrialsMat([1 3 end],:,:); auconditions  = auconditions([1 3 end]);
        end
        
        %         fulldata_viscorr(mou,ses,1:size(selectionopto,2),1:size(selectionvis,2))        = percsidecorrectvisual; %Store in matrix
        %         fulldata_visincorr(mou,ses,1:size(selectionopto,2),1:size(selectionvis,2))      = percsideincorrectvisual; %Store in matrix
        
        %if session meets criteria for performance:
        if FullRespMat(end,1,1,1) > Par.minPerf && FullRespMat(1,end,1,2) > Par.minPerf
            %% Compute d-prime for each condition:
            FullRespMat_trialn = FullRespMat .* repmat(FullnTrialsMat,1,1,1,3);
            
            for iTrial = 1:2
                for iOpto = 1:2
                    outcome         = NaN(3,3);
                    outcome(1,:)    = squeeze(FullRespMat_trialn(1,1+iTrial,iOpto,:)); %visual trials
                    outcome(2,:)    = squeeze(FullRespMat_trialn(1+iTrial,1,iOpto,:)); %audio trials
                    outcome(3,:)    = squeeze(FullRespMat_trialn(1,1,iOpto,:)); %probe trials
                    outcome         = outcome(:,[2 1 3]); %Swap response coding, visual first.
                    
                    if sum(outcome(:))>Par.sumTrialsCondition
                        [dVis_PPC(mou,ses,iTrial,iOpto),dAud_PPC(mou,ses,iTrial,iOpto),~,~] = ...
                            MOL_Fit_2ADC_Full_Session(tempsessionData,trialData,0,outcome);
                    end
                end
            end
        end
    end
end

%% Organize data structure:
%Put all sessions below each other and remove singleton dimension:
% if ndims(fulldata_aucorr) == 4
%     fulldata_aucorr     = reshape(fulldata_aucorr,[size(fulldata_aucorr,1)*size(fulldata_aucorr,2) size(fulldata_aucorr,3) size(fulldata_aucorr,4)]);
%     fulldata_auincorr   = reshape(fulldata_auincorr,[size(fulldata_auincorr,1)*size(fulldata_auincorr,2) size(fulldata_auincorr,3) size(fulldata_auincorr,4)]);
%
%     fulldata_viscorr     = reshape(fulldata_viscorr,[size(fulldata_viscorr,1)*size(fulldata_viscorr,2) size(fulldata_viscorr,3) size(fulldata_viscorr,4)]);
%     fulldata_visincorr   = reshape(fulldata_visincorr,[size(fulldata_visincorr,1)*size(fulldata_visincorr,2) size(fulldata_visincorr,3) size(fulldata_visincorr,4)]);
% end
%
% if Par.normSessions
%     for ises = 1:size(fulldata_aucorr,1)
%         fulldata_aucorr(ises,:,:) = fulldata_aucorr(ises,:,:) - repmat(fulldata_aucorr(ises,:,1),1,1,10);
%         fulldata_auincorr(ises,:,:) = fulldata_auincorr(ises,:,:) - repmat(fulldata_auincorr(ises,:,1),1,1,10);
%         fulldata_viscorr(ises,:,:) = fulldata_viscorr(ises,:,:) - repmat(fulldata_viscorr(ises,:,1),1,1,10);
%         fulldata_visincorr(ises,:,:) = fulldata_visincorr(ises,:,:) - repmat(fulldata_visincorr(ises,:,1),1,1,10);
%     end
% end

%% Show average response rate figure
if Par.showMeanFig
    
    figure; set(gcf,'color','w','units','normalized','Position', [0.1 0.1 .8 .8]);
    
    %for auditory amkeup:
    if aumultiply100; auconditions = auconditions*100; end
    if numel(auconditions)==1; audiolog = 0; end
    if audiolog %If logarithmic axis place probe -2 times interval on log scale
        auprobepos = exp(min(log(auconditions)) - 2*mean(diff(log(auconditions))));
    else auprobepos = 0;
    end
    
    %for visual:
    if vismultiply100; visconditions = visconditions*100; end
    if numel(visconditions)==1; visuallog = 0; end
    if visuallog %If logarithmic axis place probe -2 times interval on log scale
        visprobepos = exp(min(log(visconditions)) - 2*mean(diff(log(visconditions))));
    else visprobepos = 0;
    end
    
    % Construct x from position of the probe trial and conditions
    xdata_vis = [visprobepos visconditions];
    
    % Construct x from position of the probe trial and conditions
    xdata_au = [auprobepos auconditions];
    
    meancorrau          = squeeze(nanmean(fulldata_aucorr,1));
    meanincorrau        = squeeze(nanmean(fulldata_auincorr,1));
    stdcorrau          = squeeze(nanstd(fulldata_aucorr,1) / sqrt(sum(~isnan(fulldata_aucorr(:,1,1)))));
    stdincorrau        = squeeze(nanstd(fulldata_auincorr,1) / sqrt(sum(~isnan(fulldata_auincorr(:,1,1)))));
    
    subplot(1,2,1); hold all;
    xdata_au = [auprobepos 1/32 1/2];
    
    set(0,'Defaultlinelinewidth',5)
    optoColors_vis = {[0 0 1] [0.8 0.8 1] [0.4 0.4 1]};
    optoColors_aud = {[1 0 0] [1 0.8 0.8] [1 0.4 0.4]};
    hold all;
    
    for iOpto = 1:3
        Par.linehandles(iOpto) = errorbar(xdata_au,meancorrau(1:3,iOpto),stdcorrau(1:3,iOpto),'.-','Color',optoColors_aud{iOpto},'LineWidth',5,'MarkerSize',50);
        Par.linehandles(iOpto+3) = errorbar(xdata_au,meanincorrau(1:3,iOpto),stdincorrau(1:3,iOpto),'.:','Color',optoColors_vis{iOpto},'LineWidth',5,'MarkerSize',50);
    end
    
    %     for iOpto = 1:3
    %         errorbar(xdata_au,meancorrau(1:3,iOpto),stdcorrau(1:3,iOpto),['r' optomarkers{iOpto}],'MarkerSize',optomarkersizes(iOpto));
    %         errorbar(xdata_au,meancorrau(1:3,iOpto),stdcorrau(1:3,iOpto),'r-','MarkerSize',30);
    %         errorbar(xdata_au,meanincorrau(1:3,iOpto),stdincorrau(1:3,iOpto),['b' optomarkers{iOpto}],'MarkerSize',optomarkersizes(iOpto));
    %         errorbar(xdata_au,meanincorrau(1:3,iOpto),stdincorrau(1:3,iOpto),'b:','MarkerSize',15);
    %     end
    
    %Figure Make up
    xlim([auprobepos max(xdata_au)*1.02])
    if Par.normSessions
        ylim([-0.5 0.3])
    else
        ylim([0 1])
    end
    xlabel(auxaxislabel,'FontSize', 20)
    ylabel('Response Rate','FontSize', 20)
    set(gca,'FontSize',15,'YAxisLocation','left','linewidth',2)
    if audiolog;  set(gca,'XScale','log'); end
    if  strcmp(tempsessionData.auChangeUnit,'Oct')
        XTickLabels = ['Probe' num2cell(xdata_au(2:end))];
    else
        XTickLabels = ['Probe' num2cell(xdata_au(2:end))];
    end
    set(gca,'Xdir','reverse','XTick',xdata_au,'XTickLabels',XTickLabels);
    box on
    
    switch Par.PhotostimArea
        case 'V1'
            legend(Par.linehandles,{'Au - Control' 'Au - V1 Early' 'Au - V1 Late' 'Vis - Control' 'Vis - V1 Early' 'Vis - V1 Late'},'FontSize',15,'Location','NorthEast');
        case 'S1'
            legend(Par.linehandles([1:2 4:5]),{'Au - Control' 'Au - S1 Early' 'Vis - Control' 'Vis - S1 Early'},'FontSize',15,'Location','NorthEast');
    end
    legend(gca,'boxoff');
    
    subplot(1,2,2); hold all;
    
    meancorrvis          = squeeze(nanmean(fulldata_viscorr,1));
    meanincorrvis        = squeeze(nanmean(fulldata_visincorr,1));
    stdcorrvis          = squeeze(nanstd(fulldata_viscorr,1) / sqrt(sum(~isnan(fulldata_viscorr(:,1,1)))));
    stdincorrvis        = squeeze(nanstd(fulldata_visincorr,1) / sqrt(sum(~isnan(fulldata_visincorr(:,1,1)))));
    
    xdata_vis = [visprobepos 12 90];
    
    for iOpto = 1:3
        Par.linehandles(iOpto) = errorbar(xdata_vis,meancorrvis(1:3,iOpto),stdcorrvis(1:3,iOpto),'.-','Color',optoColors_vis{iOpto},'LineWidth',5,'MarkerSize',50);
        Par.linehandles(iOpto+3) = errorbar(xdata_vis,meanincorrvis(1:3,iOpto),stdincorrvis(1:3,iOpto),'.:','Color',optoColors_aud{iOpto},'LineWidth',5,'MarkerSize',50);
    end
    
    %     for iOpto = 1:3
    %         errorbar(xdata_vis,meancorrvis(1:3,iOpto),stdcorrvis(1:3,iOpto),['b' optomarkers{iOpto}],'MarkerSize',optomarkersizes(iOpto));
    %         errorbar(xdata_vis,meancorrvis(1:3,iOpto),stdcorrvis(1:3,iOpto),'b-','MarkerSize',30);
    %         errorbar(xdata_vis,meanincorrvis(1:3,iOpto),stdincorrvis(1:3,iOpto),['r' optomarkers{iOpto}],'MarkerSize',optomarkersizes(iOpto));
    %         errorbar(xdata_vis,meanincorrvis(1:3,iOpto),stdincorrvis(1:3,iOpto),'r:','MarkerSize',15);
    %     end
    
    % Make up
    xlim([visprobepos max(xdata_vis)*1.02])
    if Par.normSessions
        ylim([-0.5 0.3])
    else
        ylim([0 1])
    end
    xlabel(visxaxislabel,'FontSize', 20)
    ylabel('Response Rate','FontSize', 20)
    set(gca,'FontSize',15,'YAxisLocation','right')
    set(gca,'linewidth',2)
    if visuallog; set(gca,'XScale','log'); end
    XTickLabels = ['Probe' num2cell(round(xdata_vis(2:end)))];
    set(gca,'Xdir','normal','XTick',xdata_vis,'XTickLabels',XTickLabels);
    box on
    switch Par.PhotostimArea
        case 'V1'
            legend(Par.linehandles,{'Vis - Control' 'Vis - V1 Early' 'Vis - V1 Late' 'Au - Control' 'Au - V1 Early' 'Au - V1 Late'},'FontSize',15,'Location','NorthWest');
        case 'S1'
            legend(Par.linehandles([1:2 4:5]),{'Vis - Control' 'Vis - S1 Early' 'Au - Control' 'Au - S1 Early'},'FontSize',15,'Location','NorthWest');
    end
    legend(gca,'boxoff');
    
end

%% Dprime figure:
if Par.showDprimeFig
    %Put all sessions below each other and remove resulting singleton dimension:
    % V1
    dVis_V1_reshape             = squeeze(reshape(dVis_V1,size(dVis_V1,1)*size(dVis_V1,2),1,size(dVis_V1,3),size(dVis_V1,4)));
    dAud_V1_reshape             = squeeze(reshape(dAud_V1,size(dAud_V1,1)*size(dAud_V1,2),1,size(dAud_V1,3),size(dAud_V1,4)));
    
    dVis_Small_V1_mean          = squeeze(nanmean(dVis_V1_reshape(:,1,:))); dVis_Small_V1_mean = dVis_Small_V1_mean(~isnan(dVis_Small_V1_mean));
    dVis_Small_V1_sem           = squeeze(nanstd(dVis_V1_reshape(:,1,:))) / sqrt(sum(~isnan(dVis_V1_reshape(:,1,1)))); dVis_Small_V1_sem = dVis_Small_V1_sem(~isnan(dVis_Small_V1_sem));
    
    dAud_Small_V1_mean          = squeeze(nanmean(dAud_V1_reshape(:,1,:))); dAud_Small_V1_mean = dAud_Small_V1_mean(~isnan(dAud_Small_V1_mean));
    dAud_Small_V1_sem           = squeeze(nanstd(dAud_V1_reshape(:,1,:))) / sqrt(sum(~isnan(dAud_V1_reshape(:,1,1)))); dAud_Small_V1_sem = dAud_Small_V1_sem(~isnan(dAud_Small_V1_sem));
    
    dVis_Large_V1_mean          = squeeze(nanmean(dVis_V1_reshape(:,2,:))); dVis_Large_V1_mean = dVis_Large_V1_mean(~isnan(dVis_Large_V1_mean));
    dVis_Large_V1_sem           = squeeze(nanstd(dVis_V1_reshape(:,2,:))) / sqrt(sum(~isnan(dVis_V1_reshape(:,2,1)))); dVis_Large_V1_sem = dVis_Large_V1_sem(~isnan(dVis_Large_V1_sem));
    
    dAud_Large_V1_mean          = squeeze(nanmean(dAud_V1_reshape(:,2,:))); dAud_Large_V1_mean = dAud_Large_V1_mean(~isnan(dAud_Large_V1_mean));
    dAud_Large_V1_sem           = squeeze(nanstd(dAud_V1_reshape(:,2,:))) / sqrt(sum(~isnan(dAud_V1_reshape(:,2,1)))); dAud_Large_V1_sem = dAud_Large_V1_sem(~isnan(dAud_Large_V1_sem));
    % PPC
    dVis_PPC_reshape            = squeeze(reshape(dVis_PPC,size(dVis_PPC,1)*size(dVis_PPC,2),1,size(dVis_PPC,3),size(dVis_PPC,4)));
    dAud_PPC_reshape            = squeeze(reshape(dAud_PPC,size(dAud_PPC,1)*size(dAud_PPC,2),1,size(dAud_PPC,3),size(dAud_PPC,4)));
    
    dVis_Small_PPC_mean         = squeeze(nanmean(dVis_PPC_reshape(:,1,:))); dVis_Small_PPC_mean = dVis_Small_PPC_mean(~isnan(dVis_Small_PPC_mean));
    dVis_Small_PPC_sem          = squeeze(nanstd(dVis_PPC_reshape(:,1,:))) / sqrt(sum(~isnan(dVis_PPC_reshape(:,1,1)))); dVis_Small_PPC_sem = dVis_Small_PPC_sem(~isnan(dVis_Small_PPC_sem));
    
    dAud_Small_PPC_mean         = squeeze(nanmean(dAud_PPC_reshape(:,1,:))); dAud_Small_PPC_mean = dAud_Small_PPC_mean(~isnan(dAud_Small_PPC_mean));
    dAud_Small_PPC_sem          = squeeze(nanstd(dAud_PPC_reshape(:,1,:))) / sqrt(sum(~isnan(dAud_PPC_reshape(:,1,1)))); dAud_Small_PPC_sem = dAud_Small_PPC_sem(~isnan(dAud_Small_PPC_sem));
    
    dVis_Large_PPC_mean         = squeeze(nanmean(dVis_PPC_reshape(:,2,:))); dVis_Large_PPC_mean = dVis_Large_PPC_mean(~isnan(dVis_Large_PPC_mean));
    dVis_Large_PPC_sem          = squeeze(nanstd(dVis_PPC_reshape(:,2,:))) / sqrt(sum(~isnan(dVis_PPC_reshape(:,2,1)))); dVis_Large_PPC_sem = dVis_Large_PPC_sem(~isnan(dVis_Large_PPC_sem));
    
    dAud_Large_PPC_mean         = squeeze(nanmean(dAud_PPC_reshape(:,2,:))); dAud_Large_PPC_mean = dAud_Large_PPC_mean(~isnan(dAud_Large_PPC_mean));
    dAud_Large_PPC_sem          = squeeze(nanstd(dAud_PPC_reshape(:,2,:))) / sqrt(sum(~isnan(dAud_PPC_reshape(:,2,1)))); dAud_Large_PPC_sem = dAud_Large_PPC_sem(~isnan(dAud_Large_PPC_sem));
    
    % Make figure:
    figure; set(gcf,'color','w','units','normalized','Position', [0.1 0.1 .8 .8]); hold all;

    xpos        = 1:2;
    offset      = 2;
    errorbar(xpos, dVis_Small_V1_mean,dVis_Small_V1_sem,':ob','MarkerSize',20,'MarkerEdgeColor','blue','LineWidth',3);
    errorbar(xpos, dVis_Large_V1_mean,dVis_Large_V1_sem,'-ob','MarkerSize',20,'MarkerEdgeColor','blue','MarkerFaceColor','blue','LineWidth',3);
    xpos = xpos + offset;
    errorbar(xpos, dAud_Small_V1_mean,dAud_Small_V1_sem,':or','MarkerSize',20,'MarkerEdgeColor','red','LineWidth',3);
    errorbar(xpos, dAud_Large_V1_mean,dAud_Large_V1_sem,'-or','MarkerSize',20,'MarkerEdgeColor','red','MarkerFaceColor','red','LineWidth',3);
    
    xpos = xpos + offset;
    errorbar(xpos, dVis_Small_PPC_mean,dVis_Small_PPC_sem,':ob','MarkerSize',20,'MarkerEdgeColor','blue','LineWidth',3);
    errorbar(xpos, dVis_Large_PPC_mean,dVis_Large_PPC_sem,'-ob','MarkerSize',20,'MarkerEdgeColor','blue','MarkerFaceColor','blue','LineWidth',3);
    xpos = xpos + offset;
    errorbar(xpos, dAud_Small_PPC_mean,dAud_Small_PPC_sem,':or','MarkerSize',20,'MarkerEdgeColor','red','LineWidth',3);
    errorbar(xpos, dAud_Large_PPC_mean,dAud_Large_PPC_sem,'-or','MarkerSize',20,'MarkerEdgeColor','red','MarkerFaceColor','red','LineWidth',3);
    
    signrank(dVis_V1_reshape(~isnan(dVis_V1_reshape(:,1,1)),1,1),dVis_V1_reshape(~isnan(dVis_V1_reshape(:,1,2)),1,2))
    ranksum(dVis_V1_reshape(~isnan(dVis_V1_reshape(:,1,1)),1,1),dVis_V1_reshape(~isnan(dVis_V1_reshape(:,1,2)),1,2))
    ranksum(dVis_V1_reshape(~isnan(dVis_V1_reshape(:,2,1)),2,1),dVis_V1_reshape(~isnan(dVis_V1_reshape(:,2,2)),2,2))
%     errorbar(xpos + audoffset, dVis_V1_mean,dVis_V1_sem,':or','MarkerSize',20,'MarkerEdgeColor','red','LineWidth',3)
%     errorbar(xpos + audoffset, dVis_V1_mean,dVis_V1_sem,'-or','MarkerSize',20,'MarkerEdgeColor','red','MarkerFaceColor','red','LineWidth',3)
    
    ylim([-0.5 3])
    xlim([-.5 max(xpos)+0.5])
    ylabel('Dprime')
    
    XTickLabels = repmat({'Control' 'Opto'},1,4);
    
    set(gca,'XTick',1:max(xpos),'XTickLabels',XTickLabels);
    grid on;
    legend({'Small visual change' 'Large visual change' 'Small auditory change' 'Large auditory change'},'Location','north');
    legend boxoff
    
end


end