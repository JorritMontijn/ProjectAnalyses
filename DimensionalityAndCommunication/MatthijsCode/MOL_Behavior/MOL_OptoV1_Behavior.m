function [aupar,vispar] = MOL_OptoV1_Behavior(varargin)
%% Get input arguments:
if nargin==2
    sessionData     = varargin{1};
    trialData       = varargin{2};
else
    [Data] = MOL_GetData('K:','CHDET',{'ChangeDetectionConflict'},{'2019' '2020' '2021' '2022' '2023'},[],{'sessionData' 'trialData'});
    sessionData     = Data.sessionData;
    trialData       = Data.trialData;
    %% Remove last 20 trials:
    trialData = MOL_RemoveLastnTrials(trialData,20);
end

%% General settings:
Par.showMeanFig         = 1;
Par.showDprimeFig       = 1;
Par.normSessions        = 0;
Par.PhotostimArea       = 'V1';

Par.sumTrialsCondition  = 30;

set(0,'defaultAxesFontSize',20)
sessioncounter      = 0;

%% Initialize structure for saving output fit parameters:
vispar  = struct();
aupar   = struct();

fulldata_viscorr        = NaN(10,10,10,10); %Init matrix for storing all data
fulldata_visincorr      = NaN(10,10,10,10); %Init matrix for storing all data
fulldata_aucorr         = NaN(10,10,10,10); %Init matrix for storing all data
fulldata_auincorr       = NaN(10,10,10,10); %Init matrix for storing all data

dVis                    = NaN(10,10,10,10); %Init matrix for storing all data
dAud                    = NaN(10,10,10,10); %Init matrix for storing all data
cVis                    = NaN(10,10,10,10); %Init matrix for storing all data
cAud                    = NaN(10,10,10,10); %Init matrix for storing all data

%% Loop over mice and then sessions:
mouseids = unique(sessionData.mousename)';
for mou = 1:length(mouseids)
    
    mouseid                 = mouseids{mou};
    sesselec                = unique(sessionData.session_ID(strcmp(sessionData.mousename,mouseid)))';
    
    sesidx = strcmp(sessionData.mousename,mouseid) & sessionData.UseOpto==1 & strcmp(sessionData.PhotostimArea,Par.PhotostimArea);
    sesselec = unique(sessionData.session_ID(sesidx));
    
    for ses = 1:length(sesselec)
        sesid               = sesselec(ses);
        sessioncounter      = sessioncounter+1;         %Add one to the overall counter
        
        [tempsessionData,temptrialData] = MOL_getTempPerSes(sesid,sessionData,trialData);%Get the sessionData for each session individually:
        
        %% Get the trialtypes per session:
        trialtypes = unique(temptrialData.trialType);
        
        %% Define probe, auditory and visual trials and associated settings:
        seshasprobe = 0; seshasau = 0; seshasvis = 0;
        for tt = 1:length(trialtypes)
            switch trialtypes{tt}
                case 'P'
                    seshasprobe         = 1;
                    probetrialtype      = 'P';
                case 'Q'
                    seshasprobe         = 1;
                    probetrialtype      = 'Q';
                case 'R'
                    seshasprobe         = 1;
                    probetrialtype      = 'R';
                case 'V'
                    if strcmp(tempsessionData.strTask,'mol_det') || strcmp(tempsessionData.strTask,'mol_mod_discr')
                        seshasvis           = 1;
                        visuallog          	= 1;
                        vistrialtype        = 'V';
                        vismultiply100      = 1;
                        visxaxislabel       = 'Contrast (%)';
                        vissessionfield     = 'vecVisStimIntensities';
                        vistrialfield       = 'visualInt';
                    elseif strcmp(tempsessionData.strTask,'mol_ch_det')
                        seshasvis           = 1;
                        vismultiply100      = 0;
                        visuallog           = 1;
                        vistrialtype        = 'V';
                        vissessionfield     = 'vecOriChange';
                        vistrialfield       = 'visualOriChange';
                        visxaxislabel       = 'Change in orientation (degrees)';
                    end
                case 'X'
                    seshasvis           = 1;
                    vismultiply100      = 0;
                    visuallog           = 1;
                    vistrialtype        = 'X';
                    vissessionfield     = 'vecOriChange';
                    vistrialfield       = 'visualOriChange';
                    visxaxislabel       = 'Change in orientation (degrees)';
                case 'A'
                    if strcmp(tempsessionData.strTask,'mol_det') || strcmp(tempsessionData.strTask,'mol_mod_discr')
                        seshasau            = 1;
                        audiolog            = 0;
                        aumultiply100       = 1;
                        autrialtype         = 'A';
                        ausessionfield      = 'vecAuStimIntensities';
                        autrialfield        = 'audioInt';
                        auxaxislabel        = 'Sound level (dB)';
                    elseif strcmp(tempsessionData.strTask,'mol_ch_det')
                        seshasau            = 1;
                        aumultiply100       = 0;
                        audiolog            = 1;
                        autrialtype         = 'A';
                        if strcmp(tempsessionData.auChangeUnit,'Hz') || isempty(tempsessionData.auChangeUnit{1})
                            ausessionfield      = 'vecFreqChange';
                            autrialfield        = 'audioFreqChange';
                            auxaxislabel        = 'Change in frequency (Hz)';
                        else strcmp(tempsessionData.auChangeUnit,'Oct')
                            ausessionfield      = 'vecOctChange';
                            autrialfield        = 'audioOctChange';
                            auxaxislabel        = 'Change in octave (Oct)';
                        end
                    end
                case 'Y'
                    seshasau            = 1;
                    aumultiply100       = 0;
                    audiolog            = 1;
                    autrialtype         = 'Y';
                    if strcmp(tempsessionData.auChangeUnit,'Hz') || isempty(tempsessionData.auChangeUnit{1})
                        ausessionfield      = 'vecFreqChange';
                        autrialfield        = 'audioFreqChange';
                        auxaxislabel        = 'Change in frequency (Hz)';
                    else strcmp(tempsessionData.auChangeUnit,'Oct')
                        ausessionfield      = 'vecOctChange';
                        autrialfield        = 'audioOctChange';
                        auxaxislabel        = 'Change in octave (Oct)';
                    end
            end
        end
        
        if seshasvis
            if iscell(tempsessionData.(vissessionfield))
                visconditions   = tempsessionData.(vissessionfield){1};
                %                 visconditions    = unique(abs(temptrialData.(vistrialfield)))';
            else
                visconditions   = tempsessionData.(vissessionfield);
                %                 visconditions    = unique(abs(temptrialData.(vistrialfield)))';
            end
        end
        
        if seshasau
            if iscell(tempsessionData.(ausessionfield))
                auconditions    = tempsessionData.(ausessionfield){1};
                %                 auconditions    = unique(abs(temptrialData.(autrialfield)))';
            else
                auconditions    = tempsessionData.(ausessionfield);
                %                 auconditions    = unique(abs(temptrialData.(autrialfield)))';
            end
        end
        
        %% If session has opto:
        selectionopto = [];
        if isfield(tempsessionData,'UseOpto')
            seshasopto = tempsessionData.UseOpto;
        else         seshasopto = 0;
        end
        
        if seshasopto && isfield(temptrialData,'hasphotostim') && any(temptrialData.hasphotostim==1)
            selectionopto(:,1) = ~(temptrialData.hasphotostim==1);
            selectionopto(:,2) = (temptrialData.hasphotostim==1) & temptrialData.PostChangeOptoStart==0;
            selectionopto(:,3) = (temptrialData.hasphotostim==1) & temptrialData.PostChangeOptoStart==0.2;
        else selectionopto = ones(length(temptrialData.trialStart),1);
        end
        
        %% Selection of responses per trial:
        if sessionData.VisualLeftCorrectSide((strcmp(sessionData.mousename,mouseid)))
            responseasauditory  = strcmp(temptrialData.responseSide,'R');
            responseasvisual    = strcmp(temptrialData.responseSide,'L');
        else
            responseasauditory  = strcmp(temptrialData.responseSide,'L');
            responseasvisual    = strcmp(temptrialData.responseSide,'R');
        end
        
        %% Selection of probe trials:
        selectionprobe              = strcmp(temptrialData.trialType,probetrialtype);
        fprintf('%d Probe Trials\n',sum(selectionprobe))
        
        %% Auditory:
        % Auditory trials:
        selectionau                         = NaN(length(temptrialData.trialType),length(auconditions)+1);
        selectionau(:,1)                    = selectionprobe;
        for i = 1:length(auconditions)
            selectionau(:,i+1)        = strcmp(temptrialData.trialType,autrialtype) & ismember(abs(temptrialData.(autrialfield)),auconditions(i));
            fprintf('%d Audio Trials (%d)\n',sum(selectionau(:,i+1)),auconditions(i))
        end
        
        percsidecorrectauditory             = NaN(size(selectionau,2),size(selectionopto,2));
        percsideincorrectauditory           = NaN(size(selectionau,2),size(selectionopto,2));
        
        for iOpto = 1:size(selectionopto,2)
            for iAu = 1:size(selectionau,2)
                percsidecorrectauditory(iAu,iOpto)      = sum(selectionau(:,iAu) & selectionopto(:,iOpto) & responseasauditory) / sum(selectionau(:,iAu) & selectionopto(:,iOpto));
                percsideincorrectauditory(iAu,iOpto)    = sum(selectionau(:,iAu) & selectionopto(:,iOpto) & responseasvisual) / sum(selectionau(:,iAu) & selectionopto(:,iOpto));
            end
        end
        
        fulldata_aucorr(mou,ses,1:size(selectionopto,2),1:size(selectionau,2))         = percsidecorrectauditory; %Store in matrix
        fulldata_auincorr(mou,ses,1:size(selectionopto,2),1:size(selectionau,2))       = percsideincorrectauditory; %Store in matrix
        
        %% Visual:
        selectionvis                = NaN(length(temptrialData.trialType),length(visconditions)+1);
        selectionvis(:,1)           = selectionprobe;
        
        for i = 1:length(visconditions)
            selectionvis(:,i+1)    = strcmp(temptrialData.trialType,vistrialtype) & ismember(abs(temptrialData.(vistrialfield)),visconditions(i));
            fprintf('%d Visual Trials (%d)\n',sum(selectionvis(:,i+1)),visconditions(i))
        end
        percsidecorrectvisual   = NaN(size(selectionvis,2),size(selectionopto,2));
        percsideincorrectvisual = NaN(size(selectionvis,2),size(selectionopto,2));
        
        for iOpto = 1:size(selectionopto,2)
            for iVis = 1:size(selectionvis,2)
                percsidecorrectvisual(iVis,iOpto)    = sum(selectionvis(:,iVis) & selectionopto(:,iOpto) & responseasvisual ) / sum(selectionvis(:,iVis)  & selectionopto(:,iOpto));
                percsideincorrectvisual(iVis,iOpto)    = sum(selectionvis(:,iVis) & selectionopto(:,iOpto) & responseasauditory) / sum(selectionvis(:,iVis)  & selectionopto(:,iOpto));
            end
        end
        
        fulldata_viscorr(mou,ses,1:size(selectionopto,2),1:size(selectionvis,2))        = percsidecorrectvisual; %Store in matrix
        fulldata_visincorr(mou,ses,1:size(selectionopto,2),1:size(selectionvis,2))      = percsideincorrectvisual; %Store in matrix
        
        %% Compute d-prime for each condition:
        for iTrial = 1:2
            for iOpto = 1:3
                idx                 = (sum(selectionau(:,[1 iTrial+1],:),2) | sum(selectionvis(:,[1 iTrial+1],:),2)) & ...
                    selectionopto(:,iOpto);
                optotrialData       = struct();
                trialFields         = fieldnames(temptrialData);
                for iField = 1:length(trialFields)
                    optotrialData.(trialFields{iField}) = temptrialData.(trialFields{iField})(idx);
                end
                
                if any(selectionopto(:,iOpto)) && Par.showDprimeFig && sum(idx) > Par.sumTrialsCondition
                    [dVis(mou,ses,iTrial,iOpto),dAud(mou,ses,iTrial,iOpto),cVis(mou,ses,iTrial,iOpto),cAud(mou,ses,iTrial,iOpto)] = ...
                        MOL_Fit_2ADC_Full_Session(tempsessionData,optotrialData,0);
                end
            end
        end
    end
end

%% Organize data structure:
%Put all sessions below each other and remove singleton dimension:
if ndims(fulldata_aucorr) == 4
    fulldata_aucorr     = reshape(fulldata_aucorr,[size(fulldata_aucorr,1)*size(fulldata_aucorr,2) size(fulldata_aucorr,3) size(fulldata_aucorr,4)]);
    fulldata_auincorr   = reshape(fulldata_auincorr,[size(fulldata_auincorr,1)*size(fulldata_auincorr,2) size(fulldata_auincorr,3) size(fulldata_auincorr,4)]);
    
    fulldata_viscorr     = reshape(fulldata_viscorr,[size(fulldata_viscorr,1)*size(fulldata_viscorr,2) size(fulldata_viscorr,3) size(fulldata_viscorr,4)]);
    fulldata_visincorr   = reshape(fulldata_visincorr,[size(fulldata_visincorr,1)*size(fulldata_visincorr,2) size(fulldata_visincorr,3) size(fulldata_visincorr,4)]);
end

if Par.normSessions
    for ises = 1:size(fulldata_aucorr,1)
        fulldata_aucorr(ises,:,:) = fulldata_aucorr(ises,:,:) - repmat(fulldata_aucorr(ises,:,1),1,1,10);
        fulldata_auincorr(ises,:,:) = fulldata_auincorr(ises,:,:) - repmat(fulldata_auincorr(ises,:,1),1,1,10);
        fulldata_viscorr(ises,:,:) = fulldata_viscorr(ises,:,:) - repmat(fulldata_viscorr(ises,:,1),1,1,10);
        fulldata_visincorr(ises,:,:) = fulldata_visincorr(ises,:,:) - repmat(fulldata_visincorr(ises,:,1),1,1,10);
    end
end

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
    dVis_reshape = squeeze(reshape(dVis,size(dVis,1)*size(dVis,2),1,size(dVis,3),size(dVis,4)));
    dAud_reshape = squeeze(reshape(dAud,size(dAud,1)*size(dAud,2),1,size(dAud,3),size(dAud,4)));
    
    % Make figure:
    figure; set(gcf,'color','w','units','normalized','Position', [0.1 0.1 .8 .8]); hold all;
    y_mean      = squeeze(nanmean(dVis_reshape(:,1,:))); y_mean = y_mean(~isnan(y_mean));
    y_std       = squeeze(nanstd(dVis_reshape(:,1,:))) / sqrt(sum(~isnan(dVis_reshape(:,1,1)))); y_std = y_std(~isnan(y_std));
    xpos        = 1:length(y_mean);
    audoffset   = length(y_mean)+1;
    errorbar(xpos, y_mean,y_std,':ob','MarkerSize',20,'MarkerEdgeColor','blue','LineWidth',3);
    
    y_mean      = squeeze(nanmean(dVis_reshape(:,2,:)));  y_mean = y_mean(~isnan(y_mean));
    y_std       = squeeze(nanstd(dVis_reshape(:,2,:))) / sqrt(sum(~isnan(dVis_reshape(:,2,1)))); y_std = y_std(~isnan(y_std));
    errorbar(xpos, y_mean,y_std,'-ob','MarkerSize',20,'MarkerEdgeColor','blue','MarkerFaceColor','blue','LineWidth',3);
    
    y_mean      = squeeze(nanmean(dAud_reshape(:,1,:))); y_mean = y_mean(~isnan(y_mean));
    y_std       = squeeze(nanstd(dAud_reshape(:,1,:))) / sqrt(sum(~isnan(dAud_reshape(:,1,1))));  y_std = y_std(~isnan(y_std));
    errorbar(xpos + audoffset, y_mean,y_std,':or','MarkerSize',20,'MarkerEdgeColor','red','LineWidth',3)
    
    y_mean      = squeeze(nanmean(dAud_reshape(:,2,:))); y_mean = y_mean(~isnan(y_mean));
    y_std       = squeeze(nanstd(dAud_reshape(:,2,:))) / sqrt(sum(~isnan(dAud_reshape(:,2,1)))); y_std = y_std(~isnan(y_std));
    errorbar(xpos + audoffset, y_mean,y_std,'-or','MarkerSize',20,'MarkerEdgeColor','red','MarkerFaceColor','red','LineWidth',3)
    
    ylim([-0.5 3])
    ylabel('Dprime')
    
    switch Par.PhotostimArea
        case 'V1'
            XTickLabels = repmat({'Control' 'Full Inh' 'Late Inh'},1,2);
        case 'S1'
            XTickLabels = repmat({'Control' 'Full Inh'},1,2);
    end
    
    set(gca,'XTick',[xpos xpos+audoffset],'XTickLabels',XTickLabels);
    grid on;
    legend({'Small visual change' 'Large visual change' 'Small auditory change' 'Large auditory change'},'Location','north');
    legend boxoff
    
end


end