function MOL_2ADC_Psy(varargin)

%% Get input arguments:
if nargin==2
    sessionData     = varargin{1};
    trialData       = varargin{2};
else
%     [Data] = MOL_GetData('E:','CHDET',{'BehaviorConflict' 'BehaviorPsychophysics'},{'2019' '2020' '2021' '2022' '2023'},[],{'sessionData' 'trialData'});
    [Data] = MOL_GetData('E:','CHDET',{'BehaviorConflict' 'BehaviorPsychophysics'},{'2003' '2004' '2007' '2009' '2010' '2011' '2012' '2013'},[],{'sessionData' 'trialData'});
    sessionData     = Data.sessionData;
    trialData       = Data.trialData;
    %% Remove last 20 trials:
    trialData = MOL_RemoveLastnTrials(trialData,20);
end

%% General settings:
Par.showIndSesFig          = 0;
Par.showExAnimalFig        = 0;
Par.showMeanAllFig         = 1;
Par.showStatsFig           = 1;
set(0,'defaultAxesFontSize',20)

%% Figure Settings dependent on au protocol:
switch sessionData.auChangeUnit{1}
    case 'Hz'
        Par.auprobepos = 0.5;
        Par.auticks = [10 50 100 1000 4000];
        Par.auticklabels = ['Probe' num2cell(Par.auticks)];
        Par.auxaxislabel  = 'Change in frequency (Hz)';
        Par.auystatslabel = 'Auditory threshold (Hz)';
    case 'Oct'
        Par.auprobepos = 0.001;
        Par.auticks = [1/256 1/64 1/8 1/2];
        Par.auticklabels = {'Probe' '1/256' '1/64' '1/8' '1/2'};
        Par.auxaxislabel  = 'Change in octave';
        Par.auystatslabel = 'Auditory threshold (partial octave)';
end

% Settings:
Par.visprobepos     = 0.5;
Par.visticks        = [2 5 15 30 90];
Par.vistickslabels  = ['Probe' num2cell(Par.visticks)];
Par.visxaxislabel   = 'Change in orientation (Degrees)';
Par.visystatslabel  = 'Visual threshold (Degrees)';

Par.yticks          = [0 0.25 0.5 0.75 1];

%% Init variables to store output data:
FullParam           = NaN(8,20,20); %Init matrix to store all param
FullParamError      = NaN(8,20,20); %Init matrix to store all param
% FullVisConditions   = NaN(8,20,20); %Init matrix to store all param
% FullAuConditions   = NaN(8,20,20); %Init matrix to store all param

%% Loop over mice and then sessions:
mouseids = unique(sessionData.mousename)';
for mou = 1:length(mouseids)
    %Get the relevant sessions for this mouse:
    sesselec        = unique(sessionData.session_ID(strcmp(sessionData.mousename,mouseids{mou})))';
    
    for ses = 1:length(sesselec) %Loop over sessions for this mouse
        sesid = sesselec(ses);                      %Get sessionID for this session to subselect trials and info about this sessions:
        [tempsessionData,temptrialData]         = MOL_getTempPerSes(sesid,sessionData,trialData);%Get the sessionData for each session individually
        
        [visconditions,auconditions,FullRespMat,~] = MOL_Psy_GetTrialMat(tempsessionData,temptrialData);
        
        %Rewrite fullresponse table to contingency table for model fitting (without conflict trials)
        ctable              = NaN(3,3,numel(visconditions));
        ctable(1,:,:)       = permute(FullRespMat(2:end,1,:),[3 1 2]);              %Auditory
        ctable(2,:,:)       = permute(FullRespMat(1,2:end,:),[1 3 2]);              %Visual
        ctable(3,:,:)       = repmat(permute(FullRespMat(1,1,:),[1 3 2]),1,1,numel(visconditions));    %Probe trials
        
        %[ctable,auconditions,visconditions]     = MOL_Prep_Psy2ADC(tempsessionData,temptrialData);
        
        %         [theta_est, theta_err, LLF, ctable_mod, ctable_fit] = mADC_model_fit(ctable);
        %% Align visual and auditory intensities
        %visual and auditory conditions should be normalized such that
        %value of 1 corresponds to expected asymptotic dmax
        meanconditions = mean([visconditions / max(visconditions); auconditions / max(auconditions)],1);
        
        %% Fit psychometric 2ADC model:
        fprintf('Fitting session %d/%d, of animal %d/%d\n',ses,length(sesselec),mou,length(mouseids));
        [theta_est, theta_err, LLF, ctable_mod, ctable_fit, sivals, psyc_mate, ce] = mADC_model_fit_psyc(ctable,meanconditions,[]); %#ok<ASGLU>
        
        %% Store parameters:
        %theta_est = 8 parameters: % % 3 d' -related parameters and one 'c' parameter for each location
        % 3 d' parameters are dmax, n and s50, for each location
        % dmax = theta(1:M); n = theta(M+1:2*M); s50 = theta(2*M+1:3*M); c = theta(3*M+1:end);
        theta_est(5) = theta_est(5) * max(auconditions);
        theta_est(6) = theta_est(6) * max(visconditions);
        theta_err(5) = theta_err(5) * max(auconditions);
        theta_err(6) = theta_err(6) * max(visconditions);
        
        FullParam(:,mou,ses)           = theta_est; %store parameters
        FullParamError(:,mou,ses)     = theta_err; %store parameters
        
        %% Figure:
        if Par.showIndSesFig
            %% Generate contingency table from fitted parameters:
            [xvals_fit_au,xvals_fit_vis,ctable_fit_mat] = MOL_Gen2ADC_PsyCurve(FullParam(:,mou,ses),Par);
            
            %% Show figure
            figure; set(gcf,'color','w','units','normalized','Position', [0.1 0.1 .8 .8]);
            
            %% Auditory part:
            subplot(1,2,1); hold all; 
            %Plot individual datapoints, auditory responses first:
            Par.aulinehandles(1) = plot([Par.auprobepos auconditions],[ctable(3,1,1); squeeze(ctable(1,1,:))],'r.','MarkerSize',30);
            Par.aulinehandles(2) = plot([Par.auprobepos auconditions],[ctable(3,2,1); squeeze(ctable(1,2,:))],'b.','MarkerSize',15);
            %Plot single session fit, auditory first:
            plot(xvals_fit_au,squeeze(ctable_fit_mat(1,1,:)),'k','LineWidth',3);
            plot(xvals_fit_au,squeeze(ctable_fit_mat(1,2,:)),'k:','LineWidth',3);
            
            %% Visual figure:
            subplot(1,2,2); hold all;
            %Plot individual datapoints, visual responses first:
            Par.vislinehandles(1) = plot([Par.visprobepos visconditions],[ctable(3,2,1); squeeze(ctable(2,2,:))],'b.','MarkerSize',30);
            Par.vislinehandles(2) = plot([Par.visprobepos visconditions],[ctable(3,1,1); squeeze(ctable(2,1,:))],'r.','MarkerSize',15);
            %Plot single session fit, visual first:
            plot(xvals_fit_vis,squeeze(ctable_fit_mat(2,2,:)),'k','LineWidth',3);
            plot(xvals_fit_vis,squeeze(ctable_fit_mat(2,1,:)),'k:','LineWidth',3);
            
            %% Common figure make up:
            MOL_Psy2Sided_FigMakeup(Par);
            suptitle(['Fitted 2-alternative detection psychometric model - Example session: ' sesid]);

        end
        
    end
end

%% Figure: mean fit curve for one animal:
if Par.showExAnimalFig
    example_animal = 1;
    
    figure; set(gcf,'color','w','units','normalized','Position', [0.1 0.1 .8 .8]);
    for iSes = 1:size(FullParam,3)
        if ~isnan(FullParam(1,example_animal,iSes))
            [xvals_fit_au,xvals_fit_vis,ctable_fit_mat] = MOL_Gen2ADC_PsyCurve(FullParam(:,example_animal,iSes),Par);

            %% Auditory figure:
            subplot(1,2,1); hold all;
            %Plot single session fit, auditory responses first:
            plot(xvals_fit_au,squeeze(ctable_fit_mat(1,1,:)),'Color',[0.6 0.6 0.6],'LineWidth',3);
            plot(xvals_fit_au,squeeze(ctable_fit_mat(1,2,:)),'Color',[0.6 0.6 0.6],'LineWidth',3);
            
            %% Visual figure:
            subplot(1,2,2); hold all;
            %Plot single session fit, visual responses first:
            plot(xvals_fit_vis,squeeze(ctable_fit_mat(2,1,:)),'Color',[0.6 0.6 0.6],'LineWidth',3);
            plot(xvals_fit_vis,squeeze(ctable_fit_mat(2,2,:)),'Color',[0.6 0.6 0.6],'LineWidth',3);
            
        end
    end
    
    %% Plot the mean psychometric curve:
    [xvals_fit_au,xvals_fit_vis,ctable_fit_mat] = MOL_Gen2ADC_PsyCurve(nanmean(FullParam(:,example_animal,:),3),Par);
                
    %% Auditory figure:
    subplot(1,2,1); hold all;
    %Plot mean fit for all sessions of this animal, auditory responses first:
    Par.aulinehandles(1) = plot(xvals_fit_au,squeeze(ctable_fit_mat(1,1,:)),'r','LineWidth',5);
    Par.aulinehandles(2) = plot(xvals_fit_au,squeeze(ctable_fit_mat(1,2,:)),'b:','LineWidth',5);
            
    %% Visual figure:
    subplot(1,2,2); hold all;
    %Plot mean fit for all sessions of this animal, visual responses first:
    Par.vislinehandles(1) = plot(xvals_fit_vis,squeeze(ctable_fit_mat(2,2,:)),'b','LineWidth',5);
    Par.vislinehandles(2) = plot(xvals_fit_vis,squeeze(ctable_fit_mat(2,1,:)),'r:','LineWidth',5);

    %% Common figure make up:
    MOL_Psy2Sided_FigMakeup(Par);
    suptitle(['Fitted 2-alternative detection psychometric model - Example animal: ' mouseids{example_animal}]);
    
end

%% Figure: mean fit curve for one animal:
if Par.showMeanAllFig
    figure; set(gcf,'color','w','units','normalized','Position', [0.1 0.1 .8 .8]);
    colors = lines(length(mouseids));
    aumouselinehandles = NaN(1,length(mouseids));
    for iMou = 1:length(mouseids)
        if ~isnan(FullParam(1,iMou,1))
            [xvals_fit_au,xvals_fit_vis,ctable_fit_mat] = MOL_Gen2ADC_PsyCurve(nanmean(FullParam(:,iMou,:),3),Par);
            
            %% Auditory figure:
            subplot(1,2,1); hold all;
            %Plot mean fit for each animal, auditory responses first:
            aumouselinehandles(iMou) = plot(xvals_fit_au,squeeze(ctable_fit_mat(1,1,:)),'Color',colors(iMou,:),'LineWidth',3);
            plot(xvals_fit_au,squeeze(ctable_fit_mat(1,2,:)),'Color',colors(iMou,:),'LineWidth',3);
            
            %% Visual figure:
            subplot(1,2,2); hold all;
            %Plot mean fit for each animal, visual responses first:
            plot(xvals_fit_vis,squeeze(ctable_fit_mat(2,1,:)),'Color',colors(iMou,:),'LineWidth',3);
            plot(xvals_fit_vis,squeeze(ctable_fit_mat(2,2,:)),'Color',colors(iMou,:),'LineWidth',3);
            
        end
    end
    
    %% Plot the overall mean psychometric curve:
    overall_theta = nanmean(nanmean(FullParam(:,:,:),3),2); %Average over mou and session dimensions:
    [xvals_fit_au,xvals_fit_vis,ctable_fit_mat] = MOL_Gen2ADC_PsyCurve(overall_theta,Par);
    
    %Auditory figure:
    subplot(1,2,1); hold all;
    plot(xvals_fit_au,squeeze(ctable_fit_mat(1,1,:)),'k','LineWidth',5);
    plot(xvals_fit_au,squeeze(ctable_fit_mat(1,2,:)),'k:','LineWidth',5);
    
    %Visual figure:
    subplot(1,2,2); hold all;
    plot(xvals_fit_vis,squeeze(ctable_fit_mat(2,1,:)),'k:','LineWidth',5);
    plot(xvals_fit_vis,squeeze(ctable_fit_mat(2,2,:)),'k','LineWidth',5);
    
    %% Put legend:
    subplot(1,2,1); hold all;
    legend(aumouselinehandles,mouseids,'FontSize',15,'Location','NorthEast');
    legend(gca,'boxoff');
    
    %% Common figure make up:
    if isfield(Par,'aulinehandles')
        Par = rmfield(Par,{'aulinehandles' 'vislinehandles'});
    end
    MOL_Psy2Sided_FigMakeup(Par);
    suptitle('Psychometric model fit - all animals');
    
end


%% 
if Par.showStatsFig
    %Trim FullParam to number of animals:
    FullParam = FullParam(:,1:length(mouseids),:);
    %Compute some values:
    FullParamReshape = reshape(FullParam,size(FullParam,1),size(FullParam,2)*size(FullParam,3));
    
    %%Checks and balances: 
    Param_animalmean    = nanmean(FullParam,3);
%     Param_animalsem     = nanstd(FullParam,[],3) ./ sqrt(sum(~isnan(FullParam(1,:)))); %not necessary?
    Param_fullmean      = nanmean(FullParamReshape,2);
    Param_fullsem       = nanstd(FullParamReshape,[],2) / sqrt(sum(~isnan(FullParamReshape(1,:))));
    
    %% Dprime figure
    figure; set(gcf,'color','w','units','normalized','Position', [0.1 0.2 .6 .5]);
    subplot(1,2,1);
    plot(1:length(mouseids),squeeze(FullParam(1,:,:)),'r.','MarkerSize',20); hold all;
    plot(1:length(mouseids),Param_animalmean(1,:),'k.','MarkerSize',45); hold all;
    errorbar(length(mouseids)+1,Param_fullmean(1),Param_fullsem(1),'-or','MarkerSize',10,'MarkerEdgeColor','red','MarkerFaceColor','red')

    ylabel('d-Prime Auditory')
    set(gca, 'XTick', 1:length(mouseids),'XTickLabels', mouseids,'XTickLabelRotation',45)
    ylim([0 4])
    xlim([0.5 length(mouseids)+1.5])
    
    subplot(1,2,2);
    plot(1:length(mouseids),squeeze(FullParam(2,:,:)),'b.','MarkerSize',20); hold all;
    plot(1:length(mouseids),Param_animalmean(2,:),'k.','MarkerSize',45); hold all;
    errorbar(length(mouseids)+1,Param_fullmean(2),Param_fullsem(2),'-ob','MarkerSize',10,'MarkerEdgeColor','blue','MarkerFaceColor','blue')

    ylabel('d-Prime Visual')
    set(gca, 'XTick', 1:length(mouseids),'XTickLabels', mouseids,'XTickLabelRotation',45)
    ylim([0 4])
    xlim([0.5 length(mouseids)+1.5])
    
    %% Threshold figure:
    figure; set(gcf,'color','w','units','normalized','Position', [0.1 0.2 .6 .5]);
    subplot(1,2,1);
    plot(1:length(mouseids),squeeze(FullParam(5,:,:)),'r.','MarkerSize',20); hold all;
    plot(1:length(mouseids),Param_animalmean(5,:),'k.','MarkerSize',45); hold all;
    errorbar(length(mouseids)+1,Param_fullmean(5),Param_fullsem(5),'-or','MarkerSize',10,'MarkerEdgeColor','red','MarkerFaceColor','red')

    ylabel(Par.auystatslabel)
    set(gca, 'XTick', 1:length(mouseids),'XTickLabels', mouseids,'XTickLabelRotation',45)
    switch sessionData.auChangeUnit{1}
        case 'Hz'
            set(gca, 'YScale', 'log')
            ylim([5 6000])
            set(gca,'YTick',[10 100 1000],'YTickLabels',[10 100 1000])
        case 'Oct'
            ylim([0 0.1])
    end
    xlim([0.5 length(mouseids)+1.5])
    
    subplot(1,2,2);
    plot(1:length(mouseids),squeeze(FullParam(6,:,:)),'b.','MarkerSize',20); hold all;
    plot(1:length(mouseids),Param_animalmean(6,:),'k.','MarkerSize',45); hold all;
    errorbar(length(mouseids)+1,Param_fullmean(6),Param_fullsem(6),'-ob','MarkerSize',10,'MarkerEdgeColor','blue','MarkerFaceColor','blue')

    ylabel(Par.visystatslabel)
    set(gca, 'XTick', 1:length(mouseids),'XTickLabels', mouseids,'XTickLabelRotation',45)
    ylim([0 20])
    xlim([0.5 length(mouseids)+1.5])
    
end
