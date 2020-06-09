function [aupar,vispar] = MOL_Psy_2Sided(varargin)
%% Get input arguments:
sessionData     = varargin{1};
trialData       = varargin{2};

%% General settings:
showIndFig          = 1;
newplot             = 0;
showResFig          = 0;

set(0,'defaultAxesFontSize',20)
sessioncounter      = 0;

%% Initialize structure for saving output fit parameters:
vispar  = struct();
aupar   = struct();

fulldata_viscorr        = NaN(10,10,10,10); %Init matrix for storing all data
fulldata_visincorr      = NaN(10,10,10,10); %Init matrix for storing all data
fulldata_aucorr         = NaN(10,10,10,10); %Init matrix for storing all data
fulldata_auincorr       = NaN(10,10,10,10); %Init matrix for storing all data

%% Loop over mice and then sessions:
mouseids = unique(sessionData.mousename)';
for mou = 1:length(mouseids)
    
    mouseid                 = mouseids{mou};
    sesselec                = unique(sessionData.session_ID(strcmp(sessionData.mousename,mouseid)))';
    
    %Give title based on sessiondata information:
    if showIndFig
        figure; set(gcf,'color','w','units','normalized','Position', [0.1 0.1 .8 .8]);
        if length(sesselec)==1
            if isfield(sessionData,'Date')
                h = suptitle(sprintf('Psychometric curves for %s - %s',mouseid,sessionData.Date{strcmp(sessionData.session_ID,sesselec)}));
            elseif isfield(sessionData,'Rec_datetime')
                h = suptitle(sprintf('Psychometric curves for %s - %s',mouseid,sessionData.Rec_datetime{strcmp(sessionData.session_ID,sesselec)}));
            end
        else
            h = suptitle(sprintf('Psychometric curves for %s - %d days',mouseid,length(sesselec)));
        end
        set(h,'FontSize',20,'FontWeight','normal')
    end
    
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
%         seshasopto = 0;
        if seshasopto && isfield(temptrialData,'hasphotostim') && any(temptrialData.hasphotostim==1)
            if strcmp(tempsessionData.OptoStimMode,'V1mixed')
                selectionopto(:,1) = ~(temptrialData.hasphotostim==1);
                selectionopto(:,2) = (temptrialData.hasphotostim==1) & (temptrialData.PostChangeOptoStart==tempsessionData.PostChangeOptoStart{1}(1));
                selectionopto(:,3) = (temptrialData.hasphotostim==1) & (temptrialData.PostChangeOptoStart==tempsessionData.PostChangeOptoStart{1}(2));
            else
                selectionopto(:,1) = ~(temptrialData.hasphotostim==1);
                selectionopto(:,2) = (temptrialData.hasphotostim==1);
            end
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
        if seshasprobe
            selectionprobe      = strcmp(temptrialData.trialType,probetrialtype);
            fprintf('%d Probe Trials\n',sum(selectionprobe))
        else selectionprobe = strcmp(temptrialData.trialType,'Z'); %Nonsensical index
        end
        
        %% Auditory:
        if seshasau
            for iOpto = 1:size(selectionopto,2)
                % Auditory trials:
                selectionau         = NaN(length(temptrialData.trialType),length(auconditions)+1);
                selectionau(:,1)    = selectionprobe & selectionopto(:,iOpto);
                for i = 1:length(auconditions)
                    %                 selectionau(:,i+1)    = strcmp(temptrialData.trialType,autrialtype) & ismember(-temptrialData.(autrialfield),auconditions(i));
                    selectionau(:,i+1)    = strcmp(temptrialData.trialType,autrialtype) & ismember(abs(temptrialData.(autrialfield)),auconditions(i)) & selectionopto(:,iOpto);
                    fprintf('%d Audio Trials (%d) %s\n',sum(selectionau(:,i+1)),auconditions(i),repmat('Opto',1,iOpto-1))
                end
                
                percsidecorrectauditory = NaN(1,size(selectionau,2));
                percsideincorrectauditory = NaN(1,size(selectionau,2));
                for i = 1:size(selectionau,2)
                    percsidecorrectauditory(i)    = sum(selectionau(:,i) & responseasauditory) / sum(selectionau(:,i));
                    percsideincorrectauditory(i)    = sum(selectionau(:,i) & responseasvisual) / sum(selectionau(:,i));
                end
                
                if aumultiply100; auconditions = auconditions*100; end
                
                if numel(auconditions)==1; audiolog = 0; end
                
                if audiolog %If logarithmic axis place probe -2 times interval on log scale
                    auprobepos = exp(min(log(auconditions)) - 2*mean(diff(log(auconditions))));
                else auprobepos = 0;
                end
                
                % Construct x from position of the probe trial and conditions
                xdata_au = [auprobepos auconditions];
                
                %% Psychometric fit:
                if numel(xdata_au)>=4 %If sufficient data points: fit a psychometric curve
                    if audiolog
                        xaulogdata = log10(xdata_au);
                        [mu, sd, gr, lr, aucurve] = MOL_FitCumGauss(xaulogdata,percsidecorrectauditory);
                        aucurve(:,1)  =   10.^aucurve(:,1);
                        mu = 10.^mu;
                        sd = 10.^sd;
                    else
                        [mu, sd, gr, lr, aucurve] = MOL_FitCumGauss(xdata_au,percsidecorrectauditory);
                    end
                    aupar.mu(mou,ses) = mu;
                    aupar.sd(mou,ses) = sd;
                    aupar.gr(mou,ses) = gr;
                    aupar.lr(mou,ses) = lr;
                end
                
                %% Figure
                if showIndFig
                    subplot(1,2,1); hold all;
                    optomarkers = {'.' 'o' '*'};
                    optomarkersizes = [30 10 10];
                    plot(xdata_au,percsidecorrectauditory,['r' optomarkers{iOpto}],'MarkerSize',optomarkersizes(iOpto));
                    plot(xdata_au,percsideincorrectauditory,['b' optomarkers{iOpto}],'MarkerSize',optomarkersizes(iOpto));
                    plot(xdata_au,percsideincorrectauditory,'b:','MarkerSize',15);
                    if exist('aucurve','var')
                        plot(aucurve(:,1),aucurve(:,2),sprintf('k%s',repmat('-',iOpto,1)),'LineWidth',3);
                    else
                        plot(xdata_au,percsidecorrectauditory,'r-','MarkerSize',30); %Replaced if psychofit
                    end
                    
                    %Figure Make up
                    xlim([auprobepos max(xdata_au)*1.02])
                    ylim([0 1])
                    xlabel(auxaxislabel,'FontSize', 20)
                    ylabel('% Response Auditory','FontSize', 20)
                    set(gca,'FontSize',15,'YAxisLocation','left','linewidth',2)
                    if audiolog;  set(gca,'XScale','log'); end
%                     if strcmp(tempsessionData.auChangeUnit,'Oct') 
                        XTickLabels = ['Probe' num2cell(xdata_au(2:end))];
%                     else
%                         XTickLabels = ['Probe' num2cell(xdata_au(2:end))];
%                     end
                    set(gca,'Xdir','reverse','XTick',xdata_au,'XTickLabels',XTickLabels);
                    box on
                end
                
                fulldata_aucorr(mou,ses,iOpto,1:numel(xdata_au))         = percsidecorrectauditory; %Store in matrix
                fulldata_auincorr(mou,ses,iOpto,1:numel(xdata_au))        = percsideincorrectauditory; %Store in matrix
            end
        end
        
        %% Visual:
        if seshasvis
            for iOpto = 1:size(selectionopto,2)
                selectionvis        = NaN(length(temptrialData.trialType),length(visconditions)+1);
                selectionvis(:,1)   = selectionprobe & selectionopto(:,iOpto);
                for i = 1:length(visconditions)
                    %                 selectionvis(:,i+1)    = strcmp(temptrialData.trialType,vistrialtype) & ismember((temptrialData.(vistrialfield)),visconditions(i));
                    selectionvis(:,i+1)    = strcmp(temptrialData.trialType,vistrialtype) & ismember(abs(temptrialData.(vistrialfield)),visconditions(i)) & selectionopto(:,iOpto);
                    fprintf('%d Visual Trials (%d) %s\n',sum(selectionvis(:,i+1)),visconditions(i),repmat('Opto',1,iOpto-1))
                end
                percsidecorrectvisual   = NaN(1,size(selectionvis,2));
                percsideincorrectvisual = NaN(1,size(selectionvis,2));
                for i = 1:size(selectionvis,2)
                    percsidecorrectvisual(i)    = sum(selectionvis(:,i) & responseasvisual ) / sum(selectionvis(:,i));
                    percsideincorrectvisual(i)    = sum(selectionvis(:,i) & responseasauditory) / sum(selectionvis(:,i));
                end
                
                if vismultiply100; visconditions = visconditions*100; end
                
                if numel(visconditions)==1; visuallog = 0; end
                
                if visuallog %If logarithmic axis place probe -2 times interval on log scale
                    visprobepos = exp(min(log(visconditions)) - 2*mean(diff(log(visconditions))));
                else visprobepos = 0;
                end
                
                % Construct x from position of the probe trial and conditions
                xdata_vis = [visprobepos visconditions];
                
                %% Psychometric fit:
                if numel(xdata_vis)>=4 %If sufficient data points: fit a psychometric curve
                    if visuallog
                        xvislogdata = log10(xdata_vis);
                        [mu, sd, gr, lr, viscurve] = MOL_FitCumGauss(xvislogdata,percsidecorrectvisual);
                        viscurve(:,1)  =   10.^viscurve(:,1);
                        mu = 10.^mu;
                        sd = 10.^sd;
                    else
                        [mu, sd, gr, lr, viscurve] = MOL_FitCumGauss(xdata_vis,percsidecorrectvisual);
                    end
                    vispar.mu(mou,ses) = mu;
                    vispar.sd(mou,ses) = sd;
                    vispar.gr(mou,ses) = gr;
                    vispar.lr(mou,ses) = lr;
                end
                
                %% Figure
                if showIndFig
                    subplot(1,2,2); hold all;
                    optomarkers = {'.' 'o' '*'};
                    optomarkersizes = [30 10 10];
                    plot(xdata_vis,percsidecorrectvisual,['b' optomarkers{iOpto}],'MarkerSize',optomarkersizes(iOpto));
                    plot(xdata_vis,percsideincorrectvisual,['r' optomarkers{iOpto}],'MarkerSize',optomarkersizes(iOpto));
                    plot(xdata_vis,percsideincorrectvisual,'r:','MarkerSize',15);
                    if exist('viscurve','var')
                        plot(viscurve(:,1),viscurve(:,2),sprintf('k%s',repmat('-',iOpto,1)),'LineWidth',3);
                    else
                        plot(xdata_vis,percsidecorrectvisual,'b-','MarkerSize',30); %Replaced if psychofit
                    end
                    % Make up
                    xlim([visprobepos max(xdata_vis)*1.02])
                    ylim([0 1])
                    xlabel(visxaxislabel,'FontSize', 20)
                    ylabel('% Response Visual','FontSize', 20)
                    set(gca,'FontSize',15,'YAxisLocation','right')
                    set(gca,'linewidth',2)
                    if visuallog; set(gca,'XScale','log'); end
                    XTickLabels = ['Probe' num2cell(round(xdata_vis(2:end)))];
                    set(gca,'Xdir','normal','XTick',xdata_vis,'XTickLabels',XTickLabels);
                    box on
                end
                
                fulldata_viscorr(mou,ses,iOpto,1:numel(xdata_vis))        = percsidecorrectvisual; %Store in matrix
                fulldata_visincorr(mou,ses,iOpto,1:numel(xdata_vis))      = percsideincorrectvisual; %Store in matrix
            end
        end
    end
    
    %% Generate the mean + curve
    if newplot
        figure; set(gcf,'color','w','units','normalized','Position', [0.1 0.1 .8 .8]);
        
        %% Auditory fit parameters, take mean value:
        F=@(mu,sd,gr,lr,x) gr+(1-gr-lr)*0.5*(1+erf((x-mu)/sqrt(2*sd^2)));
        mu = mean(aupar.mu);
        sd = mean(aupar.sd);
        gr = mean(aupar.gr);
        lr = mean(aupar.lr);
        
        %% Show:
        subplot(1,2,1); hold all;
        %  Generate curve from fit
        if audiolog;
            fineX = linspace(min(xdata_au),max(xdata_au),numel(xdata_au)*1000); % Create a new xAxis with higher resolution
            aucurve = feval(F,log10(mu),log10(sd),gr,lr,log10(fineX));
            %             aucurve = 10.^(feval(F,log10(mu),log10(sd),gr,lr,log10(fineX)));
        else
            fineX = linspace(min(xdata_au),max(xdata_au),numel(xdata_au)*1000); % Create a new xAxis with higher resolution
            aucurve = feval(F,mu,sd,gr,lr,fineX);
        end
        
        aucurve = [fineX', aucurve'];
        plot(aucurve(:,1),aucurve(:,2),'k-','LineWidth',3); hold on;
        if audiolog;  set(gca,'XScale','log'); end
        xlim([auprobepos max(xdata_au)*1.02])
        ylim([0 1])
        set(gca,'Xdir','reverse','XTick',xdata_au);
        box on;
        
        %% Show values to take at specified factors away from threshold:
        factors = [-2 -1 0 1 2]; %Times std away from threshold
        for f = 1:length(factors)
            if visuallog
                authrstd = 10.^(log10(mu) + log10(sd)*factors(f));
            else
                authrstd    = mu + sd*factors(f);
            end
            over        = find(aucurve(:,1)>authrstd);
            if ~isempty(over)
                plot(authrstd,aucurve(over(1),2),'r.','MarkerSize',40)
                if  strcmp(tempsessionData.auChangeUnit,'Oct')
                    text(authrstd,aucurve(over(1),2)+0.05,sprintf('%2.0f%%@%2.2f=%0.4f',aucurve(over(1),2)*100,authrstd,round(authrstd*256)))
                else
                    text(authrstd,aucurve(over(1),2)+0.05,sprintf('%2.0f%%@%2.2f=%2.0f',aucurve(over(1),2)*100,authrstd,round(authrstd,-1)))
                end
            end
        end
        
        %         for f = 1:length(xdata_psy_au)
        %             over        = find(aucurve(:,1)>xdata_psy_au(f)-0.01);
        %             plot(xdata_psy_au(f),aucurve(over(1),2),'r.','MarkerSize',40)
        %         end
        
        
        %% Visual fit parameters, take mean value:
        mu = mean(vispar.mu);
        sd = mean(vispar.sd);
        gr = mean(vispar.gr);
        lr = mean(vispar.lr);
        
        %% Show:
        subplot(1,2,2); hold all;
        %  Generate curve from fit
        if visuallog;
            fineX = linspace(min(xdata_vis),max(xdata_vis),numel(xdata_vis)*1000); % Create a new xAxis with higher resolution
            viscurve = feval(F,log10(mu),log10(sd),gr,lr,log10(fineX));
            %             aucurve = 10.^(feval(F,log10(mu),log10(sd),gr,lr,log10(fineX)));
        else
            fineX = linspace(min(xdata_vis),max(xdata_vis),numel(xdata_vis)*1000); % Create a new xAxis with higher resolution
            viscurve = feval(F,mu,sd,gr,lr,fineX);
        end
        
        viscurve = [fineX', viscurve'];
        plot(viscurve(:,1),viscurve(:,2),'k-','LineWidth',3); hold on;
        if visuallog;  set(gca,'XScale','log'); end
        xlim([visprobepos max(xdata_vis)*1.02])
        ylim([0 1])
        set(gca,'XTick',xdata_vis);
        box on;
        
        for f = 1:length(factors)
            if visuallog
                visthrstd = 10.^(log10(mu) + log10(sd)*factors(f));
            else visthrstd = mu + sd*factors(f);
            end
            over = find(viscurve(:,1)>visthrstd);
            if ~isempty(over)
                plot(visthrstd,viscurve(over(1),2),'b.','MarkerSize',40)
                text(visthrstd,viscurve(over(1),2)+0.05,sprintf('%2.0f%%@%2.2f=%2.0f',viscurve(over(1),2)*100,visthrstd,round(visthrstd)))
            end
        end
        
        %         for f = 1:length(xdata_psy_vis)
        %             over        = find(10.^viscurve(:,1)>xdata_psy_vis(f)-0.01);
        %             plot(xdata_psy_vis(f),viscurve(over(1),2),'b.','MarkerSize',40)
        %         end
        
    end
    
end


%%

% vispar.mu(vispar.mu==0 | vispar.mu==90) = NaN;
% aupar.mu(aupar.mu==0) = NaN;

%%

if showResFig
    f = figure; set(gcf,'color','w','units','normalized','Position', [0.1 0.4 .8 .4]);
    visdata     = vispar.mu;
    audata      = aupar.mu;
    
    subplot(1,2,1);
    plot(1:length(mouseids),visdata,'b.','MarkerSize',30); hold all;
    
    mVis = nanmean(visdata(:));
    stdVis = nanstd(visdata(:)) ;%/sqrt(sum(~isnan(dVis(:))));
    plot(length(mouseids)+1,mVis,'k.','MarkerSize',30)
    errorbar(length(mouseids)+1,mVis,stdVis,stdVis,'k','LineWidth',2)
    
    ylabel('Threshold Visual')
    set(gca, 'XTick', 1:length(mouseids),'XTickLabels', mouseids,'XTickLabelRotation',45)
    ylim([0 20])
    xlim([0.5 length(mouseids)+1.5])
    
    subplot(1,2,2);
    plot(1:length(mouseids),audata,'r.','MarkerSize',30); hold all;
    
    mAud = nanmean(audata(:));
    stdAud = nanstd(audata(:)); %/sqrt(sum(~isnan(dAud(:))));
    plot(length(mouseids)+1,mAud,'k.','MarkerSize',30)
    errorbar(length(mouseids)+1,nanmean(audata(:)),stdAud,stdAud,'k','LineWidth',2)
    
    ylabel('Threshold Auditory')
    set(gca, 'XTick', 1:length(mouseids),'XTickLabels', mouseids,'XTickLabelRotation',45)
    ylim([0 40])
    xlim([0.5 length(mouseids)+1.5])
    
end


end