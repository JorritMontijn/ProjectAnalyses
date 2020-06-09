function [vis_par_series,au_par_series] = MOL_PsychoChDet(sessionData,trialData,spikeData)

set(0,'defaultAxesFontSize',20)

vis_par_series = [];
au_par_series = [];

vis_curve_series = {};
au_curve_series = {};

for mouseid = unique(sessionData.mousename)'
    
    figure; set(gcf,'color','w','units','normalized','Position', [0.1 0.1 .8 .8]);
    %Give title based on sessiondata information:
    sesselec = unique(sessionData.session_ID(strcmp(sessionData.mousename,mouseid)))';
    if length(sesselec)==1
        h = suptitle(sprintf('Psychometric curves for %s - %s',mouseid{1},sessionData.Date{1}));
    else
        h = suptitle(sprintf('Psychometric curves for %s - %d days',mouseid{1},length(sesselec)));
    end
    set(h,'FontSize',20,'FontWeight','normal')
    
    for ses = 1:length(sesselec)
        sesid = sesselec(ses);
        %Get the sessionData for each session individually:
        sesidx = sessionData.session_ID == sesid;
        sesfields = fieldnames(sessionData);
        tempsessionData = struct();
        for field = 1:length(sesfields)
            tempsessionData.(sesfields{field}) = sessionData.(sesfields{field})(sesidx);
        end
        
        %Get the trialData for each session individually:
        trialidx = trialData.session_ID == sesid;
        trialfields = fieldnames(trialData);
        temptrialData = struct();
        for field = 1:length(trialfields)
            temptrialData.(trialfields{field}) = trialData.(trialfields{field})(trialidx);
        end
        
        probetrialtype  = 'P';
        viscondfield    = 'visualOriChange';
        aucondfield     = 'audioFreqChange';
        vistrialtype    = 'X';
        autrialtype     = 'Y';
        visxaxislabel   = 'Change in orientation (degrees)';
        auxaxislabel    = 'Change in frequency (kHz)';
        
        % Probe trials:
        selectionprobe      = strcmp(temptrialData.trialType,probetrialtype);
        fprintf('%d Probe Trials\n',sum(selectionprobe))
        
        % Selection of responses trials:
        if sessionData.VisualLeftCorrectSide((strcmp(sessionData.mousename,mouseid)))
            responseasauditory = strcmp(temptrialData.responseSide,'R');
            responseasvisual = strcmp(temptrialData.responseSide,'L');
        else
            responseasauditory = strcmp(temptrialData.responseSide,'L');
            responseasvisual = strcmp(temptrialData.responseSide,'R');
        end
        
        %% Auditory:
        if isfield(tempsessionData,'vecFreqChange')
            if iscell(tempsessionData.vecFreqChange)
                auconditions    = tempsessionData.vecFreqChange{1};
            else
                auconditions    = tempsessionData.vecFreqChange;
            end
            
            % Auditory trials:
            selectionau = NaN(length(temptrialData.trialType),length(auconditions));
            for i = 1:length(auconditions)
                selectionau(:,i)    = strcmp(temptrialData.trialType,autrialtype) & ismember(abs(temptrialData.(aucondfield)),auconditions(i));
                fprintf('%d Audio Trials in cond %d\n',sum(selectionau(:,i)),i)
            end
            
            percsidecorrectproberight       = sum(selectionprobe & responseasauditory) / sum(selectionprobe);
            percsideincorrectproberight     = sum(selectionprobe & responseasvisual) / sum(selectionprobe);
            percsidecorrectright = [];
            percsideincorrectright = [];
            for i = 1:length(auconditions)
                percsidecorrectright(i)    = sum(selectionau(:,i) & responseasauditory) / sum(selectionau(:,i));
                percsideincorrectright(i)    = sum(selectionau(:,i) & responseasvisual) / sum(selectionau(:,i));
            end
            
            subplot(1,2,1)
            
            probepos = 1;

            %probepos = min(auconditions)*100-10;
            xdata_au = [probepos auconditions];
            ydata_au = [percsidecorrectproberight percsidecorrectright];
            
            plot(xdata_au,[percsidecorrectproberight percsidecorrectright],'r.','MarkerSize',30); hold on;
            plot(xdata_au,[percsidecorrectproberight percsidecorrectright],'r-','MarkerSize',30); hold on;
            plot(xdata_au,[percsideincorrectproberight percsideincorrectright],'b.','MarkerSize',15);hold on;
            plot(xdata_au,[percsideincorrectproberight percsideincorrectright],'b:','MarkerSize',15);
            
            if numel(xdata_au)>=4
                xdata_psy_au = xdata_au;
                [mu, stddev, guessrate, lapserate, curve] = MOL_FitCumGauss(xdata_au,ydata_au);
                au_par_series(1,ses) = mu;
                au_par_series(2,ses) = stddev;
                au_curve_series{ses} = curve;
                au_par_series(3,ses) = lapserate;
                au_par_series(4,ses) = guessrate;
%                              plot(curve(:,1),curve(:,2),'k-','LineWidth',3); hold on;
            else
                plot(xdata_au,[percsidecorrectproberight percsidecorrectright],'r.','MarkerSize',30); hold on;
            end
            
            % Makeup
            xlabel(auxaxislabel,'FontSize', 20)
            ylabel('% Response Auditory','FontSize', 20)
            set(gca,'FontSize',15,'YAxisLocation','left')
            set(gca,'linewidth',2)
            set(gca,'XScale','log')
            xlim([probepos max(xdata_au)])
            ylim([0 1])
            XTickLabels = ['Probe' num2cell(round(xdata_au(2:end)))];
            set(gca,'Xdir','reverse','XTick',xdata_au,'XTickLabels',XTickLabels);
            box on
        end
        
        %% Visual:
        if isfield(tempsessionData,'vecOriChange')
            if iscell(tempsessionData.vecOriChange)
                visconditions    = tempsessionData.vecOriChange{1};
            else visconditions    = tempsessionData.vecOriChange;
            end

            % Visual trials:
            selectionvis = NaN(length(temptrialData.trialType),length(visconditions));
            for i = 1:length(visconditions)
                selectionvis(:,i)    = strcmp(temptrialData.trialType,vistrialtype) & ismember(abs(temptrialData.(viscondfield)),visconditions(i));
                fprintf('%d Visual Trials in cond %d\n',sum(selectionvis(:,i)),i)
            end
            
            percsidecorrectprobeleft = sum(selectionprobe & responseasvisual ) / sum(selectionprobe);
            percsideincorrectprobeleft = sum(selectionprobe & responseasauditory) / sum(selectionprobe);
            percsidecorrectleft = [];
            percsideincorrectleft = [];
            for i = 1:length(visconditions)
                percsidecorrectleft(i)    = sum(selectionvis(:,i) & responseasvisual ) / sum(selectionvis(:,i));
                percsideincorrectleft(i)    = sum(selectionvis(:,i) & responseasauditory) / sum(selectionvis(:,i));
            end
            
            subplot(1,2,2); hold all;
            probepos = 0.1;
            xdata_vis = [probepos visconditions];
            
            plot(xdata_vis,[percsidecorrectprobeleft percsidecorrectleft],'b.','MarkerSize',30);
            plot(xdata_vis,[percsidecorrectprobeleft percsidecorrectleft],'b-','MarkerSize',30);
            plot(xdata_vis,[percsideincorrectprobeleft percsideincorrectleft],'r.','MarkerSize',15);
            plot(xdata_vis,[percsideincorrectprobeleft percsideincorrectleft],'r:','MarkerSize',15);
            
            xlogdata = log10([0.00001 visconditions]);
            
            ydata_vis = [percsidecorrectprobeleft percsidecorrectleft];
            
            if numel(xdata_vis)>=4
                xdata_psy_vis = xdata_vis;
                [mu, stddev, guessrate, lapserate, curve] = MOL_FitCumGauss(xdata_vis,ydata_vis);
%                 [mu, stddev, guessrate, lapserate, curve] = MOL_FitCumGauss(xlogdata,ydata_vis);
                vis_par_series(1,ses) = mu;
                vis_par_series(2,ses) = stddev;
                
                vis_par_series(3,ses) = lapserate;
                vis_par_series(4,ses) = guessrate;
                
                tempcurve = curve;
                tempcurve(:,1) = 10.^tempcurve(:,1);
                vis_curve_series{ses} = tempcurve;
                
%                 plot(10.^curve(:,1),curve(:,2),'k-','LineWidth',3); hold on;
%                 plot(curve(:,1),curve(:,2),'k-','LineWidth',3); hold on;
            else
                plot(xdata_vis,[percsidecorrectprobeleft percsidecorrectleft],'b.','MarkerSize',30);
            end
            
            xlabel(visxaxislabel,'FontSize', 20)
            ylabel('% Response Visual','FontSize', 20)
            set(gca,'FontSize',15,'YAxisLocation','right')
            set(gca,'linewidth',2)
            xlim([probepos max(xdata_vis)])
            ylim([0 1])
            set(gca,'XScale','log')
            XTickLabels = ['Probe' num2cell(round(xdata_vis(2:end)))];
            set(gca,'Xdir','normal','XTick',xdata_vis,'XTickLabels',XTickLabels);
            box on
        end
    end
    
    %% Generate the mean + curve
    F=@(g,l,u,v,x) g+(1-g-l)*0.5*(1+erf((x-u)/sqrt(2*v^2)));
    factors = [-2 -1 0 1 2]; %Times std away from threshold
    newplot = 0;
    if newplot
        %         figure; set(gcf,'color','w','units','normalized','Position', [0.1 0.1 .8 .8]);
        
        %% Auditory
        subplot(1,2,1); hold all;
        
        % Create a new xAxis with higher resolution
        fineX = linspace(min(xdata_psy_au),max(xdata_psy_au),numel(xdata_psy_au)*50);
        
        u = mean(au_par_series(1,:));
        v = mean(abs(au_par_series(2,:)));
        l = mean(au_par_series(3,:));
        g = mean(au_par_series(4,:));
        
        % Generate curve from fit
        aucurve = feval(F,g,l,u*100,v*100,fineX);
        aucurve = [fineX', aucurve'];
        plot(aucurve(:,1),aucurve(:,2),'k-','LineWidth',3); hold on;
        
        %         for f = 1:length(factors)
        %             authrstd    = (mean(au_par_series(1,:)) + mean(abs(au_par_series(2,:)))*factors(f))* 100;
        %             over        = find(aucurve(:,1)>authrstd);
        %             plot(authrstd,aucurve(over(1),2),'r.','MarkerSize',40)
        %         end
        
        for f = 1:length(xdata_psy_au)
            over        = find(aucurve(:,1)>xdata_psy_au(f)-0.01);
            plot(xdata_psy_au(f),aucurve(over(1),2),'r.','MarkerSize',40)
        end
        
        
        %% Visual
        % Create a new xAxis with higher resolution
        fineX = log10(linspace(min(xdata_psy_vis),max(xdata_psy_vis),numel(xdata_psy_vis)*50));
        
        subplot(1,2,2); hold all;
        
        u = mean(vis_par_series(1,:));
        v = mean(abs(vis_par_series(2,:)));
        l = mean(vis_par_series(3,:));
        g = mean(vis_par_series(4,:));
        
        % Generate curve from fit
        viscurve = feval(F,g,l,u,v,fineX);
        viscurve = [fineX', viscurve'];
        plot(10.^viscurve(:,1),viscurve(:,2),'k-','LineWidth',3); hold on;
        
        %         for f = 1:length(factors)
        %             visthrstd = 10^(mean(vis_par_series(1,:)) + mean(abs(vis_par_series(2,:)))*factors(f));
        %             over = find(10.^viscurve(:,1)>visthrstd);
        %             plot(visthrstd,viscurve(over(1),2),'b.','MarkerSize',40)
        %         end
        
        for f = 1:length(xdata_psy_vis)
            over        = find(10.^viscurve(:,1)>xdata_psy_vis(f)-0.01);
            plot(xdata_psy_vis(f),viscurve(over(1),2),'b.','MarkerSize',40)
        end
        
    end
    
end


end