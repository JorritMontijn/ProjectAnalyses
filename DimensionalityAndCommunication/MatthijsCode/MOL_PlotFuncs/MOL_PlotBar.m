function MOL_PlotBar(plotydata,colors,labels)

dostatistics = 0;

%% Make the figure window:
barfig = figure;
set(barfig,'units','normalized','Position',[0.2 0 0.8 1],'color','w')

%% Size of input data:
[rows,columns] = size(plotydata);

%% Get averaged data per session
nConditions = 0;
for rowplot = 1:rows
    for columnplot = 1:columns
        if ~isempty(plotydata{rowplot,columnplot})
            nConditions = nConditions+1;
            sessionaverage = nanmean(plotydata{rowplot,columnplot},2);
%             sessionaverage = nanmedian(plotydata{rowplot,columnplot},2);
            sessionaveragedata(1:100,nConditions) = NaN;
            sessionaveragedata(1:numel(sessionaverage),nConditions) = sessionaverage;
            boxplotlabels{nConditions} = labels{rowplot,columnplot};
            boxplotcolors(nConditions,:) = colors{rowplot,columnplot};
        end
    end
end

% sessionaveragedata(sessionaveragedata>8000) = NaN;

%% Parametric or non-parametric data?
% Testing for normality
parametric = 1;
for nCond = 1:nConditions
    if kstest(sessionaveragedata(:,nConditions))
        parametric = 0;
    end
end
parametric = 1;

%% For parametric data
if parametric
    
    % Statistics on the data
    if dostatistics
        if nConditions==2 && parametric %T-test
            [~,p,~,stats] = ttest2(sessionaveragedata(:,1),sessionaveragedata(:,2));
            stats
            p
            if p < 0.05; disp('Happiness'); end
        elseif nConditions>2 && parametric
            %Testing for homogeneity of variance
            HomoVar = 1;
            p = vartestn(sessionaveragedata,'display','off')  ;
            if p < 0.05; HomoVar = 0; end
            %Test
            if HomoVar; [p,~,stats] = anova1(sessionaveragedata); %Anova
            else [p,~,stats] = kruskalwallis(sessionaveragedata); end %Kruskal Wallis
            %Post-hoc
            if p < 0.05
                multcompare(stats,'alpha',0.05,'CType','bonferroni','display','off')
            end
        end
    end
    
    %Average the data for the bars + errorbar
    for rowplot = 1:rows
        for columnplot = 1:columns
            [nSessions,~] = size(plotydata{rowplot,columnplot});
            barmeandata{rowplot,columnplot} = nanmean(nanmean(plotydata{rowplot,columnplot},2));
            if nSessions>1
                barstddevdata{rowplot,columnplot} = nanstd(nanmean(plotydata{rowplot,columnplot},2))/sqrt(nSessions)*2;   %Plot +-2SEM;
            else SessionWithData = plotydata{rowplot,columnplot}(~all(isnan(plotydata{rowplot,columnplot}),2),:);
                barstddevdata{rowplot,columnplot} =  nanstd(SessionWithData);
            end
        end
    end
    
    %Convert
    barmeandata = cell2mat(barmeandata);
    barstddevdata = cell2mat(barstddevdata);
    
    %Create offset to plot bars on right x location
    switch rows
        case 1
            rowoffset = 0;
        case 2
            rowoffset = [-0.2 0.2];
        case 3
            rowoffset = [-0.2667 0 0.2667];
        case 4
            rowoffset = [-0.3 -0.1 0.1 0.3];
    end
    
    figure(barfig);
    hold all;
    for rowplot = 1:rows
        for columnplot = 1:columns
            bar(columnplot+rowoffset(rowplot),barmeandata(rowplot,columnplot),0.8/rows,'EdgeColor',[0 0 0],'FaceColor',colors{rowplot,columnplot});
            if ~isnan(barstddevdata(rowplot,columnplot));      
                errorbar(columnplot+rowoffset(rowplot),barmeandata(rowplot,columnplot),barstddevdata(rowplot,columnplot),'k.');
            end
        end
    end
    
    % Get the tick labels right
    XTicks = [];
    XTickLabels = {};
    for rowplot = 1:rows
        for columnplot = 1:columns
            XTicks(end+1) = columnplot+rowoffset(rowplot);
            XTickLabels{end+1} = labels{rowplot,columnplot};
        end
    end
    [XTicks,sortindex] = sort(XTicks);
    XTickLabels = XTickLabels(sortindex);
    
    % Define limits of axes
    if nSessions>1
        if min(min(barmeandata))<0
            YLim = [min(min(barmeandata-barstddevdata))*1.1 max(max(barmeandata+barstddevdata))*1.1+0.4];
        else YLim = [0 max(max(barmeandata+barstddevdata))*1.1];
        end
    else    YLim = [0 max(max(barmeandata))*1.1];
    end
    XLim=[0.5 columns+0.5];
    
    %Make up:
    set(gca,'Color',[1 1 1],'YLim',YLim,'XLim',XLim,'XTick',XTicks','XTickLabels',XTickLabels);
    set(gca,'XTickLabelRotation',60)
    
else
    %% For non-parametric data
    
    %Statistics
    if dostatistics
        if nConditions==2 %Wilcoxon signed rank sum test
            [p,h,stats] = ranksum(sessionaveragedata(:,1),sessionaveragedata(:,2))
            if p < 0.05;
            end
            
        elseif nConditions>2 %Kruskal Wallis
            
            [p,~,stats] = kruskalwallis(sessionaveragedata);
            if p < 0.05
                multcompare(stats,'alpha',0.05,'CType','bonferroni');
            end
        end
    end
    
    %Plot the boxplots
    axes(gca); hold all;
    boxplot(sessionaveragedata,'colors',boxplotcolors,'labels',boxplotlabels,'plotstyle','compact');
end

end