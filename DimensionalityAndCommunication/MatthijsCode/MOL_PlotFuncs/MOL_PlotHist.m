function MOL_PlotHist(plotydata,plotxaxisdata,colors,labels,Plot_type)

set(0,'defaultAxesFontSize',20)

%% Parameters:
smoothing = 5;      %Smooth per session over n adjacent data points with gaussian window
resolution = 20;    %With large range of values, what is resolution of binning?
                   %One bin per n data points in one condition, minimum of 5, maximum 50;
%Define minimum and maximum of data to be plotted, using quantiles
%Change quantiles used to define limits of histogram:
Xminquantile = 0.01;    %Set quantile for start x-axis, min = 0;
Xmaxquantile = 0.99;    %Set quantile for end x-axis, max = 1;

%% Make the figure window:
figure;
set(gcf,'units','normalized','Position',[0.2 0 0.8 1],'color','w'); hold all;

%% Get size of plotdata
if all(size(plotydata) == size(plotxaxisdata));
    [rows,columns]= size(plotydata); end

%% Define properties of the input data for the x-axis
switch Plot_type
    case 'XYDATA'
        DataToGetAxisProperties = plotxaxisdata;
    otherwise
        DataToGetAxisProperties = plotydata;
end

for nRow = 1:rows
    for nCol = 1:columns
        minquantiles{nRow,nCol} = min(nanmean(quantile(DataToGetAxisProperties{nRow,nCol},Xminquantile,2)));
        maxquantiles{nRow,nCol} = max(nanmean(quantile(DataToGetAxisProperties{nRow,nCol},Xmaxquantile,2)));
    end
end
minedge = min(min(cell2mat(minquantiles)))*0.9-0.2; %Use this multiplicative + error to be sure to select all relevant values
maxedge = max(max(cell2mat(maxquantiles)))*1.1+0.1;

%Determine number of unique values within the range
TotalUniqueValues = [];
for nRow = 1:rows
    for nCol = 1:columns
        NewUniqueValues = unique(DataToGetAxisProperties{nRow,nCol}(~isnan(DataToGetAxisProperties{nRow,nCol})));
        if iscolumn(NewUniqueValues); NewUniqueValues = NewUniqueValues'; end
        Uniquevalues{nRow,nCol} = NewUniqueValues(NewUniqueValues>minedge & NewUniqueValues<maxedge);
        TotalUniqueValues = unique([TotalUniqueValues NewUniqueValues]);
    end
end
nUniquevalues = length(Uniquevalues{1,1});

%Determine if input are integers
setswithvalues      = find(~cellfun(@isempty,DataToGetAxisProperties));
firstsetwithvalues  = setswithvalues(1,1);
firstvec            = DataToGetAxisProperties{1,firstsetwithvalues}(1,~isnan(DataToGetAxisProperties{1,firstsetwithvalues}(1,:)));
inputintegers       = all(firstvec==floor(firstvec));

%% Define the edges
for rowplot = 1:rows
    for columnplot = 1:columns
        if inputintegers && nUniquevalues <50
            histedges{rowplot,columnplot} = ceil(minedge)-0.5:1:ceil(maxedge)+0.5;
            plotxvalues{rowplot,columnplot} = histedges{rowplot,columnplot}(1:end-1)+0.5;
        elseif nUniquevalues<50
            if any(Uniquevalues{rowplot,columnplot})
                histedges{rowplot,columnplot}  = [Uniquevalues{rowplot,columnplot}(1) Uniquevalues{rowplot,columnplot}(2:end)-diff(Uniquevalues{rowplot,columnplot})/2 Uniquevalues{rowplot,columnplot}(end)];
                plotxvalues{rowplot,columnplot} = Uniquevalues{rowplot,columnplot};
            else
                histedges{rowplot,columnplot} = NaN;
                plotxvalues{rowplot,columnplot} = NaN;
            end
        else    %Define number of edges, more data = more edges
            n_edges = floor(numel(DataToGetAxisProperties{1}(~isnan(DataToGetAxisProperties{1})))/resolution);
            if n_edges < 5; n_edges=5; elseif n_edges>50; n_edges = 50;  end
            
            %Make edges
            histedges{rowplot,columnplot}=linspace(minedge,maxedge,n_edges);
            
            %Adapt x position of edges for plotting
            edges_d = mean(diff(histedges{rowplot,columnplot}));
            plotxvalues{rowplot,columnplot} = histedges{rowplot,columnplot}(1:end-1)+edges_d/2;
        end
    end
end

%% Plot data
YDataToPlot_max = NaN(rows,columns);
for rowplot = 1:rows
    for columnplot = 1:columns
        [nSessions,~] = size(plotydata{rowplot,columnplot});
        if nSessions
            currenthistedges = histedges{rowplot,columnplot};
            currentplotxvalues = plotxvalues{rowplot,columnplot};
            switch Plot_type
                case 'HIST'
                    for nSession = 1:nSessions
                        YDataToPlot{rowplot,columnplot}(nSession,:) = histcounts(plotydata{rowplot,columnplot}(nSession,:),currenthistedges, 'Normalization','probability');
                    end
                case 'CUM'
                    for nSession = 1:nSessions
                        YDataToPlot{rowplot,columnplot}(nSession,:) = histcounts(plotydata{rowplot,columnplot}(nSession,:),currenthistedges, 'Normalization','cdf');
                    end
                case 'HAZARD'
                    for nSession = 1:nSessions
                        hazarddata  = histcounts(plotydata{rowplot,columnplot}(nSession,:),currenthistedges, 'Normalization','probability');
                        for edg = 1:length(hazarddata)
                            hazardplotdata(edg) = hazarddata(edg)/sum(hazarddata(edg:end));
                        end
                        YDataToPlot{rowplot,columnplot}(nSession,:) = hazardplotdata;
                    end
                case 'XYDATA'
                    YDataToPlot{rowplot,columnplot} = NaN(nSessions,length(currentplotxvalues));
                    for nSession = 1:nSessions
                        for histbin = 1:length(currenthistedges)-1
                            index = plotxaxisdata{rowplot,columnplot}(nSession,:) > currenthistedges(histbin) & plotxaxisdata{rowplot,columnplot}(nSession,:) <= currenthistedges(histbin+1);
                            YDataToPlot{rowplot,columnplot}(nSession,histbin) = mean(plotydata{rowplot,columnplot}(nSession,index));
                        end
                    end
            end
            if smoothing
                for nSession = 1:nSessions
                    win=gausswin(smoothing); %convolution with gaussian
                    win=win/sum(win); %normalized
                    smoothsession          = padarray(YDataToPlot{rowplot,columnplot}(nSession,:),[0 round(length(win)/2)],'symmetric','both'); %pad the array on both sides for convolution
                    smoothsession          = conv(smoothsession,win,'valid'); %Take only the valid overlapping center of convolution
                    smoothsession          = smoothsession(1:length(YDataToPlot{rowplot,columnplot}(nSession,:))); %slight correction to get same size (edge vs convolution)
                    
                    YDataToPlot{rowplot,columnplot}(nSession,:) = smoothsession;
                end
            end
            
            YDataToPlot_average{rowplot,columnplot} = nanmean(YDataToPlot{rowplot,columnplot},1);
            YDataToPlot_max(rowplot,columnplot) = max(max(YDataToPlot_average{rowplot,columnplot}));
        end
        if nSessions == 1
            LineHandles{rowplot,columnplot} = plot(currentplotxvalues,YDataToPlot_average{rowplot,columnplot},'LineWidth',2,'Color',colors{rowplot,columnplot});
        elseif nSessions>1
            YDataToPlot_SEM{rowplot,columnplot} = nanstd(YDataToPlot{rowplot,columnplot},0,1)/sqrt(nSessions);   %Plot +-2SEM
            CurrLineHandle = shadedErrorBar(currentplotxvalues,YDataToPlot_average{rowplot,columnplot},YDataToPlot_SEM{rowplot,columnplot},{'-','LineWidth',1,'Color',colors{rowplot,columnplot}},1);
            LineHandles{rowplot,columnplot} = CurrLineHandle.mainLine;
        end
    end
end

%% Set limits of axes
maxalldata = max(max(YDataToPlot_max));
YLim = [0 maxalldata*1.15];
plotxmin = min(min(cellfun(@min,plotxvalues)));
plotxmax = max(max(cellfun(@max,plotxvalues)));
XLim=[plotxmin plotxmax];
set(gca,'Color',[1 1 1],'YLim',YLim,'XLim',XLim);

%% Add Legend
PlotLineHandles = [];
PlotLabels ={};
for rowplot = 1:rows
    for columnplot = 1:columns
        if isobject(LineHandles{rowplot,columnplot})
            PlotLineHandles(end+1) = LineHandles{rowplot,columnplot};
            PlotLabels{end+1} = labels{rowplot,columnplot};
        end
    end
end
legend(gca,PlotLineHandles,PlotLabels); %,'Parent',plot_fig_handles_right.subplot(subplot))
legend('boxoff');

end