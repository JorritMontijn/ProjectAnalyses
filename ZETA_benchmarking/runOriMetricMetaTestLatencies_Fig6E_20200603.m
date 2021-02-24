clear all;
close all;
%strPath = 'F:\Data\Results\OriMetric\Data\';
strPath = 'F:\Data\Processed\ZETA\Latencies\';
strFigPath = 'F:\Data\Results\OriMetric\';

%% load data
sDir=dir([strPath 'ZetaDataPPL3Poisson*']);
intFiles=numel(sDir);
mapC = redbluepurple(intFiles);
sLoad1=load([strPath sDir(1).name]);
intNeurons=numel(sLoad1.cellInterpT);

vecJitter = nan(1,intFiles);
vecBaseRate = nan(1,intFiles);
vecResampNum = nan(1,intFiles);
matPeakT = nan(intNeurons,intFiles);
matPeakW = nan(intNeurons,intFiles);
intC = 0;
for intFile=1:intFiles
	strFile = sDir(intFile).name;
	%get params
	vecJitter(intFile) = str2double(getFlankedBy(strFile,'Jitter','Resamp'));
	vecResampNum(intFile) = str2double(getFlankedBy(strFile,'Resamp','.mat'));
	vecBaseRate(intFile) = str2double(getFlankedBy(strFile,'PeakRate','Jitter'));
	
	%get data
	sLoad=load([strPath strFile]);
	matPeakT(:,intFile) = sLoad.vecPeakT;
	matPeakW(:,intFile) = sLoad.vecPeakW;
end
vecJitter = vecJitter*2;
matPeakT = ((matPeakT-0.1)*1000);
matLowHigh = getCI(matPeakT,1,0.05);
vecMean = mean(matPeakT,1);
matPeakW = (matPeakW*1000);
matLowHighW = getCI(matPeakW,1,0.05);
vecMeanW = mean(matPeakW,1);

%% calculate data
vecUniqueBaseRates = unique(vecBaseRate);
vecUniqueJitters = unique(vecJitter);
matC = redbluepurple(numel(vecUniqueBaseRates));
vecEdgeR = [(vecUniqueBaseRates(1)-1) (vecUniqueBaseRates(1:end-1) + diff(vecUniqueBaseRates)/2) (vecUniqueBaseRates(end)+1)];
vecEdgeJ = [(vecUniqueJitters(1)-1) (vecUniqueJitters(1:end-1) + diff(vecUniqueJitters)/2) (vecUniqueJitters(end)+1)];
matJitter = repmat(vecJitter,[intNeurons 1]);
matBaseRates = repmat(vecBaseRate,[intNeurons 1]);
vecEstMu = mean(vecMean,1);
vecEstSd = std(matPeakT,[],1)/2;
vecEstRange = (matLowHigh(:,2) - matLowHigh(:,1));

vecWidMu = mean(vecMeanW,1);
vecWidSd = std(matPeakW,[],1)/2;
vecWidRange = (matLowHighW(:,2) - matLowHighW(:,1));

%% calculate correlations
[R1,P1] = corr(matBaseRates(:),matPeakT(:));
[R2,P2] = corr(vecBaseRate(:),vecEstMu(:));
[R3,P3] = corr(vecBaseRate(:),vecEstSd(:));
[R4,P4] = corr(matJitter(:),matPeakT(:));
[R5,P5] = corr(vecJitter(:),vecEstMu(:));
[R6,P6] = corr(vecJitter(:),vecEstSd(:));

[h ,crit_p, adj_p]=fdr_bh([P1 P2 P3 P4 P5 P6]);
adj_p = [P1 P2 P3 P4 P5 P6]*6;

%% plot peak time
figure
maxfig();
subplot(2,3,1)
scatter(matBaseRates(:),matPeakT(:))%,[],label2idx(matBaseRates(:)));
xlabel('Baseline spiking rate (Hz)');
set(gca,'xscale','log')
ylabel('Error in peak latency estimate (ms)');

title(sprintf('n=%d neurons, r=%.3f,p=%.3f',numel(matJitter),R1,adj_p(1)));
%colormap(gca,redbluepurple);
fixfig

subplot(2,3,2)
hold on
%plot([0 max(vecBaseRate(:))],[0 max(vecBaseRate(:))],'k--')
%plot([0 max(vecBaseRate(:))],-[0 max(vecBaseRate(:))],'k--')
scatter(vecBaseRate(:),vecEstMu(:))%,[],label2idx(vecJitter));
hold off
set(gca,'xscale','log')
xlabel('Baseline spiking rate (Hz)');
ylabel('Peak latency estimate mean (ms)');

title(sprintf('n=%d param sets, r=%.3f,p=%.3f',numel(vecJitter),R2,adj_p(2)));
%colormap(gca,parula);
fixfig

subplot(2,3,3)
hold on
%plot(sqrt(1/12)*[0 max(vecJitter(:))],sqrt(1/12)*[0 max(vecJitter(:))],'k--')
%plot(sqrt(1/12)*[0 max(vecJitter(:))],-sqrt(1/12)*[0 max(vecJitter(:))],'k--')
scatter(vecBaseRate(:),vecEstSd(:))%,[],label2idx(vecJitter));
hold off
set(gca,'xscale','log')
xlabel('Baseline spiking rate (Hz)');
ylabel('Peak latency estimate sd (ms)');
%colormap(gca,parula);

title(sprintf('n=%d param sets, r=%.3f,p=%.3f',numel(vecJitter),R3,adj_p(3)));
fixfig

% plot jitter
subplot(2,3,4)
hold on
plot([0 max(vecJitter(:))],[0 max(vecJitter(:))],'k--')
plot([0 max(vecJitter(:))],-[0 max(vecJitter(:))],'k--')
scatter(matJitter(:),matPeakT(:))%,[],label2idx(matJitter(:)));
hold off
xlabel('Peak latency jitter range (ms)');
ylabel('Error in peak latency estimate (ms)');

title(sprintf('n=%d neurons, r=%.3f,p=%.3f',numel(matJitter),R4,adj_p(4)));
%colormap(gca,parula);
fixfig

subplot(2,3,5)
hold on
plot(sqrt(1/12)*[0 max(vecJitter(:))],sqrt(1/12)*[0 max(vecJitter(:))],'k--')
plot(sqrt(1/12)*[0 max(vecJitter(:))],-sqrt(1/12)*[0 max(vecJitter(:))],'k--')
scatter(sqrt(1/12)*(vecJitter(:)),vecEstMu(:))%,[],label2idx(vecBaseRate));
hold off
xlabel('Peak latency jitter sd (ms)');
ylabel('Peak latency estimate mean (ms)');
%colormap(gca,redbluepurple);

title(sprintf('n=%d param sets, r=%.3f,p=%.3f',numel(vecJitter),R5,adj_p(5)));
fixfig

subplot(2,3,6)
hold on
plot(sqrt(1/12)*[0 max(vecJitter(:))],sqrt(1/12)*[0 max(vecJitter(:))],'k--')
scatter(sqrt(1/12)*(vecJitter(:)),vecEstSd(:))%,[],label2idx(vecBaseRate));
hold off
xlabel('Peak latency jitter sd (ms)');
ylabel('Peak latency estimate sd (ms)');
%colormap(gca,redbluepurple);
title(sprintf('n=%d param sets, r=%.3f,p=%.3f',numel(vecJitter),R6,adj_p(6)));
fixfig

drawnow;
export_fig(sprintf('%sLatencyBenchmarkFig.tif',strFigPath));
export_fig(sprintf('%sLatencyBenchmarkFigEF.pdf',strFigPath));
print(gcf,'-dpdf', sprintf('%sLatencyBenchmarkFig.pdf',strFigPath));

%% calculate correlations
[R1,P1] = corr(matBaseRates(:),matPeakW(:));
[R2,P2] = corr(vecBaseRate(:),vecWidMu(:));
[R3,P3] = corr(vecBaseRate(:),vecWidSd(:));
[R4,P4] = corr(matJitter(:),matPeakW(:));
[R5,P5] = corr(vecJitter(:),vecWidMu(:));
[R6,P6] = corr(vecJitter(:),vecWidSd(:));

[h ,crit_p, adj_p]=fdr_bh([P1 P2 P3 P4 P5 P6]);
adj_p = [P1 P2 P3 P4 P5 P6]*6;

%% plot peak width
figure
maxfig();
subplot(2,3,1)
scatter(matBaseRates(:),abs(matPeakW(:) - matJitter(:)))%,[],label2idx(matBaseRates(:)));
xlabel('Baseline spiking rate (Hz)');
set(gca,'xscale','log')
ylabel('Error in peak width estimate (ms)');

title(sprintf('n=%d neurons, r=%.3f,p=%.3f',numel(matJitter),R1,adj_p(1)));
%colormap(gca,redbluepurple);
fixfig

subplot(2,3,2)
hold on
%plot([0 max(vecBaseRate(:))],[0 max(vecBaseRate(:))],'k--')
%plot([0 max(vecBaseRate(:))],-[0 max(vecBaseRate(:))],'k--')
scatter(vecBaseRate(:),(vecWidMu(:)))%,[],label2idx(vecJitter));
hold off
set(gca,'xscale','log')
xlabel('Baseline spiking rate (Hz)');
ylabel('Peak width estimate mean (ms)');

title(sprintf('n=%d param sets, r=%.3f,p=%.3f',numel(vecJitter),R2,adj_p(2)));
%colormap(gca,parula);
fixfig

subplot(2,3,3)
hold on
%plot(sqrt(1/12)*[0 max(vecJitter(:))],sqrt(1/12)*[0 max(vecJitter(:))],'k--')
%plot(sqrt(1/12)*[0 max(vecJitter(:))],-sqrt(1/12)*[0 max(vecJitter(:))],'k--')
scatter(vecBaseRate(:),vecWidSd(:))%,[],label2idx(vecJitter));
hold off
set(gca,'xscale','log')
xlabel('Baseline spiking rate (Hz)');
ylabel('Peak width estimate sd (ms)');
%colormap(gca,parula);

title(sprintf('n=%d param sets, r=%.3f,p=%.3f',numel(vecJitter),R3,adj_p(3)));
fixfig

% plot jitter
subplot(2,3,4)
hold on
plot([0 max(vecJitter(:))],[0 max(vecJitter(:))],'k--')
%plot([0 max(vecJitter(:))],-[0 max(vecJitter(:))],'k--')
scatter(matJitter(:),(matPeakW(:)))%,[],label2idx(matJitter(:)));
hold off
xlabel('Real peak width (ms)');
ylabel('Estimated peak width (ms)');

title(sprintf('n=%d neurons, r=%.3f,p=%.3f',numel(matJitter),R4,adj_p(4)));
%colormap(gca,parula);
fixfig

subplot(2,3,5)
hold on
plot([0 max(vecJitter(:))],[0 max(vecJitter(:))],'k--')
%plot(sqrt(1/12)*[0 max(vecJitter(:))],-sqrt(1/12)*[0 max(vecJitter(:))],'k--')
scatter((vecJitter(:)),vecWidMu(:))%,[],label2idx(vecBaseRate));
hold off
xlabel('Mean real peak width (ms)');
ylabel('Mean estimated peak width (ms)');
%colormap(gca,redbluepurple);

title(sprintf('n=%d param sets, r=%.3f,p=%.3f',numel(vecJitter),R5,adj_p(5)));
fixfig

drawnow;
export_fig(sprintf('%sLatencyBenchmarkWidthFig.tif',strFigPath));
export_fig(sprintf('%sLatencyBenchmarkWidthFigEF.pdf',strFigPath));
print(gcf,'-dpdf', sprintf('%sLatencyBenchmarkWidthFig.pdf',strFigPath));

%% plot summary
%prep vars
vecRates = unique(matBaseRates);
vecRateSteps = diff(vecRates)/2;
vecRateEdges = [(vecRates(1) - vecRateSteps(1));(vecRates(1:(end-1)) + vecRateSteps);(vecRates(end) + vecRateSteps(end))];

vecJitters = [0;unique(matJitter)];
vecJitterSteps = diff(vecJitters)/2;
vecJitterEdges = [(vecJitters(1) - vecJitterSteps(1));(vecJitters(1:(end-1)) + vecJitterSteps);(vecJitters(end) + vecJitterSteps(end))];

vecPeakTs = (-12:12)';
vecPeakTSteps = diff(vecPeakTs)/2;
vecPeakTEdges = [(vecPeakTs(1) - vecPeakTSteps(1));(vecPeakTs(1:(end-1)) + vecPeakTSteps);(vecPeakTs(end) + vecPeakTSteps(end))];

vecPeakWs = (-12:12)';
vecPeakWSteps = diff(vecPeakWs)/2;
vecPeakWEdges = [(vecPeakWs(1) - vecPeakWSteps(1));(vecPeakWs(1:(end-1)) + vecPeakWSteps);(vecPeakWs(end) + vecPeakWSteps(end))];

%vecQuantiles = [normcdf(-1,0,1) normcdf(0,0,1) normcdf(1,0,1)];
vecQuantiles = [0.25 0.5 0.75];

%plot
figure
maxfig()
subplot(3,4,1)
%heat map peak / base
[n,x,y] = histcounts2(matBaseRates(:),matPeakT(:),vecRateEdges,vecPeakTEdges);%,[],label2idx(matBaseRates(:)));
imagesc(n');axis xy;
set(gca,'ytick',[3:5:25])
set(gca,'yticklabel',vecPeakTs(get(gca,'ytick')))
set(gca,'xtick',[1:4:numel(vecRates)])
set(gca,'xticklabel',vecRates(get(gca,'xtick')))
xlabel('Background rate (Hz)');
ylabel('Error in peak est. (ms)');
colormap(parula)
title('All jitters')
fixfig;grid off;

subplot(3,4,2)
%heat map peak / jitter
[n,x,y] = histcounts2(matJitter(:),matPeakT(:),vecJitterEdges,vecPeakTEdges);%,[],label2idx(matBaseRates(:)));
imagesc(n');axis xy;
hold on
plot([find(vecJitters==2) find(vecJitters==20)],[find(vecPeakTs==1) find(vecPeakTs==10)],'w--');
plot([find(vecJitters==2) find(vecJitters==20)],[find(vecPeakTs==-1) find(vecPeakTs==-10)],'w--');
hold off
set(gca,'ytick',[3:5:25])
set(gca,'yticklabel',vecPeakTs(get(gca,'ytick')))
set(gca,'xtick',vecJitters(1:5:numel(vecJitters))+1)
set(gca,'xticklabel',vecJitters(1:5:numel(vecJitters)));
xlabel('Jitter in peak (ms)');
ylabel('Error in peak est. (ms)');
title('All background rates')
fixfig;grid off;

subplot(3,4,3)
%heat map width / base
[n,x,y] = histcounts2(matBaseRates(:),matPeakW(:)-matJitter(:),vecRateEdges,vecPeakWEdges);%,[],label2idx(matBaseRates(:)));
imagesc(n');axis xy;
set(gca,'ytick',[3:5:25])
set(gca,'yticklabel',vecPeakWs(get(gca,'ytick')))
set(gca,'xtick',[1:4:numel(vecRates)])
set(gca,'xticklabel',vecRates(get(gca,'xtick')))
xlabel('Background rate (Hz)');
ylabel('Error in width est. (ms)');
colormap(parula)
title('All jitters')
fixfig;grid off;

subplot(3,4,4)
%heat map width / jitter
[n,x,y] = histcounts2(matJitter(:),matPeakW(:),vecJitterEdges,vecJitterEdges);%,[],label2idx(matBaseRates(:)));

imagesc(n');axis xy;
hold on
plot([find(vecJitters==2) find(vecJitters==20)],[find(vecJitters==2) find(vecJitters==20)],'w--');
hold off
set(gca,'xtick',vecJitters(1:5:numel(vecJitters))+1)
set(gca,'xticklabel',vecJitters(1:5:numel(vecJitters)));
set(gca,'ytick',vecJitters(1:5:numel(vecJitters))+1)
set(gca,'yticklabel',vecJitters(1:5:numel(vecJitters)));
xlabel('Jitter in peak (ms)');
ylabel('Width est. (ms)');
colormap(parula)
title('All background rates')
fixfig;grid off;

subplot(3,4,5)
%mean plot peak / base
[vecCounts,vecMeans,vecSDs,cellVals] = makeBins(matBaseRates(:),matPeakT(:),vecRateEdges);
matCIs = cell2mat(cellfun(@quantile,cellVals,cellfill(vecQuantiles,size(cellVals)),'uniformoutput',false)');
errorbar(1:numel(vecRates),matCIs(:,2),matCIs(:,2)-matCIs(:,1),matCIs(:,2)-matCIs(:,3))
xlim([0.5 numel(vecMeans)+0.5]);
ylim([min(vecPeakTs) max(vecPeakTs)]);
set(gca,'xtick',[1:4:numel(vecRates)])
set(gca,'xticklabel',vecRates(get(gca,'xtick')))
xlabel('Background rate (Hz)');
ylabel('Error in peak est. (ms)');
title('25th, median, 75th perc.')
fixfig;

%{
subplot(3,4,6)
%mean plot peak / jitter
[vecCounts,vecMeans,vecSDs] = makeBins(matJitter(:),matPeakT(:),vecJitterEdges);
plot([find(vecJitters==2) find(vecJitters==20)],[1 10],'k--');
hold on
plot([find(vecJitters==2) find(vecJitters==20)],[-1 -10],'k--');
errorbar(1:numel(vecMeans),vecMeans,vecSDs)
hold off
xlim([0.5 numel(vecMeans)+0.5]);
ylim([min(vecPeakTs) max(vecPeakTs)]);
set(gca,'xtick',[1:6:numel(vecJitters)])
set(gca,'xticklabel',vecJitters(get(gca,'xtick')))
xlabel('Peak jitter (ms)');
ylabel('Error in peak est. (ms)');
fixfig;
%}

subplot(3,4,6)
%mean plot peak / jitter
[vecCounts,vecMeans,vecSDs,cellVals] = makeBins(matJitter(:),matPeakT(:),vecJitterEdges);
matCIs = cell2mat(cellfun(@quantile,cellVals,cellfill(vecQuantiles,size(cellVals)),'uniformoutput',false)');
errorbar(vecJitters,matCIs(:,2),matCIs(:,2)-matCIs(:,1),matCIs(:,2)-matCIs(:,3))
hold on
plot([0 20],[0 10],'k--');
plot([0 20],[0 -10],'k--');
hold off
ylim([min(vecPeakTs) max(vecPeakTs)]);
set(gca,'xtick',[0 5 10 15 20])
xlabel('Jitter in peak (ms)');
ylabel('Error in peak est. (ms)');
title('25th, median, 75th perc.')
fixfig;


subplot(3,4,7)
%mean plot width / base
[vecCounts,vecMeans,vecSDs,cellVals] = makeBins(matBaseRates(:),matPeakW(:)-matJitter(:),vecRateEdges);
matCIs = cell2mat(cellfun(@quantile,cellVals,cellfill(vecQuantiles,size(cellVals)),'uniformoutput',false)');
errorbar(1:numel(vecMeans),matCIs(:,2),matCIs(:,2)-matCIs(:,1),matCIs(:,2)-matCIs(:,3))
xlim([0.5 numel(vecMeans)+0.5]);
ylim([min(vecPeakWs) max(vecPeakWs)]);
set(gca,'xtick',[1:4:numel(vecRates)])
set(gca,'xticklabel',vecRates(get(gca,'xtick')))
xlabel('Background rate (Hz)');
ylabel('Error in width est. (ms)');
title('25th, median, 75th perc.')
fixfig;

subplot(3,4,8)
%mean plot width / jitter
[vecCounts,vecMeans,vecSDs,cellVals] = makeBins(matJitter(:),matPeakW(:),vecJitterEdges);
matCIs = cell2mat(cellfun(@quantile,cellVals,cellfill(vecQuantiles,size(cellVals)),'uniformoutput',false)');

plot([1 20],[1 20],'k--');
hold on
errorbar(vecJitters,matCIs(:,2),matCIs(:,2)-matCIs(:,1),matCIs(:,2)-matCIs(:,3))
hold off
%xlim([0.5 numel(vecMeans)+0.5]);
%ylim([0 max(vecPeakWs)]);
set(gca,'xtick',[0 5 10 15 20])
%set(gca,'xticklabel',vecJitters(get(gca,'xtick')))
xlabel('Jitter in peak (ms)');
ylabel('Width est. (ms)');
title('25th, median, 75th perc.')
fixfig;


drawnow;
export_fig(sprintf('%sLatencyBenchmarkSummaryFig.tif',strFigPath));
export_fig(sprintf('%sLatencyBenchmarkSummaryFigEF.pdf',strFigPath));
print(gcf,'-dpdf', sprintf('%sLatencyBenchmarkSummaryFig.pdf',strFigPath));

%%
figure;colormap(parula);colorbar;
export_fig(sprintf('%sColorbarParula.tif',strFigPath));
export_fig(sprintf('%sColorbarParula.pdf',strFigPath));
print(gcf,'-dpdf', sprintf('%sColorbarParula.pdf',strFigPath));
