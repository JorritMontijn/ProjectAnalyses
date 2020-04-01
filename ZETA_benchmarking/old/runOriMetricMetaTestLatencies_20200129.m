clear all;
close all;
strPath = 'D:\Data\Results\OriMetric\Data\';
strFigPath = 'D:\Data\Results\OriMetric\';

%% load data
sDir=dir([strPath 'ZetaDataPPL2Poisson*']);
intFiles=numel(sDir);
mapC = redbluepurple(intFiles);
sLoad1=load([strPath sDir(1).name]);
intNeurons=numel(sLoad1.cellInterpT);

vecJitter = nan(1,intFiles);
vecBaseRate = nan(1,intFiles);
vecResampNum = nan(1,intFiles);
matPeakT = nan(intNeurons,intFiles);
intC = 0;
for intFile=1:intFiles
	strFile = sDir(intFile).name;
	%get params
	vecJitter(intFile) = str2double(getFlankedBy(strFile,'Jitter','Resamp'));
	vecResampNum(intFile) = str2double(getFlankedBy(strFile,'Resamp','.mat'));
	vecBaseRate(intFile) = str2double(getFlankedBy(strFile,'PeakRate','Jitter'));
	
	%get data
	sLoad=load([strPath strFile]);
	[a,vecMaxT]=cellfun(@max,sLoad.cellDeriv);
	for intN=1:numel(sLoad.cellInterpT)
		matPeakT(intN,intFile) = sLoad.cellInterpT{intN}(vecMaxT(intN));
	end
end
matLowHigh = getCI(matPeakT,1,0.05);
vecMean = mean(matPeakT,1);

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

%% calculate correlations
[R1,P1] = corr(matBaseRates(:),matPeakT(:));
[R2,P2] = corr(vecBaseRate(:),vecEstMu(:));
[R3,P3] = corr(vecBaseRate(:),vecEstSd(:));
[R4,P4] = corr(matJitter(:),matPeakT(:));
[R5,P5] = corr(vecJitter(:),vecEstMu(:));
[R6,P6] = corr(vecJitter(:),vecEstSd(:));

[h ,crit_p, adj_p]=fdr_bh([P1 P2 P3 P4 P5 P6]);
adj_p = [P1 P2 P3 P4 P5 P6]*6;

%% plot
figure
maxfig();
subplot(2,3,1)
scatter(matBaseRates(:),(matPeakT(:)-0.1)*1000)%,[],label2idx(matBaseRates(:)));
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
scatter(vecBaseRate(:),(vecEstMu(:)-0.1)*1000)%,[],label2idx(vecJitter));
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
scatter(vecBaseRate(:),vecEstSd(:)*1000)%,[],label2idx(vecJitter));
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
scatter(matJitter(:),(matPeakT(:)-0.1)*1000)%,[],label2idx(matJitter(:)));
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
scatter(sqrt(1/12)*(vecJitter(:)),(vecEstMu(:)-0.1)*1000)%,[],label2idx(vecBaseRate));
hold off
xlabel('Peak latency jitter sd (ms)');
ylabel('Peak latency estimate mean (ms)');
%colormap(gca,redbluepurple);

title(sprintf('n=%d param sets, r=%.3f,p=%.3f',numel(vecJitter),R5,adj_p(5)));
fixfig

subplot(2,3,6)
hold on
plot(sqrt(1/12)*[0 max(vecJitter(:))],sqrt(1/12)*[0 max(vecJitter(:))],'k--')
scatter(sqrt(1/12)*(vecJitter(:)),vecEstSd(:)*1000)%,[],label2idx(vecBaseRate));
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