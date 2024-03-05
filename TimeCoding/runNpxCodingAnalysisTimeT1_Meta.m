%% aim
%{
https://www.nature.com/articles/s41598-021-95037-z
https://www.nature.com/articles/s42003-021-02437-y

q1: Can absence of stimulus information in initial period be explained by neurons responding only to
black/white patches in the first x ms?

q2: How precise are spike times in the natural movie repetitions?


%}
%% define qualifying areas
clear all;
strUseArea = 'Primary visual area';
if isfolder('F:\Drive\PopTimeCoding') && isfolder('F:\Data\Processed\Neuropixels\')
	strDataPath = 'F:\Data\Processed\Neuropixels\';
	strFigurePathSR = 'F:\Drive\PopTimeCoding\single_recs';
	strFigurePath = 'F:\Drive\PopTimeCoding\figures\';
	strTargetDataPath = 'F:\Drive\PopTimeCoding\data\';
else
	strDataPath = 'E:\DataPreProcessed\';
	strFigurePathSR = 'C:\Drive\PopTimeCoding\single_recs';
	strFigurePath = 'C:\Drive\PopTimeCoding\figures\';
	strTargetDataPath = 'C:\Drive\PopTimeCoding\data\';
end

%% prep
sFiles = dir([strTargetDataPath 'T1Data*.mat']);
intRecNum = numel(sFiles);
sRec1=load(fullpath(sFiles(1).folder,sFiles(1).name));

%real/shufftid/poiss/poissgain
cellTypes = {'sReal','sShuffTid','sPoiss','sPoissGain'};

%ifr distro => var
strRec1 = sRec1.strRec;
vecIFR_bins = sRec1.sReal.vecBinsIFR;
matIFR_counts = nan(intRecNum,numel(vecIFR_bins),4);
vecIFR_means = nan(intRecNum,1,4);
vecIFR_median = nan(intRecNum,1,4);
vecIFR_sd = nan(intRecNum,1,4);

%% load
for intFile=1:numel(sFiles)
	intFile
	sLoad=load(fullpath(sFiles(intFile).folder,sFiles(intFile).name));
	
	for intType=1:4
		strType = cellTypes{intType};
		sData = sLoad.(strType);
		vecIFR_bins = sData.vecBinsIFR;
		matIFR_counts(intFile,:,intType) = sData.vecBinsCount;
		vecIFR_means(intFile,:,intType) = sData.dblMeanIFR;
		vecIFR_median(intFile,:,intType) = sData.dblMedianIFR;
		vecIFR_sd(intFile,:,intType) = sData.dblSdIFR;
	end
end
cellTypes = cellfun(@(x) x(2:end),cellTypes,'uniformoutput',false);

%% plot
hFig = figure;
dblBinStep = mean(diff(vecIFR_bins));
matIFR_normcounts1 = matIFR_counts(1,:,:) ./ sum(matIFR_counts(1,:,:),2);
matCol = lines(4);
maxfig;
subplot(2,4,1)
hold on
for intType=1:4
	plot(vecIFR_bins,matIFR_counts(1,:,intType),'color',matCol(intType,:));
end
hold off
legend(cellTypes,'location','best')
title(sprintf('%s',strRec1),'interpreter','none');
xlabel('Firing rate (IFR)');
ylabel('Count');
ylim([0 4e4]);

%{
%normalize by real
subplot(2,4,2)
hold on
vecSum=nan(1,4);
for intType=1:4
	vecIFR_normcounts = ((matIFR_counts(1,:,intType) ./ sum(matIFR_counts(1,:,intType),2))*vecIFR_sd(1,1,1))/dblBinStep;
	vecIFR_normbins = (vecIFR_bins - vecIFR_means(1,1,1))/vecIFR_sd(1,1,1);
	plot(vecIFR_normbins,vecIFR_normcounts,'color',matCol(intType,:));
	vecSum(intType) = sum(diff(vecIFR_normbins).*vecIFR_normcounts(2:end));
end
hold off
legend(cellTypes)
%}

% plot
vecIFR_interpbins = -4:0.1:4;
matIFR_interpcounts = zeros(intRecNum,numel(vecIFR_interpbins),4);
matSum = nan(4,intRecNum);
for intType=1:4
	subplot(2,4,4+intType);
	hold on
	for intRec=1:intRecNum
		vecIFR_normcounts = ((matIFR_counts(intRec,:,intType) ./ sum(matIFR_counts(intRec,:,intType),2))*vecIFR_sd(intRec,1,1))/dblBinStep;
		vecIFR_normbins = (vecIFR_bins - vecIFR_means(intRec,1,1))/vecIFR_sd(intRec,1,1);
		vecIFR_interpcounts = interp1(vecIFR_normbins,vecIFR_normcounts,vecIFR_interpbins,'linear',0);
		matIFR_interpcounts(intRec,:,intType) = vecIFR_interpcounts;
		plot(vecIFR_interpbins,vecIFR_interpcounts,'color',[0.5 0.5 0.5])
		matSum(intType,intRec) = sum(diff(vecIFR_interpbins).*vecIFR_interpcounts(2:end));
	end
	
	%plot mean
	vecMean = mean(matIFR_interpcounts(:,:,intType),1);
	plot(vecIFR_interpbins,vecMean,'color',matCol(intType,:));
	xlim([-4 4]);
	ylim([0 1]);
	hold off;
	title(cellTypes{intType});
	xlabel('Norm. firing rate');
	ylabel('Probability density');
end

%compare distros for means over recs
subplot(2,4,2);cla;
hold on;
for intType=1:4
	%plot mean
	vecMean = mean(matIFR_interpcounts(:,:,intType),1);
	plot(vecIFR_interpbins,vecMean,'color',matCol(intType,:));
end
hold off;
xlim([-4 4]);
ylim([0 1])
xlabel('Norm. firing rate');
ylabel('Probability density');
title('Mean over recs','interpreter','none');
cellTypesShort = {'R','ST','P','PG'};

%compare relative variability
subplot(2,4,3);
matSds = squeeze(vecIFR_sd);
plot(1:4,matSds','color',[0.5 0.5 0.5]);
hold on;
plot(1:4,mean(matSds,1),'k')
set(gca,'xtick',1:4,'xticklabel',cellTypes);
ylabel('Variance of pop firing rate');
strTests = '';
[h,vecP_nonorm]=ttest(matSds,repmat(matSds(:,1),[1 4]));
for intType=1:4
	errorbar(intType,mean(matSds(:,intType)),std(matSds(:,intType))./sqrt(intRecNum),'x','color',matCol(intType,:));
	if intType>1
		strTests = [strTests sprintf('; %s=%.1e',cellTypesShort{intType},vecP_nonorm(intType))];
	end
end
hold off
title([strTests(3:end)]);

%compare relative variability
subplot(2,4,4);
matNormSds = matSds./repmat(matSds(:,1),[1 4]);
matX = repmat((1:4),[15 1]);
vecX = matX(:);
vecY = matNormSds(:);
scatter(vecX,vecY,[],[0.5 0.5 0.5],'marker','.');
hold on
strTests = '';
[h,vecP]=ttest(matNormSds,repmat(matNormSds(:,1),[1 4]));
for intType=1:4
	errorbar(intType-0.1,mean(matNormSds(:,intType)),std(matNormSds(:,intType))./sqrt(intRecNum),'x','color',matCol(intType,:));
	if intType>1
		strTests = [strTests sprintf('; %s=%.1e',cellTypesShort{intType},vecP(intType))];
	end
end
hold off
set(gca,'xtick',1:4,'xticklabel',cellTypes);
ylabel('Variance of pop firing rate (norm)');
title(strTests(3:end));
ylim([0.4 1.2]);
xlim([0.75 4.25]);
fixfig;

%% save
export_fig(fullpath(strFigurePath,sprintf('T1_MetaIFR.tif')));
export_fig(fullpath(strFigurePath,sprintf('T1_MetaIFR.pdf')));
		