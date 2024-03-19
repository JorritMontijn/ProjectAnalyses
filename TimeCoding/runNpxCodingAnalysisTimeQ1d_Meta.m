%% aim
%{
https://www.nature.com/articles/s41598-021-95037-z
https://www.nature.com/articles/s42003-021-02437-y

q1: Progression of orientation information over time: how does trial activity evolve? what is the function of the onset peak?
> is decoding better when matched for stimulus phase? => no.

q2: "Spike time and rate coding can be represented within a single model of spiking probability as a
function of time: rate codes are uniform over a certain period tau, while spike time codes are
temporally localized peaks"
> is this true?

q3: Rate codes do not exist; a rate code is simply a subset of spike time codes where the temporal
integration window is very large. But what about multi dim codes? Those are all rate based. Can we
formulate a multidimensional spike-time code? I.e., can we make a taxonomy of neural codes?

q4: How does information evolve over time, is initial peak indeed less tuned? Is pop activity
rhythmic? Are stimuli encoded invariant to brain state? Eg, high arousal, low arousal. Or is
stimulus manifold dynamic over time? Does manifold scale with arousal? => How does manifold depend
on binning size? What is the optimal time window?

%}
%% define qualifying areas
clear all;%close all;
cellUseAreas = {...
	'Primary visual area',...
	...'posteromedial visual area',...
	};
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

%% load data
strArea = 'V1';
cellTypes = {'Real','ShuffTid'};
cellSupraGranuInfra = {'Supragranular','Granular','Infragranular'};
for intLayer=1:numel(cellSupraGranuInfra)
for intType=1:numel(cellTypes)
for boolFixedSpikeGroupSize = false%[true false]
for dblOnset = [0 0.25]
strType = cellTypes{intType};
strLayer = cellSupraGranuInfra{intLayer};
intQuantiles = 10;
vecBinsDur = 0:0.01:1;
vecBinsDurC = vecBinsDur(2:end) - diff(vecBinsDur(1:2))/2;
vecBinsRate = 100:20:700;
vecBinsRateC = vecBinsRate(2:end) - diff(vecBinsRate(1:2))/2;
vecColGrey = [0.7 0.7 0.7];

strOnset = sprintf('%.0f',dblOnset*1000);
if boolFixedSpikeGroupSize
	sFiles = dir ([strTargetDataPath 'Q1dData*' strType '*Fixed*_Npx_DG_' strLayer '_Onset' strOnset '*']);
	strSGS = 'FixedSGS';
else
	sFiles = dir ([strTargetDataPath 'Q1dData*' strType '*Var*_Npx_DG_' strLayer '_Onset' strOnset '*']);
	strSGS = 'VarSGS';
end
[strPath,strName,strExt] = fileparts(sFiles(1).name);
cellSplit = strsplit(strName,'_');
strRecType = strjoin(cellSplit((end-2):end),',');
strRunType = [strType ',' strSGS ',' strRecType];
intRecNum = numel(sFiles);
cellQuantileDur = cell(intRecNum,intQuantiles);
cellQuantileConf = cell(intRecNum,intQuantiles);
matLatConf = nan(intRecNum,numel(vecBinsDurC));
vecSpikeGroupSize= nan(intRecNum,1);

cellQuantileRateChange = cell(intRecNum,intQuantiles);
cellQuantileConf_dR = cell(intRecNum,intQuantiles);
matRateConf = nan(intRecNum,numel(vecBinsRateC));

%% plot
figure;maxfig;
h1=subplot(2,3,1);
h2=subplot(2,3,2);
h3=subplot(2,3,3);
h4=subplot(2,3,4);
h5=subplot(2,3,5);
h6=subplot(2,3,6);


vecColChance = [0.3 0.3 0.3];
plot(h1,[0 1],([1 1]./12),'--','color',vecColChance);
hold(h1,'on')
title(h1,strType);

plot(h2,[0 100],([1 1]./12),'--','color',vecColChance);
hold(h2,'on')

hold(h3,'on')

plot(h4,[0 1],([1 1]./12),'--','color',vecColChance);
hold(h4,'on')
title(h4,strType);

plot(h5,[-300 300],([1 1]./12),'--','color',vecColChance);
hold(h5,'on')

hold(h6,'on')

for intFile=1:intRecNum
	%% load
	sLoad=load(fullpath(sFiles(intFile).folder,sFiles(intFile).name),'-mat');
	sSpikeGroup = sLoad.sSpikeGroup;
	strRec = sLoad.strRec;
	strType = sLoad.strType;
	vecOri180 = sLoad.vecOri180;
	vecOrientation = sLoad.vecOrientation;
	vecStimOffTime = sLoad.vecStimOffTime;
	vecStimOnTime = sLoad.vecStimOnTime;
	intSpikeGroupSize = sLoad.intSpikeGroupSize;
	
	%% create derived variables
	%put confidence in deciles per recording
	vecSpikeGroupDuration = [sSpikeGroup.Duration];
	vecSpikeGroupCorrect = [sSpikeGroup.Correct];
	vecSpikeGroupConfidence = [sSpikeGroup.Confidence];
	vecSpikeGroupLatency = dblOnset+[sSpikeGroup.Latency];
	vecSpikeGroupAvgRate = [sSpikeGroup.AvgRate];
	vecSpikeGroupRateChange = [sSpikeGroup.RateChange];
	
	%remove long durs
	%indRem = vecSortedDur > 0.05;
	%indRem = vecSortedLat < 0.15 | vecSortedDur > 0.05;
	indRem = [];
	vecSpikeGroupDuration(indRem) = [];
	vecSpikeGroupCorrect(indRem) = [];
	vecSpikeGroupConfidence(indRem) = [];
	vecSpikeGroupLatency(indRem) = [];
	vecSpikeGroupRateChange(indRem) = [];
	vecSpikeGroupAvgRate(indRem) = [];
	
	%sort
	[vecSortedDur,vecSort]=sort(vecSpikeGroupDuration);
	vecSortedCorr = vecSpikeGroupCorrect(vecSort);
	vecSortedConf = vecSpikeGroupConfidence(vecSort);
	vecSortedLat = vecSpikeGroupLatency(vecSort);
	
	%sort by dRate
	[vecSortedRateChange,vecSort]=sort(vecSpikeGroupRateChange);
	vecSortedCorr_dR = vecSpikeGroupCorrect(vecSort);
	vecSortedConf_dR = vecSpikeGroupConfidence(vecSort);
	vecSortedLat_dR = vecSpikeGroupLatency(vecSort);
	
	%calculate fraction correct and confidence per bin of equal size
	intBins = 10;
	intOriNum = numel(unique(vecOri180));
	intSperBin = floor(numel(vecSortedDur)/intBins);
	cellValsDur = cell(1,intBins);
	vecMeanDur = nan(1,intBins);
	vecSemDur = nan(1,intBins);
	vecMeanCorr = nan(1,intBins);
	matCiCorr = nan(2,intBins);
	cellValsCorr = cell(1,intBins);
	vecMeanConf = nan(1,intBins);
	vecSemConf = nan(1,intBins);
	cellValsConf = cell(1,intBins);
	
	cellVals_dR = cell(1,intBins);
	vecMean_dR = nan(1,intBins);
	vecSem_dR = nan(1,intBins);
	vecMeanConf_dR = nan(1,intBins);
	vecSemConf_dR = nan(1,intBins);
	cellValsConf_dR = cell(1,intBins);
	
	vecSampleGroup = zeros(size(vecSortedDur));
	for intBin=1:intBins
		intEndS = intSperBin*intBin;
		vecSamples = (intEndS-intSperBin+1):intEndS;
		
		%bin dur
		cellValsDur{intBin} = vecSortedDur(vecSamples);
		vecMeanDur(intBin) = mean(vecSortedDur(vecSamples));
		vecSemDur(intBin) = std(vecSortedDur(vecSamples))./sqrt(intSperBin);
		
		[phat,pci] = binofit(sum(vecSortedCorr(vecSamples)),intSperBin);
		vecMeanCorr(intBin) = phat;
		matCiCorr(:,intBin) = pci;
		cellValsCorr{intBin} = vecSortedCorr(vecSamples);
		
		vecMeanConf(intBin) = mean(vecSortedConf(vecSamples));
		vecSemConf(intBin) = std(vecSortedConf(vecSamples))./sqrt(intSperBin);
		cellValsConf{intBin} = vecSortedConf(vecSamples);
		vecSampleGroup(vecSamples) = intBin;
		
		%bin rate change
		cellVals_dR{intBin} = vecSortedRateChange(vecSamples);
		vecMean_dR(intBin) = mean(vecSortedRateChange(vecSamples));
		vecSem_dR(intBin) = std(vecSortedRateChange(vecSamples))./sqrt(intSperBin);

		vecMeanConf_dR(intBin) = mean(vecSortedConf_dR(vecSamples));
		vecSemConf_dR(intBin) = std(vecSortedConf_dR(vecSamples))./sqrt(intSperBin);
		cellValsConf_dR{intBin} = vecSortedConf_dR(vecSamples);
	end
	
	%conf with time
	[vecCounts,vecMeans,vecSDs]=makeBins(vecSpikeGroupLatency,vecSpikeGroupConfidence,vecBinsDur);
	plot(h1,vecBinsDurC,vecMeans,'color',vecColGrey);%,vecSDs./sqrt(vecCounts),vecSDs./sqrt(vecCounts),'color',lines(1));
	
	%conf with dur
	errorbar(h2,vecMeanDur*1000,vecMeanConf,vecSemConf,vecSemConf,vecSemDur,vecSemDur,'color',lines(1));
	
	%conf with rate
	[vecCountsR,vecMeansR,vecSDsR]=makeBins(vecSpikeGroupAvgRate,vecSpikeGroupConfidence,vecBinsRate);
	plot(h4,vecBinsRateC,vecMeansR,'color',vecColGrey);%,vecSDs./sqrt(vecCounts),vecSDs./sqrt(vecCounts),'color',lines(1));
	
	%conf with delta-rate
	errorbar(h5,vecMean_dR,vecMeanConf_dR,vecSemConf_dR,vecSemConf_dR,vecSem_dR,vecSem_dR,'color',lines(1));
	
	
	%% save data
	cellQuantileDur(intFile,:) = cellValsDur;
	cellQuantileConf(intFile,:) = cellValsConf;
	matLatConf(intFile,:) = vecMeans;
	vecSpikeGroupSize(intFile) = intSpikeGroupSize;
	cellQuantileRateChange(intFile,:) = cellVals_dR;
	cellQuantileConf_dR(intFile,:) = cellValsConf_dR;
	matRateConf(intFile,:) = vecMeansR;
end
plot(h1,vecBinsDurC,mean(matLatConf,1),'color',lines(1));%,vecSDs./sqrt(vecCounts),vecSDs./sqrt(vecCounts),'color',lines(1));
hold(h1,'off')
xlabel(h1,'Latency of n-spike block after stim onset (s)');
ylabel(h1,'Decoder confidence');
ylim(h1,[0 max(get(h1,'ylim'))]);

hold(h2,'off');
xlabel(h2,'Duration of n-spike block (ms)');
ylabel(h2,'Decoder confidence');
title(h2,strRunType);
ylim(h2,[0 max(get(h2,'ylim'))]);

plot(h4,vecBinsRateC,nanmean(matRateConf,1),'color',lines(1));%,vecSDs./sqrt(vecCounts),vecSDs./sqrt(vecCounts),'color',lines(1));
hold(h4,'off')
xlabel(h4,'Spiking rate during block (Hz)');
ylabel(h4,'Decoder confidence');
ylim(h4,[0 max(get(h4,'ylim'))]);

hold(h5,'off');
xlabel(h5,'Change in rate during block (Hz)');
ylabel(h5,'Decoder confidence');
title(h5,strRunType);
ylim(h5,[0 max(get(h5,'ylim'))]);

%%
matQuantDur = cellfun(@mean,cellQuantileDur);
matQuantConf = cellfun(@mean,cellQuantileConf);
matQuantConf = zscore(matQuantConf,[],2);
matX = repmat(1:10,[intRecNum 1]);
mdl = fitlm(matX(:),matQuantConf(:));
r=mdl.Coefficients.Estimate(2);
p=mdl.Coefficients.pValue(2);

hold(h3,'on')
plot(h3,matX',matQuantConf','color',vecColGrey);
plot(h3,mean(matX,1),mean(matQuantConf,1),'color',lines(1));
hold(h3,'off')
xlabel(h3,'Duration decile of n-spike block');
ylabel(h3,sprintf('Decoder confidence, z-scored per rec (%s)',getGreek('sigma')));
%ylim(h3,[0 max(get(h3,'ylim'))]);
title(h3,sprintf('Lin reg, r=%.3f %s/decile, p=%.3e',r,getGreek('sigma'),p));

matQuantRateChange = cellfun(@mean,cellQuantileRateChange);
matQuantConf_dR = cellfun(@mean,cellQuantileConf_dR);
matQuantConf_dR = zscore(matQuantConf_dR,[],2);
matX = repmat(1:10,[intRecNum 1]);
mdl = fitlm(matX(:),matQuantConf_dR(:));
r_dR=mdl.Coefficients.Estimate(2);
p_dR=mdl.Coefficients.pValue(2);

hold(h6,'on')
plot(h6,matX',matQuantConf_dR','color',vecColGrey);
plot(h6,mean(matX,1),mean(matQuantConf_dR,1),'color',lines(1));
hold(h6,'off')
xlabel(h6,'Rate-change decile of n-spike block');
ylabel(h6,sprintf('Decoder confidence, z-scored per rec (%s)',getGreek('sigma')));
%ylim(h6,[0 max(get(h6,'ylim'))]);
title(h6,sprintf('Lin reg, r=%.3f %s/decile, p=%.3e',r_dR,getGreek('sigma'),p_dR));

%%
fixfig;

%%
drawnow;
strType = strrep(strRunType,',','_');
export_fig(fullpath(strFigurePath,sprintf('Q1d_SpikeBlockDecoding%s.tif',strType)));
export_fig(fullpath(strFigurePath,sprintf('Q1d_SpikeBlockDecoding%s.pdf',strType)));
end
end
end
end