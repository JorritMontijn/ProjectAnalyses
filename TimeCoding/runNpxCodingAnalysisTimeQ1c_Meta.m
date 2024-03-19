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
strType = 'Real'; %'Real' 'ShuffTid'
boolFixedSpikeGroupSize = true;
intQuantiles = 10;
vecBinsDur = 0:0.01:1;
vecBinsDurC = vecBinsDur(2:end) - diff(vecBinsDur(1:2))/2;
vecColGrey = [0.7 0.7 0.7];

if boolFixedSpikeGroupSize
	sFiles = dir ([strTargetDataPath 'Q1cData*' strType '*Fixed*.mat']);
	strSGS = 'FixedSGS';
else
	sFiles = dir ([strTargetDataPath 'Q1cData*' strType '*Var*.mat']);
	strSGS = 'VarSGS';
end

intRecNum = numel(sFiles);
cellQuantileDur = cell(intRecNum,intQuantiles);
cellQuantileConf = cell(intRecNum,intQuantiles);
matLatConf = nan(intRecNum,numel(vecBinsDurC));
vecSpikeGroupSize = nan(intRecNum,1);

vecR_Num_Conf = nan(intRecNum,1);
vecR_Rate_Conf = nan(intRecNum,1);
vecR_Tune_Conf = nan(intRecNum,1);
vecR_IFR_Conf = nan(intRecNum,1);

vecR_CombConf_CellNum = nan(intRecNum,1);
vecR_CombConf_Rate = nan(intRecNum,1);
vecR_CombConf_Tune = nan(intRecNum,1);
vecR_CombConf_IFR = nan(intRecNum,1);


vecR_Num_IFR = nan(intRecNum,1);
vecR_Rate_IFR = nan(intRecNum,1);
vecR_Tune_IFR = nan(intRecNum,1);

vecR_CombIFR_Rate = nan(intRecNum,1);
vecR_CombIFR_Tune = nan(intRecNum,1);

vecR_Group_Rate_Tune = nan(intRecNum,1);

%% plot
figure;maxfig;
h1=subplot(2,3,2);
h2=subplot(2,3,1);
h3=subplot(2,3,3);
vecColChance = [0.3 0.3 0.3];
plot(h1,[0 100],([1 1]./12),'--','color',vecColChance);
hold(h1,'on')

plot(h2,[0 1],([1 1]./12),'--','color',vecColChance);
hold(h2,'on')
title(h2,sprintf('%s, %s',strType,strSGS));

%plot(h3,[0 10],([1 1]./12),'--','color',vecColChance);
hold(h3,'on')

for intFile=1:intRecNum
	%% load
	sLoad=load(fullpath(sFiles(intFile).folder,sFiles(intFile).name));
	sSpikeGroup = sLoad.sSpikeGroup;
	strRec = sLoad.strRec;
	strType = sLoad.strType;
	vecOri180 = sLoad.vecOri180;
	vecOrientation = sLoad.vecOrientation;
	vecStimOffTime = sLoad.vecStimOffTime;
	vecStimOnTime = sLoad.vecStimOnTime;
	intSpikeGroupSize = sLoad.intSpikeGroupSize;
	rSpike_IFR_Rate = sLoad.rSpike_IFR_Rate;
	pSpike_IFR_Rate = sLoad.pSpike_IFR_Rate;
	rSpike_IFR_Tune = sLoad.rSpike_IFR_Tune;
	pSpike_IFR_Tune = sLoad.pSpike_IFR_Tune;
	vecTuningPerCell = sLoad.vecTuningPerCell;
	vecRatePerCell = sLoad.vecRatePerCell;
	
	%% create derived variables
	%put confidence in deciles per recording
	vecSpikeGroupDuration = [sSpikeGroup.Duration]';
	vecSpikeGroupCorrect = [sSpikeGroup.Correct]';
	vecSpikeGroupConfidence = [sSpikeGroup.Confidence]';
	vecSpikeGroupLatency = [sSpikeGroup.Latency]';
	vecSpikeGroupNumOfCells = [sSpikeGroup.NumOfCells]';
	vecSpikeGroupAvgRateOfCells = [sSpikeGroup.AvgRateOfCells]';
	vecSpikeGroupAvgTuningOfCells = [sSpikeGroup.AvgTuningOfCells]';
	vecSpikeGroupAvgIFR = [sSpikeGroup.AvgRate]';
	
	%sort
	[vecSortedDur,vecSort]=sort(vecSpikeGroupDuration);
	vecSortedCorr = vecSpikeGroupCorrect(vecSort);
	vecSortedConf = vecSpikeGroupConfidence(vecSort);
	vecSortedLat = vecSpikeGroupLatency(vecSort);
	%remove long durs
	%indRem = vecSortedDur > 0.05;
	%indRem = vecSortedLat < 0.15 | vecSortedDur > 0.05;
	indRem = [];
	vecSortedDur(indRem) = [];
	vecSortedCorr(indRem) = [];
	vecSortedConf(indRem) = [];
	vecSortedLat(indRem) = [];
	
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
	vecSampleGroup = zeros(size(vecSortedDur));
	for intBin=1:intBins
		intEndS = intSperBin*intBin;
		vecSamples = (intEndS-intSperBin+1):intEndS;
		
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
	end
	
	%conf with dur
	errorbar(h1,vecMeanDur*1000,vecMeanConf,vecSemConf,vecSemConf,vecSemDur,vecSemDur,'color',lines(1));
	
	%conf with time
	[vecCounts,vecMeans,vecSDs]=makeBins(vecSpikeGroupLatency,vecSpikeGroupConfidence,vecBinsDur);
	plot(h2,vecBinsDurC,vecMeans,'color',vecColGrey);%,vecSDs./sqrt(vecCounts),vecSDs./sqrt(vecCounts),'color',lines(1));
	
	%%
	%predict conf
	[r_Num_Conf,p_Num_Conf]=corr(vecSpikeGroupNumOfCells,vecSpikeGroupConfidence);
	[r_Rate_Conf,p_Rate_Conf]=corr(vecSpikeGroupAvgRateOfCells,vecSpikeGroupConfidence);
	[r_Tune_Conf,p_Tune_Conf]=corr(vecSpikeGroupAvgTuningOfCells,vecSpikeGroupConfidence);
	[r_IFR_Conf,p_IFR_Conf]=corr(vecSpikeGroupAvgIFR,vecSpikeGroupConfidence);
	
	tbl = table(vecSpikeGroupNumOfCells,vecSpikeGroupAvgRateOfCells,vecSpikeGroupAvgTuningOfCells,vecSpikeGroupAvgIFR,vecSpikeGroupConfidence,...
		'VariableNames',{'CellNumber','AvgCellRate','AvgCellTuning','AvgPopIFR','Conf'});
	mdl_predconf = fitlm(tbl,'linear');
	
	%predict ifr
	[r_Num_IFR,p_Num_IFR]=corr(vecSpikeGroupNumOfCells,vecSpikeGroupAvgIFR);
	[r_Rate_IFR,p_Rate_IFR]=corr(vecSpikeGroupAvgRateOfCells,vecSpikeGroupAvgIFR);
	[r_Tune_IFR,p_Tune_IFR]=corr(vecSpikeGroupAvgTuningOfCells,vecSpikeGroupAvgIFR);
	tbl = table(vecSpikeGroupAvgRateOfCells,vecSpikeGroupAvgTuningOfCells,vecSpikeGroupAvgIFR,...
		'VariableNames',{'AvgCellRate','AvgCellTuning','AvgPopIFR'});
	mdl_predifr = fitlm(tbl,'linear');
	
	%rate vs tuning
	[rGroup_Rate_Tune,pGroup_Rate_Tune]=corr(vecTuningPerCell,vecRatePerCell);
	
	%% save
	%result: during high ifr epochs, cells with low firing rate and low tuning are also active,
	%while during low ifr epochs, only cells with high firing rate or high tuning are active
	vecR_Num_Conf(intFile) = r_Num_Conf;
	vecR_Rate_Conf(intFile) = r_Rate_Conf;
	vecR_Tune_Conf(intFile) = r_Tune_Conf;
	vecR_IFR_Conf(intFile) = r_IFR_Conf;
	
	vecR_CombConf_CellNum(intFile) = table2array(mdl_predconf.Coefficients(2,1));
	vecR_CombConf_Rate(intFile) = table2array(mdl_predconf.Coefficients(3,1));
	vecR_CombConf_Tune(intFile) = table2array(mdl_predconf.Coefficients(4,1));
	vecR_CombConf_IFR(intFile) = table2array(mdl_predconf.Coefficients(5,1));
	
	
	vecR_Num_IFR(intFile) = r_Num_IFR;
	vecR_Rate_IFR(intFile) = r_Rate_IFR;
	vecR_Tune_IFR(intFile) = r_Tune_IFR;
	
	vecR_CombIFR_Rate(intFile) = table2array(mdl_predifr.Coefficients(2,1));
	vecR_CombIFR_Tune(intFile) = table2array(mdl_predifr.Coefficients(3,1));
	
	vecR_Group_Rate_Tune(intFile) = rGroup_Rate_Tune;
	
	cellQuantileDur(intFile,:) = cellValsDur;
	cellQuantileConf(intFile,:) = cellValsConf;
	matLatConf(intFile,:) = vecMeans;
	vecSpikeGroupSize(intFile) = intSpikeGroupSize;
end
hold(h1,'off');
xlabel(h1,'Duration of n-spike block (ms)');
ylabel(h1,'Decoder confidence');
title(h1,sprintf('Deciles'));
ylim(h1,[0 max(get(h1,'ylim'))]);

plot(h2,vecBinsDurC,mean(matLatConf,1),'color',lines(1));%,vecSDs./sqrt(vecCounts),vecSDs./sqrt(vecCounts),'color',lines(1));
hold(h2,'off')
xlabel(h2,'Latency of n-spike block after stim onset (s)');
ylabel(h2,'Decoder confidence');
ylim(h2,[0 max(get(h2,'ylim'))]);


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
fixfig;

%% plot corrs
subplot(2,3,4)
hold on
vecX = ones(size(vecR_Num_Conf));
swarmchart(1*vecX,vecR_Num_Conf);
swarmchart(2*vecX,vecR_Rate_Conf);
swarmchart(3*vecX,vecR_Tune_Conf);
swarmchart(4*vecX,vecR_IFR_Conf);
 





vecR_Num_IFR
vecR_Rate_IFR
vecR_Tune_IFR


%%
return
drawnow;
export_fig(fullpath(strFigurePath,sprintf('Q1c_SpikeBlockDecoding%s%s.tif',strType,strSGS)));
export_fig(fullpath(strFigurePath,sprintf('Q1c_SpikeBlockDecoding%s%s.pdf',strType,strSGS)));
