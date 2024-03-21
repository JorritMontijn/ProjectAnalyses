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
boolSaveFig = true;
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
strType = 'ShuffTid'; %'Real' 'ShuffTid'
cellTypes = {'Real','ShuffTid'};
boolFixedSpikeGroupSize = false;
intQuantiles = 10;
intGroupOrDecile = 1;
dblRemOnset = 0.25; %remove onset period in seconds

vecBinsDur = 0:0.01:1;
vecBinsDurC = vecBinsDur(2:end) - diff(vecBinsDur(1:2))/2;
vecColGrey = [0.7 0.7 0.7];

for dblRemOnset=[0 0.25]
for intType=1:2
strType = cellTypes{intType};
if boolFixedSpikeGroupSize
	sFiles = dir ([strTargetDataPath 'Q1cData*' strType '*Fixed*.mat']);
	strSGS = 'FixedSGS';
else
	sFiles = dir ([strTargetDataPath 'Q1cData*' strType '*Var*.mat']);
	strSGS = 'VarSGS';
end
if dblRemOnset == 0
	strOnset = '';
	indWithOnset = cellfun(@(x) strcmp(x((end-7):end),['0.25.mat']), {sFiles.name}');
	sFiles(indWithOnset) = [];
else
	strOnset = sprintf('%.2f',dblRemOnset);
	indWithOnset = cellfun(@(x) strcmp(x((end-7):end),['0.25.mat']), {sFiles.name}');
	sFiles(~indWithOnset) = [];
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
h1=subplot(2,3,1);
h2=subplot(2,3,2);
h3=subplot(2,3,3);
vecColChance = [0.3 0.3 0.3];
plot(h1,[0 1],([1 1]./12),'--','color',vecColChance);
hold(h1,'on')
title(h1,sprintf('%s, %s',strType,strSGS));

plot(h2,[0 100],([1 1]./12),'--','color',vecColChance);
hold(h2,'on')

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
	%rSpike_IFR_Rate = sLoad.rSpike_IFR_Rate;
	%pSpike_IFR_Rate = sLoad.pSpike_IFR_Rate;
	%rSpike_IFR_Tune = sLoad.rSpike_IFR_Tune;
	%pSpike_IFR_Tune = sLoad.pSpike_IFR_Tune;
	vecTuningPerCell = sLoad.vecTuningPerCell;
	vecRatePerCell = sLoad.vecRatePerCell;
	
	%% create derived variables
	%put confidence in deciles per recording
	vecSpikeGroupDuration = [sSpikeGroup.Duration]';
	vecSpikeGroupCorrect = [sSpikeGroup.Correct]';
	vecSpikeGroupConfidence = [sSpikeGroup.Confidence]';
	vecSpikeGroupLatency = dblRemOnset+[sSpikeGroup.Latency]';
	vecSpikeGroupNumOfCells = [sSpikeGroup.NumOfCells]';
	vecSpikeGroupAvgRateOfCells = [sSpikeGroup.AvgRateOfCells]';
	vecSpikeGroupAvgTuningOfCells = [sSpikeGroup.AvgTuningOfCells]';
	vecSpikeGroupAvgIFR = [sSpikeGroup.AvgRate]';
	
	%sort by duration
	[vecSortedDur,vecSort]=sort(vecSpikeGroupDuration);
	vecSortedCorr = vecSpikeGroupCorrect(vecSort);
	vecSortedConf = vecSpikeGroupConfidence(vecSort);
	vecSortedLat = vecSpikeGroupLatency(vecSort);
	vecSortedNum = vecSpikeGroupNumOfCells(vecSort);
	vecSortedRate = vecSpikeGroupAvgRateOfCells(vecSort);
	vecSortedTune = vecSpikeGroupAvgTuningOfCells(vecSort);
	vecSortedIFR = vecSpikeGroupAvgIFR(vecSort);
	
	%calculate fraction correct and confidence per bin of equal size
	intQuantileNum = 10;
	[vecMeanDur,vecSemDur,vecMeanConf,vecSemConf,vecQuantile,cellValsDur,cellValsConf]=getQuantiles(vecSortedDur,vecSortedConf,intQuantileNum);
	
	%conf with dur
	errorbar(h2,vecMeanDur*1000,vecMeanConf,vecSemConf,vecSemConf,vecSemDur,vecSemDur,'color',lines(1));
	
	%conf with time
	[vecCounts,vecMeans,vecSDs]=makeBins(vecSpikeGroupLatency,vecSpikeGroupConfidence,vecBinsDur);
	plot(h1,vecBinsDurC,vecMeans,'color',vecColGrey);%,vecSDs./sqrt(vecCounts),vecSDs./sqrt(vecCounts),'color',lines(1));
	
	%% group or decile
	if intGroupOrDecile == 2
		%% decile based
		% deciles of dur
		[vecMeanDurDeciles,vecSemDurDeciles,vecMeanNum,vecSemNum]=getQuantiles(vecSpikeGroupDuration,vecSpikeGroupNumOfCells,intQuantileNum);
		[vecMeanDurDeciles,vecSemDurDeciles,vecMeanRate,vecSemRate]=getQuantiles(vecSpikeGroupDuration,vecSpikeGroupAvgRateOfCells,intQuantileNum);
		[vecMeanDurDeciles,vecSemDurDeciles,vecMeanTune,vecSemTune]=getQuantiles(vecSpikeGroupDuration,vecSpikeGroupAvgTuningOfCells,intQuantileNum);
		[vecMeanDurDeciles,vecSemDurDeciles,vecMeanIFR,vecSemIFR]=getQuantiles(vecSpikeGroupDuration,vecSpikeGroupAvgIFR,intQuantileNum);
		[vecMeanDurDeciles,vecSemDurDeciles,vecMeanConf,vecSemConf]=getQuantiles(vecSpikeGroupDuration,vecSpikeGroupConfidence,intQuantileNum);
		vecMeanDurDeciles = vecMeanDurDeciles';
		vecMeanNum = vecMeanNum';
		vecMeanRate = vecMeanRate';
		vecMeanTune = vecMeanTune';
		vecMeanIFR = vecMeanIFR';
		vecMeanConf = vecMeanConf';
		
		%predict dur
		[r_Num_Dur,p_Num_Dur]=corr(vecMeanNum,vecMeanDurDeciles);
		[r_Rate_Dur,p_Rate_Dur]=corr(vecMeanRate,vecMeanDurDeciles);
		[r_Tune_Dur,p_Tune_Dur]=corr(vecMeanTune,vecMeanDurDeciles);
		[r_IFR_Dur,p_IFR_Dur]=corr(vecMeanIFR,vecMeanDurDeciles);
		[r_Conf_Dur,p_Conf_Dur]=corr(vecMeanConf,vecMeanDurDeciles);
		
		tbl = table(vecMeanNum,vecMeanRate,vecMeanTune,vecMeanIFR,vecMeanConf,vecMeanDurDeciles,...
			'VariableNames',{'CellNumber','AvgCellRate','AvgCellTuning','AvgPopIFR','Conf','Dur'});
		mdl_preddur = fitlm(tbl,'linear');
		
		% deciles of conf
		[vecMeanConfDeciles,vecSemConfDeciles,vecMeanNum,vecSemNum]=getQuantiles(vecSpikeGroupConfidence,vecSpikeGroupNumOfCells,intQuantileNum);
		[vecMeanConfDeciles,vecSemConfDeciles,vecMeanRate,vecSemRate]=getQuantiles(vecSpikeGroupConfidence,vecSpikeGroupAvgRateOfCells,intQuantileNum);
		[vecMeanConfDeciles,vecSemConfDeciles,vecMeanTune,vecSemTune]=getQuantiles(vecSpikeGroupConfidence,vecSpikeGroupAvgTuningOfCells,intQuantileNum);
		[vecMeanConfDeciles,vecSemConfDeciles,vecMeanIFR,vecSemIFR]=getQuantiles(vecSpikeGroupConfidence,vecSpikeGroupAvgIFR,intQuantileNum);
		vecMeanConfDeciles = vecMeanConfDeciles';
		vecMeanNum = vecMeanNum';
		vecMeanRate = vecMeanRate';
		vecMeanTune = vecMeanTune';
		vecMeanIFR = vecMeanIFR';
		
		%predict conf
		[r_Num_Conf,p_Num_Conf]=corr(vecMeanNum,vecMeanConfDeciles);
		[r_Rate_Conf,p_Rate_Conf]=corr(vecMeanRate,vecMeanConfDeciles);
		[r_Tune_Conf,p_Tune_Conf]=corr(vecMeanTune,vecMeanConfDeciles);
		[r_IFR_Conf,p_IFR_Conf]=corr(vecMeanIFR,vecMeanConfDeciles);
		
		tbl = table(vecMeanNum,vecMeanRate,vecMeanTune,vecMeanIFR,vecMeanConfDeciles,...
			'VariableNames',{'CellNumber','AvgCellRate','AvgCellTuning','AvgPopIFR','Conf'});
		mdl_predconf = fitlm(tbl,'linear');
		
		% deciles of ifr
		[vecMeanIFRDeciles,vecSemIFRDeciles,vecMeanNum,vecSemNum,vecConfQ]=getQuantiles(vecSpikeGroupAvgIFR,vecSpikeGroupNumOfCells,intQuantileNum);
		[vecMeanIFRDeciles,vecSemIFRDeciles,vecMeanRate,vecSemRate]=getQuantiles(vecSpikeGroupAvgIFR,vecSpikeGroupAvgRateOfCells,intQuantileNum);
		[vecMeanIFRDeciles,vecSemIFRDeciles,vecMeanTune,vecSemTune]=getQuantiles(vecSpikeGroupAvgIFR,vecSpikeGroupAvgTuningOfCells,intQuantileNum);
		vecMeanIFRDeciles = vecMeanIFRDeciles';
		vecMeanNum = vecMeanNum';
		vecMeanRate = vecMeanRate';
		vecMeanTune = vecMeanTune';
		
		%predict ifr
		[r_Num_IFR,p_Num_IFR]=corr(vecMeanNum,vecMeanIFRDeciles);
		[r_Rate_IFR,p_Rate_IFR]=corr(vecMeanRate,vecMeanIFRDeciles);
		[r_Tune_IFR,p_Tune_IFR]=corr(vecMeanTune,vecMeanIFRDeciles);
		tbl = table(vecMeanRate,vecMeanTune,vecMeanIFRDeciles,...
			'VariableNames',{'AvgCellRate','AvgCellTuning','AvgPopIFR'});
		mdl_predifr = fitlm(tbl,'linear');
		
		%rate vs tuning
		[rGroup_Rate_Tune,pGroup_Rate_Tune]=corr(vecTuningPerCell,vecRatePerCell);
		
	else
		%% single-group based
		%predict dur
		[r_Num_Dur,p_Num_Dur]=corr(vecSpikeGroupNumOfCells,vecSpikeGroupDuration);
		[r_Rate_Dur,p_Rate_Dur]=corr(vecSpikeGroupAvgRateOfCells,vecSpikeGroupDuration);
		[r_Tune_Dur,p_Tune_Dur]=corr(vecSpikeGroupAvgTuningOfCells,vecSpikeGroupDuration);
		[r_IFR_Dur,p_IFR_Dur]=corr(vecSpikeGroupAvgIFR,vecSpikeGroupDuration);
		[r_Conf_Dur,p_Conf_Dur]=corr(vecSpikeGroupConfidence,vecSpikeGroupDuration);
		
		tbl = table(vecSpikeGroupNumOfCells,vecSpikeGroupAvgRateOfCells,vecSpikeGroupAvgTuningOfCells,vecSpikeGroupAvgIFR,vecSpikeGroupConfidence,vecSpikeGroupDuration,...
			'VariableNames',{'CellNumber','AvgCellRate','AvgCellTuning','AvgPopIFR','Conf','Dur'});
		mdl_preddur = fitlm(tbl,'linear');
		
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
		%}
	end
	%% save
	%result: during high ifr epochs, cells with low firing rate and low tuning are also active,
	%while during low ifr epochs, only cells with high firing rate or high tuning are active
	
	%dur
	vecR_Num_Dur(intFile) = r_Num_Dur;
	vecR_Rate_Dur(intFile) = r_Rate_Dur;
	vecR_Tune_Dur(intFile) = r_Tune_Dur;
	vecR_IFR_Dur(intFile) = r_IFR_Dur;
	vecR_Conf_Dur(intFile) = r_Conf_Dur;
	
	vecR_CombDur_CellNum(intFile) = table2array(mdl_preddur.Coefficients(2,1));
	vecR_CombDur_Rate(intFile) = table2array(mdl_preddur.Coefficients(3,1));
	vecR_CombDur_Tune(intFile) = table2array(mdl_preddur.Coefficients(4,1));
	vecR_CombDur_IFR(intFile) = table2array(mdl_preddur.Coefficients(5,1));
	vecR_CombDur_Conf(intFile) = table2array(mdl_preddur.Coefficients(6,1));
	
	%conf
	vecR_Num_Conf(intFile) = r_Num_Conf;
	vecR_Rate_Conf(intFile) = r_Rate_Conf;
	vecR_Tune_Conf(intFile) = r_Tune_Conf;
	vecR_IFR_Conf(intFile) = r_IFR_Conf;
	
	vecR_CombConf_CellNum(intFile) = table2array(mdl_predconf.Coefficients(2,1));
	vecR_CombConf_Rate(intFile) = table2array(mdl_predconf.Coefficients(3,1));
	vecR_CombConf_Tune(intFile) = table2array(mdl_predconf.Coefficients(4,1));
	vecR_CombConf_IFR(intFile) = table2array(mdl_predconf.Coefficients(5,1));
	
	%ifr
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
	
	%% add quantile data
	vecSpikeGroupDuration = [sSpikeGroup.Duration]';
	vecSpikeGroupCorrect = [sSpikeGroup.Correct]';
	vecSpikeGroupConfidence = [sSpikeGroup.Confidence]';
	vecSpikeGroupLatency = [sSpikeGroup.Latency]';
	vecSpikeGroupNumOfCells = [sSpikeGroup.NumOfCells]';
	vecSpikeGroupAvgRateOfCells = [sSpikeGroup.AvgRateOfCells]';
	vecSpikeGroupAvgTuningOfCells = [sSpikeGroup.AvgTuningOfCells]';
	vecSpikeGroupAvgIFR = [sSpikeGroup.AvgRate]';
	
	%by dur
	intQuantileNum = 10;
	[~,~,~,~,~,cellValsDur_Dur,cellValsConf_Dur]=getQuantiles(vecSpikeGroupDuration,vecSpikeGroupConfidence,intQuantileNum);
	[~,~,~,~,~,~,cellValsNum_Dur]=getQuantiles(vecSpikeGroupDuration,vecSpikeGroupNumOfCells,intQuantileNum);
	[~,~,~,~,~,~,cellValsRate_Dur]=getQuantiles(vecSpikeGroupDuration,vecSpikeGroupAvgRateOfCells,intQuantileNum);
	[~,~,~,~,~,~,cellValsTune_Dur]=getQuantiles(vecSpikeGroupDuration,vecSpikeGroupAvgTuningOfCells,intQuantileNum);
	[~,~,~,~,~,~,cellValsIFR_Dur]=getQuantiles(vecSpikeGroupDuration,vecSpikeGroupAvgIFR,intQuantileNum);
	cellQuantileDur_Dur(intFile,:) = cellValsDur_Dur;
	cellQuantileConf_Dur(intFile,:) = cellValsConf_Dur;
	cellQuantileNum_Dur(intFile,:) = cellValsNum_Dur;
	cellQuantileRate_Dur(intFile,:) = cellValsRate_Dur;
	cellQuantileTune_Dur(intFile,:) = cellValsTune_Dur;
	cellQuantileIFR_Dur(intFile,:) = cellValsIFR_Dur;

	%by conf
	[~,~,~,~,~,cellValsConf_Conf,cellValsNum_Conf]=getQuantiles(vecSpikeGroupConfidence,vecSpikeGroupNumOfCells,intQuantileNum);
	[~,~,~,~,~,~,cellValsRate_Conf]=getQuantiles(vecSpikeGroupConfidence,vecSpikeGroupAvgRateOfCells,intQuantileNum);
	[~,~,~,~,~,~,cellValsTune_Conf]=getQuantiles(vecSpikeGroupConfidence,vecSpikeGroupAvgTuningOfCells,intQuantileNum);
	[~,~,~,~,~,~,cellValsIFR_Conf]=getQuantiles(vecSpikeGroupConfidence,vecSpikeGroupAvgIFR,intQuantileNum);
	cellQuantileConf_Conf(intFile,:) = cellValsConf_Conf;
	cellQuantileNum_Conf(intFile,:) = cellValsNum_Conf;
	cellQuantileRate_Conf(intFile,:) = cellValsRate_Conf;
	cellQuantileTune_Conf(intFile,:) = cellValsTune_Conf;
	cellQuantileIFR_Conf(intFile,:) = cellValsIFR_Conf;
	
	%by ifr
	[~,~,~,~,~,cellValsIFR_IFR,cellValsConf_IFR]=getQuantiles(vecSpikeGroupAvgIFR,vecSpikeGroupConfidence,intQuantileNum);
	[~,~,~,~,~,~,cellValsNum_IFR]=getQuantiles(vecSpikeGroupAvgIFR,vecSpikeGroupNumOfCells,intQuantileNum);
	[~,~,~,~,~,~,cellValsRate_IFR]=getQuantiles(vecSpikeGroupAvgIFR,vecSpikeGroupAvgRateOfCells,intQuantileNum);
	[~,~,~,~,~,~,cellValsTune_IFR]=getQuantiles(vecSpikeGroupAvgIFR,vecSpikeGroupAvgTuningOfCells,intQuantileNum);
	cellQuantileConf_IFR(intFile,:) = cellValsConf_IFR;
	cellQuantileNum_IFR(intFile,:) = cellValsNum_IFR;
	cellQuantileRate_IFR(intFile,:) = cellValsRate_IFR;
	cellQuantileTune_IFR(intFile,:) = cellValsTune_IFR;
	cellQuantileIFR_IFR(intFile,:) = cellValsIFR_IFR;
end
hold(h2,'off');
xlabel(h2,'Duration of n-spike block (ms)');
ylabel(h2,'Decoder confidence');
title(h2,sprintf('Deciles'));
ylim(h2,[0 max(get(h2,'ylim'))]);

plot(h1,vecBinsDurC,mean(matLatConf,1),'color',lines(1));%,vecSDs./sqrt(vecCounts),vecSDs./sqrt(vecCounts),'color',lines(1));
hold(h1,'off')
xlabel(h1,'Latency of n-spike block after stim onset (s)');
ylabel(h1,'Decoder confidence');
ylim(h1,[0 max(get(h1,'ylim'))]);


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
%test vs 0
[h,p_Num_Dur] = ttest(vecR_Num_Dur);
[h,p_Rate_Dur] = ttest(vecR_Rate_Dur);
[h,p_Tune_Dur] = ttest(vecR_Tune_Dur);
[h,p_IFR_Dur] = ttest(vecR_IFR_Dur);
[h,p_Conf_Dur] = ttest(vecR_Conf_Dur);
[h, crit_p, adj_p_Dur]=fdr_bh([p_Num_Dur p_Rate_Dur p_Tune_Dur p_IFR_Dur p_Conf_Dur],0.05);

%plot
subplot(2,3,4);cla
hold on;
plot([0 10],[0 0],'--','color',[0.5 0.5 0.5]);
matCol = lines(5);
matCol = matCol([2 1 3 4 5],:);
matCol(3,:) = 0;
matColP = 1-((1-matCol)./2);

intRecNum = numel(vecR_Num_Dur);
vecX = ones(size(vecR_Num_Dur));
h2=swarmchart(1*vecX,vecR_Num_Dur,[],matColP(1,:),'filled');
swarmchart(3*vecX,vecR_Rate_Dur,[],matColP(2,:),'filled');
swarmchart(5*vecX,vecR_Tune_Dur,[],matColP(3,:),'filled');
swarmchart(7*vecX,vecR_IFR_Dur,[],matColP(4,:),'filled');
swarmchart(9*vecX,vecR_Conf_Dur,[],matColP(5,:),'filled');
errorbar(1,mean(vecR_Num_Dur),std(vecR_Num_Dur)./sqrt(intRecNum),'x','color',matCol(1,:));
errorbar(3,mean(vecR_Rate_Dur),std(vecR_Rate_Dur)./sqrt(intRecNum),'x','color',matCol(2,:));
errorbar(5,mean(vecR_Tune_Dur),std(vecR_Tune_Dur)./sqrt(intRecNum),'x','color',matCol(3,:));
errorbar(7,mean(vecR_IFR_Dur),std(vecR_IFR_Dur)./sqrt(intRecNum),'x','color',matCol(4,:));
errorbar(9,mean(vecR_Conf_Dur),std(vecR_Conf_Dur)./sqrt(intRecNum),'x','color',matCol(5,:));
hold off
set(gca,'xtick',[1:2:9],'xticklabel',{'# of cells','Cell Rate','Cell Tune','Pop Rate','Conf'});
xtickangle(gca,45);
ylabel('Correlation with group dur (r)');
title(sprintf('p: #=%.2e,R=%.3f,T=%.2e,IFR=%.2e,C=%.2e',adj_p_Dur))
fixfig;

%% conf
%test vs 0
[h,p_Num_Conf] = ttest(vecR_Num_Conf);
[h,p_Rate_Conf] = ttest(vecR_Rate_Conf);
[h,p_Tune_Conf] = ttest(vecR_Tune_Conf);
[h,p_IFR_Conf] = ttest(vecR_IFR_Conf);
[h, crit_p, adj_p_Conf]=fdr_bh([p_Num_Conf p_Rate_Conf p_Tune_Conf p_IFR_Conf],0.05);

%plot
subplot(2,3,5);cla
hold on;
plot([0 8],[0 0],'--','color',[0.5 0.5 0.5]);
matCol = lines(5);
matCol = matCol([2 1 3:end],:);
matCol(3,:) = 0;
matColP = 1-((1-matCol)./2);
intRecNum = numel(vecR_Num_Conf);
vecX = ones(size(vecR_Num_Conf));
h2=swarmchart(1*vecX,vecR_Num_Conf,[],matColP(1,:),'filled');
swarmchart(3*vecX,vecR_Rate_Conf,[],matColP(2,:),'filled');
swarmchart(5*vecX,vecR_Tune_Conf,[],matColP(3,:),'filled');
swarmchart(7*vecX,vecR_IFR_Conf,[],matColP(4,:),'filled');
errorbar(1,mean(vecR_Num_Conf),std(vecR_Num_Conf)./sqrt(intRecNum),'x','color',matCol(1,:));
errorbar(3,mean(vecR_Rate_Conf),std(vecR_Rate_Conf)./sqrt(intRecNum),'x','color',matCol(2,:));
errorbar(5,mean(vecR_Tune_Conf),std(vecR_Tune_Conf)./sqrt(intRecNum),'x','color',matCol(3,:));
errorbar(7,mean(vecR_IFR_Conf),std(vecR_IFR_Conf)./sqrt(intRecNum),'x','color',matCol(4,:));
hold off
set(gca,'xtick',[1:2:7],'xticklabel',{'# of cells','Cell Rate','Cell Tune','Pop Rate'});
xtickangle(gca,45);
ylabel('Correlation with confidence (r)');
title(sprintf('p: #=%.3f,R=%.3e,T=%.3e,IFR=%.3e',adj_p_Conf))
fixfig;

%% ifr
%test vs 0
[h,p_Num_IFR] = ttest(vecR_Num_IFR);
[h,p_Rate_IFR] = ttest(vecR_Rate_IFR);
[h,p_Tune_IFR] = ttest(vecR_Tune_IFR);
[h, crit_p, adj_p_IFR]=fdr_bh([p_Num_IFR p_Rate_IFR p_Tune_IFR],0.05);

%plot
subplot(2,3,6);cla
hold on;
plot([0 6],[0 0],'--','color',[0.5 0.5 0.5]);

h2=swarmchart(1*vecX,vecR_Num_IFR,[],matColP(1,:),'filled');
swarmchart(3*vecX,vecR_Rate_IFR,[],matColP(2,:),'filled');
swarmchart(5*vecX,vecR_Tune_IFR,[],matColP(3,:),'filled');

errorbar(1,mean(vecR_Num_IFR),std(vecR_Num_IFR)./sqrt(intRecNum),'x','color',matCol(1,:));
errorbar(3,mean(vecR_Rate_IFR),std(vecR_Rate_IFR)./sqrt(intRecNum),'x','color',matCol(2,:));
errorbar(5,mean(vecR_Tune_IFR),std(vecR_Tune_IFR)./sqrt(intRecNum),'x','color',matCol(3,:));

hold off
set(gca,'xtick',[1:2:5],'xticklabel',{'# of cells','Cell Rate','Cell Tune'});
xtickangle(gca,45);
ylabel('Correlation with pop rate (r)');
title(sprintf('p: #=%.3e,R=%.3f,T=%.3e',adj_p_IFR))
fixfig;
%%
if boolSaveFig
drawnow;
export_fig(fullpath(strFigurePath,sprintf('Q1c_SpikeBlockDecoding%s%s%s.tif',strType,strSGS,strOnset)));
export_fig(fullpath(strFigurePath,sprintf('Q1c_SpikeBlockDecoding%s%s%s.pdf',strType,strSGS,strOnset)));
end

%% plot deciles by dur
matQuantDur = cellfun(@mean,cellQuantileDur_Dur);
figure;maxfig;
cellVars = {'Conf','Num','Rate','Tune','IFR'};
hSummary=subplot(2,3,6);
cla(hSummary);
hold(hSummary,'on');
plot(hSummary,[0 10],[0 0],'--','color',[0.5 0.5 0.5]);
title(hSummary,sprintf('%s, Bonferroni-corrected p & CI',strType));
for i=1:5
	if i == 1
		cellY = cellQuantileConf_Dur;
	elseif i == 2
		cellY = cellQuantileNum_Dur;
	elseif i == 3
		cellY = cellQuantileRate_Dur;
	elseif i == 4
		cellY = cellQuantileTune_Dur;
	elseif i == 5
		cellY = cellQuantileIFR_Dur;
	end
	strVar = cellVars{i};
	
	subplot(2,3,i)
	matQuantY = cellfun(@mean,cellY);
	matQuantY = zscore(matQuantY,[],2);
	matX = repmat(1:10,[intRecNum 1]);
	mdl = fitlm(matX(:),matQuantY(:));
	r=mdl.Coefficients.Estimate(2);
	r_SE=mdl.Coefficients.SE(2);
	dblAlpha = 0.05/5;
	r_CI = coefCI(mdl,dblAlpha);
	r_CI = r_CI(2,:);
	p=mdl.Coefficients.pValue(2)*5;
	
	hold('on')
	plot(matX',matQuantY','color',vecColGrey);
	plot(mean(matX,1),mean(matQuantY,1),'color',matCol(i,:));
	hold('off')
	xlabel('Duration decile of n-spike block');
	ylabel(sprintf('%s, z-scored per rec (%s)',strVar,getGreek('sigma')));
	title(sprintf('Lin reg, r=%.3f %s/decile, p=%.3e',r,getGreek('sigma'),p));
	
	%plot in summary
	errorbar(hSummary,(i*2-1),r,r-r_CI(1),r-r_CI(2),'x','color',matCol(i,:));
end
hold(hSummary,'off');
set(hSummary,'xtick',[1:2:9],'xticklabel',cellVars);
xtickangle(hSummary,45);
ylabel(hSummary,sprintf('Lin reg slope, mean +/- 95 CI (%s/decile)',getGreek('sigma')));
ylim(hSummary,[-0.4 0.4]);
fixfig;

if boolSaveFig
	%%
	export_fig(fullpath(strFigurePath,sprintf('Q1c_SpikeBlockCorrWithDur%s%s%s.tif',strType,strSGS,strOnset)));
	export_fig(fullpath(strFigurePath,sprintf('Q1c_SpikeBlockCorrWithDur%s%s%s.pdf',strType,strSGS,strOnset)));
end

%% plot deciles by conf
%{
matQuantConf = cellfun(@mean,cellQuantileConf_Conf);
figure;maxfig;
hSummary=subplot(2,3,6);
cla(hSummary);
hold(hSummary,'on');
plot(hSummary,[0 10],[0 0],'--','color',[0.5 0.5 0.5]);
title(hSummary,sprintf('%s, Bonferroni-corrected p & CI',strType));
for i=1:5
	if i == 1
		continue;
	elseif i == 2
		cellY = cellQuantileNum_Conf;
	elseif i == 3
		cellY = cellQuantileRate_Conf;
	elseif i == 4
		cellY = cellQuantileTune_Conf;
	elseif i == 5
		cellY = cellQuantileIFR_Conf;
	end
	strVar = cellVars{i};
	
	subplot(2,3,i)
	matQuantY = cellfun(@mean,cellY);
	matQuantY = zscore(matQuantY,[],2);
	matX = repmat(1:10,[intRecNum 1]);
	mdl = fitlm(matX(:),matQuantY(:));
	r=mdl.Coefficients.Estimate(2);
	r_SE=mdl.Coefficients.SE(2);
	dblAlpha = 0.05/4;
	r_CI = coefCI(mdl,dblAlpha);
	r_CI = r_CI(2,:);
	p=mdl.Coefficients.pValue(2)*4;
	
	hold('on')
	plot(matX',matQuantY','color',vecColGrey);
	plot(mean(matX,1),mean(matQuantY,1),'color',matCol(i,:));
	hold('off')
	xlabel('Confidence decile of n-spike block');
	ylabel(sprintf('%s, z-scored per rec (%s)',strVar,getGreek('sigma')));
	title(sprintf('Lin reg, r=%.3f %s/decile, p=%.3e',r,getGreek('sigma'),p));
	
	%plot in summary
	errorbar(hSummary,(i*2-1),r,r-r_CI(1),r-r_CI(2),'x','color',matCol(i,:));
end
hold(hSummary,'off');
set(hSummary,'xtick',[1:2:9],'xticklabel',cellVars);
xtickangle(hSummary,45);
ylabel(hSummary,sprintf('Lin reg slope, mean +/- 95 CI (%s/decile)',getGreek('sigma')));
ylim(hSummary,[-0.4 0.4]);
fixfig;

%%
if boolSaveFig
	export_fig(fullpath(strFigurePath,sprintf('Q1c_SpikeBlockCorrWithConf%s%s%s.tif',strType,strSGS,strOnset)));
	export_fig(fullpath(strFigurePath,sprintf('Q1c_SpikeBlockCorrWithConf%s%s%s.pdf',strType,strSGS,strOnset)));
end
%}
%% plot deciles by ifr
matQuantIFR = cellfun(@mean,cellQuantileIFR_IFR);
figure;maxfig;
hSummary=subplot(2,3,6);
cla(hSummary);
hold(hSummary,'on');
plot(hSummary,[0 10],[0 0],'--','color',[0.5 0.5 0.5]);
title(hSummary,sprintf('%s, Bonferroni-corrected p & CI',strType));
for i=1:4
	if i == 1
		cellY = cellQuantileConf_IFR;
	elseif i == 2
		cellY = cellQuantileNum_IFR;
	elseif i == 3
		cellY = cellQuantileRate_IFR;
	elseif i == 4
		cellY = cellQuantileTune_IFR;
	end
	strVar = cellVars{i};
	
	subplot(2,3,i)
	matQuantY = cellfun(@mean,cellY);
	matQuantY = zscore(matQuantY,[],2);
	matX = repmat(1:10,[intRecNum 1]);
	mdl = fitlm(matX(:),matQuantY(:));
	r=mdl.Coefficients.Estimate(2);
	r_SE=mdl.Coefficients.SE(2);
	dblAlpha = 0.05/4;
	r_CI = coefCI(mdl,dblAlpha);
	r_CI = r_CI(2,:);
	p=mdl.Coefficients.pValue(2)*4;
	
	hold('on')
	plot(matX',matQuantY','color',vecColGrey);
	plot(mean(matX,1),mean(matQuantY,1),'color',matCol(i,:));
	hold('off')
	xlabel('IFR decile of n-spike block');
	ylabel(sprintf('%s, z-scored per rec (%s)',strVar,getGreek('sigma')));
	title(sprintf('Lin reg, r=%.3f %s/decile, p=%.3e',r,getGreek('sigma'),p));
	
	%plot in summary
	errorbar(hSummary,(i*2-1),r,r-r_CI(1),r-r_CI(2),'x','color',matCol(i,:));
end
hold(hSummary,'off');
set(hSummary,'xtick',[1:2:9],'xticklabel',cellVars);
xtickangle(hSummary,45);
ylabel(hSummary,sprintf('Lin reg slope, mean +/- 95 CI (%s/decile)',getGreek('sigma')));
ylim(hSummary,[-0.4 0.4]);
fixfig;

%%
if boolSaveFig
	export_fig(fullpath(strFigurePath,sprintf('Q1c_SpikeBlockCorrWithIFR%s%s%s.tif',strType,strSGS,strOnset)));
	export_fig(fullpath(strFigurePath,sprintf('Q1c_SpikeBlockCorrWithIFR%s%s%s.pdf',strType,strSGS,strOnset)));
end
end
end