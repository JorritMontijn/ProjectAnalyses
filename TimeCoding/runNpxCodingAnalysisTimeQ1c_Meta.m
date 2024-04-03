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
cellRunTypes = {'RecTopo','SimDG18'};
vecMaxT = [1 0.5];
vecStimNum = [12 9];
vecRemOnsets = [0.25 0.125];
intRunType = 2;
strRunType = cellRunTypes{intRunType}; %Sim or ABI or Npx?
dblMaxT = vecMaxT(intRunType);
intStimNum = vecStimNum(intRunType);

boolFixedSpikeGroupSize = false;
intQuantiles = 10;
intGroupOrDecile = 1;
intOnsetType = 2;

vecBinsDur = 0:0.01:1;
vecBinsDurC = vecBinsDur(2:end) - diff(vecBinsDur(1:2))/2;
vecColGrey = [0.7 0.7 0.7];

for boolRemOnset=true%[false true]
	if boolRemOnset
		dblRemOnset = vecRemOnsets(intRunType);
	else
		dblRemOnset = 0;
	end
for intType=1:2
strType = cellTypes{intType};
if boolFixedSpikeGroupSize
	sFiles = dir ([strTargetDataPath 'Q1cData_' strRunType '*' strType '*Fixed*.mat']);
	strSGS = 'FixedSGS';
else
	sFiles = dir ([strTargetDataPath 'Q1cData_' strRunType '*' strType '*Var*.mat']);
	strSGS = 'VarSGS';
end
if dblRemOnset == 0
	strOnset = '';
	indWithOnset = cellfun(@(x) strcmp(x((end-7):end),[strOnset '.mat']), {sFiles.name}');
	sFiles(indWithOnset) = [];
else
	strOnset = sprintf('%.2f',dblRemOnset);
	indWithOnset = cellfun(@(x) strcmp(x((end-7):end),[strOnset '.mat']), {sFiles.name}');
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
plot(h1,[0 dblMaxT],([1 1]./intStimNum),'--','color',vecColChance);
hold(h1,'on')
title(h1,sprintf('%s, %s',strType,strSGS));

plot(h2,[0 100],([1 1]./intStimNum),'--','color',vecColChance);
hold(h2,'on')

%plot(h3,[0 10],([1 1]./intStimNum),'--','color',vecColChance);
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

	vecSpikeGroupFractionInterneurons = [sSpikeGroup.FractionInterneurons]';
	vecSpikeGroupAvgPrefDistToStim = [sSpikeGroup.AvgPrefDistToStim]';
	vecSpikeGroupAvgBandwidthOfCells = [sSpikeGroup.AvgBandwidthOfCells]';

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

		[vecMeanDurDeciles,vecSemDurDeciles,vecMeanFracInt,vecSemFracInt]=getQuantiles(vecSpikeGroupDuration,vecSpikeGroupFractionInterneurons,intQuantileNum);
		[vecMeanDurDeciles,vecSemDurDeciles,vecMeanPrefDist,vecSemPrefDist]=getQuantiles(vecSpikeGroupDuration,vecSpikeGroupAvgPrefDistToStim,intQuantileNum);
		[vecMeanDurDeciles,vecSemDurDeciles,vecMeanBandw,vecSemBandw]=getQuantiles(vecSpikeGroupDuration,vecSpikeGroupAvgBandwidthOfCells,intQuantileNum);
		vecMeanDurDeciles = vecMeanDurDeciles';
		vecMeanNum = vecMeanNum';
		vecMeanRate = vecMeanRate';
		vecMeanTune = vecMeanTune';
		vecMeanIFR = vecMeanIFR';
		vecMeanConf = vecMeanConf';
		vecMeanFracInt = vecMeanFracInt';
		vecMeanPrefDist = vecMeanPrefDist';
		vecMeanBandw = vecMeanBandw';

		%predict dur
		[r_Num_Dur,p_Num_Dur]=corr(vecMeanNum,vecMeanDurDeciles);
		[r_Rate_Dur,p_Rate_Dur]=corr(vecMeanRate,vecMeanDurDeciles);
		[r_Tune_Dur,p_Tune_Dur]=corr(vecMeanTune,vecMeanDurDeciles);
		[r_IFR_Dur,p_IFR_Dur]=corr(vecMeanIFR,vecMeanDurDeciles);
		[r_Conf_Dur,p_Conf_Dur]=corr(vecMeanConf,vecMeanDurDeciles);

		[r_FrInt_Dur,p_FrInt_Dur]=corr(vecMeanFracInt,vecMeanDurDeciles);
		[r_PrDi_Dur,p_PrDi_Dur]=corr(vecMeanPrefDist,vecMeanDurDeciles);
		[r_BndW_Dur,p_BndW_Dur]=corr(vecMeanBandw,vecMeanDurDeciles);

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
		[r_FrInt_Dur,p_FrInt_Dur]=corr(vecSpikeGroupFractionInterneurons,vecSpikeGroupDuration);
		[r_PrDi_Dur,p_PrDi_Dur]=corr(vecSpikeGroupAvgPrefDistToStim,vecSpikeGroupDuration);
		[r_BndW_Dur,p_BndW_Dur]=corr(vecSpikeGroupAvgBandwidthOfCells,vecSpikeGroupDuration);

		tbl = table(vecSpikeGroupNumOfCells,vecSpikeGroupAvgRateOfCells,vecSpikeGroupAvgTuningOfCells,vecSpikeGroupAvgIFR,vecSpikeGroupConfidence,vecSpikeGroupDuration,...
			'VariableNames',{'CellNumber','AvgCellRate','AvgCellTuning','AvgPopIFR','Conf','Dur'});
		mdl_preddur = fitlm(tbl,'linear');

		%predict conf
		[r_Num_Conf,p_Num_Conf]=corr(vecSpikeGroupNumOfCells,vecSpikeGroupConfidence);
		[r_Rate_Conf,p_Rate_Conf]=corr(vecSpikeGroupAvgRateOfCells,vecSpikeGroupConfidence);
		[r_Tune_Conf,p_Tune_Conf]=corr(vecSpikeGroupAvgTuningOfCells,vecSpikeGroupConfidence);
		[r_IFR_Conf,p_IFR_Conf]=corr(vecSpikeGroupAvgIFR,vecSpikeGroupConfidence);
		[r_FrInt_Conf,p_FrInt_Conf]=corr(vecSpikeGroupFractionInterneurons,vecSpikeGroupConfidence);
		[r_PrDi_Conf,p_PrDi_Conf]=corr(vecSpikeGroupAvgPrefDistToStim,vecSpikeGroupConfidence);
		[r_BndW_Conf,p_BndW_Conf]=corr(vecSpikeGroupAvgBandwidthOfCells,vecSpikeGroupConfidence);

		tbl = table(vecSpikeGroupNumOfCells,vecSpikeGroupAvgRateOfCells,vecSpikeGroupAvgTuningOfCells,vecSpikeGroupAvgIFR,vecSpikeGroupConfidence,...
			'VariableNames',{'CellNumber','AvgCellRate','AvgCellTuning','AvgPopIFR','Conf'});
		mdl_predconf = fitlm(tbl,'linear');

		%predict ifr
		[r_Num_IFR,p_Num_IFR]=corr(vecSpikeGroupNumOfCells,vecSpikeGroupAvgIFR);
		[r_Rate_IFR,p_Rate_IFR]=corr(vecSpikeGroupAvgRateOfCells,vecSpikeGroupAvgIFR);
		[r_Tune_IFR,p_Tune_IFR]=corr(vecSpikeGroupAvgTuningOfCells,vecSpikeGroupAvgIFR);
		[r_Conf_IFR,p_Conf_IFR]=corr(vecSpikeGroupConfidence,vecSpikeGroupAvgIFR);
		[r_FrInt_IFR,p_FrInt_IFR]=corr(vecSpikeGroupFractionInterneurons,vecSpikeGroupAvgIFR);
		[r_PrDi_IFR,p_PrDi_IFR]=corr(vecSpikeGroupAvgPrefDistToStim,vecSpikeGroupAvgIFR);
		[r_BndW_IFR,p_BndW_IFR]=corr(vecSpikeGroupAvgBandwidthOfCells,vecSpikeGroupAvgIFR);

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
	vecR_FrIn_Dur(intFile) = r_FrInt_Dur;
	vecR_PrDi_Dur(intFile) = r_PrDi_Dur;
	vecR_Bndw_Dur(intFile) = r_BndW_Dur;

	% 	vecR_CombDur_CellNum(intFile) = table2array(mdl_preddur.Coefficients(2,1));
	% 	vecR_CombDur_Rate(intFile) = table2array(mdl_preddur.Coefficients(3,1));
	% 	vecR_CombDur_Tune(intFile) = table2array(mdl_preddur.Coefficients(4,1));
	% 	vecR_CombDur_IFR(intFile) = table2array(mdl_preddur.Coefficients(5,1));
	% 	vecR_CombDur_Conf(intFile) = table2array(mdl_preddur.Coefficients(6,1));

	%conf
	vecR_Num_Conf(intFile) = r_Num_Conf;
	vecR_Rate_Conf(intFile) = r_Rate_Conf;
	vecR_Tune_Conf(intFile) = r_Tune_Conf;
	vecR_IFR_Conf(intFile) = r_IFR_Conf;
	vecR_FrIn_Conf(intFile) = r_FrInt_Conf;
	vecR_PrDi_Conf(intFile) = r_PrDi_Conf;
	vecR_Bndw_Conf(intFile) = r_BndW_Conf;

	% 	vecR_CombConf_CellNum(intFile) = table2array(mdl_predconf.Coefficients(2,1));
	% 	vecR_CombConf_Rate(intFile) = table2array(mdl_predconf.Coefficients(3,1));
	% 	vecR_CombConf_Tune(intFile) = table2array(mdl_predconf.Coefficients(4,1));
	% 	vecR_CombConf_IFR(intFile) = table2array(mdl_predconf.Coefficients(5,1));

	%ifr
	vecR_Num_IFR(intFile) = r_Num_IFR;
	vecR_Rate_IFR(intFile) = r_Rate_IFR;
	vecR_Tune_IFR(intFile) = r_Tune_IFR;
	vecR_Conf_IFR(intFile) = r_Conf_IFR;
	vecR_FrIn_IFR(intFile) = r_FrInt_IFR;
	vecR_PrDi_IFR(intFile) = r_PrDi_IFR;
	vecR_Bndw_IFR(intFile) = r_BndW_IFR;

	%vecR_CombIFR_Rate(intFile) = table2array(mdl_predifr.Coefficients(2,1));
	%vecR_CombIFR_Tune(intFile) = table2array(mdl_predifr.Coefficients(3,1));

	%vecR_Group_Rate_Tune(intFile) = rGroup_Rate_Tune;

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
	vecSpikeGroupFractionInterneurons= [sSpikeGroup.FractionInterneurons]';
	vecSpikeGroupAvgPrefDistToStim = [sSpikeGroup.AvgPrefDistToStim]';
	vecSpikeGroupAvgBandwidthOfCells = [sSpikeGroup.AvgBandwidthOfCells]';

	%by dur
	intQuantileNum = 10;
	[~,~,~,~,~,cellValsDur_Dur,cellValsConf_Dur]=getQuantiles(vecSpikeGroupDuration,vecSpikeGroupConfidence,intQuantileNum);
	[~,~,~,~,~,~,cellValsNum_Dur]=getQuantiles(vecSpikeGroupDuration,vecSpikeGroupNumOfCells,intQuantileNum);
	[~,~,~,~,~,~,cellValsRate_Dur]=getQuantiles(vecSpikeGroupDuration,vecSpikeGroupAvgRateOfCells,intQuantileNum);
	[~,~,~,~,~,~,cellValsTune_Dur]=getQuantiles(vecSpikeGroupDuration,vecSpikeGroupAvgTuningOfCells,intQuantileNum);
	[~,~,~,~,~,~,cellValsIFR_Dur]=getQuantiles(vecSpikeGroupDuration,vecSpikeGroupAvgIFR,intQuantileNum);
	[~,~,~,~,~,~,cellValsFrIn_Dur]=getQuantiles(vecSpikeGroupDuration,vecSpikeGroupFractionInterneurons,intQuantileNum);
	[~,~,~,~,~,~,cellValsPrDi_Dur]=getQuantiles(vecSpikeGroupDuration,vecSpikeGroupAvgPrefDistToStim,intQuantileNum);
	[~,~,~,~,~,~,cellValsBndw_Dur]=getQuantiles(vecSpikeGroupDuration,vecSpikeGroupAvgBandwidthOfCells,intQuantileNum);
	cellQuantileDur_Dur(intFile,:) = cellValsDur_Dur;
	cellQuantileConf_Dur(intFile,:) = cellValsConf_Dur;
	cellQuantileNum_Dur(intFile,:) = cellValsNum_Dur;
	cellQuantileRate_Dur(intFile,:) = cellValsRate_Dur;
	cellQuantileTune_Dur(intFile,:) = cellValsTune_Dur;
	cellQuantileIFR_Dur(intFile,:) = cellValsIFR_Dur;
	cellQuantileFrIn_Dur(intFile,:) = cellValsFrIn_Dur;
	cellQuantilePrDi_Dur(intFile,:) = cellValsPrDi_Dur;
	cellQuantileBndw_Dur(intFile,:) = cellValsBndw_Dur;

	%by conf
	[~,~,~,~,~,cellValsConf_Conf,cellValsNum_Conf]=getQuantiles(vecSpikeGroupConfidence,vecSpikeGroupNumOfCells,intQuantileNum);
	[~,~,~,~,~,~,cellValsRate_Conf]=getQuantiles(vecSpikeGroupConfidence,vecSpikeGroupAvgRateOfCells,intQuantileNum);
	[~,~,~,~,~,~,cellValsTune_Conf]=getQuantiles(vecSpikeGroupConfidence,vecSpikeGroupAvgTuningOfCells,intQuantileNum);
	[~,~,~,~,~,~,cellValsIFR_Conf]=getQuantiles(vecSpikeGroupConfidence,vecSpikeGroupAvgIFR,intQuantileNum);
	[~,~,~,~,~,~,cellValsFrIn_Conf]=getQuantiles(vecSpikeGroupConfidence,vecSpikeGroupFractionInterneurons,intQuantileNum);
	[~,~,~,~,~,~,cellValsPrDi_Conf]=getQuantiles(vecSpikeGroupConfidence,vecSpikeGroupAvgPrefDistToStim,intQuantileNum);
	[~,~,~,~,~,~,cellValsBndw_Conf]=getQuantiles(vecSpikeGroupConfidence,vecSpikeGroupAvgBandwidthOfCells,intQuantileNum);
	cellQuantileConf_Conf(intFile,:) = cellValsConf_Conf;
	cellQuantileNum_Conf(intFile,:) = cellValsNum_Conf;
	cellQuantileRate_Conf(intFile,:) = cellValsRate_Conf;
	cellQuantileTune_Conf(intFile,:) = cellValsTune_Conf;
	cellQuantileIFR_Conf(intFile,:) = cellValsIFR_Conf;
	cellQuantileFrIn_Conf(intFile,:) = cellValsFrIn_Conf;
	cellQuantilePrDi_Conf(intFile,:) = cellValsPrDi_Conf;
	cellQuantileBndw_Conf(intFile,:) = cellValsBndw_Conf;

	%by ifr
	[~,~,~,~,~,cellValsIFR_IFR,cellValsConf_IFR]=getQuantiles(vecSpikeGroupAvgIFR,vecSpikeGroupConfidence,intQuantileNum);
	[~,~,~,~,~,~,cellValsNum_IFR]=getQuantiles(vecSpikeGroupAvgIFR,vecSpikeGroupNumOfCells,intQuantileNum);
	[~,~,~,~,~,~,cellValsRate_IFR]=getQuantiles(vecSpikeGroupAvgIFR,vecSpikeGroupAvgRateOfCells,intQuantileNum);
	[~,~,~,~,~,~,cellValsTune_IFR]=getQuantiles(vecSpikeGroupAvgIFR,vecSpikeGroupAvgTuningOfCells,intQuantileNum);
	[~,~,~,~,~,~,cellValsFrIn_IFR]=getQuantiles(vecSpikeGroupAvgIFR,vecSpikeGroupFractionInterneurons,intQuantileNum);
	[~,~,~,~,~,~,cellValsPrDi_IFR]=getQuantiles(vecSpikeGroupAvgIFR,vecSpikeGroupAvgPrefDistToStim,intQuantileNum);
	[~,~,~,~,~,~,cellValsBndw_IFR]=getQuantiles(vecSpikeGroupAvgIFR,vecSpikeGroupAvgBandwidthOfCells,intQuantileNum);
	cellQuantileConf_IFR(intFile,:) = cellValsConf_IFR;
	cellQuantileNum_IFR(intFile,:) = cellValsNum_IFR;
	cellQuantileRate_IFR(intFile,:) = cellValsRate_IFR;
	cellQuantileTune_IFR(intFile,:) = cellValsTune_IFR;
	cellQuantileIFR_IFR(intFile,:) = cellValsIFR_IFR;
	cellQuantileFrIn_IFR(intFile,:) = cellValsFrIn_IFR;
	cellQuantilePrDi_IFR(intFile,:) = cellValsPrDi_IFR;
	cellQuantileBndw_IFR(intFile,:) = cellValsBndw_IFR;

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
cellSuffices = {'Dur','Conf','IFR'};
cellVars = {'Num','Rate','Tune','IFR','Conf','FrIn','PrDi','Bndw'};
matP = nan(3,numel(cellVars));
for i=1:3
	%transform variable names and calculate p
	strSuffix = cellSuffices{i};
	for intVar=1:numel(cellVars)
		strVar = cellVars{intVar};
		if strcmp(strVar,strSuffix)
			eval(['vecR_' strVar ' = nan(size(vecR_' cellVars{1} '_' strSuffix '));']);
		else
			eval(['[h,matP(i,intVar)] = ttest(vecR_' strVar '_' strSuffix ');']);
			eval(['vecR_' strVar ' = vecR_' strVar '_' strSuffix ';']);
		end
	end
	vecP=matP(i,:);
	[h, crit_p, adj_p]=fdr_bh(vecP,0.05);

	%% plot
	subplot(2,3,3+i);cla
	hold on;
	plot([0 15],[0 0],'--','color',[0.5 0.5 0.5]);
	matCol = lines(8);
	matCol = matCol([2 1 3:size(matCol,1)],:);
	matCol(3,:) = 0;
	matColP = 1-((1-matCol)./2);

	intRecNum = numel(vecR_Num);
	vecX = ones(size(vecR_Num));
	h2=swarmchart(1*vecX,vecR_Num,[],matColP(1,:),'filled');
	swarmchart(3*vecX,vecR_Rate,[],matColP(2,:),'filled');
	swarmchart(5*vecX,vecR_Tune,[],matColP(3,:),'filled');
	swarmchart(7*vecX,vecR_IFR,[],matColP(4,:),'filled');
	swarmchart(9*vecX,vecR_Conf,[],matColP(5,:),'filled');

	swarmchart(11*vecX,vecR_FrIn,[],matColP(6,:),'filled');
	swarmchart(13*vecX,vecR_PrDi,[],matColP(7,:),'filled');
	swarmchart(15*vecX,vecR_Bndw,[],matColP(8,:),'filled');

	errorbar(1,mean(vecR_Num),std(vecR_Num)./sqrt(intRecNum),'x','color',matCol(1,:));
	errorbar(3,mean(vecR_Rate),std(vecR_Rate)./sqrt(intRecNum),'x','color',matCol(2,:));
	errorbar(5,mean(vecR_Tune),std(vecR_Tune)./sqrt(intRecNum),'x','color',matCol(3,:));
	errorbar(7,mean(vecR_IFR),std(vecR_IFR)./sqrt(intRecNum),'x','color',matCol(4,:));
	errorbar(9,mean(vecR_Conf),std(vecR_Conf)./sqrt(intRecNum),'x','color',matCol(5,:));

	errorbar(11,mean(vecR_FrIn),std(vecR_FrIn)./sqrt(intRecNum),'x','color',matCol(6,:));
	errorbar(13,mean(vecR_PrDi),std(vecR_PrDi)./sqrt(intRecNum),'x','color',matCol(7,:));
	errorbar(15,mean(vecR_Bndw),std(vecR_Bndw)./sqrt(intRecNum),'x','color',matCol(8,:));

	hold off
	set(gca,'xtick',[1:2:15],'xticklabel',{'# of cells','Cell Rate','Cell Tune','Pop Rate','Conf',...
		'% InterN','Pref-d','BndW'});
	xtickangle(gca,45);
	ylabel(['Correlation with group ' strSuffix ' (r)']);
	title(sprintf('p: #=%.2e,R=%.3f,T=%.2e,IFR=%.2e,\nC=%.2e,FrI=%.2e,PrDi=%.2e,B=%.2e',adj_p))
end
fixfig;

%%
if boolSaveFig
	drawnow;
	export_fig(fullpath(strFigurePath,sprintf('Q1c_SpikeBlockDecoding%s%s%s%s.tif',strRunType,strType,strSGS,strOnset)));
	export_fig(fullpath(strFigurePath,sprintf('Q1c_SpikeBlockDecoding%s%s%s%s.pdf',strRunType,strType,strSGS,strOnset)));
end

%% plot deciles per var
%test vs 0
cellSuffices = {'Dur','Conf','IFR'};
cellVars = {'Num','Rate','Tune','IFR','Conf','FrIn','PrDi','Bndw'};
matP = nan(3,numel(cellVars));
for i=1:3
	%prep plot
	figure;maxfig;
	hSummary=subplot(2,5,10);
	cla(hSummary);
	hold(hSummary,'on');
	plot(hSummary,[0 (2*numel(cellVars))-1],[0 0],'--','color',[0.5 0.5 0.5]);
	title(hSummary,sprintf('%s-%s, Bonferroni-corrected p & CI',strRunType,strType));

	%transform variable names and calculate p
	strSuffix = cellSuffices{i};
	for intVar=1:numel(cellVars)
		strVar = cellVars{intVar};
		eval(['cellY = cellQuantile' strVar '_' strSuffix ';']);


		subplot(2,5,intVar)
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
		plot(mean(matX,1),mean(matQuantY,1),'color',matCol(intVar,:));
		hold('off')
		xlabel([strSuffix ' decile of n-spike block']);
		ylabel(sprintf('%s, z-scored per rec (%s)',strVar,getGreek('sigma')));
		title(sprintf('OLS, y=%.2f %s/x, p=%.2e',r,getGreek('sigma'),p));

		%plot in summary
		errorbar(hSummary,(intVar*2-1),r,r-r_CI(1),r-r_CI(2),'x','color',matCol(intVar,:));
	end
	hold(hSummary,'off');
	set(hSummary,'xtick',[1:2:(2*numel(cellVars))],'xticklabel',cellVars);
	xtickangle(hSummary,45);
	ylabel(hSummary,sprintf('Lin reg slope, mean +/- 95 CI (%s/decile)',getGreek('sigma')));
	ylim(hSummary,[-0.4 0.4]);
	fixfig;

	if boolSaveFig
		%%
		export_fig(fullpath(strFigurePath,sprintf('Q1c_SpikeBlockCorr%s%s%s%s%s.tif',strType,strRunType,strSuffix,strSGS,strOnset)));
		export_fig(fullpath(strFigurePath,sprintf('Q1c_SpikeBlockCorr%s%s%s%s%s.pdf',strType,strRunType,strSuffix,strSGS,strOnset)));
	end
end
end
end