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
clear all;close all;
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
sFiles = dir ([strTargetDataPath 'Q4Data*.mat']);
strArea = 'V1';
intRecNum = numel(sFiles);

%% pre-allocate
vecR2_M = nan(1,intRecNum);
vecR2_G = nan(1,intRecNum);
vecR2_R = nan(1,intRecNum);
vecR2_adjusted_M = nan(1,intRecNum);
vecR2_adjusted_G = nan(1,intRecNum);
vecR2_adjusted_R = nan(1,intRecNum);
vecCorrMeanPupil = nan(1,intRecNum);
vecCorrGainPupil = nan(1,intRecNum);
indRem = false(1,intRecNum);
for intFile=1:intRecNum
	%% load
	load(fullpath(sFiles(intFile).folder,sFiles(intFile).name));
	strRecFile = sFiles(intFile).name;
	strRecFile = strRecFile((1+numel('Q4Data_Rec')):(16+numel('Q4Data_Rec')));
	%save
	if dblP_M < 0.05
	vecR2_M(intFile) = dblR2_M;
	vecR2_G(intFile) = dblR2_G;
	vecR2_R(intFile) = dblR2_R;
	vecR2_adjusted_M(intFile) = dblR2_adjusted_M;
	vecR2_adjusted_G(intFile) = dblR2_adjusted_G;
	vecR2_adjusted_R(intFile) = dblR2_adjusted_R;
	vecCorrMeanPupil(intFile) = dblCorrMeanPupil;
	vecCorrGainPupil(intFile) = dblCorrGainPupil;
	else
		indRem(intFile) = true;
	end
end
vecR2_M(indRem) = [];
vecR2_G(indRem) = [];
vecR2_R(indRem) = [];
vecR2_adjusted_M(indRem) = [];
vecR2_adjusted_G(indRem) = [];
vecR2_adjusted_R(indRem) = [];
vecCorrMeanPupil(indRem) = [];
vecCorrGainPupil(indRem) = [];

%%
figure;maxfig;
subplot(2,3,1)
dblMeanR2_M = mean(vecR2_M);
dblMeanR2_G = mean(vecR2_G);
dblMeanR2_R = mean(vecR2_R);
errorbar(1:3,[dblMeanR2_M dblMeanR2_G dblMeanR2_R],[std(vecR2_M) std(vecR2_G) std(vecR2_R)]./sqrt(intRecNum),'x');
set(gca,'xtick',1:3,'xticklabel',{'Mean','Gain','Reg'});
ylabel('Uncorrected R^2');
xlim([0.5 3.5]);

subplot(2,3,2)
dblMeanR2_adjusted_M = mean(vecR2_adjusted_M);
dblMeanR2_adjusted_G = mean(vecR2_adjusted_G);
dblMeanR2_adjusted_R = mean(vecR2_adjusted_R);
errorbar(1:3,[dblMeanR2_adjusted_M dblMeanR2_adjusted_G dblMeanR2_adjusted_R],[std(vecR2_adjusted_M) std(vecR2_adjusted_G) std(vecR2_adjusted_R)]./sqrt(intRecNum),'x');
set(gca,'xtick',1:3,'xticklabel',{'Mean','Gain','Reg'});
ylabel('Adjusted R^2');
xlim([0.5 3.5]);

subplot(2,3,3)
dblMeanR2_adjusted_M = mean(vecR2_adjusted_M);
dblMeanR2_adjusted_G = mean(vecR2_adjusted_G);
dblMeanR2_adjusted_R = mean(vecR2_adjusted_R);
errorbar(1:2,[mean(vecCorrMeanPupil) mean(vecCorrGainPupil)],[std(vecCorrMeanPupil) std(vecCorrGainPupil)]./sqrt(intRecNum),'x');
set(gca,'xtick',1:3,'xticklabel',{'Mean','Gain'});
ylabel('Correlation with pupil size');
xlim([0.5 2.5]);
fixfig;

drawnow;
export_fig(fullpath(strFigurePath,sprintf('Q4_VariabilityPrediction.tif')));
export_fig(fullpath(strFigurePath,sprintf('Q4_VariabilityPrediction.pdf')));
