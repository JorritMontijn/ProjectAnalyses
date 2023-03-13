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
boolHome = true;
if boolHome
	strDataPath = 'F:\Data\Processed\Neuropixels\';
	strFigurePathSR = 'F:\Drive\PopTimeCoding\single_recs';
	strFigurePath = 'F:\Drive\PopTimeCoding\figures\';
	strTargetDataPath = 'F:\Drive\PopTimeCoding\data\';
else
	strFigurePath = 'D:\Data\Results\PopTimeCoding\figures\';
	strTargetDataPath = 'D:\Data\Results\PopTimeCoding\data\';
end

%% load data
sFiles = dir ([strTargetDataPath 'Q3Data*.mat']);
strArea = 'V1';
intRecNum = numel(sFiles);

%% pre-allocate
vecBins = -25:2.5:25;
matCounts = [];
vecAggRelVarMean = [];
vecAggRelVarNorm = [];
matAggPercDiff = [];
cellRec = {};

for intFile=1:intRecNum
	%% load
	load(fullpath(sFiles(intFile).folder,sFiles(intFile).name));
	strRecFile = sFiles(intFile).name;
	strRecFile = strRecFile((1+numel('Q3Data_Rec')):(16+numel('Q3Data_Rec')));
	
	
	vecPercDiff = 100*((vecRelVarNorm./vecRelVarMean)-1);
	vecBinsC = vecBins(2:end)-median(diff(vecBins))/2;
	vecCounts =histcounts(vecPercDiff,vecBins);
	
	matAggPercDiff(end+1,:) = vecPercDiff;
	matCounts(end+1,:) = vecCounts;
	cellRec(end+1) = {strRecFile};
	vecAggRelVarMean(end+1) = dblRelVarMean;
	vecAggRelVarNorm(end+1) = dblRelVarNorm;
end

subplot(2,3,1)
bar(vecBinsC,sum(matCounts,1),'hist');
xlabel('Variability of pop act. as Euclidian vs mean distance ((L2/L1) - 1) (%)');
title('Trial-to-trial variability of pop act per stim type');
ylabel('Number of stimulus types');
fixfig;
		
vecPercDiff2 = 100*((vecAggRelVarNorm./vecAggRelVarMean)-1);
vecCounts2 = histcounts(vecPercDiff2,vecBins);
subplot(2,3,2)
bar(vecBinsC,vecCounts2,'hist');
xlabel('Variability of pop act. as Euclidian vs mean distance ((L2/L1) - 1) (%)');
title('Trial-to-trial variability of pop act per stim type');
ylabel('Number of stimulus types');
fixfig;