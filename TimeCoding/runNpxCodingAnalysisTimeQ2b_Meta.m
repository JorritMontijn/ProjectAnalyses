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

%% define parameters
%single rec plots
boolSingleRecPlots = false;

%% load data
dblRemOnset = 0; %remove onset period in seconds; 0.125 for sim, 0.25 for npx
%% onset string
if dblRemOnset == 0
	strOnset = '';
else
	strOnset = sprintf('%.2f',dblRemOnset);
end
sFiles = dir ([strTargetDataPath 'Q2bData*' strOnset '.mat']);
intRecNum = numel(sFiles);

%% pre-allocate

cellTypes =  {'Real','Poiss','ShuffTid','Shuff','PoissGain','Uniform'};
for intFile=1:intRecNum
	%% load
	sLoad=load(fullpath(sFiles(intFile).folder,sFiles(intFile).name));
	strRec = getFlankedBy(sFiles(intFile).name,'Q2bData','_g0_t0');
	
end

%% plot mean over recordings
%error plot matSlopes

