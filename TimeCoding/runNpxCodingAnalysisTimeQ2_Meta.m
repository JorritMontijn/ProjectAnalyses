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
sFiles = dir ([strTargetDataPath 'Q2Data*.mat']);
strArea = 'V1';
intRecNum = numel(sFiles);
for intFile=1:intRecNum
	%% load
	load(fullpath(sFiles(intFile).folder,sFiles(intFile).name));
	strRec = sFiles(intFile).name;
	strRec = strRec((1+numel('Q2Data_Rec')):(16+numel('Q2Data_Rec')));
	
	%% create derived variables
	%get data
	vecAllISI = cell2vec(cellISI_perTrial);
	cellIFR_Reduced = cellfun(@(x) x(2:(end-2)),cellIFR_perTrial,'UniformOutput',false);
	vecAllIFR = cell2vec(cellIFR_Reduced);
	
	%% save data
	dblStartT
	strRec
	cellIFR_perTrial
	cellTimeIFR_perTrial
	cellISI_perTrial
	vecRperTrial
	vecSperTrial
	vecHperTrial
	vecLperTrial
	vecMperTrial
	vecPopSparseness
	vecPopSparsenessS
	vecSperTrial_S
	vecHperTrial_S
	vecLperTrial_S
	vecMperTrial_S
	vecSperTrial_SN
	vecHperTrial_SN
	vecLperTrial_SN
	vecMperTrial_SN
	vecSperTrial_SS
	vecHperTrial_SS
	vecLperTrial_SS
	vecMperTrial_SS
	
	
	%% single plot 1
	figure
	subplot(2,3,1)
	dblBinSizeISI = 1/1000;
	vecBinEdgesISI = 0:dblBinSizeISI:0.05;
	vecBinCentersISI = vecBinEdgesISI(2:end)-dblBinSizeISI/2;
	dblLambda = 1./mean(vecAllISI);
	vecExpPdf = dblLambda.*exp(-dblLambda.*vecBinCentersISI);
	vecCounts = histcounts(vecAllISI,vecBinEdgesISI);
	hold on
	plot(vecBinCentersISI*1000,vecCounts./sum(vecCounts(:)))
	plot(vecBinCentersISI*1000,vecExpPdf./sum(vecExpPdf(:)))
	hold off
	set(gca,'yscale','log');
	xlabel('Inter-spike interval (ms)');
	ylabel('Normalized count (n)');
	legend({'Observed','Theory (Exponential)'});
	title('Population spiking dynamics');
	fixfig;
	
	
	subplot(2,3,2)
	%vars
	dblMaxT = 1000*0.05;
	dblStep = 1000*0.002;
	vecBinD = 0:dblStep:dblMaxT;
	vecBinD_c = vecBinD(2:end)-dblStep/2;
	
	%real
	vecD1 = 1000*vecAllISI(1:(end-1));
	vecD2 = 1000*vecAllISI(2:end);
	[r,p]=corr(vecD1,vecD2);
	[vecCounts,vecMeans,vecSDs,cellVals,cellIDs] = makeBins(vecD1,vecD2,vecBinD);
	
	%renewal process
	vecAllISI_R = exprnd(1/dblLambda,size(vecAllISI));
	vecD1_R = 1000*vecAllISI_R(1:(end-1));
	vecD2_R = 1000*vecAllISI_R(2:end);
	[r,p]=corr(vecD1,vecD2);
	[vecCounts_R,vecMeans_R,vecSDs_R,cellVals_R,cellIDs_R] = makeBins(vecD1_R,vecD2_R,vecBinD);
	
	%shuffled
	vecAllISI_S = vecAllISI(randperm(numel(vecAllISI)));
	vecD1_S = 1000*vecAllISI_S(1:(end-1));
	vecD2_S = 1000*vecAllISI_S(2:end);
	[rS,pS]=corr(vecD1,vecD2);
	[vecCounts_S,vecMeans_S,vecSDs_S,cellVals_S,cellIDs_S] = makeBins(vecD1_S,vecD2_S,vecBinD);
	
	%plot
	errorbar(vecBinD_c,vecMeans,vecSDs./sqrt(vecCounts));
	hold on
	errorbar(vecBinD_c,vecMeans_R-0.02,vecSDs_R./sqrt(vecCounts_R));
	errorbar(vecBinD_c,vecMeans_S+0.02,vecSDs_S./sqrt(vecCounts_S));
	
	xlabel('ISI spikes i,i+1 (ms)');
	ylabel('ISI spikes i+1,i+2 (ms)');
	title(sprintf('ISI correlation r(d(i,j),d(i+1,j+1)), r=%.3f, p=%.3f',r,p));
	fixfig;
	
	legend({'Observed','Theory (Exponential)','Shuffled'},'location','best');
	
	%% single plot 2
	
	
end

%% plot mean over recordings
