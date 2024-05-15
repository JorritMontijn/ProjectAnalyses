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
sFiles = dir ([strTargetDataPath 'Q1Data*.mat']);
boolSinglePlots = false;
strArea = 'V1';
intRecNum = numel(sFiles);
intBinNum = 300;
vecDecP = nan(1,intRecNum);
matDecPerf = nan(intBinNum,intRecNum);
matRateRaw = nan(intBinNum,intRecNum);
matAcrossTimeDecoderAgg = nan(intBinNum,intBinNum,intRecNum);
for intFile=1:intRecNum
	%% load
	load(fullpath(sFiles(intFile).folder,sFiles(intFile).name));
	strRec = sFiles(intFile).name;
	strRec = strRec((1+numel('Q1Data_Rec')):(16+numel('Q1Data_Rec')));
	
	%% create derived variables
	[vecTrialTypeIdx,vecUnique,vecPriorDistribution,cellSelect,vecRepetition] = val2idx(vecOri180);
	dblChance = 1/numel(vecPriorDistribution);
	indUseForTest = vecStimTime>0.05 & vecStimTime < 0.15;
	vecDec = vecDecConf;
	[h,dblP] = ttest(vecDecConf(indUseForTest),dblChance);
	vecDecP(intFile) = dblP;
	matAcrossTimeDecoder(tril(true(size(matAcrossTimeDecoder))) & triu(true(size(matAcrossTimeDecoder)))) = vecDec;
	
	if boolSinglePlots
	%% plot
	figure;maxfig;
	subplot(2,3,1)
	plot(vecStimTime,vecSpikesPerBin./dblBinWidth)
	title('Spikes per bin')
	xlabel('Time after onset (s)');
	ylabel('Binned population spiking rate (Hz)');
	fixfig;
	
	subplot(2,3,2)
	plot(vecStimTime,vecDec(:,end)');
	hold on
	plot([vecStimTime(1) vecStimTime(end)],[dblChance dblChance],'--','color',[0.5 0.5 0.5]);
	hold off
	title(sprintf('Dec perf %s; %s',strRec,strArea),'interpreter','none')
	xlabel('Time after onset (s)');
	ylabel('Fraction correct decoded');
	fixfig;
	
	subplot(2,3,3)
	plot(vecStimTime,(vecDec(:,end)./dblChance)'./(vecSpikesPerBin./dblBinWidth))
	title('Dec perf / spike')
	xlabel('Time after onset (s)');
	ylabel('Performance/spike');
	fixfig;
	
	hS=subplot(2,3,4);
	cMap=colormap(hS,circcol);
	hB=colorbar;
	set(gca,'clim',[min(vecStimTime) max(vecStimTime)]);
	h=cline([vecSpikesPerBin(:); vecSpikesPerBin(1)]./dblBinWidth,[vecDec(:,end); vecDec(1,end)],[vecStimTime(:); vecStimTime(1)]);
	set(h,'LineWidth',2);
	hold on
	plot([min(vecSpikesPerBin) max(vecSpikesPerBin)]./dblBinWidth,[dblChance dblChance],'--','color',[0.5 0.5 0.5]);
	hold off
	hB.Label.String = 'Time (s)';
	xlabel('Binned population spiking rate (Hz)');
	ylabel('Decoding performance');
	title(sprintf('Decoding>chance p-value=%.3e',dblP));
	fixfig;
	
	hS2=subplot(2,3,5);
	colormap(hS2,parula);
	imagesc(vecStimTime,vecStimTime,matAcrossTimeDecoder);
	hB2=colorbar;
	hB2.Label.String = 'Decoding performance';
	axis xy
	xlabel('Training bin');
	ylabel('Testing bin');
	fixfig;grid off
	
	drawnow;
	export_fig(fullpath(strFigurePathSR,sprintf('Q1%s.tif',strRec)));
	export_fig(fullpath(strFigurePathSR,sprintf('Q1%s.pdf',strRec)));
	end
	%% save data
	matDecPerf(:,intFile) = vecDec;
	matRateRaw(:,intFile) = (vecSpikesPerBin./dblBinWidth);
	matAcrossTimeDecoderAgg(:,:,intFile) = matAcrossTimeDecoder;
end

%% plot mean over recordings
vecCorrP = bonf_holm(vecDecP);
vecAvgPerf=mean(matDecPerf(indUseForTest,:),1);
indUseRecs = vecCorrP<0.01 & vecAvgPerf>0.09;

matRate = matRateRaw(:,indUseRecs);
matPerf = matDecPerf(:,indUseRecs);
matConf = matAcrossTimeDecoderAgg(:,:,indUseRecs);
matConfOverTime = nan(size(matPerf));
intB = size(matConfOverTime,1);
for i=1:size(matConf,3)
	matSq = matAcrossTimeDecoderAgg(:,:,i);
	matConfOverTime(:,i) = matSq(diag(diag(true(intB))));
end
vecMeanRate = mean(matRate,2);
vecMeanPerf = mean(matPerf,2);
vecMeanConf = mean(matConfOverTime,2);
matMeanConf = mean(matConf,3);
intUseRecs = size(matRate,2);

figure;maxfig;
%avg rate/decoding/confusion matrix


subplot(2,3,1)
plot(vecStimTime,vecMeanRate);
title(sprintf('Average spiking per %.1f ms bin, n=%d recs',dblBinWidth*1000,intUseRecs))
xlabel('Time after onset (s)');
ylabel('Binned population spiking rate (Hz)');
fixfig;

subplot(2,3,2)
plot(vecStimTime,vecMeanPerf,'color','k');
hold on
plot([vecStimTime(1) vecStimTime(end)],[dblChance dblChance],'--','color',[0.5 0.5 0.5]);
hold off
title(sprintf('Dec perf'),'interpreter','none')
xlabel('Time after onset (s)');
ylabel('Decoder accuracy');
fixfig;


hS2=subplot(2,3,3);
colormap(hS2,parula);
imagesc(vecStimTime,vecStimTime,matMeanConf);
hB2=colorbar;
hB2.Label.String = 'Decoder confidence';
axis xy
xlabel('Training bin');
ylabel('Testing bin');
fixfig;grid off

% pairwise plots
hS=subplot(2,3,4);
cMap=colormap(hS,parula);
hB=colorbar;
set(gca,'clim',[min(vecStimTime) max(vecStimTime)]);
hold on

for intRec=1:intUseRecs
	%vecH(intRec)=plot([matRateZ(:,intRec); matRateZ(1,intRec)],[matPerfZ(:,intRec); matPerfZ(1,intRec)],'color',[0.5 0.5 0.5]);
end
h=cline([vecMeanRate],[vecMeanPerf],[vecStimTime(:)],[vecStimTime(:)]);
plot([min(vecMeanRate) max(vecMeanRate)],[dblChance dblChance],'--','color',[0.3 0.3 0.3]);
set(h,'LineWidth',2);
hold off
hB.Label.String = 'Time (s)';
xlabel('Mean population spiking rate (Hz)');
ylabel('Mean decoding performance');
title('Population dynamics of spiking rate&decoding');
fixfig;
%set(vecH,'LineWidth',1);


%get baseline
indBase = vecStimTime<0;
indPeak = vecStimTime>0 & vecStimTime<0.1;
vecBaseRate = mean(matRate(indBase,:),1);
matPeakRateZ = bsxfun(@rdivide,bsxfun(@minus,matRate(indPeak,:),vecBaseRate),std(matRate));
matPeakPerfZ = bsxfun(@rdivide,bsxfun(@minus,matPerf(indPeak,:),dblChance),std(matPerf));
%matPeakPerf = matPerfZ(indPeak,:);
vecPeakT = vecStimTime(indPeak);
intUsePoints = size(matPeakRateZ,1);
hS=subplot(2,3,5);
cMap=colormap(hS,parula(intUsePoints));
hB=colorbar;
set(gca,'clim',[min(vecPeakT) max(vecPeakT)]);
hold on
vecH = [];
h=cline(mean(matPeakRateZ,2),[mean(matPeakPerfZ,2)],[vecPeakT(:)],[vecPeakT(:)]);
vecRateP = nan(1,intUsePoints);
vecPerfP = nan(1,intUsePoints);
vecDiffP = nan(1,intUsePoints);
for intP=1:intUsePoints
	errorbar(mean(matPeakRateZ(intP,:),2),mean(matPeakPerfZ(intP,:),2),...
		std(matPeakPerfZ(intP,:),[],2)./sqrt(intUseRecs),std(matPeakPerfZ(intP,:),[],2)./sqrt(intUseRecs),...
		std(matPeakRateZ(intP,:),[],2)./sqrt(intUseRecs),std(matPeakRateZ(intP,:),[],2)./sqrt(intUseRecs),'color',cMap(intP,:));
	
	[h,pR]=ttest(matPeakRateZ(intP,:));
	[h,pP]=ttest(matPeakPerfZ(intP,:));
	[h,pD]=ttest(matPeakRateZ(intP,:),matPeakPerfZ(intP,:));
	vecRateP(intP) = pR;
	vecPerfP(intP) = pP;
	vecDiffP(intP) = pD;
end
hB.Label.String = 'Time (s)';
xlabel('Norm. pop. spiking rate (sd above base rate)');
ylabel('Norm. decoding perf. (sd above chance)');
title('Normalized rate/perf of onset (\mu +/- SEM)');
fixfig;

vecRateSigma = -norminv(vecRateP/2);
vecPerfSigma = -norminv(vecPerfP/2);
vecDiffSigma = -norminv(vecDiffP/2);

%real vals
%vecRateSigma = imnorm(mean(matRate(indPeak,:),2));
%vecPerfSigma = imnorm(mean(matPerf(indPeak,:),2));
%vecDiffSigma = imnorm(vecDiffSigma);

%matConfOverTime

subplot(2,3,6);
hold on
colororder(gca,{'k'});
plot([min(vecPeakT) max(vecPeakT)],-norminv(0.05/2)*[1 1],'--','color',[0.5 0.5 0.5]);
plot(vecPeakT,vecRateSigma,'color',lines(1));
plot(vecPeakT,vecPerfSigma,'color','k');
plot(vecPeakT,vecDiffSigma,'color',[0.8 0 0]);
ylabel('Significance (sigma)');
xlabel('Time (s)');
legend({'\alpha=0.05','\DeltaRate (vs base)','\DeltaDecoding (vs chance)','\DeltaRate vs \DeltaDecoding'},'location','best');
title('Spiking response precedes stimulus information');
vecSigmaTicks = get(gca,'ytick');
vecPvalTicks = (1-normcdf(vecSigmaTicks))*2;
fixfig;
yyaxis right
ylabel('Significance (p-value)');
ylim([min(vecSigmaTicks) max(vecSigmaTicks)]);
set(gca,'ytick',vecSigmaTicks,'yticklabel',cellfun(@sprintf,cellfill('%.1e',size(vecPvalTicks)),vec2cell(vecPvalTicks),'UniformOutput',false));
fixfig;

drawnow;
export_fig(fullpath(strFigurePath,sprintf('Q1_DynamicsSpikingDecoding.tif')));
export_fig(fullpath(strFigurePath,sprintf('Q1_DynamicsSpikingDecoding.pdf')));
	