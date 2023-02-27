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
	strFigurePathSR = 'F:\Drive\PopTimeCoding\single_recs';
	strFigurePath = 'F:\Drive\PopTimeCoding\figures\';
	strTargetDataPath = 'F:\Drive\PopTimeCoding\data\';
else
	strFigurePath = 'D:\Data\Results\PopTimeCoding\figures\';
	strTargetDataPath = 'D:\Data\Results\PopTimeCoding\data\';
end

%% load data
sFiles = dir ([strTargetDataPath 'Q1Data*.mat']);
strArea = 'V1';
intRecNum = numel(sFiles);
intBinNum = 240;
vecDecP = nan(1,intRecNum);
matDecPerf = nan(intBinNum,intRecNum);
matRate = nan(intBinNum,intRecNum);
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
	[h,dblP] = ttest(vecDecPerf(indUseForTest),dblChance);
	vecDecP(intFile) = dblP;
	matAcrossTimeDecoder(tril(true(size(matAcrossTimeDecoder))) & triu(true(size(matAcrossTimeDecoder)))) = vecDecPerf;
	
	%% plot
	figure;maxfig;
	subplot(2,3,1)
	plot(vecStimTime,vecSpikesPerBin./dblBinWidth)
	title('Spikes per bin')
	xlabel('Time after onset (s)');
	ylabel('Binned population spiking rate (Hz)');
	fixfig;
	
	subplot(2,3,2)
	plot(vecStimTime,vecDecPerf(:,end)');
	hold on
	plot([vecStimTime(1) vecStimTime(end)],[dblChance dblChance],'--','color',[0.5 0.5 0.5]);
	hold off
	title(sprintf('Dec perf %s; %s',strRec,strArea),'interpreter','none')
	xlabel('Time after onset (s)');
	ylabel('Fraction correct decoded');
	fixfig;
	
	subplot(2,3,3)
	plot(vecStimTime,(vecDecPerf(:,end)./dblChance)'./(vecSpikesPerBin./dblBinWidth))
	title('Dec perf / spike')
	xlabel('Time after onset (s)');
	ylabel('Performance/spike');
	fixfig;
	
	hS=subplot(2,3,4);
	cMap=colormap(hS,circcol);
	hB=colorbar;
	set(gca,'clim',[min(vecStimTime) max(vecStimTime)]);
	h=cline([vecSpikesPerBin(:); vecSpikesPerBin(1)]./dblBinWidth,[vecDecPerf(:,end); vecDecPerf(1,end)],[vecStimTime(:); vecStimTime(1)]);
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
	
	%% save data
	matDecPerf(:,intFile) = vecDecPerf;
	matRate(:,intFile) = (vecSpikesPerBin./dblBinWidth);
	matAcrossTimeDecoderAgg(:,:,intFile) = matAcrossTimeDecoder;
end

%% plot mean over recordings
vecCorrP = bonf_holm(vecDecP);
indUseRecs = vecCorrP<0.01;
matRateZ = matRate(:,indUseRecs);
matPerfZ = matDecPerf(:,indUseRecs);
matConfZ = matAcrossTimeDecoderAgg(:,:,indUseRecs);
vecMeanRateZ = mean(matRateZ,2);
vecMeanPerfZ = mean(matPerfZ,2);
matMeanConfZ = mean(matConfZ,3);
intUseRecs = size(matRateZ,2);

figure;maxfig;
%avg rate/decoding/confusion matrix


subplot(2,3,1)
plot(vecStimTime,vecMeanRateZ);
title(sprintf('Average spiking per %.1f ms bin, n=%d recs',dblBinWidth*1000,intUseRecs))
xlabel('Time after onset (s)');
ylabel('Binned population spiking rate (Hz)');
fixfig;

subplot(2,3,2)
plot(vecStimTime,vecMeanPerfZ,'color','k');
hold on
plot([vecStimTime(1) vecStimTime(end)],[dblChance dblChance],'--','color',[0.5 0.5 0.5]);
hold off
title(sprintf('Dec perf'),'interpreter','none')
xlabel('Time after onset (s)');
ylabel('Fraction correct decoded');
fixfig;


hS2=subplot(2,3,3);
colormap(hS2,parula);
imagesc(vecStimTime,vecStimTime,matMeanConfZ);
hB2=colorbar;
hB2.Label.String = 'Decoding performance';
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
h=cline([vecMeanRateZ],[vecMeanPerfZ],[vecStimTime(:)],[vecStimTime(:)]);
plot([min(vecMeanRateZ) max(vecMeanRateZ)],[dblChance dblChance],'--','color',[0.3 0.3 0.3]);
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
vecBaseRate = mean(matRateZ(indBase,:),1);
matPeakRate = bsxfun(@rdivide,bsxfun(@minus,matRateZ(indPeak,:),vecBaseRate),std(matRateZ));
matPeakPerf = bsxfun(@rdivide,bsxfun(@minus,matPerfZ(indPeak,:),dblChance),std(matPerfZ));
%matPeakPerf = matPerfZ(indPeak,:);
vecPeakT = vecStimTime(indPeak);
intUsePoints = size(matPeakRate,1);
hS=subplot(2,3,5);
cMap=colormap(hS,parula(intUsePoints));
hB=colorbar;
set(gca,'clim',[min(vecPeakT) max(vecPeakT)]);
hold on
vecH = [];
h=cline(mean(matPeakRate,2),[mean(matPeakPerf,2)],[vecPeakT(:)],[vecPeakT(:)]);
vecRateP = nan(1,intUsePoints);
vecPerfP = nan(1,intUsePoints);
vecDiffP = nan(1,intUsePoints);
for intP=1:intUsePoints
	errorbar(mean(matPeakRate(intP,:),2),mean(matPeakPerf(intP,:),2),...
		std(matPeakPerf(intP,:),[],2)./sqrt(intUseRecs),std(matPeakPerf(intP,:),[],2)./sqrt(intUseRecs),...
		std(matPeakRate(intP,:),[],2)./sqrt(intUseRecs),std(matPeakRate(intP,:),[],2)./sqrt(intUseRecs),'color',cMap(intP,:));
	
	[h,pR]=ttest(matPeakRate(intP,:));
	[h,pP]=ttest(matPeakPerf(intP,:));
	[h,pD]=ttest(matPeakRate(intP,:),matPeakPerf(intP,:));
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
	