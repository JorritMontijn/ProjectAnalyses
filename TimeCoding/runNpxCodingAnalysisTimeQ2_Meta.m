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

%% define parameters
%ISI distro
dblBinSizeISI = 1/1000;
vecBinEdgesISI = 0:dblBinSizeISI:0.05;
vecBinCentersISI = vecBinEdgesISI(2:end)-dblBinSizeISI/2;

%temp corr
dblMaxT = 1000*0.05;
dblStep = 1000*0.002;
vecBinD = 0:dblStep:dblMaxT;
vecBinD_c = vecBinD(2:end)-dblStep/2;

%quantiles
intNumQ = 19;
vecQuantiles = linspace(0,1,intNumQ+2);
vecQuantiles([1 end]) = [];

%mean/sd/cv
dblActBinW = 50;
vecActBins = 0:dblActBinW:750;
vecActBinsC = vecActBins(2:end)-dblActBinW/2;

%% load data
sFiles = dir ([strTargetDataPath 'Q2Data*.mat']);
strArea = 'V1';
intRecNum = numel(sFiles);

%% pre-allocate
matISINormCounts = nan(numel(vecBinCentersISI),intRecNum);
matISIExpPdf = nan(numel(vecBinCentersISI),intRecNum);
matISICounts = nan(numel(vecBinCentersISI),intRecNum);
matISICorrs = nan(numel(vecBinD_c),3,intRecNum);
matISICorrCounts = nan(numel(vecBinD_c),3,intRecNum);
matRatioQ = nan(intNumQ,intRecNum);
matRealQ = nan(intNumQ,intRecNum);
matExpQ = nan(intNumQ,intRecNum);

cellAllMperTrial = cell(1,intRecNum);
cellAllSperTrial = cell(1,intRecNum);
cellAllMperTrial_S = cell(1,intRecNum);
cellAllSperTrial_S = cell(1,intRecNum);
matCountsCV = nan(numel(vecActBinsC),intRecNum);
matMeansSd = nan(numel(vecActBinsC),intRecNum);
matMeansCV = nan(numel(vecActBinsC),intRecNum);
matCountsCV_S = nan(numel(vecActBinsC),intRecNum);
matMeansSd_S = nan(numel(vecActBinsC),intRecNum);
matMeansCV_S = nan(numel(vecActBinsC),intRecNum);
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
	dblMaxDur = 1-dblStartT;
	%{
	dblStartT
	strRec
	cellIFR_perTrial
	cellTimeIFR_perTrial
	cellISI_perTrial
	cellIFR_perTrial_S
	cellIFR_perTrial_SN
	cellIFR_perTrial_SS
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
	%}
	
	%% single plot 1
	figure;maxfig;
	subplot(2,3,1)
	dblLambda = 1./mean(vecAllISI);
	vecExpPdf = dblLambda.*exp(-dblLambda.*vecBinCentersISI);
	vecCountsISI = histcounts(vecAllISI,vecBinEdgesISI);
	hold on
	plot(vecBinCentersISI*1000,vecCountsISI./sum(vecCountsISI(:)))
	plot(vecBinCentersISI*1000,vecExpPdf./sum(vecExpPdf(:)))
	hold off
	set(gca,'yscale','log');
	xlabel('Inter-spike interval (ms)');
	ylabel('Normalized count (n)');
	legend({'Observed','Theory (Exponential)'});
	title('Population spiking dynamics');
	fixfig;
	
	
	subplot(2,3,2)
	
	%real
	vecD1 = 1000*vecAllISI(1:(end-1));
	vecD2 = 1000*vecAllISI(2:end);
	[r,p]=corr(vecD1,vecD2);
	[vecCounts_D,vecMeans_D,vecSDs_D,cellVals_D,cellIDs_D] = makeBins(vecD1,vecD2,vecBinD);
	
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
	errorbar(vecBinD_c,vecMeans_D,vecSDs_D./sqrt(vecCounts_D));
	hold on
	errorbar(vecBinD_c,vecMeans_R-0.02,vecSDs_R./sqrt(vecCounts_R));
	errorbar(vecBinD_c,vecMeans_S+0.02,vecSDs_S./sqrt(vecCounts_S));
	
	xlabel('ISI spikes i,i+1 (ms)');
	ylabel('ISI spikes i+1,i+2 (ms)');
	title(sprintf('ISI temporal correlation, r=%.3f, p=%.3f',r,p));
	fixfig;
	
	legend({'Observed','Theory (Exponential)','Shuffled'},'location','best');
	
	% single plot 2
	%plot distribution of firing rates vs shuffled
	
	intTrials = numel(cellISI_perTrial);
	matISIQ_Real = nan(intNumQ,intTrials);
	matISIQ_Exp = nan(intNumQ,intTrials);
	%get quantiles for each trial
	cellGenISI = cell(1,intTrials);
	for i=1:intTrials
		vecRealISI = cellISI_perTrial{i};
		if numel(vecRealISI) < round(numel(vecQuantiles)+1)
			continue;
		end
		
		dblLambda = numel(vecRealISI)/dblMaxDur;
		
		vecGenISI = exprnd(1/dblLambda,[1 round(dblMaxDur*dblLambda*2)]);
		vecT = cumsum(vecGenISI);
		while max(vecT) < dblMaxDur
			vecGenISI = cat(2,vecGenISI,exprnd(1/dblLambda,[1 round(dblMaxDur*dblLambda)]));
			vecT = cumsum(vecGenISI);
		end
		vecGenISI = vecGenISI(vecT<dblMaxDur);
		cellGenISI{i} = vecGenISI;
		
		%sort ISIs
		vecSortRealISI = sort(vecRealISI);
		vecSortGenISI = sort(vecGenISI);
		if numel(vecSortGenISI) < round(numel(vecQuantiles)+1)
			continue;
		end
		
		matISIQ_Real(:,i) = vecSortRealISI(round(vecQuantiles*numel(vecSortRealISI)));
		matISIQ_Exp(:,i) = vecSortGenISI(round(vecQuantiles*numel(vecSortGenISI)));
		
	end
	vecAllGenISI = cell2vec(cellGenISI);
	
	%take lowest and highest 5%
	vecMeanQ_Real = nanmean(matISIQ_Real,2);
	vecMeanQ_Exp = nanmean(matISIQ_Exp,2);
	
	subplot(2,3,3)
	plot(vecQuantiles,vecMeanQ_Real./vecMeanQ_Exp);
	hold on
	plot(vecQuantiles([1 end]),[1 1],'--','color',[0.5 0.5 0.5]);
	hold off
	ylabel('ISI Ratio Real/Exponential');
	xlabel('Quantile of ISIs');
	title('Real dynamics show highly variable rates (avg over trials)');
	legend({'Real ISI distribution','Exponential process'},'location','best');
	fixfig;
	
	subplot(2,3,4)
	plot(vecQuantiles,1000*vecMeanQ_Real);
	hold on
	plot(vecQuantiles,1000*vecMeanQ_Exp);
	hold off
	ylabel('ISI (ms)');
	xlabel('Quantile of ISIs');
	title('Source distros');
	legend({'Real ISI distribution','Exponential process'},'location','best');
	fixfig;
	
	%mean and sd over time (IFR)
	subplot(2,3,5)
	[vecCountsSd,vecMeansSd,vecSDsSd] = makeBins(vecMperTrial(:)',vecSperTrial(:)',vecActBins);
	indPlotBins = vecCountsSd>10;
	
	[vecCountsSd_S,vecMeansSd_S,vecSDsSd_S] = makeBins(vecMperTrial_S,vecSperTrial_S,vecActBins);
	indPlotBins_S = vecCountsSd_S>10;
	
	hold on
	scatter(vecMperTrial_S,vecSperTrial_S(:),[],[0.5 0.5 0.5],'.');
	errorbar(vecActBinsC(indPlotBins_S),vecMeansSd_S(indPlotBins_S),vecSDsSd_S(indPlotBins_S)./sqrt(vecCountsSd_S(indPlotBins_S)),'color',[0.5 0.5 0.5])
	scatter(vecMperTrial,vecSperTrial(:),[],[0.5 0.5 1],'.');
	errorbar(vecActBinsC(indPlotBins),vecMeansSd(indPlotBins),vecSDsSd(indPlotBins)./sqrt(vecCountsSd(indPlotBins)),'color',[0 0 1])
	xlabel('Mean population activity (Hz)')
	ylabel('Sd of population activity (Hz)')
	fixfig;
	
	
	subplot(2,3,6)
	vecCVperTrial = vecSperTrial./vecMperTrial;
	[vecCountsCV,vecMeansCV,vecSDsCV] = makeBins(vecMperTrial,vecCVperTrial,vecActBins);
	vecCVperTrial_S = vecSperTrial./vecMperTrial;
	[vecCountsCV_S,vecMeansCV_S,vecSDsCV_S] = makeBins(vecMperTrial_S,vecCVperTrial_S,vecActBins);
	hold on
	scatter(vecMperTrial_S,vecCVperTrial_S(:),[],[0.5 0.5 0.5],'.');
	errorbar(vecActBinsC(indPlotBins_S),vecMeansCV_S(indPlotBins_S),vecSDsCV_S(indPlotBins_S)./sqrt(vecCountsCV_S(indPlotBins_S)),'color',[0.5 0.5 0.5])
	scatter(vecMperTrial,vecCVperTrial(:),[],[0.5 0.5 1],'.');
	errorbar(vecActBinsC(indPlotBins),vecMeansCV(indPlotBins),vecSDsCV(indPlotBins)./sqrt(vecCountsCV(indPlotBins)),'color',[0 0 1])
	xlabel('Mean population activity (Hz)')
	ylabel('CV of population activity (Hz)')
	fixfig;
	
	% 	%sparseness vs mean rate
	% 	r1=corr(vecLperTrial(:),vecPopSparseness(:));
	% 	r2=corr(vecMperTrial(:),vecPopSparseness(:));
	% 	[rSpH,pSpH]=nancorr(vecMperTrial(:),vecPopSparseness(:));
	% 	dblActBinW = 50;
	% 	vecActBinsH = 0:dblActBinW:1700;
	% 	vecActBinsHC = vecActBinsH(2:end)-dblActBinW/2;
	% 	[vecCounts2,vecMeans2,vecSDs2] = makeBins(vecMperTrial,vecPopSparseness,vecActBinsH);
	% 	indPlotBins2 = vecCounts2>10;
	% 	mdl = fitlm(vecMperTrial,vecPopSparseness);
	% 	ci = coefCI(mdl);
	% 	[ypred,yci] = predict(mdl,vecActBinsHC(indPlotBins2)');
	% 	subplot(2,3,5)
	% 	hold on;
	% 	scatter(vecMperTrial,vecPopSparseness,[],1-(1-lines(1))*(2/3),'.');
	% 	%errorbar(vecActBinsHC(indPlotBins2),vecMeans2(indPlotBins2),vecSDs2(indPlotBins2)./sqrt(vecCounts2(indPlotBins2)),'color',lines(1));
	% 	plot(vecActBinsHC(indPlotBins2),yci(:,1),'--','color',lines(1));
	% 	plot(vecActBinsHC(indPlotBins2),ypred,'color',lines(1));
	% 	plot(vecActBinsHC(indPlotBins2),yci(:,2),'--','color',lines(1));
	%
	%
	% 	r1=corr(vecLperTrial_SN(:),vecPopSparsenessS(:));
	% 	r2=corr(vecMperTrial_SN(:),vecPopSparsenessS(:));
	% 	[rSpHS,pSpHS]=nancorr(vecMperTrial_SN(:),vecPopSparsenessS(:));
	% 	vecActBinsH = 0:dblActBinW:1700;
	% 	vecActBinsHC = vecActBinsH(2:end)-dblActBinW/2;
	% 	[vecCounts2,vecMeans2,vecSDs2] = makeBins(vecMperTrial_SN,vecPopSparsenessS,vecActBinsH);
	% 	indPlotBins2 = vecCounts2>10;
	% 	mdl = fitlm(vecMperTrial_SN,vecPopSparsenessS);
	% 	ci = coefCI(mdl);
	% 	[ypred,yci] = predict(mdl,vecActBinsHC(indPlotBins2)');
	%
	% 	hold on;
	% 	scatter(vecMperTrial_SN,vecPopSparsenessS,[],[0.7 0.7 0.7],'.');
	% 	%errorbar(vecActBinsHC(indPlotBins2),vecMeans2(indPlotBins2),vecSDs2(indPlotBins2)./sqrt(vecCounts2(indPlotBins2)),'color',lines(1));
	% 	plot(vecActBinsHC(indPlotBins2),yci(:,1),'--','color',[0.5 0.5 0.5]);
	% 	plot(vecActBinsHC(indPlotBins2),ypred,'color',[0.5 0.5 0.5]);
	% 	plot(vecActBinsHC(indPlotBins2),yci(:,2),'--','color',[0.5 0.5 0.5]);
	%
	% 	hold off;
	% 	ylabel('Population sparseness per trial');
	% 	xlabel('Mean pop. firing rate per trial (Hz) ');
	% 	title(sprintf('Corr(sparse, rate); r=%.3f, p=%.1e; shuff:r=%.3f, p=%.1e',rSpH,pSpH,rSpHS,pSpHS));
	% 	fixfig;
	
	
	%hold on
	%plot(vecQuantiles,vecMeanQ_Exp);
	%hold off
	% 	%transform to quantiles
	% 	vecSortRealISI = sort(vecAllISI);
	% 	vecSortGenISI = sort(vecAllGenISI);
	% 	vecIdxRealQ = round(linspace(1,numel(vecSortRealISI),intNumQ+2));
	% 	vecIdxGenQ = round(linspace(1,numel(vecSortGenISI),intNumQ+2));
	% 	vecIdxRealQ([1 end]) = [];
	% 	vecIdxGenQ([1 end]) = [];
	% 	vecRealQ = vecSortRealISI(vecIdxRealQ);
	% 	vecGenQ = vecSortGenISI(vecIdxGenQ);
	% 	%%{
	% 	subplot(2,3,3)
	% 	scatter(vecGenQ,vecRealQ)
	% 	set(gca,'xscale','log','yscale','log');
	% 	vecLimX = get(gca,'xlim');
	% 	vecLimY = get(gca,'ylim');
	% 	vecLim = [min([vecLimX vecLimY]) max([vecLimX vecLimY])];
	% 	%vecLim = [3.3*1e-5 max([vecLimX vecLimY])];
	% 	hold on
	% 	plot(vecLim,vecLim,'k');
	% 	hold off
	% 	xlim(vecLim);
	% 	ylim(vecLim);
	% 	xlabel('Exponential ISIs (s)');
	% 	ylabel('Real ISIs (s)');
	% 	fixfig;
	% 	%%}
	
	%% save data
	vecLogExpPdf = log10(vecExpPdf./sum(vecExpPdf(:)));
	vecLogCounts = log10(vecCountsISI./sum(vecCountsISI(:)));
	vecNormCounts = (vecLogCounts - min(vecLogExpPdf))./range(vecLogExpPdf);
	
	%subplot(2,3,4)
	%plot(vecBinCentersISI*1000,vecNormCounts);
	%hold on
	%plot(vecBinCentersISI*1000,linspace(1,0,numel(vecBinCentersISI*1000)));
	%hold off;
	
	matISINormCounts(:,intFile) = vecNormCounts;
	matISIExpPdf(:,intFile) = vecExpPdf;
	matISICounts(:,intFile) = vecCountsISI;
	matISICorrs(:,:,intFile) = cat(2,vecMeans_D,vecMeans_R,vecMeans_S);
	matISICorrCounts(:,:,intFile) = cat(2,vecCounts_D,vecCounts_R,vecCounts_S);
	matRatioQ(:,intFile) = vecMeanQ_Real./vecMeanQ_Exp;
	matRealQ(:,intFile) = vecMeanQ_Real;
	matExpQ(:,intFile) = vecMeanQ_Exp;
	%vecCorrSparseness(intFile) = rSpH;
	%vecCorrSparsenessS(intFile) = rSpHS;
	
	matCountsCV(:,intFile) = vecCountsSd;
	matMeansSd(:,intFile) = vecMeansSd;
	matMeansCV(:,intFile) = vecMeansCV;
	matCountsCV_S(:,intFile) = vecCountsSd_S;
	matMeansSd_S(:,intFile) = vecMeansSd_S;
	matMeansCV_S(:,intFile) = vecMeansCV_S;
	
	cellAllMperTrial{intFile} = vecMperTrial;
	cellAllSperTrial{intFile} = vecSperTrial;
	cellAllMperTrial_S{intFile} = vecMperTrial_S;
	cellAllSperTrial_S{intFile} = vecSperTrial_S;
end

%% plot mean over recordings
%{
subplot(2,3,4)
matC = lines(2);
%plot(vecBinCentersISI*1000,bsxfun(@rdivide,matISINormCounts,mean(matISINormCounts)),'color',matC(1,:));
plot(vecBinCentersISI*1000,matISINormCounts,'color',matC(1,:));
hold on
plot(vecBinCentersISI([1 end])*1000,[1 0],'color',matC(2,:));
hold off
%}
figure;maxfig;
subplot(2,3,1)
matMeanISIcorrs = mean(matISICorrs,3);
matSdISIcorrs = std(matISICorrs,[],3);
hold on
for intVar=1:3
	errorbar(vecBinD_c,matMeanISIcorrs(:,intVar),matSdISIcorrs(:,intVar)./sqrt(intRecNum));
end
hold off
xlabel('ISI spikes i,i+1 (ms)');
ylabel('ISI spikes i+1,i+2 (ms)');
title(sprintf('Mean +/- SEM over recs'));
fixfig;

legend({'Observed','Theory (Exponential)','Shuffled'},'location','best');

fixfig;

subplot(2,3,2)
errorbar(vecQuantiles,mean(matRealQ*1000,2),std(matRealQ*1000,[],2)./sqrt(intRecNum));
hold on
errorbar(vecQuantiles,mean(matExpQ*1000,2),std(matExpQ*1000,[],2)./sqrt(intRecNum));
hold off
ylabel('ISI (ms)');
xlabel('Quantile of ISIs');
title('Mean +/- SEM over recs');
legend({'Real ISI distribution','Exponential process'},'location','best');
fixfig;

subplot(2,3,3)
errorbar(vecQuantiles,mean(matRatioQ,2),std(matRatioQ,[],2)./sqrt(intRecNum));
hold on
plot(vecQuantiles([1 end]),[1 1],'--','color',[0.5 0.5 0.5]);
hold off
ylabel('ISI Ratio Real/Exponential');
xlabel('Quantile of ISIs');
title('Mean +/- SEM over recs');
legend({'Real ISI distribution','Exponential process'},'location','best');
fixfig;
%
% subplot(2,3,4);
% vecBins = -1:0.1:1;
% vecSparseCounts = histcounts(vecCorrSparseness,vecBins);
% vecSparseCountsS = histcounts(vecCorrSparsenessS,vecBins);
% vecPlotBins = vecBins(2:end) - median(diff(vecBins))/2;
% plot(vecPlotBins,vecSparseCounts);
% hold on
% plot(vecPlotBins,vecSparseCountsS,'color',[0.5 0.5 0.5]);
% hold off
% fixfig;

%generate data for exponential processes
dblLambda = 1;
dblDuration = 1;
vecPopSizes = 50:50:max(vecActBins);
intPopNum = numel(vecPopSizes);
intTrials = numel(vecRperTrial);
matMperTrial_G = nan(intTrials,intPopNum);
matSperTrial_G = nan(intTrials,intPopNum);
for intPopSizeIdx = 1:intPopNum
	intPopSize = vecPopSizes(intPopSizeIdx)
	
	matSpikeCounts = nan(intPopSize,intTrials);
	for intTrial=1:intTrials
		cellSpikeT = cell(1,intPopSize);
		for intNeuron=1:intPopSize
		
			vecISI = exprnd(dblLambda,[1 100]);
			vecSpikeT = cumsum(vecISI(2:end))-vecISI(1);
			cellSpikeT{intNeuron} = vecSpikeT(vecSpikeT>0 & vecSpikeT<dblDuration);
		end
		
		%combine spikes from whole population & calc IFR
		vecSpikeT = sort(cell2vec(cellSpikeT));
		[vecTimeIFRG,vecIFRG] = getIFR(vecSpikeT,0,dblDuration,[],[],[],0);
		matMperTrial_G(intTrial,intPopSizeIdx) = mean(vecIFRG);
		matSperTrial_G(intTrial,intPopSizeIdx) = std(vecIFRG);
	end
end


%% mean and sd scale together: CV is constant
subplot(2,3,4);cla;
%theory
vecMperTrial_G = matMperTrial_G(:);
vecSperTrial_G = matSperTrial_G(:);
vecBeta_G = polyfit(vecMperTrial_G,vecSperTrial_G,1);
dblSlope_G = ((vecMperTrial_G' * vecMperTrial_G) \ vecMperTrial_G') * vecSperTrial_G;
vecCVperTrial_G = vecSperTrial_G./vecMperTrial_G;
[vecCountsCV_G,vecMeansCV_G,vecSDsCV_G] = makeBins(vecMperTrial_G,vecCVperTrial_G,vecActBins);
[vecCountsSd_G,vecMeansSd_G,vecSDsSd_G] = makeBins(vecMperTrial_G,vecSperTrial_G,vecActBins);
indPlotBins_G = vecCountsSd_G>10;

%shuffled
matMeansSdRem_S = matMeansSd_S;
matMeansSdRem_S(matCountsCV_S<10)=nan;
matMeansCVRem_S = matMeansCV_S;
matMeansCVRem_S(matCountsCV_S<10)=nan;
vecX_S = cell2vec(cellAllMperTrial_S);
vecY_S = cell2vec(cellAllSperTrial_S);
dblSlope_S = ((vecX_S' * vecX_S) \ vecX_S') * vecY_S;

matC = lines(3);
%real
matMeansSdRem = matMeansSd;
matMeansSdRem(matCountsCV<10)=nan;
matMeansCVRem = matMeansCV;
matMeansCVRem(matCountsCV<10)=nan;
vecX = cell2vec(cellAllMperTrial);
vecY = cell2vec(cellAllSperTrial);
dblSlope = ((vecX' * vecX) \ vecX') * vecY;


hold on
plot(vecActBinsC,dblSlope*vecActBinsC,'--','color',matC(1,:));
plot(vecActBinsC,dblSlope_S*vecActBinsC,'--','color',matC(3,:));
plot(vecActBinsC,dblSlope_G*vecActBinsC,'--','color',matC(2,:));
errorbar(vecActBinsC,nanmean(matMeansSdRem,2),nanstd(matMeansSdRem,[],2)./sqrt(intRecNum),'linestyle','none','color',matC(1,:));
errorbar(vecActBinsC,nanmean(matMeansSdRem_S,2),nanstd(matMeansSdRem_S,[],2)./sqrt(intRecNum),'linestyle','none','color',matC(3,:));
errorbar(vecActBinsC,nanmean(vecMeansSd_G,2),nanstd(vecSDsSd_G,[],2)./sqrt(intRecNum),'linestyle','none','color',matC(2,:));
hold off
legend({'Neural data','Shuffled','Theory (Exponential Process)'},'location','best')
xlabel('Mean of spiking rate (Hz)');
ylabel('Sd of spiking rate (Hz)');
title(sprintf('Constant coefficient of variation (CV) = %.3f',dblSlope));
fixfig;

h=subplot(2,3,5);delete(h);subplot(2,3,5);
vecDeviationM = (nanmean(matMeansCVRem,2) ./ dblSlope)-1;
vecDeviationSEM = nanstd(matMeansCVRem,[],2) ./ dblSlope;
vecDeviationM_S = (nanmean(matMeansCVRem_S,2) ./ dblSlope_S)-1;
vecDeviationSEM_S = nanstd(matMeansCVRem_S,[],2) ./ dblSlope_S;
vecDeviationM_G = (nanmean(vecMeansCV_G,2) ./ dblSlope_G)-1;
vecDeviationSEM_G = nanstd(matMeansCVRem,[],2) ./ dblSlope_G;


hold on
plot(vecActBinsC,100*vecDeviationM,'color',matC(1,:));
plot(vecActBinsC,100*vecDeviationM_S,'color',matC(3,:));
plot(vecActBinsC,100*vecDeviationM_G,'color',matC(2,:));
hold off
%ylim([0 1]);
xlabel('Mean of spiking rate (Hz)');
ylabel('Deviation from constant CV (%)');
title(sprintf('How best to plot this?'));
fixfig;

%% plot mean vs euclidian


drawnow;
export_fig(fullpath(strFigurePath,sprintf('Q2_ConstantCV.tif')));
export_fig(fullpath(strFigurePath,sprintf('Q2_ConstantCV.pdf')));
	