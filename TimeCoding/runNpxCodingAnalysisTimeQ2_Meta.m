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
dblRemOnset = 0; %remove onset period in seconds; 0.125 for sim, 0.25 for npx
%% onset string
if dblRemOnset == 0
	strOnset = '';
else
	strOnset = sprintf('%.2f',dblRemOnset);
end
sFiles = dir ([strTargetDataPath 'Q2Data*' strOnset '.mat']);
intRecNum = numel(sFiles);

%% pre-allocate
matISINormCounts = nan(numel(vecBinCentersISI),intRecNum);
matISIExpPdf = nan(numel(vecBinCentersISI),intRecNum);
matISICounts = nan(numel(vecBinCentersISI),intRecNum);
matISICorrs = nan(numel(vecBinD_c),4,intRecNum);
matISICorrCounts = nan(numel(vecBinD_c),4,intRecNum);
matRatioQ = nan(intNumQ,intRecNum);
matRealQ = nan(intNumQ,intRecNum);
matExpQ = nan(intNumQ,intRecNum);

cellAllMperTrial = cell(1,intRecNum);
cellAllSperTrial = cell(1,intRecNum);
matCountsCV = nan(numel(vecActBinsC),intRecNum);
matMeansSd = nan(numel(vecActBinsC),intRecNum);
matMeansCV = nan(numel(vecActBinsC),intRecNum);
cellTypes =  {'Real','Poiss','ShuffTid','Shuff','PoissGain','Uniform'};
for intFile=1:intRecNum
	%% load
	load(fullpath(sFiles(intFile).folder,sFiles(intFile).name));
	strRec = getFlankedBy(sFiles(intFile).name,'Q2Data','_g0_t0');
	
	%% unpack
	%real
	intReal = 1;
	strType = sAggData(intReal).strType;
	cellISI_perTrial = sAggData(intReal).cellISI_perTrial;
	cellIFR_perTrial = sAggData(intReal).cellIFR_perTrial;
	
	matResp = sAggData(intReal).matResp;
	vecSperTrial = sAggData(intReal).vecSperTrial;
	vecHperTrial = sAggData(intReal).vecHperTrial;
	vecLperTrial = sAggData(intReal).vecLperTrial;
	vecMperTrial = sAggData(intReal).vecMperTrial;
	
	vecPopSparseness = sAggData(intReal).vecPopSparseness;
	
	%shufftid
	intS = 3;
	strType_S = sAggData(intS).strType;
	cellISI_perTrial_S = sAggData(intS).cellISI_perTrial;
	vecSperTrial_S = sAggData(intS).vecSperTrial;
	vecHperTrial_S = sAggData(intS).vecHperTrial;
	vecLperTrial_S = sAggData(intS).vecLperTrial;
	vecMperTrial_S = sAggData(intS).vecMperTrial;
	vecPopSparseness_S = sAggData(intS).vecPopSparseness;
	
	%poisson
	intP = 2;
	strType_P = sAggData(intP).strType;
	cellISI_perTrial_P = sAggData(intP).cellISI_perTrial;
	vecSperTrial_P = sAggData(intP).vecSperTrial;
	vecHperTrial_P = sAggData(intP).vecHperTrial;
	vecLperTrial_P = sAggData(intP).vecLperTrial;
	vecMperTrial_P = sAggData(intP).vecMperTrial;
	vecPopSparseness_P = sAggData(intP).vecPopSparseness;
	
	%uniform
	intU = 6;
	strType_U = sAggData(intU).strType;
	cellISI_perTrial_U = sAggData(intU).cellISI_perTrial;
	vecSperTrial_U = sAggData(intU).vecSperTrial;
	vecHperTrial_U = sAggData(intU).vecHperTrial;
	vecLperTrial_U = sAggData(intU).vecLperTrial;
	vecMperTrial_U = sAggData(intU).vecMperTrial;
	vecPopSparseness_U = sAggData(intU).vecPopSparseness;
	
	%% create derived variables
	%get data
	vecAllISI = cell2vec(cellISI_perTrial);
	vecAllISI_P = cell2vec(cellISI_perTrial_P);
	vecAllISI_S = cell2vec(cellISI_perTrial_S);
	vecAllISI_U = cell2vec(cellISI_perTrial_U);
	cellIFR_Reduced = cellfun(@(x) x(2:(end-2)),cellIFR_perTrial,'UniformOutput',false);
	vecAllIFR = cell2vec(cellIFR_Reduced);
	dblMaxDur = 1;
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
	
	%% run calculations
	dblLambda = 1./mean(vecAllISI);
	vecExpPdf = dblLambda.*exp(-dblLambda.*vecBinCentersISI);
	vecCountsISI = histcounts(vecAllISI,vecBinEdgesISI);
	vecCountsISI_U = histcounts(vecAllISI_U,vecBinEdgesISI);
	vecCountsISI_S = histcounts(vecAllISI_S,vecBinEdgesISI);
	vecCountsISI_P = histcounts(vecAllISI_P,vecBinEdgesISI);
	
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
	
	%trial-shuffled
	vecD1_U = 1000*vecAllISI_U(1:(end-1));
	vecD2_U = 1000*vecAllISI_U(2:end);
	[rTS,pTS]=corr(vecD1_U,vecD2_U);
	[vecCounts_U,vecMeans_U,vecSDs_U,cellVals_U,cellIDs_U] = makeBins(vecD1_U,vecD2_U,vecBinD);
	
	%ISI-shuffled
	vecD1_S = 1000*vecAllISI_S(1:(end-1));
	vecD2_S = 1000*vecAllISI_S(2:end);
	[rS,pS]=corr(vecD1,vecD2);
	[vecCounts_S,vecMeans_S,vecSDs_S,cellVals_S,cellIDs_S] = makeBins(vecD1_S,vecD2_S,vecBinD);
	
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
	
	%make bins
	[vecCountsSd,vecMeansSd,vecSDsSd] = makeBins(vecMperTrial(:)',vecSperTrial(:)',vecActBins);
	indPlotBins = vecCountsSd>10;
	
	[vecCountsSd_S,vecMeansSd_S,vecSDsSd_S] = makeBins(vecMperTrial_S,vecSperTrial_S,vecActBins);
	indPlotBins_S = vecCountsSd_S>10;
	
	[vecCountsSd_U,vecMeansSd_U,vecSDsSd_U] = makeBins(vecMperTrial_U,vecSperTrial_U,vecActBins);
	indPlotBins_U = vecCountsSd_U>10;
	
	[vecCountsSd_P,vecMeansSd_P,vecSDsSd_P] = makeBins(vecMperTrial_P,vecSperTrial_P,vecActBins);
	indPlotBins_P = vecCountsSd_P>10;
	
	%CVs
	vecCVperTrial = vecSperTrial./vecMperTrial;
	[vecCountsCV,vecMeansCV,vecSDsCV] = makeBins(vecMperTrial,vecCVperTrial,vecActBins);
	vecCVperTrial_S = vecSperTrial_S./vecMperTrial_S;
	[vecCountsCV_S,vecMeansCV_S,vecSDsCV_S] = makeBins(vecMperTrial_S,vecCVperTrial_S,vecActBins);
	vecCVperTrial_U = vecSperTrial_U./vecMperTrial_U;
	[vecCountsCV_U,vecMeansCV_U,vecSDsCV_U] = makeBins(vecMperTrial_U,vecCVperTrial_U,vecActBins);
	vecCVperTrial_P = vecSperTrial_P./vecMperTrial_P;
	[vecCountsCV_P,vecMeansCV_P,vecSDsCV_P] = makeBins(vecMperTrial_P,vecCVperTrial_P,vecActBins);
	
	%fit per recording
	intK = 1;
	dblSlope = ((vecMperTrial' * vecMperTrial) \ vecMperTrial') * vecSperTrial;
	vecFitY = vecMperTrial*dblSlope;
	[dblR2,dblSS_tot,dblSS_res,dblT,dblP,dblR2_adjusted,dblR2_SE] = getR2(vecSperTrial,vecFitY,intK);
	
	dblSlope_S = ((vecMperTrial_S' * vecMperTrial_S) \ vecMperTrial_S') * vecSperTrial_S;
	vecFitY_S = vecMperTrial_S*dblSlope_S;
	[dblR2_S,dblSS_tot_S,dblSS_res_S,dblT_S,dblP_S,dblR2_adjusted_S,dblR2_SE_S] = getR2(vecSperTrial_S,vecFitY_S,intK);
	
	dblSlope_U = ((vecMperTrial_U' * vecMperTrial_U) \ vecMperTrial_U') * vecSperTrial_U;
	vecFitY_U = vecMperTrial_U*dblSlope_U;
	[dblR2_U,dblSS_tot_U,dblSS_res_U,dblT_U,dblP_U,dblR2_adjusted_U,dblR2_SE_U] = getR2(vecSperTrial_U,vecFitY_U,intK);
	
	dblSlope_P = ((vecMperTrial_P' * vecMperTrial_P) \ vecMperTrial_P') * vecSperTrial_P;
	vecFitY_P = vecMperTrial_P*dblSlope_P;
	[dblR2_P,dblSS_tot_P,dblSS_res_P,dblT_P,dblP_P,dblR2_adjusted_P,dblR2_SE_P] = getR2(vecSperTrial_P,vecFitY_P,intK);
	
	%% save data
	vecLogExpPdf = log10(vecExpPdf./sum(vecExpPdf(:)));
	vecLogCounts = log10(vecCountsISI./sum(vecCountsISI(:)));
	vecNormCounts = (vecLogCounts - min(vecLogExpPdf))./range(vecLogExpPdf);
	
	matSlopes(1,intFile) = dblSlope;
	matSlopes(2,intFile) = dblSlope_S;
	matSlopes(3,intFile) = dblSlope_U;
	matSlopes(4,intFile) = dblSlope_P;
	
	matR2(1,intFile) = dblR2;
	matR2(2,intFile) = dblR2_S;
	matR2(3,intFile) = dblR2_U;
	matR2(4,intFile) = dblR2_P;
	
	matISINormCounts(:,intFile) = vecNormCounts;
	matISIExpPdf(:,intFile) = vecExpPdf;
	matISICounts(:,intFile) = vecCountsISI;
	matISICorrs(:,:,intFile) = cat(2,vecMeans_D,vecMeans_R,vecMeans_U,vecMeans_S);
	matISICorrCounts(:,:,intFile) = cat(2,vecCounts_D,vecCounts_R,vecCounts_U,vecCounts_S);
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
	
	matCountsCV_U(:,intFile) = vecCountsSd_U;
	matMeansSd_U(:,intFile) = vecMeansSd_U;
	matMeansCV_U(:,intFile) = vecMeansCV_U;
	
	matCountsCV_P(:,intFile) = vecCountsSd_P;
	matMeansSd_P(:,intFile) = vecMeansSd_P;
	matMeansCV_P(:,intFile) = vecMeansCV_P;
	
	cellAllMperTrial{intFile} = vecMperTrial;
	cellAllSperTrial{intFile} = vecSperTrial;
	cellAllMperTrial_S{intFile} = vecMperTrial_S;
	cellAllSperTrial_S{intFile} = vecSperTrial_S;
	cellAllMperTrial_U{intFile} = vecMperTrial_U;
	cellAllSperTrial_U{intFile} = vecSperTrial_U;
	cellAllMperTrial_P{intFile} = vecMperTrial_P;
	cellAllSperTrial_P{intFile} = vecSperTrial_P;
	
	if boolSingleRecPlots
		%% single plot 1
		figure;maxfig;
		subplot(2,3,1)
		hold on
		plot(vecBinCentersISI*1000,vecCountsISI./sum(vecCountsISI(:)))
		plot(vecBinCentersISI*1000,vecCountsISI_S./sum(vecCountsISI_S(:)))
		plot(vecBinCentersISI*1000,vecCountsISI_U./sum(vecCountsISI_U(:)))
		%plot(vecBinCentersISI*1000,vecExpPdf./sum(vecExpPdf(:)))
		plot(vecBinCentersISI*1000,vecCountsISI_P./sum(vecCountsISI_P(:)))
		
		
		hold off
		set(gca,'xscale','log');
		xlabel('Inter-spike interval (ms)');
		ylabel('Normalized count (n)');
		legend({strType,strType_S,strType_U,strType_P},'location','best');
		title('Population spiking dynamics');
		
		
		subplot(2,3,2)
		%plot
		errorbar(vecBinD_c,vecMeans_D,vecSDs_D./sqrt(vecCounts_D));
		hold on
		errorbar(vecBinD_c,vecMeans_U+0.02,vecSDs_U./sqrt(vecCounts_U));
		errorbar(vecBinD_c,vecMeans_S+0.02,vecSDs_S./sqrt(vecCounts_S));
		errorbar(vecBinD_c,vecMeans_R-0.02,vecSDs_R./sqrt(vecCounts_R));
		
		xlabel('ISI spikes i,i+1 (ms)');
		ylabel('ISI spikes i+1,i+2 (ms)');
		title(sprintf('ISI temporal correlation, r=%.3f, p=%.3f',r,p));
		legend({strType,strType_S,strType_U,strType_P},'location','best');
		
		% single plot 2
		%plot distribution of firing rates vs shuffled
		
		subplot(2,3,3)
		plot(vecQuantiles,vecMeanQ_Real./vecMeanQ_Exp);
		hold on
		plot(vecQuantiles([1 end]),[1 1],'--','color',[0.5 0.5 0.5]);
		hold off
		ylabel('ISI Ratio Real/Exponential');
		xlabel('Quantile of ISIs');
		title('Real dynamics show highly variable rates (avg over trials)');
		legend({'Real ISI distribution','Exponential process'},'location','best');
		
		subplot(2,3,4)
		plot(vecQuantiles,1000*vecMeanQ_Real);
		hold on
		plot(vecQuantiles,1000*vecMeanQ_Exp);
		hold off
		ylabel('ISI (ms)');
		xlabel('Quantile of ISIs');
		title('Source distros');
		legend({'Real ISI distribution','Exponential process'},'location','best');
		
		%mean and sd over time (IFR)
		subplot(2,3,5)
		
		hold on
		scatter(vecMperTrial_S,vecSperTrial_S(:),[],[0.5 0.5 0.5],'.');
		errorbar(vecActBinsC(indPlotBins_S),vecMeansSd_S(indPlotBins_S),vecSDsSd_S(indPlotBins_S)./sqrt(vecCountsSd_S(indPlotBins_S)),'color',[0.5 0.5 0.5])
		scatter(vecMperTrial,vecSperTrial(:),[],[0.5 0.5 1],'.');
		errorbar(vecActBinsC(indPlotBins),vecMeansSd(indPlotBins),vecSDsSd(indPlotBins)./sqrt(vecCountsSd(indPlotBins)),'color',[0 0 1])
		xlabel('Mean population activity (Hz)')
		ylabel('Sd of population activity (Hz)')
		
		
		subplot(2,3,6)
		hold on
		scatter(vecMperTrial_S,vecCVperTrial_S(:),[],[0.5 0.5 0.5],'.');
		errorbar(vecActBinsC(indPlotBins_S),vecMeansCV_S(indPlotBins_S),vecSDsCV_S(indPlotBins_S)./sqrt(vecCountsCV_S(indPlotBins_S)),'color',[0.5 0.5 0.5])
		scatter(vecMperTrial,vecCVperTrial(:),[],[0.5 0.5 1],'.');
		errorbar(vecActBinsC(indPlotBins),vecMeansCV(indPlotBins),vecSDsCV(indPlotBins)./sqrt(vecCountsCV(indPlotBins)),'color',[0 0 1])
		xlabel('Mean population activity (Hz)')
		ylabel('CV of population activity (Hz)')
		fixfig;
		
		%save fig
		export_fig(fullpath(strFigurePathSR,sprintf('Q2C_ISIDistros%s_%s.tif',strRec,strOnset)));
		export_fig(fullpath(strFigurePathSR,sprintf('Q2C_ISIDistros%s_%s.pdf',strRec,strOnset)))
		
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
		
		%% plot real vs shuffled
		cellMarkers = {'x','o'};
		
		figure;maxfig;
		[rML,pML]=corr(vecMperTrial(:),vecLperTrial(:));
		[rMH,pMH]=corr(vecMperTrial,vecHperTrial(:));
		[rHL,pHL]=corr(vecLperTrial(:),vecHperTrial(:));
		
		subplot(2,3,1)
		[vecCountsML,vecMeansML,vecSDsML] = makeBins(vecMperTrial,vecLperTrial,vecActBins);
		[vecCountsMH,vecMeansMH,vecSDsMH] = makeBins(vecMperTrial,vecHperTrial,vecActBins);
		[vecCountsML_S,vecMeansML_S,vecSDsML_S] = makeBins(vecMperTrial_S,vecLperTrial_S,vecActBins);
		[vecCountsMH_S,vecMeansMH_S,vecSDsMH_S] = makeBins(vecMperTrial_S,vecHperTrial_S,vecActBins);hold on
		indPlotBins = vecCountsMH>10;
		
		scatter(cat(1,vecMperTrial_S,vecMperTrial_S),cat(1,vecLperTrial_S(:),vecHperTrial_S(:)),[],[0.7 0.7 0.7],'.');
		scatter(vecMperTrial,vecLperTrial(:),[],[0.5 0.5 1],'.');
		scatter(vecMperTrial,vecHperTrial(:),[],[1 0.5 0.5],'.');
		
		errorbar(vecActBinsC(indPlotBins),vecMeansML_S(indPlotBins),vecSDsML_S(indPlotBins)./sqrt(vecCountsML_S(indPlotBins)),'color',[0.5 0.5 0.5])
		errorbar(vecActBinsC(indPlotBins),vecMeansMH_S(indPlotBins),vecSDsMH_S(indPlotBins)./sqrt(vecCountsMH_S(indPlotBins)),'color',[0.5 0.5 0.5])
		errorbar(vecActBinsC(indPlotBins),vecMeansML(indPlotBins),vecSDsML(indPlotBins)./sqrt(vecCountsML(indPlotBins)),'color',[0 0 1])
		errorbar(vecActBinsC(indPlotBins),vecMeansMH(indPlotBins),vecSDsMH(indPlotBins)./sqrt(vecCountsMH(indPlotBins)),'color',[1 0 0])
		
		hold off
		fixfig;
		xlabel('Mean firing rate during trial (Hz)');
		ylabel('FR in upper/lower 5% (Hz)');
		title(['Real rate fluctuations vs ' strType_S])
		fixfig;
		
		subplot(2,3,2)
		dblStep = 0.05;
		vecBinE = -1:dblStep:1.5;
		vecBinC = vecBinE(2:end)-dblStep/2;
		vecCountsLow = histcounts((vecLperTrial-vecLperTrial_S)./vecLperTrial_S,vecBinE);
		vecCountsHigh = histcounts((vecHperTrial-vecHperTrial_S)./vecHperTrial_S,vecBinE);
		plot(vecBinC*100,vecCountsLow,'color',[0 0 1]);
		hold on
		plot(vecBinC*100,vecCountsHigh,'color',[1 0 0]);
		hold off
		xlabel(['% change in FR over ' strType_S]);
		legend({'Lowest 5%','Highest 5%'});
		ylabel('Number of trials (count)');
		vecL = (vecLperTrial-vecLperTrial_S)./vecLperTrial_S;
		vecH = (vecHperTrial-vecHperTrial_S)./vecHperTrial_S;
		[h,pL]=ttest(vecL);
		[h,pH]=ttest(vecH);
		title(sprintf('Low, avg=%.1f%%, p=%.1e; high, avg=+%.1f%%, p=%.1e',mean(vecL)*100,pL,mean(vecH)*100,pH));
		fixfig;
		
		%real
		r1=corr(vecLperTrial(:),vecPopSparseness(:));
		r2=corr(vecMperTrial(:),vecPopSparseness(:));
		[rSpH,pSpH]=corr(vecMperTrial(:),vecPopSparseness(:));
		vecActBinsH = 0:dblActBinW:1700;
		vecActBinsHC = vecActBinsH(2:end)-dblActBinW/2;
		[vecCounts2,vecMeans2,vecSDs2] = makeBins(vecMperTrial,vecPopSparseness,vecActBinsH);
		%shuff
		r1_S=corr(vecLperTrial_S(:),vecPopSparseness_S(:));
		r2_S=corr(vecMperTrial_S(:),vecPopSparseness_S(:));
		[rSpH_S,pSpH_S]=corr(vecMperTrial_S(:),vecPopSparseness_S(:));
		vecActBinsH_S = 0:dblActBinW:1700;
		vecActBinsHC_S = vecActBinsH_S(2:end)-dblActBinW/2;
		[vecCounts2_S,vecMeans2_S,vecSDs2_S] = makeBins(vecMperTrial_S,vecPopSparseness_S,vecActBinsH_S);
		
		%fit
		indPlotBins2 = vecCounts2>10 | vecCounts2_S>10;
		
		mdl = fitlm(vecMperTrial,vecPopSparseness);
		ci = coefCI(mdl);
		[ypred,yci] = predict(mdl,vecActBinsHC(indPlotBins2)');
		
		mdl_S = fitlm(vecMperTrial_S,vecPopSparseness_S);
		ci_S = coefCI(mdl_S);
		[ypred_S,yci_S] = predict(mdl_S,vecActBinsHC_S(indPlotBins2)');
		
		subplot(2,3,3)
		hold on;
		scatter(vecMperTrial,vecPopSparseness,[],1-(1-lines(1))*(2/3),'.');
		scatter(vecMperTrial_S,vecPopSparseness_S,[],[0.7 0.7 0.7],'.');
		plot(vecActBinsHC(indPlotBins2),yci(:,1),'--','color',lines(1));
		plot(vecActBinsHC(indPlotBins2),ypred,'color',lines(1));
		plot(vecActBinsHC(indPlotBins2),yci(:,2),'--','color',lines(1));
		plot(vecActBinsHC_S(indPlotBins2),yci_S(:,1),'--','color',[0.5 0.5 0.5]);
		plot(vecActBinsHC_S(indPlotBins2),ypred_S,'color',[0.5 0.5 0.5]);
		plot(vecActBinsHC_S(indPlotBins2),yci_S(:,2),'--','color',[0.5 0.5 0.5]);
		hold off;
		hold off;
		ylabel('Pop. sparseness per trial');
		xlabel('Mean pop. firing rate per trial (Hz) ');
		title(sprintf('Real=%.2f, p=%.1e; %s=%.2f, p=%.1e',rSpH,pSpH,strType,rSpH_S,pSpH_S));
		fixfig;
		
		%split population in highest 50/lowest 50
		intUseUpperCells = min(sum(matResp>0,1));
		[intNeurons,intTrialNum] = size(matResp);
		vecHighAct = nan(intTrialNum,1);
		vecLowAct =  nan(intTrialNum,1);
		vecQuantiles3 = [1/3 1/2 2/3];
		vecQuantileIdx = round(vecQuantiles3*intNeurons);
		matQuantileAct = nan(intTrialNum,numel(vecQuantileIdx));
		vecMeanOfActiveCells = nan(intTrialNum,1);
		for intTrial=1:intTrialNum
			vecR = sort(matResp(:,intTrial));
			matQuantileAct(intTrial,:) = vecR(vecQuantileIdx);
			vecMeanOfActiveCells(intTrial) = mean(vecR((end-intUseUpperCells+1):end));
		end
		[vecHsorted,vecReorder]=sort(vecHperTrial);
		
		%save fig
		export_fig(fullpath(strFigurePathSR,sprintf('Q2D_QuantileDeviations%s_%s.tif',strRec,strOnset)));
		export_fig(fullpath(strFigurePathSR,sprintf('Q2D_QuantileDeviations%s_%s.pdf',strRec,strOnset)))
			error do fit per recording, then plot fits
	end
end

%% plot mean over recordings
%error plot matSlopes

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
for intVar=1:4
	errorbar(vecBinD_c,matMeanISIcorrs(:,intVar),matSdISIcorrs(:,intVar)./sqrt(intRecNum));
end
hold off
xlabel('ISI spikes i,i+1 (ms)');
ylabel('ISI spikes i+1,i+2 (ms)');
title(sprintf('Mean +/- SEM over recs'));
fixfig;

legend({strType,strType_S,strType_U,strType_P},'location','best');
		
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

% mean and sd scale together: CV is constant
%do analysis where you shuffle trials independently for each neuron to test if blue curve is population effect
%error calculate fit on binned points
intK = 1;
vecBinX = vecActBinsC';
subplot(2,3,4);cla;

%poisson
matMeansSdRem_P = matMeansSd_P;
matMeansSdRem_P(matCountsCV_P<10)=nan;
matMeansCVRem_P = matMeansCV_P;
matMeansCVRem_P(matCountsCV_P<10)=nan;
vecMeansSd_P = nanmean(matMeansSdRem_P,2);
indKeep_P = ~isnan(vecMeansSd_P);
vecMeansSd_P = vecMeansSd_P(indKeep_P);
vecX_P = vecBinX(indKeep_P);
dblSlopeBins_P = ((vecX_P' * vecX_P) \ vecX_P') * vecMeansSd_P;
vecFitY_P = vecX_P*dblSlopeBins_P;
[dblR2_P,dblSS_tot_P,dblSS_res_P,dblT_P,dblP_P,dblR2_adjusted_P,dblR2_SE_P] = getR2(vecMeansSd_P,vecFitY_P,intK);

%neuron-wise ISI-shuffled
matMeansSdRem_S = matMeansSd_S;
matMeansSdRem_S(matCountsCV_S<10)=nan;
matMeansCVRem_S = matMeansCV_S;
matMeansCVRem_S(matCountsCV_S<10)=nan;
vecMeansSd_S = nanmean(matMeansSdRem_S,2);
indKeep_S = ~isnan(vecMeansSd_S);
vecMeansSd_S = vecMeansSd_S(indKeep_S);
vecX_S = vecBinX(indKeep_S);
dblSlopeBins_S = ((vecX_S' * vecX_S) \ vecX_S') * vecMeansSd_S;
vecFitY_S = vecX_S*dblSlopeBins_S;
[dblR2_S,dblSS_tot_S,dblSS_res_S,dblT_S,dblP_S,dblR2_adjusted_S,dblR2_SE_S] = getR2(vecMeansSd_S,vecFitY_S,intK);

%trial-shuffled
matMeansSdRem_U = matMeansSd_U;
matMeansSdRem_U(matCountsCV_U<10)=nan;
matMeansCVRem_U = matMeansCV_U;
matMeansCVRem_U(matCountsCV_U<10)=nan;
vecMeansSd_U = nanmean(matMeansSdRem_U,2);
indKeep_U = ~isnan(vecMeansSd_U);
vecMeansSd_U = vecMeansSd_U(indKeep_U);
vecX_U = vecBinX(indKeep_U);
dblSlopeBins_U = ((vecX_U' * vecX_U) \ vecX_U') * vecMeansSd_U;
vecFitY_U = vecX_U*dblSlopeBins_U;
[dblR2_U,dblSS_tot_U,dblSS_res_U,dblT_U,dblP_U,dblR2_adjusted_U,dblR2_SE_U] = getR2(vecMeansSd_U,vecFitY_U,intK);

%real
matMeansSdRem = matMeansSd;
matMeansSdRem(matCountsCV<10)=nan;
matMeansCVRem = matMeansCV;
matMeansCVRem(matCountsCV<10)=nan;
vecMeansSd = nanmean(matMeansSdRem,2);
indKeep = ~isnan(vecMeansSd);
vecMeansSd = vecMeansSd(indKeep);
vecX = vecBinX(indKeep);
dblSlopeBins = ((vecX' * vecX) \ vecX') * vecMeansSd;
vecFitY = vecX*dblSlopeBins;
[dblR2,dblSS_tot,dblSS_res,dblT,dblP,dblR2_adjusted,dblR2_SE] = getR2(vecMeansSd,vecFitY,intK);

%plot
matC = lines(4);
hold on
errorbar(vecActBinsC,nanmean(matMeansSdRem,2),nanstd(matMeansSdRem,[],2)./sqrt(intRecNum),'linestyle','none','color',matC(1,:));
errorbar(vecActBinsC,nanmean(matMeansSdRem_S,2),nanstd(matMeansSdRem_S,[],2)./sqrt(intRecNum),'linestyle','none','color',matC(2,:));
errorbar(vecActBinsC,nanmean(matMeansSdRem_U,2),nanstd(matMeansSdRem_U,[],2)./sqrt(intRecNum),'linestyle','none','color',matC(3,:));
errorbar(vecActBinsC,nanmean(matMeansSdRem_P,2),nanstd(matMeansSdRem_P,[],2)./sqrt(intRecNum),'linestyle','none','color',matC(4,:));
hold off
legend({strType,strType_S,strType_U,strType_P},'location','best');
xlabel('Mean of spiking rate (Hz)');
ylabel('Sd of spiking rate (Hz)');
title(sprintf('Constant coefficient of variation (CV) = %.3f',dblSlopeBins));
fixfig;

h=subplot(2,3,5);delete(h);subplot(2,3,5);
hold on;
errorbar(1,dblR2,dblR2_SE,'x','color',matC(1,:),'capsize',20);
errorbar(2,dblR2_S,dblR2_SE_S,'x','color',matC(2,:),'capsize',20);
errorbar(3,dblR2_U,dblR2_SE_U,'x','color',matC(3,:),'capsize',20);
errorbar(4,dblR2_P,dblR2_SE_P,'x','color',matC(4,:),'capsize',20);
xlim([0 5]);
set(gca,'xtick',1:4,'xticklabel',{strType,strType_S,strType_U,strType_P});
ylabel('Linearity of sd/mean (R^2)');
fixfig;

subplot(2,3,6)
matSlopesNorm = matSlopes./matSlopes(4,:);
hold on;
vecX = ones(size(matSlopesNorm(1,:)));
swarmchart(vecX,matSlopesNorm(1,:),[30],matC(1,:),'marker','o','JitterWidth',0.2);
errorbar(1,mean(matSlopesNorm(1,:),2),std(matSlopesNorm(1,:),[],2)./sqrt(intRecNum),'x','color',matC(1,:),'capsize',20);
swarmchart(2*vecX,matSlopesNorm(2,:),[],matC(2,:),'marker','o','JitterWidth',0.2);
errorbar(2,mean(matSlopesNorm(2,:),2),std(matSlopesNorm(2,:),[],2)./sqrt(intRecNum),'x','color',matC(2,:),'capsize',20);
swarmchart(3*vecX,matSlopesNorm(3,:),[],matC(3,:),'marker','o','JitterWidth',0.2);
errorbar(3,mean(matSlopesNorm(3,:),2),std(matSlopesNorm(3,:),[],2)./sqrt(intRecNum),'x','color',matC(3,:),'capsize',20);
swarmchart(4*vecX,matSlopesNorm(4,:),[],matC(4,:),'marker','o','JitterWidth',0.2);
errorbar(4,mean(matSlopesNorm(4,:),2),std(matSlopesNorm(4,:),[],2)./sqrt(intRecNum),'x','color',matC(4,:),'capsize',20);
xlim([0 5]);
set(gca,'xtick',1:4,'xticklabel',{strType,strType_S,strType_U,strType_P});
ylabel('Slope of sd/mean');
fixfig;

drawnow;
export_fig(fullpath(strFigurePath,sprintf('Q2_ConstantCV.tif')));
export_fig(fullpath(strFigurePath,sprintf('Q2_ConstantCV.pdf')));
