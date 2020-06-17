%% description
%{
    - The lack of a 2nd bump in the visual miss trials was expected. I find
    it very interesting that there may also be a difference between hits
    and misses for what pertains the first bump. Do you think it would be
    possible to look at B1_hit / B1_miss alignment? Of course we could also
    try a simpler approach (just look at whether the trajectories in a
    latent space differ between hits and misses).

When you look at the absolute numbers (in the non-normalized figures) in
figures 3B/4B, there actually doesn't seem to be that much difference in B1
dynamics between hits and misses: it's mostly the presence/absence of B2
that changes the relative size of the B1 trajectories. But of course the
B1-B2 are different logistic regressors, so we can't say for certain based
on these graphs.     
I think the right analysis to answer your question here would be to do a
reduced rank regression like before to compare the performance of
B1_hit(B1_hit)->B1_hit with B1_hit(B1_miss)->B1_hit (and the other way
around). This would allow us to quantify the difference in neural subspaces
between B1 hits and misses.
%}

%% clear
clear all;

%% analysis
%is B1 activity hit-specific?
%predict B1_hit(B1_hit)->B1_hit and compare with B1_hit(B1_miss)->B1_hit => x% generalizable
%is B1 hit-specificity larger than B2 hit-specificity?
%predict B2_hit(B2_hit)->B2_hit and compare with B2_hit(B2_miss)->B2_hit  => y% generalizable
%x < y? p=?

%% parameters
%On which timestamp to align as t=0
strAlignOn = 'Change';
intFullIters = 1;
intResamplings = 10;
boolShuffle = false;
dblLambda = 500;
boolSavePlots = true;

%% HEADER, load data
MOL_Header;

%% Remove last 20 trials:
trialData = MOL_RemoveLastnTrials(trialData,20);

%% Filter good neurons:
spikeData = MOL_filterNeurons(sessionData,trialData,spikeData);

%% split data by session
[cellSessions,ia,vecSessionIdx] = unique(spikeData.session_ID);

%% split by area
[cellAreas,ia,vecAreaIdx] = unique(spikeData.area);

%% Get significant first and second bump responses:
%[sessionData,trialData,spikeData] = MOL_calc_signresp(sessionData,trialData,spikeData);

%% make populations
[cellPopID,ia,vecPopulationIdx] = unique(vecSessionIdx*10 + vecAreaIdx);

%% define V1 populations
indV1=vecAreaIdx==find(strcmp(cellAreas,'V1'));
vecPopsV1 = unique(vecPopulationIdx(indV1));
intPops = numel(vecPopsV1);

%% pre-allocate output
vecCellsPerPop = nan(1,intPops);
vecFracCorr = nan(1,intPops);
vecTrialNum = nan(1,intPops);
vecHitTrialNum = nan(1,intPops);
vecMissTrialNum = nan(1,intPops);

cellM_Ep1R2_hh2h = cell(1,intPops);
cellM_Ep1R2_hm2h = cell(1,intPops);
cellM_Ep2R2_hh2h = cell(1,intPops);
cellM_Ep2R2_hm2h = cell(1,intPops);

cellM_yy2y = cell(1,intPops);
cellM_yy2x = cell(1,intPops);
cellM_yRRR2x = cell(1,intPops);
cellM_yx2x = cell(1,intPops);

vecRemoveSessions = [3 6 9]; %due to bad performance

%%  analyze
dblTimeConversion = 10^6;
for intPopIdx=1:intPops
	%% prepare, retrieve and transform data
	%retrieve cells
	intPopulation = vecPopsV1(intPopIdx);
	indUseCells = vecPopulationIdx==intPopulation;
	intCells = sum(indUseCells);
	intCutOffCellNr = 10;
	vecCellsPerPop(intPopIdx) = intCells;
	if intCells < intCutOffCellNr
		fprintf('Pop %d, Number of cells is %d, which is under %d, skipping... [%s]\n',intPopIdx,intCells,intCutOffCellNr,getTime);
		continue;
	end
	
	%get spikes and transform to seconds
	cellSpikes = cellfun(@rdivide,spikeData.ts(indUseCells),cellfill(dblTimeConversion,[intCells 1]),'uniformoutput',false);
	strSesID = unique(spikeData.session_ID(indUseCells));
	if numel(strSesID)~=1,error([mfilename ':SessionNotUnique'],'Something went wrong with population selection...');end
	
	%get trial data
	indSesTrials = strcmp(trialData.session_ID,strSesID);
	intSesTrials = sum(indSesTrials);
	
	vecTrialType = trialData.trialType(indSesTrials);
	vecTrialVisOriPre = trialData.visualOriPreChange(indSesTrials);
	vecTrialVisOriChange = trialData.visualOriChange(indSesTrials);
	vecTrialAudFreqPre = trialData.audioFreqPreChange(indSesTrials);
	vecTrialAudFreqChange = trialData.audioFreqChange(indSesTrials);
	vecTrialCorrect = trialData.correctResponse(indSesTrials);
	vecTrialRespSide = trialData.responseSide(indSesTrials);
	
	vecTrialLickSecs = cellfun(@rdivide,trialData.lickTime(indSesTrials),cellfill(dblTimeConversion,[intSesTrials 1]),'uniformoutput',false);
	
	vecTrialStartSecs = trialData.trialStart(indSesTrials)./dblTimeConversion;
	vecTrialChangeSecs = trialData.stimChange(indSesTrials)./dblTimeConversion;
	vecTrialEndSecs = trialData.trialEnd(indSesTrials)./dblTimeConversion;
	vecTrialStartITISecs = trialData.itiStart(indSesTrials)./dblTimeConversion;
	vecTrialEndITISecs = trialData.itiEnd(indSesTrials)./dblTimeConversion;
	vecTrialRespSecs = trialData.responseLatency(indSesTrials)./dblTimeConversion;
	
	indKeepTrials = abs(vecTrialVisOriChange) == 90;
	dblCutOff = 0;
	dblFracCorr = sum(vecTrialCorrect(indKeepTrials))/sum(indKeepTrials);
	vecFracCorr(intPopIdx) = dblFracCorr;
	if dblFracCorr < dblCutOff%|| ismember(intPopIdx,vecRemoveSessions)
		fprintf('Pop %d, Correct response proportion is %.3f, which is under %.3f, skipping... [%s]\n',intPopIdx,dblFracCorr,dblCutOff,getTime);
		continue;
	end
	%message
	fprintf('Analyzing pop %d/%d; %d neurons, %d trials [%s]\n',intPopIdx,intPops,intCells,sum(indKeepTrials),getTime);
	vecTrialNum(intPopIdx) = sum(indKeepTrials);
	
	%define 1st vs 2nd bump windows
	if strcmp(strAlignOn,'Change')
		vecT0 = vecTrialChangeSecs(indKeepTrials);
	end
	
	%build vectors
	intTrials = numel(vecT0);
	vecStartBaseT = vecT0 + dblStartBaseT;
	vecStopBaseT = vecT0 + dblStopBaseT;
	vecStartEp1T = vecT0 + dblStartEp1T;
	vecStopEp1T = vecT0 + dblStopEp1T;
	vecStartEp2T = vecT0 + dblStartEp2T;
	vecStopEp2T = vecT0 + dblStopEp2T;
	vecEdges = sort(cat(1,vecStartBaseT,vecStopBaseT,vecStartEp1T,vecStopEp1T,vecStartEp2T,vecStopEp2T));
	vecBaseEdges = sort(cat(1,vecStartBaseT,vecStopBaseT));
	vecEp1Edges = sort(cat(1,vecStartEp1T,vecStopEp1T));
	vecEp2Edges = sort(cat(1,vecStartEp2T,vecStopEp2T));
	%calculate response matrix
	matRespBase = nan(intCells,intTrials);
	matRespEp1 = nan(intCells,intTrials);
	matRespEp2 = nan(intCells,intTrials);
	for intNeuron=1:intCells
		vecBaseR = histcounts(cellSpikes{intNeuron},vecBaseEdges);
		matRespBase(intNeuron,:) = vecBaseR(1:2:end);
		vecEp1R = histcounts(cellSpikes{intNeuron},vecEp1Edges);
		matRespEp1(intNeuron,:) = vecEp1R(1:2:end);
		vecEp2R = histcounts(cellSpikes{intNeuron},vecEp2Edges);
		matRespEp2(intNeuron,:) = vecEp2R(1:2:end);
	end
	
	%% sub-select trials
	%is B1 activity hit-specific?
%predict B1_hit(B1_hit)->B1_hit and compare with B1_hit(B1_miss)->B1_hit => x% generalizable
%is B1 hit-specificity larger than B2 hit-specificity?
%predict B2_hit(B2_hit)->B2_hit and compare with B2_hit(B2_miss)->B2_hit  => y% generalizable
%x < y? p=?

%sub-select trials that are large change, visual, and hit
	%sub-select trials that are large change, visual, and miss
	%split into B1 and B2
	indHits = vecTrialCorrect(indKeepTrials)==1;
	indMiss = isnan(vecTrialRespSecs(indKeepTrials));
	vecHitTrialNum(intPopIdx) = sum(indHits);
	vecMissTrialNum(intPopIdx) = sum(indMiss);

	%analysis 1
	matRespEp1Hits = matRespEp1(:,indHits); %B1
	matRespEp1Miss = matRespEp1(:,indMiss); %B2
	
	%analysis 2
	matRespEp2Hits = matRespEp2(:,indHits); %B1
	matRespEp2Miss = matRespEp2(:,indMiss); %B2
	
		%% analyze
		%B1_hit(B1_hit)->B1_hit and compare with B1_hit(B1_miss)->B1_hit
		%B1_hit(B1_hit)->B1_hit: getRegInSpaceCV(B1_hit_train,PC(B1_hit_train),B1_hit_train,B1_hit_test,B1_hit_test)
		%B1_hit(B1_miss)->B1_hit: getRegInSpaceCV(B1_hit_train,PC(B1_miss_train),B1_hit_train,B1_hit_test,B1_hit_test)
		
		%REPEAT FOR:
		%B2_hit(B2_hit)->B2_hit: getRegInSpaceCV(B2_hit_train,PC(B2_hit_train),B1_hit_train,B1_hit_test,B1_hit_test)
		%B2_hit(B2_miss)->B2_hit: getRegInSpaceCV(B2_hit_train,PC(B2_miss_train),B2_hit_train,B2_hit_test,B2_hit_test)
		[matEp1R2_hh2h,matEp1R2_hm2h] = doDimCompareSpaceCV_HitsMiss(matRespEp1Hits,matRespEp1Miss,intResamplings,dblLambda,boolShuffle);
		[matEp2R2_hh2h,matEp2R2_hm2h] = doDimCompareSpaceCV_HitsMiss(matRespEp2Hits,matRespEp2Miss,intResamplings,dblLambda,boolShuffle);
		
		%% save output
		cellM_Ep1R2_hh2h{1,intPopIdx} = matEp1R2_hh2h;
		cellM_Ep1R2_hm2h{1,intPopIdx} = matEp1R2_hm2h;
		cellM_Ep2R2_hh2h{1,intPopIdx} = matEp2R2_hh2h;
		cellM_Ep2R2_hm2h{1,intPopIdx} = matEp2R2_hm2h;
end
warning('on','MATLAB:nearlySingularMatrix')

%% process data
cellM2_Ep1R2_hh2h = cell(1,intPops);
cellM2_Ep1R2_hm2h = cell(1,intPops);
cellM2_Ep2R2_hh2h = cell(1,intPops);
cellM2_Ep2R2_hm2h = cell(1,intPops);

cellS2_Ep1R2_hh2h = cell(1,intPops);
cellS2_Ep1R2_hm2h = cell(1,intPops);
cellS2_Ep2R2_hh2h = cell(1,intPops);
cellS2_Ep2R2_hm2h = cell(1,intPops);

% AUCs
matAUC_Ep1R2_hh2h = nan(intResamplings,intPops);
matAUC_Ep1R2_hm2h = nan(intResamplings,intPops);
matAUC_Ep2R2_hh2h = nan(intResamplings,intPops);
matAUC_Ep2R2_hm2h = nan(intResamplings,intPops);

for intPopIdx=1:intPops
	for intResampling=1:intResamplings
		if isempty(cellM_Ep1R2_hh2h{1,intPopIdx}),continue;end
		vecEp1R2_hh2h = cellM_Ep1R2_hh2h{1,intPopIdx}(intResampling,:);
		matAUC_Ep1R2_hh2h(intResampling,intPopIdx) = sum((vecEp1R2_hh2h)/vecEp1R2_hh2h(end)) / numel(vecEp1R2_hh2h);
		
		vecEp1R2_hm2h = cellM_Ep1R2_hm2h{1,intPopIdx}(intResampling,:);
		matAUC_Ep1R2_hm2h(intResampling,intPopIdx) = sum((vecEp1R2_hm2h)/vecEp1R2_hm2h(end)) / numel(vecEp1R2_hm2h);
		
		vecEp2R2_hh2h = cellM_Ep2R2_hh2h{1,intPopIdx}(intResampling,:);
		matAUC_Ep2R2_hh2h(intResampling,intPopIdx) = sum((vecEp2R2_hh2h)/vecEp2R2_hh2h(end)) / numel(vecEp2R2_hh2h);
		
		vecEp2R2_hm2h = cellM_Ep2R2_hm2h{1,intPopIdx}(intResampling,:);
		matAUC_Ep2R2_hm2h(intResampling,intPopIdx) = sum((vecEp2R2_hm2h)/vecEp2R2_hm2h(end)) / numel(vecEp2R2_hm2h);
	end
	
	%mean curves
	%process data
	cellM2_Ep1R2_hh2h{intPopIdx} = mean(cell2mat(cellM_Ep1R2_hh2h(:,intPopIdx)),1);
	cellS2_Ep1R2_hh2h{intPopIdx} = std(cell2mat(cellM_Ep1R2_hh2h(:,intPopIdx)),[],1)./sqrt(intResamplings);
	
	cellM2_Ep1R2_hm2h{intPopIdx} = mean(cell2mat(cellM_Ep1R2_hm2h(:,intPopIdx)),1);
	cellS2_Ep1R2_hm2h{intPopIdx} = std(cell2mat(cellM_Ep1R2_hm2h(:,intPopIdx)),[],1)./sqrt(intResamplings);
	
	cellM2_Ep2R2_hh2h{intPopIdx} = mean(cell2mat(cellM_Ep2R2_hh2h(:,intPopIdx)),1);
	cellS2_Ep2R2_hh2h{intPopIdx} = std(cell2mat(cellM_Ep2R2_hh2h(:,intPopIdx)),[],1)./sqrt(intResamplings);
	
	cellM2_Ep2R2_hm2h{intPopIdx} = mean(cell2mat(cellM_Ep2R2_hm2h(:,intPopIdx)),1);
	cellS2_Ep2R2_hm2h{intPopIdx} = std(cell2mat(cellM_Ep2R2_hm2h(:,intPopIdx)),[],1)./sqrt(intResamplings);
end

%% plot
%prep plot
%close all;
figure
vecAxPtr=nan(1,4);
for intPlot=1:4
	vecAxPtr(intPlot) = subplot(2,2,intPlot);
	hold on;
end
matColors = [0.1 0.1 0.8;...
	0.1 0.5 0.8;...
	0.5 0.1 0.8;...
	0.4 0.4 0.8;...
	0.8 0.1 0.1;...
	0.8 0.5 0.1;...
	0.8 0.1 0.5;...
	0.8 0.4 0.4];

%labels
cellLegend = {'A) Epoch1 Hit(Hit)-Hit',...
	'B) Epoch1 Hit(Miss)-Hit',...
	'C) Epoch2 Hit(Hit)-Hit',...
	'D) Epoch2 Hit(Miss)-Hit'...
	};

cellLabelsX = {['B1 hit-space, sorted by Hit-' getGreek('lambda','lower')],...
	['B1 hit-space, sorted by Miss-' getGreek('lambda','lower')],...
	['B2 hit-space, sorted by Hit-' getGreek('lambda','lower')],...
	['B2 hit-space, sorted by Miss-' getGreek('lambda','lower')],...
	};

cellLabelsY = {['Fraction of B1 hit-activity (R^2)'],...
	['Fraction of B1 hit-activity (R^2)'],...
	['Fraction of B1 hit-activity (R^2)'],...
	['Fraction of B1 hit-activity (R^2)'],...
	};

% plot with error bars
%{
for intPopIdx=1:intPops
	vecX = (1:numel(cellM2_xx2x{intPopIdx}))';%1:intMaxDim;%linspace(0,1,intMaxDim);
	errorbar(vecAxPtr(1),vecX,cellM2_xx2x{intPopIdx},cellS2_xx2x{intPopIdx},'Color',matColors(1,:));
	errorbar(vecAxPtr(2),vecX,cellM2_xx2y{intPopIdx},cellS2_xx2y{intPopIdx},'Color',matColors(2,:));
	errorbar(vecAxPtr(3),vecX,cellM2_xRRR2y{intPopIdx},cellS2_xRRR2y{intPopIdx},'Color',matColors(3,:));
	errorbar(vecAxPtr(4),vecX,cellM2_xy2y{intPopIdx},cellS2_xy2y{intPopIdx},'Color',matColors(4,:));
	
	errorbar(vecAxPtr(5),vecX,cellM2_yy2y{intPopIdx},cellS2_yy2y{intPopIdx},'Color',matColors(5,:));
	errorbar(vecAxPtr(6),vecX,cellM2_yy2x{intPopIdx},cellS2_yy2x{intPopIdx},'Color',matColors(6,:));
	errorbar(vecAxPtr(7),vecX,cellM2_yRRR2x{intPopIdx},cellS2_yRRR2x{intPopIdx},'Color',matColors(7,:));
	errorbar(vecAxPtr(8),vecX,cellM2_yx2x{intPopIdx},cellS2_yx2x{intPopIdx},'Color',matColors(8,:));
end
%}
% plot only means
for intPopIdx=1:intPops
	vecX = (1:numel(cellM2_Ep1R2_hh2h{intPopIdx}))';%1:intMaxDim;%linspace(0,1,intMaxDim);
	plot(vecAxPtr(1),vecX,cellM2_Ep1R2_hh2h{intPopIdx},'Color',matColors(1,:));
	plot(vecAxPtr(2),vecX,cellM2_Ep1R2_hm2h{intPopIdx},'Color',matColors(2,:));
	plot(vecAxPtr(3),vecX,cellM2_Ep2R2_hh2h{intPopIdx},'Color',matColors(3,:));
	plot(vecAxPtr(4),vecX,cellM2_Ep2R2_hm2h{intPopIdx},'Color',matColors(4,:));
end

% finish plots
for intPlot=1:4
	axes(vecAxPtr(intPlot));
	hold off
	%ylim([0 max()]);
	title(cellLegend{intPlot});
	fixfig;
	xlabel(cellLabelsX{intPlot});
	ylabel(cellLabelsY{intPlot});
	if mod(intPlot,4) == 1
		ylim([0 0.3]);
	else
		ylim([0 0.3]);
	end
end

%save
maxfig();
if boolSavePlots
	strFigFile = ['DimCom5ST_HitMiss' getDate];
	export_fig([strFigDir strFigFile '.tif']);
	export_fig([strFigDir strFigFile '.pdf']);
end

%% is B1 activity hit-specific?
%predict B1_hit(B1_hit)->B1_hit and compare with B1_hit(B1_miss)->B1_hit => x% generalizable
%is B1 hit-specificity larger than B2 hit-specificity?
%predict B2_hit(B2_hit)->B2_hit and compare with B2_hit(B2_miss)->B2_hit  => y% generalizable
%x < y? p=?

% summary figure
figure
matHitSpecificityB1 = matAUC_Ep1R2_hh2h - matAUC_Ep1R2_hm2h;
matHitSpecificityB2 = matAUC_Ep2R2_hh2h - matAUC_Ep2R2_hm2h;

intUsedPops = sum(any(~isnan(matAUC_Ep1R2_hh2h),1));

vecHitSpec1 = mean(matHitSpecificityB1);
vecHitSpec2 = mean(matHitSpecificityB2);
vecHitSpecDiff = vecHitSpec2 - vecHitSpec1;

dblR_corr = nancorr(vecHitSpecDiff',vecFracCorr');
dblR_trials= nancorr(vecHitSpecDiff',vecTrialNum');
dblR_cells = nancorr(vecHitSpecDiff',vecCellsPerPop');
dblR_hits = nancorr(vecHitSpecDiff',vecHitTrialNum');
dblR_miss = nancorr(vecHitSpecDiff',vecMissTrialNum');

[h,p1]=ttest(vecHitSpecDiff)
subplot(2,2,1)
plot([1.1 1.9],[vecHitSpec1' vecHitSpec2'],'color',[0.5 0.5 0.5]);
hold on
errorbar([1.05 1.95],[nanmean(vecHitSpec1) nanmean(vecHitSpec2)],[nanstd(vecHitSpec1) nanstd(vecHitSpec2)]./sqrt(intUsedPops),'b');
hold off
ylabel('d(AUC) between hit/miss subspaces');
title(sprintf('(B1 dAUC) vs (B2 dAUC), p=%.3f',p1));
set(gca,'xtick',[1.1 1.9],'xticklabel',{'B1','B2'});
fixfig;

vecHitSpecDiff2 = (vecHitSpec2 ./ vecHitSpec1) - 1;
[h,p2]=ttest(vecHitSpecDiff2);
subplot(2,2,2)
bplot(vecHitSpecDiff2,'linewidth',35);
set(gca,'xtick',[])
ylabel('Hit-specificity index');
title(sprintf('index=((B1 dAUC) / (B2 dAUC)) - 1, p=%.3f',p2));
fixfig;
maxfig;


if boolSavePlots
	strFigFile = ['DimCom5ST_Generalizability' getDate];
	export_fig([strFigDir strFigFile '.tif']);
	export_fig([strFigDir strFigFile '.pdf']);
end



