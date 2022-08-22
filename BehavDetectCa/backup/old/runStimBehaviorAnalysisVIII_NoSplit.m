%function vecPredictorAssemblyZscores = doPlotActivityDistributions(cellContrastDataAD,sMetaDataAD,sPC,sParams)

%if ~exist('sParams','var'),sParams = struct;end
%if isfield(sParams,'intStartSes'),intStartSes=sParams.intStartSes;else intStartSes=1;end
%if isfield(sParams,'intStopSes'),intStopSes=sParams.intStopSes;else intStopSes=length(cellContrastDataAD);end
close all
boolUseNeuropilSubtraction = false;
global vecAllH;
vecAllH = [];
global intPopCounter;
intPopCounter = 0;
for intMouse=1:8;%[1 2 4 5 6]
	close all;
	clearvars -except intMouse boolUseNeuropilSubtraction vecAllH intPopCounter
	%get block data
	if ~exist('cellMultiSes','var')
		if intMouse == 1
			strSes = '20140207';
		elseif intMouse == 2
			strSes = '20140314';
		elseif intMouse == 3
			strSes = '20140425';
		elseif intMouse == 4
			strSes = '20140507';
		elseif intMouse == 5
			strSes = '20140530';
		elseif intMouse == 6
			strSes = '20140604';
		elseif intMouse == 7
			strSes = '20140711';
		elseif intMouse == 8
			strSes = '20140715';
		end
		if boolUseNeuropilSubtraction
			strSes = ['NPS' strSes];
		end
		load(['D:\Data\Results\stimdetection\dataPreProAggregate' strSes '.mat']);
	end
	
	%vecBlock = [1 1 1 2 2 2 2 2]; %define whether neurons are in same population or not
	vecBlockTypes = unique(vecBlock);
	intNumBlocks = length(vecBlockTypes);
	%vecNeuronNum = zeros(1,intNumBlocks);
	%cellKeepList = cell(1,intNumBlocks);
	%#ok<*ASGLU>
	%#ok<*AGROW>
	clear sLoad sSesAggregate ses
	sParams.strFigDir = ['D:\Data\Results\stimdetection' filesep strSes filesep];
	sParams.boolSavePlots = true;
	sParams.boolSaveData = true;
	strOldDir = cd(sParams.strFigDir);
	
	%change name for no split
	strSes = ['NS_Extra_' strSes];
	
	%remove last block for mouse 8
	if intMouse==8
		%remove last trial
		vecRem = true(size(cellMultiSes{1}.structStim.Orientation));
		vecRem((end-(48*1)+1):end) = false;
		cellMultiSes{1}.structStim = remel(cellMultiSes{1}.structStim,vecRem);
	end
	cellMultiSource = cellMultiSes;
	for intPopulation = vecBlockTypes
		close all;
		fprintf('Processing mouse %d/%d [pop %d] [%s]\n',intMouse,8,intPopulation,getTime);
		cellMultiSes = cellMultiSource;
		%[indTuned,vecNeuronPrefStim,vecOSI,cellAggregate] = getTunedStimDetectionNeurons(cellSes(vecBlock==intPopulation));
		[indTuned,vecNeuronPrefStim,vecOSI,cellAggregate] = getTunedStimDetectionNeurons(cellMultiSes{intPopulation});
		structStim = cellMultiSes{intPopulation}.structStim;
		structStim.vecTrialRepetition = ceil((1:(length(structStim.Orientation)))/48);
		
		% remove trials with reaction time <100ms
		dblRemSecs = 0.15;
		indTooFast = (structStim.vecTrialRespSecs-structStim.SecsOn)<dblRemSecs & structStim.vecTrialResponse == 1;
		fprintf('Removed %d trials with responses <%dms\n',sum(indTooFast),round(dblRemSecs*1000));
		%structStim.vecTrialResponse(indTooFast) = 0;
		structStim = remel(structStim,~indTooFast);
		%structStim.FrameOff = structStim.FrameOn+1;
		
		%remove trials that were too slow
		dblRemSecs = 3;
		indTooSlow = (structStim.vecTrialRespSecs-structStim.SecsOn)>dblRemSecs;
		fprintf('Removed %d trials with responses >%dms\n',sum(indTooSlow),round(dblRemSecs*1000));
		structStim.vecTrialResponse(indTooSlow) = 0;
		structStim.vecTrialRespSecs(indTooSlow) = nan;
		%structStim.FrameOff = structStim.FrameOn + 10;
		
		%take opposite directions as the same
		structStim.Orientation = mod(structStim.Orientation,180);
		vecOrientations = unique(structStim.Orientation);
		cellMultiSes{intPopulation}.structStim = structStim;
		
		%get orientation-based trial selection vectors
		sTypesOri = getStimulusTypes(structStim,{'Orientation'});
		cellSelectOri = getSelectionVectors(structStim,sTypesOri);
		
		%remove non-tuned neurons
		cellMultiSes{intPopulation}.neuron(~indTuned) = [];
		vecNeuronPrefStim(~indTuned) = [];
		intNeurons = numel(cellMultiSes{intPopulation}.neuron);
		intOris = length(vecOrientations);
		
		%% get pre-formatted data variables
		%get neuronal responses per trial
		[matTrialResponse,cellSelectContrasts,vecTrialDur] = getTrialResponseData(cellMultiSes{intPopulation},structStim);
		intContrasts = length(cellSelectContrasts);
		
		% get stimulus data
		[vecStimOris,dummy]=find(bsxfun(@eq,structStim.Orientation,vecOrientations'));
		vecStimOris = vecStimOris';
		
		[vecStimContrasts,dummy]=find(bsxfun(@eq,structStim.Contrast,unique(structStim.Contrast)'));
		vecStimContrasts = vecStimContrasts';
		vecTrialRepetition = structStim.vecTrialRepetition;
		
		vecRTs = cellMultiSes{intPopulation}.structStim.vecTrialRespSecs - cellMultiSes{intPopulation}.structStim.SecsOn; %get RTs for hit trials and selected contrasts
		indSelectRespTrials = ~isnan(vecRTs);
		intTrials = length(indSelectRespTrials);
		
		%% pre-analyses
		%normalize
		matRespNormPerNeuron = zscore(matTrialResponse,[],2);
		
		%get R
		matR = nan(intNeurons,intNeurons,intTrials);
		for intNeuron1=1:intNeurons
			for intNeuron2=intNeuron1:intNeurons
				matR(intNeuron1,intNeuron2,:) = matRespNormPerNeuron(intNeuron1,:).*matRespNormPerNeuron(intNeuron2,:);
			end
		end
		
		%get measures per trial
		matSelectAllR = false(size(matR));
		vecR = nan(intTrials,1);
		vecR_SD = nan(intTrials,1);
		matSelect = tril(true(intNeurons,intNeurons),-1);
		vecHeterogeneity = nan(intTrials,1);
		vecHetRawAct = nan(intTrials,1);
		vecHetMultiDim = nan(intTrials,1);
		vecHetMultiDim2 = nan(intTrials,1);
		vecActPrefPop = nan(intTrials,1);
		vecActWholePop = nan(intTrials,1);
		vecZActAll = nan(intTrials,1);
		vecActZPrefPop = nan(intTrials,1);
		for intTrial=1:intTrials
			vecZ = matRespNormPerNeuron(:,intTrial);
			matZ1 = repmat(vecZ,[1 intNeurons]);
			matZ2 = repmat(vecZ',[intNeurons 1]);
			matDist = abs(matZ1 - matZ2);
			vecHeterogeneity(intTrial) = mean(matDist(matSelect));
			%vecHetMultiDim(intTrial) = sqrt(sum(abs(vecZ-mean(vecZ))))/length(vecZ);
			vecHetMultiDim(intTrial) = sqrt(sum((vecZ-mean(vecZ)).^2))/length(vecZ);
			
			vecAct = matTrialResponse(:,intTrial);
			matA1 = repmat(vecAct,[1 intNeurons]);
			matA2 = repmat(vecAct',[intNeurons 1]);
			matDistA = abs(matA1 - matA2);
			vecHetRawAct(intTrial) = mean(matDistA(matSelect));
			
			vecHetMultiDim2(intTrial) = sqrt(sum((vecAct-mean(vecAct)).^2));
			
			%get pref pop
			%intTrialOri = [];
			vecPrefNeurons = vecNeuronPrefStim == vecStimOris(intTrial);
			vecActPrefPop(intTrial) = mean(vecAct(vecPrefNeurons));
			vecActWholePop(intTrial) = mean(vecAct);
			vecZActAll(intTrial) = mean(matRespNormPerNeuron(:,intTrial));
			vecActZPrefPop(intTrial) = mean(matRespNormPerNeuron(vecPrefNeurons,intTrial));
			
			matSelectR = matSelectAllR;
			matSelectR(:,:,intTrial) = triu(true(intNeurons),1);
			vecR(intTrial) = mean(matR(matSelectR));
			vecR_SD(intTrial) = std(matR(matSelectR));
		end
		
		%% RT dependence
		vecRTs = cellMultiSes{intPopulation}.structStim.vecTrialRespSecs - cellMultiSes{intPopulation}.structStim.SecsOn; %get RTs for hit trials and selected contrasts
		indSelectRespTrials = ~isnan(vecRTs);
		vecRT = vecRTs(indSelectRespTrials);
		[vecRTsSorted,vecRTSortIndex] = sort(vecRT,'ascend');
		
		%% check bimodality of miss trials
		%select trials
		intPopCounter = intPopCounter + 1;
		vecAllH{intPopCounter} = [];
		for intC=2:5
			indSelect = cellSelectContrasts{intC} & ~indSelectRespTrials;
			subplot(2,2,intC-1)
			vecAllH{intPopCounter} = [vecAllH{intPopCounter}; zscore(vecHeterogeneity(indSelect))];
			histx(vecHeterogeneity(indSelect))
		
		end
		figure
		histx(vecAllH{intPopCounter})
		continue
		
		%% 
		vecAllH_Agg = cell2mat(vecAllH')
		[y,vecX] = histx(vecAllH_Agg)
		
		%extend
		intN = numel(vecX);
		dblStep = vecX(2)-vecX(1);
		vecX = vecX(1)-dblStep*3:dblStep:vecX(end)+dblStep*3;
		y = [zeros(1,(numel(vecX)-intN)/2) y zeros(1,(numel(vecX)-intN)/2)];
		
		%fit single gauss
		p0Uni = [mean(vecAllH_Agg) std(vecAllH_Agg) max(y)];
		lb = [-inf eps eps];%[mu1,sigma1,peak1]
		ub = [inf inf inf];
		
		[vecParamsUni,resnorm,residual,exitflag] = curvefitfun(@getGaussian,p0Uni,vecX,y,lb,ub);

		vecGaussian = getGaussian(vecParamsUni,vecX)
		
		%fit double gauss
		p0 = [mean(vecAllH_Agg) std(vecAllH_Agg) max(y) mean(vecAllH_Agg) std(vecAllH_Agg) max(y)];
		lb = [-inf eps eps -inf eps eps];%[mu1,sigma1,peak1,mu2,sigma2,peak2]
		ub = [inf inf inf inf inf inf];
		
		options.MaxFunEvals = 1000;
		[vecParams,resnorm,residual,exitflag] = curvefitfun(@getBimodalGaussCF,p0,vecX,y,lb,ub,options);

		
		figure
		bar(vecX,y)
		hold on
		plot(vecX,getGaussian(vecParamsUni,vecX),'Color',[1 1 0])
		plot(vecX,getBimodalGaussCF(vecParams,vecX),'b')
		plot(vecX,normpdf(vecX,vecParams(1),vecParams(2))*vecParams(3),'g')
		plot(vecX,normpdf(vecX,vecParams(4),vecParams(5))*vecParams(6),'r')
		xlim([-4 4])
		
		
		
		[dip, p_value, xlow,xup]=HartigansDipSignifTest(vecAllH_Agg,1000)

		
		
		[dip, p_value, xlow,xup]=HartigansDipSignifTest([randn(1,10000)*2+1 randn(1,1000)*2+2],1000)

		
		%% perform correlation sliding window analysis
		%{
		vecScales = [1]*25.4;
		vecFreqs = wavscal2frq(vecScales,'cmor1-1',1);
		intScales=length(vecScales);
		intPairs = (intNeurons*(intNeurons-1))/2;
		intT = length(cellMultiSes{intPopulation}.neuron(1).dFoF);
		matSlidingCorr=zeros(intScales,intT,intPairs);
		intPairCounter = 0;
		fprintf('Performing sliding window correlation analysis, please wait...\n');
		for intNeuron1=1:(intNeurons-1)
			for intNeuron2=(intNeuron1+1):intNeurons
				intPairCounter = intPairCounter + 1;
				%if mod(intPairCounter,100) == 0,fprintf('Pair %d/%d [%s]\n',intPairCounter,intPairs,getTime);end
				%matSlidingCorr(:,:,intPairCounter) = semblance(1:intTrials,,,intUseWaves);
				
				y1 = cellMultiSes{intPopulation}.neuron(intNeuron1).dFoF;
				y2 = cellMultiSes{intPopulation}.neuron(intNeuron2).dFoF;
				
				y1=y1-mean(y1(:)); y2=y2-mean(y2(:));
				
				c1=wavcwt(y1,vecScales,'cmor1-1');
				c2=wavcwt(y2,vecScales,'cmor1-1');
				ctc=c1.*conj(c2);                  % Cross wavelet transform amplitude
				spt=atan2(imag(ctc),real(ctc));
				matSlidingCorr(:,:,intPairCounter) = cos(spt);
			end
		end
		%paper used:
		% Cooper, G.R.J., and Cowan, D.R., 2008.
		% Comparing Time Series using Wavelet Based Semblance Analysis
		% Computers & Geosciences v.34(2) p.95-102.
		
		%get mean+SD of sliding window correlations
		matMeanCorr = mean(matSlidingCorr,3);
		matSDCorr = std(matSlidingCorr,[],3);
		
		%transform to trials
		vecWavCorrMean = nan(1,intTrials);
		vecWavCorrSD = nan(1,intTrials);
		for intTrial=1:intTrials
			vecTrialFrames = cellMultiSes{intPopulation}.structStim.FrameOn(intTrial):(cellMultiSes{intPopulation}.structStim.FrameOn(intTrial)+vecTrialDur(intTrial));
			vecWavCorrMean(intTrial) = mean(mean(matMeanCorr(:,vecTrialFrames),2),1);
			vecWavCorrSD(intTrial) = mean(mean(matSDCorr(:,vecTrialFrames),2),1);
		end
		
		%mean sliding corr
		vecWCMHit = vecWavCorrMean(indSelectRespTrials);
		vecWCMRT = vecWCMHit(vecRTSortIndex);
		
		%sd sliding corr
		vecWCSHit = vecWavCorrSD(indSelectRespTrials);
		vecWCSRT = vecWCSHit(vecRTSortIndex);
		
		
		%plot
		if sParams.boolSavePlots
			figure
			subplot(2,2,1);
			scatter(vecRTsSorted,vecWCMRT,'kx')
			
			%perform regressions
			sStatsC=regstats(vecWCMRT,vecRTsSorted,'linear',{'beta','rsquare','tstat'});
			vecX = get(gca,'XLim');
			vecY = polyval(sStatsC.beta([2 1]),vecX);
			hold on
			plot(vecX,vecY,'r')
			hold off
			xlim([0 3])
			title(sprintf('Mean sliding window corr; Lin reg: slope=%.3f; p=%.3f; intercept=%.3f; p=%.3f; R^2=%.3f',...
				sStatsC.beta(2),sStatsC.tstat.pval(2),sStatsC.beta(1),sStatsC.tstat.pval(1),sStatsC.rsquare))
			
			xlabel('Reaction Time (s)')
			ylabel('Mean sliding window corr')
			
			subplot(2,2,2);
			scatter(vecRTsSorted,vecWCSRT,'kx')
			
			%perform regressions
			sStatsC=regstats(vecWCSRT,vecRTsSorted,'linear',{'beta','rsquare','tstat'});
			vecX = get(gca,'XLim');
			vecY = polyval(sStatsC.beta([2 1]),vecX);
			hold on
			plot(vecX,vecY,'r')
			hold off
			xlim([0 3])
			title(sprintf('SD sliding window corr; Lin reg: slope=%.3f; p=%.3f; intercept=%.3f; p=%.3f; R^2=%.3f',...
				sStatsC.beta(2),sStatsC.tstat.pval(2),sStatsC.beta(1),sStatsC.tstat.pval(1),sStatsC.rsquare))
			
			xlabel('Reaction Time (s)')
			ylabel('SD sliding window corr')
			
			subplot(2,2,3);
			scatter(vecWCMRT,vecWCSRT,'kx')
			xlabel('Mean sliding window correlation')
			ylabel('SD sliding window correlation')
			
			drawnow;
			jFig = get(handle(gcf), 'JavaFrame');
			jFig.setMaximized(true);
			figure(gcf);
			drawnow;
			strFig = sprintf('%sagg_RT_dependence_semblance_pop%d_raw',strSes,intPopulation);
			export_fig([strFig '.tif']);
			export_fig([strFig '.pdf']);
		end
		
		%% other RT measures
		%heterogeneity
		vecHetHit = vecHeterogeneity(indSelectRespTrials);
		vecHetRT = vecHetHit(vecRTSortIndex);
		
		%Z-scored activity
		matActZ = matRespNormPerNeuron(:,indSelectRespTrials)';
		
		matZRTSorted = matActZ(vecRTSortIndex,:);
		vecZRT = mean(matZRTSorted,2);
		
		%dF/F
		matActdFoF = matTrialResponse(:,indSelectRespTrials)';
		
		matARTSorted = matActdFoF(vecRTSortIndex,:);
		vecART = mean(matARTSorted,2);
		
		%variance
		vecVRT = var(matARTSorted,[],2);
		
		%sparseness
		vecSparseness = kurtosis(matARTSorted,[],2)-3;
		
		%pearson-like similarity metric
		vecRHit = vecR(indSelectRespTrials);
		vecCorrLike = vecRHit(vecRTSortIndex);
		
		%spread in pearson-like similarity metric
		vecRSDHit = vecR_SD(indSelectRespTrials);
		vecCorrLikeSpread = vecRSDHit(vecRTSortIndex);
		
		%dF/F pref pop
		vecActPPHit = vecActPrefPop(indSelectRespTrials);
		vecAct_PP = vecActPPHit(vecRTSortIndex);
		
		%z-scored dF/F pref pop
		vecActZPPHit = vecActZPrefPop(indSelectRespTrials);
		vecActZ_PP = vecActZPPHit(vecRTSortIndex);
		
		
		%save
		cellSaveRTDependency{1,intPopulation} = vecRTsSorted';
		cellSaveRTDependency{2,intPopulation} = vecHetRT;
		cellSaveRTDependency{3,intPopulation} = vecZRT;
		cellSaveRTDependency{4,intPopulation} = vecART;
		cellSaveRTDependency{5,intPopulation} = vecVRT;
		cellSaveRTDependency{6,intPopulation} = vecSparseness;
		cellSaveRTDependency{7,intPopulation} = vecCorrLike;
		cellSaveRTDependency{8,intPopulation} = vecCorrLikeSpread;
		cellSaveRTDependency{9,intPopulation} = vecAct_PP;
		cellSaveRTDependency{10,intPopulation} = vecActZ_PP;
		cellSaveRTDependency{11,intPopulation} = vecWCMRT'; %mean
		cellSaveRTDependency{12,intPopulation} = vecWCSRT'; %sd
		
		%plot
		if sParams.boolSavePlots
			figure
			subplot(3,3,1);
			scatter(vecRTsSorted,vecHetRT,'kx')
			
			%perform regressions
			sStatsC=regstats(vecHetRT,vecRTsSorted,'linear',{'beta','rsquare','tstat'});
			vecX = get(gca,'XLim');
			vecY = polyval(sStatsC.beta([2 1]),vecX);
			hold on
			plot(vecX,vecY,'r')
			hold off
			xlim([0 3])
			title(sprintf('Pop heterogeneity; Lin reg: slope=%.3f; p=%.3f; intercept=%.3f; p=%.3f; R^2=%.3f',...
				sStatsC.beta(2),sStatsC.tstat.pval(2),sStatsC.beta(1),sStatsC.tstat.pval(1),sStatsC.rsquare))
			
			xlabel('Reaction Time (s)')
			ylabel('Mean population activity dissimilarity')
			
			%z-scored activity
			subplot(3,3,2);
			scatter(vecRTsSorted,vecZRT,'kx')
			
			%perform regressions
			sStatsC=regstats(vecZRT,vecRTsSorted,'linear',{'beta','rsquare','tstat'});
			vecX = get(gca,'XLim');
			vecY = polyval(sStatsC.beta([2 1]),vecX);
			hold on
			plot(vecX,vecY,'r')
			hold off
			xlim([0 3])
			title(sprintf('Z-scored act; Lin reg: slope=%.3f; p=%.3f; intercept=%.3f; p=%.3f; R^2=%.3f',...
				sStatsC.beta(2),sStatsC.tstat.pval(2),sStatsC.beta(1),sStatsC.tstat.pval(1),sStatsC.rsquare))
			
			xlabel('Reaction Time (s)')
			ylabel('Mean z-scored population activity')
			
			%dF/F activity
			subplot(3,3,3);
			scatter(vecRTsSorted,vecART,'kx')
			
			%perform regressions
			sStatsC=regstats(vecART,vecRTsSorted,'linear',{'beta','rsquare','tstat'});
			vecX = get(gca,'XLim');
			vecY = polyval(sStatsC.beta([2 1]),vecX);
			hold on
			plot(vecX,vecY,'r')
			hold off
			xlim([0 3])
			title(sprintf('dF/F0; Lin reg: slope=%.3f; p=%.3f; intercept=%.3f; p=%.3f; R^2=%.3f',...
				sStatsC.beta(2),sStatsC.tstat.pval(2),sStatsC.beta(1),sStatsC.tstat.pval(1),sStatsC.rsquare))
			
			xlabel('Reaction Time (s)')
			ylabel('Mean dF/F population activity')
			
			%dF/F activity
			subplot(3,3,4);
			scatter(vecRTsSorted,vecVRT,'kx')
			
			%perform regressions
			sStatsC=regstats(vecVRT,vecRTsSorted,'linear',{'beta','rsquare','tstat'});
			vecX = get(gca,'XLim');
			vecY = polyval(sStatsC.beta([2 1]),vecX);
			hold on
			plot(vecX,vecY,'r')
			hold off
			xlim([0 3])
			title(sprintf('Variance; Lin reg: slope=%.3f; p=%.3f; intercept=%.3f; p=%.3f; R^2=%.3f',...
				sStatsC.beta(2),sStatsC.tstat.pval(2),sStatsC.beta(1),sStatsC.tstat.pval(1),sStatsC.rsquare))
			
			xlabel('Reaction Time (s)')
			ylabel('Variance of population activity')
			
			
			% sparseness
			subplot(3,3,5);
			scatter(vecRTsSorted,vecSparseness,'kx')
			
			%perform regressions
			sStatsC=regstats(vecSparseness,vecRTsSorted,'linear',{'beta','rsquare','tstat'});
			vecX = get(gca,'XLim');
			vecY = polyval(sStatsC.beta([2 1]),vecX);
			hold on
			plot(vecX,vecY,'r')
			hold off
			xlim([0 3])
			title(sprintf('Sparseness; Lin reg: slope=%.3f; p=%.3f; intercept=%.3f; p=%.3f; R^2=%.3f',...
				sStatsC.beta(2),sStatsC.tstat.pval(2),sStatsC.beta(1),sStatsC.tstat.pval(1),sStatsC.rsquare))
			
			xlabel('Reaction Time (s)')
			ylabel('Population Sparseness')
			
			% corr-like
			subplot(3,3,6);
			scatter(vecRTsSorted,vecCorrLike,'kx')
			
			%perform regressions
			sStatsC=regstats(vecCorrLike,vecRTsSorted,'linear',{'beta','rsquare','tstat'});
			vecX = get(gca,'XLim');
			vecY = polyval(sStatsC.beta([2 1]),vecX);
			hold on
			plot(vecX,vecY,'r')
			hold off
			xlim([0 3])
			title(sprintf('Pearson-like; Lin reg: slope=%.3f; p=%.3f; intercept=%.3f; p=%.3f; R^2=%.3f',...
				sStatsC.beta(2),sStatsC.tstat.pval(2),sStatsC.beta(1),sStatsC.tstat.pval(1),sStatsC.rsquare))
			
			xlabel('Reaction Time (s)')
			ylabel('Pearson-like')
			
			% corr-like spread
			subplot(3,3,7);
			scatter(vecRTsSorted,vecCorrLikeSpread,'kx')
			
			%perform regressions
			sStatsC=regstats(vecCorrLikeSpread,vecRTsSorted,'linear',{'beta','rsquare','tstat'});
			vecX = get(gca,'XLim');
			vecY = polyval(sStatsC.beta([2 1]),vecX);
			hold on
			plot(vecX,vecY,'r')
			hold off
			xlim([0 3])
			title(sprintf('CorrLikeSpread; Lin reg: slope=%.3f; p=%.3f; intercept=%.3f; p=%.3f; R^2=%.3f',...
				sStatsC.beta(2),sStatsC.tstat.pval(2),sStatsC.beta(1),sStatsC.tstat.pval(1),sStatsC.rsquare))
			
			xlabel('Reaction Time (s)')
			ylabel('Pearson-like sd')
			
			% pref pop dF/F
			subplot(3,3,8);
			scatter(vecRTsSorted,vecAct_PP,'kx')
			
			%perform regressions
			sStatsC=regstats(vecAct_PP,vecRTsSorted,'linear',{'beta','rsquare','tstat'});
			vecX = get(gca,'XLim');
			vecY = polyval(sStatsC.beta([2 1]),vecX);
			hold on
			plot(vecX,vecY,'r')
			hold off
			xlim([0 3])
			title(sprintf('PrefPop Act; Lin reg: slope=%.3f; p=%.3f; intercept=%.3f; p=%.3f; R^2=%.3f',...
				sStatsC.beta(2),sStatsC.tstat.pval(2),sStatsC.beta(1),sStatsC.tstat.pval(1),sStatsC.rsquare))
			
			xlabel('Reaction Time (s)')
			ylabel('dF/F pref-pop')
			
			% pref pop actZ
			subplot(3,3,9);
			scatter(vecRTsSorted,vecActZ_PP,'kx')
			
			%perform regressions
			sStatsC=regstats(vecActZ_PP,vecRTsSorted,'linear',{'beta','rsquare','tstat'});
			vecX = get(gca,'XLim');
			vecY = polyval(sStatsC.beta([2 1]),vecX);
			hold on
			plot(vecX,vecY,'r')
			hold off
			xlim([0 3])
			title(sprintf('PrefPop ActZ; Lin reg: slope=%.3f; p=%.3f; intercept=%.3f; p=%.3f; R^2=%.3f',...
				sStatsC.beta(2),sStatsC.tstat.pval(2),sStatsC.beta(1),sStatsC.tstat.pval(1),sStatsC.rsquare))
			
			xlabel('Reaction Time (s)')
			ylabel('Z-scored dF/F pref-pop')
			
			
			drawnow;
			jFig = get(handle(gcf), 'JavaFrame');
			jFig.setMaximized(true);
			figure(gcf);
			drawnow;
			strFig = sprintf('%sagg_RT_dependence_pop%d_raw',strSes,intPopulation);
			export_fig([strFig '.tif']);
			export_fig([strFig '.pdf']);
		end
		
		%}
		
		%% analyze orthogonal vs parallel properties
		figure
		subplot(2,2,1)
		%scatter3(vecZActAll,vecHeterogeneity,vecHetMultiDim)
		%xlabel('Z-scored dF/F0');
		%ylabel('Heterogeneity');
		%zlabel('Multi-dimensional heterogeneity');
		
		vecMeanActResp = nan(1,6);
		vecSEActResp = nan(1,6);
		vecMeanHResp = nan(1,6);
		vecSEHResp = nan(1,6);
		vecMeanActNoResp = nan(1,6);
		vecSEActNoResp = nan(1,6);
		vecMeanHNoResp = nan(1,6);
		vecSEHNoResp = nan(1,6);
		vecHitMissSeparabilityH = nan(1,6);
		vecHitMissSeparabilityA = nan(1,6);
		for intC=1:6
			vecHet = vecHeterogeneity;%vecHetMultiDim;%vecHeterogeneity
			vecAct = vecActWholePop;%vecZActAll;%vecActPrefPop
			subplot(2,3,intC)
			indC = cellSelectContrasts{intC};
			h= scatter(vecAct(indC)',vecHet(indC)','CDATA',1-double(indSelectRespTrials(indC)));
			colormap(redgreen);
			xlabel('Mean activity (dF/F0)');
			ylabel('Heterogeneity');
			
			%get data
			vecRespA = vecAct(indC & indSelectRespTrials);
			vecNoRespA = vecAct(indC & ~indSelectRespTrials);
			vecRespH = vecHet(indC & indSelectRespTrials);
			vecNoRespH = vecHet(indC & ~indSelectRespTrials);
			
			%mean+sem
			vecMeanActResp(intC) = mean(vecRespA);
			vecSEActResp(intC) = std(vecRespA)/sqrt(length(vecRespA));
			vecMeanHResp(intC) =  mean(vecRespH);
			vecSEHResp(intC) = std(vecRespH)/sqrt(length(vecRespH));
			vecMeanActNoResp(intC) =  mean(vecNoRespA);
			vecSEActNoResp(intC) = std(vecNoRespA)/sqrt(length(vecNoRespA));
			vecMeanHNoResp(intC) =  mean(vecNoRespH);
			vecSEHNoResp(intC) = std(vecNoRespH)/sqrt(length(vecNoRespH));
			
			%separability
			indTheseTrialsResp = indSelectRespTrials(indC);
			
			%het
			vecTheseTrialsHet =  vecHet(indC);
			if numel(unique(indTheseTrialsResp))==1
				AUC=1;
				optrocpt=[1 1];
				h= scatter(vecAct(indC)',vecHet(indC)');
			else
				[X,Y,T,AUC,optrocpt] = perfcurve(indTheseTrialsResp',vecTheseTrialsHet,true);
				if intC==3
					XH = X;
					YH = Y;
				end
			end
			%vecHitMissSeparabilityH(intC) = abs(0.5-AUC)+0.5;%0.5+(optrocpt(2)-optrocpt(1))/2;%AUC
			vecHitMissSeparabilityH(intC) = AUC;%0.5+(optrocpt(2)-optrocpt(1))/2;%AUC
			
			%df/f
			vecTheseTrialsAct =  vecAct(indC);
			if numel(unique(indTheseTrialsResp))==1,AUC=1;else
				[X,Y,T,AUC,optrocpt] = perfcurve(indTheseTrialsResp',vecTheseTrialsAct,true);
			end
			%vecHitMissSeparabilityA(intC) = abs(0.5-AUC)+0.5;%0.5+(optrocpt(2)-optrocpt(1))/2;%AUC
			vecHitMissSeparabilityA(intC) = AUC;%0.5+(optrocpt(2)-optrocpt(1))/2;%AUC
			
			title(sprintf('Separability; Act: %.3f; Het: %.3f',vecHitMissSeparabilityA(intC),vecHitMissSeparabilityH(intC)));
		end
		%continue
		cellSaveSeparability(intPopulation,1) = mean(vecHitMissSeparabilityH(2:5));
		cellSaveSeparability(intPopulation,2) = mean(vecHitMissSeparabilityA(2:5));
		
		%{
		figure
		subplot(2,2,1)
		vecPlotC=[0.2 0.5 2 8 32 100];
		errorfill(vecPlotC,vecMeanActNoResp,vecSEActNoResp,[1 0 0],[1 0.5 0.5]);
		hold on
		errorfill(vecPlotC,vecMeanActResp,vecSEActResp,[0 1 0],[0.5 1 0.5]);
		hold off
		ylabel('dF/F0')
		xlabel('Stimulus contrast (%)')
		set(gca,'xscale','log');
		
		subplot(2,2,2)
		vecPlotC=[0.2 0.5 2 8 32 100];
		errorfill(vecPlotC,vecMeanHNoResp,vecSEHNoResp,[1 0 0],[1 0.5 0.5]);
		hold on
		errorfill(vecPlotC,vecMeanHResp,vecSEHResp,[0 1 0],[0.5 1 0.5]);
		hold off
		ylabel('Heterogeneity')
		xlabel('Stimulus contrast (%)')
		set(gca,'xscale','log');
		%}
		%{
		%% calculate mean hit/miss for different measures and save to output variable (in Cohen's d)
		matCohensD_HitMiss = nan(11,6);%[measures x contrasts]
		for intC=1:6
			indC = cellSelectContrasts{intC};
			%heterogeneity
			if sum(indSelectRespTrials&indC) > 1 && sum(~indSelectRespTrials&indC) > 1
				vecHetHit = vecHeterogeneity(indSelectRespTrials&indC);
				vecHetMiss = vecHeterogeneity(~indSelectRespTrials&indC);
				dblCohensD_Het = getCohensD(vecHetHit,vecHetMiss);
				matCohensD_HitMiss(1,intC) = dblCohensD_Het;
				
				%z-scored dF/F
				vecZActWholePopHit = mean(matRespNormPerNeuron(:,indSelectRespTrials&indC),1);
				vecZActWholePopMiss = mean(matRespNormPerNeuron(:,~indSelectRespTrials&indC),1);
				dblCohensD_WPZ = getCohensD(vecZActWholePopHit,vecZActWholePopMiss);
				matCohensD_HitMiss(2,intC) = dblCohensD_WPZ;
				
				%dF/F
				vecActWholePopHit = mean(matTrialResponse(:,indSelectRespTrials&indC),1);
				vecActWholePopMiss = mean(matTrialResponse(:,~indSelectRespTrials&indC),1);
				dblCohensD_WPA = getCohensD(vecActWholePopHit,vecActWholePopMiss);
				matCohensD_HitMiss(3,intC) = dblCohensD_WPA;
				
				%variance
				vecVarHit = var(matTrialResponse(:,indSelectRespTrials&indC),[],1);
				vecVarMiss = var(matTrialResponse(:,~indSelectRespTrials&indC),[],1);
				dblCohensD_Var = getCohensD(vecVarHit,vecVarMiss);
				matCohensD_HitMiss(4,intC) = dblCohensD_Var;
				
				%sparseness
				vecSparsenessHit = kurtosis(matTrialResponse(:,indSelectRespTrials&indC),[],1)-3;
				vecSparsenessMiss = kurtosis(matTrialResponse(:,~indSelectRespTrials&indC),[],1)-3;
				dblCohensD_Sparse = getCohensD(vecSparsenessHit,vecSparsenessMiss);
				matCohensD_HitMiss(5,intC) = dblCohensD_Sparse;
				
				%pearson-like similarity metric
				vecRHit = vecR(indSelectRespTrials&indC);
				vecRMiss = vecR(~indSelectRespTrials&indC);
				dblCohensD_R = getCohensD(vecRHit,vecRMiss);
				matCohensD_HitMiss(6,intC) = dblCohensD_R;
				
				%spread in pearson-like similarity metric
				vecRSDHit = vecR_SD(indSelectRespTrials&indC);
				vecRSDMiss = vecR_SD(~indSelectRespTrials&indC);
				dblCohensD_HRSD = getCohensD(vecRSDHit,vecRSDMiss);
				matCohensD_HitMiss(7,intC) = dblCohensD_HRSD;
				
				%dF/F pref pop
				vecActPPHit = vecActPrefPop(indSelectRespTrials&indC);
				vecActPPMiss = vecActPrefPop(~indSelectRespTrials&indC);
				dblCohensD_APP = getCohensD(vecActPPHit,vecActPPMiss);
				matCohensD_HitMiss(8,intC) = dblCohensD_APP;
				
				%z-scored dF/F pref pop
				vecActZPPHit = vecActZPrefPop(indSelectRespTrials&indC);
				vecActZPPMiss = vecActZPrefPop(~indSelectRespTrials&indC);
				dblCohensD_ZPP = getCohensD(vecActZPPHit,vecActZPPMiss);
				matCohensD_HitMiss(9,intC) = dblCohensD_ZPP;
				
				%mean sliding corr
				vecWCMHit = vecWavCorrMean(indSelectRespTrials&indC);
				vecWCMMiss = vecWavCorrMean(~indSelectRespTrials&indC);
				dblCohensD_WCM = getCohensD(vecWCMHit,vecWCMMiss);
				matCohensD_HitMiss(10,intC) = dblCohensD_WCM;
				
				%sd sliding corr
				vecWCSHit = vecWavCorrSD(indSelectRespTrials&indC);
				vecWCSMiss = vecWavCorrSD(~indSelectRespTrials&indC);
				dblCohensD_WCS = getCohensD(vecWCSHit,vecWCSMiss);
				matCohensD_HitMiss(11,intC) = dblCohensD_WCS;
			end
		end
		vecMeanD = nanmean(matCohensD_HitMiss(:,2:5),2);
		vecSD_D = nanstd(matCohensD_HitMiss(:,2:5),[],2);
		
		%save data
		cellSaveCohensD_HitMiss{intPopulation} = matCohensD_HitMiss;
		%}
		%{
		%% calculate inter-trial correlations
		%loop through contrasts & stim types for trial reordering
		intStimType = 0;
		vecStimType = nan(1,length(cellSelectOri{1}));
		for intStim=unique(vecNeuronPrefStim)
			for intContrastIndex=2:length(cellSelectContrasts)
				intStimType = intStimType + 1;
				%trials
				indSelectTrials = cellSelectOri{intStim} & cellSelectContrasts{intContrastIndex};
				vecStimType(indSelectTrials) = intStimType;
			end
		end
		vecOriType = mod(vecStimType,4);
		vecOriType(vecOriType==0) = 4;
		
		%source matrices
		%matTrialResponse(intNeuron,intTrial);
		%matRespNormPerContrastintNeuron,intTrial);
		%normalize per contrast
		matRespNormPerContrast = nan(size(matTrialResponse));
		for intContrastIndex=1:length(cellSelectContrasts)
			vecSelectContrastTrials = cellSelectContrasts{intContrastIndex};
			matRespNormPerContrast(:,vecSelectContrastTrials) = zscore(matTrialResponse(:,vecSelectContrastTrials),[],2);
		end
		
		%for whole pop or only preferred pop
		if sParams.boolSavePlots,hPopCorr = figure;end
		for intPopType=1:2
			if intPopType == 1
				strPopType = 'PrefNeurons';
			elseif intPopType == 2
				strPopType = 'NonPrefNeurons';
			end
			
			%% get real data
			for intPrefType=1:4
				indSelect = vecNeuronPrefStim' == intPrefType;
				if intPopType == 2
					indSelect = ~indSelect;
				end
				%strPopType
				%sum(indSelect)
				
				%calculate trial-by-trial correlations of population activity
				if sum(indSelect) > 0
					cellCorr{intPrefType}=corr(matTrialResponse(indSelect,:));
				else
					cellCorr{intPrefType} = nan(intTrials,intTrials);
				end
			end
			matCorr = nan(size(cellCorr{1}));
			for intOriType=vecOriType
				if isnan(intOriType),continue;end
				%get trials
				indSelectTrials = vecOriType == intOriType;
				%assign
				matCorr(indSelectTrials,indSelectTrials) = cellCorr{intPrefType}(indSelectTrials,indSelectTrials);
			end
			
			%% get shuffled data
			intIters = 100;
			for intIter=1:intIters
				if mod(intIter,10) == 0,fprintf('Shuffling iter %d/%d [%s]\n',intIter,intIters,getTime);end
				
				%prep variables
				cellCorr = [];
				
				%shuffle data
				matTrialRespShuffled = matTrialResponse;
				for intNeuron=1:intNeurons
					vecAllTrials = 1:length(vecOriType);
					for intOriType=vecOriType
						if isnan(intOriType),continue;end
						%get trials
						indSelectTrials = vecOriType == intOriType;
						
						%assign
						vecTheseOriTrials=find(indSelectTrials);
						
						vecAllTrials(indSelectTrials) = vecTheseOriTrials(randperm(sum(indSelectTrials)));
					end
					matTrialRespShuffled(intNeuron,:) = matTrialResponse(intNeuron,vecAllTrials);
				end
				
				%get correlations
				for intPrefType=1:4
					indSelect = vecNeuronPrefStim' == intPrefType;
					if intPopType == 2
						indSelect = ~indSelect;
					end
					%strPopType
					%sum(indSelect)
					
					%calculate trial-by-trial correlations of population activity
					if sum(indSelect) > 0
						cellCorr{intPrefType}=corr(matTrialRespShuffled(indSelect,:));
					else
						cellCorr{intPrefType} = nan(intTrials,intTrials);
					end
				end
				if intIter == 1
					matCorrShuffled = nan([size(cellCorr{1}) intIters]);
				end
				for intOriType=vecOriType
					if isnan(intOriType),continue;end
					%get trials
					indSelectTrials = vecOriType == intOriType;
					%assign
					matCorrShuffled(indSelectTrials,indSelectTrials,intIter) = cellCorr{intPrefType}(indSelectTrials,indSelectTrials);
				end
			end
			matCorrMeanShuffled = mean(matCorrShuffled,3);
			
			%% get means
			%compute absolute trial distance
			matTrialDist = abs(repmat(1:intTrials,[intTrials 1])-repmat((1:intTrials)',[1 intTrials]));
			%matSelect = tril(true(size(matTrialDist)),-1);
			indMiss = cellMultiSes{intPopulation}.structStim.vecTrialResponse == 0;
			indHit = cellMultiSes{intPopulation}.structStim.vecTrialResponse == 1;
			indFast = cellMultiSes{intPopulation}.structStim.vecTrialRespSecs <= nanmedian(cellMultiSes{intPopulation}.structStim.vecTrialRespSecs);
			indSlow = cellMultiSes{intPopulation}.structStim.vecTrialRespSecs > nanmedian(cellMultiSes{intPopulation}.structStim.vecTrialRespSecs);
			
			%separate across-trial correlation correlations for different trial types (1-4)
			%also split for hit/miss trials
			intHitTrials = sum(indHit);
			intMissTrials = sum(indMiss);
			intFastTrials = sum(indFast);
			intSlowTrials = sum(indSlow);
			vecCorrHit = [];
			vecCorrMiss = [];
			vecCorrFast = [];
			vecCorrSlow = [];
			vecCorrHitShuffled = [];
			vecCorrMissShuffled = [];
			vecCorrFastShuffled = [];
			vecCorrSlowShuffled = [];
			for intOriType=vecOriType
				%% hit; select trials
				if isnan(intOriType),continue;end
				indSelectHitTrials = vecOriType == intOriType & indHit;
				if sum(indSelectHitTrials) > 2
					%get trials
					matCorrSubHit = matCorr(indSelectHitTrials,indSelectHitTrials);
					matCorrSubHitShuffled = matCorrMeanShuffled(indSelectHitTrials,indSelectHitTrials);
					
					%subselection matrix hit
					matSubSelectHit = tril(true(size(matCorrSubHit)),-1);
					
					%subselect
					vecCorrHit = [vecCorrHit matCorrSubHit(matSubSelectHit)'];
					vecCorrHitShuffled = [vecCorrHitShuffled matCorrSubHitShuffled(matSubSelectHit)'];
				end
				
				%% miss; select trials
				indSelectMissTrials = vecOriType == intOriType & indMiss;
				if sum(indSelectMissTrials) > 2
					%get trials
					matCorrSubMiss = matCorr(indSelectMissTrials,indSelectMissTrials);
					matCorrSubMissShuffled = matCorrMeanShuffled(indSelectMissTrials,indSelectMissTrials);
					
					%subselection matrix hit
					matSubSelectMiss = tril(true(size(matCorrSubMiss)),-1);
					
					%subselect
					vecCorrMiss = [vecCorrMiss matCorrSubMiss(matSubSelectMiss)'];
					vecCorrMissShuffled = [vecCorrMissShuffled matCorrSubMissShuffled(matSubSelectMiss)'];
				end
				
				%% fast; select trials
				indSelectFastTrials = vecOriType == intOriType & indFast;
				if sum(indSelectFastTrials) > 2
					%get trials
					matCorrSubFast = matCorr(indSelectFastTrials,indSelectFastTrials);
					matCorrSubFastShuffled = matCorrMeanShuffled(indSelectFastTrials,indSelectFastTrials);
					
					%subselection matrix hit
					matSubSelectFast = tril(true(size(matCorrSubFast)),-1);
					
					%subselect
					vecCorrFast = [vecCorrFast matCorrSubFast(matSubSelectFast)'];
					vecCorrFastShuffled = [vecCorrFastShuffled matCorrSubFastShuffled(matSubSelectFast)'];
				end
				
				%% slow; select trials
				indSelectSlowTrials = vecOriType == intOriType & indSlow;
				if sum(indSelectSlowTrials) > 2
					%get trials
					matCorrSubSlow = matCorr(indSelectSlowTrials,indSelectSlowTrials);
					matCorrSubSlowShuffled = matCorrMeanShuffled(indSelectSlowTrials,indSelectSlowTrials);
					
					%subselection matrix hit
					matSubSelectSlow = tril(true(size(matCorrSubSlow)),-1);
					
					%subselect
					vecCorrSlow = [vecCorrSlow matCorrSubSlow(matSubSelectSlow)'];
					vecCorrSlowShuffled = [vecCorrSlowShuffled matCorrSubSlowShuffled(matSubSelectSlow)'];
				end
			end
			dblMeanCorrHit = mean(vecCorrHit);
			dblMeanCorrMiss = mean(vecCorrMiss);
			dblMeanCorrFast = mean(vecCorrFast);
			dblMeanCorrSlow = mean(vecCorrSlow);
			dblErrCorrHit = std(vecCorrHit)/sqrt(length(vecCorrHit));
			dblErrCorrMiss = mean(vecCorrMiss)/sqrt(length(vecCorrMiss));
			dblErrCorrFast = mean(vecCorrFast)/sqrt(length(vecCorrFast));
			dblErrCorrSlow = mean(vecCorrSlow)/sqrt(length(vecCorrSlow));
			dblShuffledMeanCorrHit = mean(vecCorrHitShuffled);
			dblShuffledMeanCorrMiss = mean(vecCorrMissShuffled);
			dblShuffledMeanCorrFast = mean(vecCorrFastShuffled);
			dblShuffledMeanCorrSlow = mean(vecCorrSlowShuffled);
			dblShuffledErrCorrHit = std(vecCorrHitShuffled)/sqrt(length(vecCorrHitShuffled));
			dblShuffledErrCorrMiss = std(vecCorrMissShuffled)/sqrt(length(vecCorrMissShuffled));
			dblShuffledErrCorrFast = std(vecCorrFastShuffled)/sqrt(length(vecCorrFastShuffled));
			dblShuffledErrCorrSlow = std(vecCorrSlowShuffled)/sqrt(length(vecCorrSlowShuffled));
			
			if intPopType==1 %pref
				vecMeanCorrPref = [dblMeanCorrMiss dblMeanCorrSlow dblMeanCorrFast];
				vecMeanCorrPrefShuffled = [dblShuffledMeanCorrMiss dblShuffledMeanCorrSlow dblShuffledMeanCorrFast];
			else
				vecMeanCorrNonPref = [dblMeanCorrMiss dblMeanCorrSlow dblMeanCorrFast];
				vecMeanCorrNonPrefShuffled = [dblShuffledMeanCorrMiss dblShuffledMeanCorrSlow dblShuffledMeanCorrFast];
			end
			[h,pHM] = ttest2(vecCorrHit,vecCorrMiss);
			[h,pMF] = ttest2(vecCorrMiss,vecCorrFast);
			[h,pMS] = ttest2(vecCorrMiss,vecCorrSlow);
			[h,pFS] = ttest2(vecCorrFast,vecCorrSlow);
			
			%plot
			if sParams.boolSavePlots
				subplot(2,2,1+(intPopType-1)*2)
				matPlotCorr = matCorr(~isnan(vecOriType),~isnan(vecOriType));
				matPlotCorr(diag(true(1,length(matPlotCorr)))) = nan;
				vecLim = [-1 1];%[-nanmax(abs(matPlotCorr(:))) nanmax(abs(matPlotCorr(:)))];
				imagesc(matPlotCorr,vecLim);
				matC = redbluepurple(128);
				nancolorbar(matPlotCorr,vecLim,matC);
				title(sprintf('Inter-trial correlation of %s population activity (block %d)',strPopType,intPopulation))
				xlabel('Trial')
				ylabel('Trial')
				
				subplot(2,2,2+(intPopType-1)*2)
				
				errorbar(1:4,[dblMeanCorrMiss dblMeanCorrHit dblMeanCorrFast dblMeanCorrSlow],[dblErrCorrMiss dblErrCorrHit dblErrCorrFast dblErrCorrSlow],'Linestyle','none','Marker','x','Color','b');
				hold on
				errorbar(1:4,[dblShuffledMeanCorrMiss dblShuffledMeanCorrHit dblShuffledMeanCorrFast dblShuffledMeanCorrSlow],[dblShuffledErrCorrMiss dblShuffledErrCorrHit dblShuffledErrCorrFast dblShuffledErrCorrSlow],'Linestyle','none','Marker','x','Color','k');
				hold off
				set(gca,'XTick',1:4,'XTickLabel',{'Miss','Hit-All','Hit-Fast','Hit-Slow'})
				title(sprintf('Mean inter-trial correlation of %s population activity\nT-tests; Hit-miss,p=%.3f; Miss-Fast,p=%.3f; Miss-Slow,p=%.3f; Fast-Slow,p=%.3f',strPopType,pHM,pMF,pMS,pFS))
				ylabel('Pearson correlation of population activity')
			end
		end
		if sParams.boolSavePlots
			drawnow;
			jFig = get(handle(gcf), 'JavaFrame');
			jFig.setMaximized(true);
			figure(gcf);
			drawnow;
			strFig = sprintf('%sintertrialcorr_pop%d_raw',strSes,intPopulation);
			export_fig([strFig '.tif']);
			export_fig([strFig '.pdf']);
		end
		%%%% OUTPUT
		cellSaveITC{intPopulation,1} = [vecMeanCorrNonPref;vecMeanCorrPref];
		cellSaveITC{intPopulation,2} = [vecMeanCorrNonPrefShuffled;vecMeanCorrPrefShuffled]; %shuffled
		
		%% test for multidimensional non-uniformity
		%calculate distribution of pairwise inter-point distances in
		%multidimensional space
		%}
		%type (norm or raw)
		intType = 2;
		indTestTrials = vecStimContrasts ~= 1 & vecStimContrasts ~= 6;
		indTestRespTrials = indSelectRespTrials(indTestTrials);
		if intType == 1
			matTestResps = matTrialResponse(:,indTestTrials);
			%binning vector
			dblStep = 0.02; %dF/F: 0.02
			vecBins = 0:dblStep:1.4;
			strLabelX = 'Distance (dF/F0)';
		else
			matTestResps = matRespNormPerNeuron(:,indTestTrials);
			%binning vector
			dblStep = 0.4; %dF/F: 0.02
			vecBins = 0:dblStep:30;
			strLabelX = 'Distance (z-score units)';
		end
		intTestTrials = size(matTestResps,2);
		matPairwiseDistances = nan(intTestTrials,intTestTrials);
		for intTrial=1:intTestTrials
			matPairwiseDistances(intTrial,:) = sqrt(sum(bsxfun(@minus,matTestResps,matTestResps(:,intTrial)).^2,1));
		end
		matSelectTrials = tril(true(intTestTrials,intTestTrials),-1);
		matSelectHitTrials = matSelectTrials;
		matSelectHitTrials(~indTestRespTrials,~indTestRespTrials) = false;
		matSelectMissTrials = matSelectTrials;
		matSelectMissTrials(indTestRespTrials,indTestRespTrials) = false;
		
		
		%subplot(2,2,1)
		%imagesc(matPairwiseDistances);colormap(hot);colorbar;
		%title('Pairwise (trials) Euclidian distances');
		%xlabel('Trial')
		%ylabel('Trial')
		
		
		vecDistroAll = hist(matPairwiseDistances(matSelectTrials),vecBins);
		vecDistroAll = vecDistroAll/sum(vecDistroAll);
		vecDistroHit = hist(matPairwiseDistances(matSelectHitTrials),vecBins);
		vecDistroHit = vecDistroHit/sum(vecDistroHit);
		vecDistroMiss = hist(matPairwiseDistances(matSelectMissTrials),vecBins);
		vecDistroMiss = vecDistroMiss/sum(vecDistroMiss);
		
		if sParams.boolSavePlots
			figure
			subplot(2,2,1);
			stairs(vecBins,vecDistroAll,'Color','k');
			hold on
			stairs(vecBins,vecDistroHit,'Color','g');
			stairs(vecBins,vecDistroMiss,'Color','r');
			hold off
			legend({'all','hit','miss'});
			xlabel(strLabelX);
			ylabel('Normalized count of trial-pairs');
			title(sprintf('Mean distances, all=%.3f, hit=%.3f, miss=%.3f',mean(matPairwiseDistances(matSelectTrials)),mean(matPairwiseDistances(matSelectHitTrials)),mean(matPairwiseDistances(matSelectMissTrials))));
		end
		
		%for each point (i.e., trial), reflect/mirror that neuron through
		%the multidimensional minimal-heterogeneity axis (i.e., mean
		%response gradient)
		matMirrorDistroAll = nan(intTestTrials,length(vecBins));
		matMirrorDistroHit = nan(intTestTrials,length(vecBins));
		matMirrorDistroMiss = nan(intTestTrials,length(vecBins));
		for intTrial = 1:intTestTrials;
			vecPoint = matTestResps(:,intTrial);
			dblMean = mean(vecPoint);
			vecRel = vecPoint - dblMean;
			vecMirrored = -vecRel + dblMean;
			matMirrorResps = matTestResps;
			matMirrorResps(:,intTrial) = vecMirrored;
			
			%recalculate the distribution of pairwise inter-point distances in
			%multidimensional space
			matPairwiseDistancesMirrored = matPairwiseDistances;
			vecMirrorDistances = sqrt(sum(bsxfun(@minus,matMirrorResps,matMirrorResps(:,intTrial)).^2,1));
			matPairwiseDistancesMirrored(intTrial,:) = vecMirrorDistances;
			matPairwiseDistancesMirrored(:,intTrial) = vecMirrorDistances;
			
			vecMirrorDistroAll = hist(matPairwiseDistancesMirrored(matSelectTrials),vecBins);
			vecMirrorDistroAll = vecMirrorDistroAll/sum(vecMirrorDistroAll);
			vecMirrorDistroHit = hist(matPairwiseDistancesMirrored(matSelectHitTrials),vecBins);
			vecMirrorDistroHit = vecMirrorDistroHit/sum(vecMirrorDistroHit);
			vecMirrorDistroMiss = hist(matPairwiseDistancesMirrored(matSelectMissTrials),vecBins);
			vecMirrorDistroMiss = vecMirrorDistroMiss/sum(vecMirrorDistroMiss);
			
			matMirrorDistroAll(intTrial,:) = vecMirrorDistroAll-vecDistroAll;
			matMirrorDistroHit(intTrial,:) = vecMirrorDistroHit-vecDistroHit;
			matMirrorDistroMiss(intTrial,:) = vecMirrorDistroMiss-vecDistroMiss;
			
		end
		%compare this distribution to the original; if changes occur, the
		%population response is asymmetricaly distributed around the axis,
		%and specific multidimensional response clusters that correspond to
		%hit-trials are present. if no changes occur, the population
		%response is symmetric, and heterogeneity is the main factor
		%determining stimulus detection
		
		
		vecMirrorDistroAll = mean(matMirrorDistroAll,1);
		%vecMirrorDistroAll = vecMirrorDistroAll/sum(vecMirrorDistroAll);
		vecMirrorDistroHit = mean(matMirrorDistroHit,1);
		%vecMirrorDistroHit = vecMirrorDistroHit/sum(vecMirrorDistroHit);
		vecMirrorDistroMiss = mean(matMirrorDistroMiss,1);
		%vecMirrorDistroMiss = vecMirrorDistroMiss/sum(vecMirrorDistroMiss);
		
		if sParams.boolSavePlots
			subplot(2,2,2);
			%stairs(vecBins,vecDistroAll,'Color','k');
			hold on
			stairs(vecBins,vecMirrorDistroAll,'Color','b');
			hold off
			xlabel(strLabelX);
			ylabel('Normalized count of trial-pairs');
			title(sprintf('Mean distance All trials; Mirrored-Raw'));%, data=%.3f, mirrored=%.3f',mean(matPairwiseDistances(matSelectTrials)),mean(matPairwiseDistances(matSelectHitTrials)),mean(matPairwiseDistances(matSelectMissTrials))));
			
			subplot(2,2,3);
			%stairs(vecBins,vecDistroHit,'Color','k');
			hold on
			stairs(vecBins,vecMirrorDistroHit,'Color','b');
			hold off
			xlabel(strLabelX);
			ylabel('Normalized count of trial-pairs');
			title(sprintf('Mean distance Hit trials; Mirrored-Raw'));%, data=%.3f, mirrored=%.3f',mean(matPairwiseDistances(matSelectTrials)),mean(matPairwiseDistances(matSelectHitTrials)),mean(matPairwiseDistances(matSelectMissTrials))));
			
			subplot(2,2,4);
			%stairs(vecBins,vecDistroMiss,'Color','k');
			hold on
			stairs(vecBins,vecMirrorDistroMiss,'Color','b');
			hold off
			xlabel(strLabelX);
			ylabel('Normalized count of trial-pairs');
			title(sprintf('Mean distance Miss trials; Mirrored-Raw'));%, data=%.3f, mirrored=%.3f',mean(matPairwiseDistances(matSelectTrials)),mean(matPairwiseDistances(matSelectHitTrials)),mean(matPairwiseDistances(matSelectMissTrials))));
			
			drawnow;
			jFig = get(handle(gcf), 'JavaFrame');
			jFig.setMaximized(true);
			figure(gcf);
			drawnow;
			strFig = sprintf('%sdistance_multidimensional_neural_response_space_pop%d',strSes,intPopulation);
			export_fig([strFig '.tif']);
			export_fig([strFig '.pdf']);
		end
		
		%save data
		cellSaveMultiDimDistanceDistros(1,:,intPopulation) = vecDistroAll;
		cellSaveMultiDimDistanceDistros(2,:,intPopulation) = vecDistroHit;
		cellSaveMultiDimDistanceDistros(3,:,intPopulation) = vecDistroMiss;
		cellSaveMultiDimDistanceDistros(4,:,intPopulation) = vecMirrorDistroAll;
		cellSaveMultiDimDistanceDistros(5,:,intPopulation) = vecMirrorDistroHit;
		cellSaveMultiDimDistanceDistros(6,:,intPopulation) = vecMirrorDistroMiss;
		%}
		%% compare normal heterogeneity with multi-dim heterogeneity
		%% hit/miss
		matCohensD2_HitMiss = nan(4,6);%[measures x contrasts]
		for intC=1:6
			indC = cellSelectContrasts{intC};
			
			if sum(indSelectRespTrials&indC) > 1 && sum(~indSelectRespTrials&indC) > 1
				%dF/F
				vecActWholePopHit = mean(matTrialResponse(:,indSelectRespTrials&indC),1);
				vecActWholePopMiss = mean(matTrialResponse(:,~indSelectRespTrials&indC),1);
				dblCohensD_WPA = getCohensD(vecActWholePopHit,vecActWholePopMiss);
				matCohensD2_HitMiss(1,intC) = dblCohensD_WPA;
				
				%heterogeneity
				vecHetHit = vecHeterogeneity(indSelectRespTrials&indC);
				vecHetMiss = vecHeterogeneity(~indSelectRespTrials&indC);
				dblCohensD_Het = getCohensD(vecHetHit,vecHetMiss);
				matCohensD2_HitMiss(2,intC) = dblCohensD_Het;
				
				%heterogeneity MD
				vecHetMDHit = vecHetMultiDim(indSelectRespTrials&indC);
				vecHetMDMiss = vecHetMultiDim(~indSelectRespTrials&indC);
				dblCohensD_HetMD = getCohensD(vecHetMDHit,vecHetMDMiss);
				matCohensD2_HitMiss(3,intC) = dblCohensD_HetMD;
				
				%heterogeneity MD2
				vecHetMD2Hit = vecHetMultiDim2(indSelectRespTrials&indC);
				vecHetMD2Miss = vecHetMultiDim2(~indSelectRespTrials&indC);
				dblCohensD_HetMD2 = getCohensD(vecHetMD2Hit,vecHetMD2Miss);
				matCohensD2_HitMiss(4,intC) = dblCohensD_HetMD2;
				
			end
		end
		vecMeanD = nanmean(matCohensD2_HitMiss(:,2:5),2);
		vecSD_D = nanstd(matCohensD2_HitMiss(:,2:5),[],2);
		
		%save data
		cellSaveCohensD_HitMiss2{intPopulation} = matCohensD2_HitMiss;
		
		%plot
		if sParams.boolSavePlots
			figure
			errorbar(1:4,vecMeanD,vecSD_D/sqrt(4),'x');
			set(gca,'xtick',1:4,'xticklabel',{'dF/F0','Het','Het MD','Het MD non-z'});
			ylabel('Cohen''s D')
			xlabel('Measure')
			
			%save figure
			drawnow;
			strFig = sprintf('%sHet_types_CohensD_pop%d',strSes,intPopulation);
			export_fig([strFig '.tif']);
			export_fig([strFig '.pdf']);
		end
		
		
		%% RT
		%heterogeneity
		vecHetHit = vecHeterogeneity(indSelectRespTrials);
		vecHetRT = vecHetHit(vecRTSortIndex);
		
		%heterogeneity
		vecHetHitMD = vecHetMultiDim(indSelectRespTrials); %z-score based
		vecHetMDRT = vecHetHitMD(vecRTSortIndex);
		
		%heterogeneity
		vecHetHitMD2 = vecHetMultiDim2(indSelectRespTrials); %non z-score based
		vecHetMD2RT = vecHetHitMD2(vecRTSortIndex);
		
		%dF/F0
		vecActdFoF = mean(matTrialResponse(:,indSelectRespTrials),1)';
		vecART = vecActdFoF(vecRTSortIndex);
		
		%plot
		if sParams.boolSavePlots
			figure
			for intPlot=1:4
				if intPlot == 1
					vecData = vecART;
					strLabel = 'Mean dF/F0';
				elseif intPlot == 2
					vecData = vecHetRT;
					strLabel = 'Heterogeneity';
				elseif intPlot == 3
					vecData = vecHetMDRT;
					strLabel = 'Het_MD';
				elseif intPlot == 4
					vecData = vecHetMD2RT;
					strLabel = 'Het_MD2';
				end
				subplot(2,2,intPlot);
				scatter(vecRTsSorted,vecData,'kx')
				
				%perform regressions
				sStatsC=regstats(vecData,vecRTsSorted,'linear',{'beta','rsquare','tstat'});
				vecX = get(gca,'XLim');
				vecY = polyval(sStatsC.beta([2 1]),vecX);
				hold on
				plot(vecX,vecY,'r')
				hold off
				xlim([0 3])
				title(sprintf('%s; Lin reg: slope=%.3f; p=%.3f; intercept=%.3f; p=%.3f; R^2=%.3f',...
					strLabel,sStatsC.beta(2),sStatsC.tstat.pval(2),sStatsC.beta(1),sStatsC.tstat.pval(1),sStatsC.rsquare))
				
				xlabel('Reaction Time (s)')
				ylabel(strLabel)
			end
			%save figure
			drawnow;
			jFig = get(handle(gcf), 'JavaFrame');
			jFig.setMaximized(true);
			figure(gcf);
			drawnow;
			strFig = sprintf('%sHet_types_RT_dep_pop%d',strSes,intPopulation);
			export_fig([strFig '.tif']);
			export_fig([strFig '.pdf']);
		end
		
		%save
		cellSaveRTDependency2{1,intPopulation} = vecRTsSorted';
		cellSaveRTDependency2{2,intPopulation} = vecART;
		cellSaveRTDependency2{3,intPopulation} = vecHetRT;
		cellSaveRTDependency2{4,intPopulation} = vecHetMDRT;
		cellSaveRTDependency2{5,intPopulation} = vecHetMD2RT;
		
		%% remove mean, heterogeneity, none or both
		matTestRespsNorm = matRespNormPerNeuron(:,indTestTrials);
		matMeanRemovedTestResps = bsxfun(@minus,matTestRespsNorm,mean(matTestRespsNorm,1));
		matHeteroRemovedTestResps = nan(size(matTestRespsNorm));
		vecHetMultiDimHetRem = nan(1,size(matTestRespsNorm,2));
		for intTrial=1:size(matTestRespsNorm,2)
			%get data
			vecA = matTestRespsNorm(:,intTrial);
			
			%remove heterogeneity
			vecHetRem = vecA/sqrt(sum((vecA-mean(vecA)).^2));
			vecHetRem = vecHetRem + (mean(vecA)-mean(vecHetRem));
			matHeteroRemovedTestResps(:,intTrial) = vecHetRem;
			vecHetMultiDimHetRem(intTrial) = sqrt(sum((vecHetRem-mean(vecHetRem)).^2));
		end
		%remove both
		matHeteroMeanRemovedTestResps = bsxfun(@minus,matHeteroRemovedTestResps,mean(matHeteroRemovedTestResps,1));
		
		if sParams.boolSavePlots
			figure
			colormap(hot)
			subplot(2,2,1)
			imagesc(matTestRespsNorm)
			title('Original data')
			xlabel('Trial')
			ylabel('Neuron')
			
			subplot(2,2,2)
			imagesc(matMeanRemovedTestResps)
			title('Mean removed')
			xlabel('Trial')
			ylabel('Neuron')
			
			subplot(2,2,3)
			imagesc(matHeteroRemovedTestResps)
			title('Heterogeneity removed')
			xlabel('Trial')
			ylabel('Neuron')
			
			subplot(2,2,4)
			imagesc(matHeteroMeanRemovedTestResps)
			title('Mean+heterogeneity removed')
			xlabel('Trial')
			ylabel('Neuron')
			
			%save figure
			drawnow;
			jFig = get(handle(gcf), 'JavaFrame');
			jFig.setMaximized(true);
			figure(gcf);
			drawnow;
			strFig = sprintf('%sneural_responses_het_mean_removed_pop%d',strSes,intPopulation);
			export_fig([strFig '.tif']);
			export_fig([strFig '.pdf']);
		end
		
		%% decode hit/miss
		%raw data (mean+heterogeneity+multidim)
		[dblPerformanceRawML,vecDecodedIndexCV,matPosteriorProbabilityCV] = doCrossValidatedDecodingML(matTestRespsNorm',indTestRespTrials,1);
		[dummy,vecCI_Raw] = binofit(dblPerformanceRawML*intTestTrials,intTestTrials);
		
		%heterogeneity removed (mean+multidim)
		[dblPerformanceNoHetML,vecDecodedIndexCV,matPosteriorProbabilityCV] = doCrossValidatedDecodingML(matHeteroRemovedTestResps',indTestRespTrials,1);
		[dummy,vecCI_NoHet] = binofit(dblPerformanceNoHetML*intTestTrials,intTestTrials);
		
		%mean removed (heterogeneity+multidim)
		[dblPerformanceNoMeanML,vecDecodedIndexCV,matPosteriorProbabilityCV] = doCrossValidatedDecodingML(matMeanRemovedTestResps',indTestRespTrials,1);
		[dummy,vecCI_NoMean] = binofit(dblPerformanceNoMeanML*intTestTrials,intTestTrials);
		
		%mean+heterogeneity removed (only multidim)
		[dblPerformanceNoHetMeanML,vecDecodedIndexCV,matPosteriorProbabilityCV] = doCrossValidatedDecodingML(matHeteroMeanRemovedTestResps',indTestRespTrials,1);
		[dummy,vecCI_NoHetMean] = binofit(dblPerformanceNoHetMeanML*intTestTrials,intTestTrials);
		
		vecMean = [dblPerformanceRawML dblPerformanceNoMeanML dblPerformanceNoHetML dblPerformanceNoHetMeanML];
		if sParams.boolSavePlots
			figure
			errorbar(1:4,vecMean,([vecCI_Raw(1) vecCI_NoMean(1) vecCI_NoHet(1) vecCI_NoHetMean(1)]-vecMean)/sqrt(intTestTrials),([vecCI_Raw(2) vecCI_NoMean(2) vecCI_NoHet(2) vecCI_NoHetMean(2)]-vecMean)/sqrt(intTestTrials),'xb');
			ylim([0.5 0.7])
			set(gca,'xtick',1:4,'xticklabel',{'Original Data','Mean removed','Het removed','Both removed'});
			ylabel('Hit/miss decoding accuracy');
			
			%save figure
			drawnow;
			strFig = sprintf('%sdec_perf_neural_responses_het_mean_removed_pop%d',strSes,intPopulation);
			export_fig([strFig '.tif']);
			export_fig([strFig '.pdf']);
		end
		
		%save data
		cellSaveHitMissDecoding(:,intPopulation) = vecMean;
		
		%% stimulus presence decoding
		%decode 0% and 100% stimuli for presence
		%afterwards, split for hit/miss trials
		
		%define variables
		vecContrasts = [0.2 0.5 2 8 32 100];
		indNoStim=cellSelectContrasts{1};
		indStim=cellSelectContrasts{end};
		intNumTrials = length(indNoStim);
		indTrials = true(1,intNumTrials);
		vecStimType = cellMultiSes{intPopulation}.structStim.Contrast;
		matDetect = nan(3,6,2); %[upper/mean/lower] x [contrasts] x [hit/miss]
		vecStimDecoded = zeros(1,intNumTrials);
		matPostProb = zeros(2,intNumTrials);
		vecRatioProb = zeros(1,intNumTrials);
		vecOriTrials = (cellMultiSes{intPopulation}.structStim.Orientation/45)+1;
		for intTrial = 1:intNumTrials
			%make likelihood indices for trial
			indLikelihoodStim = indStim;
			indLikelihoodNoStim = indNoStim;
			indLikelihoodStim(intTrial) = false;
			indLikelihoodNoStim(intTrial) = false;
			
			%loop through neurons
			intNeuronCounter = 0;
			vecPrefPop = 1:length(vecNeuronPrefStim);%find(vecNeuronPrefStim == vecOriTrials(intTrial));
			vecStimP = nan(1,length(vecPrefPop));
			vecNoStimP = nan(1,length(vecPrefPop));
			for intNeuron=vecPrefPop %take only pref pop
				intNeuronCounter = intNeuronCounter + 1;
				%build likelihood
				vecLikeStimAct = matTrialResponse(intNeuron,indLikelihoodStim);
				vecLikeNoStimAct = matTrialResponse(intNeuron,indLikelihoodNoStim);
				
				%get data
				dblAct = matTrialResponse(intNeuron,intTrial);
				
				%do decoding
				dblPostStimTemp = normpdf(dblAct,mean(vecLikeStimAct),std(vecLikeStimAct));
				dblPostNoStimTemp = normpdf(dblAct,mean(vecLikeNoStimAct),std(vecLikeNoStimAct));
				
				%put in temp output
				vecStimP(intNeuronCounter) = dblPostStimTemp;
				vecNoStimP(intNeuronCounter) = dblPostNoStimTemp;
			end
			
			%get population decoding
			dblPostNoStim = prod(vecNoStimP);
			dblPostStim = prod(vecStimP);
			
			%put in output
			matPostProb(:,intTrial) = [dblPostNoStim dblPostStim];
			[dummy,vecStimDecoded(intTrial)] = max([dblPostNoStim dblPostStim]);
			vecRatioProb(:,intTrial) = dblPostStim/(dblPostNoStim+dblPostStim);
		end
		
		%decoded stim presence probability + 95% CI
		for intContrast=1:length(vecContrasts)
			
			vecResp = vecStimDecoded(cellSelectContrasts{intContrast} & indSelectRespTrials) == 2;
			vecNoResp = vecStimDecoded(cellSelectContrasts{intContrast} & ~indSelectRespTrials) == 2;
			
			[dblP,dblCI] = binofit(sum(vecResp),length(vecResp));
			matDetect(1,intContrast,1) = dblCI(2); %upper
			matDetect(2,intContrast,1) = dblP; %mean
			matDetect(3,intContrast,1) = dblCI(1); %lower
			
			
			[dblP,dblCI] = binofit(sum(vecNoResp),length(vecNoResp));
			matDetect(1,intContrast,2) = dblCI(2); %upper
			matDetect(2,intContrast,2) = dblP; %mean
			matDetect(3,intContrast,2) = dblCI(1); %lower
			
			%decoding all trials
			[dblP,dblCI] = binofit(sum(vecStimDecoded(cellSelectContrasts{intContrast}) == 2),sum(cellSelectContrasts{intContrast}));
			matDetect(1,intContrast,3) = dblCI(2); %upper
			matDetect(2,intContrast,3) = dblP; %mean
			matDetect(3,intContrast,3) = dblCI(1); %lower
			
			%behavior all trials
			[dblP,dblCI] = binofit(sum(indSelectRespTrials(cellSelectContrasts{intContrast}) == 1),sum(cellSelectContrasts{intContrast}));
			matDetect(1,intContrast,4) = dblCI(2); %upper
			matDetect(2,intContrast,4) = dblP; %mean
			matDetect(3,intContrast,4) = dblCI(1); %lower
			
			
		end
		
		%fit curves
		vecMeanDetect = squeeze(matDetect(2,:,3)); %[contrasts] x [hit/miss]
		
		dblMinY = vecMeanDetect(1);
		dblMaxY = vecMeanDetect(6);
		vecX = 2:5;
		vecDecY = vecMeanDetect(vecX);
		vecBehavY = squeeze(matDetect(2,vecX,4));
		
		lb = [1,eps];
		ub = [6,inf];
		p0 = [4,1];
		[vecParams,resnorm,residual,exitflag] = curvefitfun(@(p0,vecX) getCumGauss(p0,vecX,[dblMinY dblMaxY]),p0,vecX,vecDecY,lb,ub);
		vecFitY = getCumGauss(vecParams,vecX,[dblMinY dblMaxY]);
		
		%calculate R^2 of behavior by decoding
		dblSS_res = sum((vecBehavY - vecDecY).^2);
		dblSS_tot = sum((vecBehavY - mean(vecBehavY)).^2);
		dblR2 = 1 - (dblSS_res/dblSS_tot);
		
		%plot
		if sParams.boolSavePlots
			hDecResp = figure;
			hold on;
			
			vecWindow = 1:length(vecContrasts);
			vecWindowSelect = vecWindow(1):vecWindow(end);
			intWL = length(vecWindowSelect);
			vecWindowInv = intWL:-1:1;
			vecWindowPlotInv = vecWindow(end):-1:vecWindow(1);
			vecX = [vecContrasts(vecWindowSelect) vecContrasts(vecWindowPlotInv)];
			vecLineX = vecContrasts(vecWindowSelect);
			
			for intDetect=[0 1]
				if intDetect == 1
					vecMeanTrace = matDetect(2,:,1);
					vecMinTrace = matDetect(3,:,1);
					vecMaxTrace = matDetect(1,:,1);
					vecColorFill = [0.7 1.0 0.7];
					vecColorLine = [0 1 0];
				else
					vecMeanTrace = matDetect(2,:,2);
					vecMinTrace = matDetect(3,:,2);
					vecMaxTrace = matDetect(1,:,2);
					vecColorLine = [1 0 0];
					vecColorFill = [1 0.7 0.7];
				end
				vecY = [vecMinTrace vecMaxTrace(vecWindowInv)];
				
				%plot
				hold on
				errorbar(vecLineX,vecMeanTrace,abs(vecMinTrace-vecMeanTrace),vecMaxTrace-vecMeanTrace,'Color',vecColorLine);
				hold off
			end
			set(gca,'XScale','log','YScale','linear')
			set(gca,'XTick',vecLineX,'XTickLabel',[0 vecLineX(2:end)])
			title('Stimulus presence decoding with 95% CI')
			grid on
			xlabel('Contrast')
			ylabel('Decoded stimulus presence')
			xlim([min(vecContrasts(vecWindow))-eps max(vecContrasts(vecWindow))])
			ylim([0 1])
			legend({'Miss','Hit'},'Location','Best')
			
			%save figure
			drawnow;
			strFig = sprintf('%s_decodeStimPresence%d_raw',strSes,intPopulation);
			export_fig([strFig '.tif']);
			export_fig([strFig '.pdf']);
		end
		
		%save data
		matMeanDetect = squeeze(matDetect(2,:,:)); %[contrasts] x [hit/miss]
		
		vecDecodedIndexCV = (vecStimDecoded-1);
		
		%perform chi-square
		intHitCorr = sum(vecDecodedIndexCV==indSelectRespTrials&indSelectRespTrials==1);
		intHitIncorr = sum(vecDecodedIndexCV~=indSelectRespTrials&indSelectRespTrials==1);
		intMissCorr = sum(vecDecodedIndexCV==indSelectRespTrials&indSelectRespTrials==0);
		intMissIncorr = sum(vecDecodedIndexCV~=indSelectRespTrials&indSelectRespTrials==0);
		
		%[table,chi2,p,labels] = crosstab(vecDecodedIndexCV,indSelectRespTrials);
		cellSaveMatDetect{intPopulation,1} = matDetect; %[up/mean/low] x [contrasts] x [hit/miss/all/behavior]
		cellSaveMatDetect{intPopulation,2} = vecDecodedIndexCV;
		cellSaveMatDetect{intPopulation,3} = indSelectRespTrials;
		cellSaveMatDetect{intPopulation,4} = vecStimType; %contrasts
	end
	
	%% save data structures
	if sParams.boolSaveData
		vecClock = fix(clock);
		strDate = num2str(vecClock(1:3));
		strFile = ['data_aggregate' strSes '_' strrep(strDate,'    ','_')];
		%save(['D:\Data\Results\stimdetection\' strFile],'-v7.3');
		%if exist(['D:\Data\Results\stimdetection\' strFile '.mat'],'file')
		%	load(['D:\Data\Results\stimdetection\' strFile]);
		%end
		cellVariables=who('cellSave*');
		if isempty(cellVariables)
			warning([mfilename ':NoVariablesToSave'],'No variables found to be saved to disk');
		else
			save(['D:\Data\Results\stimdetection\' strFile],'cellSave*','-v7.3');
		end
		cd(strOldDir);
	end
end

%{
figures to save and take mean over animals
- cellSavePreStim: dF/F detect/no-detect preceding stimulus frames [5 (c) x 2 (d/nd) cell array] x n (animals) with every cell = vector of values
- cellSaveBehavDetect: behavioral detection [5 (c) x n (animals)] cell array with every cell = vector of detect/no detect values (1/0)
- cellSaveBehavRT: behavioral RT [5 (c) x n (animals)] cell array with every cell = vector of RTs
- cellSaveBaseSubtrStimAct: dF/F detect/no-detect baseline-subtracted activity [5 (c) x 2 (d/nd) cell array] x n (animals) with every cell = vector of values
- cellSaveDuringStim: dF/F detect/no-detect during stimulus frames [5 (c) x 2 (d/nd) cell array] x n (animals) with every cell = vector of values
- cellSaveTraceAct: dF/F trace detect/no-detect pref/non-pref [2 (pref/no-pref) x 2 (detect/no-detect) x 6 (c)] x n (animals) with every cell = [Nreps x Nframes] matrix of trace dF/F values
-  AD: cellSaveMatrices: blob matrices in [1 x 4 (detect/no-detect/diff/significance)] cell array with [5 (c) x 31 (bins)] matrix
-  AD: cellSaveHCAENeuronConsistency: [1 x i (neuronal population)] cell array with vector of shuffled correlation values, where first value is actual data
-  AD: cellSaveSignalCorrs: signal correlations [3 (within-HCAE, between, within-non-HCAE) x n (animals/blocks)] with every cell = vector of pairwise correlation values
-  AD: cellSaveNoiseCorrs: noise correlations [3 (within-HCAE, between, within-non-HCAE) x n (animals/blocks)] with every cell = vector of pairwise correlation values
-  AD: cellSaveHCAE: normalized increase in dF/F for detection trials [5 (c) x n (animals/blocks)] cell array with each cell = vector of neurons with normalized dF/F increase for detection trials
-  AD: cellSaveSignalCorrsBiDir: signal correlations [3 (within-HCAE, between, within-non-HCAE) x n (animals/blocks)] with every cell = vector of pairwise correlation values
-  AD: cellSaveNoiseCorrsBiDir: noise correlations [3 (within-HCAE, between, within-non-HCAE) x n (animals/blocks)] with every cell = vector of pairwise correlation values
-  AD: cellSaveNormActDissim: within-group z-scored activation dissimilarity [2 (high/low HCAE) x 2 (detect/no-detect)] x n (animals/blocks) with every cell = vector of dissimilarity values
-  AD: cellSaveDissimCorrITD: inter-trial-distance dependence of assembly consistency [2 (high/low HCAE) x 2 (hits/misses) x n (animals/blocks)] with every cell = [p x 2] matrix with intertrial distance and correlation values for each trial pair
	- extra graph of above with means comparison of hits/misses high/low HCAE

%}