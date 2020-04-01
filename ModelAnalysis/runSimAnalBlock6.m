%% analyze input strength dependency for dimensionality of pop responses

%% initialize
%close all;
clearvars;
intType = 3;

	close all
	clearvars -except intType
vecRunAreas = [3 1];

boolLoadSync = false;
boolLoadLGN = false;
boolLoadSpikeTimes = false;
boolLoadTrialSplit = false;
boolDoSplitAnalysis = false;
%vecRunSims = [-1 -2 11:14];

	
if intType == 1
	
	vecRunSims = [21 23];
	intUseIndepTrials = 4000;
	intUseRepetitionMax = intUseIndepTrials;
	vecRunSizesX = [40:40:200];
	vecRunSizesY = [40:40:200];
elseif intType == 2
	vecRunSims = [21 23];
	intUseIndepTrials = 4000;
	intUseRepetitionMax = intUseIndepTrials;
	vecRunSizesX = [40:40:200];
	vecRunSizesY = 30*ones(size(vecRunSizesX));
elseif intType == 3
	vecRunSims = [21 23];
	intUseIndepTrials = 4000;
	intUseRepetitionMax = intUseIndepTrials;
	vecRunSizesY = [40:40:200];
	vecRunSizesX = 200*ones(size(vecRunSizesY));
end
for intLoadSim=vecRunSims
	close all;
	clearvars -except vecRunSims intLoadSim bool* intUseRepetitionMax intUseIndepTrials vecRunAreas intSizeIdx vecRunSizesY vecRunSizesX
	boolLoad = true;
		
	%% get simulation name [strSimulation] from [intLoadSim]
	loadSim;
	
	%% RUN: #header
	strFigDir = 'D:\Data\Results\Block6\';
	if boolLoad
		runAnalysisHeader;
	end
	%return
	%save(['DataForRex_' strSimulation '.mat'],'matModelResp','vecCellArea','vecTrialStimType');
	%return
	
	cellIn = strsplit(strSimulation,'_');
	strFiller = cellIn{1};
	strType = cell2mat(cellIn(2:(end-1)));
	strDate = cellIn{end};
	strTag = [strType '_' strDate];
	if boolDoSplitAnalysis
		strTag = ['TS_' strTag];
		strType = ['TS_' strType];
	end
	
	
	%% set parameters
	for intSizeIdx=1:numel(vecRunSizesX)
	
	boolShowThree = true; %show only pred, rand, full CV
	boolIncludeShuff = false;
	boolDoRegularization = true;
	intFoldK = 10;
	if ~exist('intUseRepetitionMax','var'),intUseRepetitionMax = 4000;end
	intSizeX = vecRunSizesX(intSizeIdx);
	intSizeY = vecRunSizesY(intSizeIdx);
	boolSaveFigs = true;
	boolSaveData = false;
	dblLambdaInfo = 0;
	dblLambdaPred = 0;
	dblNoiseLevel = 0;
	intMaxDimAnal = 30;
	%intWithinArea = 1; %1, within V1; 2, within V2; 3, across V1-V2;
	intIters = 5;
	
	if ~boolDoRegularization
		strTag = ['NoR_' strTag];
		strType = ['NoR_' strType];
	end
	
	%% calculate which parameter to use
	strSizeXY = sprintf('X%dY%d',intSizeX,intSizeY);
	strSizeT = strcat('T',num2str(intUseRepetitionMax));
	vecRanges = range(matStimTypeCombos(2:end,:),2);
	[d,intUseParam]=max(vecRanges); %exclude orientation
	intUseParam = intUseParam + 1;
	if all(vecRanges==0), intUseParam = 0;vecParamStimType = [0];strParam='None';
		vecParamIdx = ones(1,size(matStimTypeCombos,2));
		vecParamVals = 1;
		intParamNum=1;
		mapC = redbluepurple(intParamNum);
	else
		if intUseParam == 1,vecParamStimType = vecStimTypeOris;strParam='Ori';
		elseif intUseParam == 2,vecParamStimType = vecStimTypeSFs;strParam='SF';
		elseif intUseParam == 3,vecParamStimType = vecStimTypeTFs;strParam='TF';
		elseif intUseParam == 4,vecParamStimType = vecStimTypeContrasts;strParam='Contrast';
		elseif intUseParam == 5,vecParamStimType = vecStimTypeLuminance;strParam='Luminance';
		end
		vecParamIdx = matStimTypeCombos(intUseParam,:);
		vecParamVals = unique(vecParamStimType);
		intParamNum=numel(vecParamVals);
		mapC = redbluepurple(intParamNum);
	end
	cellStrArea = {'V1 from V1','V2 from V2','V2 from V1'};
	
	%% build comparison matrix for stim types
	vecUniqueStimTypes = unique(vecTrialStimType);
	intMultiComp = false;
	intStimTypes = numel(vecUniqueStimTypes);
	if ~exist('matCompareTypes','var')
		if sum(range(matStimTypeCombos,2)>0)>1
			intMultiComp = true;
			matCompareTypes = [ones(1,size(matStimTypeCombos,2)-1)' (2:size(matStimTypeCombos,2))'];
			matCompareTypes = [1 2; 1 3; 1 4];
		else
			
			matCompareTypes = [(1:intStimTypes)' circshift((1:intStimTypes)',-1)];
		end
	end
	
	%% get data
	if boolLoadTrialSplit
		fprintf('Transforming data [%s]\n',getTime);
		[vecTrialStimTypeUnsplit,matRespUnsplit,vecTrialStimTypeSplit,matRespSplitNorm,matRespSplit,matSplitTrialIdx] = ...
			prepSplitData(matModelRespTS3,vecTrialStimType);
		
		if boolDoSplitAnalysis
			matData = matRespSplitNorm;
			matDataFull = matRespSplit;
			vecTrialStimType = vecTrialStimTypeSplit;
		else
			matData = matRespUnsplit;
			vecTrialStimType = vecTrialStimTypeUnsplit;
		end
	else
		matData = double(matModelResp);
		vecTrialStimType = vecTrialStimType;
	end
	%[~,~,dblInfoCheck]=getSeparation(matData,vecTrialStimType,false,range(vecStimTypeOris([1 2])))
	
	%% run analyses on different areas
	for intWithinArea=vecRunAreas
		%msg
		strPredArea = cellStrArea{intWithinArea};
		fprintf('Starting area %d; predicting %s [X=%d,Y=%d] [%s]\n',intWithinArea,strPredArea,intSizeX,intSizeY,getTime);
		
		%% pre-allocate variables dim dep
		hFigDD = figure;
		%full screen
		drawnow;
		jFig = get(handle(gcf), 'JavaFrame');
		jFig.setMaximized(true);
		figure(gcf);
		drawnow;
		
		%% start
		vecUseStimTypes = unique(round(vecTrialStimType));
		indUseTrials = ismember(vecTrialStimType,vecUseStimTypes);
		
		%% set default parameters
		dblDiffTheta = range(vecStimTypeOris(vecUseStimTypes([1 2])));
		
		%check if no spikes
		if (sum(sum(matData(1:intCellsV1,vecTrialStimType==vecUseStimTypes(1)))) + sum(sum(matData(1:intCellsV1,vecTrialStimType==vecUseStimTypes(end))))) == 0
			boolAllZero = true;else boolAllZero = false;end
		
		%% dependence of dimensionality of stimulus intensity
		if ~boolAllZero
			%% get data splits
			%set parameters
			sParamsAnalSplit=struct;
			sParamsAnalSplit.intSizeX = intSizeX;
			sParamsAnalSplit.intSizeY = intSizeY;
			sParamsAnalSplit.intResamplings = intIters;
			sParamsAnalSplit.vecCellArea = vecCellArea;
			sParamsAnalSplit.intWithinArea = intWithinArea;
			sParamsAnalSplit.vecUseStimTypes = vecUseStimTypes;
			sParamsAnalSplit.intMaxReps = intUseRepetitionMax;
			%get splits
			[cellMatX,cellNeuronsX,cellMatY,cellNeuronsY,cellTrials] = doDimDataSplits(matData,vecTrialStimType,sParamsAnalSplit);
			fprintf('  Created %dx%d data splits [%s]\n',size(cellMatX),getTime);
			
			%% predict with reduced rank regression, normal
			if boolDoSplitAnalysis
				sFullData = struct;
				sFullData.matData = matDataFull;
				sFullData.cellNeuronsX = cellNeuronsX;
				sFullData.cellNeuronsY = cellNeuronsY;
				sFullData.cellTrials = cellTrials;
				
				sPredAndI = doDimPredAndI_RR(cellMatX,cellMatY,dblLambdaInfo,dblDiffTheta,sFullData,intFoldK,boolIncludeShuff,boolDoRegularization,matCompareTypes);
				
			else
				sPredAndI = doDimPredAndI_RR(cellMatX,cellMatY,dblLambdaInfo,dblDiffTheta,[],intFoldK,boolIncludeShuff,boolDoRegularization,matCompareTypes);
			end
			
			if intMultiComp
				return
				%%
				if ~exist('sPredAndI2','var')
					sPredAndI2 = sPredAndI;
				end
				sPredAndI = sPredAndI2;
				intAnalyzeComp = 3;
				cellFields = fieldnames(sPredAndI);
				for intField=1:numel(cellFields)
					strField = cellFields{intField};
					if strcmpi(strField(1:3),'mat')
						intDims = ndims(sPredAndI.(strField));
						strExpR = ['sPredAndI.' strField '=sPredAndI.' strField '(:,' num2str(intAnalyzeComp) repmat(',:',[1 intDims-2]) ');'];
						eval(strExpR);
					end
				end
			end
			
			%% get dimensionality with factor analysis
			%[matSourceExplainedVarLL,matSourcePopDimensionality,matSourcePopLogLikeTest] = doDimFactorAnal(cellMatX,intMaxDimAnal);
			%fprintf('  C%s; Mean number of dominant dimensions in source: <%s\b> [%s]\n',strC,sprintf('%.2f ',xmean(matSourcePopDimensionality,1)),getTime);
			%[matTargetExplainedVarLL,matTargetPopDimensionality,matTargetPopLogLikeTest] = doDimFactorAnal(cellMatY,intMaxDimAnal);
			%fprintf('  C%s; Mean number of dominant dimensions in target: <%s\b> [%s]\n',strC,sprintf('%.2f ',xmean(matTargetPopDimensionality,1)),getTime);
			
			%% extract
			cellShuffStr = {'none','neurons','reps','xTrials'};
			vecRank = sPredAndI.vecRank;
			intMaxShuff = size(sPredAndI.matPredFullR2_CV,4);
			for intShuffIdx=1:intMaxShuff
				strShuff = cellShuffStr{intShuffIdx};
				
				matPredFullR2_CV = nanmean(sPredAndI.matPredFullR2_CV(:,:,:,intShuffIdx),3);
				matPredDimDepR2_CV = nanmean(sPredAndI.matPredDimDepR2_CV(:,:,:,:,intShuffIdx),4);
				matPredDimDepR2_NonCV = nanmean(sPredAndI.matPredDimDepR2_NonCV(:,:,:,:,intShuffIdx),4);
				matPredDimDepR2_Rem =  nanmean(sPredAndI.matPredDimDepR2_Rem(:,:,:,:,intShuffIdx),4);
				
				matPredDimDepR2Rand_CV = nanmean(sPredAndI.matPredDimDepR2Rand_CV(:,:,:,:,:,intShuffIdx),5);
				matPredDimDepR2Rand_NonCV = nanmean(sPredAndI.matPredDimDepR2Rand_NonCV(:,:,:,:,:,intShuffIdx),5);
				matPredDimDepR2RandRem = nanmean(sPredAndI.matPredDimDepR2RandRem(:,:,:,:,:,intShuffIdx),5);
				matFisherDimDepProj = nanmean(sPredAndI.matFisherDimDepProj(:,:,:,:,intShuffIdx),4);
				matFisherDimDepOrth = nanmean(sPredAndI.matFisherDimDepOrth(:,:,:,:,intShuffIdx),4);
				matFisherDimDepTotal = nanmean(sPredAndI.matFisherDimDepTotal(:,:,:,:,intShuffIdx),4);
				
				matFisherDimDepTotalRand = nanmean(sPredAndI.matFisherDimDepTotalRand(:,:,:,:,:,intShuffIdx),5);
				matFisherDimDepProjRand = nanmean(sPredAndI.matFisherDimDepProjRand(:,:,:,:,:,intShuffIdx),5);
				matFisherDimDepOrthRand = nanmean(sPredAndI.matFisherDimDepOrthRand(:,:,:,:,:,intShuffIdx),5);
				
				matFisherDimDepConsecProj = sPredAndI.matFisherDimDepConsecProj(:,:,:,intShuffIdx);
				matFisherDimDepConsecOrth = sPredAndI.matFisherDimDepConsecOrth(:,:,:,intShuffIdx);
				matFisherDimRedundancyPerc = sPredAndI.matFisherDimRedundancyPerc(:,:,:,intShuffIdx);
				matNoisePredR2ConsecProj = sPredAndI.matNoisePredR2ConsecProj(:,:,:,intShuffIdx);
				
				matFisherDimDepProjT = nanmean(sPredAndI.matFisherDimDepProjT(:,:,:,:,intShuffIdx),4);
				matFisherDimDepOrthT = nanmean(sPredAndI.matFisherDimDepOrthT(:,:,:,:,intShuffIdx),4);
				matFisherDimDepTotalT = nanmean(sPredAndI.matFisherDimDepTotalT(:,:,:,:,intShuffIdx),4);
				
				matFisherDimDepTotalRandT = nanmean(sPredAndI.matFisherDimDepTotalRandT(:,:,:,:,:,intShuffIdx),5);
				matFisherDimDepProjRandT = nanmean(sPredAndI.matFisherDimDepProjRandT(:,:,:,:,:,intShuffIdx),5);
				matFisherDimDepOrthRandT = nanmean(sPredAndI.matFisherDimDepOrthRandT(:,:,:,:,:,intShuffIdx),5);
				
				
				%transform
				matPredDimDepR2_CV_Avg = squeeze(nanmean(matPredDimDepR2_CV,2));
				matPredDimDepR2_NonCV_Avg = squeeze(nanmean(matPredDimDepR2_NonCV,2));
				
				matPredDimDepR2RemAvg = squeeze(nanmean(matPredDimDepR2_Rem,2));
				matPredDimDepR2Rand_CV_Avg = squeeze(nanmean(nanmean(matPredDimDepR2Rand_CV,4),2));
				matPredDimDepR2Rand_NonCV_Avg = squeeze(nanmean(nanmean(matPredDimDepR2Rand_NonCV,4),2));
				matPredDimDepR2RandRemAvg = squeeze(nanmean(nanmean(matPredDimDepR2RandRem,4),2));
				matFisherDimDepProjAvg = squeeze(nanmean(matFisherDimDepProj,2));
				matFisherDimDepOrthAvg = squeeze(nanmean(matFisherDimDepOrth,2));
				matFisherDimDepTotalAvg = squeeze(nanmean(matFisherDimDepTotal,2));
				matFisherDimDepProjRandAvg = squeeze(nanmean(nanmean(matFisherDimDepProjRand,4),2));
				matFisherDimDepOrthRandAvg = squeeze(nanmean(nanmean(matFisherDimDepOrthRand,4),2));
				matFisherDimDepTotalRandAvg = squeeze(nanmean(nanmean(matFisherDimDepTotalRand,4),2));
				
				matFisherDimDepConsecProjAvg = squeeze(nanmean(matFisherDimDepConsecProj,2));
				matFisherDimDepConsecOrthAvg = squeeze(nanmean(matFisherDimDepConsecOrth,2));
				%matFisherDimDepConsecProjAvg2 = squeeze(nanmean(matFisherDimDepConsecProj2,2));
				%matFisherDimDepConsecOrthAvg2 = squeeze(nanmean(matFisherDimDepConsecOrth2,2));
				matFisherDimDepRedundancyAvg = squeeze(nanmean(matFisherDimRedundancyPerc,2));
				matPredConsecProjR2Avg = squeeze(nanmean(matNoisePredR2ConsecProj,2));
				
				matFisherDimDepProjAvgT = squeeze(nanmean(matFisherDimDepProjT,2));
				matFisherDimDepOrthAvgT = squeeze(nanmean(matFisherDimDepOrthT,2));
				matFisherDimDepTotalAvgT = squeeze(nanmean(matFisherDimDepTotalT,2));
				matFisherDimDepProjRandAvgT = squeeze(nanmean(nanmean(matFisherDimDepProjRandT,4),2));
				matFisherDimDepOrthRandAvgT = squeeze(nanmean(nanmean(matFisherDimDepOrthRandT,4),2));
				matFisherDimDepTotalRandAvgT = squeeze(nanmean(nanmean(matFisherDimDepTotalRandT,4),2));
				
				%% predict with full pop
				fprintf('  Mean full population predictions: <%s\b> [%s]\n',sprintf('%.4f ',xmean(matPredFullR2_CV,1)),getTime);
				%msg
				
				%% do plotting
				if ishandle(hFigDD),figure(hFigDD);end
				clf;
				subplot(2,3,1)
				matFullPredCopied = repmat(mean(matPredFullR2_CV,2),[1 numel(vecRank)]);
				hold on
				errorbar(vecRank,mean(matFullPredCopied),std(matFullPredCopied)/sqrt(intIters),'k-');
				errorbar(vecRank,mean(matPredDimDepR2_CV_Avg),std(matPredDimDepR2_CV_Avg)/sqrt(intIters),'b-');
				errorbar(vecRank,mean(matPredDimDepR2Rand_CV_Avg),std(matPredDimDepR2Rand_CV_Avg)/sqrt(intIters),'r-');
				if ~boolShowThree
					errorbar(vecRank,mean(matPredDimDepR2_NonCV_Avg),std(matPredDimDepR2_NonCV_Avg)/sqrt(intIters),'b--');
					errorbar(vecRank,mean(matPredDimDepR2Rand_NonCV_Avg),std(matPredDimDepR2Rand_NonCV_Avg)/sqrt(intIters),'r--');
					%errorbar(vecRank,mean(matPredDimDepR2RandRemAvg),std(matPredDimDepR2RandRemAvg)/sqrt(intIters),'r--');
					errorbar(vecRank,mean(matPredDimDepR2RemAvg),std(matPredDimDepR2RemAvg)/sqrt(intIters),'c-');
				end
				hold off
				title(sprintf('Pred %s; shuff %s; %s,%s',strPredArea,strShuff,strType,strSizeT),'Interpreter','none');
				ylabel('Noise predictability (R^2)');
				xlabel('Dimensionality of subspace');
				xlim([min(vecRank)-1 max(vecRank+1)]);
				if max(get(gca,'ylim')) < 0, dblMinY = min(get(gca,'ylim'));else dblMinY=0;end
				ylim([dblMinY max(get(gca,'ylim'))]);
				fixfig;
				h=legend({'Full model CV','Predictive subsp CV','Random subsp CV','Pred subsp train','Random Train','Orth to pred subsp'});
				set(h,'Location','Best');
				
				subplot(2,3,2)
				%matFullMinusPred = matFisherDimDepTotalAvg - matFisherDimDepProjAvg;
				hold on
				errorbar(vecRank,mean(matFisherDimDepTotalAvg),std(matFisherDimDepTotalAvg)/sqrt(intIters),'k-');
				errorbar(vecRank,mean(matFisherDimDepProjAvg),std(matFisherDimDepProjAvg)/sqrt(intIters),'b-');
				errorbar(vecRank,mean(matFisherDimDepProjRandAvg),std(matFisherDimDepProjRandAvg)/sqrt(intIters),'r-');
				if ~boolShowThree
					errorbar(vecRank,mean(matFisherDimDepOrthRandAvg),std(matFisherDimDepOrthRandAvg)/sqrt(intIters),'m-');
					errorbar(vecRank,mean(matFisherDimDepOrthAvg),std(matFisherDimDepOrthAvg)/sqrt(intIters),'c-');
					%errorbar(vecRank,mean(matFullMinusPred),std(matFullMinusPred)/sqrt(intIters),'m-');
				end
				hold off
				%h2=legend({'Total','Orth to Rand Subsp','Orth to Pred Subsp','Pred Subspace','Rand Subspace','(Full) - (Pred Subsp)'});
				%h2=legend({'Total','Orth to Rand Subsp','Orth to Pred Subsp','Pred Subspace','Rand Subspace'});
				%set(h2,'Location','west');
				title(sprintf('Info in source subsp [N=%d,T=%d] (+/-SEM)',intSizeX,size(cellMatX{1},1)));
				ylabel('Fisher information (d''^2/d{\theta}^2)');
				xlabel('Dimensionality of subspace');
				xlim([min(vecRank)-1 max(vecRank+1)]);
				if max(get(gca,'ylim')) < 0, dblMinY = min(get(gca,'ylim'));else dblMinY=0;end
				ylim([dblMinY max(get(gca,'ylim'))]);
				fixfig;
				
				
				subplot(2,3,3)
				hold on
				errorbar(vecRank,mean(matFisherDimDepTotalAvgT),std(matFisherDimDepTotalAvgT)/sqrt(intIters),'k-');
				%errorbar(vecRank,mean(matFisherDimDepOrthRandAvgT),std(matFisherDimDepOrthRandAvgT)/sqrt(intIters),'r--');
				%errorbar(vecRank,mean(matFisherDimDepOrthAvgT),std(matFisherDimDepOrthAvgT)/sqrt(intIters),'b--');
				errorbar(vecRank,mean(matFisherDimDepProjAvgT),std(matFisherDimDepProjAvgT)/sqrt(intIters),'b-');
				errorbar(vecRank,mean(matFisherDimDepProjRandAvgT),std(matFisherDimDepProjRandAvgT)/sqrt(intIters),'r-');
				hold off
				%h=legend({'Total','Orth to Rand Subsp','Orth to Pred Subsp','Pred Subspace','Rand Subspace'});
				%set(h,'Location','Best');
				title(sprintf('Info in target subsp [N=%d,T=%d] (+/-SEM)',intSizeY,size(cellMatX{1},1)));
				ylabel('Fisher information (d''^2/d{\theta}^2)');
				xlabel('Dimensionality of subspace');
				xlim([min(vecRank)-1 max(vecRank+1)]);
				try
					if max(get(gca,'ylim')) < 0, dblMinY = min(get(gca,'ylim'));else dblMinY=0;end
					ylim([dblMinY max(get(gca,'ylim'))]);
				catch
					
				end
				fixfig;
				%{
				subplot(2,3,4)
				matPercR2 = (matPredDimDepR2Avg./matFullPredCopied)*100;
				matPercI_PredSubsp = (matFisherDimDepProjAvg./matFisherDimDepTotalAvg)*100;
				matPercI_RandSubsp = (matFisherDimDepProjRandAvg./matFisherDimDepTotalAvg)*100;
				hold on
				errorbar(vecRank,mean(matPercR2),std(matPercR2)/sqrt(intIters),'g-');
				errorbar(vecRank,mean(matPercI_PredSubsp),std(matPercI_PredSubsp)/sqrt(intIters),'b-');
				errorbar(vecRank,mean(matPercI_RandSubsp),std(matPercI_RandSubsp)/sqrt(intIters),'r-');
				hold off
				h=legend({'R^2 prediction','Fisher I in Pred Subsp','Fisher I in Rand Subsp'});
				set(h,'Location','Best');
				title('Prediction and Fisher info in % of full model');
				ylabel('% of full model');
				xlabel('Dimensionality of subspace');
				xlim([min(vecRank)-1 max(vecRank+1)]);
				ylim([0 110]);
				fixfig;
				%}
				%{
				subplot(2,3,6)
				matPredPlusOrth = matFisherDimDepOrthAvg + matFisherDimDepProjAvg;
				hold on
				errorbar(vecRank,mean(matFisherDimDepTotalAvg),std(matFisherDimDepTotalAvg)/sqrt(intIters),'k-');
				errorbar(vecRank,mean(matPredPlusOrth),std(matPredPlusOrth)/sqrt(intIters),'b-');
				hold off
				title('Redundancy Fisher info');
				ylabel('Fisher information (d''^2)');
				xlabel('Dimensionality of subspace');
				xlim([min(vecRank)-1 max(vecRank+1)]);
				ylim([0 max(get(gca,'ylim'))]);
				fixfig;
				h=legend({'I(Total)','I(Pred subsp) + I(Orth)'});
				set(h,'Location','Best');
				
				subplot(2,3,4)
				hold on
				matPredConsecProjR2AvgCumSum = matPredConsecProjR2Avg;
				errorbar(vecRank,mean(matFullPredCopied),std(matFullPredCopied)/sqrt(intIters),'k-');
				errorbar(vecRank,mean(matPredConsecProjR2AvgCumSum),std(matPredConsecProjR2AvgCumSum)/sqrt(intIters),'b-');
				errorbar(vecRank,mean(matPredDimDepR2RandAvg),std(matPredDimDepR2RandAvg)/sqrt(intIters),'r-');
				hold off
				title('Prediction using consecutive coding dimensions');
				ylabel('Cumulative noise predictability (R^2)');
				xlabel('Optimal dimension #');
				xlim([0 10]);
				ylim([0 max(get(gca,'ylim'))]);
				fixfig;
				h=legend({'Full model','\Sigma^{-1} * f''','Random subspace'});
				set(h,'Location','Best');
				
				subplot(2,3,5)
				hold on
				errorbar(vecRank,mean(matFisherDimDepTotalAvg),std(matFisherDimDepTotalAvg)/sqrt(intIters),'k-');
				errorbar(vecRank,mean(matFisherDimDepConsecProjAvg),std(matFisherDimDepConsecProjAvg)/sqrt(intIters),'b-');
				%errorbar(vecRank,mean(matFisherDimDepConsecProjAvg2),std(matFisherDimDepConsecProjAvg2)/sqrt(intIters),'b--');
				hold off
				title('Fisher I with consecutive coding dimensions');
				ylabel('Fisher information (d''^2)');
				xlabel('Optimal dimension #');
				xlim([0 10]);
				ylim([0 max(get(gca,'ylim'))]);
				fixfig;
				h=legend({'Total','Consec Proj'});
				set(h,'Location','Best');
				%}
				%save figure
				if boolSaveFigs
					%figure(hFigDD);
					drawnow;
					jFig = get(handle(gcf), 'JavaFrame');
					jFig.setMaximized(true);
					figure(gcf);
					drawnow;
					export_fig([strFigDir  'Block6Subsp_Area' num2str(intWithinArea) '_' strSizeXY strSizeT strTag '_Shuff' strShuff '.tif']);
					export_fig([strFigDir  'Block6Subsp_Area' num2str(intWithinArea) '_' strSizeXY strSizeT strTag '_Shuff' strShuff '.pdf']);
				end
			end
		else
			
		end
		
		
	end
end
end
%end
