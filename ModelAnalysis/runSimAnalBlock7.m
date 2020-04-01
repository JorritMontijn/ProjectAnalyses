%% merge adjacent stimulus classes for pred subspace & covariance calculation, see if diff corrs show up

%% initialize
%close all;
clearvars;
intType = 2;
for intType =[1 2 3]
	close all
	clearvars -except intType
vecRunAreas = [1];

boolLoadSync = false;
boolLoadLGN = false;
boolLoadSpikeTimes = false;
boolLoadTrialSplit = false;
boolDoSplitAnalysis = false;
%vecRunSims = [-1 -2 11:14];

if intType == 1
	vecRunSims = [-21 -5 -3 -2 -1 11];
	intSizeX = 80;
	intSizeY = 30;
	intUseIndepTrials = 4000;
	intUseRepetitionMax = intUseIndepTrials;
elseif intType == 2
	vecRunSims = [-11 -10];
	intUseIndepTrials = 900;
	intUseRepetitionMax = intUseIndepTrials;
	intSizeX = 40;
	intSizeY = 10;
elseif intType == 3
	vecRunSims = [-7 -6 -4];
	intUseIndepTrials = 2500;
	intUseRepetitionMax = intUseIndepTrials;
	intSizeX = 40;
	intSizeY = 15;
end
for intLoadSim=vecRunSims
	clearvars -except vecRunSims boolDoRegularization intLoadSim bool* intUseRepetitionMax intUseIndepTrials vecRunAreas intSizeX intSizeY
	boolLoad = true;
		
	%% get simulation name [strSimulation] from [intLoadSim]
	loadSim;
	
	%% RUN: #header
	strBlockNr = getFlankedBy(mfilename,'Block','');
	strBlockNr = strBlockNr(1);
	strFigDir = ['D:\Data\Results\Block' strBlockNr '\'];
	strDataDir = ['D:\Data\Results\Data' strBlockNr '\'];
	if isempty(strBlockNr),error;end
	
	if boolLoad
		runAnalysisHeader;
	end
	
	
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
	boolShowThree = true; %show only pred, rand, full CV
	boolIncludeShuff = false;
	if ~exist('boolDoRegularization','var'),boolDoRegularization = true;end
	intFoldK = 10;
	if ~exist('intUseRepetitionMax','var'),intUseRepetitionMax = 4000;end
	if ~exist('intSizeX','var'),intSizeX = 80;end
	if ~exist('intSizeY','var'),intSizeY = 30;end
	boolSaveFigs = true;
	boolSaveData = true;
	dblLambdaInfo = 0;
	dblLambdaPred = 0;
	dblNoiseLevel = 0;
	intMaxDimAnal = 30;
	%intWithinArea = 1; %1, within V1; 2, within V2; 3, across V1-V2;
	intIters = 10;
	
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
		fprintf('Starting area %d; predicting %s [%s]\n',intWithinArea,strPredArea,getTime);
		
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
		dblDiffTheta = range(vecStimTypeOris(matCompareTypes(1,:)));
		
		%check if no spikes
		if (sum(sum(matData(1:intCellsV1,vecTrialStimType==vecUseStimTypes(1)))) + sum(sum(matData(1:intCellsV1,vecTrialStimType==vecUseStimTypes(end))))) == 0
			boolAllZero = true;else boolAllZero = false;end
		indRemUnits = false(intCellsV1,1);
		for intStimType=(vecUseStimTypes(:)')
			indRemUnits(range(matData(:,vecTrialStimType==intStimType),2)==0)=true;
		end
		matData(indRemUnits,:) = [];
		if exist('matDataFull','var')
			matDataFull(indRemUnits,:) = [];
		end
		vecCellArea(indRemUnits) = [];
		
		%% dependence of dimensionality of stimulus intensity
		if ~boolAllZero
			%% get data splits
			%set parameters
			sDataParams=struct;
			sDataParams.intSizeX = intSizeX;
			sDataParams.intSizeY = intSizeY;
			sDataParams.intResamplings = intIters;
			sDataParams.vecCellArea = vecCellArea;
			sDataParams.intWithinArea = intWithinArea;
			sDataParams.vecUseStimTypes = 1:intStimTypes;
			sDataParams.intMaxReps = intUseRepetitionMax;
			%get splits
			[cellMatX,cellNeuronsX,cellMatY,cellNeuronsY,cellTrials] = doDimDataSplits(matData,vecTrialStimType,sDataParams);
			fprintf('  Created %dx%d data splits [%s]\n',size(cellMatX),getTime);
			warning off;
			
			%% predict with reduced rank regression, normal
			if ~boolDoSplitAnalysis
				matDataFull = [];
			end
			
			sPredAndI = doDimPredAndI_RR_Merge(cellNeuronsX,cellNeuronsY,cellTrials,matData,matDataFull,dblDiffTheta,intFoldK,boolIncludeShuff,boolDoRegularization,matCompareTypes);
			
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
				%ylim([dblMinY max(get(gca,'ylim'))]);
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
				title(sprintf('Merged source[N=%d,T=%d] (+/-SEM)',intSizeX,size(cellMatX{1},1)));
				ylabel('Fisher information (d''^2/d{\theta}^2)');
				xlabel('Dimensionality of subspace');
				xlim([min(vecRank)-1 max(vecRank+1)]);
				if max(get(gca,'ylim')) < 0, dblMinY = min(get(gca,'ylim'));else dblMinY=0;end
				%ylim([dblMinY max(get(gca,'ylim'))]);
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
				title(sprintf('Merged target [N=%d,T=%d] (+/-SEM)',intSizeY,size(cellMatX{1},1)));
				ylabel('Fisher information (d''^2/d{\theta}^2)');
				xlabel('Dimensionality of subspace');
				xlim([min(vecRank)-1 max(vecRank+1)]);
				try
					if max(get(gca,'ylim')) < 0, dblMinY = min(get(gca,'ylim'));else dblMinY=0;end
					%ylim([dblMinY max(get(gca,'ylim'))]);
				catch
					
				end
				fixfig;
				
				%save figure
				if boolSaveFigs
					%figure(hFigDD);
					drawnow;
					jFig = get(handle(gcf), 'JavaFrame');
					jFig.setMaximized(true);
					figure(gcf);
					drawnow;
					export_fig([strFigDir  'Block' strBlockNr 'SubspMergeCov_Area' num2str(intWithinArea) '_' strSizeXY strSizeT strTag '_Shuff' strShuff '.tif']);
					export_fig([strFigDir  'Block' strBlockNr 'SubspMergeCov_Area' num2str(intWithinArea) '_' strSizeXY strSizeT strTag '_Shuff' strShuff '.pdf']);
				end
				
				
				%% save data
				if boolSaveData
					sSave = struct;
					%ID, params
					sSave.strName = strType;
					sSave.vecRank = vecRank;
					sSave.sDataParams = sDataParams;
					%R^2
					sSave.matFullR2 = matFullPredCopied;
					sSave.matPredR2 = matPredDimDepR2_CV_Avg;
					sSave.matRandR2 = matPredDimDepR2Rand_CV_Avg;
					%Fisher I
					sSave.matFullI = matFisherDimDepTotalAvg;
					sSave.matPredI = matFisherDimDepProjAvg;
					sSave.matRandI = matFisherDimDepProjRandAvg;
					
					%save
					save([strDataDir 'Block' strBlockNr 'Subsp_Area' num2str(intWithinArea) '_' strSizeXY strSizeT strTag '_Shuff' strShuff],'-struct','sSave');
				end
			end
		else
			
		end
		
		
	end
end
end
