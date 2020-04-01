%% analyze input strength dependency for dimensionality of pop responses
%% initialize
close all;
clearvars;
intType = 3;
vecRunAreas = 1;%[1 2];

boolLoadSync = false;
boolLoadLGN = false;
boolLoadSpikeTimes = false;
%vecRunSims = [-1 -2 11:14];
boolLoadTrialSplit = false;
boolDoSplitAnalysis = false;

if intType == 1
	boolLoadTrialSplit = false;
	%vecRunSims = [-21 -11 -10 -7:-1 21:26 31:33];
	%vecRunSims = [-21 -11 -10 -7:-1];
	vecRunSims = 2;%[-21];
	intUseRepetitionMax = 800;
elseif intType == 2
	boolLoadTrialSplit = false;
	vecRunSims = [112];
	intUseRepetitionMax = 400;
elseif intType == 3
	boolLoadTrialSplit = false;
	%vecRunSims = [-1 -2];
	vecRunSims = [100 110 150 1 2 102 104 106 108];
	intUseIndepTrials = 4000;
	intUseRepetitionMax = 1000;%intUseIndepTrials*10;
elseif intType == 4
	boolLoadTrialSplit = true;
	%vecRunSims = [-1 -2];
	vecRunSims = [113 114];
	intUseIndepTrials = 400;
	intUseRepetitionMax = intUseIndepTrials*10;
end
for intLoadSim=vecRunSims
	clearvars -except vecRunSims intLoadSim bool* intUseRepetitionMax intUseIndepTrials vecRunAreas
	boolLoad = true;
	
	%% get simulation name [strSimulation] from [intLoadSim]
	loadSim;
	
	%% RUN: #header
	strFigDir = 'D:\Data\Results\Block1\';
	if boolLoad
		runAnalysisHeader;
	end
	
	%% build comparison matrix for stim types
	vecUniqueStimTypes = unique(vecTrialStimType);
	intStimTypes = numel(vecUniqueStimTypes);
	if ~exist('matCompareTypes','var')
		if sum(range(matStimTypeCombos,2)>0)>1
			intMultiComp = true;
			matCompareTypes = [ones(1,size(matStimTypeCombos,2)-1)' (2:size(matStimTypeCombos,2))'];
			matCompareTypes = [1 2; 1 3; 1 4];
		else
			intMultiComp = false;
			matCompareTypes = [(1:intStimTypes)' circshift((1:intStimTypes)',-1)];
		end
	end
	cellIn = strsplit(strSimulation,'_');
	strFiller = cellIn{1};
	strType = cellIn{2};
	strDate = cellIn{3};
	strTag = [strType '_' strDate];
	
	%% set parameters
	
	dblDiffTheta = range(vecStimTypeOris([1 2]));
	strTag = [sprintf('dTheta%d',round(dblDiffTheta)) strTag];
	boolBiasCorrection = true;
	boolSaveFigs = true;
	boolSaveData = true;
	dblLambdaInfo = 1;
	dblLambdaPred = 0;
	boolAnalI = true;
	boolAnalPred = false;
	boolAnalPopM = false;
	intMaxDimAnal = 30;
	%intWithinArea = 1; %1, within V1; 2, within V2; 3, across V1-V2;
	intIters = 50;
	if ~boolBiasCorrection
		strTag = ['NBC_' strTag];
	end
	if boolLoadTrialSplit
		strTag = ['TS_' strTag];
	end
	
	%% calculate which parameter to use
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
	
	%% transform model resp
	if boolLoadTrialSplit
		vecTrialStimType = vecTrialStimType(vecTrialIdxTS);
		matModelRespP = double(matModelRespTS);
	else
		matModelRespP = double(matModelResp);
	end
	intNeurons = size(matModelRespP,1);
	intTrials = size(matModelRespP,2);
	
	if exist('matModelRespLGN_ON','var')
		matModelRespLGN_ON = double(matModelRespLGN_ON);
		matModelRespLGN_OFF = double(matModelRespLGN_OFF);
		intNeuronsLGN = size(matModelRespLGN_ON,1);
	end
	for intWithinArea=vecRunAreas
		%msg
		fprintf('Starting area %d [%s]\n',intWithinArea,getTime);
		
		%% pre-allocate
		cellPredictionsR2_SingleNeurons=cell(1,intParamNum);
		cellPredictionsR2=cell(1,intParamNum);
		cellPredDimDepR2=cell(1,intParamNum);
		cellPredictiveDimensions=cell(1,intParamNum);
		cellTargetPopDimensionality=cell(1,intParamNum);
		cellR2RemPredDim=cell(1,intParamNum);
		cellR2UseDomDim=cell(1,intParamNum);
		cellI=cell(1,intParamNum);
		cellI_shuff=cell(1,intParamNum);
		cellI_NCV_=cell(1,intParamNum);
		cellI_NCV_shuff=cell(1,intParamNum);
		
		matNoiseCorrsS1 = [];
		matNoiseCorrsS2 = [];
		matNC_Diffs = [];
		
		%% pre-allocate variables Fisher I
		if boolAnalI && (intWithinArea == 1 || intWithinArea == 2)
			intUseN = min(sum(vecCellArea==intWithinArea),round(intTrials/intStimTypes));
			intUseN = 110;
			if intUseN < 30
				continue;
			elseif intUseN < 200
				vecGroupSizes = [10:10:(roundi(intUseN-10,-1,'floor'))];
			elseif intUseN < 600
				vecGroupSizes = [10:10:40 50:50:(roundi(intUseN-10,-2,'floor'))];
			elseif intUseN < 2000
				vecGroupSizes = [10:10:40 50:50:500 600:100:(roundi(intUseN-10,-2,'floor'))];
			else
				vecGroupSizes = [10:10:40 50:50:500 600:100:1000 1500:500:(roundi(intUseN-10,-3,'floor'))];
			end
			vecGroupSizes(vecGroupSizes>intUseN) = [];
			matI_LR = nan(intParamNum,intIters,numel(vecGroupSizes));
			matI_shuff = nan(intParamNum,intIters,numel(vecGroupSizes));
			matI_Direct = nan(intParamNum,intIters,numel(vecGroupSizes));
			matI_Direct_shuff = nan(intParamNum,intIters,numel(vecGroupSizes));
			matI_NCV = nan(intParamNum,intIters,numel(vecGroupSizes));
			matI_NCV_shuff = nan(intParamNum,intIters,numel(vecGroupSizes));
			
			
			
			for intCompareIdx = 1:size(matCompareTypes,1)
				vecUseStimTypes = matCompareTypes(intCompareIdx,:);
				
				%check if no spikes
				if (sum(sum(matModelRespP(1:intCellsV1,vecTrialStimType==vecUseStimTypes(1)))) + sum(sum(matModelRespP(1:intCellsV1,vecTrialStimType==vecUseStimTypes(end))))) == 0
					boolAllZero = true;else boolAllZero = false;end
				
				%% calculate information in cortex
				%set parameters
				sParamsAnal = struct;
				sParamsAnal.vecUseStimTypes = vecUseStimTypes;
				sParamsAnal.dblLambda = dblLambdaInfo;
				sParamsAnal.boolDirectI = true;
				sParamsAnal.boolLogRegI = false;
				sParamsAnal.intIters = intIters;
				sParamsAnal.vecGroupSizes = vecGroupSizes;
				sParamsAnal.boolVerbose = false;
				sParamsAnal.dblDiffTheta = dblDiffTheta;
				sParamsAnal.boolBiasCorrection = boolBiasCorrection;
				if intWithinArea == 1
					vecCells = 1:intCellsV1;
				elseif intWithinArea == 2
					vecCells = (intCellsV1+1):size(matModelRespP,1);
				end
				if boolAllZero
					%get output
					matI_LR(intCompareIdx,:) = 0;
					matI_shuff(intCompareIdx,:) = 0;
				else
					if ~isinf(intUseRepetitionMax)
						vecTrials1 = find(vecTrialStimType==vecUseStimTypes(1));
						vecTrials1 = vecTrials1(randperm(numel(vecTrials1),intUseRepetitionMax));
						vecTrials2 = find(vecTrialStimType==vecUseStimTypes(2));
						vecTrials2 = vecTrials2(randperm(numel(vecTrials2),intUseRepetitionMax));
						vecUseTrials = sort([vecTrials1 vecTrials2],'ascend');
					else
						vecUseTrials = 1:size(matModelRespP,2);
					end
					matData = matModelRespP(vecCells,vecUseTrials);
					vecUseTrialStimType = vecTrialStimType(vecUseTrials);
					
					%do analysis
					[matFisherFull,sOut] = doFisherAnal(matData,vecUseTrialStimType,sParamsAnal);
					if sParamsAnal.boolLogRegI
						%get output
						cellI{intCompareIdx} = nanmean(sOut.matI_LogReg_bc_CV,3); %mean across folds
						cellI_shuff{intCompareIdx} = nanmean(sOut.matI_LogReg_bc_CV_shuff,3); %mean across folds
						matI_LR(intCompareIdx,:,:) = nanmean(sOut.matI_LogReg_bc_CV,3); %mean across folds and resmplings
						matI_shuff(intCompareIdx,:,:) = nanmean(sOut.matI_LogReg_bc_CV_shuff,3); %mean across folds and resmplings
					end
					if sParamsAnal.boolDirectI
						matI_NCV(intCompareIdx,:,:) = nanmean(sOut.matI_LogReg_bc,3); %mean across folds and resmplings
						matI_NCV_shuff(intCompareIdx,:,:) = nanmean(sOut.matI_LogReg_bc_shuff,3); %mean across folds and resmplings
						
						matI_Direct(intCompareIdx,:,:) = nanmean(sOut.matI_Direct_bc,3); %mean across folds and resmplings
						matI_Direct_shuff(intCompareIdx,:,:) = nanmean(sOut.matI_Direct_bc_shuff,3); %mean across folds and resmplings
					end
					cellI_NCV{intCompareIdx} = nanmean(sOut.matI_LogReg_bc,3); %mean across folds
					cellI_NCV_shuff{intCompareIdx} = nanmean(sOut.matI_LogReg_bc_shuff,3); %mean across folds
				end
				
				
			end
			%plot
			strPopSize = [num2str(intUseN) ' neurons; '];
			hFigI = figure;
			hAxI1 = subplot(2,2,1);ylabel('Bias-corrected Fisher information (d''^2)');xlabel('Population size (# of neurons)');
			title(['Unshuffled; ' strPopSize sprintf('ori diff: %.1fdeg; ',dblDiffTheta) strType],'Interpreter','none');hold on;
			hAxI2 = subplot(2,2,2);ylabel('Bias-corrected Fisher information (d''^2)');xlabel('Population size (# of neurons)');
			title(['Shuffled; area ' num2str(intWithinArea) '; ' strTag],'Interpreter','none');hold on;
			%full screen
			drawnow;
			jFig = get(handle(gcf), 'JavaFrame');
			jFig.setMaximized(true);
			figure(gcf);
			drawnow;
			% plot information
			matI_DirectPlot = squeeze(nanmean(matI_Direct,1));
			matI_shuffPlot = squeeze(nanmean(matI_Direct_shuff,1));
			
			%errorbar(hAxI1,vecGroupSizes,nanmean(matI_NCV(intCompareIdx,:,:),2),nanstd(matI(intCompareIdx,:,:),[],2),'Linestyle',':','Color',mapC(intCompareIdx,:));
			errorbar(hAxI1,vecGroupSizes,nanmean(matI_DirectPlot,1),nanstd(matI_DirectPlot,[],1));
			%errorbar(hAxI1,vecGroupSizes,nanmean(matI(intCompareIdx,:,:),2),nanstd(matI(intCompareIdx,:,:),[],2),'Color',mapC(intCompareIdx,:));
			
			%errorbar(hAxI2,vecGroupSizes,nanmean(matI_NCV_shuff(intCompareIdx,:,:),2),nanstd(matI(intCompareIdx,:,:),[],2),'Linestyle',':','Color',mapC(intCompareIdx,:));
			errorbar(hAxI2,vecGroupSizes,nanmean(matI_shuffPlot,1),nanstd(matI_shuffPlot,[],1));
			%errorbar(hAxI2,vecGroupSizes,nanmean(matI_Direct_shuff(intCompareIdx,:,:),2),nanstd(matI_Direct_shuff(intCompareIdx,:,:),[],2),'Linestyle','--','Color',mapC(intCompareIdx,:));
			drawnow;
			
			%set props
			axes(hAxI1);fixfig;ylim([0 (max(get(gca,'ylim')))]);
			axes(hAxI2);fixfig;ylim([0 (max(get(gca,'ylim')))]);
			drawnow;
			
			%% plot if # of trials < # of neurons
			if intUseRepetitionMax < max(vecGroupSizes)
			hAxI3 = subplot(2,2,3);ylabel('Bias-corrected Fisher information (d''^2)');xlabel('Population size (# of neurons)');
			title(['Unshuffled; ' strPopSize sprintf('ori diff: %.1fdeg; ',dblDiffTheta) strType],'Interpreter','none');hold on;
			hAxI4 = subplot(2,2,4);ylabel('Bias-corrected Fisher information (d''^2)');xlabel('Population size (# of neurons)');
			title(['Shuffled; area ' num2str(intWithinArea) '; ' strTag],'Interpreter','none');hold on;
			
			idxSub = intUseRepetitionMax > vecGroupSizes;
			errorbar(hAxI3,vecGroupSizes(idxSub),nanmean(matI_DirectPlot(:,idxSub),1),nanstd(matI_DirectPlot(:,idxSub),[],1));
			
			errorbar(hAxI4,vecGroupSizes(idxSub),nanmean(matI_DirectPlot(:,idxSub),1),nanstd(matI_shuffPlot(:,idxSub),[],1));
			drawnow;
			
			%set props
			axes(hAxI3);fixfig;ylim([0 (max(get(gca,'ylim')))]);
			axes(hAxI4);fixfig;ylim([0 (max(get(gca,'ylim')))]);
			drawnow;
			end
			%save figure
			if boolSaveFigs
				figure(hFigI);drawnow;
				strN = num2str(round(max(vecGroupSizes)));
				strSizeT = strcat('T',num2str(intUseRepetitionMax));
				export_fig([strFigDir  'AnalBlock1Info_N' strN strSizeT '_Area' num2str(intWithinArea) '_' strTag '.tif']);
				export_fig([strFigDir  'AnalBlock1Info_N' strN strSizeT '_Area' num2str(intWithinArea) '_' strTag '.pdf']);
			end
		end
		
		
		
		%% save data
		if boolSaveData
			strDataFile = ['D:\Data\Processed\V1_LIFmodel\Simulation_' strSimulation '_AB1_Area' num2str(intWithinArea) '.mat'];
			if exist(strDataFile,'file')
				sSave = load(strDataFile);
			else
				sSave = struct;
			end
			if boolAnalPred
				sSave.strParam = strParam;
				sSave.matPredictionsR2_SingleNeurons=matPredictionsR2_SingleNeurons;
				sSave.matPredictionsR2=matPredictionsR2;
				sSave.matPredDimDepR2=matPredDimDepR2;
				sSave.matPredictiveDimensions=matPredictiveDimensions;
				sSave.matTargetPopDimensionality=matTargetPopDimensionality;
				sSave.matR2RemPredDim=matR2RemPredDim;
				sSave.matR2UseDomDim=matR2UseDomDim;
				
				sSave.cellPredictionsR2_SingleNeurons=cellPredictionsR2_SingleNeurons;
				sSave.cellPredictionsR2=cellPredictionsR2;
				sSave.cellPredDimDepR2=cellPredDimDepR2;
				sSave.cellPredictiveDimensions=cellPredictiveDimensions;
				sSave.cellTargetPopDimensionality=cellTargetPopDimensionality;
				sSave.cellR2RemPredDim=cellR2RemPredDim;
				sSave.cellR2UseDomDim=cellR2UseDomDim;
				sSave.sParamsAnalSplit=sParamsAnalSplit;
			end
			if boolAnalI
				sSave.matI=matI_LR;
				sSave.matI_shuff=matI_shuff;
				sSave.cellI=cellI;
				sSave.cellI_shuff=cellI_shuff;
				sSave.sParamsAnal=sParamsAnal;
			end
			if boolAnalPopM
				%gather population-level metrics
				sSave.vecPopMax = vecPopMax;
				sSave.vecPopStDev =  vecPopStDev;
				sSave.vecPopMean = vecPopMean;
				sSave.vecPopHet = vecPopHet;
				
				sSave.matNoiseCorrsS1 = matNoiseCorrsS1;
				sSave.matNoiseCorrsS2 = matNoiseCorrsS2;
			end
			save(strDataFile,'-struct','sSave');
		end
	end
end
