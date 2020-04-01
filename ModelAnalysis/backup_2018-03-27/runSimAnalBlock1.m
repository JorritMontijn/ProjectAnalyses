%% analyze input strength dependency for dimensionality of pop responses

%% initialize
clearvars;
boolLoad = false;
if ~exist('matModelResp','var')
	clearvars -except bool*;%close all;
	intLoadSim = 42; %% SET WHICH SIM TO LOAD HERE; 41 vs 42; 2 vs. 5
	boolLoad = true;
end
boolLoadTrialSplit = true;
boolLoadSync = false;
boolLoadLGN = false;
boolLoadSpikeTimes = false;
vecRunSims = [-1 -2 11 12 13 14];
for intLoadSim=vecRunSims
	%if boolLoad
	clearvars -except vecRunSims intLoadSim bool*
	boolLoad = true;
	matUseStimTypes = [1 2];
	if intLoadSim == -1 && boolLoad
		strSimulation = 'xAreaExperiment_106r001p26-t_2018-02-21';
	elseif intLoadSim == -2 && boolLoad
		strSimulation = 'xAreaExperiment_107l003p143-t_2018-02-21';
		
	elseif intLoadSim == 11 && boolLoad
		strSimulation = 'xAreaDistributed_Ret32Col180Ori2Att0_2018-03-14';
	elseif intLoadSim == 12 && boolLoad
		strSimulation = 'xAreaDistributed_Ret32Col180Ori2Noise5_2018-03-14';
	elseif intLoadSim == 13 && boolLoad
		strSimulation = 'xAreaDistributed_Ret32Col180Ori8Noise0_2018-03-14';
	elseif intLoadSim == 14 && boolLoad
		strSimulation = 'xAreaDistributed_Ret32Col180Ori8Noise5_2018-03-23';
	end
	
	%% RUN: #header
	strFigDir = 'D:\Data\Results\Block1\';
	if boolLoad
		runAnalysisHeader;
	end
	
	for intDiffIdx = 1:size(matUseStimTypes,1)
		
		cellIn = strsplit(strSimulation,'_');
		strFiller = cellIn{1};
		strType = cellIn{2};
		strDate = cellIn{3};
		strTag = [strType '_' strDate];
		
		%% set parameters
		vecUseStimTypes = matUseStimTypes(intDiffIdx,:);
		dblDiffTheta = range(vecStimTypeOris(vecUseStimTypes));
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
		for intWithinArea=[1 2 3]
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
				intUseN = sum(vecCellArea==intWithinArea);
				vecGroupSizes = [2.^(0:10)];
				vecGroupSizes(vecGroupSizes>intUseN) = [];
				matI = nan(intParamNum,intIters,numel(vecGroupSizes));
				matI_shuff = nan(intParamNum,intIters,numel(vecGroupSizes));
				matI_Direct = nan(intParamNum,intIters,numel(vecGroupSizes));
				matI_Direct_shuff = nan(intParamNum,intIters,numel(vecGroupSizes));
				matI_NCV = nan(intParamNum,intIters,numel(vecGroupSizes));
				matI_NCV_shuff = nan(intParamNum,intIters,numel(vecGroupSizes));
						
				hFigI = figure;
				hAxI1 = subplot(2,2,1);ylabel('Bias-corrected Fisher information (d''^2)');xlabel('Population size (# of neurons)');
				title(['Unshuffled; ' strType]);hold on;
				hAxI2 = subplot(2,2,2);ylabel('Bias-corrected Fisher information (d''^2)');xlabel('Population size (# of neurons)');
				title(['Shuffled; area ' num2str(intWithinArea) sprintf('; ori diff: %.1fdeg',dblDiffTheta)]);hold on;
				%full screen
				drawnow;
				jFig = get(handle(gcf), 'JavaFrame');
				jFig.setMaximized(true);
				figure(gcf);
				drawnow;
			end
			%% pre-allocate variables dim dep
			if boolAnalPred
				hFigDD = figure;
				hAxD1 = subplot(2,3,1);
				xlabel('Number of predictive dimensions');
				ylabel('Predictive performance');
				title(['Semedo Fig 4;' strType]);
				hold on;
				
				hAxD2 = subplot(2,3,2);
				plot([0 10],[0 10],'k--');
				hold on;
				xlabel('Target population dimensionality')
				ylabel('Number of predictive dimensions')
				
				hAxD3 = subplot(2,3,3);
				xlabel('Number of dominant dimensions');
				ylabel('Predictive performance');
				hold on;
				%full screen
				drawnow;
				jFig = get(handle(gcf), 'JavaFrame');
				jFig.setMaximized(true);
				figure(gcf);
				drawnow;
			end
			%% pre-allocate variables mean
			if boolAnalPopM && (intWithinArea == 1 || intWithinArea == 2)
				%#ok<*SAGROW> %set tag to disable pre-allocation warning, because matlab is apparently retarded
				intUseN = sum(vecCellArea==intWithinArea);
				matNoiseCorrsS1 = nan(intUseN,intUseN,intParamNum);
				matNoiseCorrsS2 = nan(intUseN,intUseN,intParamNum);
				matNC_Diffs = nan(intUseN,intUseN,intParamNum);
				matTril = tril(true(intUseN),-1);
				intDistroVals = sum(matTril(:));
				matNC_Distro1 = nan(intDistroVals,intParamNum);
				matNC_Distro2 = nan(intDistroVals,intParamNum);
				
				matAct = matModelRespP(vecCellArea==intWithinArea,:);
				vecActivity = squeeze(mean(matAct,1));
				
				%if ~exist('cellHet','var') || numel(cellHet) < intWithinArea
				%	[vecHeterogeneity,vecActivity] = calcMatRespHeteroGen(matAct);
				%	cellHet{intWithinArea} = vecHeterogeneity;
				%	cellMean{intWithinArea} = vecActivity;
				%	save(['D:\Data\Processed\V1_LIFmodel\Simulation_' strSimulation '_prepro2.mat'],'cellHet','cellMean');
				%else
				%	vecHeterogeneity = cellHet{intWithinArea};
				%	vecActivity = cellMean{intWithinArea};
				%end
				vecPopMax = nan(1,intParamNum);
				vecPopStDev = nan(1,intParamNum);
				vecPopMean = nan(1,intParamNum);
				vecPopHet = nan(1,intParamNum);
				
				
				%prepare figures
				%hFigPM1 = figure;
				
				%full screen
				drawnow;
				jFig = get(handle(gcf), 'JavaFrame');
				jFig.setMaximized(true);
				figure(gcf);
				drawnow;
				
				%fig 2
				hFigPM = figure;
				hAxM1 = subplot(2,3,1);
				hold on;
				hAxM2 = subplot(2,3,2);
				hold on;
				hAxM3 = subplot(2,3,3);
				hold on;
				hAxM4 = subplot(2,3,4);
				hold on;
				
				%full screen
				drawnow;
				jFig = get(handle(gcf), 'JavaFrame');
				jFig.setMaximized(true);
				figure(gcf);
				drawnow;
				
			end
			
			%% compare different stimulus intensity levels; from blank white noise to structured stimulus
			for intStimTypeIdx=1:intParamNum
				
				%% start
				dblParamVal = vecParamVals(intStimTypeIdx);
				cellC{intStimTypeIdx} = sprintf('%03d',dblParamVal);
				strC = cellC{intStimTypeIdx};
				fprintf('Starting %s %s (%d/%d) [%s]\n',strParam,strC,intStimTypeIdx,intParamNum,getTime);
				if ~exist('vecUseStimTypes','var')
					vecUseStimTypes = find(sort(vecParamIdx,'ascend')==intStimTypeIdx);
				end
				
				%check if no spikes
				if (sum(sum(matModelRespP(1:intCellsV1,vecTrialStimType==vecUseStimTypes(1)))) + sum(sum(matModelRespP(1:intCellsV1,vecTrialStimType==vecUseStimTypes(end))))) == 0
					boolAllZero = true;else boolAllZero = false;end
				
				%% dependence of Fisher information on stimulus intensity
				if boolAnalI && (intWithinArea == 1 || intWithinArea == 2)
					%msg
					fprintf('  C%s; running Fisher information analysis [%s]\n',strC,getTime);
					
					%% calculate Fisher information on input LGN units
					if exist('matModelRespLGN_ON','var')
						%select trials
						indUseTrials = ismember(vecTrialStimType,vecUseStimTypes);
						vecUseTrialStimType = vecTrialStimType(indUseTrials);
						dblDiffTheta = range(vecStimTypeOris(vecUseStimTypes));
			
						%get data
						matModelRespLGN = cat(1,matModelRespLGN_ON(:,indUseTrials),matModelRespLGN_OFF(:,indUseTrials));
						[dblPredA,matPredA,dblD2,dblD2mat,dblD2_diag] = getSeparation(matModelRespLGN,vecUseTrialStimType,0,dblDiffTheta);
						
						intGroupSize = size(matModelRespLGN,1);
						intTrials12 = size(matModelRespLGN,2); %check if this one
						dblSubFac =(2*intGroupSize)/(intTrials12*(range(vecTrialOris).^2));
						dblProdFacRaw = ((2*intTrials12-intGroupSize-3)/(2*intTrials12-2));
						dblInputI = dblD2*dblProdFacRaw-dblSubFac;
						if ~boolBiasCorrection,dblInputI = dblD2;end %not bias-corrected
						
						%%
						%{
						%calc log reg
						sParamsInfoLGN = struct;
						sParamsInfoLGN.vecGroupSizes = intGroupSize;
						sParamsInfoLGN.vecUseStimTypes = vecUseStimTypes;
						sParamsInfoLGN.dblLambda = dblLambdaInfo;
						sParamsInfoLGN.boolDirectI = true;
						sParamsInfoLGN.intIters = 3;%intIters;
						sParamsInfoLGN.boolVerbose = false;
						sParamsInfoLGN.dblDiffTheta = range(vecTrialOris(vecUseStimTypes));
						sParamsInfoLGN.boolBiasCorrection = boolBiasCorrection;
						[matFisherFull,sAgg] = doFisherAnal(matModelRespLGN,vecUseTrialStimType,sParamsInfoLGN);
						
						dblI_CV = mean(sAgg.matI_LogReg_bc_CV(:));
						dblI = mean(sAgg.matI_LogReg_bc(:));
						dblI_Dir = mean(sAgg.matI_Direct_bc(:));
						
						
						%drifting: CV, -0.1; non-CV, 2.6
						%}
						%return
					else
						if isfield(sData,'vecTrialOriNoise')
							dblNoise = 0;%mean(sData.vecTrialOriNoise(:));
							dblInputI = 0;%(range(vecTrialOris) / sqrt(0.5*(dblNoise^2 + dblNoise^2))).^2;
						else
							dblNoise = 0;
							dblInputI = 0;
						end
					end
					% HIERO!!
					
					%% calculate information in cortex
					%set parameters
					sParamsAnal = struct;
					sParamsAnal.vecUseStimTypes = vecUseStimTypes;
					sParamsAnal.dblLambda = dblLambdaInfo;
					sParamsAnal.boolDirectI = true;
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
						matI(intStimTypeIdx,:) = 0;
						matI_shuff(intStimTypeIdx,:) = 0;
					else
						
						%do analysis
						[matFisherFull,sOut] = doFisherAnal(matModelRespP(vecCells,:),vecTrialStimType,sParamsAnal);
						
						%get output
						cellI{intStimTypeIdx} = nanmean(sOut.matI_LogReg_bc_CV,3); %mean across folds
						cellI_shuff{intStimTypeIdx} = nanmean(sOut.matI_LogReg_bc_CV_shuff,3); %mean across folds
						matI(intStimTypeIdx,:,:) = nanmean(sOut.matI_LogReg_bc_CV,3); %mean across folds and resmplings
						matI_shuff(intStimTypeIdx,:,:) = nanmean(sOut.matI_LogReg_bc_CV_shuff,3); %mean across folds and resmplings
						matI_NCV(intStimTypeIdx,:,:) = nanmean(sOut.matI_LogReg_bc,3); %mean across folds and resmplings
						matI_NCV_shuff(intStimTypeIdx,:,:) = nanmean(sOut.matI_LogReg_bc_shuff,3); %mean across folds and resmplings
						
						matI_Direct(intStimTypeIdx,:,:) = nanmean(sOut.matI_Direct_bc,3); %mean across folds and resmplings
						matI_Direct_shuff(intStimTypeIdx,:,:) = nanmean(sOut.matI_Direct_bc_shuff,3); %mean across folds and resmplings
						
						cellI_NCV{intStimTypeIdx} = nanmean(sOut.matI_LogReg_bc,3); %mean across folds
						cellI_NCV_shuff{intStimTypeIdx} = nanmean(sOut.matI_LogReg_bc_shuff,3); %mean across folds
					end
					
					
					% plot information
					errorbar(hAxI1,vecGroupSizes,nanmean(matI_NCV(intStimTypeIdx,:,:),2),nanstd(matI(intStimTypeIdx,:,:),[],2),'Linestyle',':','Color',mapC(intStimTypeIdx,:));
					errorbar(hAxI1,vecGroupSizes,nanmean(matI_Direct(intStimTypeIdx,:,:),2),nanstd(matI_Direct(intStimTypeIdx,:,:),[],2),'Linestyle','--','Color',mapC(intStimTypeIdx,:));
					errorbar(hAxI1,vecGroupSizes,nanmean(matI(intStimTypeIdx,:,:),2),nanstd(matI(intStimTypeIdx,:,:),[],2),'Color',mapC(intStimTypeIdx,:));
					plot(hAxI1,vecGroupSizes,ones(size(vecGroupSizes))*dblInputI,'k--');
					legend({'LR','Direct bc','LR CV'})
					
					errorbar(hAxI2,vecGroupSizes,nanmean(matI_NCV_shuff(intStimTypeIdx,:,:),2),nanstd(matI(intStimTypeIdx,:,:),[],2),'Linestyle',':','Color',mapC(intStimTypeIdx,:));
					errorbar(hAxI2,vecGroupSizes,nanmean(matI_shuff(intStimTypeIdx,:,:),2),nanstd(matI_shuff(intStimTypeIdx,:,:),[],2),'Color',mapC(intStimTypeIdx,:));
					errorbar(hAxI2,vecGroupSizes,nanmean(matI_Direct_shuff(intStimTypeIdx,:,:),2),nanstd(matI_Direct_shuff(intStimTypeIdx,:,:),[],2),'Linestyle','--','Color',mapC(intStimTypeIdx,:));
					drawnow;
				end
				
				%% dependence of dimensionality of stimulus intensity
				if boolAnalPred
					if ~boolAllZero
						%% get data splits
						%set parameters
						sParamsAnalSplit=struct;
						sParamsAnalSplit.intSizeX = 110;
						sParamsAnalSplit.intSizeY = 30;
						sParamsAnalSplit.intResamplings = intIters;
						sParamsAnalSplit.vecCellArea = vecCellArea;
						sParamsAnalSplit.intWithinArea = intWithinArea;
						sParamsAnalSplit.vecUseStimTypes = vecUseStimTypes;
						%get splits
						[cellMatX,cellNeuronsX,cellMatY,cellNeuronsY] = doDimDataSplits(matModelRespP,vecTrialStimType,sParamsAnalSplit);
						fprintf('  C%s; Created %dx%d data splits [%s]\n',strC,size(cellMatX),getTime);
						
						%% predict
						%% predict with single neuron
						matPredictionsR2_SingleNeurons = doDimPredSingle(cellMatX,cellMatY,dblLambdaPred);
						fprintf('  C%s; Mean single neuron predictions: <%s\b> [%s]\n',strC,sprintf('%.4f ',xmean(xmean(matPredictionsR2_SingleNeurons,3),1)),getTime);
						
						%% predict with full pop
						matPredictionsR2 = doDimPredFull(cellMatX,cellMatY,dblLambdaPred);
						fprintf('  C%s; Mean full population predictions: <%s\b> [%s]\n',strC,sprintf('%.4f ',xmean(matPredictionsR2,1)),getTime);
						
						
						%% predict with RRR
						matPredDimDepR2 = doDimPredRR(cellMatX,cellMatY);
						matPercOfTotPred = bsxfun(@rdivide,matPredDimDepR2,matPredictionsR2);
						matPredictiveDimensions = sum(matPercOfTotPred<0.95,3);
						fprintf('  C%s; Mean number of predictive dimensions: <%s\b> [%s]\n',strC,sprintf('%.2f ',xmean(matPredictiveDimensions,1)),getTime);
						
						%% get dimensionality with factor analysis
						%[matSourcePopLogLikeTest,matSourcePopDimensionality] = doDimFactorAnal(cellMatX,15);
						[matTargetPopLogLikeTest,matTargetPopDimensionality] = doDimFactorAnal(cellMatY,intMaxDimAnal);
						fprintf('  C%s; Mean number of dominant dimensions: <%s\b> [%s]\n',strC,sprintf('%.2f ',xmean(matTargetPopDimensionality,1)),getTime);
						
						%% predict with removing predictive dimensions
						matR2RemPredDim = doDimPredRemPred(cellMatX,cellMatY,matPredictiveDimensions,matPredictionsR2,dblLambdaPred);
						fprintf('  C%s; Prediction before removal of pred dims: <%s\b>; After removal (d=%d): <%s\b> [%s]\n',strC,sprintf('%.4f ',xmean(matR2RemPredDim(:,:,1),1)),max(matPredictiveDimensions(:)),sprintf('%.4f ',xmean(matR2RemPredDim(:,:,end),1)),getTime);
						
						%% predict with removing dominant dimensions
						matR2UseDomDim = doDimPredRemDom(cellMatX,cellMatY,intMaxDimAnal,dblLambdaPred);
						fprintf('  C%s; Prediction using first dominant dim: <%s\b>; Using d=%d: <%s\b> [%s]\n',strC,sprintf('%.4f ',xmean(matR2UseDomDim(:,:,1),1)),intMaxDimAnal,sprintf('%.4f ',xmean(matR2UseDomDim(:,:,end),1)),getTime);
						
					else
						matPredictionsR2_SingleNeurons=0;
						matPredictionsR2=0;
						matPredDimDepR2=zeros(1,1,intMaxDimAnal);
						matPredictiveDimensions=0;
						matTargetPopDimensionality=0;
						matR2RemPredDim=0;
						matR2UseDomDim=zeros(1,1,intMaxDimAnal);
						%matI
						%matI_shuff
						%sParamsAnal
						%sParamsAnalSplit
					end
					%% do plotting
					%fig 4
					axes(hAxD1);
					vecPredDD = squeeze(nanmean(nanmean(matPredDimDepR2,2),1));
					plot(1:length(vecPredDD),vecPredDD,'Color',mapC(intStimTypeIdx,:));
					scatter(length(vecPredDD)-0.5,nanmean(matPredictionsR2(:)),80,mapC(intStimTypeIdx,:),'x');
					%if intContrastIdx==1,text(1,mean(matPredictionsR2(:))+0.02,'Full prediction','Color',mapC(intContrastIdx,:),'FontSize',14,'Rotation',45);end
					
					%fig 5
					axes(hAxD2);
					dblDomDim = mean(matTargetPopDimensionality(:));
					dblPredDim = mean(matPredictiveDimensions(:));
					scatter(dblDomDim,dblPredDim,80,mapC(intStimTypeIdx,:),'x');
					
					%fig 7
					axes(hAxD3);
					vecPlot = (1:intMaxDimAnal);
					vecDomDD =  squeeze(nanmean(nanmean(matR2UseDomDim,2),1));
					%plot(vecPlot,vecPredDD(vecPlot)','Color',mapC(intStimTypeIdx,:),'Marker','x');
					plot(vecPlot,vecDomDD(vecPlot)','Color',mapC(intStimTypeIdx,:));%,'Marker','o');
					drawnow;
					
					%% save data
					cellPredictionsR2_SingleNeurons{intStimTypeIdx}=matPredictionsR2_SingleNeurons;
					cellPredictionsR2{intStimTypeIdx}=matPredictionsR2;
					cellPredDimDepR2{intStimTypeIdx}=matPredDimDepR2;
					cellPredictiveDimensions{intStimTypeIdx}=matPredictiveDimensions;
					cellTargetPopDimensionality{intStimTypeIdx}=matTargetPopDimensionality;
					cellR2RemPredDim{intStimTypeIdx}=matR2RemPredDim;
					cellR2UseDomDim{intStimTypeIdx}=matR2UseDomDim;
					
				end
				if boolAnalPopM && (intWithinArea == 1 || intWithinArea == 2)
					%% prep
					%vecUseStimTypes = [1 2]+20
					indS1 = vecTrialStimType == vecUseStimTypes(1);
					indS2 = vecTrialStimType == vecUseStimTypes(end);
					indBoth = indS1 | indS2;
					matThisAct = matAct(:,indBoth);
					
					%gather population-level metrics
					vecPopMax(intStimTypeIdx) = mean(max(matAct(:,indBoth),[],1));
					vecPopStDev(intStimTypeIdx) = mean(xstd(matAct(:,indBoth),1));
					vecPopMean(intStimTypeIdx) = mean(vecActivity(indBoth));
					%vecPopHet(intStimTypeIdx) = mean(vecHeterogeneity(indBoth));
					
					matNoiseCorr1 = corr(matAct(:,indS1)');
					matNoiseCorr2 = corr(matAct(:,indS2)');
					matNoiseDiff = (matNoiseCorr2-matNoiseCorr1);
					
					%save
					matNoiseCorrsS1(:,:,intStimTypeIdx) = matNoiseCorr1;
					matNoiseCorrsS2(:,:,intStimTypeIdx) = matNoiseCorr2;
					matNC_Distro1(:,intStimTypeIdx) = matNoiseCorr1(matTril);
					matNC_Distro2(:,intStimTypeIdx) = matNoiseCorr2(matTril);
					
					vecNC_StDev(intStimTypeIdx) = (std(matNC_Distro1(:,intStimTypeIdx)) + std(matNC_Distro2(:,intStimTypeIdx)))/2;
					vecNC_Mean(intStimTypeIdx) = (mean(matNC_Distro1(:,intStimTypeIdx)) + mean(matNC_Distro2(:,intStimTypeIdx)))/2;
					
					matNC_Diffs(:,:,intStimTypeIdx) = matNoiseDiff./vecNC_StDev(intStimTypeIdx);
					
					%calculate noise correlations as a function of tuning similariy
					indOriTuned = ~isnan(vecPrefOri(vecCellArea==intWithinArea));
					if ~isempty(indOriTuned) && ~all(~indOriTuned)
						vecNeuronsPrefs = vecPrefOri(indOriTuned)*2;
						matPrefOriDists = roundi(abs(bsxfun(@circ_dist,vecNeuronsPrefs,vecNeuronsPrefs')),5);
						vecDiffVals = unique(matPrefOriDists);
						matPrefOriDists(diag(true(1,numel(vecNeuronsPrefs)),0)) = nan;
						
						intVals = numel(vecDiffVals);
						vecMeanCorr = nan(1,intVals);
						vecSDCorr = nan(1,intVals);
						vecCountCorr = nan(1,intVals);
						for intDiff=1:intVals
							dblDiffVal = vecDiffVals(intDiff);
							vecCorrVals = matNoiseCorrsS1(matPrefOriDists==dblDiffVal);
							vecMeanCorr(intDiff) = mean(vecCorrVals);
							vecSDCorr(intDiff) = std(vecCorrVals);
							vecCountCorr(intDiff) = numel(vecCorrVals);
						end
						axes(hAxM4)
						errorbar(vecDiffVals,vecMeanCorr,vecSDCorr./sqrt(vecCountCorr));;
					end
					%%
					%{
					figure(hFigPM1)
					subplot(3,4,intStimTypeIdx)
					imagesc((matNoiseCorr1),[-1 1])
					colorbar
					colormap(redblue)
					%}
					%%
					[N,EDGES] = histcounts(matNC_Distro1(:,intStimTypeIdx));
					figure(hFigPM)
					plot(hAxM1,EDGES,[0 (N/sum(N))./diff(EDGES)],'Color',mapC(intStimTypeIdx,:));
					
					
				end
				
				
			end
			if boolAnalI && (intWithinArea == 1 || intWithinArea == 2)
				%set props
				axes(hAxI1);fixfig;ylim([0 (max(get(gca,'ylim')))]);
				axes(hAxI2);fixfig;ylim([0 (max(get(gca,'ylim')))]);
				drawnow;
				
				%save figure
				if boolSaveFigs
					figure(hFigI);drawnow;
					export_fig([strFigDir  'AnalBlock1Info_Area' num2str(intWithinArea) '_' strTag '.tif']);
					export_fig([strFigDir  'AnalBlock1Info_Area' num2str(intWithinArea) '_' strTag '.pdf']);
				end
			end
			if boolAnalPred
				%set props
				axes(hAxD1);fixfig;
				axes(hAxD2);fixfig;legend([{''};cellC(:)],'Location','northwest');
				axes(hAxD3);fixfig;legend({'Predictive','Dominant'},'Location','Best')
				drawnow;
				
				%save figure
				if boolSaveFigs
					figure(hFigDD);drawnow;
					export_fig([strFigDir  'AnalBlock1Pred_Area' num2str(intWithinArea) '_' strTag '.tif']);
					export_fig([strFigDir  'AnalBlock1Pred_Area' num2str(intWithinArea) '_' strTag '.pdf']);
				end
			end
			if boolAnalPopM && (intWithinArea == 1 || intWithinArea == 2)
				%set props
				axes(hAxM1);
				ylabel('Normalized density');
				xlabel('Noise correlation (r)');fixfig;
				
				axes(hAxM2);
				plot(hAxM2,vecParamVals,vecNC_StDev,'kx');
				hold on
				plot(hAxM2,vecParamVals,vecNC_Mean,'rx');
				if numel(vecParamVals) > 1,xlim([min(vecParamVals) max(vecParamVals)]);end
				hold off
				legend('SD','Mean')
				xlabel('Stimulation Intensity')
				ylabel('Noise correlation (r)')
				fixfig;
				
				axes(hAxM3);
				plot(vecParamVals,vecPopMax,'x');
				xlabel('Stimulation Intensity')
				ylabel('Maximum firing rate (Hz)')
				fixfig;
				
				axes(hAxM4);
				xlabel('Difference in pref ori')
				ylabel('Noise correlation')
				fixfig;
				
				subplot(2,3,5)
				plot(vecParamVals,vecPopStDev,'kx');
				hold on
				plot(vecParamVals,vecPopMean,'rx');
				hold off
				legend('SD','Mean')
				xlabel('Stimulation Intensity')
				ylabel('Population firing rate (Hz)')
				ylim([0 max(get(gca,'ylim'))]);
				fixfig;
				
				subplot(2,3,6)
				vecCV = vecPopStDev./vecPopMean;
				vecCV(abs(vecCV)>10)=nan;
				plot(vecParamVals,vecCV,'x');
				xlabel('Stimulation Intensity')
				ylabel('Population CV')
				fixfig;
				
				
				
				
				%drawnow;
				
				%save figure
				if boolSaveFigs
					figure(hFigPM);drawnow;
					export_fig([strFigDir  'AnalBlock1PopM_Area' num2str(intWithinArea) '_' strTag '.tif']);
					export_fig([strFigDir  'AnalBlock1PopM_Area' num2str(intWithinArea) '_' strTag '.pdf']);
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
					sSave.matI=matI;
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
end
