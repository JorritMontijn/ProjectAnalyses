%% analyze input strength dependency for dimensionality of pop responses

%% initialize
close all;
clearvars;
boolLoad = false;
if ~exist('matModelResp','var')
	clearvars -except bool*;%close all;
	intLoadSim = 2; %% SET WHICH SIM TO LOAD HERE; 41 vs 42; 2 vs. 5
	boolLoad = true;
end

vecRunSims = 22;
for intLoadSim=vecRunSims
	clearvars -except vecRunSims intLoadSim
	boolLoad = true;
	if intLoadSim == 11 && boolLoad
		strSimulation = 'xAreaDistributed_OriDrift18Att0_2017-11-16';
	elseif intLoadSim == 12 && boolLoad
		strSimulation = 'xAreaDistributed_OriDrift18Att1_2017-11-16';
		
	elseif intLoadSim == 21 && boolLoad
		strSimulation = 'xAreaDistributed_Dyn5OriDrift2Att0_2018-01-10'; %new connectivity
	elseif intLoadSim == 211 && boolLoad
		boolLoadLGN = true;
		strSimulation = 'xAlsoLGN_Distrib_Dyn5OriDrift2Att0_2018-01-10'; %new connectivityelseif intLoadSim == 22 && boolLoad
	
	elseif intLoadSim == 22 && boolLoad
		strSimulation = 'xAreaDistributed_Noise0Drift2Att0_2018-01-16';
	elseif intLoadSim == 221 && boolLoad
		boolLoadLGN = true;
		strSimulation = 'xAlsoLGN_Distrib_Noise0Drift2Att0_2018-01-16';
	end
	
	%% RUN: #header
	strFigDir = 'D:\Data\Results\Block4\';
	if boolLoad
		runModelHeader;
	end
	cellIn = strsplit(strSimulation,'_');
	strFiller = cellIn{1};
	strType = cellIn{2};
	strDate = cellIn{3};
	strTag = [strType strDate];
	
	%% set parameters
	boolSaveFigs = true;
	boolSaveData = true;
	dblCutOff = 0.90;
	dblLambdaInfo = 1;
	dblLambdaPred = 0;
	dblNoiseLevel = 0;
	boolAnalInfoSub = true;
	intMaxDimAnal = 30;
	%intWithinArea = 1; %1, within V1; 2, within V2; 3, across V1-V2;
	intIters = 30;
	
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
	cellStrArea = {'V1 from V1','V2 from V2','V2 from V1'};
	
	%% add noise to improve numerical stability
	matModelResp = double(matModelResp);
	intNeurons = size(matModelResp,1);
	intTrials = size(matModelResp,2);
	vecNeuronSD = xstd(matModelResp,2);
	matNoise = normrnd(zeros(intNeurons,intTrials),repmat(vecNeuronSD*dblNoiseLevel,[1 intTrials]));
	matModelRespP = matModelResp + matNoise;
	%matModelRespP = zscore(matModelRespP,[],1);
	
	for intWithinArea=[1 2 3]
		%msg
		strPredArea = cellStrArea{intWithinArea};
		fprintf('Starting area %d; predicting %s [%s]\n',intWithinArea,strPredArea,getTime);
				
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
		
		
		%% pre-allocate variables dim dep
		if boolAnalInfoSub
			hFigDD = figure;
			%full screen
			drawnow;
			jFig = get(handle(gcf), 'JavaFrame');
			jFig.setMaximized(true);
			figure(gcf);
			drawnow;
		end
		
		%% set default parameters
		sParamsAnal = struct;
		sParamsAnal.dblLambda = dblLambdaInfo;
		sParamsAnal.boolDirectI = false;
		sParamsAnal.boolVerbose = false;
		sParamsAnal.dblDiffTheta = range(vecTrialOris);
		
		%% compare different stimulus intensity levels; from blank white noise to structured stimulus
		for intStimTypeIdx=1:intParamNum
			
			%% start
			dblParamVal = vecParamVals(intStimTypeIdx);
			cellC{intStimTypeIdx} = sprintf('%03d',dblParamVal);
			strC = cellC{intStimTypeIdx};
			fprintf('Starting %s %s (%d/%d) [%s]\n',strParam,strC,intStimTypeIdx,intParamNum,getTime);
			vecUseStimTypes = find(sort(vecParamIdx,'ascend')==intStimTypeIdx);
			vecUseStimTypes = [1 2];
			
			%check if no spikes
			if (sum(sum(matModelRespP(1:intCellsV1,vecTrialStimType==vecUseStimTypes(1)))) + sum(sum(matModelRespP(1:intCellsV1,vecTrialStimType==vecUseStimTypes(end))))) == 0
				boolAllZero = true;else boolAllZero = false;end
			
			%% dependence of dimensionality of stimulus intensity
			if boolAnalInfoSub
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
					
					%[matI_LogReg_bc_CV,sOut] = doFisherFull(cellMatX)
					
					%{
					%% get dimensionality with factor analysis
					[matSourceExplainedVarLL,matSourcePopDimensionality,matSourcePopLogLikeTest] = doDimFactorAnal(cellMatX,dblCutOff,intMaxDimAnal);
					fprintf('  C%s; Mean number of dominant dimensions in source: <%s\b> [%s]\n',strC,sprintf('%.2f ',xmean(matSourcePopDimensionality,1)),getTime);
					[matTargetExplainedVarLL,matTargetPopDimensionality,matTargetPopLogLikeTest] = doDimFactorAnal(cellMatY,dblCutOff,intMaxDimAnal);
					fprintf('  C%s; Mean number of dominant dimensions in target: <%s\b> [%s]\n',strC,sprintf('%.2f ',xmean(matTargetPopDimensionality,1)),getTime);
					%}	
				
					%% predict with reduced rank regression
					sPredAndI = doDimPredAndI_RR(cellMatX,cellMatY);
					
					%% extract
					vecRank = sPredAndI.vecRank;
					matPredFullR2 = sPredAndI.matPredFullR2;
					matPredDimDepR2 = sPredAndI.matPredDimDepR2;
					matPredDimDepR2Rem = sPredAndI.matPredDimDepR2Rem;
					matPredDimDepR2Rand = sPredAndI.matPredDimDepR2Rand;
					matFisherDimDepProj = sPredAndI.matFisherDimDepProj;
					matFisherDimDepOrth = sPredAndI.matFisherDimDepOrth;
					matFisherDimDepTotal = sPredAndI.matFisherDimDepTotal;
					
					matFisherDimDepTotalRand = sPredAndI.matFisherDimDepTotalRand;
					matFisherDimDepProjRand = sPredAndI.matFisherDimDepProjRand;
					matFisherDimDepOrthRand = sPredAndI.matFisherDimDepOrthRand;
					
					matFisherDimDepConsecProj = sPredAndI.matFisherDimDepConsecProj;
					matFisherDimDepConsecOrth = sPredAndI.matFisherDimDepConsecOrth;
					matFisherDimDepConsecProj2 = sPredAndI.matFisherDimDepConsecProj2;
					matFisherDimDepConsecOrth2 = sPredAndI.matFisherDimDepConsecOrth2;
					matFisherDimRedundancyPerc = sPredAndI.matFisherDimRedundancyPerc;
					matNoisePredR2ConsecProj = sPredAndI.matNoisePredR2ConsecProj;
					
					matFisherDimDepProjT = sPredAndI.matFisherDimDepProjT;
					matFisherDimDepOrthT = sPredAndI.matFisherDimDepOrthT;
					matFisherDimDepTotalT = sPredAndI.matFisherDimDepTotalT;
					
					matFisherDimDepTotalRandT = sPredAndI.matFisherDimDepTotalRandT;
					matFisherDimDepProjRandT = sPredAndI.matFisherDimDepProjRandT;
					matFisherDimDepOrthRandT = sPredAndI.matFisherDimDepOrthRandT;
					
					
					%transform
					matPredDimDepR2Avg = squeeze(mean(matPredDimDepR2,2));
					matPredDimDepR2RemAvg = squeeze(mean(matPredDimDepR2Rem,2));
					matPredDimDepR2RandAvg = squeeze(mean(mean(matPredDimDepR2Rand,4),2));
					matFisherDimDepProjAvg = squeeze(mean(matFisherDimDepProj,2));
					matFisherDimDepOrthAvg = squeeze(mean(matFisherDimDepOrth,2));
					matFisherDimDepTotalAvg = squeeze(mean(matFisherDimDepTotal,2));
					matFisherDimDepProjRandAvg = squeeze(mean(mean(matFisherDimDepProjRand,4),2));
					matFisherDimDepOrthRandAvg = squeeze(mean(mean(matFisherDimDepOrthRand,4),2));
					matFisherDimDepTotalRandAvg = squeeze(mean(mean(matFisherDimDepTotalRand,4),2));
					
					matFisherDimDepConsecProjAvg = squeeze(mean(matFisherDimDepConsecProj,2));
					matFisherDimDepConsecOrthAvg = squeeze(mean(matFisherDimDepConsecOrth,2));
					%matFisherDimDepConsecProjAvg2 = squeeze(mean(matFisherDimDepConsecProj2,2));
					%matFisherDimDepConsecOrthAvg2 = squeeze(mean(matFisherDimDepConsecOrth2,2));
					matFisherDimDepRedundancyAvg = squeeze(mean(matFisherDimRedundancyPerc,2));
					matPredConsecProjR2Avg = squeeze(mean(matNoisePredR2ConsecProj,2));
					
					matFisherDimDepProjAvgT = squeeze(mean(matFisherDimDepProjT,2));
					matFisherDimDepOrthAvgT = squeeze(mean(matFisherDimDepOrthT,2));
					matFisherDimDepTotalAvgT = squeeze(mean(matFisherDimDepTotalT,2));
					matFisherDimDepProjRandAvgT = squeeze(mean(mean(matFisherDimDepProjRandT,4),2));
					matFisherDimDepOrthRandAvgT = squeeze(mean(mean(matFisherDimDepOrthRandT,4),2));
					matFisherDimDepTotalRandAvgT = squeeze(mean(mean(matFisherDimDepTotalRandT,4),2));
					
					%% predict with full pop
					fprintf('  C%s; Mean full population predictions: <%s\b> [%s]\n',strC,sprintf('%.4f ',xmean(matPredFullR2,1)),getTime);
					%msg
					
					matPredictiveDimensions = sum(bsxfun(@rdivide,matPredDimDepR2,matFullPredictionsR2)<dblCutOff,3)+1;
					fprintf('  C%s; Mean number of predictive dimensions: <%s\b> [%s]\n',strC,sprintf('%.2f ',xmean(matPredictiveDimensions,1)),getTime);
					
					%% predict with removing predictive dimensions
					matR2RemPredDim = doDimPredRemPred(cellMatX,cellMatY,matPredictiveDimensions,matFullPredictionsR2,dblLambdaPred);
					fprintf('  C%s; Prediction before removal of pred dims: <%s\b>; After removal (d=%d): <%s\b> [%s]\n',strC,sprintf('%.4f ',xmean(matR2RemPredDim(:,:,1),1)),max(matPredictiveDimensions(:)),sprintf('%.4f ',xmean(matR2RemPredDim(:,:,end),1)),getTime);
					
					%% predict with removing dominant dimensions
					matR2UseDomDim = doDimPredRemDom(cellMatX,cellMatY,intMaxDimAnal,dblLambdaPred);
					fprintf('  C%s; Prediction using first dominant dim: <%s\b>; Using d=%d: <%s\b> [%s]\n',strC,sprintf('%.4f ',xmean(matR2UseDomDim(:,:,1),1)),intMaxDimAnal,sprintf('%.4f ',xmean(matR2UseDomDim(:,:,end),1)),getTime);
				
					%% do plotting
					if ishandle(hFigDD),figure(hFigDD);end
					clf;
					subplot(2,3,1)
					matFullPredCopied = repmat(mean(matPredFullR2,2),[1 numel(vecRank)]);
					hold on
					errorbar(vecRank,mean(matFullPredCopied),std(matFullPredCopied)/sqrt(intIters),'k-');
					errorbar(vecRank,mean(matPredDimDepR2Avg),std(matPredDimDepR2Avg)/sqrt(intIters),'b-');
					errorbar(vecRank,mean(matPredDimDepR2RemAvg),std(matPredDimDepR2RemAvg)/sqrt(intIters),'b--');
					errorbar(vecRank,mean(matPredDimDepR2RandAvg),std(matPredDimDepR2RandAvg)/sqrt(intIters),'r-');
					hold off
					h=legend({'Full model','Predictive subspace','Orth to pred subsp','Random subspace'});
					set(h,'Location','Best');
					title(sprintf('Predictability of %s activity (mean +/- SEM)',strPredArea));
					ylabel('Noise predictability (R^2)');
					xlabel('Dimensionality of subspace');
					xlim([min(vecRank)-1 max(vecRank+1)]);
					ylim([0 max(get(gca,'ylim'))]);
					fixfig; 
					
					subplot(2,3,2)
					%matFullMinusPred = matFisherDimDepTotalAvg - matFisherDimDepProjAvg;
					hold on
					errorbar(vecRank,mean(matFisherDimDepTotalAvg),std(matFisherDimDepTotalAvg)/sqrt(intIters),'k-');
					errorbar(vecRank,mean(matFisherDimDepOrthRandAvg),std(matFisherDimDepOrthRandAvg)/sqrt(intIters),'r--');
					errorbar(vecRank,mean(matFisherDimDepOrthAvg),std(matFisherDimDepOrthAvg)/sqrt(intIters),'b--');
					errorbar(vecRank,mean(matFisherDimDepProjAvg),std(matFisherDimDepProjAvg)/sqrt(intIters),'b-');
					errorbar(vecRank,mean(matFisherDimDepProjRandAvg),std(matFisherDimDepProjRandAvg)/sqrt(intIters),'r-');
					%errorbar(vecRank,mean(matFullMinusPred),std(matFullMinusPred)/sqrt(intIters),'m-');
					hold off
					%h2=legend({'Total','Orth to Rand Subsp','Orth to Pred Subsp','Pred Subspace','Rand Subspace','(Full) - (Pred Subsp)'});
					h2=legend({'Total','Orth to Rand Subsp','Orth to Pred Subsp','Pred Subspace','Rand Subspace'});
					set(h2,'Location','west');
					title('Fisher I in pred subsp of source pop (mean +/- SEM)');
					ylabel('Fisher information (d''^2)');
					xlabel('Dimensionality of subspace');
					xlim([min(vecRank)-1 max(vecRank+1)]);
					ylim([0 max(get(gca,'ylim'))]);
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
					title('Fisher I in pred subsp of target pop (mean +/- SEM)');
					ylabel('Fisher information (d''^2)');
					xlabel('Dimensionality of subspace');
					xlim([min(vecRank)-1 max(vecRank+1)]);
					ylim([0 max(get(gca,'ylim'))]);
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
					subplot(2,3,6)
					matPredPlusOrth = matFisherDimDepOrthAvg + matFisherDimDepProjAvg;
					hold on
					errorbar(vecRank,mean(matFisherDimDepTotalAvg),std(matFisherDimDepTotalAvg)/sqrt(intIters),'k-');
					errorbar(vecRank,mean(matPredPlusOrth),std(matPredPlusOrth)/sqrt(intIters),'b-');
					hold off
					h=legend({'I(Total)','I(Pred subsp) + I(Orth)'});
					set(h,'Location','Best');
					title('Redundancy Fisher info');
					ylabel('Fisher information (d''^2)');
					xlabel('Dimensionality of subspace');
					xlim([min(vecRank)-1 max(vecRank+1)]);
					ylim([0 max(get(gca,'ylim'))]);
					fixfig;
					
					subplot(2,3,4)
					hold on
					matPredConsecProjR2AvgCumSum = matPredConsecProjR2Avg;
					errorbar(vecRank,mean(matFullPredCopied),std(matFullPredCopied)/sqrt(intIters),'k-');
					errorbar(vecRank,mean(matPredConsecProjR2AvgCumSum),std(matPredConsecProjR2AvgCumSum)/sqrt(intIters),'b-');
					errorbar(vecRank,mean(matPredDimDepR2RandAvg),std(matPredDimDepR2RandAvg)/sqrt(intIters),'r-');
					hold off
					h=legend({'Full model','\Sigma^{-1} * f''','Random subspace'});
					set(h,'Location','Best');
					title('Prediction using consecutive coding dimensions');
					ylabel('Cumulative noise predictability (R^2)');
					xlabel('Optimal dimension #');
					xlim([0 10]);
					ylim([0 max(get(gca,'ylim'))]);
					fixfig;
					
					subplot(2,3,5)
					hold on
					errorbar(vecRank,mean(matFisherDimDepTotalAvg),std(matFisherDimDepTotalAvg)/sqrt(intIters),'k-');
					errorbar(vecRank,mean(matFisherDimDepConsecProjAvg),std(matFisherDimDepConsecProjAvg)/sqrt(intIters),'b-');
					%errorbar(vecRank,mean(matFisherDimDepConsecProjAvg2),std(matFisherDimDepConsecProjAvg2)/sqrt(intIters),'b--');
					hold off
					h=legend({'Total','Consec Proj'});
					set(h,'Location','Best');
					title('Fisher I with consecutive coding dimensions');
					ylabel('Fisher information (d''^2)');
					xlabel('Optimal dimension #');
					xlim([0 10]);
					ylim([0 max(get(gca,'ylim'))]);
					fixfig;
					
					%save figure
					if boolSaveFigs
						figure(hFigDD);drawnow;
						export_fig([strFigDir  'AnalBlock4Subsp_Area' num2str(intWithinArea) '_' strTag '.tif']);
						export_fig([strFigDir  'AnalBlock4Subsp_Area' num2str(intWithinArea) '_' strTag '.pdf']);
					end
				else
					
				end
				
			end
		end
		
		%% save data
		if boolSaveData
			strDataFile = ['D:\Data\Processed\V1_LIFmodel\Simulation_' strSimulation '_AB4_Area' num2str(intWithinArea) '.mat'];
			if exist(strDataFile,'file')
				sSave = load(strDataFile);
			else
				sSave = struct;
			end
			if boolAnalInfoSub
				sSave.strParam = strParam;
				sSave.cellMatX = cellMatX;
				sSave.cellMatY = cellMatY;
				sSave.sParamsAnalSplit = sParamsAnalSplit;
				sSave.matFullPredictionsR2 = matPredFullR2;
				sSave.sPredAndI = sPredAndI;
			end
			
			save(strDataFile,'-struct','sSave');
		end
	end
end
