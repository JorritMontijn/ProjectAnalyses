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
boolLoadSync = false;
boolLoadLGN = false;
boolLoadSpikeTimes = true;
vecRunSims = [35];
for intLoadSim=vecRunSims
	clearvars -except vecRunSims intLoadSim bool*
	boolLoad = true;
	if intLoadSim == 10 && boolLoad
		strSimulation = 'xAreaDistributed_Ret16Noise0Ori2Drift2Att0_2018-01-30';
	elseif intLoadSim == 11 && boolLoad
		strSimulation = 'xAreaDistributed_Ret16Noise1Ori2Drift2Att0_2018-02-05';
	elseif intLoadSim == 13 && boolLoad
		strSimulation = 'xAreaDistributed_Ret16Noise3Ori2Drift2Att0_2018-02-05';
	elseif intLoadSim == 15 && boolLoad
		strSimulation = 'xAreaDistributed_Ret16Noise5Ori2Drift2Att0_2018-02-05';
		
		
		
	elseif intLoadSim == 20 && boolLoad
		strSimulation = 'xAreaDistributed_Ret32Noise0Ori2Drift2Att0_2018-01-29'; %new connectivity
	elseif intLoadSim == 25 && boolLoad
		strSimulation = 'xAreaDistributed_Ret32Noise5Ori2Drift2Att0_2018-02-01';
		
	elseif intLoadSim == 30 && boolLoad
		strSimulation = '';
	elseif intLoadSim == 35 && boolLoad
		strSimulation = 'xAreaDistributed_Ret64Noise5Ori2Drift2Att0_2018-02-08';
	end
	
	%% RUN: #header
	strFigDir = 'D:\Data\Results\Block4\';
	if boolLoad
		runAnalysisHeader;
	end
	cellIn = strsplit(strSimulation,'_');
	strFiller = cellIn{1};
	strType = cellIn{2};
	strDate = cellIn{3};
	strTag = [strType strDate];
	
	%% set parameters
	boolSaveFigs = true;
	boolSaveData = true;
	dblLambdaInfo = 1;
	dblLambdaPred = 0;
	dblNoiseLevel = 0;
	boolAnalInfoSub = true;
	intMaxDimAnal = 30;
	%intWithinArea = 1; %1, within V1; 2, within V2; 3, across V1-V2;
	intIters = 50;
	
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
		
		
		%% compare different stimulus intensity levels; from blank white noise to structured stimulus
		for intStimTypeIdx=1:intParamNum
			
			%% start
			dblParamVal = vecParamVals(intStimTypeIdx);
			cellC{intStimTypeIdx} = sprintf('%03d',dblParamVal);
			strC = cellC{intStimTypeIdx};
			fprintf('Starting %s %s (%d/%d) [%s]\n',strParam,strC,intStimTypeIdx,intParamNum,getTime);
			%vecUseStimTypes = find(sort(vecParamIdx,'ascend')==intStimTypeIdx);
			vecUseStimTypes = [1 2];
			indUseTrials = ismember(vecTrialStimType,vecUseStimTypes);
			
			%% set default parameters
			dblDiffTheta = range(vecTrialOris(vecUseStimTypes));
			sParamsAnal = struct;
			sParamsAnal.dblLambda = dblLambdaInfo;
			sParamsAnal.boolDirectI = false;
			sParamsAnal.boolVerbose = false;
			sParamsAnal.dblDiffTheta = dblDiffTheta;

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
					
					%%
					intCoincN = 30;
					vecNeurons = randperm(intCellsV1,intCoincN);
					cellSpikeTimes = cellSpikeTimesCortex(vecNeurons);
					vecNumObs = round(2.^[6:14]);
					[sMetricsOut, vecNumObs] = testCoincTrialDep(cellSpikeTimes,vecTrialStimType,vecStimStartSecs,vecStimStopSecs,vecNumObs,dblDiffTheta);

					%% predict and I with RRR; sync-based
					intCoincN = 30;
					vecNeurons = randperm(intCellsV1,intCoincN);
					cellSpikeTimes = cellSpikeTimesCortex(vecNeurons);
					intUseTrials = (intCoincN*2)+500;
					dblStopT = vecStimStopSecs(intUseTrials);
					for intN=1:intCoincN
						cellSpikeTimes{intN} = cellSpikeTimes{intN}(cellSpikeTimes{intN} < dblStopT);
					end
					
					%transform spike times to integers
					boolUseInt = true;
					hTic = tic;
					fprintf('Preparing coinc cort calc... [%s]\n',getTime);
					if boolUseInt && ~all(isint(cellSpikeTimes{1}))
						cellSpikeTimesInt = doSpikeDoubleToInt(cellSpikeTimes,dblStepSize,0);
						vecStimStartInt = doSpikeDoubleToInt(vecStimStartSecs,dblStepSize,0);
						vecStimStopInt = doSpikeDoubleToInt(vecStimStopSecs,dblStepSize,0);
					else
						cellSpikeTimesInt = cellSpikeTimes;
						vecStimStartInt = vecStimStartSecs;
						vecStimStopInt = vecStimStopSecs;
					end
					
					matAggMeanCoincCort = nan(intCoincN,intUseTrials);
					matAggDiagCoincCort = nan(intCoincN,intUseTrials);
					matAggAllCoincCort = nan(intCoincN,intCoincN,intUseTrials);
					parfor intTrial=1:intUseTrials
						vecWindow = [vecStimStartInt(intTrial) vecStimStopInt(intTrial)];
						matCoincidence = calcCoincidence(cellSpikeTimesInt,vecWindow);
						matAggMeanCoincCort(:,intTrial) = mean(matCoincidence);
						matAggDiagCoincCort(:,intTrial) = diag(matCoincidence);
						matAggAllCoincCort(:,:,intTrial) = matCoincidence;
						%if toc(hTic) > 5 || intTrial==1
						%	hTic = tic;
						%	fprintf('Coinc cort calc; Now at trial %d/%d [%s]\n',intTrial,intTrials,getTime);
						%end
					end
					toc(hTic);
					
					%% get bias correction & prep
					vecUseTrialStimType = vecTrialStimType(1:intUseTrials);
					matMuAct = matModelRespP(vecNeurons,1:intUseTrials);
					matMeanCoincCort = matAggMeanCoincCort(:,1:intUseTrials);
					matDiagCoincCort = matAggDiagCoincCort(:,1:intUseTrials);
					matAllCoincCort = matAggAllCoincCort(:,1:intUseTrials);
					intRepetitions = intUseTrials/2;
					
					%% 30 predictors
					intPredictors = 30;
					dblSubFac =(2*intPredictors)/(intRepetitions*(dblDiffTheta.^2));
					dblProdFacRaw = ((2*intRepetitions-intPredictors-3)/(2*intRepetitions-2));
					[dblPredA,matPredA,dblI_Rate,dblImat,dblI_diag] = getSeparation(matMuAct,vecUseTrialStimType,0,dblDiffTheta);
					[dblPredA,matPredA,dblI_PopCoinc,dblImat,dblI_diag] = getSeparation(matMeanCoincCort,vecUseTrialStimType,0,dblDiffTheta);
					[dblPredA,matPredA,dblI_SelfCoinc,dblImat,dblI_diag] = getSeparation(matDiagCoincCort,vecUseTrialStimType,0,dblDiffTheta);
					dblI_Rate_bc = dblI_Rate*dblProdFacRaw-dblSubFac;
					dblI_PopCoinc_bc = dblI_PopCoinc*dblProdFacRaw-dblSubFac;
					dblI_SelfCoinc_bc = dblI_SelfCoinc*dblProdFacRaw-dblSubFac;
					
					%% 30 predictors, composites
					matBurstiness = matDiagCoincCort./matMuAct;
					[dblPredA,matPredA,dblI_Burst,dblImat,dblI_diag] = getSeparation(matBurstiness,vecUseTrialStimType,0,dblDiffTheta);
					
					matSynchronicity = matMeanCoincCort./matMuAct;
					[dblPredA,matPredA,dblI_Sync,dblImat,dblI_diag] = getSeparation(matSynchronicity,vecUseTrialStimType,0,dblDiffTheta);
					
					matSimultaneity = matMeanCoincCort./matDiagCoincCort;
					[dblPredA,matPredA,dblI_Simult,dblImat,dblI_diag] = getSeparation(matSimultaneity,vecUseTrialStimType,0,dblDiffTheta);
					dblI_Burst_bc = dblI_Burst*dblProdFacRaw-dblSubFac;
					dblI_Sync_bc = dblI_Sync*dblProdFacRaw-dblSubFac;
					dblI_Simult_bc = dblI_Simult*dblProdFacRaw-dblSubFac;
					
					%% 60 predictors
					intPredictors = 60;
					dblSubFac =(2*intPredictors)/(intRepetitions*(dblDiffTheta.^2));
					dblProdFacRaw = ((2*intRepetitions-intPredictors-3)/(2*intRepetitions-2));
					
					matRateAndPopCoinc = cat(1,matMeanCoincCort,matMuAct);
					[dblPredA,matPredA,dblI_RateAndPopCoinc,dblImat,dblI_diag] = getSeparation(matRateAndPopCoinc,vecUseTrialStimType,0,dblDiffTheta);
					
					matRateAndSelfCoinc = cat(1,matDiagCoincCort,matMuAct);
					[dblPredA,matPredA,dblI_RateAndSelfCoinc,dblImat,dblI_diag] = getSeparation(matRateAndSelfCoinc,vecUseTrialStimType,0,dblDiffTheta);
					
					matSelfAndPopCoinc = cat(1,matDiagCoincCort,matMeanCoincCort);
					[dblPredA,matPredA,dblI_SelfAndPopCoinc,dblImat,dblI_diag] = getSeparation(matSelfAndPopCoinc,vecUseTrialStimType,0,dblDiffTheta);
					
					dblI_RateAndPop_bc = dblI_RateAndPopCoinc*dblProdFacRaw-dblSubFac;
					dblI_RateAndSelfCoinc_bc = dblI_RateAndSelfCoinc*dblProdFacRaw-dblSubFac;
					dblI_SelfAndPopCoinc_bc = dblI_SelfAndPopCoinc*dblProdFacRaw-dblSubFac;
					
					%% 90 predictors
					intPredictors = 90;
					dblSubFac =(2*intPredictors)/(intRepetitions*(dblDiffTheta.^2));
					dblProdFacRaw = ((2*intRepetitions-intPredictors-3)/(2*intRepetitions-2));
					
					matComposite = cat(1,matDiagCoincCort,matMeanCoincCort,matMuAct);
					[dblPredA,matPredA,dblI_Composite,dblImat,dblI_diag] = getSeparation(matComposite,vecUseTrialStimType,0,dblDiffTheta);
					
					%is burstiness redundant with diag/mu separate?
					matTriple = cat(1,matBurstiness,matDiagCoincCort,matMuAct);
					[dblPredA,matPredA,dblI_Triple,dblImat,dblI_diag] = getSeparation(matTriple,vecUseTrialStimType,0,dblDiffTheta);
					
					dblI_Composite_bc = dblI_Composite*dblProdFacRaw-dblSubFac;
					dblI_Triple_bc = dblI_Triple*dblProdFacRaw-dblSubFac;
					
					
					%% all coincidence information
					matSelect = tril(true([intCoincN intCoincN]),0);
					intPredictors = sum(matSelect(:));
					matLinCoincCort = nan(intPredictors,intUseTrials);
					for intTrial=1:intUseTrials
						matTemp = matAllCoincCort(:,:,intTrial);
						matLinCoincCort(:,intTrial) = matTemp(matSelect);
					end
					matAll = cat(1,matMuAct,matLinCoincCort);
					[dblPredA,matPredA,dblI_All,dblImat,dblI_diag] = getSeparation(matAll,vecUseTrialStimType,0,dblDiffTheta);
					
					
					%% do plotting
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
						export_fig([strFigDir  'AnalBlock2Subsp_Area' num2str(intWithinArea) '_' strTag '.tif']);
						export_fig([strFigDir  'AnalBlock2Subsp_Area' num2str(intWithinArea) '_' strTag '.pdf']);
					end
				else
					
				end
				
			end
		end
		
		%% save data
		if boolSaveData
			strDataFile = ['D:\Data\Processed\V1_LIFmodel\Simulation_' strSimulation '_AB2_Area' num2str(intWithinArea) '.mat'];
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

					
					
