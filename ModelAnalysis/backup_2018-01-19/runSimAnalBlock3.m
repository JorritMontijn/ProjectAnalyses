
%% initialize
clearvars;
boolLoad = false;
if ~exist('matModelResp','var')
	clearvars -except bool*;%close all;
	intLoadSim = 24; %% SET WHICH SIM TO LOAD HERE; 41 vs 42; 2 vs. 5
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
	strFigDir = 'D:\Data\Results\Block3\';
	if boolLoad
		runModelHeader;
	end
	cellIn = strsplit(strSimulation,'_');
	strFiller = cellIn{1};
	strType = cellIn{2};
	strDate = cellIn{3};
	strTag = [strStimType strType strDate];
	
	%% set parameters
	boolSaveFigs = true;
	boolSaveData = false;
	dblCutOff = 0.90;
	dblLambdaInfo = 1;
	dblLambdaPred = 0;
	dblNoiseLevel = 0;
	boolAnalInfoSpike = true;
	%boolAnalPred = true;
	%boolAnalPopM = true;
	
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
	
	%% add noise to improve numerical stability
	matModelResp = double(matModelResp);
	intNeurons = size(matModelResp,1);
	intTrials = size(matModelResp,2);
	vecNeuronSD = xstd(matModelResp,2);
	matNoise = normrnd(zeros(intNeurons,intTrials),repmat(vecNeuronSD*dblNoiseLevel,[1 intTrials]));
	matModelRespP = matModelResp + matNoise;
	%matModelRespP = zscore(matModelRespP,[],1);
	cellArea = {'V1','V2','V1-V2'};
	h1 = figure;
	for intWithinArea=[1 2]
		%% pre-allocate
		vecMeanInfoPerSpike = nan(1,intParamNum);
		vecSDInfoPerSpike = nan(1,intParamNum);
		vecMeanInfoPerSpikeShuff = nan(1,intParamNum);
		vecSDInfoPerSpikeShuff = nan(1,intParamNum);
			
		%% compare different stimulus intensity levels; from blank white noise to structured stimulus
		for intStimTypeIdx=1:intParamNum
			%% start
			dblParamVal = vecParamVals(intStimTypeIdx);
			cellC{intStimTypeIdx} = sprintf('%03d',dblParamVal);
			strC = cellC{intStimTypeIdx};
			fprintf('Starting %s %s (%d/%d) [%s]\n',strParam,strC,intStimTypeIdx,intParamNum,getTime);
			vecUseStimTypes = find(sort(vecParamIdx,'ascend')==intStimTypeIdx);
			
			%check if no spikes
			if (sum(sum(matModelRespP(1:intCellsV1,vecTrialStimType==vecUseStimTypes(1)))) + sum(sum(matModelRespP(1:intCellsV1,vecTrialStimType==vecUseStimTypes(end))))) == 0
				boolAllZero = true;else boolAllZero = false;end
			
			%set parameters
			sParams = struct;
			sParams.vecUseStimTypes = vecUseStimTypes;
			
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
				matData = matModelRespP(vecCells,:);
				[dblPerformance,vecDecodedIndex,matPosteriorProbability] = doCrossValidatedDecodingLR(matData,vecTrialStimType,dblLambdaInfo);
				dblPerformance
				
				matData1 = matModelRespP(1:110,:)';
				vecTrSt1 = vecTrialStimType;
				
				matData2 = [cellMatX{1,1}; cellMatX{1,2}];
				vecTrSt2 = [ones(1,size(cellMatX{1,1},1)) 2*ones(1,size(cellMatX{1,2},1))];
				[dblPredA1,matPredA,dblDprimeSquaredOffDiagonal,matDprimeSquared,matDprimeSquared_diagonal] = getSeparation(matData1,vecTrSt1,1);
				[dblPredA2,matPredA,dblDprimeSquaredOffDiagonal,matDprimeSquared,matDprimeSquared_diagonal] = getSeparation(matData2,vecTrSt2,1);
				return
				%% get Fisher info with full pop
				sParamsFish=struct;
				
				sParamsFish.intSizeY = 0;
				sParamsFish.intResamplings = 2;
				sParamsFish.vecCellArea = vecCellArea;
				sParamsFish.intWithinArea = intWithinArea;
				sParamsFish.vecUseStimTypes = vecUseStimTypes;
				sParamsFish.dblDiffTheta = range(vecTrialOris);
				sParamsFish.dblLambda = dblLambdaInfo;
			
				vecGroupSizes = ([2.^[3:9] intCellsV1]);
				intGroups = numel(vecGroupSizes);
				vecInfoMeans = nan(1,intGroups);
				vecInfoSDs = nan(1,intGroups);
				for intGroupSizeIdx=1:intGroups
					dblGroupSize = vecGroupSizes(intGroupSizeIdx)
					sParamsFish.intSizeX = dblGroupSize;
					[matFisherFull,sOut] = doFisherAnal(matModelRespP,vecTrialStimType,sParamsFish);
					vecI = xmean(matFisherFull,3);
					vecInfoMeans(intGroupSizeIdx) = mean(vecI);
					vecInfoSDs(intGroupSizeIdx) = std(vecI);
				end
				
				subplot(2,2,intWithinArea);
				errorbar(vecGroupSizes,vecInfoMeans,vecInfoSDs);
				xlabel('Number of cells');
				ylabel('Fisher information');
				title(cellArea{intWithinArea});
				
			end
		end
		
		
		%% save data
		if boolSaveData
			strDataFile = ['D:\Data\Processed\V1_LIFmodel\Simulation_' strSimulation '_AB3_Area' num2str(intWithinArea) '.mat'];
			if exist(strDataFile,'file')
				sSave = load(strDataFile);
			else
				sSave = struct;
			end
			if boolAnalInfoSpike
				
			end
			
			save(strDataFile,'-struct','sSave');
		end
	end
	%full screen
	figure(h1);
	drawnow;
	jFig = get(handle(gcf), 'JavaFrame');
	jFig.setMaximized(true);
	figure(gcf);
	drawnow;
	if boolSaveFigs
		%save figure
		export_fig([strFigDir  'AnalBlock3Info_Area' num2str(intWithinArea) '_' strTag '.tif']);
		export_fig([strFigDir  'AnalBlock3Info_Area' num2str(intWithinArea) '_' strTag '.pdf']);
	end
	
end
