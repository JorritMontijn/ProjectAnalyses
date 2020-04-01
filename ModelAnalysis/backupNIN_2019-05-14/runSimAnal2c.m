
%% initialize
clearvars;
boolLoad = false;
if ~exist('matModelResp','var')
	clearvars -except bool*;%close all;
	intLoadSim = 24; %% SET WHICH SIM TO LOAD HERE; 41 vs 42; 2 vs. 5
	boolLoad = true;
end

vecRunSims = [21 24 31 34];
for intLoadSim=vecRunSims
	clearvars -except vecRunSims intLoadSim
	boolLoad = true;
	if intLoadSim == 11 && boolLoad
		strSimulation = 'xAreaDistributed_OriFull_2017-06-15';
	elseif intLoadSim == 12 && boolLoad
		strSimulation = 'xAreaDistributed_ContExcOnly_2017-06-08'; %contrast
	elseif intLoadSim == 13 && boolLoad
		strSimulation = 'xAreaDistributed_ContInhOnly_2017-06-12'; %contrast
	elseif intLoadSim == 14 && boolLoad
		strSimulation = 'xAreaDistributed_ContNone_2017-06-12'; %contrast
		
	elseif intLoadSim == 21 && boolLoad
		strSimulation = 'xAreaDistributed_LumFull_2017-06-26'; %luminance full
	elseif intLoadSim == 22 && boolLoad
		strSimulation = ''; %lum exc
	elseif intLoadSim == 23 && boolLoad
		strSimulation = ''; %lum inh
	elseif intLoadSim == 24 && boolLoad
		strSimulation = 'xAreaDistributed_LumNone_2017-06-26'; %lum none
		
		
	elseif intLoadSim == 31 && boolLoad
		strSimulation = 'xAreaDistributed_ContFull_2017-06-26'; %contrast full
	elseif intLoadSim == 32 && boolLoad
		strSimulation = ''; %contrast exc
	elseif intLoadSim == 33 && boolLoad
		strSimulation = ''; %contrast inh
	elseif intLoadSim == 34 && boolLoad
		strSimulation = 'xAreaDistributed_ContNone_2017-06-26'; %contrast none
		
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
	boolSaveData = true;
	dblCutOff = 0.90;
	dblLambdaInfo = 1;
	dblLambdaPred = 0;
	dblNoiseLevel = 0;
	boolAnalInfoSpike = true;
	%boolAnalPred = true;
	%boolAnalPopM = true;
	
	%% calculate which parameter to use
	[d,intUseParam]=max(range(matStimTypeCombos,2));
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
	
	%% add noise to improve numerical stability
	matModelResp = double(matModelResp);
	intNeurons = size(matModelResp,1);
	intTrials = size(matModelResp,2);
	vecNeuronSD = xstd(matModelResp,2);
	matNoise = normrnd(zeros(intNeurons,intTrials),repmat(vecNeuronSD*dblNoiseLevel,[1 intTrials]));
	matModelRespP = matModelResp + matNoise;
	%matModelRespP = zscore(matModelRespP,[],1);
	cellArea = {'V1','V2','V1-V2'};
	figure
	for intWithinArea=[1 2]
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
			sParams.dblLambda = dblLambdaInfo;
			sParams.dblDiffTheta = range(vecTrialOris);
			
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
				sOut = doFisherAnalTrials(matData,vecTrialStimType,sParams);
				
				%get output
				vecA = sOut.matA_LogReg(:);
				vecA(vecA==0) = nan;
				vecA_shuff = sOut.matA_LogReg_shuff(:);
				vecA_shuff(vecA_shuff==0) = nan;
				vecMeanInfoPerSpike(intStimTypeIdx) = mean(sOut.matI_LogReg_bc_CV(:)./vecA);
				vecSDInfoPerSpike(intStimTypeIdx) = std(sOut.matI_LogReg_bc_CV(:)./vecA);
				vecMeanInfoPerSpikeShuff(intStimTypeIdx) = mean(sOut.matI_LogReg_bc_CV_shuff(:)./vecA_shuff);
				vecSDInfoPerSpikeShuff(intStimTypeIdx) = std(sOut.matI_LogReg_bc_CV_shuff(:)./vecA_shuff);
				cellA{intStimTypeIdx} = vecA;
				cellI{intStimTypeIdx} = sOut.matI_LogReg_bc_CV(:);
				cellA_shuff{intStimTypeIdx} = vecA_shuff;
				cellI_shuff{intStimTypeIdx} = sOut.matI_LogReg_bc_CV_shuff(:);
			end
		end
		subplot(2,2,intWithinArea)
		errorbar(vecParamVals,vecMeanInfoPerSpike,vecSDInfoPerSpike,'b');
		hold on
		errorbar(vecParamVals,vecMeanInfoPerSpikeShuff,vecSDInfoPerSpikeShuff,'r');
		hold off
		xlabel(strParam);
		ylabel('Information per spike');
		legend({'Unshuffled',' Shuffled'});
		strTitle = [cellArea{intWithinArea} '  ' strType];
		title(strTitle);
		xlim([-1 101]);
		ylim([0 0.8])
		fixfig;
		
		subplot(2,2,intWithinArea+2)
		scatter(cell2mat(cellA'),cell2mat(cellI'),'b')
		hold on
		scatter(cell2mat(cellA_shuff'),cell2mat(cellI_shuff'),'r')
		hold off
		xlabel('# of spikes in single trial');
		ylabel('Information on stimulus orientation');
		strTitle = [cellArea{intWithinArea} '  ' strType];
		title(strTitle);
		xlim([0 8000]);
		if intWithinArea==1
			ylim([0 3000]);
		else
			ylim([0 600]);
		end
		fixfig;
		
		%% save data
		if boolSaveData
			strDataFile = ['D:\Data\Processed\V1_LIFmodel\Simulation_' strSimulation '_AB2_Area' num2str(intWithinArea) '.mat'];
			if exist(strDataFile,'file')
				sSave = load(strDataFile);
			else
				sSave = struct;
			end
			if boolAnalInfoSpike
				sSave.strParam = strParam;
				sSave.vecMeanInfoPerSpike=vecMeanInfoPerSpike;
				sSave.vecSDInfoPerSpike=vecSDInfoPerSpike;
				sSave.vecMeanInfoPerSpikeShuff=vecMeanInfoPerSpikeShuff;
				sSave.vecSDInfoPerSpikeShuff=vecSDInfoPerSpikeShuff;
				sSave.cellA=cellA;
				sSave.cellI=cellI;
				sSave.cellA_shuff=cellA_shuff;
				sSave.cellI_shuff=cellI_shuff;
				
				sSave.strArea = cellArea{intWithinArea};
				sSave.strType = strType;
				sSave.vecParamVals = vecParamVals;
				sSave.sParams=sParams;
			end
			
			save(strDataFile,'-struct','sSave');
		end
	end
	
	%full screen
	drawnow;
	jFig = get(handle(gcf), 'JavaFrame');
	jFig.setMaximized(true);
	figure(gcf);
	drawnow;
	if boolSaveFigs
		%save figure
		export_fig([strFigDir  'AnalBlock2InfoSpike_Area' num2str(intWithinArea) '_' strTag '.tif']);
		export_fig([strFigDir  'AnalBlock2InfoSpike_Area' num2str(intWithinArea) '_' strTag '.pdf']);
	end
end
