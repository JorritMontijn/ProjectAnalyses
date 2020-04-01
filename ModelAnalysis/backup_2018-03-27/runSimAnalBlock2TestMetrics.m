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
vecRunSims = [-1 -2];
for intLoadSim=vecRunSims
	clearvars -except vecRunSims intLoadSim bool*
	boolLoad = true;
	if intLoadSim == -1 && boolLoad
		strSimulation = 'xAreaExperiment_106r001p26-t_2018-02-21';
	elseif intLoadSim == -2 && boolLoad
		strSimulation = 'xAreaExperiment_107l003p143-t_2018-02-21';
		
	elseif intLoadSim == 10 && boolLoad
		strSimulation = 'xAreaDistributed_Ret16Noise0Ori2Drift2Att0_2018-01-30';
	elseif intLoadSim == 11 && boolLoad
		strSimulation = 'xAreaDistributed_Ret16Noise1Ori2Drift2Att0_2018-02-05';
	elseif intLoadSim == 13 && boolLoad
		strSimulation = 'xAreaDistributed_Ret16Noise3Ori2Drift2Att0_2018-02-05';
	elseif intLoadSim == 15 && boolLoad
		strSimulation = 'xAreaDistributed_Ret16Noise5Ori2Drift2Att0_2018-02-05';
		
		
		
	elseif intLoadSim == 20 && boolLoad
		strSimulation = 'xAreaDistributed_Ret32Noise0Ori2Drift2Att0_2018-01-29'; %new connectivity
	elseif intLoadSim == 21 && boolLoad
		strSimulation = 'xAreaDistributed_Ret32Noise1Ori2Drift2Att0_2018-01-29';
	elseif intLoadSim == 23 && boolLoad
		strSimulation = 'xAreaDistributed_Ret32Noise3Ori2Drift2Att0_2018-02-12';
	elseif intLoadSim == 25 && boolLoad
		strSimulation = 'xAreaDistributed_Ret32Noise5Ori2Drift2Att0_2018-02-01';
		
	elseif intLoadSim == 30 && boolLoad
		strSimulation = 'xAreaDistributed_Ret64Noise0Ori2Drift2Att0_2018-02-12';
	elseif intLoadSim == 35 && boolLoad
		strSimulation = 'xAreaDistributed_Ret64Noise5Ori2Drift2Att0_2018-02-08';
		
	elseif intLoadSim == 40 && boolLoad
		strSimulation = 'xAreaDistributed_Ret32Col600Noise0Ori2Drift2Att0_2018-02-12';
	end
	
	%% RUN: #header
	strFigDir = 'D:\Data\Results\Block2\';
	if boolLoad
		runAnalysisHeader;
	end
	return
	cellIn = strsplit(strSimulation,'_');
	strFiller = cellIn{1};
	strType = cellIn{2};
	strDate = cellIn{3};
	strTag = [strType strDate];
	
	%% set parameters
	boolSaveFigs = false;
	boolSaveData = true;
	dblLambdaInfo = 1;
	dblLambdaPred = 0;
	dblNoiseLevel = 0;
	boolAnalInfoSub = true;
	intMaxDimAnal = 30;
	%intWithinArea = 1; %1, within V1; 2, within V2; 3, across V1-V2;
	intIters = 5;
	intStimTypes = numel(vecStimTypeOris);
	intStimPairs = intStimTypes;
	if intStimPairs == 2,intStimPairs=1;end
	
	%% run loop
	cellStrArea = {'V1','V2',''};
	for intWithinArea=[1 2]
		for intCoincN=[30 110]
			for intStimTypeIdx=1%:intStimTypes
				%% start
				vecUseStimTypes = mod([intStimTypeIdx intStimTypeIdx+1],intStimTypes);
				vecUseStimTypes(vecUseStimTypes==0)=intStimTypes;
				%vecUseStimTypes = 1:intStimTypes;
				indUseTrials = ismember(vecTrialStimType,vecUseStimTypes);
				
				%msg
				strArea = cellStrArea{intWithinArea};
				fprintf('Starting area %d (%s) with pop size %d, stims [%s\b] [%s]\n',intWithinArea,strArea,intCoincN,sprintf('%d ',vecUseStimTypes),getTime);
				
				%check if data exists
				strDataFile = ['D:\Data\Processed\V1_LIFmodel\Simulation_' strSimulation '_AB2_Area' num2str(intWithinArea) '.mat'];
				%if exist(strDataFile,'file')
				%	sSave = load(strDataFile);
				%else
					sSave = struct;
				%end
				if isfield(sSave,'vecNumObs') && numel(sSave.vecNumObs) >= intCoincN && ~isempty(sSave.vecNumObs{intCoincN})
					fprintf('Data for area %d (%s) with pop size %d for data set [%s] already exists; skipping.. [%s]\n',...
						intWithinArea,strArea,intCoincN,strSimulation,getTime);
				%	continue;
				end
				
				
				%% pre-allocate variables dim dep
				hFigDD = figure;
				%full screen
				drawnow;
				jFig = get(handle(gcf), 'JavaFrame');
				jFig.setMaximized(true);
				figure(gcf);
				drawnow;
				
				%% select one pair of orientations
				vecUseTrialStimType = vecTrialStimType(indUseTrials);
				vecUseStimStartSecs = vecStimStartSecs(indUseTrials);
				vecUseStimStopSecs = vecStimStopSecs(indUseTrials);
				dblDiffTheta = range(vecTrialOris(vecUseStimTypes));
				%%
				%vecTrialStimType = vecUseTrialStimType
				%vecStimStartSecs = vecUseStimStartSecs
				%vecStimStopSecs = vecUseStimStopSecs
				
				%% pre-alloc
				intUseTrials = numel(vecUseTrialStimType);
				vecNumObs = round(2.^[6:14]);
				vecNumObs(vecNumObs > intUseTrials) = [];
				intNumObs = numel(vecNumObs);
				matI_Rate = nan(intIters,intNumObs);
				matI_PopCoinc = nan(intIters,intNumObs);
				matI_SelfCoinc = nan(intIters,intNumObs);
				matI_Burst = nan(intIters,intNumObs);
				matI_Sync = nan(intIters,intNumObs);
				matI_Simult = nan(intIters,intNumObs);
				cellNeurons = cell(intIters);
				%%
				for intResampling=1:intIters
					fprintf('Running resampling %d/%d for pop size %d [%s]\n',intResampling,intIters,intCoincN,getTime);
					if intWithinArea==1
						vecNeurons = randperm(intCellsV1,intCoincN);
					elseif intWithinArea==2
						vecNeurons = randperm(intCellsV2,intCoincN)+intCellsV1;
					end
					cellSpikeTimes = cellSpikeTimesCortex(vecNeurons);
					[sMetricsOut, vecNumObs] = testCoincTrialDep(cellSpikeTimes,vecUseTrialStimType,vecUseStimStartSecs,vecUseStimStopSecs,vecNumObs,dblDiffTheta);
					matI_Rate(intResampling,:) = sMetricsOut.vecI_Rate_bc;
					matI_PopCoinc(intResampling,:) = sMetricsOut.vecI_PopCoinc_bc;
					matI_SelfCoinc(intResampling,:) = sMetricsOut.vecI_SelfCoinc_bc;
					matI_Burst(intResampling,:) = sMetricsOut.vecI_Burst_bc;
					matI_Sync(intResampling,:) = sMetricsOut.vecI_Sync_bc;
					matI_Simult(intResampling,:) = sMetricsOut.vecI_Simult_bc;
					cellNeurons{intResampling} = vecNeurons;
				end
				sCoincMetrics.matI_Rate = matI_Rate;
				sCoincMetrics.matI_PopCoinc = matI_PopCoinc;
				sCoincMetrics.matI_SelfCoinc = matI_SelfCoinc;
				sCoincMetrics.matI_Burst = matI_Burst;
				sCoincMetrics.matI_Sync = matI_Sync;
				sCoincMetrics.matI_Simult = matI_Simult;
				sCoincMetrics.cellNeurons = cellNeurons;
				
				%%
				hold all
				errorbar(vecNumObs,xmean(matI_Rate,1),xstd(matI_Rate,1)/sqrt(intIters))
				errorbar(vecNumObs,xmean(matI_PopCoinc,1),xstd(matI_PopCoinc,1)/sqrt(intIters))
				errorbar(vecNumObs,xmean(matI_SelfCoinc,1),xstd(matI_SelfCoinc,1)/sqrt(intIters))
				errorbar(vecNumObs,xmean(matI_Burst,1),xstd(matI_Burst,1)/sqrt(intIters))
				errorbar(vecNumObs,xmean(matI_Sync,1),xstd(matI_Sync,1)/sqrt(intIters))
				errorbar(vecNumObs,xmean(matI_Simult,1),xstd(matI_Simult,1)/sqrt(intIters))
				fixfig
				xlabel('Number of trials');
				ylabel('Bias-corr Fisher I');
				legend({'Rate','Pop coinc','Self coinc','Burst','Sync','Simult'})
				title(sprintf('%s %s; %d predictors, %d resamplings; Burst=SelfC/Rate; Sync=PopC/Rate, Simult=PopC/SelfC',strArea,strType,intCoincN,intIters));
				
				%save figure
				if boolSaveFigs
					figure(hFigDD);drawnow;
					export_fig([strFigDir  'CoincidenceMetricsTestPop' num2str(intCoincN) '_Area' num2str(intWithinArea) '_' strTag '.tif']);
					export_fig([strFigDir  'CoincidenceMetricsTestPop' num2str(intCoincN) '_Area' num2str(intWithinArea) '_' strTag '.pdf']);
				end
				
				%% save data
				if boolSaveData
					%check if data exists
					strDataFile = ['D:\Data\Processed\V1_LIFmodel\Simulation_' strSimulation '_AB2_Area' num2str(intWithinArea) '.mat'];
					if exist(strDataFile,'file')
						sSave = load(strDataFile);
					else
						sSave = struct;
					end
					sSave.vecNumObs{intCoincN,intStimTypeIdx} = vecNumObs;
					sSave.sCoincMetrics{intCoincN,intStimTypeIdx} = sCoincMetrics;
					
					save(strDataFile,'-struct','sSave');
				end
			end
		end
	end
end


