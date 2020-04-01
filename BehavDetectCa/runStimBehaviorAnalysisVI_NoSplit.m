%function vecPredictorAssemblyZscores = doPlotActivityDistributions(cellContrastDataAD,sMetaDataAD,sPC,sParams)

%if ~exist('sParams','var'),sParams = struct;end
%if isfield(sParams,'intStartSes'),intStartSes=sParams.intStartSes;else intStartSes=1;end
%if isfield(sParams,'intStopSes'),intStopSes=sParams.intStopSes;else intStopSes=length(cellContrastDataAD);end
for intMouse=1:8
	close all
	clearvars -except intMouse
	%get block data
	if ~exist('cellMultiSes','var')
		if intMouse == 1
			strSes = '20140207';
		elseif intMouse == 2
			strSes = '20140314';
		elseif intMouse == 3
			strSes = '20140425';
		elseif intMouse == -1
			%return; %exclude?; bad behavior, weird signals
			strSes = '20140430';
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
	strSes = ['NS_Supp2_' strSes];
	
	for intPopulation = vecBlockTypes
		%get neuronal tuning
		if intMouse==8
			%remove last trial
			vecRem = true(size(cellMultiSes{1}.structStim.Orientation));
			vecRem((end-47):end) = false;
			cellMultiSes{1}.structStim = remel(cellMultiSes{1}.structStim,vecRem);
			%recalc dfof
			%cellMultiSes{1} = doRecalcdFoF(cellMultiSes{1},3);
		end
		
		%[indTuned,vecNeuronPrefStim,vecOSI,cellAggregate] = getTunedStimDetectionNeurons(cellSes(vecBlock==intPopulation));
		[indTuned,vecNeuronPrefStim,vecOSI,cellAggregate] = getTunedStimDetectionNeurons(cellMultiSes{intPopulation});
		structStim = cellMultiSes{intPopulation}.structStim;
		
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
		%cellMultiSes{intPopulation}.structStim = structStim;
		
		%% eye movements/blinks
		%get data
		vecAreaZ = cellMultiSes{intPopulation}.sEyeTracking.vecArea_Z;
		vecRoundZ = cellMultiSes{intPopulation}.sEyeTracking.vecRoundness_Z;
		vecPosX_Z = cellMultiSes{intPopulation}.sEyeTracking.vecPosX_Z;
		vecPosY_Z = cellMultiSes{intPopulation}.sEyeTracking.vecPosY_Z;
		vecPupilLuminance_Z = cellMultiSes{intPopulation}.sEyeTracking.vecPupilLuminance_Z;
		
		%remove unrealistic epochs
		vecDiffA = abs(diff(vecAreaZ)) > 1;
		vecRem = [false vecDiffA] | [vecDiffA false];
		
		vecRoundZ(vecRem) = nan;
		vecPosX_Z(vecRem) = nan;
		vecPosY_Z(vecRem) = nan;
		
		%get blinks and saccades
		vecBlinks = vecPupilLuminance_Z < -5; %-4
		vecSaccades = abs(vecRoundZ) < 3 & (abs(vecPosY_Z) > 5 | abs(vecPosX_Z) > 5);

		%get trials
		[dummy,vecTrialsBlinks] = getStimAtFrame(cellMultiSes{intPopulation},find(vecBlinks));
		[dummy,vecTrialsSaccades] = getStimAtFrame(cellMultiSes{intPopulation},find(vecSaccades));
		
		vecRemTrials = unique([vecTrialsBlinks vecTrialsSaccades]);
		vecKeepTrials = 1:length(cellMultiSes{intPopulation}.structStim.Orientation);
		vecKeepTrials(vecRemTrials) = [];
		fprintf('Session %s pop %d; %d/%d trials removed [%.1f%%]\n',cellMultiSes{intPopulation}.session,intPopulation,length(vecRemTrials),length(cellMultiSes{intPopulation}.structStim.Orientation),(length(vecRemTrials)/length(cellMultiSes{intPopulation}.structStim.Orientation))*100);
		
		%remove trials
		structStim = structfun(@(x) ( x(vecKeepTrials) ), cellMultiSes{intPopulation}.structStim, 'UniformOutput', false);
		cellMultiSes{intPopulation}.structStim = structStim;
		
		%% get pre-formatted data variables
		%get neuronal responses per trial
		[matTrialResponse,cellSelectContrasts] = getTrialResponseData(cellMultiSes{intPopulation},structStim);
		intContrasts = length(cellSelectContrasts);
		intTrials = length(structStim.Orientation);
		
		% get stimulus data
		[vecStimOris,dummy]=find(bsxfun(@eq,structStim.Orientation,vecOrientations'));
		vecStimOris = vecStimOris';
		
		[vecStimContrasts,dummy]=find(bsxfun(@eq,structStim.Contrast,unique(structStim.Contrast)'));
		vecStimContrasts = vecStimContrasts';
		
		indStimResp = structStim.vecTrialResponse==1;
		
		%% calc heterogeneity & mean pupil size
		%normalize
		matRespNormPerNeuron = zscore(matTrialResponse,[],2);
		
		matSelect = tril(true(intNeurons,intNeurons),-1);
		vecHeterogeneity = nan(intTrials,1);
		vecPupilSize = nan(intTrials,1);
		for intTrial=1:intTrials
			% STANDARD
			%perform calculation for all neurons
			vecActivity = matRespNormPerNeuron(:,intTrial);
			matDistAll = abs(bsxfun(@minus,vecActivity,vecActivity'));
			
			%save data
			vecHeterogeneity(intTrial) = mean(matDistAll(matSelect));
			
			%pupil size
			intStartFrame = structStim.FrameOn(intTrial)-round(0.5*cellMultiSes{intPopulation}.samplingFreq);
			intStopFrame = structStim.FrameOn(intTrial)-1;
			vecPupilSize(intTrial) = nanmean(vecAreaZ(intStartFrame:intStopFrame));
		end
		
		%% eye-tracking stuff
		%pref pop act miss still/move + hit still/move
		%get other data
		vecOriTrials = (cellMultiSes{intPopulation}.structStim.Orientation/45)+1;
		indSelectResp = cellMultiSes{intPopulation}.structStim.vecTrialResponse;
		
		%pre-allocate
		vecContrasts = unique(cellMultiSes{1}.structStim.Contrast)*100;
		vecContrasts(1) = 0.2;
		intContrasts = length(vecContrasts);
		intContrastCounter = 0;
		for intContrastIndex=1:intContrasts
			%get contrast
			dblContrast = vecContrasts(intContrastIndex);
			intContrastCounter = intContrastCounter + 1;
			
			%overall
			intTrialCounter = 0;
			vecTrials = find(cellSelectContrasts{intContrastCounter} & indSelectResp);
			vecActHit = nan(1,length(vecTrials));
			for intTrial=vecTrials
				intTrialCounter = intTrialCounter + 1;
				indPrefPop = vecNeuronPrefStim == vecOriTrials(intTrial);
				vecActHit(intTrialCounter) = mean(matTrialResponse(indPrefPop,intTrial));
			end
			
			intTrialCounter = 0;
			vecTrials = find(cellSelectContrasts{intContrastCounter} & ~indSelectResp);
			vecActMiss = nan(1,length(vecTrials));
			for intTrial=vecTrials
				intTrialCounter = intTrialCounter + 1;
				indPrefPop = vecNeuronPrefStim == vecOriTrials(intTrial);
				vecActMiss(intTrialCounter) = mean(matTrialResponse(indPrefPop,intTrial));
			end
			
			%heterogen
			vecHetHit = vecHeterogeneity(cellSelectContrasts{intContrastCounter} & indSelectResp)';
			vecHetMiss = vecHeterogeneity(cellSelectContrasts{intContrastCounter} & ~indSelectResp)';
			
			%pupil size
			vecPupHit = vecPupilSize(cellSelectContrasts{intContrastCounter} & indSelectResp)';
			vecPupMiss = vecPupilSize(cellSelectContrasts{intContrastCounter} & ~indSelectResp)';
			
			%put in output
			cellSaveMatContAct{intPopulation}(intContrastCounter,1) = mean(vecActHit); %overall hit
			cellSaveMatContAct{intPopulation}(intContrastCounter,2) = mean(vecActMiss); %overall miss
			cellSaveMatContAct{intPopulation}(intContrastCounter,3) = mean(vecHetHit); %overall het hit
			cellSaveMatContAct{intPopulation}(intContrastCounter,4) = mean(vecHetMiss); %overall het miss
			cellSaveMatContAct{intPopulation}(intContrastCounter,5) = nanmean(vecPupHit); %overall pup hit
			cellSaveMatContAct{intPopulation}(intContrastCounter,6) = nanmean(vecPupMiss); %overall pup miss
			cellSaveMatContAct{intPopulation}(intContrastCounter,7) = getCohensD(vecActHit,vecActMiss); %overall MES dFoF
			cellSaveMatContAct{intPopulation}(intContrastCounter,8) = getCohensD(vecHetHit,vecHetMiss); %overall MES Het
		end
		
		%% remove different quintiles
		vecHetAll = nan(1,intTrials);
		matHetQ = nan(5,intTrials);
		for intTrial=1:intTrials
			%% STANDARD
			%perform calculation for all neurons
			%perform calculation for all neurons
			vecActivity = matRespNormPerNeuron(:,intTrial);
			matDistAll = abs(bsxfun(@minus,vecActivity,vecActivity'));
			
			%save data as vectors in matrix
			matSelectAll = tril(true(size(matDistAll)),-1);
			vecHetAll(intTrial) = mean(matDistAll(matSelectAll));
			
			%% REMOVE quintiles
			%do calculation
			dblSize = 0.2;
			intRemNeurons = round(intNeurons*dblSize);
			vecNeurons = 1:intNeurons;
			for intQuintile = 1:5
				intStartNeuron = floor(intNeurons*((intQuintile-1)*dblSize))+1;
				vecNeuronsQ = vecNeurons;
				vecNeuronsQ(intStartNeuron:(intStartNeuron+intRemNeurons-1)) = [];
				
				vecActivityQ = vecActivity(vecNeuronsQ);
				matDist = abs(bsxfun(@minus,vecActivityQ,vecActivityQ'));
				
				%save data as vectors in matrix
				matSelect = tril(true(size(matDist)),-1);
				matHetQ(intQuintile,intTrial) = mean(matDist(matSelect));
			end
		end
		
		%save data
		cellSaveHetQs{intPopulation} = matHetQ; 
		cellSaveHetAll{intPopulation} = vecHetAll;
		cellSavePupAll{intPopulation} = vecPupilSize;
		cellSaveRespTimeTrials{intPopulation} = cellMultiSes{intPopulation}.structStim.vecTrialRespSecs-cellMultiSes{intPopulation}.structStim.SecsOn; 
		cellSaveContrastTrials{intPopulation} = cellMultiSes{intPopulation}.structStim.Contrast; 
		
		%% brain state analysis
		%take blocks of time of good/bad performance, then look at neural
		%correlates
		% => canned
		
		%% investigate multiplicative gain for 33% HCN from miss to hit
		% => canned
		
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
		save(['D:\Data\Results\stimdetection\' strFile],'cellSave*','-v7.3');
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