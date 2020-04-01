%function vecPredictorAssemblyZscores = doPlotActivityDistributions(cellContrastDataAD,sMetaDataAD,sPC,sParams)

%if ~exist('sParams','var'),sParams = struct;end
%if isfield(sParams,'intStartSes'),intStartSes=sParams.intStartSes;else intStartSes=1;end
%if isfield(sParams,'intStopSes'),intStopSes=sParams.intStopSes;else intStopSes=length(cellContrastDataAD);end
clear all;
%close all;
for intMouse=5
	%close all
	clearvars -except intMouse
	intUseNeuropilSubtraction = -1; %[-1=no prepro,0=none,1=pre,2=post]
	boolEvents = false;
	boolExcludeLocomotor = false;
	strAnalyzeType = 'neuron';
	
	%% load data
	loadSes;
	if boolEvents
		strSes = ['Ev' strSes];
	end
	
	%% run populations
	for intPopulation = vecBlockTypes
		%% GET RESPONSES AFTER COMBINING BLOCKS OF SAME POPULATION INTO ONE
		%msg
		fprintf('Now processing %s [pop %d], analyzing %ss [%s]\n',strSes,intPopulation,strAnalyzeType,getTime);
		
		%recalc dfof
		if intUseNeuropilSubtraction == 2
			dblNeuropilSubtractionFactor = [];
			cellMultiSes{intPopulation} = doRecalcdFoF(cellMultiSes{intPopulation},5,[],strAnalyzeType,[],[],dblNeuropilSubtractionFactor);
		elseif intUseNeuropilSubtraction == 0
			cellMultiSes{intPopulation} = doRecalcdFoF(cellMultiSes{intPopulation},3,[],strAnalyzeType);
		else
			%do nothing
		end
		
		%replace neurons with astrocytes
		sObject = cellMultiSes{intPopulation}.(strAnalyzeType);
		boolOnlyPresence = true;
		
		%% run header
		runAstroHeader;
		intTotFrameNr = numel(sObject(1).dFoF);
		dblSampFreq = cellMultiSes{intPopulation}.samplingFreq;
		
		%% get pre-formatted data variables
		%get neuronal responses per trial
		[matTrialResponse,cellSelectContrasts] = getTrialResponseData(cellMultiSes{intPopulation},structStim);
		intContrasts = length(cellSelectContrasts);
		
		% get stimulus data
		[vecStimOris,dummy]=find(bsxfun(@eq,structStim.Orientation,vecOrientations'));
		vecStimOris = vecStimOris';
		
		[vecStimContrasts,dummy]=find(bsxfun(@eq,structStim.Contrast,unique(structStim.Contrast)'));
		vecStimContrasts = vecStimContrasts';
		
		indStimResp = structStim.vecTrialResponse==1;
		
		%% run analysis
		intObjects = numel(sObject);
		for intObject=1:intObjects
			vecdFoF = sObject(intObject).dFoF-median(sObject(intObject).dFoF);
			if boolEvents
				indAbove0 = vecdFoF > 0;
				vecStarts = indAbove0 - [false indAbove0(1:(end-1))];
				intEpochStarts = sum(vecStarts==1);
				vecEpochStarts = find(vecStarts==1);
				vecEpochStops = find(vecStarts==-1);
				if numel(vecEpochStops) < numel(vecEpochStarts),vecEpochStops(end+1) = numel(indAbove0);end
				vecEpochDurs = vecEpochStops - vecEpochStarts;
				%remove all durations under 500ms
				indDelete = vecEpochDurs < (dblSampFreq*0.5);
				vecEpochStarts(indDelete) = [];
				vecEpochStops(indDelete) = [];
				vecEpochDurs(indDelete) = [];
				vecAUC = nan(size(vecEpochDurs));
				vecR2 = nan(size(vecEpochDurs));
				
				%% fit
				for intEpoch=1:numel(vecEpochStarts)
					dblStart = vecEpochStarts(intEpoch)-1;
					dblStop = vecEpochStops(intEpoch);
					vecY = vecdFoF(dblStart:dblStop);
					%vecY = vecY - min(vecY);
					vecX = 1:numel(vecY);
					vecP0 = [3 1 1 1]; %alpha, beta, y-scale, x-scale
					vecFitP = curvefitfun(@getGamma,vecP0,1:numel(vecY),vecY);
					vecFitY = getGamma(vecFitP,vecX);
					vecAUC(intEpoch) = sum(vecY);
					vecR2(intEpoch) =  getR2(vecY,vecFitY);
					
					%plot(vecX,vecY);
					%hold on
					%plot(vecX,vecFitY,'--');
					%hold off
					%title(sprintf('Epoch %d; R^2=%.3f; AUC=%.3f,dur=%d',intEpoch,vecR2(intEpoch),vecAUC(intEpoch),vecEpochDurs(intEpoch)));
					%pause;
				end
				
				%% remove R^2 <0.5
				indDelete = vecR2 < 0.7;
				vecEpochStarts(indDelete) = [];
				vecEpochStops(indDelete) = [];
				vecEpochDurs(indDelete) = [];
				vecAUC(indDelete) = [];
				vecR2(indDelete) = [];
				%add events to vector
				vecEvents = zeros(1,intTotFrameNr);
				for intEpoch=1:numel(vecEpochDurs)
					vecEvents(vecEpochStarts(intEpoch):vecEpochStops(intEpoch)) = 1;
				end
				%get data
				vecData = vecEvents;
			else
				%get data
				vecData = vecdFoF;
			end
			
			%% do analysis
			%prep
			cellContsMiss = cell(1,intContrasts);
			cellContsHit = cell(1,intContrasts);
			vecPEP = round(-3*dblSampFreq):round(5*dblSampFreq);
			%go through contrasts
			%clf;
			for intContrast=1:intContrasts
				
				%select trials
				indContrastTrials = cellSelectContrasts{intContrast};
				vecHits = find(indContrastTrials & indStimResp);
				vecMisses = find(indContrastTrials & ~indStimResp);
				
				%prep data
				matHitPEP = nan(numel(vecHits),numel(vecPEP));
				matMissPEP = nan(numel(vecMisses),numel(vecPEP));
				
				%go through hits & misses
				for intHit=1:numel(vecHits)
					intTrial = vecHits(intHit);
					vecSelect = vecPEP+structStim.FrameOn(intTrial);
					indKeep = ~(vecSelect<1 | vecSelect>intTotFrameNr);
					matHitPEP(intHit,indKeep) = vecData(vecSelect(indKeep));
				end
				for intMiss=1:numel(vecMisses)
					intTrial = vecMisses(intMiss);
					vecSelect = vecPEP+structStim.FrameOn(intTrial);
					indKeep = ~(vecSelect<1 | vecSelect>intTotFrameNr);
					matMissPEP(intMiss,indKeep) = vecData(vecSelect(indKeep));
				end
				%save data
				cellSaveRespPEP{intPopulation}{intObject,intContrast,1} = matMissPEP;
				cellSaveRespPEP{intPopulation}{intObject,intContrast,2} = matHitPEP;
				%plot
				%{
			subplot(2,3,intContrast)
			errorfill(vecPEP/dblSampFreq,nanmean(matHitPEP,1),nanstd(matHitPEP,[],1)/sqrt(numel(vecHits)),[1 0 0],[1 0.7 0.7]);
			hold on
			errorfill(vecPEP/dblSampFreq,nanmean(matMissPEP,1),nanstd(matMissPEP,[],1)/sqrt(numel(vecMisses)));
				%}
			end
			%title(sprintf('%s %d/%d [%s (pop %d)]',strAnalyzeType,intObject,intObjects,strSes,intPopulation));
			%drawnow;
			
			%save
			
		end
		%%
		figure;
		for intContrast=1:intContrasts
			cellMiss = cellSaveRespPEP{intPopulation}(:,intContrast,1);
			matMeanMiss = cell2mat(cellfun(@nanmean,cellMiss,cellfill(1,size(cellMiss)),'UniformOutput',false));
			cellSaveMeanMiss{intPopulation}(:,:,intContrast) = matMeanMiss;
			
			cellHit = cellSaveRespPEP{intPopulation}(:,intContrast,2);
			matMeanHit = cell2mat(cellfun(@nanmean,cellHit,cellfill(1,size(cellHit)),'UniformOutput',false));
			cellSaveMeanHit{intPopulation}(:,:,intContrast) = matMeanHit;
			
			%plot
			subplot(2,3,intContrast)
			errorfill(vecPEP/dblSampFreq,nanmean(matMeanHit,1),nanstd(matMeanHit,[],1)/sqrt(numel(vecHits)),[0 1 0],[0.7 1 0.7]);
			hold on
			errorfill(vecPEP/dblSampFreq,nanmean(matMeanMiss,1),nanstd(matMeanMiss,[],1)/sqrt(numel(vecMisses)),[1 0 0],[1 0.7 0.7]);
		end
		drawnow;
	end
	%% save data structures
	if sParams.boolSaveData && ~isempty(whos('cellSave*'))
		vecClock = fix(clock);
		strDate = num2str(vecClock(1:3));
		strFile = ['data_aggregate_' strAnalyzeType strSes '_' strrep(strDate,'    ','_')];
		save(['D:\Data\ResultsAstroAnalysis\' strFile],'cellSave*','-v7.3');
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