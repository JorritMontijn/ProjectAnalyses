%function doStimDetectDecoding
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%get block data
clear all;
for intMouse=1:8
	clearvars -except intMouse
	if intMouse == 1
		strSes = '20140207';
	elseif intMouse == 2
		strSes = '20140314';
	elseif intMouse == 3
		strSes = '20140425';
	elseif intMouse == -1 %exclude
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
	load(['D:\Data\Results\stimdetection\dataRawPre_aggregate' strSes '.mat']);
	
	%parameters
	dblBaselineSecs = 0.1;
	intIterMax = 250;
	cellSaveDecoding = struct;
	
	%% GET RESPONSES AFTER COMBINING BLOCKS OF SAME POPULATION INTO ONE
	for intPopulation = vecBlockTypes
		
		%take opposite directions as the same
		cellMultiSes{intPopulation}.structStim.Orientation = mod(cellMultiSes{intPopulation}.structStim.Orientation,180);
		vecOrientations = unique(cellMultiSes{intPopulation}.structStim.Orientation);
		intNeurons = numel(cellMultiSes{intPopulation}.neuron);
		
		%get neuronal responses per trial
		[matTrialResponse,cellSelectContrasts] = getTrialResponseData(cellMultiSes{intPopulation},cellMultiSes{intPopulation}.structStim);
		intContrasts = length(cellSelectContrasts);
		
		%get orientation-based trial selection vectors
		sTypesOri = getStimulusTypes(cellMultiSes{intPopulation},{'Orientation'});
		cellSelectOri = getSelectionVectors(cellMultiSes{intPopulation}.structStim,sTypesOri);
		
		%make normalized activity matrix per contrast
		matTrialNormResponse = zeros(size(matTrialResponse));
		for intContrast=1:intContrasts
			matTrialNormResponse(:,cellSelectContrasts{intContrast}) = imnorm(matTrialResponse(:,cellSelectContrasts{intContrast}));
		end
		
		%calculate preferred stimulus orientation for all neurons
		vecOriPrefContrasts = 4:6; % 8% - 100%
		matPref = nan(length(vecOriPrefContrasts),intNeurons);
		intCounter = 0;
		for intContrast = vecOriPrefContrasts
			intCounter = intCounter + 1;
			matResp = matTrialNormResponse(:,cellSelectContrasts{intContrast});
			structStimC{intContrast} = cellMultiSes{intPopulation}.structStim;
			cellFields = fieldnames(cellMultiSes{intPopulation}.structStim);
			for intField=1:length(cellFields)
				strField = cellFields{intField};
				structStimC{intContrast}.(strField) = structStimC{intContrast}.(strField)(cellSelectContrasts{intContrast});
			end
			cellSelect = getSelectionVectors(structStimC{intContrast},sTypesOri);
			sTuning{intContrast} = calcTuningRespMat(matResp,cellSelect,vecOrientations);
			matPref(intCounter,:) = sTuning{intContrast}.vecPrefIndex;
		end
		vecNeuronPrefStim = nan(1,intNeurons);
		for intOri=1:length(vecOrientations);
			vecNeuronPrefStim(sum(matPref == intOri,1) > 1) = intOri;
		end
		
		%remove non-tuned neurons
		cellMultiSes{intPopulation}.neuron(isnan(vecNeuronPrefStim)) = [];
		vecNeuronPrefStim(isnan(vecNeuronPrefStim)) = [];
		intNeurons = numel(cellMultiSes{intPopulation}.neuron);
		intTrials = length(cellMultiSes{intPopulation}.structStim.Orientation);
		
		%get updated response matrices
		%get neuronal responses per trial
		[matTrialResponse,cellSelectContrasts] = getTrialResponseData(cellMultiSes{intPopulation},cellMultiSes{intPopulation}.structStim);
		
		%make normalized activity matrix per contrast
		matTrialNormResponse = zeros(size(matTrialResponse));
		for intContrast=1:intContrasts
			matTrialNormResponse(:,cellSelectContrasts{intContrast}) = imnorm(matTrialResponse(:,cellSelectContrasts{intContrast}));
		end
		
		
		%% calculate hit-correlated activity enhancement
		%pre-allocate
		matHCAE = nan(intNeurons,length(cellSelectContrasts));
		indHitTrials = logical(cellMultiSes{intPopulation}.structStim.vecTrialResponse);
		for intContrast=1:intContrasts
			indSelectContrast = cellSelectContrasts{intContrast};
			%select pref stim trials per neuron
			for intNeuron=1:intNeurons
				intPrefStim = vecNeuronPrefStim(intNeuron);
				indPrefStimTrials = cellSelectOri{intPrefStim};
				indPrefHitTrials = indPrefStimTrials & indHitTrials & indSelectContrast;
				indPrefMissTrials = indPrefStimTrials & ~indHitTrials & indSelectContrast;
				
				%get normalized response to calculate hit correlated activity enhancement (range [-1 1])
				vecHitNormResp = matTrialNormResponse(intNeuron,indPrefHitTrials);
				vecMissNormResp = matTrialNormResponse(intNeuron,indPrefMissTrials);
				dblMeanHitResp = mean(vecHitNormResp);if dblMeanHitResp < 0,dblMeanHitResp = 0;end
				dblMeanMissResp = mean(vecMissNormResp);if dblMeanMissResp < 0,dblMeanMissResp = 0;end
				matHCAE(intNeuron,intContrast) = (dblMeanHitResp-dblMeanMissResp)/(dblMeanHitResp+dblMeanMissResp);
			end
		end
		
		%set nans to 0
		matHCAE(isnan(matHCAE)) = 0;
		vecStimDetectActInc = mean(matHCAE(:,2:5),2);
		
		%define correlated/uncorrelated neurons
		dblFrac = 1/3;
		intNumNeurons = round(intNeurons * dblFrac);
		[vecHCAE,vecNeuronsHCAE] = findmax(vecStimDetectActInc,intNumNeurons);
		[vecHUAE,vecNeuronsHAAE] = findmin(vecStimDetectActInc,intNumNeurons);
		indNeuronsHCAE = false(size(vecStimDetectActInc));
		indNeuronsHCAE(vecNeuronsHCAE) = true;
		indNeuronsHAAE = false(size(vecStimDetectActInc));
		indNeuronsHAAE(vecNeuronsHAAE) = true;
		
		indSelectHitCorrelatedNeurons = false(size(vecStimDetectActInc));
		indSelectHitCorrelatedNeurons(vecNeuronsHCAE) = true;
		indSelectHitUncorrelatedNeurons = false(size(vecStimDetectActInc));
		indSelectHitUncorrelatedNeurons(~indNeuronsHCAE & ~indNeuronsHAAE) = true;
		
		%rename to ses
		ses = cellMultiSes{intPopulation};
		
		% get stim list
		sTypes = getStimulusTypes(ses);
		cellSelect = getSelectionVectors(ses.structStim,sTypes);
		intContrastNum = sTypes.vecNumTypes(2);
		intOriNum = sTypes.vecNumTypes(1);
		sTypesC = getStimulusTypes(ses,{'Contrast'});
		cellSelectC = getSelectionVectors(ses.structStim,sTypesC);
		sTypesO = getStimulusTypes(ses,{'Orientation'});
		cellSelectO = getSelectionVectors(ses.structStim,sTypesO);
		
		% get hit/miss trials
		vecHits = ses.structStim.vecTrialResponse;
		
		%set base parameters
		verbose = 1;
		sParams = struct;
		sParams.sTypes = sTypesO;
		sParams.cellSelect = cellSelectO;
		sParams.verbose = 0;
		sParams.intPreBaselineRemoval = round(ses.samplingFreq*dblBaselineSecs);
		indLikelihoodTemplate = false(1,length(cellSelectO));
		
		%pre-allocate output
		matHD_Orientation = nan(intTrials,2);
		matHD_Contrast = nan(intTrials,2);
		matHD_ReactionTime = nan(intTrials,2);
		matHD_DecodingAccuracy = nan(intTrials,2);
		matHD_HeterogeneityPre = nan(intTrials,2);
		matHD_ActivityPre = nan(intTrials,2);
		matHD_HeterogeneityDuring = nan(intTrials,2);
		matHD_ActivityDuring = nan(intTrials,2);
		
		for intNeuronPopulationType = 1:2
			if intNeuronPopulationType == 1
				%hit-correlated neurons
				vecIncludeCells = find(indSelectHitCorrelatedNeurons);
			elseif intNeuronPopulationType == 2
				%non-hit-correlated neurons
				vecIncludeCells = find(indSelectHitUncorrelatedNeurons);
			end
			
			%pre-alloc output
			vecCorrect = false(1,intTrials);
			cellVecPost = cell(1,intTrials);
			vecDecodedStimType = nan(1,intTrials);
			vecStimType = nan(1,intTrials);
		
			%set neuronal populations
			sParams.vecIncludeCells = vecIncludeCells;
			
			
			%decode
			
			%split data set per contrast; perform decoding of 4
			%orientations for each contrast level
			%=> output is for each trial; its orientation, contrast, RT,
			%decoding accuracy (correct/incorrect), pre-stim heterogeneity
			%and during stim heterogeneity for both HCN and HUN
			for intContrast=1:intContrastNum
				%select trials
				indContrastSelect = cellSelectC{intContrast};
				indSelectTheseTrials = indContrastSelect;
				vecOutC = nan(1,length(indSelectTheseTrials));
				
				%leave-one-repetition-out cross validation
				intNumRepetitions = sum(indContrastSelect) / intOriNum;
				for intRepetition = 1:intNumRepetitions
					%select subset of trials to build likelihood
					vecToBeDecodedTrials = nan(1,intOriNum);
					for intOriType=1:intOriNum
						vecSubSelect = find(cellSelectO{intOriType}&indSelectTheseTrials,intRepetition);
						vecToBeDecodedTrials(intOriType) = vecSubSelect(end);
					end
					
					%make likelihood vector
					vecTheseLikelihoodTrials = indContrastSelect;
					vecTheseLikelihoodTrials(vecToBeDecodedTrials) = false;
					sParams.vecLikelihoodTrials = vecTheseLikelihoodTrials;
					sParams.vecDecodeTrials = vecToBeDecodedTrials;
					
					%do decoding
					sOutSD = doStimDecoding(ses,sParams);
					
					%put in output
					cellVecPost(vecToBeDecodedTrials) = sOutSD.cellVecPost(vecToBeDecodedTrials);
					vecDecodedStimType(vecToBeDecodedTrials) = sOutSD.vecDecodedStimType(vecToBeDecodedTrials);
					vecStimType(vecToBeDecodedTrials) = sOutSD.vecStimType(vecToBeDecodedTrials);
					vecOutC(vecToBeDecodedTrials) = vecDecodedStimType(vecToBeDecodedTrials) == vecStimType(vecToBeDecodedTrials);
					
					%message
					for intTrial = vecToBeDecodedTrials
						intDecodedStimType = vecDecodedStimType(intTrial);
						intStimType = vecStimType(intTrial);
						if intStimType == intDecodedStimType
							strErr = '';
							vecCorrect(intTrial) = true;
						else
							strErr = '; Incorrect';
							vecCorrect(intTrial) = false;
						end
						%if verbose,fprintf('Decoded stimulus %d of %d; decoded type=%d, actual type=%d%s\n',intTrial,intTrials,vecDecodedStimType(intTrial),vecStimType(intTrial),strErr);end
					end
				end
				intTrialsC = sum(indSelectTheseTrials);
				intCorr = nansum(vecOutC);
				if verbose,fprintf('\nDecoding performance for C=%d: %d of %d [%.0f%%] correct\n',intContrast,intCorr,intTrialsC,(intCorr/intTrialsC)*100);end
			end
			
			%calculate heterogeneity
			%matResp [neuron x time]
			[matTrialResponse,cellSelectContrasts] = getTrialResponseData(ses,ses.structStim);
			[vecHeterogeneity,vecActivity] = calcMatRespHeteroGen(matTrialResponse(vecIncludeCells,:));
			
			structParams.intStopOffset = 0;
			structParams.intStartOffset = -round(ses.samplingFreq);
			matPreResp = getNeuronResponse(ses,1:numel(ses.neuron),1:length(ses.structStim.FrameOff),structParams);
			[vecHeteroPre,vecActivityPre] = calcMatRespHeteroGen(matPreResp(vecIncludeCells,:));
			
			%put in output
			matHD_Orientation(:,intNeuronPopulationType) = ses.structStim.Orientation;
			matHD_Contrast(:,intNeuronPopulationType) = ses.structStim.Contrast;
			matHD_ReactionTime(:,intNeuronPopulationType) = ses.structStim.vecTrialRespSecs;
			matHD_DecodingAccuracy(:,intNeuronPopulationType) = vecCorrect;
			matHD_HeterogeneityPre(:,intNeuronPopulationType) = vecHeteroPre;
			matHD_ActivityPre(:,intNeuronPopulationType) = vecActivityPre;
			matHD_HeterogeneityDuring(:,intNeuronPopulationType) = vecHeterogeneity;
			matHD_ActivityDuring(:,intNeuronPopulationType) = vecActivity;
		end
		%save data
		cellSaveDecoding2(intPopulation).matHD_Orientation = matHD_Orientation;
		cellSaveDecoding2(intPopulation).matHD_Contrast = matHD_Contrast;
		cellSaveDecoding2(intPopulation).matHD_ReactionTime = matHD_ReactionTime;
		cellSaveDecoding2(intPopulation).matHD_DecodingAccuracy = matHD_DecodingAccuracy;
		cellSaveDecoding2(intPopulation).matHD_HeterogeneityPre = matHD_HeterogeneityPre;
		cellSaveDecoding2(intPopulation).matHD_ActivityPre = matHD_ActivityPre;
		cellSaveDecoding2(intPopulation).matHD_HeterogeneityDuring = matHD_HeterogeneityDuring;
		cellSaveDecoding2(intPopulation).matHD_ActivityDuring = matHD_ActivityDuring;
	end
	
	%% save data structure
	vecClock = fix(clock);
	strDate = num2str(vecClock(1:3));
	strRecs = num2str(vecRecordings);
	strFile = ['data_aggregateSD2_' strSes '_' strrep(strDate,'    ','_') '_' strRecs(~isspace(strRecs))];
	%save(['D:\Data\Results\stimdetection\' strFile],'-v7.3');
	%if exist(['D:\Data\Results\stimdetection\' strFile '.mat'],'file')
	%	load(['D:\Data\Results\stimdetection\' strFile]);
	%end
	save(['D:\Data\Results\stimdetection\' strFile],'cellSaveDecoding2','-v7.3');
	cd(strOldDir);
end