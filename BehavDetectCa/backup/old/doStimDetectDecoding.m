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
		
		% get hit/miss trials
		vecHits = ses.structStim.vecTrialResponse;
		
		%set base parameters
		sParams = struct;
		sParams.sTypes = sTypes;
		sParams.cellSelect = cellSelect;
		sParams.verbose = 0;
		sParams.intPreBaselineRemoval = round(ses.samplingFreq*dblBaselineSecs);
		
		for intNeuronPopulationType = 1:2
			%pre-allocate output
			matRawBootstrappedDecodingOutput = nan(intOriNum,intContrastNum,2,2);
			matMeanBootstrappedDecodingOutput = nan(intContrastNum,2,2);
			matMixBootstrappedDecodingOutput = nan(intIterMax,intContrastNum,2);
			
			if intNeuronPopulationType == 1
				%hit-correlated neurons
				vecIncludeCells = find(indSelectHitCorrelatedNeurons);
			elseif intNeuronPopulationType == 2
				%non-hit-correlated neurons
				vecIncludeCells = find(indSelectHitUncorrelatedNeurons);
			end
			
			%set neuronal populations
			sParams.vecIncludeCells = vecIncludeCells;
			
			%loop for hit/miss
			matDecodedOutput = nan(sTypes.vecNumTypes(1),sTypes.vecNumTypes(2),2,2); %[ori] x [contrast] x [hit/miss] x [likelihood-hit/miss]
			matHitMissTrials = nan(sTypes.vecNumTypes(1),sTypes.vecNumTypes(2),2); %number of hits/misses per stim type
			vecHitMiss = [1 0];
			for intHitLikelihood = 1:2
				%make likelihood
				indLikelihoodSelect = vecHits == vecHitMiss(intHitLikelihood);
				sParams.vecLikelihoodTrials = indLikelihoodSelect;
				
				for intContrast=1:intContrastNum
					%get contrast
					dblContrast = sTypes.matTypes(2,(intContrast-1)*intOriNum + 1);
					
					%select trials
					indContrastSelect = sTypes.matTypes(2,:) == dblContrast;
					for intOrientation = 1:intOriNum
						dblOri = sTypes.matTypes(1,intOrientation);
						
						%select trials
						indOriSelect = sTypes.matTypes(1,:) == dblOri;
						indSelectTrialType = indContrastSelect & indOriSelect;
						indTrialsOfType = cellSelect{indSelectTrialType};
						
						
						intHitNr = sum(indTrialsOfType & vecHits == 1);
						intMissNr = sum(indTrialsOfType & vecHits == 0);
						matHitMissTrials(intOrientation,intContrast,1) = intHitNr;
						matHitMissTrials(intOrientation,intContrast,2) = intMissNr;
						
						%put in structure
						sParams.vecDecodeTrials = indTrialsOfType;
						
						%decode
						sOut = doStimDecoding(ses,sParams);
						%split hit/miss trials
						
						vecHitDecodedStimType = sOut.vecDecodedStimType(indTrialsOfType & vecHits);
						vecHitStimType = sOut.vecStimType(indTrialsOfType & vecHits);
						
						vecMissDecodedStimType = sOut.vecDecodedStimType(indTrialsOfType & ~vecHits);
						vecMissStimType = sOut.vecStimType(indTrialsOfType & ~vecHits);
						
						dblHitCorrect = sum(vecHitDecodedStimType == vecHitStimType)/length(vecHitDecodedStimType);
						dblMissCorrect = sum(vecMissDecodedStimType == vecMissStimType)/length(vecMissDecodedStimType);
						matDecodedOutput(intOrientation,intContrast,1,intHitLikelihood) = dblHitCorrect; %hit trials
						matDecodedOutput(intOrientation,intContrast,2,intHitLikelihood) = dblMissCorrect; %miss trials
					end
				end
			end
			
			%put in output
			matRawBootstrappedDecodingOutput(:,:,:,:) = matDecodedOutput;
			matOut = nanmean(matDecodedOutput,1);
			matMeanBootstrappedDecodingOutput(:,:,:) = matOut;
			
			%% decode with mixed likelihood
			%select all trial types with at least 2 hits and 2 misses
			[vecOri,vecContrast] = find(matHitMissTrials(:,:,1)>1 & matHitMissTrials(:,:,2)>1);
			
			for intIter=1:intIterMax
				%msg
				fprintf('Now processing %s; iteration %d of %d for population %d; time is %s\n',strSes,intIter,intIterMax,intNeuronPopulationType,getTime);
				
				%create likelihood & which trials to decode
				vecLikelihoodTrials = false(1,length(vecHits));
				for intTrialType=1:length(vecOri)
					intContrast = vecContrast(intTrialType);
					intOrientation = vecOri(intTrialType);
					
					%get contrast & ori
					dblContrast = sTypes.matTypes(2,(intContrast-1)*intOriNum + 1);
					dblOri = sTypes.matTypes(1,intOrientation);
					
					%select trials
					indContrastSelect = sTypes.matTypes(2,:) == dblContrast;
					indOriSelect = sTypes.matTypes(1,:) == dblOri;
					
					indSelectTrialType = indContrastSelect & indOriSelect;
					indTrialsOfType = cellSelect{indSelectTrialType};
					
					vecHitOfType = find(indTrialsOfType & vecHits == 1);
					vecMissOfType = find(indTrialsOfType & vecHits == 0);
					
					%select likelihood trials
					intHitNr = matHitMissTrials(intOrientation,intContrast,1);
					intMissNr = matHitMissTrials(intOrientation,intContrast,2);
					intSelectTrials = min(intHitNr,intMissNr);
					
					vecHitSelect = randperm(intHitNr,intSelectTrials);
					vecMissSelect = randperm(intMissNr,intSelectTrials);
					
					vecHitLikelihood = vecHitOfType(vecHitSelect);
					vecMissLikelihood = vecMissOfType(vecMissSelect);
					
					
					vecLikelihoodTrials(vecHitLikelihood) = true;
					vecLikelihoodTrials(vecMissLikelihood) = true;
					
				end
				
				%pre-allocate output
				matMixDecodedOutput = nan(sTypes.vecNumTypes(1),sTypes.vecNumTypes(2),2); %[ori] x [contrast] x [hit/miss]
				sParams.vecLikelihoodTrials = vecLikelihoodTrials;
				sParams.vecDecodeTrials = 1:length(vecLikelihoodTrials);
				
				%decode
				sOut = doStimDecoding(ses,sParams);
				
				%put hit/miss trials in output matrix
				for intTrialType=1:length(vecOri)
					intContrast = vecContrast(intTrialType);
					intOrientation = vecOri(intTrialType);
					
					%get contrast & ori
					dblContrast = sTypes.matTypes(2,(intContrast-1)*intOriNum + 1);
					dblOri = sTypes.matTypes(1,intOrientation);
					
					%select trials
					indContrastSelect = sTypes.matTypes(2,:) == dblContrast;
					indOriSelect = sTypes.matTypes(1,:) == dblOri;
					
					indSelectTrialType = indContrastSelect & indOriSelect;
					indTrialsOfType = cellSelect{indSelectTrialType};
					
					vecHitOfType = find(indTrialsOfType & vecHits == 1);
					vecMissOfType = find(indTrialsOfType & vecHits == 0);
					
					%split hit/miss trials
					vecHitDecodedStimType = sOut.vecDecodedStimType(indTrialsOfType & vecHits & vecLikelihoodTrials);
					vecHitStimType = sOut.vecStimType(indTrialsOfType & vecHits & vecLikelihoodTrials);
					
					vecMissDecodedStimType = sOut.vecDecodedStimType(indTrialsOfType & ~vecHits & vecLikelihoodTrials);
					vecMissStimType = sOut.vecStimType(indTrialsOfType & ~vecHits & vecLikelihoodTrials);
					
					dblHitCorrect = sum(vecHitDecodedStimType == vecHitStimType)/length(vecHitDecodedStimType);
					dblMissCorrect = sum(vecMissDecodedStimType == vecMissStimType)/length(vecMissDecodedStimType);
					matMixDecodedOutput(intOrientation,intContrast,1) = dblHitCorrect; %hit trials
					matMixDecodedOutput(intOrientation,intContrast,2) = dblMissCorrect; %miss trials
				end
				matMixOut = nanmean(matMixDecodedOutput,1);
				matMixBootstrappedDecodingOutput(intIter,:,:) = matMixOut;
			end
			cellSaveDecoding(intNeuronPopulationType,intPopulation).matMixBootstrappedDecodingOutput = matMixBootstrappedDecodingOutput;
			cellSaveDecoding(intNeuronPopulationType,intPopulation).matMeanBootstrappedDecodingOutput = matMeanBootstrappedDecodingOutput;
			cellSaveDecoding(intNeuronPopulationType,intPopulation).matRawBootstrappedDecodingOutput = matRawBootstrappedDecodingOutput;
		end
	end
	
	%% save data structure
	vecClock = fix(clock);
	strDate = num2str(vecClock(1:3));
	strRecs = num2str(vecRecordings);
	strFile = ['data_aggregateSD' strSes '_' strrep(strDate,'    ','_') '_' strRecs(~isspace(strRecs))];
	%save(['D:\Data\Results\stimdetection\' strFile],'-v7.3');
	%if exist(['D:\Data\Results\stimdetection\' strFile '.mat'],'file')
	%	load(['D:\Data\Results\stimdetection\' strFile]);
	%end
	save(['D:\Data\Results\stimdetection\' strFile],'cellSave*','-v7.3');
	cd(strOldDir);
end