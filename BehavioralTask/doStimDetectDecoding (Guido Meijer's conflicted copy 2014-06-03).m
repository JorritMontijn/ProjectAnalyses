%function doStimDetectDecoding
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%parameters
dblBaselineSecs = 0.1;
intIterMax = 100;
cellSaveDecoding = struct;

%input
ses = cellMultiSes{end};
ses.structStim.Orientation = mod(ses.structStim.Orientation,180);

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
	matRawBootstrappedDecodingOutput = nan(intIterMax,intOriNum,intContrastNum,2,2);
	matMeanBootstrappedDecodingOutput = nan(intIterMax,intContrastNum,2,2);
	matMixBootstrappedDecodingOutput = nan(intIterMax,intContrastNum,2);
	
	for intIter=1:intIterMax
		%msg
		fprintf('Now processing iteration %d of %d for population %d; time is %s\n',intIter,intIterMax,intNeuronPopulationType,getTime);
		if intNeuronPopulationType == 1
			%blob neurons
			vecIncludeCells = find(indSelectHitCorrelatedNeurons);
		elseif intNeuronPopulationType == 2
			%random set of neurons same size as blob
			vecIncludeCells = randperm(numel(ses.neuron),sum(indSelectHitCorrelatedNeurons));
		end
		
		%set neuronal populations
		sParams.vecIncludeCells = vecIncludeCells;
		
		%loop for hit/miss
		matDecodedOutput = nan(sTypes.vecNumTypes(1),sTypes.vecNumTypes(2),2,2); %[ori] x [contrast] x [hit/miss] x [likelihood-hit/miss]
		matHitMissTrials = nan(sTypes.vecNumTypes(1),sTypes.vecNumTypes(2),2);
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
		matRawBootstrappedDecodingOutput(intIter,:,:,:,:) = matDecodedOutput;
		matOut = nanmean(matDecodedOutput,1);
		matMeanBootstrappedDecodingOutput(intIter,:,:,:) = matOut;
		
		%% decode with mixed likelihood
		%select all trial types with at least 2 hits and 2 misses
		[vecOri,vecContrast] = find(matHitMissTrials(:,:,1)>1 & matHitMissTrials(:,:,2)>1);
		
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
	cellSaveDecoding(intNeuronPopulationType).matMixBootstrappedDecodingOutput = matMixBootstrappedDecodingOutput;
	cellSaveDecoding(intNeuronPopulationType).matMeanBootstrappedDecodingOutput = matMeanBootstrappedDecodingOutput;
	cellSaveDecoding(intNeuronPopulationType).matRawBootstrappedDecodingOutput = matRawBootstrappedDecodingOutput;
end

%% save data structure
vecClock = fix(clock);
strDate = num2str(vecClock(1:3));
strRecs = num2str(vecRecordings);
strFile = ['data_aggregateAD' strSes '_' strrep(strDate,'    ','_') '_' strRecs(~isspace(strRecs))];
%save(['D:\Data\Results\stimdetection\' strFile],'-v7.3');
if exist(['D:\Data\Results\stimdetection\' strFile '.mat'],'file')
	load(['D:\Data\Results\stimdetection\' strFile]);
end
save(['D:\Data\Results\stimdetection\' strFile],'cellSave*','-v7.3');
cd(strOldDir);