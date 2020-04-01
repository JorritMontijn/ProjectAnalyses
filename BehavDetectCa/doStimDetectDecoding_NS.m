%function doStimDetectDecoding
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%get block data
clear all;
for intMouse=2:8
	clearvars -except intMouse
	if intMouse == 1
		strSes = '20140207';
	elseif intMouse == 2
		strSes = '20140314';
	elseif intMouse == 3
		strSes = '20140425';
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
	
	%parameters
	boolExcludeLocomotor = false;
	dblBaselineSecs = 0;
	intIterMax = 100;
	cellSaveDecoding = struct;
	
	%% GET RESPONSES AFTER COMBINING BLOCKS OF SAME POPULATION INTO ONE
	for intPopulation = 1:numel(cellMultiSes)
		%% GET RESPONSES AFTER COMBINING BLOCKS OF SAME POPULATION INTO ONE
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
		%structStim.FrameOff = structStim.FrameOn+1;
		
		%remove locomotion trials
		if boolExcludeLocomotor
			indLoco = false(1,numel(structStim.Orientation));
			for intTrial=1:numel(structStim.Orientation)
				vecMoveTimes = structStim.cellMoveSecs{intTrial};
				indLoco(intTrial) = sum(vecMoveTimes > structStim.SecsOn(intTrial) & vecMoveTimes < structStim.SecsOff(intTrial)) > 0;
			end
			structStim = remel(structStim,~indLoco);
			fprintf('Removed %d trials with locomotion\n',sum(indLoco));
		end
		
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
		intTrials = length(structStim.Orientation);
		intOris = length(vecOrientations);
		%cellMultiSes{intPopulation}.structStim = structStim;
		
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
		
		
		
		%% calculate hit-correlated activity enhancement
		%pre-allocate
		indHitTrials = logical(cellMultiSes{intPopulation}.structStim.vecTrialResponse);
		
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
		if dblBaselineSecs == 0,sParams.intPreBaselineRemoval = [];else sParams.intPreBaselineRemoval = round(ses.samplingFreq*dblBaselineSecs);end
		
		%pre-allocate output
		matRawBootstrappedDecodingOutput = nan(intOriNum,intContrastNum,2,2);
		matMeanBootstrappedDecodingOutput = nan(intContrastNum,2,2);
		matMixBootstrappedDecodingOutput = nan(intIterMax,intContrastNum,2);
		
		vecIncludeCells = true(1,numel(ses.neuron));
		
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
			fprintf('Now processing %s; iteration %d of %d for population %d; time is %s\n',strSes,intIter,intIterMax,intPopulation,getTime);
			
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
		cellSaveDecoding(intPopulation).matMixBootstrappedDecodingOutput = matMixBootstrappedDecodingOutput;
		cellSaveDecoding(intPopulation).matMeanBootstrappedDecodingOutput = matMeanBootstrappedDecodingOutput;
		cellSaveDecoding(intPopulation).matRawBootstrappedDecodingOutput = matRawBootstrappedDecodingOutput;
	end
	
	%% save data structure
	vecClock = fix(clock);
	strDate = num2str(vecClock(1:3));
	strFile = ['data_aggregateNSSD' strSes '_' strrep(strDate,'    ','_')];
	%save(['D:\Data\Results\stimdetection\' strFile],'-v7.3');
	%if exist(['D:\Data\Results\stimdetection\' strFile '.mat'],'file')
	%	load(['D:\Data\Results\stimdetection\' strFile]);
	%end
	save(['D:\Data\Results\stimdetection\' strFile],'cellSaveDecoding','-v7.3');
end