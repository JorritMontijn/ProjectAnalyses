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
	cellSaveDecoding = struct;
	
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
		structStim.FrameOff = structStim.FrameOn;
		
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
		if dblBaselineSecs == 0,sParams.intPreBaselineRemoval=[];else sParams.intPreBaselineRemoval = round(ses.samplingFreq*dblBaselineSecs);end
		indLikelihoodTemplate = false(1,length(cellSelectO));
		
		%pre-allocate output
		matHD_Orientation = nan(intTrials,1);
		matHD_Contrast = nan(intTrials,1);
		matHD_ReactionTime = nan(intTrials,1);
		matHD_DecodingAccuracy = nan(intTrials,1);
		matHD_HeterogeneityPre = nan(intTrials,1);
		matHD_ActivityPre = nan(intTrials,1);
		matHD_HeterogeneityDuring = nan(intTrials,1);
		matHD_ActivityDuring = nan(intTrials,1);
		
		vecIncludeCells = true(1,numel(ses.neuron));
		
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
			vecTrials = find(indSelectTheseTrials);
			
			%leave-one-repetition-out cross validation
			for intTrial=1:length(vecTrials)
				intRealTrial = vecTrials(intTrial);
				
				%make likelihood vector
				vecTheseLikelihoodTrials = indContrastSelect;
				vecTheseLikelihoodTrials(intRealTrial) = false;
				sParams.vecLikelihoodTrials = vecTheseLikelihoodTrials;
				sParams.vecDecodeTrials = intRealTrial;
				
				%do decoding
				sOutSD = doStimDecoding(ses,sParams);
				
				%put in output
				cellVecPost(intRealTrial) = sOutSD.cellVecPost(intRealTrial);
				vecDecodedStimType(intRealTrial) = sOutSD.vecDecodedStimType(intRealTrial);
				vecStimType(intRealTrial) = sOutSD.vecStimType(intRealTrial);
				vecOutC(intTrial) = vecDecodedStimType(intRealTrial) == vecStimType(intRealTrial);
				vecCorrect(intRealTrial) = vecDecodedStimType(intRealTrial) == vecStimType(intRealTrial);
			end
			intTrialsC = sum(indSelectTheseTrials);
			intCorr = nansum(vecOutC);
			if verbose,fprintf('\nDecoding performance for C=%d: %d of %d [%.0f%%] correct\n',intContrast,intCorr,intTrialsC,(intCorr/intTrialsC)*100);end
			
			
			%calculate heterogeneity
			%matResp [neuron x time]
			[matTrialResponse,cellSelectContrasts] = getTrialResponseData(ses,ses.structStim);
			[vecHeterogeneity,vecActivity] = calcMatRespHeteroGen(matTrialResponse(vecIncludeCells,:));
			
			structParams.intStopOffset = 0;
			structParams.intStartOffset = -round(ses.samplingFreq);
			matPreResp = getNeuronResponse(ses,1:numel(ses.neuron),1:length(ses.structStim.FrameOff),structParams);
			[vecHeteroPre,vecActivityPre] = calcMatRespHeteroGen(matPreResp(vecIncludeCells,:));
			
			%put in output
			matHD_Orientation(:) = ses.structStim.Orientation;
			matHD_Contrast(:) = ses.structStim.Contrast;
			matHD_ReactionTime(:) = ses.structStim.vecTrialRespSecs;
			matHD_DecodingAccuracy(:) = vecCorrect;
			matHD_HeterogeneityPre(:) = vecHeteroPre;
			matHD_ActivityPre(:) = vecActivityPre;
			matHD_HeterogeneityDuring(:) = vecHeterogeneity;
			matHD_ActivityDuring(:) = vecActivity;
			
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
		strFile = ['data_aggregateNSSD2_' strSes '_' strrep(strDate,'    ','_')];
		%save(['D:\Data\Results\stimdetection\' strFile],'-v7.3');
		%if exist(['D:\Data\Results\stimdetection\' strFile '.mat'],'file')
		%	load(['D:\Data\Results\stimdetection\' strFile]);
		%end
		save(['D:\Data\Results\stimdetection\' strFile],'cellSaveDecoding2','-v7.3');
	end
end