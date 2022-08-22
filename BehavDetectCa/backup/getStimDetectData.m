function [sContrast,sMetaData] = getStimDetectData(ses)
	%UNTITLED3 Summary of this function goes here
	%   Detailed explanation goes here
	
	
	%get common variables
	intNeurons = numel(ses.neuron);
	
	% get indexing vectors for unique stimulus combinations
	sTypes = getStimulusTypes(ses);
	cellSelect = getSelectionVectors(ses.structStim,sTypes);
	
	%remove empty types
	intTypes = length(cellSelect);
	vecKeep = true(1,intTypes);
	for intType = 1:intTypes
		if sum(cellSelect{intType}) == 0
			vecKeep(intType) = false;
		end
	end
	cellSelect = cellSelect(vecKeep);
	sTypes = getStimulusTypes(ses);
	sTypes.matTypes = sTypes.matTypes(:,vecKeep);
	
	%transform structure-based data to raw dFoF matrix
	intFrames = length(ses.neuron(1).dFoF);
	matActivity = zeros(intNeurons,intFrames);
	for intNeuron=1:intNeurons
		matActivity(intNeuron,:) = ses.neuron(intNeuron).dFoF;
	end
	%{
	%get pref stim per neuron
	structStimCorrs = calcStimCorrs(ses);
	matSignalCorrs = structStimCorrs.matSignalCorrs;
	matNoiseCorrs = structStimCorrs.matNoiseCorrs;
	matStimResponse = structStimCorrs.matStimResponse;
	matSignalResponse = structStimCorrs.matSignalResponse;
	[vecStimAct,vecPrefStim]=max(matSignalResponse,[],1);
	%}
	
	%general selection vectors
	vecCorrResp = logical(ses.structStim.vecTrialResponse);
	cellFieldsC{1} = 'Contrast';
	sTypesC = getStimulusTypes(ses,cellFieldsC);
	cellSelectC = getSelectionVectors(ses.structStim,sTypesC);
	vecContrasts = nan(1,length(cellSelectC));
	
	%select only full contrast stimuli
	vecTrialsFullC = cellSelectC{end} | cellSelectC{end-1};
	cellSubFields = fieldnames(ses.structStim);
	sesFullC = ses;
	for intField=1:length(cellSubFields)
		strField = cellSubFields{intField};
		sesFullC.structStim.(strField) = sesFullC.structStim.(strField)(vecTrialsFullC);
	end
	
	%get pref stim per neuron
	structStimCorrsFullC = calcStimCorrs(sesFullC);
	matSignalResponse = structStimCorrsFullC.matSignalResponse;
	[vecStimAct,vecPrefStim]=max(matSignalResponse,[],1);
	structStimCorrs = calcStimCorrs(ses);
	
	%create structures
	sContrast = struct;
	sMetaData = struct;
		
	%loop through contrasts
	for intContrastIndex=1:length(cellSelectC)
		vecSelectC = cellSelectC{intContrastIndex};
		dblContrast = sTypesC.matTypes(intContrastIndex);
		vecContrasts(intContrastIndex) = dblContrast;
		vecCorrRespC = logical(ses.structStim.vecTrialResponse) & vecSelectC;
		
		%pre-allocate output matrices
		vecWindow = round([-3*ses.samplingFreq 5*ses.samplingFreq]);
		intWindowFrames = vecWindow(2) - vecWindow(1) + 1;
		intPreAllocRespN = intNeurons * sum(vecCorrRespC);
		intPreAllocNoRespN = intNeurons * sum(~vecCorrRespC);
		matPR = nan(intPreAllocRespN,intWindowFrames);
		matNPR = nan(intPreAllocRespN,intWindowFrames);
		matPN = nan(intPreAllocNoRespN,intWindowFrames);
		matNPN = nan(intPreAllocNoRespN,intWindowFrames);
		intCounterPR = 1;
		intCounterNPR = 1;
		intCounterPN = 1;
		intCounterNPN = 1;
		
		%plot per type
		for intStimType=1:intTypes
			vecPrefNeurons = vecPrefStim==intStimType;
			intPrefNeurons = sum(vecPrefNeurons);
			vecNonPrefNeurons = ~vecPrefNeurons;
			intNonPrefNeurons = sum(vecNonPrefNeurons);
			if intPrefNeurons>0
				
				%correct trials
				%get selection vectors
				vecRespTrials = structStimCorrs.cellSelect{intStimType} & vecCorrResp & vecSelectC;
				vecRespOn = ses.structStim.FrameOn(vecRespTrials) + vecWindow(1);
				vecRespOff = ses.structStim.FrameOn(vecRespTrials) + vecWindow(2);
				
				for intTrial=1:length(vecRespOn)
					%pref + corr
					matPR(intCounterPR:(intCounterPR+intPrefNeurons-1),:) = matActivity(vecPrefNeurons,vecRespOn(intTrial):vecRespOff(intTrial));
					intCounterPR = intCounterPR + intPrefNeurons;
					
					%non-pref + corr
					matNPR(intCounterNPR:(intCounterNPR+intNonPrefNeurons-1),:) = matActivity(vecNonPrefNeurons,vecRespOn(intTrial):vecRespOff(intTrial));
					intCounterNPR = intCounterNPR + intNonPrefNeurons;
				end
				
				
				%incorrect trials
				%get selection vector
				vecNoRespTrials = structStimCorrs.cellSelect{intStimType} & ~vecCorrResp & vecSelectC;
				vecNoRespOn = ses.structStim.FrameOn(vecNoRespTrials) + vecWindow(1);
				vecNoRespOff = ses.structStim.FrameOn(vecNoRespTrials) + vecWindow(2);
				
				for intTrial=1:length(vecNoRespOn)
					%pref + incorr
					matPN(intCounterPN:(intCounterPN+intPrefNeurons-1),:) = matActivity(vecPrefNeurons,vecNoRespOn(intTrial):vecNoRespOff(intTrial));
					intCounterPN = intCounterPN + intPrefNeurons;
					
					%non-pref + incorr
					matNPN(intCounterNPN:(intCounterNPN+intNonPrefNeurons-1),:) = matActivity(vecNonPrefNeurons,vecNoRespOn(intTrial):vecNoRespOff(intTrial));
					intCounterNPN = intCounterNPN + intNonPrefNeurons;
				end
			end
		end
		%{
		%remove trailing nans
		[row,col]=find(isnan(matPR));
		matPR = matPR(1:(row-1),:);
		[row,col]=find(isnan(matNPR));
		matNPR = matNPR(1:(row-1),:);
		[row,col]=find(isnan(matPN));
		matPN = matPN(1:(row-1),:);
		[row,col]=find(isnan(matNPN));
		matNPN = matNPN(1:(row-1),:);
		%}
		%add data to structure
		sContrast(intContrastIndex).matPR = matPR;
		sContrast(intContrastIndex).matNPR = matNPR;
		sContrast(intContrastIndex).matPN = matPN;
		sContrast(intContrastIndex).matNPN = matNPN;
	end
	%add metadata to structure
	sMetaData.vecContrasts = vecContrasts;
	sMetaData.structStimCorrs = structStimCorrs;
	sMetaData.vecCorrResp = vecCorrResp;
	sMetaData.sTypesC = sTypesC;
	sMetaData.cellSelectC = cellSelectC;
	sMetaData.vecWindow = vecWindow;
	sMetaData.intWindowFrames = intWindowFrames;
end

