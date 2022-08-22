function [sContrast,sMetaData] = getStimDetectDataResps(ses)
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
	
	%general selection vectors
	vecResp = logical(ses.structStim.vecTrialResponse);
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
	
	%pre-allocate no stim variables
	intPreAllocNoStimN = intNeurons * length(ses.structStim.cellRespPulses) * 30;
	vecWindow = round([-3*ses.samplingFreq 5*ses.samplingFreq]);
	vecRemWindow = round([-5*ses.samplingFreq 5*ses.samplingFreq]);
	intWindowFrames = vecWindow(2) - vecWindow(1) + 1;
	matNS = nan(intPreAllocNoStimN,intWindowFrames); %no stimulus
	intCounterNS = 1;
	
	%get resps outside stimulus
	for intTrial=1:length(ses.structStim.cellRespPulses)
		vecRespPulses = ses.structStim.cellRespPulses{intTrial};
		intStimOn = ses.structStim.FrameOn(intTrial);
		intStimOff = ses.structStim.FrameOff(intTrial);
		intRemStart = intStimOn + vecRemWindow(1);
		intRemStop = intStimOff + vecRemWindow(2);
		vecSelect = vecRespPulses > intRemStop | vecRespPulses < intRemStart;
		vecRespPulses = vecRespPulses(vecSelect);
		
		%loop through valid responses
		for intRespPulse=vecRespPulses
			if vecWindow(1) + intRespPulse < 1, continue;end
			matNS(intCounterNS:(intCounterNS+intNeurons-1),:) = matActivity(:,(intRespPulse + vecWindow(1)):(intRespPulse + vecWindow(2)));
			intCounterNS = intCounterNS + 1;
		end
	end
	
	%remove trailing nans
	[row,col]=find(isnan(matNS));
	matNS = matNS(1:(row-1),:);
	
	%loop through contrasts
	for intContrastIndex=1:length(cellSelectC)
		vecSelectC = cellSelectC{intContrastIndex};
		dblContrast = sTypesC.matTypes(intContrastIndex);
		vecContrasts(intContrastIndex) = dblContrast;
		vecCorrRespC = logical(ses.structStim.vecTrialResponse) & vecSelectC;
		
		%pre-allocate output matrices
		intPreAllocStimN = intNeurons * sum(vecCorrRespC);
		
		matPR = nan(intPreAllocStimN,intWindowFrames); %preferred pop stim
		matNPR = nan(intPreAllocStimN,intWindowFrames); %non-preferred pop stim
		
		intCounterPR = 1;
		intCounterNPR = 1;
		
		
		%plot per type
		for intStimType=1:intTypes
			vecPrefNeurons = vecPrefStim==intStimType;
			intPrefNeurons = sum(vecPrefNeurons);
			vecNonPrefNeurons = ~vecPrefNeurons;
			intNonPrefNeurons = sum(vecNonPrefNeurons);
			if intPrefNeurons>0
				
				%correct trials
				%get selection vectors
				vecRespTrials = structStimCorrs.cellSelect{intStimType} & vecResp & vecSelectC;
				vecRespOn = ses.structStim.vecTrialRespPulses(vecRespTrials) + vecWindow(1);
				vecRespOff = ses.structStim.vecTrialRespPulses(vecRespTrials) + vecWindow(2);
				
				for intTrial=1:length(vecRespOn)
					%pref + corr
					matPR(intCounterPR:(intCounterPR+intPrefNeurons-1),:) = matActivity(vecPrefNeurons,vecRespOn(intTrial):vecRespOff(intTrial));
					intCounterPR = intCounterPR + intPrefNeurons;
					
					%non-pref + corr
					matNPR(intCounterNPR:(intCounterNPR+intNonPrefNeurons-1),:) = matActivity(vecNonPrefNeurons,vecRespOn(intTrial):vecRespOff(intTrial));
					intCounterNPR = intCounterNPR + intNonPrefNeurons;
				end
			end
		end
		
		%remove trailing nans
		[row,col]=find(isnan(matPR));
		matPR = matPR(1:(row-1),:);
		[row,col]=find(isnan(matNPR));
		matNPR = matNPR(1:(row-1),:);

		
		%add data to structure
		sContrast(intContrastIndex).matPR = matPR;
		sContrast(intContrastIndex).matNPR = matNPR;
		sContrast(intContrastIndex).matNS = matNS;
	end
	%add metadata to structure
	sMetaData.vecContrasts = vecContrasts;
	sMetaData.structStimCorrs = structStimCorrs;
	sMetaData.vecCorrResp = vecResp;
	sMetaData.sTypesC = sTypesC;
	sMetaData.cellSelectC = cellSelectC;
	sMetaData.vecWindow = vecWindow;
	sMetaData.intWindowFrames = intWindowFrames;
end

