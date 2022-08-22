function [sContrast,sMetaData] = getStimDetectDataPS(ses,sStimAggregate)
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
	vecCorrResp = logical(ses.structStim.vecTrialResponse);
	cellFieldsC{1} = 'Contrast';
	sTypesC = getStimulusTypes(ses,cellFieldsC);
	cellSelectC = getSelectionVectors(ses.structStim,sTypesC);
	vecContrasts = nan(1,length(cellSelectC));
	intNumContrasts = length(vecContrasts);
	
	
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
	
	
	%check for stim aggregate; if present, then calculate mean number of
	%frames for detected stimuli
	if exist('sStimAggregate','var') && isstruct(sStimAggregate)
		cellFieldsC_SA{1} = 'Contrast';
		sTypesC_SA = getStimulusTypes(ses,cellFieldsC_SA);
		cellSelectC_SA = getSelectionVectors(sStimAggregate,sTypesC_SA);
		
		vecMeanStimDur = nan(1,intNumContrasts);
		%get responded trials
		vecResponded = logical(sStimAggregate.vecTrialResponse);
		for intContrastIndex=1:intNumContrasts
			%get target trials
			vecSelectTrials = vecResponded & cellSelectC_SA{intContrastIndex};
			
			%calculate mean stim duration
			vecMeanStimDur(intContrastIndex) = round(mean(sStimAggregate.FrameOff(vecSelectTrials)-sStimAggregate.FrameOn(vecSelectTrials)));
		end
	else
		vecMeanStimDur = nan(1,intNumContrasts);
	end
	
	
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
		intPreAllocRespN = intNeurons * sum(vecCorrRespC);
		intPreAllocNoRespN = intNeurons * sum(~vecCorrRespC);
		vecPR = nan(1,intPreAllocRespN);
		vecNPR = nan(1,intPreAllocRespN);
		vecPN = nan(1,intPreAllocNoRespN);
		vecNPN = nan(1,intPreAllocNoRespN);
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
				vecRespOn = ses.structStim.FrameOn(vecRespTrials);
				dblBaselineSecs = 1;
				vecBaseRespStart = vecRespOn - round(ses.samplingFreq*dblBaselineSecs) - 1;
				vecBaseRespStop = vecRespOn - 1;
				
				for intTrial=1:length(vecRespOn)
					%pref + corr
					vecPrefAct = mean(matActivity(vecPrefNeurons,vecBaseRespStart(intTrial):vecBaseRespStop(intTrial)),2);
					vecPR(intCounterPR:(intCounterPR+intPrefNeurons-1)) = vecPrefAct(:);
					intCounterPR = intCounterPR + intPrefNeurons;
					
					%non-pref + corr
					vecNonPrefAct = mean(matActivity(vecNonPrefNeurons,vecBaseRespStart(intTrial):vecBaseRespStop(intTrial)),2);
					vecNPR(intCounterNPR:(intCounterNPR+intNonPrefNeurons-1)) = vecNonPrefAct(:);
					intCounterNPR = intCounterNPR + intNonPrefNeurons;
				end
				
				
				%incorrect trials
				%get selection vector
				vecNoRespTrials = structStimCorrs.cellSelect{intStimType} & ~vecCorrResp & vecSelectC;
				vecNoRespOn = ses.structStim.FrameOn(vecNoRespTrials);
				
				dblBaselineSecs = 1;
				vecBaseNoRespStart = vecNoRespOn - round(ses.samplingFreq*dblBaselineSecs) - 1;
				vecBaseNoRespStop = vecNoRespOn - 1;
				
				for intTrial=1:length(vecNoRespOn)
					%pref + incorr
					vecPrefAct = mean(matActivity(vecPrefNeurons,vecBaseNoRespStart(intTrial):vecBaseNoRespStop(intTrial)),2);
					vecPN(intCounterPN:(intCounterPN+intPrefNeurons-1)) = vecPrefAct(:);
					intCounterPN = intCounterPN + intPrefNeurons;
					
					%non-pref + incorr
					vecNonPrefAct = mean(matActivity(vecNonPrefNeurons,vecBaseNoRespStart(intTrial):vecBaseNoRespStop(intTrial)),2);
					vecNPN(intCounterNPN:(intCounterNPN+intNonPrefNeurons-1)) = vecNonPrefAct(:);
					intCounterNPN = intCounterNPN + intNonPrefNeurons;
				end
			end
		end
		
		%remove trailing nans
		[indNan]=find(isnan(vecPR),1);
		vecPR = vecPR(1:(indNan-1));
		[indNan]=find(isnan(vecNPR),1);
		vecNPR = vecNPR(1:(indNan-1));
		[indNan]=find(isnan(vecPN),1);
		vecPN = vecPN(1:(indNan-1));
		[indNan]=find(isnan(vecNPN),1);
		vecNPN = vecNPN(1:(indNan-1));
		
		%add data to structure
		sContrast(intContrastIndex).vecPR = vecPR;
		sContrast(intContrastIndex).vecNPR = vecNPR;
		sContrast(intContrastIndex).vecPN = vecPN;
		sContrast(intContrastIndex).vecNPN = vecNPN;
	end
	%add metadata to structure
	sMetaData.vecContrasts = vecContrasts;
	sMetaData.structStimCorrs = structStimCorrs;
	sMetaData.vecCorrResp = vecCorrResp;
	sMetaData.sTypesC = sTypesC;
	sMetaData.cellSelectC = cellSelectC;
end

