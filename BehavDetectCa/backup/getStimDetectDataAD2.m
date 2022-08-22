function [sContrast,sMetaData] = getStimDetectDataAD2(ses,sStimAggregate)
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
	vecContrastRef = [0 0.005 0.02 0.08 0.32 1];
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
	
	%pre-allocate
	cellContrastPR = cell(1,6);
	cellContrastNPR = cell(1,6);
	cellContrastPN = cell(1,6);
	cellContrastNPN = cell(1,6);
	cellContrastPR_index = cell(1,6);
	cellContrastNPR_index = cell(1,6);
	cellContrastPN_index = cell(1,6);
	cellContrastNPN_index = cell(1,6);
	
	%loop through contrasts
	for intContrastIndex=1:length(cellSelectC)
		vecSelectC = cellSelectC{intContrastIndex};
		dblContrast = sTypesC.matTypes(intContrastIndex);
		intContrastRefIndex = find(vecContrastRef==dblContrast);
		vecContrasts(intContrastRefIndex) = dblContrast;
		vecCorrRespC = logical(ses.structStim.vecTrialResponse) & vecSelectC;
		
		%pre-allocate output matrices
		intPreAllocRespN = intNeurons * sum(vecCorrRespC);
		intPreAllocNoRespN = intNeurons * sum(~vecCorrRespC);
		vecPR = nan(intPreAllocRespN);
		vecPR_index = nan(intPreAllocRespN);
		vecNPR = nan(intPreAllocRespN);
		vecNPR_index = nan(intPreAllocRespN);
		vecPN = nan(intPreAllocNoRespN);
		vecPN_index = nan(intPreAllocNoRespN);
		vecNPN = nan(intPreAllocNoRespN);
		vecNPN_index = nan(intPreAllocNoRespN);
		intCounterPR = 1;
		intCounterNPR = 1;
		intCounterPN = 1;
		intCounterNPN = 1;
		
		%plot per type
		for intStimType=1:intTypes
			indPrefNeurons = vecPrefStim==intStimType;
			intPrefNeurons = sum(indPrefNeurons);
			vecPrefNeurons = find(indPrefNeurons);
			indNonPrefNeurons = ~indPrefNeurons;
			intNonPrefNeurons = sum(indNonPrefNeurons);
			vecNonPrefNeurons = find(indNonPrefNeurons);
			if intPrefNeurons>0
				
				%correct trials
				%get selection vectors
				vecRespTrials = structStimCorrs.cellSelect{intStimType} & vecCorrResp & vecSelectC;
				vecRespOn = ses.structStim.FrameOn(vecRespTrials);
				vecRespOff = ses.structStim.FrameOff(vecRespTrials);
				
				for intTrial=1:length(vecRespOn)
					%pref + corr
					vecPrefAct = mean(matActivity(indPrefNeurons,vecRespOn(intTrial):vecRespOff(intTrial)),2);
					vecPR(intCounterPR:(intCounterPR+intPrefNeurons-1)) = vecPrefAct(:);
					vecPR_index(intCounterPR:(intCounterPR+intPrefNeurons-1)) = vecPrefNeurons;
					intCounterPR = intCounterPR + intPrefNeurons;
					
					%non-pref + corr
					vecNonPrefAct = mean(matActivity(indNonPrefNeurons,vecRespOn(intTrial):vecRespOff(intTrial)),2);
					vecNPR(intCounterNPR:(intCounterNPR+intNonPrefNeurons-1)) = vecNonPrefAct(:);
					vecNPR_index(intCounterNPR:(intCounterNPR+intNonPrefNeurons-1)) = vecNonPrefNeurons;
					intCounterNPR = intCounterNPR + intNonPrefNeurons;
				end
				
				%incorrect trials
				%get selection vector
				vecNoRespTrials = structStimCorrs.cellSelect{intStimType} & ~vecCorrResp & vecSelectC;
				vecNoRespOn = ses.structStim.FrameOn(vecNoRespTrials);
				
				if isnan(vecMeanStimDur(intContrastIndex))
					vecNoRespOff = ses.structStim.FrameOff(vecNoRespTrials);
				else
					vecNoRespOff = vecNoRespOn + vecMeanStimDur(intContrastIndex);
				end
				
				for intTrial=1:length(vecNoRespOn)
					%pref + incorr
					vecPrefAct = mean(matActivity(indPrefNeurons,vecNoRespOn(intTrial):vecNoRespOff(intTrial)),2);
					vecPN(intCounterPN:(intCounterPN+intPrefNeurons-1)) = vecPrefAct(:);
					vecPN_index(intCounterPN:(intCounterPN+intPrefNeurons-1)) = vecPrefNeurons;
					intCounterPN = intCounterPN + intPrefNeurons;
					
					%non-pref + incorr
					vecNonPrefAct = mean(matActivity(indNonPrefNeurons,vecNoRespOn(intTrial):vecNoRespOff(intTrial)),2);
					vecNPN(intCounterNPN:(intCounterNPN+intNonPrefNeurons-1)) = vecNonPrefAct(:);
					vecNPN_index(intCounterNPN:(intCounterNPN+intNonPrefNeurons-1)) = vecNonPrefNeurons;
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
		[indNan]=find(isnan(vecPR_index),1);
		vecPR_index = vecPR_index(1:(indNan-1));
		[indNan]=find(isnan(vecNPR_index),1);
		vecNPR_index = vecNPR_index(1:(indNan-1));
		[indNan]=find(isnan(vecPN_index),1);
		vecPN_index = vecPN_index(1:(indNan-1));
		[indNan]=find(isnan(vecNPN_index),1);
		vecNPN_index = vecNPN_index(1:(indNan-1));
		
		%add to output
		cellContrastPR{intContrastRefIndex} = vecPR;
		cellContrastNPR{intContrastRefIndex} = vecNPR;
		cellContrastPN{intContrastRefIndex} = vecPN;
		cellContrastNPN{intContrastRefIndex} = vecNPN;
		cellContrastPR_index{intContrastRefIndex} = vecPR_index;
		cellContrastNPR_index{intContrastRefIndex} = vecNPR_index;
		cellContrastPN_index{intContrastRefIndex} = vecPN_index;
		cellContrastNPN_index{intContrastRefIndex} = vecNPN_index;
	end
	
	%add data to structure
	sContrast.cellContrastPR = cellContrastPR;
	sContrast.cellContrastNPR = cellContrastNPR;
	sContrast.cellContrastPN = cellContrastPN;
	sContrast.cellContrastNPN = cellContrastNPN;
	sContrast.cellContrastPR_index = cellContrastPR_index;
	sContrast.cellContrastNPR_index = cellContrastNPR_index;
	sContrast.cellContrastPN_index = cellContrastPN_index;
	sContrast.cellContrastNPN_index = cellContrastNPN_index;
	
	%add metadata to structure
	sMetaData.vecContrasts = vecContrasts;
	sMetaData.structStimCorrs = structStimCorrs;
	sMetaData.vecCorrResp = vecCorrResp;
	sMetaData.sTypesC = sTypesC;
	sMetaData.cellSelectC = cellSelectC;
	sMetaData.vecPrefStim = vecPrefStim;
end

