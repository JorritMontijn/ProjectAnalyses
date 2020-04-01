function [sContrast,sMetaData,sDecoding] = getStimDetectDataSD(ses,sStimAggregate)
	%UNTITLED3 Summary of this function goes here
	%   Detailed explanation goes here
	
	%get common variables
	intNeurons = numel(ses.neuron);
	
	% get indexing vectors for unique stimulus combinations
	sTypes = getStimulusTypes(ses,{'Orientation'});
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
	sTypes.matTypes = sTypes.matTypes(:,vecKeep);
	
	%general selection vectors
	vecCorrResp = logical(ses.structStim.vecTrialResponse);
	cellFieldsC{1} = 'Contrast';
	sTypesC = getStimulusTypes(ses,cellFieldsC);
	cellSelectC = getSelectionVectors(ses.structStim,sTypesC);
	vecContrasts = nan(1,length(cellSelectC));
	intNumContrasts = length(vecContrasts);
	
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
	
	
	%alter structStim so that all non-response trials have the same number
	%of frames as the mean of all response trials
	%loop through contrasts
	for intContrastIndex=1:length(cellSelectC)
		vecSelectC = cellSelectC{intContrastIndex};
		
		%get mean duration of responded trials for this contrast
		intMeanDur = vecMeanStimDur(intContrastIndex);
		
		%edit frame off vector
		ses.structStim.FrameOff(vecSelectC & ~vecCorrResp) = ses.structStim.FrameOn(vecSelectC & ~vecCorrResp) + intMeanDur - 1;
	end
	
	%select all other contrasts
	vecLikelihoodTrials = cellSelectC{end-1} | cellSelectC{end};
	sParams.vecLikelihoodTrials = vecLikelihoodTrials;
		
		
	%put in params structure
	sParams.cellSelect = cellSelect;
	sParams.sTypes = sTypes;
	sParams.verbose = true;
	
	%do decoding
	sDecoding = doStimDecoding(ses,sParams);
	
	%loop through contrasts and assign to field
	for intContrastIndex=1:length(cellSelectC)
		%select all other contrasts
		%vecLikelihoodTrials = ~cellSelectC{intContrastIndex};
		%sParams.vecLikelihoodTrials = vecLikelihoodTrials;
	
		%add metadata
		vecSelectC = cellSelectC{intContrastIndex};
		dblContrast = sTypesC.matTypes(intContrastIndex);
		vecContrasts(intContrastIndex) = dblContrast;
		
		%do decoding
		%sDecoding = doStimDecoding(ses,sParams);
	
		%get selection vectors
		vecRespTrials = vecCorrResp & vecSelectC;
		vecNoRespTrials = ~vecCorrResp & vecSelectC;
		
		%put in output
		sContrast(intContrastIndex).vecDecodedStimTypeR = sDecoding.vecDecodedStimType(vecRespTrials);
		sContrast(intContrastIndex).vecStimTypeR = sDecoding.vecStimType(vecRespTrials);
		sContrast(intContrastIndex).vecDecodedStimTypeN = sDecoding.vecDecodedStimType(vecNoRespTrials);
		sContrast(intContrastIndex).vecStimTypeN = sDecoding.vecStimType(vecNoRespTrials);
	end
	
	%add metadata to structure
	sMetaData.vecContrasts = vecContrasts;
	sMetaData.vecCorrResp = vecCorrResp;
	sMetaData.sTypesC = sTypesC;
	sMetaData.cellSelectC = cellSelectC;
end

