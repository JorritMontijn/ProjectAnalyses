function [matRespNSR,vecStimTypes,vecUnique] = getStimulusResponses(matResp,vecStimTypeList)
	%[Neuron x Stimulus x Repetition]
	
	%neuron x trial
   intNeurons = size(matResp,1);
   intTrials = size(matResp,2);
   [vecStimTypes,vecUnique] = label2idx(vecStimTypeList);
   vecStimTypes = vecStimTypes(:)';
   intStimTypes = numel(unique(vecStimTypes));
   intRepetitions = sum(vecStimTypes==vecStimTypes(1));
   if (intStimTypes * intRepetitions) ~= intTrials
	   error([mfilename ':InconsistentNumberOfRepetitions'],'Number of repetitions is inconsistent');
   end
   
   matRespNSR = nan(intNeurons,intStimTypes,intRepetitions);
   for intStimType=1:intStimTypes
	   matRespNSR(:,intStimType,:) = matResp(:,vecStimTypes==intStimType);
   end
   
end
