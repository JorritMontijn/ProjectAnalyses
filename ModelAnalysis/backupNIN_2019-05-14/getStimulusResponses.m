function [matRespNSR,vecStimTypes,vecUnique] = getStimulusResponses(matResp,vecStimTypeList)
	%[Neuron x Stimulus x Repetition]
	
	%neuron x trial
   intNeurons = size(matResp,1);
   intTrials = size(matResp,2);
   [vecStimTypes,vecUnique,vecCounts] = label2idx(vecStimTypeList);
   vecStimTypes = vecStimTypes(:)';
   intStimTypes = numel(unique(vecStimTypes));
   intRepetitions = sum(vecStimTypes==vecStimTypes(1));
   if (intStimTypes * intRepetitions) ~= intTrials
	   warning([mfilename ':InconsistentNumberOfRepetitions'],'Number of repetitions is inconsistent');
	   intRepetitions = max(vecCounts);
   end
   
   matRespNSR = nan(intNeurons,intStimTypes,intRepetitions);
   
	   for intStimType=1:intStimTypes
	    indSelect = vecStimTypes==intStimType;
		intThisReps = sum(indSelect);
		for intNeuron = 1:intNeurons
	  
	   matRespNSR(intNeuron,intStimType,1:intThisReps) = matResp(:,indSelect);
   end
   
end
