function matStimResponse = doMatRespTransform(matResponse,cellSelect)
	%doMatRespTransform Transforms raw response matrix to [stimtype x stimrep x neuron]
	%Syntax: matStimResponse = doMatRespTransform(matResponse,cellSelect)
	
	% loop through stimulus types to get nr of reps
	intMaxReps = 0;
	intStimTypes = length(cellSelect);
	for intS=1:intStimTypes
		intMaxReps = max(intMaxReps,sum(cellSelect{intS}));
	end
	
	% retrieve responses
	intNeurons = size(matResponse,1);
	matStimResponse = nan(intStimTypes,intMaxReps,intNeurons);
	for intStimIndex=1:intStimTypes
		vecSelect = cellSelect{intStimIndex};
		matThisStim = matResponse(:,vecSelect);
		intReps = size(matThisStim,2);
		if intReps < intMaxReps
			matThisStim = [matThisStim nan(intNeurons,intMaxReps-intReps)];
		end
		matStimResponse(intStimIndex,:,:) = shiftdim(matThisStim,1);
	end
end

