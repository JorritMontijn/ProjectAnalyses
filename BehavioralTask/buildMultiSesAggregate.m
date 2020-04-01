function sSesAggregate = buildMultiSesAggregate(ses,sSesAggregate)
	%UNTITLED2 Summary of this function goes here
	%   Detailed explanation goes here
	
	if ~exist('sSesAggregate','var') || isempty(sSesAggregate)
		sSesAggregate = ses;
	else
		%get starting values from previous session
		intLastFrame = length(sSesAggregate.neuron(1).dFoF);
		dblLastSec = (sSesAggregate.structStim.SecsOff(end)/sSesAggregate.structStim.FrameOff(end))*length(sSesAggregate.neuron(1).dFoF);
		
		%aggregate structStim
		cellFields = fieldnames(ses.structStim);
		intNumTrials = length(ses.structStim.FrameOn);
		for intField=1:length(cellFields)
			strField = cellFields{intField};
			if ~isfield(sSesAggregate.structStim,strField),continue;end
			if ~isempty(strfind(strField,'Secs'))
				if iscell(sSesAggregate.structStim.(strField))
					for intExtraTrial=1:intNumTrials
						sSesAggregate.structStim.(strField){end+1} = ses.structStim.(strField){intExtraTrial}+dblLastSec;
					end
				else
					sSesAggregate.structStim.(strField)((end+1):(end+intNumTrials)) = ses.structStim.(strField)+dblLastSec;
				end
			elseif ~isempty(strfind(strField,'Pulse')) || strcmp(strField,'FrameOn') || strcmp(strField,'FrameOff')
				if iscell(sSesAggregate.structStim.(strField))
					for intExtraTrial=1:intNumTrials
						sSesAggregate.structStim.(strField){end+1} = ses.structStim.(strField){intExtraTrial}+intLastFrame;
					end
				else
					sSesAggregate.structStim.(strField)((end+1):(end+intNumTrials)) = ses.structStim.(strField)+intLastFrame;
				end
			else
				sSesAggregate.structStim.(strField)((end+1):(end+intNumTrials)) = ses.structStim.(strField);
			end
		end
		
		%aggregate neuron
		for intNeuron=1:numel(ses.neuron)
			sSesAggregate.neuron(intNeuron).F = [sSesAggregate.neuron(intNeuron).F ses.neuron(intNeuron).F];
			sSesAggregate.neuron(intNeuron).npF = [sSesAggregate.neuron(intNeuron).npF ses.neuron(intNeuron).npF];
			sSesAggregate.neuron(intNeuron).dFoF = [sSesAggregate.neuron(intNeuron).dFoF ses.neuron(intNeuron).dFoF];
			sSesAggregate.neuron(intNeuron).vecSpikes = [sSesAggregate.neuron(intNeuron).vecSpikes ses.neuron(intNeuron).vecSpikes];
		end
		
		%aggregate eye-tracking
		if isfield(ses,'sEyeTracking')
			%aggregate sEyeTracking
			cellFields = fieldnames(sSesAggregate.sEyeTracking);
			for intField=1:length(cellFields)
				strField = cellFields{intField};
				if strcmp(strField(1:3),'vec') || strcmp(strField(1:3),'ind')
					sSesAggregate.sEyeTracking.(strField) = [sSesAggregate.sEyeTracking.(strField) ses.sEyeTracking.(strField)];
				end
			end
		end
	end
end


