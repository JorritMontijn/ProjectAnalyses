function sSesAggregate = buildMultiSesAggregate(ses,sSesAggregate)
	%UNTITLED2 Summary of this function goes here
	%   Detailed explanation goes here
	%%
	if ~exist('sSesAggregate','var') || isempty(sSesAggregate)
		sSesAggregate = ses;
	else
		%get starting values from previous session
		intLastFrame = length(sSesAggregate.neuron(1).dFoF);
		if isfield(sSesAggregate.structStim,'SecsOff')
			dblLastSec = (sSesAggregate.structStim.SecsOff(end)/sSesAggregate.structStim.FrameOff(end))*length(sSesAggregate.neuron(1).dFoF);
		end
		
		%aggregate structStim
		cellFields = fieldnames(ses.structStim);
		intNumTrials = length(ses.structStim.FrameOn);
		for intField=1:length(cellFields)
			strField = cellFields{intField};
			if ~isfield(sSesAggregate.structStim,strField),continue;end
			if ~isempty(strfind(strField,'Secs')) %#ok<*STREMP>
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
		
		%aggregate object timeseries
		cellFields=fieldnames(ses);
		cellTimeseriesFields = {'F','npF','dFoF','vecSpikes'};
		for intField=1:numel(cellFields)
			strObjectField = cellFields{intField};
			for intTimeseriesField=1:numel(cellTimeseriesFields)
				strFieldTS = cellTimeseriesFields{intTimeseriesField};
				
				if isfield(ses.(strObjectField),strFieldTS)
					for intNeuron=1:numel(ses.(strObjectField))
						sSesAggregate.(strObjectField)(intNeuron).(strFieldTS) = [sSesAggregate.(strObjectField)(intNeuron).(strFieldTS) ses.(strObjectField)(intNeuron).(strFieldTS)];
					end
				end
			end
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


