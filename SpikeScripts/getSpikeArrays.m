function [cellSpikes,structStimSpikes] = getSpikeArrays(ses)
	%getSpikeArrays Returns time-stamped spikes and structStim from ses
	%   [cellSpikes,structStimSpikes] = getSpikeArrays(ses)
	
	%pre-allocate
	intNeurons = numel(ses.neuron);
	cellSpikes = cell(1,intNeurons);
	dblSamplingFreq = ses.structStim.FrameOff(end)/ses.structStim.SecsOff(end);
	
	%loop through neurons to get spikes
	for intNeuron=1:intNeurons
		vecSpikes = ses.neuron(intNeuron).vecSpikes;
		while max(vecSpikes) > 0
			vecSpikeIndices = find(vecSpikes>0);
			cellSpikes{intNeuron} = [cellSpikes{intNeuron} vecSpikeIndices/dblSamplingFreq];
			vecSpikes(vecSpikeIndices) = vecSpikes(vecSpikeIndices) - 1;
		end
		cellSpikes{intNeuron} = sort(cellSpikes{intNeuron},'ascend');
	end
	
	%transform structstim
	if nargout > 1
		structStimSpikes = ses.structStim;
		structStimSpikes.TimeOn = ses.structStim.FrameOn/dblSamplingFreq;
		structStimSpikes.TimeOff = ses.structStim.FrameOff/dblSamplingFreq;
		
		%check for behavior fields (anything with 'Pulses')
		cellFields = fieldnames(structStimSpikes)';
		vecMatch=cellfun(@numel,strfind(cellFields,'Pulses'));
		for intField=find(vecMatch)
			strField = cellFields{intField};
			strNewField = strrep(strField,'Pulses','Time');
			varOldData = structStimSpikes.(strField);
			if iscell(varOldData)
				structStimSpikes.(strNewField) = cellfun(@(x) x/dblSamplingFreq,varOldData,'UniformOutput', false);
			else
				structStimSpikes.(strNewField) = varOldData/dblSamplingFreq;
			end
		end
	end
end
