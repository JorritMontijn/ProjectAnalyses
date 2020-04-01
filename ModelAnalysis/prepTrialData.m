function [matModelRespNormTS3,matModelRespTS3] = prepTrialData(matModelRespTS3,vecTrialStimType)
	%UNTITLED Summary of this function goes here
	%   Detailed explanation goes here
	
	
	%% transform model resp
	%msg
	fprintf('Transforming data; removing residuals... [%s]\n',getTime);
	
	%data
	matModelRespTS3 = int8(matModelRespTS3); %[neurons * bins * trials]; data
	matModelRespNormTS3 = zeros(size(matModelRespTS3),'int8');
	
	%build constants
	vecUniqueStimTypes = unique(vecTrialStimType);
	intStimTypes = numel(vecUniqueStimTypes);
	[intNeurons,intBinsPerTrial,intTrials] = size(matModelRespTS3);
	for intStimTypeIdx=1:intStimTypes
		fprintf('  Processing stim type %d/%d... [%s]\n',intStimTypeIdx,intStimTypes,getTime);
		intStimType = vecUniqueStimTypes(intStimTypeIdx);
		indTrials = vecTrialStimType==intStimType;
		for intNeuron = 1:intNeurons
			for intBin=1:intBinsPerTrial
				matModelRespNormTS3(intNeuron,intBin,indTrials) = matModelRespTS3(intNeuron,intBin,indTrials) - mean(matModelRespTS3(intNeuron,intBin,indTrials));
			end
		end
	end
end

