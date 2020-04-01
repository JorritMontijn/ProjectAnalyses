function structOut = calcStimCorrsRespMat(matResp,cellSelect)
	%calcStimCorrs Calculates signal+noise correlations from stimulus
	%presentations
	%syntax: structOut = calcStimCorrs(ses)
	%	input:
	%	- ses: structure containing session information (prepro output)
	%
	%	output:
	%	- structOut: a structure containing the following fields:
	%		- matSignalCorrs: matrix containing signal correlations
	%		- matNoiseCorrs: matrix containing noise correlations
	%		- matStimResponse: matrix containing mean response of neuron by stimulus presentation
	%		- matSignalResponse: matrix containing mean response of neuron by stimulus type
	%		- sTypes: structure containing data on stimulus types of the
	%			different separated classes
	%
	%Dependencies:
	% - getStimulusTypes.m
	% - getSelectionVectors.m
	% - getNeuronResponse.m
	%
	%	Version history:
	%	1.0 - May 31 2013
	%	Created by Jorrit Montijn
	%	2.0 - July 25 2013
	%	Modified to work with arbitrary stimulus categories based on
	%	getStimulusTypes/getSelectionVectors/getNeuronResponse triple combo
	%	3.0 - May 26 2015
	%	Modified from calcStimCorrs() to accept response matrices
	
	%% pre-allocate response matrix
	intTypes = length(cellSelect);
	intReps = sum(cellSelect{1});
	intNeurons = size(matResp,1);
	intTrials = size(matResp,2);
	matStimResponse = nan(intTypes,intReps,intNeurons);
	
	%% retrieve responses
	for intType=1:intTypes
		vecSelect = cellSelect{intType};
		if sum(vecSelect) > 0
			matStimResponse(intType,1:sum(vecSelect),:) = shiftdim(matResp(:,vecSelect),1);
		end
	end
	
	%% calculate signal correlations
	matSignalResponse = squeeze(mean(matStimResponse,2));%nan(intStims,intNeurons);
	
	%get mean per stim type
	matSignalCorrs = corr(matSignalResponse);
	
	%set autocorrelation at 0
	matSignalCorrs(diag(true(1,intNeurons))) = 0;
	
	
	%% do same for noise correlations
	matAllNoise = nan(intNeurons,intNeurons,intTypes);
	for intType=1:intTypes
		%calculate noise correlations per stimulus type
		matNoiseResponse = squeeze(matStimResponse(intType,:,:));
		matAllNoise(:,:,intType) = corr(matNoiseResponse);
	end
	%then average
	matNoiseCorrs = mean(matAllNoise,3);
	
	%set autocorrelation to 0
	matNoiseCorrs(diag(true(1,intNeurons))) = 0;
	
	%% put in output
	structOut.matAllNoise = matAllNoise;
	structOut.matSignalCorrs = matSignalCorrs;
	structOut.matNoiseCorrs = matNoiseCorrs;
	structOut.matStimResponse = matStimResponse;
	structOut.matSignalResponse = matSignalResponse;
	structOut.cellSelect = cellSelect;
end
