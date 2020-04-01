function matResp = getNeuronResponseRespMat(matRespHiRes,structStim,vecNeurons,vecStims,structParams)
	%getNeuronResponseRespMat Retrieves neuronal response for certain stimuli
	%	Syntax: matResp = getNeuronResponseRespMat(matRespHiRes,structStim,vecNeurons,vecStims,structParams)
	%   Input:
	%	- matRespHiRes, oversampled data matrix
	%	- structStim, structure containing stimulus information
	%	- vecNeurons, vector of which neurons to include
	%	- vecStims, vector of which stimuli to include [-1 returns response
	%		outside stimulus presentations]; works well with cellSelect{}
	%		output vector (output from getSelectionVectors)
	%	- structParams, structure with optional parameters
	%
	%	Output: 
	%	- matResp, 2D matrix containing neuronal response per stimulus per neuron:
	%		matResp(intNeuron,intStimPres)
	%
	%	Version history:
	%	1.0 - July 22 2013
	%	Created by Jorrit Montijn
	%	2.0 - May 19 2014
	%	Added support for preceding baseline subtraction [by JM]
	
	%check inputs
	[intNeurons,intSizeT] = size(matRespHiRes);
	intTrials = length(structStim.FrameOn);
	if nargin < 5
		structParams = struct;
	end
	if nargin < 4 || isempty(vecStims)
		vecStims = 1:intTrials;
	end
	if nargin < 3 || isempty(vecNeurons)
		vecNeurons = 1:intNeurons;
	end
	
	%select frames
	if vecStims == -1
		%baseline
		vecStartFrames = structStim.FrameOff;
		vecStopFrames = [structStim.FrameOn(2:end) intSizeT];
	else
		%stimuli
		vecStartFrames = structStim.FrameOn(vecStims);
		vecStopFrames = structStim.FrameOff(vecStims);
	end
	
	%check if frame subset selection is requested
	if isfield(structParams,'intStopOffset')
		vecStopFrames = vecStartFrames + structParams.intStopOffset;
	end
	if isfield(structParams,'intStartOffset')
		vecStartFrames = vecStartFrames + structParams.intStartOffset;
	end
	if isfield(structParams,'intPreBaselineRemoval')
		intPreBaselineRemoval = structParams.intPreBaselineRemoval;%dblBaselineSecs
	else
		intPreBaselineRemoval = [];
	end
	
	%retrieve data
	if islogical(vecStims)
		intRepetitions = sum(vecStims);
	else
		intRepetitions = numel(vecStims);
	end
	%check if vector or scalar
	boolVec = length(vecNeurons) > 1;
	if boolVec
		if size(vecNeurons,2) == 1,vecNeurons = vecNeurons';end
		matResp = nan(max(vecNeurons),intRepetitions);
	else matResp = nan(1,intRepetitions);
	end
	
	%go through stims
	for intStimPres=1:length(vecStartFrames)
		intStartFrame = vecStartFrames(intStimPres);
		intStopFrame = vecStopFrames(intStimPres);
		if boolVec
			for intNeuron=vecNeurons
				if ~isempty(intPreBaselineRemoval),dblBaseline = mean(matRespHiRes(intNeuron,(intStartFrame-intPreBaselineRemoval):(intStartFrame-1)));
				else dblBaseline = 0;end
				matResp(intNeuron,intStimPres) = mean(matRespHiRes(intNeuron,intStartFrame:intStopFrame))-dblBaseline;
			end; 
		else
			if ~isempty(intPreBaselineRemoval),dblBaseline = mean(matRespHiRes(vecNeurons,(intStartFrame-intPreBaselineRemoval):(intStartFrame-1)));
			else dblBaseline = 0;end
			matResp(1,intStimPres) = mean(matRespHiRes(vecNeurons,intStartFrame:intStopFrame)) - dblBaseline;
		end
	end
end

