function matLikelihood = buildStimDecodingLikelihoodRespMat(matResp,structIn)
	%buildStimDecodingLikelihood Build Bayesian likelihood for supplied session
	%	Syntax: structOut = buildDecodingLikelihood(structIn)
	%   Input: structIn, a structure containing the following fields:
	%	- ses, session data
	%	- sTypes, output of getStimulusTypes()
	%	- cellSelect, output of getSelectionVectors()
	%	Optional fields are:
	%	- vecIncludeCells, vector of cells to include in analysis; [Default: all]
	%	- vecIgnoreTrial, vector trials to ignore; [Default: none]
	%	Output: 
	%	- matLikelihood, 3D matrix containing distribution parameters of likelihood function:
	%		matLikelihood(intNeuron,intStimType,intParameter)
	%		Parameter data for Gaussian is:
	%		- intParameter == 1; mean response
	%		- intParameter == 2; std response
	%
	%Dependencies:
	% - getNeuronResponse.m
	%
	%	Version history:
	%	1.0 - July 22 2013
	%	Created by Jorrit Montijn
	%	2.0 - May 19 2014
	%	Fixed likelihood & added pre-baseline removal [by JM]
	
	%% define input variables
	sTypes = structIn.sTypes;
	matTypes = sTypes.matTypes;
	cellSelect = structIn.cellSelect;
	intNeurons = size(matResp,1);
	intTrials = size(matResp,2);
	if ~isfield(structIn,'vecIncludeCells') || isempty(structIn.vecIncludeCells), vecIncludeCells = 1:intNeurons;else vecIncludeCells = structIn.vecIncludeCells;end
	if ~isfield(structIn,'vecLikelihoodTrials') || isempty(structIn.vecLikelihoodTrials), vecLikelihoodTrials = [];else vecLikelihoodTrials = structIn.vecLikelihoodTrials;end
	if ~isfield(structIn,'intPreBaselineRemoval') || isempty(structIn.intPreBaselineRemoval), intPreBaselineRemoval = [];else intPreBaselineRemoval = structIn.intPreBaselineRemoval;end
	
	
	%% pre-allocate variables
	%approximate every stimulus response separately with normal distribution
	intParams = 2;
	intTypes = size(matTypes,2);
	matLikelihood = nan(intNeurons,intTypes,intParams);
	
	%% build likelihood lookup table
	% pre-allocate
	if isempty(vecLikelihoodTrials)
		vecLikelihoodTrials = true(1,intTrials);
	elseif ~islogical(vecLikelihoodTrials)
		vecTemp = false(1,intTrials);
		vecTemp(vecLikelihoodTrials) = true;
		vecLikelihoodTrials = vecTemp;
	elseif length(vecLikelihoodTrials) ~= intTrials
		error([mfilename ':IncorrectLength'],'Length of vecLikelihoodTrials [%d] is incompatible with number of stimuli [%d]',length(vecLikelihoodTrials),intTrials);
	end
	vecStimTypeIndex = 1:intTypes;
	
	% loop through stim types
	for intStimType=vecStimTypeIndex
		vecStims = cellSelect{intStimType} & vecLikelihoodTrials;
		matRespTemp = matResp(vecIncludeCells,vecStims);
		matLikelihood(vecIncludeCells,intStimType,1) = mean(matRespTemp,2);
		matLikelihood(vecIncludeCells,intStimType,2) = std(matRespTemp,0,2);
	end
end

