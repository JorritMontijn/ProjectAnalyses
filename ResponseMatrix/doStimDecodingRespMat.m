function sOut = doStimDecodingRespMat(matResp,structStim,sParams)
	%doStimDecodingRespMat Performs stimulus decoding of arbitrary stimulus
	%						types using getStimulusTypes.m
	%   Syntax: sOut = doStimDecodingRespMat(matResp,structStim,sParams)
	%	Input:
	%	- matResp; response matrix per trial (spikes/mean dF/F)
	%	- structStim; structure containing trial information
	%	- sParams; structure containing additional inputs. Fields are:
	%		- sTypes
	%		- cellSelect
	%		- vecIncludeCells
	%		- verbose
	%	Output:
	%	- sOut; output structure, containing fields:
	%		- cellVecPost; cell for each stimulus containing vector of
	%			stimulus posterior probabilities
	%		- vecDecodedStimType; vector containing decoded stimulus type
	%			for all stimuli (ML)
	%		- vecStimType; vector containing actual stimulus types
	%
	%Dependencies:
	% - getStimulusTypes.m
	% - getSelectionVectors.m
	% - getNeuronResponse.m
	% - buildStimDecodingLikelihood.m
	% - buildStimDecodingPosterior.m
	%
	%	Version history:
	%	1.0 - November 19 2014
	%	Created by Jorrit Montijn, based on doStimDecoding.m [by JM]
	
	%pre-alloc
	sOut = struct;
	if nargin == 1, sParams = struct;end
	
	% get stim list
	if ~isfield(sParams,'sTypes'),sTypes = getStimulusTypes(structStim);else sTypes = sParams.sTypes;end
	
	% get indexing vectors for unique stimulus combinations
	if ~isfield(sParams,'cellSelect'),cellSelect = getSelectionVectors(structStim,sTypes);else cellSelect = sParams.cellSelect;end
	%get other inputs
	if ~isfield(sParams,'vecIncludeCells'),vecNeurons = 1:size(matResp,1);else vecNeurons = sParams.vecIncludeCells;end
	if ~isfield(sParams,'vecLikelihoodTrials'),vecLikelihoodTrials = 1:(length(structStim.FrameOn)/2);else vecLikelihoodTrials = sParams.vecLikelihoodTrials;end
	if ~isfield(sParams,'verbose'),verbose = 1;else verbose = sParams.verbose;end
	if ~isfield(sParams,'intPreBaselineRemoval'),intPreBaselineRemoval = [];else intPreBaselineRemoval = sParams.intPreBaselineRemoval;end
	
	%build likelihood
	structIn.sTypes = sTypes;
	structIn.cellSelect = cellSelect;
	structIn.vecLikelihoodTrials = vecLikelihoodTrials;
	structIn.intPreBaselineRemoval = intPreBaselineRemoval;
	matLikelihood = buildStimDecodingLikelihoodRespMat(matResp,structIn);
	[row1,col1]=find(matLikelihood(:,:,1)==0);
	[row2,col2]=find(matLikelihood(:,:,2)==0);
	vecRemovedNeurons = unique([row1;row2]);
	vecNeurons = vecNeurons(~ismember(vecNeurons,vecRemovedNeurons));
	
	%get additional data
	intStims = length(cellSelect{1});
	vecStimTypes = nan(1,intStims);
	for intType=1:numel(cellSelect)
		vecStimTypes(cellSelect{intType}) = intType;
	end
	vecCorrect = nan(1,intStims);
	vecStims = 1:intStims;
	if ~isempty(vecLikelihoodTrials)
		if islogical(vecLikelihoodTrials)
			vecStims = vecStims(~vecLikelihoodTrials);
		elseif length(vecLikelihoodTrials) ~= length(vecStims)
			vecStims = vecStims(~ismember(vecStims,vecLikelihoodTrials));
		end
	end
	if isfield(sParams,'vecDecodeTrials'),vecStims = sParams.vecDecodeTrials;end
	if islogical(vecStims)
		vecDummy = 1:intStims;
		vecStims = vecDummy(vecStims);
	end
	
	%pre-allocate output
	cellVecPost = cell(1,intStims);
	vecDecodedStimType = nan(1,intStims);
	vecStimType = nan(1,intStims);
	
	%loop through stimuli
	for intStim=vecStims
		%get pop resp for this stimulus
		matRespTemp = matResp(:,intStim);
		if isempty(matRespTemp),continue;end
		
		%posterior vars
		sInPost.matLikelihood = matLikelihood;
		sInPost.vecIncludeCells = vecNeurons;
		
		%get posterior
		structPosterior = buildStimDecodingPosteriorRespMat(matRespTemp,sInPost);
		matPost = structPosterior.matPost;
		vecMeanPost = structPosterior.vecMeanPost;
		
		
		%% plot
		[dummy,intDecodedStimType] = max(vecMeanPost);
		
		%% output
		intStimType = vecStimTypes(intStim);
		if intStimType == intDecodedStimType
			strErr = '';
			vecCorrect(intStim) = true;
		else
			strErr = '; Incorrect';
			vecCorrect(intStim) = false;
		end
		cellVecPost{intStim} = vecMeanPost;
		vecDecodedStimType(intStim) = intDecodedStimType;
		vecStimType(intStim) = intStimType;
		
		%% message
		if verbose == 1,fprintf('Decoded stimulus %d of %d; decoded type=%d, actual type=%d%s\n',intStim,intStims,intDecodedStimType,intStimType,strErr);end
	end
	
	%put in output
	dblPercCorr = (nansum(vecCorrect)/length(vecStims))*100;
	sOut.dblPercCorr = dblPercCorr;
	sOut.cellVecPost = cellVecPost;
	sOut.vecDecodedStimType = vecDecodedStimType;
	sOut.vecStimType = vecStimType;
	sOut.vecRemovedNeurons = vecRemovedNeurons;
	if verbose ~= 0,fprintf('\nDecoding performance: %d of %d [%.0f%%] correct\n',nansum(vecCorrect),length(vecStims),dblPercCorr);end
end

