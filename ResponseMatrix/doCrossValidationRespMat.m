function sMetaOutSD = doCrossValidationRespMat(matResp,structStim,sParams)
	%doCrossValidationRespMat Performs automatic cross-validated ML decoding
	%   Syntax: sMetaOutSD = doCrossValidationRespMat(matResp,structStim,sParams)
	
	%get inputs
	if ~exist('sParams','var'),sParams=struct;end
	if ~isfield(sParams,'verbose'),sParams.verbose = 0;end
	if ~isfield(sParams,'cellSelect') || ~isfield(sParams,'sTypes')
		sTypes = getStimulusTypes(structStim,{'Orientation'});
		cellSelect = getSelectionVectors(structStim,sTypes);
	else
		sTypes = sParams.sTypes;
		cellSelect = sParams.cellSelect;
	end
	
	%get data
	verbose = sParams.verbose;
	sParams.verbose = 0;
	intNumTypes = length(cellSelect);
	intNumTrials = numel(cellSelect{1});
	intNumRepetitions = sum(cellSelect{1});
	vecTrials = true(1,intNumTrials);
	
	%pre-alloc output
	vecCorrect = false(1,intNumTrials);
	cellVecPost = cell(1,intNumTrials);
	vecDecodedStimType = nan(1,intNumTrials);
	vecStimType = nan(1,intNumTrials);
	
	for intRepetition = 1:intNumRepetitions
		%select subset of trials to build likelihood
		vecToBeDecodedTrials = nan(1,intNumTypes);
		for intType=1:intNumTypes
			vecSubSelect = find(cellSelect{intType},intRepetition);
			vecToBeDecodedTrials(intType) = vecSubSelect(end);
		end
		
		%make likelihood vector
		vecTheseLikelihoodTrials = vecTrials;
		vecTheseLikelihoodTrials(vecToBeDecodedTrials) = false;
		sParams.vecLikelihoodTrials = vecTheseLikelihoodTrials;
		sParams.cellSelect = cellSelect;
		sParams.sTypes = sTypes;
		
		%do decoding
		sOutSD = doStimDecodingRespMat(matResp,structStim,sParams);
		
		%put in output
		cellVecPost(vecToBeDecodedTrials) = sOutSD.cellVecPost(vecToBeDecodedTrials);
		vecDecodedStimType(vecToBeDecodedTrials) = sOutSD.vecDecodedStimType(vecToBeDecodedTrials);
		vecStimType(vecToBeDecodedTrials) = sOutSD.vecStimType(vecToBeDecodedTrials);
		
		
		%message
		if verbose
			for intTrial = vecToBeDecodedTrials
				intDecodedStimType = vecDecodedStimType(intTrial);
				intStimType = vecStimType(intTrial);
				if intStimType == intDecodedStimType
					strErr = '';
					vecCorrect(intTrial) = true;
				else
					strErr = '; Incorrect';
					vecCorrect(intTrial) = false;
				end
				fprintf('Decoded stimulus %d of %d; decoded type=%d, actual type=%d%s\n',intTrial,intNumTrials,vecDecodedStimType(intTrial),vecStimType(intTrial),strErr);
			end
		end
	end
	
	%save to structure
	sMetaOutSD.cellVecPost = cellVecPost;
	sMetaOutSD.vecDecodedStimType = vecDecodedStimType;
	sMetaOutSD.vecStimType = vecStimType;
	
	%message
	if verbose,fprintf('\nDecoding performance: %d of %d [%.0f%%] correct\n',nansum(vecCorrect),intNumTrials,(nansum(vecCorrect)/intNumTrials)*100);end
end

