function structPosterior = buildStimDecodingPosteriorRespMat(matResp,sInPost)
	%buildDecodingPosterior Build Bayesian likelihood for supplied session
	%	Syntax: structPosterior = buildDecodingPosterior(sInPost)
	%   Input: sInPost, a structure containing the following fields:
	%	- structLikelihood, structure containing likelihood data
	%	- data
	%	- vecIncludeCells, vector of cells to include in analysis; [Default: all]
	%	Output:
	
	
	%% define input variables
	intNeurons = size(matResp,1);
	if ~isfield(sInPost,'vecIncludeCells') || isempty(sInPost.vecIncludeCells), vecIncludeCells = 1:intNeurons;else vecIncludeCells = sInPost.vecIncludeCells;end
	matLikelihood = sInPost.matLikelihood;
	intStimTypes = size(matLikelihood,2);
	matPost = ones(intNeurons,intStimTypes);
	
	%% build posterior probability
	for intNeuron=vecIncludeCells
		thisActivity = matResp(intNeuron);
		
		for intStimType=1:intStimTypes
			%get mu and sigma
			mu = matLikelihood(intNeuron,intStimType,1);
			sigma = matLikelihood(intNeuron,intStimType,2);
			
			%calc probability
			if mu == 0 || sigma == 0
				P_ori_given_dFoF = nan;
			else
				P_ori_given_dFoF = normpdf(thisActivity,mu,sigma);
			end
			
			%put in matrix
			matPost(intNeuron,intStimType) = P_ori_given_dFoF;
		end
	end
	
	%% put into output structure
	structPosterior = struct;
	structPosterior.matPost = matPost;
	structPosterior.vecMeanPost = nanprod(matPost,1);
	structPosterior.vecMeanPostM = nanmean(matPost,1);
end

