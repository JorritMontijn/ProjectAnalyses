function [dblPerformance,vecDecodedIndexCV,matMahalDistsCV,dblMeanErrorDegs,matConfusion] = doCrossValidatedDecodingMD(matData,vecTrialTypes,intTypeCV)
	%UNTITLED9 Summary of this function goes here
	%   Detailed explanation goes here
	
	%check which kind of cross-validation
	if nargin < 3 || isempty(intTypeCV) || ~(intTypeCV < 3)
		intTypeCV = 1;
	end
	intVerbose = 0;
	
	%get number of trials
	vecTrialTypes = vecTrialTypes(:);
	intTrials = numel(vecTrialTypes);
	
	%check if matData is [trial x neuron] or [neuron x trial]
	if size(matData,1) == intTrials && size(matData,2) == intTrials
		%number of neurons and trials is the same
		warning([mfilename ':SameNeuronsTrials'],'Number of neurons and trials is identical; please double check the proper orientation of [intNeurons x intTrials]');
	elseif size(matData,1) == intTrials
		%rotate
		matData = matData';
	elseif size(matData,2) == intTrials
		%size is correct
	else
		error([mfilename ':SameNeuronsTrials'],'Size of matData and vecTrialTypes do not match');
	end
	intNeurons = size(matData,1);
	vecUniqueTrialTypes = unique(vecTrialTypes);
	intStimTypes = length(vecUniqueTrialTypes);
	vecTrialTypeIdx = label2idx(vecTrialTypes);
	
	%pre-allocate output
	matMahalDistsCV = zeros(intTrials,intStimTypes);
	
	%cross-validate
	if intTypeCV == 0
		%no CV
		for intStimType=1:intStimTypes
			vecMu = xmean(matData(:,vecTrialTypeIdx==intStimType),2);
			matCovar = cov(matData(:,vecTrialTypeIdx==intStimType)');
			matCovarInv = inv(matCovar);
			vecMahal=getMahal(matData,vecMu,matCovarInv);
			matMahalDistsCV(:,intStimType) = vecMahal;
		end
	elseif intTypeCV == 1
		%get distances
		for intStimType=1:intStimTypes
			if mod(intStimType,10) == 0 && intVerbose > 0,fprintf('Preparing stimulus %d/%d [%s]\n',intStimType,intStimTypes,getTime);pause(eps);end
			vecMu = xmean(matData(:,vecTrialTypeIdx==intStimType),2);
			matCovar = cov(matData(:,vecTrialTypeIdx==intStimType)');
			matCovarInv = inv(matCovar);
			vecMahal=getMahal(matData,vecMu,matCovarInv);
			matMahalDistsCV(:,intStimType) = vecMahal;
		end
		
		%leave one out
		for intLeaveOut=1:intTrials
			%get info on to-be-left-out trial
			indSelect = ~isnan(vecTrialTypeIdx);
			indSelect(intLeaveOut) = false;
			intTypeCVTrial = vecTrialTypeIdx(intLeaveOut);
			
			%calc CV mean
			vecMuCV = xmean(matData(:,(vecTrialTypeIdx==intTypeCVTrial)&indSelect),2);
			
			%calc CV covar
			matCovarInvCV = inv(cov(matData(:,(vecTrialTypeIdx==intTypeCVTrial)&indSelect)'));
			
			%get CV mahal dists
			matMahalDistsCV(intLeaveOut,intTypeCVTrial) = getMahal(matData(:,intLeaveOut),vecMuCV,matCovarInvCV);
			
			%msg
			if mod(intLeaveOut,100) == 0,pause(eps);end
			if mod(intLeaveOut,1000) == 0 && intVerbose > 0,fprintf('Decoding; now at trial %d/%d [%s]\n',intLeaveOut,intTrials,getTime);end
		end
	elseif intTypeCV == 2
		%remove repetition
		intRepNum = intTrials/intStimTypes;
		if round(intRepNum) ~= intRepNum,error([mfilename ':IncompleteRepetitions'],'Number of repetitions is not an integer');end
		intTrial = 0;
		for intRep=1:intRepNum
			
			%remove trials
			intTrialStopRep = intRep*intStimTypes;
			intTrialStartRep = intTrialStopRep - intStimTypes + 1;
			indSelect = true(1,intTrials);
			indSelect(intTrialStartRep:intTrialStopRep) = false;
			matThisData = matData(:,indSelect);
			vecThisTrialType = vecTrialTypeIdx(indSelect);
			
			%recalculate covariance matrix
			for intStimType=1:intStimTypes
				%get distances
				vecMu = xmean(matThisData(:,vecThisTrialType==intStimType),2);
				matCovarInv = inv(cov(matThisData(:,vecThisTrialType==intStimType)'));
				vecMahal=getMahal(matData(:,intTrialStartRep:intTrialStopRep),vecMu,matCovarInv);
				matMahalDistsCV(intTrialStartRep:intTrialStopRep,intStimType) = vecMahal;
				
				%msg
				intTrial = intTrial + 1;
				if mod(intTrial,100) == 0,pause(eps);end
				if mod(intTrial,1000) == 0 && intVerbose > 0,fprintf('Decoding; now at trial %d/%d [%s]\n',intTrial,intTrials,getTime);end
			end
		end
	end
	
	%output
	[dummy,vecDecodedIndexCV]=min(matMahalDistsCV,[],2);
	dblPerformance=sum(vecDecodedIndexCV==vecTrialTypeIdx)/intTrials;
	
	%error
	if nargout > 4
		vecDecodedValuesCV = vecUniqueTrialTypes(vecDecodedIndexCV);
		dblMeanErrorRads = mean(abs(circ_dist(vecDecodedValuesCV,vecTrialTypes)));
		dblMeanErrorDegs = rad2ang(dblMeanErrorRads);
	end
	
	
	%confusion matrix;
	if nargout > 5
		matConfusion = getFillGrid(zeros(intStimTypes),vecDecodedIndexCV,vecTrialTypeIdx,ones(intTrials,1))/intReps;
		%imagesc(matConfusion,[0 1]);colormap(hot);axis xy;colorbar;
	end
end

