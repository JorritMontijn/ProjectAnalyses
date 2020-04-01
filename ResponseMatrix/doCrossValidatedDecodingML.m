function [dblPerformanceCV,vecDecodedIndexCV,matPosteriorProbabilityCV] = doCrossValidatedDecodingML(matData,vecTrialTypes,intTypeCV)
	%UNTITLED9 Summary of this function goes here
	%   Detailed explanation goes here
	
	%check which kind of cross-validation
	if nargin < 3 || isempty(intTypeCV) || ~(intTypeCV < 3)
		intTypeCV = 1;
	end
	
	%get info
	vecTrialTypes = vecTrialTypes(:);
	intTrials = size(matData,1);
	intNeurons = size(matData,2);
	vecUniqueTrialTypes = unique(vecTrialTypes);
	intStimTypes = length(vecUniqueTrialTypes);
	
	%pre-allocate output
	matLikelihood = nan(intNeurons,intStimTypes,2); %mean, sd
	matPosteriorProbabilityCV = zeros(intTrials,intStimTypes,intNeurons);
	
	%cross-validate
	if intTypeCV == 1
		%get distances
		for intStimType=1:intStimTypes
			intStimTypeNumber = vecUniqueTrialTypes(intStimType);
			if mod(intStimType,10) == 0,fprintf('Preparing stimulus %d/%d [%s]\n',intStimType,intStimTypes,getTime);pause(eps);end
			vecMu = mean(matData(vecTrialTypes==intStimTypeNumber,:),1); %mean response this trial type per neuron
			vecSD = std(matData(vecTrialTypes==intStimTypeNumber,:),[],1); %sd response this trial type per neuron
			
			%put data in likelihood parameter matrix
			matLikelihood(:,intStimType,1) = vecMu;
			matLikelihood(:,intStimType,2) = vecSD;
			
			%calculate non-CV posterior
			for intTrial = 1:intTrials
				vecThisActivity = matData(intTrial,:);
				for intNeuron = 1:intNeurons
					%get mu and sigma
					dblMu = matLikelihood(intNeuron,intStimType,1);
					dblSigma = matLikelihood(intNeuron,intStimType,2);
					
					%calc probability
					if dblMu == 0 || dblSigma == 0
						P_ori_given_dFoF = nan;
					else
						P_ori_given_dFoF = normpdf(vecThisActivity(intNeuron),dblMu,dblSigma);
					end
					
					%put in matrix
					matPosteriorProbabilityCV(intTrial,intStimType,intNeuron) = P_ori_given_dFoF;
				end
			end
		end
		
		%leave one out
		for intLeaveOut=1:intTrials
			%get info on to-be-left-out trial
			indSelect = ~isnan(vecTrialTypes);
			indSelect(intLeaveOut) = false;
			intTypeCVTrial = vecTrialTypes(intLeaveOut);
			intTypeCVTrialNumber = find(vecUniqueTrialTypes==intTypeCVTrial);
			
			%calc CV mean + sd
			vecMuCV = mean(matData((vecTrialTypes==intTypeCVTrial)&indSelect,:),1);
			vecSDCV = std(matData((vecTrialTypes==intTypeCVTrial)&indSelect,:),[],1);
			
			%calc CV parameters
			matLikelihoodCV = matLikelihood;
			matLikelihoodCV(:,intTypeCVTrialNumber,1) = vecMuCV;
			matLikelihoodCV(:,intTypeCVTrialNumber,2) = vecSDCV;
			
			%get posterior probability
			vecThisActivity = matData(intLeaveOut,:);
			for intNeuron = 1:intNeurons
				%get mu and sigma
				dblMu = matLikelihoodCV(intNeuron,intTypeCVTrialNumber,1);
				dblSigma = matLikelihoodCV(intNeuron,intTypeCVTrialNumber,2);
				
				%calc probability
				if dblMu == 0 || dblSigma == 0
					P_ori_given_dFoF = nan;
				else
					P_ori_given_dFoF = normpdf(vecThisActivity(intNeuron),dblMu,dblSigma);
				end
				
				%put in matrix
				matPosteriorProbabilityCV(intLeaveOut,intTypeCVTrialNumber,intNeuron) = P_ori_given_dFoF;
			end
			
			%msg
			if mod(intLeaveOut,100) == 0,pause(eps);end
			if mod(intLeaveOut,1000) == 0,fprintf('Decoding; now at trial %d/%d [%s]\n',intLeaveOut,intTrials,getTime);end
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
			matThisData = matData(indSelect,:);
			vecThisTrialType = vecTrialTypes(indSelect);
			
			%recalculate covariance matrix
			for intStimType=1:intStimTypes
				%get distances
				vecMu = mean(matThisData(vecThisTrialType==intStimType,:),1);
				matCovar = cov(matThisData(vecThisTrialType==intStimType,:));
				matCovarInv = inv(matCovar);
				vecMahal=getMahal(matData(intTrialStartRep:intTrialStopRep,:),vecMu,matCovarInv);
				matPosteriorProbabilityCV(intTrialStartRep:intTrialStopRep,intStimType) = vecMahal;
				
				%msg
				intTrial = intTrial + 1;
				if mod(intTrial,100) == 0,pause(eps);end
				if mod(intTrial,1000) == 0,fprintf('Decoding; now at trial %d/%d [%s]\n',intTrial,intTrials,getTime);end
			end
		end
	end
	
	%calculate output
	[dummy,vecDecodedIndexCV]=max(prod(matPosteriorProbabilityCV,3),[],2);
	vecDecodedIndexCV = vecUniqueTrialTypes(vecDecodedIndexCV);
	dblPerformanceCV = sum(vecDecodedIndexCV == vecTrialTypes)/length(vecDecodedIndexCV);
end

