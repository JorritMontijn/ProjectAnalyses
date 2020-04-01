function [dblPredA,matPredA,dblDprimeSquaredOffDiagonal,matDprimeSquared,matDprimeSquared_diagonal] = getSeparation(matData,vecTrialTypes,boolLinear)
	%getSeparation Syntax:
	%   [dblPredA,matPredA,dblDprimeSquaredOffDiagonal,matDprimeSquared,matDprimeSquared_diagonal] = getSeparation(matData,vecTrialTypes,boolLinear)
	
	%check which kind of cross-validation
	if nargin < 3 || isempty(boolLinear) || ~(boolLinear < 2)
		boolLinear = true;
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
	vecUniqueTrialTypes = unique(vecTrialTypes);
	intStimTypes = length(vecUniqueTrialTypes);
	vecTrialTypeIdx = label2idx(vecTrialTypes);
	
	%calculate D prime squared
	matDprimeSquared = nan(intStimTypes,intStimTypes);
	matDprimeSquared_diag = nan(intStimTypes,intStimTypes);
	if boolLinear
		%linear case; use only the diagonal
		for intStimType1=1:(intStimTypes-1)
			%get mean/std for stimulus type 1
			vecMu1 = xmean(matData(:,vecTrialTypeIdx==intStimType1),2);
			vecVar1 = xvar(matData(:,vecTrialTypeIdx==intStimType1),2);
			for intStimType2 = (intStimType1+1):intStimTypes
				%get mean/std for stimulus type 2
				vecMu2 = xmean(matData(:,vecTrialTypeIdx==intStimType2),2);
				vecVar2 = xvar(matData(:,vecTrialTypeIdx==intStimType2),2);
				
				%get separation of means, and pooled variance
				vecDeltaMu = vecMu1 - vecMu2;
				vecVarPooled = (vecVar1 + vecVar2)/2;
				
				%build diagonal covariance matrix
				matCovar = diag(vecVarPooled);
				matCovarInv = diag(1./vecVarPooled);
				
				%calculate d'^2
				dblDprimeSquared = vecDeltaMu'*matCovarInv*vecDeltaMu;
				dblDprimeSquared_diag = (vecDeltaMu'*matCovarInv*vecDeltaMu).^2 / (vecDeltaMu'*matCovarInv*matCovar*matCovarInv*vecDeltaMu);
				
				%assign to output
				matDprimeSquared(intStimType1,intStimType2) = dblDprimeSquared;
				matDprimeSquared(intStimType2,intStimType1) = dblDprimeSquared;
				
				%assign to output
				matDprimeSquared_diag(intStimType1,intStimType2) = dblDprimeSquared_diag;
				matDprimeSquared_diag(intStimType2,intStimType1) = dblDprimeSquared_diag;
			end
		end
	else
		%multivariate case
		for intStimType1=1:(intStimTypes-1)
			%get mean for stimulus 1
			vecMu1 = xmean(matData(:,vecTrialTypeIdx==intStimType1),2);
			
			%get covariance matrix for stimulus 1
			matCov1 = cov(matData(:,vecTrialTypeIdx==intStimType1)');
			%matCovInv1 = inv(matCov1);
			for intStimType2 = (intStimType1+1):intStimTypes
				%get mean for stimulus 2
				vecMu2 = mean(matData(:,vecTrialTypeIdx==intStimType2),2);
				
				%get covariance matrix for stimulus 2
				matCov2 = cov(matData(:,vecTrialTypeIdx==intStimType2)');
				%matCovInv2 = inv(matCov2);
				
				%get separation of means
				vecDeltaMu = vecMu1 - vecMu2;
				
				%get mean covariance matrix
				matCovar = (matCov1 + matCov2)/2;
				matCovarInv = inv(matCovar);
				
				%calculate d'^2
				dblDprimeSquared = vecDeltaMu'*matCovarInv*vecDeltaMu;
				dblDprimeSquared_diag = (vecDeltaMu'*matCovarInv*vecDeltaMu).^2 / (vecDeltaMu'*matCovarInv*matCovar*matCovarInv*vecDeltaMu);
				
				%assign to output
				matDprimeSquared(intStimType1,intStimType2) = dblDprimeSquared;
				matDprimeSquared(intStimType2,intStimType1) = dblDprimeSquared;
				
				%assign to output
				matDprimeSquared_diag(intStimType1,intStimType2) = dblDprimeSquared_diag;
				matDprimeSquared_diag(intStimType2,intStimType1) = dblDprimeSquared_diag;
			end
		end
	end
	
	%calculate predicted accuracy
	vecDprimeSquaredOffDiag = [diag(matDprimeSquared,-1); matDprimeSquared(1,end)];
	vecDprimeSquaredDiag_OffDiag = [diag(matDprimeSquared_diag,-1); matDprimeSquared_diag(1,end)];
	matPredA = 0.5+0.5*erf((sqrt(matDprimeSquared)/2)/sqrt(2));
	%dblPredA = mean((1/intStimTypes) + nansum((1/intStimTypes)* ((matPredA-0.5)*2)));
	%dblPredA = xmean((1/intStimTypes) + nansum((1/intStimTypes)* matPredA),2);
	dblDprimeSquaredOffDiagonal = xmean(vecDprimeSquaredOffDiag,1);
	
	matDprimeSquared_diagonal = xmean(vecDprimeSquaredDiag_OffDiag,1);
	
	dblPredA = nansum((1/intStimTypes)*erf((sqrt(vecDprimeSquaredOffDiag)/2)/sqrt(2)));
	dblPredA = dblPredA + (1-dblPredA)*(1/intStimTypes);
end

