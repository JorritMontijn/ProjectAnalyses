function matMeanRateRand = getTRandData(matMeanRateIn,vecStimIdx,strType)
	%getTRandData Randomize trial data
	%   matMeanRateRand = getTRandData(matMeanRateIn,vecStimIdx,strType)
	
	%get props
	[vecStimIdx,vecStimTypes] = val2idx(vecStimIdx);
	intStimNum = numel(vecStimTypes);
	[intNumN,intNumT] = size(matMeanRateIn);
	
	%get population mean
	vecOldMean = mean(matMeanRateIn,1);
	vecPopMeanFactor = vecOldMean ./ mean(vecOldMean);
	vecOldSd = std(matMeanRateIn,[],1);
	vecPopSdFactor = vecOldSd ./ mean(vecOldSd);
	
	%run
	matMeanRateRand = matMeanRateIn;
	if strcmp(strType,'Real')
		%do nothing
	else
		for intN=1:intNumN
			for intStim=1:intStimNum
				vecUseT = find(vecStimIdx==intStim);
				if strcmp(strType,'TShuff')
					%shuffle spikes
					matMeanRateRand(intN,vecUseT) = matMeanRateIn(intN,vecUseT(randperm(numel(vecUseT))));
				elseif strcmp(strType,'TPoiss')
					%generate spikes
					dblMean = mean(matMeanRateIn(intN,vecUseT));
					vecRates = poissrnd(dblMean,size(vecUseT));
					matMeanRateRand(intN,vecUseT) = vecRates;
				elseif strcmp(strType,'TUniStretch')
					%shuffle spikes & compensate for pop-rate change later
					matMeanRateRand(intN,vecUseT) = matMeanRateIn(intN,vecUseT(randperm(numel(vecUseT))));
				else
					%error
				end
			end
			
			%change for each neuron
			if strcmp(strType,'TSaturating')
				%smoothly saturating poisson
				vecR = matMeanRateRand(intN,:);
				
				%logistic slope is 0.5; both k and L increase slope
				vecOldHz = vecR;
				dblSatStart = mean(vecR)/2;
				indSat = vecR > dblSatStart;
				
				fLogistic = @(x,x0,L) x0 + L*2* (-0.5+1./(1+exp(-(2/L)*(x-x0))));
				x = vecR(indSat);
				x0 = dblSatStart;
				L = dblSatStart + 2*sqrt(dblSatStart);
				vecNewSat = fLogistic(x,x0,L);
				
				vecR(indSat) = vecNewSat;
				matMeanRateRand(intN,:) = vecR;
				
				%scatter(vecOldHz,vecR);
			end
		end
		if strcmp(strType,'TUniStretch')
			%unistretch
			vecNewMean = mean(matMeanRateRand,1);
			vecCompensateBy = vecOldMean./vecNewMean;
			matMeanRateRand = bsxfun(@times,matMeanRateRand,vecCompensateBy);
		elseif strcmp(strType,'TSdFixed')
			%fixed sd, scaling tuning
			
			%remove mean
			matMeanRateRand = bsxfun(@minus,matMeanRateRand,vecOldMean);
			
			%make sd uniform
			matMeanRateRand = bsxfun(@rdivide,matMeanRateRand,vecOldSd);
			
			%add mean back in
			matMeanRateRand = bsxfun(@plus,matMeanRateRand,vecOldMean);
		elseif strcmp(strType,'TSdScaling')
			%scaling sd as pop mean
			
			%remove mean
			matMeanRateRand = bsxfun(@minus,matMeanRateRand,vecOldMean);
			
			%make sd uniform
			matMeanRateRand = bsxfun(@rdivide,matMeanRateRand,vecOldSd);
			
			%multiply mean
			matMeanRateRand = bsxfun(@times,matMeanRateRand,vecOldMean); %basically recapitulates real data
			
			%add mean back in
			matMeanRateRand = bsxfun(@plus,matMeanRateRand,vecOldMean);
		elseif strcmp(strType,'TSdLinear')
			%linear scaling sd
			
			%remove mean
			matMeanRateRand = bsxfun(@minus,matMeanRateRand,vecOldMean);
			
			%make sd uniform
			matMeanRateRand = bsxfun(@rdivide,matMeanRateRand,vecOldSd);
			
			%multiply sd
			matMeanRateRand = bsxfun(@times,matMeanRateRand,vecPopMeanFactor);
			
			%add mean back in
			matMeanRateRand = bsxfun(@plus,matMeanRateRand,vecOldMean);
		elseif strcmp(strType,'TSdQuad')
			%linear scaling variance (sd quadratic)
			
			%remove mean
			matMeanRateRand = bsxfun(@minus,matMeanRateRand,vecOldMean);
			
			%make sd uniform
			matMeanRateRand = bsxfun(@rdivide,matMeanRateRand,vecOldSd);
			
			%multiply sd
			%matMeanRate = bsxfun(@times,matMeanRate,vecOldMean); %basically recapitulates real data
			matMeanRateRand = bsxfun(@times,matMeanRateRand,vecPopMeanFactor.^2);
			
			%add mean back in
			matMeanRateRand = bsxfun(@plus,matMeanRateRand,vecOldMean);
		elseif strcmp(strType,'TSdCube')
			%sd cubic
			
			%remove mean
			matMeanRateRand = bsxfun(@minus,matMeanRateRand,vecOldMean);
			
			%make sd uniform
			matMeanRateRand = bsxfun(@rdivide,matMeanRateRand,vecOldSd);
			
			%multiply sd
			%matMeanRate = bsxfun(@times,matMeanRate,vecOldMean); %basically recapitulates real data
			matMeanRateRand = bsxfun(@times,matMeanRateRand,vecPopMeanFactor.^3);
			
			%add mean back in
			matMeanRateRand = bsxfun(@plus,matMeanRateRand,vecOldMean);
			
		end
	end
	%ensure positive rates
	matMeanRateRand(matMeanRateRand<0.01)=0.01;
	
end

