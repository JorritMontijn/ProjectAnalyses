function [vecSpikeTimes,dblPrefOri] = getCorrelatedSpikeTrains(vecMu,matSigma,dblBinSizeSecs,dblTotDur)
	%getGeneratedData Generates neural data using von Mises tuning curves
	%    matResp = getGeneratedSpikingData(intN,intRep,dblKappa,dblDistDprime)
	%
	%Uses https://github.com/mackelab/CorBinian
	%
	%Version History:
	%2019-09-23 Created getGeneratedSpikingData function [by Jorrit Montijn]
	
	%% generate vecR & matSigma
	dblBinSizeSecs = 1e-3;
	dblTotDur=100;
	intN = 31;
	vecMu = 2*(abs(randn(1,intN))+3);
	matSigma1 = randn(intN)/5;
	matSigma2 = matSigma1';
	matSigma2(tril(true(intN))) = 0;
	matSigma1(triu(true(intN))) = 0;
	matSigma = matSigma1 + matSigma2 + diag(vecMu.^2);
	matSigma=matGenS
	
	%% using Macke et al. (2009), Neural Computation
	%get number of bins
	intTotBins = round(dblTotDur/dblBinSizeSecs);
	intDefBins = round(10/dblBinSizeSecs);
	intGenBins = intTotBins - intDefBins;
	
	%transform spiking rate into spiking probability per bin
	%vecR = (vecMu*intTotBins)*dblBinSizeSecs;
	vecR = vecMu*dblBinSizeSecs;
	dblCorrFactor = mean(vecR(:)) ./ mean(diag(matSigma));
	%matSigma = matSigma * dblCorrFactor;
	
	%% calculate complement gamma to r for dichotimized gaussian
	vecGamma = norminv(vecR);
	vecR_SanityCheck = normcdf(vecGamma,0,1); %must be identical to vecR
	if ~all((vecR_SanityCheck - vecR) < eps),error([mfilename ':ValueError'],'Unsupported values supplied; vecR is not identical to vecR_SanityCheck');end
	
	%% calculate complement Lambda to matSigma for dichotimized gaussian
	[samples1,gammas,Lambda,joints2D,pmfs,hists,supports]=DGPoisson(vecR,matSigma,intDefBins);
	
	%% generate additional samples
	%if intGenBins > 0
	%	[samples2,hists]=SampleDGAnyMarginal(gammas,Lambda,supports,intGenBins);
	%end
	%matSamples = cat(1,samples1',samples2);
	matSamples = samples1';
	
	%% turn into spike times
	%     gammas: the discretization thresholds, as described in the paper. When
	%       sampling. The k-th dimension of the output random variable is f if e.g.
	%       supports{k}(1)=f and gammas{k}(f) <= U(k) <= gammas{k}(f+1)
	vecGenR = (sum(matSamples,1)*dblBinSizeSecs)/dblTotDur;
	matGenS = cov(matSamples);
	
	figure
	subplot(2,2,1)
	scatter(vecR,vecGenR);
	title('R_req (x) vs R_gen (y)')
	colorbar
	
	subplot(2,2,2)
	imagesc(matSigma)
	title('Sigma_req')
	colorbar
	
	subplot(2,2,3)
	imagesc(matGenS)
	title('Sigma_gen')
	colorbar
	
	subplot(2,2,4)
	imagesc(matSigma - matGenS)
	title('Sigma_req - Sigma_gen')
	colorbar
end

