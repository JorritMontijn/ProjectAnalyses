
dblBaseRate = 1;
vecRates = [5 0 5];
vecDurs = [1 0.5 1.5];
dblStimDur = sum(vecDurs);
intTrials = 10;
vecEventTimes = (3*dblStimDur):dblStimDur:(dblStimDur*(intTrials+2.9));
dblEndT = vecEventTimes(end)+4*dblStimDur;
vecSpikeTimes = getGeneratedMultiPhasicR(dblBaseRate,vecRates,vecDurs,vecEventTimes,dblEndT);
dblZetaP = zetatest(vecSpikeTimes,vecEventTimes);

%%
[vecT,vecMapP,matQ,vecMapZ,matMapZ]  = getLatency(vecSpikeTimes,vecEventTimes);
vecMaxZ = max(matMapZ,[],2);
dblE = mean(vecMaxZ);
dblV = var(vecMaxZ);
[vecRealP,vecRealZ,dblMode,dblBeta] = getGumbel(dblE,dblV,vecMapZ);

%calculate statistical significance using empirical quantiles
vecJitZ = matMapZ(:);
vecJitZ(isnan(vecJitZ)) = [];
vecJitZ = sort(vecJitZ);
%define p-value
vecQuantiles = nan(size(vecMapZ));
for i=1:numel(vecMapZ)
	dblZ = vecMapZ(i);
	if dblZ < min(vecJitZ) || isnan(dblZ)
		idx = 0;
	elseif dblZ > max(vecJitZ) || isinf(dblZ) || numel(vecJitZ) < 3
		idx = numel(vecJitZ);
	else
		idx = findfirst(dblZ,vecJitZ);
	end
	vecQuantiles(i) = 1 - (idx/(1+numel(vecJitZ)));
end

		
%transform to output z-score
vecP = -norminv(vecQuantiles/2);

%%
figure;maxfig;
subplot(2,3,1)
plotRaster(vecSpikeTimes,vecEventTimes,dblStimDur);
xlabel('Time after onset (s)');
ylabel('Trial #');

subplot(2,3,2)
vecISI = sort(diff(sort(vecSpikeTimes)));
plot(vecISI,(1:numel(vecISI))/numel(vecISI))
xlabel('ISI (s)');
ylabel('Quantile (cdf)');

subplot(2,3,3)
imagesc(vecT,1:numel(vecEventTimes),matQ);
xlabel('Time after onset (s)');
ylabel('Trial #');
h=colorbar;
title('ISI Quantile of observing this cessation period');

subplot(2,3,4);
plot(vecT,matMapZ,'color', [0.5 0.5 0.5 0.1]);
hold on
%plot(vecT,mean(matMapZ)-2*std(matMapZ),'--','color', [0.5 0.5 0.5]);
%plot(vecT,mean(matMapZ)+2*std(matMapZ),'--','color', [0.5 0.5 0.5]);
plot(vecT,vecMapZ,'color',lines(1));
hold off
xlabel('Time after onset (s)');
ylabel(sprintf('Significance of trial-averaged ISI quantile (%s)',getGreek('sigma')));
title('Grey=jittered, blue=real');

subplot(2,3,5)
plot(vecT,vecRealZ);
xlabel('Time after onset (s)');
ylabel(sprintf('Jitter-normalized significance (%s)',getGreek('sigma')));
title('Normalization using max of jitters with Gumbel');

subplot(2,3,6);
[vecTime,vecIFR] = getIFR(vecSpikeTimes,vecEventTimes,dblStimDur);
plot(vecTime,vecIFR);
xlabel('Time after onset (s)');
ylabel('IFR (Hz)');
title(sprintf('Zeta-p = %.6f; min p of ISI-ZETA=%.6f',dblZetaP,2-2*normcdf(max(vecRealZ))));
fixfig;