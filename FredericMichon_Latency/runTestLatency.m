
dblBaseRate = 1;
vecRates = [5 0 5];
vecDurs = [1 0.5 1.5];
dblStimDur = sum(vecDurs);
intTrials = 10;
vecEventTimes = 0:dblStimDur:(dblStimDur*9.9);
dblEndT = vecEventTimes(end)+dblStimDur;
vecSpikeTimes = getGeneratedMultiPhasicR(dblBaseRate,vecRates,vecDurs,vecEventTimes,dblEndT);

dblZetaP = zetatest(vecSpikeTimes,vecEventTimes,[],[],3);

%%
[vecT,vecMapP,matQ,vecMapZ,matMapZ]  = getLatency(vecSpikeTimes,vecEventTimes);
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
vecMean = mean(matQ);
vecSEM = std(matQ)./sqrt(size(matQ,1));
subplot(2,3,5);cla;
imagesc(vecT,1:numel(vecEventTimes),matQ);

subplot(2,3,3);cla;
plot(vecT,vecP)


subplot(2,3,2);
plot(vecWindowBinCenters,mean(matMapZ),'color', [0.5 0.5 0.5]);
hold on
plot(vecWindowBinCenters,mean(matMapZ)-2*std(matMapZ),'--','color', [0.5 0.5 0.5]);
plot(vecWindowBinCenters,mean(matMapZ)+2*std(matMapZ),'--','color', [0.5 0.5 0.5]);
plot(vecWindowBinCenters,vecMapZ,'color',lines(1));
hold off