clear all;

%% generate random data
%events
dblStartEpoch = 0;
dblStopEpoch = 30;
dblEpochDur = dblStopEpoch - dblStartEpoch;
intTotSpikeNum = 20000;

%random events
vecAllSpikeTime = sort((rand(1,intTotSpikeNum)*dblEpochDur)+dblStartEpoch);

%% run
%binned counts
dblBinDur = (0.5/1000);
vecBins=(dblStartEpoch:dblBinDur:dblStopEpoch)';
vecIFR = (histcounts(vecAllSpikeTime,vecBins)./dblBinDur)';
vecTime = vecBins(2:end)-dblBinDur/2;

%find peaks
[vecPeakHeight,vecPeakLocs,w,p] = findpeaks(vecIFR);

%get pop events
vecPopEventTimes = vecTime(vecPeakLocs);
vecRandPopEventTimes = linspace(min(vecPopEventTimes)+1,max(vecPopEventTimes)-1,numel(vecPopEventTimes));
intEvents = numel(vecPopEventTimes);

% what does peak look like? plot PSTH
%raw peaks
dblStep = (1/30000);
vecEventBins = (-0.01+dblStep/2):dblStep:(0.01-dblStep/2);
intWindowSize = numel(vecEventBins)-1;
vecBinDur = diff(vecEventBins);
vecEventBinsC = vecEventBins(1:(end-1)) + vecBinDur/2;

%go through event loop
matPET = nan(intEvents,intWindowSize);
matPETR = nan(intEvents,intWindowSize);
for intEvent=1:intEvents
	%retrieve target entries
	vecTheseEdges = vecEventBins + vecPopEventTimes(intEvent);
	[vecCounts,edges] = histcounts(vecAllSpikeTime,vecTheseEdges);
	matPET(intEvent,:) = vecCounts./vecBinDur;
	
	%retrieve rand target entries
	vecTheseEdgesR = vecEventBins + vecRandPopEventTimes(intEvent);
	[vecCountsR,edges] = histcounts(vecAllSpikeTime,vecTheseEdgesR);
	matPETR(intEvent,:) = vecCountsR./vecBinDur;
end
%normal peak times
vecEventBinsC = vecEventBinsC*1000;%ms
matPopRate = sum(matPET,3);
vecMean = mean(matPopRate,1);
vecSEM = std(matPopRate,[],1)./sqrt(size(matPopRate,1));
vecMean(vecEventBinsC==0)=nan;

%random peak times
matPopRateR = sum(matPETR,3);
vecMeanR = mean(matPopRateR,1);
vecMeanR(vecEventBinsC==0)=nan;
vecSEMR = std(matPopRateR,[],1)./sqrt(size(matPopRateR,1));

%% plot
figure
plot(vecEventBinsC,vecMeanR,'color','k')
hold on
plot(vecEventBinsC,vecMean,'color',lines(1))
hold off
ylim([0 3000]);
ylabel('Pop rate (Hz)');
xlabel('Time after pop event (ms)');
drawnow;
