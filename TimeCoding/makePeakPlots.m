%first run runNpxCodingAnalysisTimeT1_Agg, then run the code below
dblStart  =2050;
dblStop = 2052;
indData = vecStimTime_Real>dblStart & vecStimTime_Real < dblStop;
vecT = vecStimTime_Real(indData);
vecV = vecStimIFR_Real(indData);


%prep data
indPeaks = vecPopEventTimes_Real>dblStart & vecPopEventTimes_Real < dblStop;
vecP_T = vecPopEventTimes_Real(indPeaks);
vecP_L = vecPopEventLocsIFR_Real(indPeaks)-find(indData,1)+1;
vecP_V = vecV(vecP_L);
[matPeakDomain,indKeepPeaks] = mergepeaks(vecT,vecV,vecP_L);

%plot
figure;maxfig;
plot(vecT,vecV);
hold on
scatter(vecP_T,vecP_V,'ro')
scatter(vecP_T(indKeepPeaks),vecP_V(indKeepPeaks),'xb');
vecHeight = linspace(0,0.1,sum(indKeepPeaks));
for intPeak=1:sum(indKeepPeaks)
	h=plot(vecT([matPeakDomain(intPeak,2) matPeakDomain(intPeak,3)]),vecV([matPeakDomain(intPeak,2) matPeakDomain(intPeak,3)]),'linewidth',2);
end
xlabel('Time (s)');
ylabel('Norm IFR');
title('Example peaks - red = pre, blue post-merging; color = peak domain');
fixfig(gca,[],1);

export_fig(fullpath(strFigurePath,sprintf('T0_ExamplePeaks_%s.tif',strRec)));
export_fig(fullpath(strFigurePath,sprintf('T0_ExamplePeaks_%s.pdf',strRec)));
		