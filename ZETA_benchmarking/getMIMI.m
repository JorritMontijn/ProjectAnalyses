function [dblMIMI_P,vecLatencies,sMIMI,sRate] = getMIMI(vecSpikeTimes,matEventTimes,dblUseMaxDur,intPlot,intLatencyPeaks,vecRestrictRange,boolVerbose,intCoeffsL1,intCoeffsG1)
	
	%% prep data
	%ensure orientation
	vecSpikeTimes = vecSpikeTimes(:);
	
	%calculate stim/base difference?
	if size(matEventTimes,2) > 2
		matEventTimes = matEventTimes';
	end
	
	%trial dur
	if ~exist('dblUseMaxDur','var') || isempty(dblUseMaxDur)
		dblUseMaxDur = min(diff(matEventTimes(:,1)));
	end
	
	%get boolPlot
	if ~exist('intPlot','var') || isempty(intPlot)
		intPlot = 0;
	end
	
	%get intLatencyPeaks
	if ~exist('intLatencyPeaks','var') || isempty(intLatencyPeaks)
		if nargout > 1
			intLatencyPeaks = 2;
		else
			intLatencyPeaks = 0;
		end
	end
	
	%get boolPlot
	if ~exist('vecRestrictRange','var') || isempty(vecRestrictRange)
		vecRestrictRange = [-inf inf];
	end
	
	%get boolVerbose
	if ~exist('boolVerbose','var') || isempty(boolVerbose)
		boolVerbose = false;
	end
	
	%get coeffs
	if ~exist('intCoeffsL1','var') || isempty(intCoeffsL1)
		intCoeffsL1 = 16;
	end
	if ~exist('intCoeffsG1','var') || isempty(intCoeffsG1)
		intCoeffsG1 = 16;
	end
	
	%% build onset/offset vectors
	vecEventStarts = matEventTimes(:,1);
	
	%% check inputs
	vecLatencies = [];
	if numel(vecSpikeTimes) < 3 || numel(vecEventStarts) < 3
		dblMIMI_P = 1;
		sMIMI = struct;
		sMIMI.dblMIMI_P = dblMIMI_P;
		sMIMI.vecSpikeT = [];
		sMIMI.dblUseMaxDur = [];
		sMIMI.vecLatencyVals = [];
		sRate = [];
		warning([mfilename ':InsufficientSamples'],'Insufficient samples to calculate MIMI');
		
		%build placeholder outputs
		if numel(vecLatencies) < intLatencyPeaks
			vecLatencies(end+1:intLatencyPeaks) = nan;
		end
		
		return
	end
	
	%% fit MIMI
	sMIMI_fit = fitMIMI(intCoeffsL1,intCoeffsG1,vecSpikeTimes,vecEventStarts,dblUseMaxDur,boolVerbose);
	
	%
	intTrials = numel(vecEventStarts);
	intBins = numel(sMIMI_fit.vecY);
	dblLambda = mean(sMIMI_fit.vecY);
	vecH0 = mean(poissrnd(dblLambda,[intTrials intBins]),1);
	vecQ0 = sort(zscore(vecH0));
	vecQ1 = sort(zscore(sMIMI_fit.vecFitY));
	
	%get KS p
	[h,dblMIMI_P,ksstat] = kstest2(sMIMI_fit.vecFitY,vecH0);
	
	
	%% plot
	if intPlot > 1
		%make maximized figure
		figure
		drawnow;
		jFig = get(handle(gcf), 'JavaFrame');
		jFig.setMaximized(true);
		figure(gcf);
		drawnow;
		
		if intPlot > 2
			subplot(2,3,1)
			plotRaster(vecSpikeTimes,vecEventStarts(:,1),dblUseMaxDur,10000);
			xlabel('Time from event (s)');
			ylabel('Trial #');
			title('Spike raster plot');
			fixfig;
			grid off;
		end
		
		%plot
		%% plot
		hPeak=subplot(2,3,2);
		plot(sMIMI_fit.vecX,sMIMI_fit.vecY);
		hold on
		plot(sMIMI_fit.vecX,sMIMI_fit.vecFitY,'k--');
		hold off
		title(sprintf('Mean spiking over trials'));
		xlabel('Time from event (s)');
		ylabel('Mean spiking rate (Hz)');
		fixfig
		legend('Data','MIMI model fit');
		
		%plot mean rates
		subplot(2,3,3)
		scatter(sort(sMIMI_fit.vecY),sort(sMIMI_fit.vecFitY))
		vecLim = [min([get(gca,'xlim') get(gca,'ylim')]) max([get(gca,'xlim') get(gca,'ylim')])];
		xlim(vecLim);ylim(vecLim);
		hold on
		plot(vecLim,vecLim,'--','color',[0.5 0.5 0.5])
		hold off
		xlabel('Real spiking rate (Hz)');
		ylabel('Fitted spiking rate (Hz)');
		fixfig
		
		%% transform to probabilities from poisson process
		subplot(2,3,4)
		histx(sMIMI_fit.vecFitY)
		title(sprintf('mean fit y=%.3f,sd fit y=%.3f',mean(sMIMI_fit.vecFitY),std(sMIMI_fit.vecFitY)))
		xlabel('Fitted spiking rate (Hz)');
		ylabel('Count (bins)');
		fixfig
		
		%% are fitted values deviating from sum of poissons?
		subplot(2,3,5)
		
		hold on
		scatter(vecQ0,vecQ1)
		vecLim = [min([get(gca,'xlim') get(gca,'ylim')]) max([get(gca,'xlim') get(gca,'ylim')])];
		xlim(vecLim);ylim(vecLim);
		plot(vecLim,vecLim,'--','color',[0.5 0.5 0.5]);
		hold off
		title(sprintf('K-S test,p=%.3f',dblMIMI_P))
		xlabel('Theoretical quantile');
		ylabel('Observed quantile');
		fixfig;
		
	end
	
	%% calculate MSD statistics
	if ~isempty(sMIMI_fit) && intLatencyPeaks > 0
		%get MSD peak
		[dblPeakRate,dblPeakTime,dblPeakWidth,vecPeakStartStop,intPeakLoc,vecPeakStartStopIdx] = getPeak(sMIMI_fit.vecFitY,sMIMI_fit.vecX,vecRestrictRange);
		sRate.dblPeakRate = dblPeakRate;
		sRate.dblPeakTime = dblPeakTime;
		sRate.dblPeakWidth = dblPeakWidth;
		sRate.vecPeakStartStop = vecPeakStartStop;
		sRate.intPeakLoc = intPeakLoc;
		sRate.vecPeakStartStopIdx = vecPeakStartStopIdx;
		
		
		if ~isnan(dblPeakTime)
			%assign array data
			if intLatencyPeaks > 1
				%get onset
				[dblOnset,dblOnsetVal] = getOnset(sMIMI_fit.vecFitY,sMIMI_fit.vecX,dblPeakTime,vecRestrictRange);
				sRate.dblOnset = dblOnset;
				vecLatencies = [dblPeakTime dblOnset];
				vecLatencyVals = [sMIMI_fit.vecFitY(intPeakLoc) dblOnsetVal];
			else
				sRate.dblOnset = [nan];
				vecLatencies = [dblPeakTime];
				vecLatencyVals = [sMIMI_fit.vecFitY(intPeakLoc)];
			end
			vecLatencies = vecLatencies(1:intLatencyPeaks);
			vecLatencyVals = vecLatencyVals(1:intLatencyPeaks);
			if intPlot > 0
				axes(hPeak);
				hold on
				if intLatencyPeaks > 1
					title(sprintf('Pk=%.0fms,On=%.2fms',dblPeakTime*1000,dblOnset*1000));
				else
					title(sprintf('Pk=%.0fms',dblPeakTime*1000));
				end
				hold off
				fixfig;
				
				if intPlot > 3
					vecHandles = get(gcf,'children');
					ptrFirstSubplot = vecHandles(find(contains(get(vecHandles,'type'),'axes'),1,'last'));
					axes(ptrFirstSubplot);
					vecY = get(gca,'ylim');
					hold on;
					if intLatencyPeaks > 1,plot(dblOnset*[1 1],vecY,'r--');end
					plot(dblPeakTime*[1 1],vecY,'g--');
					hold off
				end
			end
		else
			%placeholder peak data
			sRate.dblOnset = [nan];
			vecLatencies = [nan nan];
			vecLatencyVals = [nan nan];
		end
	else
		vecLatencies = [];
		vecLatencyVals = [];
	end
	
	%check number of latencies
	if numel(vecLatencies) < intLatencyPeaks
		vecLatencies(end+1:intLatencyPeaks) = nan;
		vecLatencyVals(end+1:intLatencyPeaks) = nan;
	end
	
	%% build optional output structure
	if nargout > 2
		sMIMI = struct;
		sMIMI.dblMIMI_P = dblMIMI_P;
		sMIMI.vecSpikeT = sMIMI_fit.vecX;
		sMIMI.dblUseMaxDur = dblUseMaxDur;
		sMIMI.vecLatencyVals = vecLatencyVals;
		sMIMI = catstruct(sMIMI,sMIMI_fit);
		
	end
end

