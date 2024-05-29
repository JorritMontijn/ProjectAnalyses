strRunStim = 'DG';
runHeaderPopTimeCoding;
%%
intRec=1;%1
runRecPrepNpx;
intNumN = intRespN;
%%
% make plot
dblStartT = 2000;%vecStimOnTime(1);
dblPlotDur = 500;%vecStimOffTime(end) - vecStimOnTime(1);
dblBinW = 0.1;
vecEdgesT = dblStartT:dblBinW:(dblStartT+dblPlotDur);
vecCenterT = vecEdgesT(2:end)-dblBinW/2;
intNumT = numel(vecCenterT);
matPlot = nan(intNumN,intNumT);
for intN=1:intNumN
	matPlot(intN,:) = (0+histcounts(cellSpikeTimes{intN},vecEdgesT));
end
vecMeanR = mean(matPlot,2)/dblBinW;
indRem = vecMeanR<0;
intPlotN = sum(~indRem);
matPlot(indRem,:)=[];
% plot

% filter
dblScale = 500;
dblFactor = (dblPlotDur/dblBinW)/dblScale;
vecFilt = normpdf((-2*dblFactor):(2*dblFactor),0,dblFactor);
matFilt = vecFilt'*vecFilt;
matFilt = matFilt./sum(matFilt(:));
matPlotFilt = imfilt(matPlot,vecFilt,0);
vecMinR =  min(matPlotFilt,[],2);
vecMaxR =  max(matPlotFilt,[],2);
matPlotFilt = (matPlotFilt-vecMinR)./(vecMaxR-vecMinR);

%sort
vecNewOrder = clustsort(matPlot');

figure

imagesc(vecCenterT,1:intPlotN,matPlotFilt(vecNewOrder,:));
colormap(flipud(bone));
xlabel('Time (s)')
ylabel('Neuron #');
title('Norm. act (Hz)');
colorbar
fixfig;

%% make lognormal plot 
vecOut = logmvnrnd(4,4,10000);
vecE = 0:0.5:15;
vecC = histcounts(vecOut,vecE);
stairs(vecE(2:end)-diff(vecE(1:2))/2,vecC);
