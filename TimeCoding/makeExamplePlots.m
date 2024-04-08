runHeaderPopTimeCoding;
%%
intRec=1;%1
runRecPrepNpx;
intNumN = intRespN;
%%
% make plot
dblStartT = 2000;%vecStimOnTime(1);
dblPlotDur = 1500;%vecStimOffTime(end) - vecStimOnTime(1);
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
dblScale = (dblPlotDur/dblBinW)/500;
vecFilt = normpdf((-2*dblScale):(2*dblScale),0,dblScale);
matFilt = vecFilt'*vecFilt;
matFilt = matFilt./sum(matFilt(:));
matPlotFilt = imfilt(matPlot,vecFilt,0);
vecMinR =  min(matPlotFilt,[],2);
vecMaxR =  max(matPlotFilt,[],2);
matPlotFilt = (matPlotFilt-vecMinR)./(vecMaxR-vecMinR);

%cluster
matCorr = corr(matPlotFilt');
vecMeanCorr = mean(matCorr,1);
vecNewOrder = nan(1,intPlotN);
[dummy,intN]=max(vecMeanCorr);
vecNewOrder(1) = intN;
matCorr(intN,intN) = nan;
i=1;

for i=2:intPlotN
	%%
	if mod(i,2)==0
		j=intPlotN-ceil(i/2)+1;
	else
		j=ceil(i/2);
	end
	intOldN = intN;
	vecR = matCorr(intOldN,:);
	matCorr(intOldN,:) = inf;
	matCorr(:,intOldN) = inf;
	[dummy,intN] = min(vecR);
	vecNewOrder(j) = intN;
	matCorr(intN,intN) = nan;
end
matCorr = corr(matPlot');
[vecNewOrder,matCorr] = clustsort(matPlot);
vecNewOrder2=vecNewOrder;
for i=1:numel(vecNewOrder)
	%vecNewOrder2(i) = find(vecNewOrder==i);
end
matCorr(isnan(matCorr))=0;

figure
h1=subplot(2,3,1);
imagesc(vecCenterT,1:intPlotN,-matPlotFilt);
h1.Colormap = bone;

h2=subplot(2,3,2);
imagesc(vecCenterT,1:intPlotN,matCorr,[-1 1]);
h2.Colormap = redblue;


h3=subplot(2,3,4);
imagesc(vecCenterT,1:intPlotN,-matPlotFilt(vecNewOrder2,:));
h3.Colormap = bone;

h4=subplot(2,3,5);
imagesc(vecCenterT,1:intPlotN,matCorr(vecNewOrder2,vecNewOrder2),[-1 1]);
h4.Colormap = redblue;


%cluster
