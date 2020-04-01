clear all;
close all;
strPath = 'D:\Data\Results\OriMetric\Data\';
strFigPath = 'D:\Data\Results\OriMetric\';
cellUniqueAreas = {'V1','SC','Poisson','Retina','CaNM','CaDG'};
strArea = 'Retina'; %V1, SC, Retina, Poisson, GCaMP
intRunType = 1;

%set var
if intRunType == 1
	strRand = '';
elseif intRunType ==2
	strRand = '-Rand';
end
strRunType = [strArea strRand];
%% load data
%ZetaDataMSDRetinaResamp50.mat
sDir=dir([strPath 'ZetaDataMSD' strRunType 'Resamp50*']);
intFiles=numel(sDir);

%load data
strFile = sDir(1).name;
sLoad=load([strPath strFile]);
cellMSDs = sLoad.cellDeriv;
cellT = sLoad.cellInterpT;

%supersample
dblSuperRate = 0.1/1000;
dblStartT = 0;
dblEndT = 2;
vecSuperT = dblStartT:dblSuperRate:dblEndT;
intSamples = numel(vecSuperT);
intNeurons = numel(cellMSDs);
matMSDs = nan(intSamples,intNeurons);

for intNeuron=1:intNeurons
	[vecOldT,vecUniqueEntries,ic] = unique(cellT{intNeuron});
	vecOldT = [dblStartT;vecOldT;dblEndT];
	
	vecOldMSD = [cellMSDs{intNeuron}(1);cellMSDs{intNeuron}(vecUniqueEntries);cellMSDs{intNeuron}(end)];
	
	matMSDs(:,intNeuron) = interp1(vecOldT,vecOldMSD,vecSuperT);
end

%% calculate correlations
intN=50;
matCorr = corr(matMSDs);
matDist = matMSDs';
figure
for intClustNum=2:26
matLinkage = linkage(matDist,'ward');
vecClusterID = cluster(matLinkage,'maxclust',intClustNum);
[vecSortedClusts,vecReorder]=sort(vecClusterID);

matCorrReordered=matCorr(vecReorder,vecReorder);
subplot(5,5,intClustNum-1)
imagesc(matCorrReordered)
drawnow;
end

%%
%6 or 17 clusters
%close all
matDist = 1-matCorr;
%matDist = matMSDs';
for intClustNum=6%[6 17]
matLinkage = linkage(matDist,'ward');
vecClusterID = cluster(matLinkage,'maxclust',intClustNum);
[vecSortedClusts,vecReorder]=sort(vecClusterID);
matCorrReordered=matCorr(vecReorder,vecReorder);
figure
imagesc(matCorrReordered,[-1 1]);colormap(redblue)
vecEdges = find(diff(vecSortedClusts)) + 0.5;
hold on
for intC=1:numel(vecEdges)
	
	plot([1 intNeurons],vecEdges(intC)*[1 1],'g--')
end
hold off
drawnow;
end

[v,m]=fastSilhouette(matDist,vecClusterID)
plot(v)
%%
%sort
[dummy,vecReorder] = sort(vecResampNum,'ascend');
matZeta = matZeta(vecReorder,:);
matTmaxZ = matTmaxZ(vecReorder,:);
matNumSpikes = matNumSpikes(vecReorder,:);
matHzP = matHzP(vecReorder,:);
matZP = matZP(vecReorder,:);

vecR = vecR(vecReorder);
vecP = vecP(vecReorder);
vecResampNum = vecResampNum(vecReorder);
vecUpperTwoSd = vecUpperTwoSd(vecReorder);
vecFivePerc = vecFivePerc(vecReorder);
%fit exp decay
vecP0 = [vecFivePerc(1) 1 vecFivePerc(end)];
[vecParams,RESNORM,RESIDUAL,EXITFLAG] = curvefitfun(@(p,x) p(3) + p(1)*exp(-p(2)*x), vecP0, vecResampNum, vecFivePerc,[0 0 0],[10000 10 10000]);
vecPredY = vecParams(3) + vecParams(1)*exp(-vecParams(2)*vecResampNum);
intThreshResamp = vecResampNum(find(vecFivePerc<vecParams(3)*1.01,1));
intUseResamp = max(vecResampNum);%intThreshResamp

%calc significance
vecHzP = matHzP(vecResampNum==intUseResamp,:);
vecZP = matZP(vecResampNum==intUseResamp,:);
vecColor = (vecHzP < 0.05) + 2*(vecZP < 0.05) ;
if ~contains(strRunType,'Rand')
	mapC = [0.5 0.5 0.5; 1 0 0; 0 1 0; 0 0 1];
else
	mapC = [0.5 0.5 0.5; 0 1 0; 1 0 0; 0 0 1];
end
indRem = ~ismember((1:size(mapC,1))-1,vecColor);
while indRem(end)
	mapC(end,:) = [];
	indRem(end) = [];
end

%plot
%plot(vecResampNum,vecUpperTwoSd,'k--');
plot(vecResampNum,vecPredY,'k--');
dblUseCutOff = mean(vecUpperTwoSd(vecResampNum==intUseResamp));
cellLegend(end+1) = {'2 Sd'};
hold off
xlabel('Number of resamplings');
ylabel(sprintf('Maximum z-score'));
title(sprintf('%s; Z-asymptote=%.3f (resamp=%d)',strRunType,vecParams(3),intThreshResamp));
fixfig
if intRunType == 1
	%ylim([0 50]);
elseif intRunType == 3
	%ylim([0 50]);
else
	ylim([0 10]);
end
%set(gca,'xscale','log')
%xlim([1 1000])


%plot number of neurons
subplot(2,2,2)
scatter(matNumSpikes(vecResampNum==intUseResamp,:),matZeta(vecResampNum==intUseResamp,:),[],vecColor);
colormap(mapC);
[r,p]=nancorr(matNumSpikes(vecResampNum==intUseResamp,:)',matZeta(vecResampNum==intUseResamp,:)');
set(gca,'xscale','log')
xlabel('Number of spikes per neuron (count)');
ylabel(sprintf('Maximum z-score'));
title(sprintf('resamp=%d, corr(spike nr,z)=%.3f, p=%.3f',intUseResamp,r,p));
fixfig

vecBins = 0:(1/18):1;
subplot(4,2,5);
scatter(matTmaxZ(vecResampNum==intUseResamp,:),matZeta(vecResampNum==intUseResamp,:),[],vecColor);
colormap(mapC);
ylabel('Z-score');
dblPercZeta = (sum(matZP(vecResampNum==intUseResamp,:)<0.05)/numel(matZP(vecResampNum==intUseResamp,:)))*100;
dblPercRate = (sum(matHzP(vecResampNum==intUseResamp,:)<0.05)/numel(matHzP(vecResampNum==intUseResamp,:)))*100;
title(sprintf('Significant cells; (zeta): %.1f%%; (rate): %.1f%%',dblPercZeta,dblPercRate));
fixfig;

subplot(4,2,7);
vecTmaxZ = matTmaxZ(vecResampNum==intUseResamp,:);
vecSub = vecTmaxZ(matZeta(vecResampNum==intUseResamp,:)>2);
histogram(vecTmaxZ,vecBins);
xlabel('Time in trial (fraction)')
ylabel('Number of neurons');
fixfig;


subplot(2,2,4)
scatter(vecHzP,vecZP,[],vecColor);
colormap(mapC);
xlabel('P-value rate-based t-test');
ylabel('P-value Zeta');
%set(gca,'xscale','log','yscale','log')
xlim([0 1]);
ylim([0 1]);
fixfig;
maxfig();
drawnow;
export_fig(sprintf('%sSummaryFig%s.tif',strFigPath,strRunType));
print(gcf,'-dpdf', sprintf('%sSummaryFig%s.pdf',strFigPath,strRunType));