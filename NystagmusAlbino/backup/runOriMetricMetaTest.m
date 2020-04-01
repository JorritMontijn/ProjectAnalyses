clear all;
close all;
strPath = 'D:\Data\ResultsOriMetric\';
cellSub = {'Subsamp0.2','Subsamp0.4','Subsamp0.6','Subsamp0.8',''};
for intSub=5
sDir=dir([strPath 'VisResp3DataRand' cellSub{intSub} 'Resamp*']);

intFiles=numel(sDir);
mapC = redbluepurple(intFiles);

figure
hold on
matRandVisZ =[];
matNumSpikes = [];
matTmaxZ = [];
intC = 0;
for intFile=1:intFiles
	strFile = sDir(intFile).name;
	intResampNum = str2double(getFlankedBy(strFile,'Resamp','.mat'));
	
	
	if intResampNum < 51
		intC = intC + 1;
		sLoad=load([strPath strFile]);
		
		matRandVisZ(intC,:) = abs(sLoad.vecRandVisZ);
		matNumSpikes(intC,:) = sLoad.vecNumSpikes;
		%scatter(sLoad.vecNumSpikes,abs(sLoad.vecRandVisZ),[],mapC(intFile,:))
		
		[vecR(intC),vecP(intC)] = corr(sLoad.vecNumSpikes',abs(sLoad.vecRandVisZ)');
		vecResampNum(intC) = intResampNum;
		[cellLegend, vecOneSd, vecTwoSd] = bplot(abs(sLoad.vecRandVisZ),intResampNum,'outliers','width',3);
		vecUpperTwoSd(intC) = vecTwoSd(2);
		
		for intN=1:numel(sLoad.cellRandVisZ)
			[dummy,intIdx]=max(sLoad.cellRandVisZ{intN});
			matTmaxZ(intC,intN)=intIdx/numel(sLoad.cellRandVisZ{intN});
		end
	end
	
end
%sort
[dummy,vecReorder] = sort(vecResampNum,'ascend');
matRandVisZ = matRandVisZ(vecReorder,:);
matTmaxZ = matTmaxZ(vecReorder,:);
matNumSpikes = matNumSpikes(vecReorder,:);
vecR = vecR(vecReorder);
vecP = vecP(vecReorder);
vecResampNum = vecResampNum(vecReorder);
vecUpperTwoSd = vecUpperTwoSd(vecReorder);
plot(vecResampNum,vecUpperTwoSd,'k--');
dblUseCutOff = mean(vecUpperTwoSd(vecResampNum>20));
cellLegend(end+1) = {'2 Sd'};
hold off
xlabel('Number of resamplings');
ylabel(sprintf('Maximum z-score'));
title(sprintf('Mean z for resamp>25=%.3f',dblUseCutOff));
fixfig
ylim([0 30]);
%set(gca,'xscale','log')
%xlim([1 1000])

end

%plot number of neurons
if intSub == 5
figure
intUseResamp = 25;
scatter(matNumSpikes(vecResampNum==intUseResamp,:),matRandVisZ(vecResampNum==intUseResamp,:));
[r,p]=corr(matNumSpikes(vecResampNum==intUseResamp,:)',matRandVisZ(vecResampNum==intUseResamp,:)');
set(gca,'xscale','log')
xlabel('Number of spikes per neuron (count)');
ylabel(sprintf('Maximum z-score'));
title(sprintf('resamp=%d, corr(spike nr,z)=%.3f, p=%.3f',intUseResamp,r,p));
fixfig
end

figure
vecBins = 0:0.001:1;
subplot(2,1,1);
scatter(matTmaxZ(:),matRandVisZ(:));
subplot(2,1,2);
histogram(matTmaxZ(:),vecBins);
