clear all;
%close all;
strDisk = 'F:';
strPath = [strDisk '\Data\Processed\ZETA\Latencies\'];
strFigPath = [strDisk '\Data\Results\ZETA\Latencies\'];

%% load data
sDir=dir([strPath 'ZetaDataPPL2bPoisson*']);
intFiles=numel(sDir);
mapC = redbluepurple(intFiles);
sLoad=load([strPath sDir(1).name]);
[intBins,intNeurons]=size(sLoad.matBinLatencies);

vecJitter = nan(1,intFiles);
vecBaseRate = nan(1,intFiles);
matLatencies = nan(intNeurons,intFiles);
matRealLatencies = nan(intNeurons,intFiles);
matBinLatencies = nan(intBins,intNeurons,intFiles);
intC = 0;
for intFile=1:intFiles
	strFile = sDir(intFile).name;
	
	%get data
	sLoad=load([strPath strFile]);
	
	%get params
	vecJitter(intFile) = sLoad.dblJitter;
	vecBaseRate(intFile) = sLoad.dblBaseRate;
	matLatencies(:,intFile) = sLoad.vecLatencies - sLoad.vecRealLatencies;
	matRealLatencies(:,intFile) = sLoad.vecRealLatencies;
	matBinLatencies(:,:,intFile) = sLoad.matBinLatencies - repmat(sLoad.vecRealLatencies,[intBins 1]);
	vecBinDurs = sLoad.vecBinDurs;
end
%vecJitter = vecJitter*2;

%%
figure
subplot(2,3,1)
[vecReordered,vecJitterIdx] = sort(vecJitter);
matError = squeeze(mean(abs(matBinLatencies(:,:,vecJitterIdx)),2));
[vecVal,vecIdx]=min(matError,[],1);
imagesc(log(matError)');
hold on
scatter(vecIdx,1:numel(vecIdx),'rx')
hold off
set(gca,'xtick',get(gca,'xtick')-1);
set(gca,'xticklabel',cellfun(@sprintf,cellfill('%.2f',size(get(gca,'xtick'))),vec2cell(1000*vecBinDurs(get(gca,'xtick'))),'uniformoutput',false))
colorbar
axis xy
ylabel('Jitter (ms');
xlabel('Bin size (ms)');
fixfig;
grid off;

subplot(2,3,2)
matPeakEst = squeeze(matBinLatencies(5,:,vecJitterIdx));
vecMeanEstBin = geomean(abs(matPeakEst),1);

vecMeanEstZ = geomean(abs(matLatencies(:,vecJitterIdx)),1);

plot(vecReordered,vecMeanEstBin,'k')
hold on
plot(vecReordered,vecMeanEstZ,'b')
hold off
xlabel('Jitter (ms)');
ylabel('L1 error (s)');
title(['Bin sizes: ',sprintf('%.2f; ',1000*vecBinDurs)])
fixfig;

%
mapC = redblackgreen(numel(vecBinDurs));
hs=subplot(2,3,3)
colormap(hs,mapC);
hold on
for intBin=1:numel(vecBinDurs)
	plot(vecReordered,squeeze(geomean(abs(matBinLatencies(intBin,:,vecJitterIdx)),2)),'color',mapC(intBin,:))
end
plot(vecReordered,vecMeanEstZ,'b')
hold off
h=colorbar;
title(h,'Bin sizes')
set(gca,'yscale','log')
xlabel('Jitter (ms)');
ylabel('Mean L1 error (s)');
fixfig;
maxfig;
return
strFileOut = 'Fig6J';
export_fig([strFigPath strFileOut '.tif']);
export_fig([strFigPath strFileOut '.pdf']);

median(matLatencies(:,vecJitterIdx),1)


%% calculate data
intUseData = vecJitter==max(vecJitter) & vecBaseRate == max(vecBaseRate);

vecTheseLatencies = matLatencies(:,intUseData);
vecTheseLatencies = abs(vecTheseLatencies - 0.1);
matTheseBinLatencies = matBinLatencies(:,:,intUseData);
matTheseBinLatencies = abs(matTheseBinLatencies - 0.1);
dblM = geomean(vecTheseLatencies);

dblPercError=normcdf(-1);
vecLowHigh = (sqrt(pi/2)*(getCI(vecTheseLatencies,1,dblPercError,1)))./sqrt(size(intNeurons,2));

vecMedians = geomean(matTheseBinLatencies,2);
matLowHigh = (sqrt(pi/2)*(getCI(matTheseBinLatencies,2,dblPercError,1)))./sqrt(size(intNeurons,2));
vecMSE = mean((matTheseBinLatencies - 0.1).^2,2);
figure
errorbar(vecBinDurs([1 end]),dblM*[1 1],[1 1]*(vecLowHigh(1)),[1 1]*(vecLowHigh(2)));
hold on
errorbar(vecBinDurs,vecMedians,(matLowHigh(:,1)),(matLowHigh(:,2)));
hold on
set(gca,'xscale','log')
set(gca,'yscale','log')
xlabel('Bin size (s)')
ylabel('Peak estimate (s)')
fixfig;


matDiff= matTheseBinLatencies - repmat(vecTheseLatencies',[size(matTheseBinLatencies,1) 1])

figure
errorbar(vecBinDurs,mean(matDiff,2),std(matDiff,[],2)./sqrt(intNeurons))
set(gca,'xscale','log')
xlabel('Bin size (s)')
ylabel('Peak estimate (s)')
fixfig;

