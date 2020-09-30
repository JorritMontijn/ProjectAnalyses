clear all;
cellFiles = {'Delier_20191015_002_Split1',...
	'Bauke_20200825_002_Split1',...
	'Bryu_20200825_002_Split1',...
	'Just_20200825_003_Split1',...
	'Bauke_20200828_004_Split1',...
	'Bryu_20200828_002_Split1',...
	'Just_20200828_002_Split1'};

boolZ = 1;
boolOnlyPyr = true;
intRemUnder = 0;

%% aggregate data
vecMisMatZetaP = [];
vecMismatchLoc = [];
vecMismatchT = [];
vecMismatchTrials = [];
vecMMNum = [];
vecVisLocZetaP = [];
vecVisTimZetaP = [];
vecInterneuron = [];

vecCorr_MM_Loc = [];
vecCorr_MM_Tim = [];
vecCorr_Loc_Tim = [];

strDisk = 'F:';
strDataPath = [strDisk '\Data\Results\ZETA\VirtCorr\'];
strDataPathRaw = [strDisk '\Data\Processed\VirtualTunnel\'];

for intFile=1:numel(cellFiles)
	strDataFileRaw = [cellFiles{intFile} 'PreProSpikes2.mat'];
	strDataFile = [cellFiles{intFile} 'ProcDatav5b.mat'];
	sLoad = load([strDataPath strDataFile]);
	sLoad2 = load([strDataPathRaw strDataFileRaw]);
	if numel(sLoad.vecMismatchT(:)) < intRemUnder,continue;end
	
	vecMismatchLoc = cat(1,vecMismatchLoc,sLoad.vecMismatchLoc(:));
	vecMismatchT = cat(1,vecMismatchT,sLoad.vecMismatchT(:));
	vecMismatchTrials = cat(1,vecMismatchTrials,sLoad.vecMismatchTrials(:));
	vecInterneuron = cat(1,vecInterneuron,cell2vec({sLoad2.sInfo.rois.red}));
	vecMMNum = cat(1,vecMMNum,numel(sLoad.vecMismatchT(:)));
	%vecTrialNum = cat(1,vecMMNum,numel(sLoad.vecMismatchT(:)));
	
	vecCorr_MM_Loc = cat(1,vecCorr_MM_Loc,corr(sLoad.vecVisLocZetaP(:),sLoad.vecMisMatZetaP(:)));
	vecCorr_MM_Tim = cat(1,vecCorr_MM_Tim,corr(sLoad.vecMisMatZetaP(:),sLoad.vecVisTimZetaP(:)));
	vecCorr_Loc_Tim = cat(1,vecCorr_Loc_Tim,corr(sLoad.vecVisLocZetaP(:),sLoad.vecVisTimZetaP(:)));
boolZ
	if boolZ == 1
		vecMisMatZetaP = cat(1,vecMisMatZetaP,nanzscore(sLoad.vecMisMatZetaP(:)));
		vecVisLocZetaP = cat(1,vecVisLocZetaP,nanzscore(sLoad.vecVisLocZetaP(:)));
		vecVisTimZetaP = cat(1,vecVisTimZetaP,nanzscore(sLoad.vecVisTimZetaP(:)));
	elseif boolZ == 2
		vecMisMatZetaP = cat(1,vecMisMatZetaP,zscore(sLoad.vecMisMatZetaP(:)));
		vecVisLocZetaP = cat(1,vecVisLocZetaP,zscore(sLoad.vecVisLocZetaP(:)));
		vecVisTimZetaP = cat(1,vecVisTimZetaP,zscore(sLoad.vecVisTimZetaP(:)));
	else
		vecMisMatZetaP = cat(1,vecMisMatZetaP,sLoad.vecMisMatZetaP(:));
		vecVisLocZetaP = cat(1,vecVisLocZetaP,sLoad.vecVisLocZetaP(:));
		vecVisTimZetaP = cat(1,vecVisTimZetaP,sLoad.vecVisTimZetaP(:));
	end
end

vecMisMatZetaP = vecMisMatZetaP';
vecMismatchLoc = vecMismatchLoc';
vecMismatchT = vecMismatchT';
vecMismatchTrials = vecMismatchTrials';
vecVisLocZetaP = vecVisLocZetaP';
vecVisTimZetaP = vecVisTimZetaP';
vecInterneuron = logical(vecInterneuron)';

%% remove interneurons
if boolOnlyPyr
vecMisMatZetaP(vecInterneuron) = [];
vecVisLocZetaP(vecInterneuron) = [];
vecVisTimZetaP(vecInterneuron) = [];
end

%% remove nans
vecIsnan = isnan(vecMisMatZetaP) | isnan(vecVisLocZetaP) | isnan(vecVisTimZetaP);
vecMisMatZetaP(vecIsnan) = [];
vecVisLocZetaP(vecIsnan) = [];
vecVisTimZetaP(vecIsnan) = [];

%% calc distros
[vecVisLocSorted,vecReorderVisLoc] = sort(vecVisLocZetaP,'ascend');
vecVL_cdf = (1:numel(vecVisLocSorted))/numel(vecVisLocSorted);
vecVL_pdf = diff([0 vecVL_cdf]) ./ diff([0 vecVisLocSorted]);

[vecMisMatSorted,vecReorderMisMat] = sort(vecMisMatZetaP,'ascend');
vecMM_cdf = (1:numel(vecMisMatSorted))/numel(vecMisMatSorted);
vecMM_pdf = diff([0 vecVL_cdf]) ./ diff([0 vecMisMatSorted]);

%% plot
[h,dblCorrP_MM_Loc] = ttest(vecCorr_MM_Loc);
[h,dblCorrP_MM_Tim] = ttest(vecCorr_MM_Tim);
[h,dblCorrP_Tim_Loc] = ttest(vecCorr_Loc_Tim);
vecCorrP = [dblCorrP_Tim_Loc dblCorrP_MM_Tim dblCorrP_MM_Loc];
close all;pause(0.5);drawnow;
hFigAgg = figure;
hFigSupp = figure;drawnow;

%% run comparisons
for intComp=1:3
	if intComp == 1
	vecData1 = vecVisLocZetaP;
	vecData2 = vecVisTimZetaP;
	strVar1 = 'Location';
	strVar2 = 'Time';
	elseif intComp == 2
	vecData1 = vecVisTimZetaP;
	vecData2 = vecMisMatZetaP;
	strVar1 = 'Time';
	strVar2 = 'Mismatch';
	
	elseif intComp == 3
	vecData1 = vecVisLocZetaP;
	vecData2 = vecMisMatZetaP;
	strVar1 = 'Location';
	strVar2 = 'Mismatch';
	
	end
	
	
vecLim = [-4 4];
%hF=figure
%subplot(2,3,1)
%plot([vecVisLocSorted],[vecVL_cdf]);

%subplot(2,3,2)
%histogram(vecVisLocSorted)

%subplot(2,3,3)
[bandwidthVL,densityVL,xmesh,cdfVL]=kde(vecData1,2^12,vecLim(1),vecLim(2));
%kde(vecData1,2^12,vecLim(1),vecLim(2));

%subplot(2,3,4)
%plot([vecMisMatSorted],[vecMM_cdf]);

%subplot(2,3,5)
%histogram(vecMisMatSorted)

%subplot(2,3,6)
[bandwidthMM,densityMM,ymesh,cdfMM]=kde(vecData2,2^12,vecLim(1),vecLim(2));
%kde(vecData2,2^12,vecLim(1),vecLim(2));
%drawnow;close(hF);pause(0.1);
%% convolve 2D points with VL/MM bandwidths to construct real distro
dblStepX = median(diff(xmesh));
vecEdgeX = [(xmesh(1) - dblStepX) xmesh] + dblStepX/2;

dblStepY = median(diff(ymesh));
vecEdgeY = [(ymesh(1) - dblStepY) ymesh] + dblStepY/2;
[matSignal] = histcounts2(vecData1,vecData2,vecEdgeX,vecEdgeY);

vecFiltMeshY = linspace(-bandwidthMM*5,bandwidthMM*5,(bandwidthMM*10)/dblStepX);
vecFilterY = normpdf(vecFiltMeshY,0,bandwidthMM);

vecFiltMeshX = linspace(-bandwidthVL*5,bandwidthVL*5,(bandwidthVL*10)/dblStepY);
vecFilterX = normpdf(vecFiltMeshX,0,bandwidthVL);

matFilter = vecFilterY' * vecFilterX;
matDensReal = imfiltreflectpad(matSignal,matFilter)';

%% cut mesh
intPoints = numel(vecVisLocSorted);
h1=subplot(3,4,1+4*(intComp-1));
matDens2D = densityMM * densityVL';
matDensReal = (matDensReal./sum(matDensReal(:)))*intPoints;
matDens2D = (matDens2D./sum(matDens2D(:)))*intPoints;
%matDens2D = imnorm(matDens2D);%(matDens2D ./ sum(matDens2D(:)))+eps;
imagesc(xmesh,ymesh,matDens2D);
colormap(h1,redwhite);
axis xy
hold on
colorbar
hold off;
xlabel([strVar1 '-mod (z-score)']);
ylabel([strVar2 '-mod (z-score)']);
%set(gca,'xscale','log','yscale','log')
title(sprintf('Independent null-distro'));
fixfig;

h2=subplot(3,4,2+4*(intComp-1));
%matDensReal = imnorm(matDensReal);%
%matDensReal = (matDensReal ./ sum(matDensReal(:)))+eps;
imagesc(xmesh,ymesh,matDensReal);
colormap(h2,redwhite);
axis xy
hold on
scatter(vecData1,vecData2,'k.');
colorbar
hold off;
xlabel([strVar1 '-mod (z-score)']);
ylabel([strVar2 '-mod (z-score)']);
%set(gca,'xscale','log','yscale','log')
title(sprintf('Real joint distribution'));
fixfig;

h3=subplot(3,4,3+4*(intComp-1));
matDiff = (matDensReal-matDens2D);
imagesc(xmesh,ymesh,matDensReal-matDens2D,max(abs(matDiff(:)))*[-1 1]);
colormap(h3,redblue);
axis xy
hold on
colorbar
hold off;
xlabel([strVar1 '-mod (z-score)']);
ylabel([strVar2 '-mod (z-score)']);
%set(gca,'xscale','log','yscale','log')
title(sprintf('%sDensity',getGreek('Delta')));
fixfig;


%{
h4=subplot(2,3,4)
scatter(vecData1,vecData2,'k.');
hold off;
xlabel([strVar1 '-modulation (z-score)']);
ylabel([strVar2 '-modulation (z-score)']);
%set(gca,'xscale','log','yscale','log')
title(sprintf('r(%s,%s)=%.3f,p=%.3f',strVar1,strVar2,rLM,pLM));
fixfig;
%[rLM,pLM]=corr(log10(vecData1)',log10(vecData2'));
[rLM,pLM]=corr(vecData1',vecData2');
vecX = linspace(min(vecData1),max(vecData1),100)';
mdl = fitlm(vecData1',vecData2');                                  % Fit Data
vecB = mdl.Coefficients.Estimate;                      % Coefficients
[vecFitY,vecFitY_CI] = predict(mdl, vecX);
%errorfill(vecX,vecFitY,vecFitY-vecFitY_CI(:,1),vecFitY_CI(:,2)-vecFitY);
plot(vecX,vecFitY,'b');
plot(vecX,vecFitY_CI(:,1),'b');
plot(vecX,vecFitY_CI(:,2),'b');
scatter(vecData1,vecData2,'k.');
hold off;
xlabel([strVar1 '-modulation (z-score)']);
ylabel([strVar2 '-modulation (z-score)']);
%set(gca,'xscale','log','yscale','log')
title(sprintf('r(%s,%s)=%.3f,p=%.3f',strVar1,strVar2,rLM,pLM));
fixfig;
%}

%% plot sampling from KDE distro
%x: loc
%y: mm
%{
figure
for intR=1:9
subplot(3,3,intR)
vecVisLocSampled = bandwidthVL*randn(size(vecData1)) + vecData1(randperm(numel(vecData1)));
vecMisMatSampled = bandwidthMM*randn(size(vecData2)) + vecData2(randperm(numel(vecData2)));
scatter(vecVisLocSampled,vecMisMatSampled,'k.');
xlim(vecLim);ylim(vecLim);fixfig;

title(sprintf('r=%.3f',corr(vecVisLocSampled',vecMisMatSampled')))
end
%}

%% do actual bootstrapping of correlation coefficient
%{
intIters = 10000;
vecRandR = nan(1,intIters);
for intIter=1:10000
	vecVisLocSampled = bandwidthVL*randn(size(vecData1)) + vecData1(randperm(numel(vecData1)));
	vecMisMatSampled = bandwidthMM*randn(size(vecData2)) + vecData2(randperm(numel(vecData2)));
	vecRandR(intIter) = corr(vecVisLocSampled',vecMisMatSampled');
end

subplot(2,3,5)
histogram(vecRandR);
hold on
plot(rLM*[1 1],[0 max(get(gca,'ylim'))],'b');
hold off

dblZ_VM = (rLM-mean(vecRandR)) / std(vecRandR);
dblP_VM = (1-normcdf(abs(dblZ_VM)))*2;
fixfig;
title(sprintf('Z=%.3f, p=%.3f',dblZ_VM,dblP_VM))
%}
%% do bootstrapping of upper quantile count
dblMeanVisLoc = mean(vecData1);
dblMeanMisMat = mean(vecData2);
intIters = 10000;
vecCountsURQ = nan(1,intIters);
for intIter=1:intIters
	vecVisLocSampled = bandwidthVL*randn(size(vecData1)) + vecData1(randperm(numel(vecData1)));
	vecMisMatSampled = bandwidthMM*randn(size(vecData2)) + vecData2(randperm(numel(vecData2)));
	vecCountsURQ(intIter) = sum((vecVisLocSampled > dblMeanVisLoc) & (vecMisMatSampled > dblMeanMisMat));
end
intRealURQ = sum((vecData1 > dblMeanVisLoc) & (vecData2 > dblMeanMisMat)); 

subplot(3,4,4+4*(intComp-1))
histogram(vecCountsURQ);drawnow;
hold on
plot(intRealURQ*[1 1],[0 max(get(gca,'ylim'))],'b');
hold off
xlim([120 200]);
xlabel(['Number of points in URQ' newline ' ']);

dblZ_VM_URQ = (intRealURQ-mean(vecCountsURQ)) / std(vecCountsURQ);
dblP_VM_URQ = (1-normcdf(abs(dblZ_VM_URQ)))*2;
fixfig;
title(sprintf('URQ-count, Z=%.3f, p=%.3f',dblZ_VM_URQ,dblP_VM_URQ))

%% aggregate summary figure
figure(hFigAgg);drawnow;
hAgg=subplot(2,3,intComp);

[rLM,pLM]=corr(vecData1',vecData2');
matDiff = (matDensReal-matDens2D);
imagesc(xmesh,ymesh,matDensReal-matDens2D,max(abs(matDiff(:)))*[-1 1]);
hold on;
colormap(hAgg,redblue);
axis xy
colorbar
scatter(vecData1,vecData2,'k.');
hold off;
xlabel([strVar1 '-modulation (z-score)']);
ylabel([strVar2 '-modulation (z-score)']);
xlim(vecLim);ylim(vecLim);
%set(gca,'xscale','log','yscale','log')
title(sprintf('r(%s,%s)=%.3f,URQ-p=%.3f,corr-p=%.3f',strVar1,strVar2,rLM,dblP_VM_URQ,vecCorrP(intComp)));
fixfig;
figure(hFigSupp);
fprintf('%s-%s, URQ-p=%.3e; R-p=%.3e\n',strVar1,strVar2,dblP_VM_URQ,vecCorrP(intComp));
end
maxfig(hFigSupp,0.9);
maxfig(hFigAgg);
return
%% save figures
figure(hFigSupp);
drawnow;
export_fig([strDataPath 'SuppFigVirtCorr.tif']);
export_fig([strDataPath 'SuppFigVirtCorr.pdf']);


figure(hFigAgg);
drawnow;
export_fig([strDataPath 'VirtCorrAggKDE.tif']);
export_fig([strDataPath 'VirtCorrAggKDE.pdf']);
