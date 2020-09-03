cellFiles = {'Bryu_20200825_002_Split1PreProSpikes2',...
	'Delier_20191015_002_Split1PreProSpikes2',...
	'Bauke_20200825_002_Split1PreProSpikes2',...
	'Just_20200825_003_Split1PreProSpikes2'};

%% aggregate data
vecMisMatZetaP = [];
vecMismatchLoc = [];
vecMismatchT = [];
vecMismatchTrials = [];
vecVisLocZetaP = [];
vecVisTimZetaP = [];
vecInterneuron = [];
strDisk = 'F:';
strDataPath = [strDisk '\Data\Results\ZETA\VirtCorr\'];
strDataPathRaw = [strDisk '\Data\Processed\VirtualTunnel\'];

for intFile=1:numel(cellFiles)
	
	strDataFileRaw = [cellFiles{intFile} '.mat'];
	strDataFile = [cellFiles{intFile} 'ProcDatav3.mat'];
	sLoad = load([strDataPath strDataFile]);
	sLoad2 = load([strDataPathRaw strDataFileRaw]);
	
	vecMisMatZetaP = cat(1,vecMisMatZetaP,zscore(sLoad.vecMisMatZetaP(:)));
	vecMismatchLoc = cat(1,vecMismatchLoc,sLoad.vecMismatchLoc(:));
	vecMismatchT = cat(1,vecMismatchT,sLoad.vecMismatchT(:));
	vecMismatchTrials = cat(1,vecMismatchTrials,sLoad.vecMismatchTrials(:));
	vecVisLocZetaP = cat(1,vecVisLocZetaP,zscore(sLoad.vecVisLocZetaP(:)));
	vecVisTimZetaP = cat(1,vecVisTimZetaP,zscore(sLoad.vecVisTimZetaP(:)));
	vecInterneuron = cat(1,vecInterneuron,cell2vec({sLoad2.sInfo.rois.red}));
end

vecMisMatZetaP = vecMisMatZetaP';
vecMismatchLoc = vecMismatchLoc';
vecMismatchT = vecMismatchT';
vecMismatchTrials = vecMismatchTrials';
vecVisLocZetaP = vecVisLocZetaP';
vecVisTimZetaP = vecVisTimZetaP';
vecInterneuron = logical(vecInterneuron)';

%% remove interneurons
vecMisMatZetaP(vecInterneuron) = [];
vecVisLocZetaP(vecInterneuron) = [];
vecVisTimZetaP(vecInterneuron) = [];

%% calc distros
[vecVisLocSorted,vecReorderVisLoc] = sort(vecVisLocZetaP,'ascend');
vecVL_cdf = (1:numel(vecVisLocSorted))/numel(vecVisLocSorted);
vecVL_pdf = diff([0 vecVL_cdf]) ./ diff([0 vecVisLocSorted]);

[vecMisMatSorted,vecReorderMisMat] = sort(vecMisMatZetaP,'ascend');
vecMM_cdf = (1:numel(vecMisMatSorted))/numel(vecMisMatSorted);
vecMM_pdf = diff([0 vecVL_cdf]) ./ diff([0 vecMisMatSorted]);

%% plot
vecLim = [-4 4];
figure
subplot(2,3,1)
plot([vecVisLocSorted],[vecVL_cdf]);

subplot(2,3,2)
histogram(vecVisLocSorted)

subplot(2,3,3)
[bandwidthVL,densityVL,xmesh,cdfVL]=kde(vecVisLocZetaP,2^12,vecLim(1),vecLim(2));
kde(vecVisLocZetaP,2^12,vecLim(1),vecLim(2));

subplot(2,3,4)
plot([vecMisMatSorted],[vecMM_cdf]);

subplot(2,3,5)
histogram(vecMisMatSorted)

subplot(2,3,6)
[bandwidthMM,densityMM,ymesh,cdfMM]=kde(vecMisMatZetaP,2^12,vecLim(1),vecLim(2));
kde(vecMisMatZetaP,2^12,vecLim(1),vecLim(2));

%% cut mesh

figure
matDens2D = densityMM * densityVL';
imagesc(xmesh,ymesh,matDens2D);
colormap(redwhite);
axis xy
hold on

%[rLM,pLM]=corr(log10(vecVisLocZetaP)',log10(vecMisMatZetaP'));
[rLM,pLM]=corr(vecVisLocZetaP',vecMisMatZetaP');
vecX = linspace(min(vecVisLocZetaP),max(vecVisLocZetaP),100)';
mdl = fitlm(vecVisLocZetaP',vecMisMatZetaP');                                  % Fit Data
vecB = mdl.Coefficients.Estimate;                      % Coefficients
[vecFitY,vecFitY_CI] = predict(mdl, vecX);
%errorfill(vecX,vecFitY,vecFitY-vecFitY_CI(:,1),vecFitY_CI(:,2)-vecFitY);
plot(vecX,vecFitY,'b');
plot(vecX,vecFitY_CI(:,1),'b');
plot(vecX,vecFitY_CI(:,2),'b');
scatter(vecVisLocZetaP,vecMisMatZetaP,'k.');
hold off;
xlabel('Location-modulation (z-score)');
ylabel('Mismatch-modulation (z-score)');
%set(gca,'xscale','log','yscale','log')
title(sprintf('r(Loc,MisM)=%.3f,p=%.3f',rLM,pLM));
fixfig;

%% plot sampling from KDE distro
%x: loc
%y: mm
figure
for intR=1:9
subplot(3,3,intR)
vecVisLocSampled = bandwidthVL*randn(size(vecVisLocZetaP)) + vecVisLocZetaP(randperm(numel(vecVisLocZetaP)));
vecMisMatSampled = bandwidthMM*randn(size(vecMisMatZetaP)) + vecMisMatZetaP(randperm(numel(vecMisMatZetaP)));
scatter(vecVisLocSampled,vecMisMatSampled,'k.');
xlim(vecLim);ylim(vecLim);fixfig;

title(sprintf('r=%.3f',corr(vecVisLocSampled',vecMisMatSampled')))
end

%% do actual bootstrapping of correlation coefficient
intIters = 10000;
vecRandR = nan(1,intIters);
for intIter=1:10000
	vecVisLocSampled = bandwidthVL*randn(size(vecVisLocZetaP)) + vecVisLocZetaP(randperm(numel(vecVisLocZetaP)));
	vecMisMatSampled = bandwidthMM*randn(size(vecMisMatZetaP)) + vecMisMatZetaP(randperm(numel(vecMisMatZetaP)));
	vecRandR(intIter) = corr(vecVisLocSampled',vecMisMatSampled');
end

figure
histogram(vecRandR);
hold on
plot(rLM*[1 1],[0 max(get(gca,'ylim'))],'b');
hold off

dblZ_VM = (rLM-mean(vecRandR)) / std(vecRandR);
dblP_VM = (1-normcdf(abs(dblZ_VM)))*2;
fixfig;
title(sprintf('Z=%.3f, p=%.3f',dblZ_VM,dblP_VM))

%% do bootstrapping of upper quantile count
dblMeanVisLoc = mean(vecVisLocZetaP);
dblMeanMisMat = mean(vecMisMatZetaP);
intIters = 10000;
vecCountsURQ = nan(1,intIters);
for intIter=1:10000
	vecVisLocSampled = bandwidthVL*randn(size(vecVisLocZetaP)) + vecVisLocZetaP(randperm(numel(vecVisLocZetaP)));
	vecMisMatSampled = bandwidthMM*randn(size(vecMisMatZetaP)) + vecMisMatZetaP(randperm(numel(vecMisMatZetaP)));
	vecCountsURQ(intIter) = sum((vecVisLocSampled > dblMeanVisLoc) & (vecMisMatSampled > dblMeanMisMat));
end
intRealURQ = sum((vecVisLocZetaP > dblMeanVisLoc) & (vecMisMatZetaP > dblMeanMisMat)); 

figure
histogram(vecCountsURQ);
hold on
plot(intRealURQ*[1 1],[0 max(get(gca,'ylim'))],'b');
hold off

dblZ_VM_URQ = (intRealURQ-mean(vecCountsURQ)) / std(vecCountsURQ);
dblP_VM_URQ = (1-normcdf(abs(dblZ_VM_URQ)))*2;
fixfig;
title(sprintf('URQ, Z=%.3f, p=%.3f',dblZ_VM_URQ,dblP_VM_URQ))
