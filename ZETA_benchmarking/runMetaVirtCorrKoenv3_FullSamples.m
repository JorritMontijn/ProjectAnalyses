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

%% plot
figure;
subplot(2,3,1)
scatter(vecVisLocZetaP,vecVisTimZetaP,[],vecMisMatZetaP)
%scatter(vecVisLocZetaP,vecVisTimZetaP,[],log10(vecMisMatZetaP))
colormap('bluepurplered');
xlabel('Location-modulation (z-score)');
ylabel('Time-modulation (z-score)');
h=colorbar;
h.Label.String = 'Mismatch-modulation (z-score)';
fixfig;
title(sprintf('%d neurons,%d mismatch trials',numel(vecVisLocZetaP),numel(vecMismatchT)));
%set(gca,'xscale','log','yscale','log')

%[rLT,pLT]=corr(log10(vecVisLocZetaP)',log10(vecVisTimZetaP'));
[rLT,pLT]=corr(vecVisLocZetaP',vecVisTimZetaP');
subplot(2,3,2)
vecX = linspace(min(vecVisLocZetaP),max(vecVisLocZetaP),100)';
mdl = fitlm(vecVisLocZetaP',vecVisTimZetaP');                                  % Fit Data
vecB = mdl.Coefficients.Estimate;                      % Coefficients
[vecFitY,vecFitY_CI] = predict(mdl, vecX);
%errorfill(vecX,vecFitY,vecFitY-vecFitY_CI(:,1),vecFitY_CI(:,2)-vecFitY);
plot(vecX,vecFitY,'b');
hold on
plot(vecX,vecFitY_CI(:,1),'b');
plot(vecX,vecFitY_CI(:,2),'b');
scatter(vecVisLocZetaP,vecVisTimZetaP,'k')
hold off
xlabel('Location-modulation (z-score)');
ylabel('Time-modulation (z-score)');
%set(gca,'xscale','log','yscale','log')
title(sprintf('r(Loc,Time)=%.3f,p=%.3f',rLT,pLT));
fixfig;

%[rLM,pLM]=corr(log10(vecVisLocZetaP)',log10(vecMisMatZetaP'));
[rLM,pLM]=corr(vecVisLocZetaP',vecMisMatZetaP');
subplot(2,3,4)
vecX = linspace(min(vecVisLocZetaP),max(vecVisLocZetaP),100)';
mdl = fitlm(vecVisLocZetaP',vecMisMatZetaP');                                  % Fit Data
vecB = mdl.Coefficients.Estimate;                      % Coefficients
[vecFitY,vecFitY_CI] = predict(mdl, vecX);
%errorfill(vecX,vecFitY,vecFitY-vecFitY_CI(:,1),vecFitY_CI(:,2)-vecFitY);
plot(vecX,vecFitY,'b');
hold on
plot(vecX,vecFitY_CI(:,1),'b');
plot(vecX,vecFitY_CI(:,2),'b');
scatter(vecVisLocZetaP,vecMisMatZetaP,'k');
hold off;
xlabel('Location-modulation (z-score)');
ylabel('Mismatch-modulation (z-score)');
%set(gca,'xscale','log','yscale','log')
title(sprintf('r(Loc,MisM)=%.3f,p=%.3f',rLM,pLM));
fixfig;

%[rTM,pTM]=corr(log10(vecVisTimZetaP)',log10(vecMisMatZetaP'));
[rTM,pTM]=corr(vecVisTimZetaP',vecMisMatZetaP');
subplot(2,3,5)
vecX = linspace(min(vecVisTimZetaP),max(vecVisTimZetaP),100)';
mdl = fitlm(vecVisTimZetaP',vecMisMatZetaP');                                  % Fit Data
vecB = mdl.Coefficients.Estimate;                      % Coefficients
[vecFitY,vecFitY_CI] = predict(mdl, vecX);
%errorfill(vecX,vecFitY,vecFitY-vecFitY_CI(:,1),vecFitY_CI(:,2)-vecFitY);
plot(vecX,vecFitY,'b');
xlim([0 4]);ylim([0 4]);
hold on
plot(vecX,vecFitY_CI(:,1),'b');
plot(vecX,vecFitY_CI(:,2),'b');
scatter(vecVisTimZetaP,vecMisMatZetaP,'k');
xlabel('Time-modulation (z-score)');
ylabel('Mismatch-modulation (z-score)');
%set(gca,'xscale','log','yscale','log')
title(sprintf('r(Time,MisM)=%.3f,p=%.3f',rTM,pTM));
fixfig;
maxfig;
normaxes('xy');

%{
subplot(2,3,3)

[h,pML]=ttest(vecCorr_MM_Loc);
[h,pMT]=ttest(vecCorr_MM_Tim);
[h,pLT]=ttest(vecCorr_Loc_Tim);

hold on
bplot(vecCorr_MM_Loc,1);
bplot(vecCorr_MM_Tim,2);
bplot(vecCorr_Loc_Tim,3);
%}

return
%%
if boolZ
	strZ = 'z-scored';
else
	strZ = 'raw';
end
if boolOnlyPyr
	strPyr = 'PyrOnly';
else
	strPyr = 'AllCells';
end
strFigFile = sprintf('VirtCorrAgg_FullSamples_%s_%s',strPyr,strZ);
export_fig([strFigFile '.tif'])
export_fig([strFigFile '.pdf'])
