clear all;
%close all;
strPath = 'D:\Data\Processed\TraceZeta\';
strFigPath = 'D:\Data\Results\TraceZeta\';

cellRunRand = {...
	'',...Rand 1
	'-Rand',...Rand 2
	};


%% prep
strArea = 'Primary visual';%cellUniqueAreas{intArea}; %V1, SC, Retina, Poisson, GCaMP
strStim = 'RunDriftingGratings';
matZetaP = [];
matMeanP = [];

%% load data
sDirAll=dir([strPath '*.mat']);
vecRand = contains({sDirAll.name},'Rand');
sDirReal = sDirAll(~vecRand);
sDirRand = sDirAll(vecRand);

for intRandType=1:2
	if intRandType == 1
		sDir = sDirReal;
	else
		sDir = sDirRand;
	end
	intFiles=numel(sDir);
	
	for intFile=1:intFiles
		strFile = sDir(intFile).name;
		sLoad=load([strPath strFile]);
		cellMeanP{intRandType}{intFile} = sLoad.vecMeanP;
		cellZetaP{intRandType}{intFile} = sLoad.vecZetaP;
	end
end

%% check if recording is above chance
vecResp = cellfun(@(x) sum(x<0.05),cellMeanP{1});
vecTot = cellfun(@numel,cellMeanP{1});
vecRespR = vecResp./vecTot;
pBino=bonf_holm(myBinomTest(vecResp,vecTot,0.05,'two'));
indRemRecs = pBino>0.05;
cellMeanP{1}(indRemRecs) = [];
cellMeanP{2}(indRemRecs) = [];
cellZetaP{1}(indRemRecs) = [];
cellZetaP{2}(indRemRecs) = [];

vecRandMeanP = cell2vec(cellMeanP{2});
vecRealMeanP = cell2vec(cellMeanP{1});
vecRandZetaP = cell2vec(cellZetaP{2});
vecRealZetaP = cell2vec(cellZetaP{1});

matMeanP = cat(2,vecRealMeanP,vecRandMeanP)';
matZetaP = cat(2,vecRealZetaP,vecRandZetaP)';
matMeanZ = cat(2,norminv(1-vecRealMeanP/2),norminv(1-vecRandMeanP/2))';
matZetaZ = cat(2,norminv(1-vecRealZetaP/2),norminv(1-vecRandZetaP/2))';

%remove nans
indRem = any(isnan(matZetaP),1) | any(isnan(matMeanP),1);
matMeanP(:,indRem)=[];
matZetaP(:,indRem)=[];
matMeanZ(:,indRem)=[];
matZetaZ(:,indRem)=[];

%% plot
matMeanP(matMeanZ(:)==0) = 1e-29;
matZetaP(matZetaP(:)==0) = 1e-29;
h1 =subplot(2,3,1);
matC = [0.5 0.5 0.5;...
	0 0.8 0;...
	0.8 0 0;...
	0 0 0.8];
vecColor1 = 1 + (matZetaP(1,:) < 0.05 & matMeanP(1,:) > 0.05) + 2*(matZetaP(1,:) > 0.05 & matMeanP(1,:) < 0.05) + 3*(matZetaP(1,:) < 0.05 & matMeanP(1,:) < 0.05);
scatter(matMeanZ(1,:),matZetaZ(1,:),100,vecColor1,'.');
colormap(h1,matC(1:max(vecColor1),:));
%xlim([0 1]);ylim([0 1]);
xlabel('Z-statistic mean-based t-test (\Phi^-^1(1-p/2))')
ylabel('ZETA (\zeta_c)')
title(sprintf('A) Inclusion at %s=0.05: %s=%.3f, %s=%.3f',getGreek('alpha'),getGreek('zeta'),sum(matZetaP(1,:)<0.05)/numel(matZetaP(1,:)),getGreek('mu'),sum(matMeanP(1,:)<0.05)/numel(matMeanP(1,:))))
%set(gca,'xscale','log','yscale','log');
fixfig;

h2=subplot(2,3,2)
vecColor2 = 1 + 1*(matZetaP(2,:) > 0.05 & matMeanP(2,:) < 0.05) + 2*(matZetaP(2,:) < 0.05 & matMeanP(2,:) > 0.05) + 3*(matZetaP(2,:) < 0.05 & matMeanP(2,:) < 0.05);
scatter(matMeanZ(2,:),matZetaZ(2,:),100,vecColor1,'.');
colormap(h2,matC(1:max(vecColor1),:));
%xlim([0 1]);ylim([0 1]);
xlabel('Z-statistic mean-based t-test (\Phi^-^1(1-p/2))')
ylabel('ZETA (\zeta_c)')
title(sprintf('B) False alarms at %s=0.05: %s=%.3f, %s=%.3f',getGreek('alpha'),getGreek('zeta'),sum(matZetaP(2,:)<0.05)/numel(matZetaP(2,:)),getGreek('mu'),sum(matMeanP(2,:)<0.05)/numel(matMeanP(2,:))))
%set(gca,'xscale','log','yscale','log');
fixfig;maxfig;

%% plot ROC
cellColor = {lines(1),'k'};
%vecH(intResampNpx) = subplot(4,3,intResampNpx);
subplot(2,3,3)
maxfig;
hold on;

for intTest=1:2
	if intTest == 1
		matData = matZetaP;
	else
		matData = matMeanP;
	end
	intCells = size(matData,2);
	vecBothData = cat(2,matData(1,:),matData(2,:));
	vecBothLabels = cat(2,zeros(size(matData(1,:))),ones(size(matData(1,:))));
	vecThresholds = sort(vecBothData);
	vecRealP = matData(1,:);
	vecShuffP = matData(2,:);
	
	vecTP = sum(vecRealP<vecThresholds',2)/intCells;
	vecFP = sum(vecShuffP<vecThresholds',2)/intCells;
	
	plot(vecFP,vecTP,'Color',cellColor{intTest});
	
	[dblAUC,Aci] = auc(cat(1,vecBothLabels,vecBothData)');
	vecAUC(intTest) = dblAUC;
end
hold off;
xlabel('False positive fraction');
ylabel('Inclusion fraction');
title(sprintf('C) ROC analysis'));
legend({sprintf('ZETA-test, AUC=%.2f',vecAUC(1)),sprintf('t-test, AUC=%.2f',vecAUC(2))},'location','best','interpreter','none');
fixfig;


%% save
drawnow;
export_fig(fullpath(strFigPath,'TraceZeta.tif'));
export_fig(fullpath(strFigPath,'TraceZeta.pdf'));