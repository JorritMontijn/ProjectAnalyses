%% aim
%{
how variable are the pop rates?

%}

%% set parameters
clear all;%close all;
cellRunTypes = {'RecTopo','SimDG18'};
intRunType = 1; %topo or sim
intOnsetType = 0; %normal (0) or rem onset (1)
cellTypes = {'Real','ShuffTid','PoissGain','Poiss'};%{'Real','Uniform','ShuffTid','PoissGain','Shuff','Poiss'}

%% define qualifying areas
boolSaveFig = true;
cellUseAreas = {...
	'Primary visual area',...
	...'posteromedial visual area',...
	};
if isfolder('F:\Drive\PopTimeCoding') && isfolder('F:\Data\Processed\Neuropixels\')
	strDataPath = 'F:\Data\Processed\Neuropixels\';
	strFigurePathSR = 'F:\Drive\PopTimeCoding\single_recs';
	strFigurePath = 'F:\Drive\PopTimeCoding\figures\';
	strTargetDataPath = 'F:\Drive\PopTimeCoding\data\';
else
	strDataPath = 'E:\DataPreProcessed\';
	strFigurePathSR = 'C:\Drive\PopTimeCoding\single_recs';
	strFigurePath = 'C:\Drive\PopTimeCoding\figures\';
	strTargetDataPath = 'C:\Drive\PopTimeCoding\data\';
end

%% load data
strArea = 'V1';
strRunType = cellRunTypes{intRunType}; %Sim or ABI or Npx?
vecRemOnsets = [0.25 0.125];
if intOnsetType == 1
	dblRemOnset = vecRemOnsets(intRunType);
else
	dblRemOnset = 0;
end
sFiles = dir ([strTargetDataPath 'D1Data_' strRunType '*.mat']);
strSelectOnset = sprintf('%.2f',dblRemOnset);

indWithOnset = cellfun(@(x) strcmp(x((end-7):end),[strSelectOnset '.mat']), {sFiles.name}');
if dblRemOnset == 0
	indWithOnset = cellfun(@(x) strcmp(x((end-7):end),[strSelectOnset '.mat']), {sFiles.name}');
	sFiles(indWithOnset) = [];
else
	sFiles(~indWithOnset) = [];
end
strOnset = sprintf('%.2f',dblRemOnset);
intRecNum = numel(sFiles);

matLogIntegral = nan(6,intRecNum);
for intFile=1:intRecNum
	%% load
	sLoad=load(fullpath(sFiles(intFile).folder,sFiles(intFile).name));
	matLogIntegral(:,intFile) = sLoad.vecLogIntegral;
	cellLoadedTypes = sLoad.cellTypes;
end

%use supplied type order
vecReorder = nan(size(cellTypes));
for i=1:numel(vecReorder)
	vecReorder(i) = find(strcmp(cellTypes{i},cellLoadedTypes));
end

%sort
matLogIntegral = matLogIntegral(vecReorder,:);

%% plot
figure;maxfig;
subplot(2,3,1)
intVarNum = size(matLogIntegral,1);
vecMean = mean(matLogIntegral,2);
plot(matLogIntegral,'color',[0.7 0.7 0.7]);
hold on
plot(vecMean,'color',lines(1));
set(gca,'xtick',1:numel(cellTypes),'xticklabel',cellTypes);
ylabel('Pop rate variability (s^2)');

%plot normalized
subplot(2,3,2)
intPoisson = find(strcmp(cellTypes,'Poiss'));
matNormLogIntegral = 100*(matLogIntegral./matLogIntegral(intPoisson,:));
vecMean = mean(matNormLogIntegral,2);
plot(matNormLogIntegral,'color',[0.7 0.7 0.7]);
hold on
plot(vecMean,'color',lines(1));
set(gca,'xtick',1:numel(cellTypes),'xticklabel',cellTypes);
ylabel('Normalized pop rate variability (%)');
title('Relative to 100%=exponential process (Poisson)');

%test
matP = nan(intVarNum,intVarNum);
for i=1:intVarNum
	for j=1:i
		[h,p]=ttest(matNormLogIntegral(i,:),matNormLogIntegral(j,:));
		matP(i,j) = p;
	end
	%matP(i,i) = 1;
end
intComps = (intVarNum^2-intVarNum)/2;
matP = matP*intComps;
subplot(2,3,3);
h = heatmap(matP);
h.CellLabelFormat = '%.3e';
h.ColorLimits = [0 0.05];
h.XDisplayLabels = cellTypes;
h.YDisplayLabels = cellTypes;
h.MissingDataColor = 'w';
colorbar;
matCol=flipud(parula);
matCol(end,:) = [0 0 0];
colormap(matCol);
title('p-values of Bonferroni-corrected paired t-tests');

%% plot mean+/-sem
subplot(2,3,4)
intVarNum = size(matLogIntegral,1);

vecMean = mean(matLogIntegral,2);
vecSem = std(matLogIntegral,[],2)./sqrt(intRecNum);
errorbar(vecMean,vecSem,'linestyle','none','marker','x');

set(gca,'xtick',1:numel(cellTypes),'xticklabel',cellTypes);
ylabel('Pop rate variability (s^2)');
xlim([0.5 intVarNum+0.5]);
%ylim([0 0.12]);

%plot normalized
subplot(2,3,5);cla;
intPoisson = find(strcmp(cellTypes,'Poiss'));
matNormLogIntegral = 100*(matLogIntegral./matLogIntegral(intPoisson,:));
vecMean = mean(matNormLogIntegral,2);
vecSem = std(matNormLogIntegral,[],2)./sqrt(intRecNum);
errorbar(vecMean(1:(intVarNum-1)),vecSem(1:(intVarNum-1)),'linestyle','none','CapSize',40);
hold on
bar((1:(intVarNum-1)),vecMean(1:(intVarNum-1)),0.5,'facecolor',lines(1));
set(gca,'xtick',1:(intVarNum-1),'xticklabel',cellTypes(1:(intVarNum-1)));
ylabel('Normalized pop rate variability (%)');
title('Relative to 100%=exponential process (Poisson)');
xlim([0.5 (intVarNum-1)+0.5]);
ylim([100 160]);




fixfig;
%%
if boolSaveFig
	%%
	export_fig(fullpath(strFigurePath,sprintf('D1_Summary%s%s.tif',strRunType,strOnset)));
	export_fig(fullpath(strFigurePath,sprintf('D1_Summary%s%s.pdf',strRunType,strOnset)));
end