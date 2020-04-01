%% perform meta plots

%% initialize
clearvars;clc;

%% set vars
strSourceDir = 'D:\Data\Processed\V1_LIFmodel\';
strFigDir = 'D:\Data\Results\Block2\';
boolSaveFigs = true;

%% load directory contents
sFiles = dir(strSourceDir);
intFiles = numel(sFiles);
strProcessingID = '';
vecProcessed = false(1,intFiles);
cellArea = {'V1','V2'};
cellStim = {'Contrast','Luminance'};

%% loop through files
matMeanInfoPerSpike = nan(2,11,2,2); %full/none x int x cont/lum x V1/V2
matSDInfoPerSpike = nan(2,11,2,2); %full/none x int x cont/lum x V1/V2
matMeanInfoPerSpikeShuff = nan(2,11,2,2); %full/none x int x cont/lum x V1/V2
matSDInfoPerSpikeShuff = nan(2,11,2,2); %full/none x int x cont/lum x V1/V2

matSubMeanInfoPerSpike = nan(2,11,2,2); %full/none x int x cont/lum x V1/V2
matSubSDInfoPerSpike = nan(2,11,2,2); %full/none x int x cont/lum x V1/V2
			
cellInfoPerSpike = cell(2,11,2,2); %full/none x int x cont/lum x V1/V2
for intFile = [intFiles 1:(intFiles-1)]
	strFile = sFiles(intFile).name;
	if length(strFile) > 4 && strcmpi(strFile(end-3:end),'.mat') && ~vecProcessed(intFile)
		strSimRun = getFlankedBy(strFile,'xAreaDistributed_','_201');
		strAnalBlock = getFlankedBy(strFile,'_AB','_Area');
		strArea = getFlankedBy(strFile,'_Area','.mat');
		intArea = str2double(strArea);
		if strcmp(strAnalBlock,'2') && intArea >= 0
			%% load
			sLoad = load([strSourceDir strFile]);
			intAreaIdx = intArea;
			
			%get cont/lum
			strParam = sLoad.strParam;
			if strcmp(strParam,'Contrast')
				intStimType = 1;
			elseif strcmp(strParam,'Luminance')
				intStimType = 2;
			end
			
			%get full/none
			strConn = strSimRun((end-3):end);
			if strcmp(strConn,'Full')
				intConnType = 1;
			elseif strcmp(strConn,'None')
				intConnType = 2;
			end
			
			%get V1/V2
			strArea = sLoad.strArea;
			intArea = find(strcmp(strArea,cellArea));
			
			%msg
			fprintf('Loaded file %s [%s] (%d/%d)\n',strFile,sLoad.strParam,intFile,intFiles);
			
			%% assign data
			matMeanInfoPerSpike(intConnType,:,intStimType,intArea) = sLoad.vecMeanInfoPerSpike;
			matSDInfoPerSpike(intConnType,:,intStimType,intArea) = sLoad.vecSDInfoPerSpike;
			matMeanInfoPerSpikeShuff(intConnType,:,intStimType,intArea) = sLoad.vecMeanInfoPerSpikeShuff;
			matSDInfoPerSpikeShuff(intConnType,:,intStimType,intArea) = sLoad.vecSDInfoPerSpikeShuff;
			
			matSubMeanInfoPerSpike(intConnType,:,intStimType,intArea) = sLoad.vecSubMeanInfoPerSpike;
			matSubSDInfoPerSpike(intConnType,:,intStimType,intArea) = sLoad.vecSubSDInfoPerSpike;
			
			for intIdx=1:numel(sLoad.cellA)
				vecIPS = sLoad.cellSubI{intIdx}./sLoad.cellSubA{intIdx};
				vecIPS(isnan(vecIPS)) = [];
				cellInfoPerSpike{intConnType,intIdx,intStimType,intArea} = vecIPS;
			end
		end
	end
end

%% plot summaries
intStimIdx = size(cellInfoPerSpike,2);
intStimType = 1;
for intStimIdx=1:size(cellInfoPerSpike,2)
figure
%full screen
drawnow;
jFig = get(handle(gcf), 'JavaFrame');
jFig.setMaximized(true);
figure(gcf);
drawnow;

%
cellIPSV1 = squeeze(cellInfoPerSpike(:,intStimIdx,intStimType,1));
cellIPSV2 = squeeze(cellInfoPerSpike(:,intStimIdx,intStimType,2));

vecFullV1 = cat(1,cellIPSV1{1,:});
vecNoneV1 = cat(1,cellIPSV1{2,:});
dblStepV1 = 0.025;
vecBinsV1 = 0:dblStepV1:0.6;
vecBinsV1Edges = [(vecBinsV1 - dblStepV1/2) vecBinsV1(end) + dblStepV1/2];
vecCountsFullV1 = histcounts(vecFullV1,vecBinsV1Edges);
vecCountsNoneV1 = histcounts(vecNoneV1,vecBinsV1Edges);
dblPercIncV1 = 100*((mean(vecFullV1)/mean(vecNoneV1))-1);
[h,pV1]=ttest2(vecNoneV1,vecFullV1);

vecFullV2 = cat(1,cellIPSV2{1,:});
vecNoneV2 = cat(1,cellIPSV2{2,:});
dblStepV2 = 0.01;
vecBinsV2 = 0:dblStepV2:0.2;
vecBinsV2Edges = [(vecBinsV2 - dblStepV2/2) vecBinsV2(end) + dblStepV2/2];
vecCountsFullV2 = histcounts(vecFullV2,vecBinsV2Edges);
vecCountsNoneV2 = histcounts(vecNoneV2,vecBinsV2Edges);
dblPercIncV2 = 100*((mean(vecFullV2)/mean(vecNoneV2))-1);
[h,pV2]=ttest2(vecNoneV2,vecFullV2);

subplot(2,2,1)
plot(vecBinsV1,vecCountsFullV1,'r');
hold on
plot(vecBinsV1,vecCountsNoneV1,'b');
hold off
title(sprintf('V1, Full cont/lum subsampled to 1000 spikes, %.1f%% more I/spike',dblPercIncV1))
legend({'Full connectivity','No lateral conns'});
xlabel('Information per spike (d''^2)');
ylabel('Number of trials (count)');
fixfig

subplot(2,2,2)
plot(vecBinsV2,vecCountsFullV2,'r');
hold on
plot(vecBinsV2,vecCountsNoneV2,'b');
hold off
title(sprintf('V2, Full cont/lum subsampled to 1000 spikes, %.1f%% more I/spike',dblPercIncV2))
legend({'Full connectivity','No lateral conns'});
xlabel('Information per spike (d''^2)');
ylabel('Number of trials (count)');
fixfig
end
return

%% save figure
if boolSaveFigs
	drawnow;
	export_fig([strFigDir 'AnalBlock2Meta2_' getDate '.tif']);
	export_fig([strFigDir 'AnalBlock2Meta2_' getDate '.pdf']);
end

%% plot summaries
figure
%full screen
drawnow;
jFig = get(handle(gcf), 'JavaFrame');
jFig.setMaximized(true);
figure(gcf);
drawnow;

%%
vecParamVals = sLoad.vecParamVals;
for intStimType=1:2
vecFracInfoIncFullV1 = (matMeanInfoPerSpike(1,:,intStimType,1) - matMeanInfoPerSpike(2,:,intStimType,1)) ./ matMeanInfoPerSpike(2,:,intStimType,1);
vecFracInfoIncFullV1SD = matSDInfoPerSpike(1,:,intStimType,1) ./ matMeanInfoPerSpike(2,:,intStimType,1);

vecFracInfoIncFullV2 = (matMeanInfoPerSpike(1,:,intStimType,2) - matMeanInfoPerSpike(2,:,intStimType,2)) ./ matMeanInfoPerSpike(2,:,intStimType,2);
vecFracInfoIncFullV2SD = matSDInfoPerSpike(1,:,intStimType,2) ./ matMeanInfoPerSpike(2,:,intStimType,2);

vecY1 = vecFracInfoIncFullV1(2:end)*100;
vecY2 = vecFracInfoIncFullV2(2:end)*100;
[h,p]=ttest(vecY1,vecY2);

subplot(2,2,intStimType)
plot(vecParamVals(2:end)-1,vecY1,'xr-')
hold on
plot(vecParamVals(2:end)+1,vecY2,'xb-')
plot([0 120],[0 0],'k--')
errorbar(110,nanmean(vecY1),nanstd(vecY1),'xr')
errorbar(110,nanmean(vecY2),nanstd(vecY2),'xb')
hold off
set(gca,'xtick',vecParamVals)
ylim([-100 100])
ylabel(sprintf('Change in information/spike\ndue to lateral connectivity (%%)'))
xlabel(['Stimulus ' cellStim{intStimType} ' (%)'])
legend({'V1','V2'})
title('Unshuffled')
fixfig
end

%
vecParamVals = sLoad.vecParamVals;
for intStimType=1:2
vecFracInfoIncFullV1Shuff = (matMeanInfoPerSpikeShuff(1,:,intStimType,1) - matMeanInfoPerSpikeShuff(2,:,intStimType,1)) ./ matMeanInfoPerSpikeShuff(2,:,intStimType,1);
vecFracInfoIncFullV1SDShuff = matSDInfoPerSpikeShuff(1,:,intStimType,1) ./ matMeanInfoPerSpikeShuff(2,:,intStimType,1);

vecFracInfoIncFullV2Shuff = (matMeanInfoPerSpikeShuff(1,:,intStimType,2) - matMeanInfoPerSpikeShuff(2,:,intStimType,2)) ./ matMeanInfoPerSpikeShuff(2,:,intStimType,2);
vecFracInfoIncFullV2SDShuff = matSDInfoPerSpikeShuff(1,:,intStimType,2) ./ matMeanInfoPerSpikeShuff(2,:,intStimType,2);

vecY1 = vecFracInfoIncFullV1Shuff(2:end)*100;
vecY2 = vecFracInfoIncFullV2Shuff(2:end)*100;
[h,p]=ttest(vecY1,vecY2);

subplot(2,2,intStimType+2)
plot(vecParamVals(2:end)-1,vecY1,'xr-')
hold on
plot(vecParamVals(2:end)+1,vecY2,'xb-')
plot([0 120],[0 0],'k--')
errorbar(110,nanmean(vecY1),nanstd(vecY1),'xr')
errorbar(110,nanmean(vecY2),nanstd(vecY2),'xb')
hold off
set(gca,'xtick',vecParamVals)
ylim([-100 100])
ylabel(sprintf('Change in information/spike\ndue to lateral connectivity (%%)'))
xlabel(['Stimulus ' cellStim{intStimType} ' (%)'])
legend({'V1','V2'})
title('Shuffled')
fixfig
end

%% save figure
if boolSaveFigs
	drawnow;
	export_fig([strFigDir 'AnalBlock2Meta1_' getDate '.tif']);
	export_fig([strFigDir 'AnalBlock2Meta1_' getDate '.pdf']);
end

