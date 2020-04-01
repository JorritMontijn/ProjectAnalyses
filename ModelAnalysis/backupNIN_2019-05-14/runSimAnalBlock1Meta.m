%% perform meta plots

%% initialize
clearvars;clc;

%% set vars
strSourceDir = 'D:\Data\Processed\V1_LIFmodel\';
strFigDir = 'D:\Data\Results\Block1\';
boolSaveFigs = true;

%% load directory contents
sFiles = dir(strSourceDir);
intFiles = numel(sFiles);
strProcessingID = '';
vecProcessed = false(1,intFiles);
cellArea = {'V1','V2','V1-V2'};
cellStim = {'Contrast','Luminance'};

%% loop through files
matMeanPredDim = nan(2,11,2,3); %full/none x int x cont/lum x V1/V2
matSDPredDim = nan(2,11,2,3); %full/none x int x cont/lum x V1/V2
%matMeanPredDimShuff = nan(2,11,2,3); %full/none x int x cont/lum x V1/V2
%matSDPredDimShuff = nan(2,11,2,3); %full/none x int x cont/lum x V1/V2

%% loop through files
for intFile = [intFiles 1:(intFiles-1)]
	strFile = sFiles(intFile).name;
	if length(strFile) > 4 && strcmpi(strFile(end-3:end),'.mat') && ~vecProcessed(intFile)
		strSimRun = getFlankedBy(strFile,'xAreaDistributed_','_201');
		strAnalBlock = getFlankedBy(strFile,'_AB','_Area');
		strArea = getFlankedBy(strFile,'_Area','.mat');
		intArea = str2double(strArea);
		if strcmp(strAnalBlock,'1') && intArea >= 0
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
			strArea = cellArea{intArea};
			
			%msg
			fprintf('Loaded file %s [%s] (%d/%d)\n',strFile,sLoad.strParam,intFile,intFiles);
			
			
			
			%% prep
			cellPredDimDepR2 = sLoad.cellPredDimDepR2;
			intIntensities = numel(cellPredDimDepR2);
			mapC = redbluepurple(intIntensities);
			vecIntensities = linspace(0,100,intIntensities);
			
			for intIntensity=1:intIntensities
				%predictability dimensionality
				vecPD = sLoad.cellPredictiveDimensions{intIntensity}(:);
				matMeanPredDim(intConnType,intIntensity,intStimType,intArea) = mean(vecPD);
				matSDPredDim(intConnType,intIntensity,intStimType,intArea) = max([std(vecPD)  mean(vecPD)/100]);
				
			end
		end
	end
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
vecParamVals = linspace(0,100,11);
for intStimType=1:2
	vecFracInfoIncFullV1 = (matMeanPredDim(1,:,intStimType,1) - matMeanPredDim(2,:,intStimType,1)) ./ matMeanPredDim(2,:,intStimType,1);
	vecFracInfoIncFullV1SD = matSDPredDim(1,:,intStimType,1) ./ matMeanPredDim(2,:,intStimType,1);
	
	vecFracInfoIncFullV2 = (matMeanPredDim(1,:,intStimType,2) - matMeanPredDim(2,:,intStimType,2)) ./ matMeanPredDim(2,:,intStimType,2);
	vecFracInfoIncFullV2SD = matSDPredDim(1,:,intStimType,2) ./ matMeanPredDim(2,:,intStimType,2);
	
	vecFracInfoIncFullV12 = (matMeanPredDim(1,:,intStimType,3) - matMeanPredDim(2,:,intStimType,3)) ./ matMeanPredDim(2,:,intStimType,3);
	vecFracInfoIncFullV12SD = matSDPredDim(1,:,intStimType,3) ./ matMeanPredDim(2,:,intStimType,3);
	
	vecY1 = vecFracInfoIncFullV1(2:end)*100;
	vecY2 = vecFracInfoIncFullV2(2:end)*100;
	vecY3 = vecFracInfoIncFullV12(2:end)*100;
	[h,p]=ttest(vecY1,vecY2);
	
	subplot(2,2,intStimType)
	plot(vecParamVals(2:end)-1,vecY1,'xr-')
	hold on
	plot(vecParamVals(2:end)+1,vecY2,'xb-')
	plot(vecParamVals(2:end)+1,vecY3,'xm-')
	plot([0 120],[0 0],'k--')
	errorbar(110,nanmean(vecY1),nanstd(vecY1),'xr')
	errorbar(110,nanmean(vecY2),nanstd(vecY2),'xb')
	errorbar(110,nanmean(vecY3),nanstd(vecY3),'xm')
	hold off
	set(gca,'xtick',vecParamVals)
	ylim([-100 100])
	ylabel(sprintf('Change in noise correlation complexity\ndue to lateral connectivity (%%)'))
	xlabel(['Stimulus ' cellStim{intStimType} ' (%)'])
	legend({'V1','V2','V1-V2'})
	title('Unshuffled')
	fixfig
end

%% save figure
if boolSaveFigs
	drawnow;
	export_fig([strFigDir 'AnalBlock1Meta3_' getDate '.tif']);
	export_fig([strFigDir 'AnalBlock1Meta3_' getDate '.pdf']);
end
