%% exploratory analysis, no proper controls

%% load data
strDataPath = 'E:\DataPreProcessed';
sFiles = dir(fullpath(strDataPath,'*_AP.mat'));
if ~exist('sExp','var') || isempty(sExp)
	sExp = [];
	for intFile=1:numel(sFiles)
		fprintf('Loading %d/%d: %s [%s]\n',intFile,numel(sFiles),sFiles(intFile).name,getTime);
		sLoad = load(fullpath(sFiles(intFile).folder,sFiles(intFile).name));
		if ~isfield(sLoad.sAP,'sPupil') || isempty(sLoad.sAP.sPupil),continue;end
		if isempty(sExp)
			sExp = sLoad.sAP;
		else
			sExp(end+1) = sLoad.sAP;
		end
	end
end

%MP_20200115 eye tracking remove last stimulus (gunk in eye)
cellUseForEyeTrackingMP = {'20191120','20191121','20191122','20191210','20191211','20191212','20191213','20191216','20191217','20200116','20200116R02'}; %don't forget to set high vid lum as blinks
cellUseForEyeTrackingMA = {'20210212','20210215','20210218','20210220','20210225','20210301'};
cellUseForEyeTracking = cat(2,cellUseForEyeTrackingMA,cellUseForEyeTrackingMP);
strTargetPath = 'D:\Data\Results\AlbinoProject';

%best rec BL6: 20191216B5 (rec 17)
%best rec DBA: 20210212B2 (rec 5)

%% define area categories
%cortex
cellUseAreas = [];
cellUseAreas{1} = {'Primary visual','Posteromedial visual','anteromedial visual'};
%NOT
cellUseAreas{2} = {'nucleus of the optic tract'};
cellAreaGroups = {'Vis. ctx','NOT'};
cellAreaGroupsAbbr = {'Ctx','NOT'};
cellSubjectGroups = {'BL6','DBA'};

vecColAlb = [0.9 0 0];
vecColBl6 = lines(1);

%% pre-allocate
intStimNr = 24;
vecAggPerfWt = nan(0,numel(cellUseAreas),1);
matAggPerfWt = nan(0,numel(cellUseAreas),intStimNr);
vecAggPerfAlb = nan(0,numel(cellUseAreas),1);
matAggPerfAlb = nan(0,numel(cellUseAreas),intStimNr);
matAggConfusionWt = nan(0,numel(cellUseAreas),intStimNr,intStimNr);
matAggConfusionAlb = nan(0,numel(cellUseAreas),intStimNr,intStimNr);

cellAggPrefWt = cell(0,numel(cellUseAreas));
cellAggPrefAlb = cell(0,numel(cellUseAreas));

cellAggMeanRespWt = cellfill(nan(0,intStimNr),[1 3]);
cellAggMeanRespAlb = cellfill(nan(0,intStimNr),[1 3]);
				
%% run
cellNameAP = arrayfun(@(x) x.sJson.file_preproAP,sExp,'uniformoutput',false);
cellExperiment = arrayfun(@(x) x.sJson.experiment,sExp,'uniformoutput',false);
cellRemove = {};%{'RecMA5_2021-03-01R01_g0_t0'};
indRemRecs = contains(cellExperiment,cellRemove);
indRemRecs2 = ~contains(cellNameAP,cellUseForEyeTracking);
cellSubjectType = arrayfun(@(x) x.sJson.subjecttype,sExp,'uniformoutput',false);
dblAverageMouseHeadTiltInSetup = -15;
for intSubType=1:2
	intPopCounter = 0;
	if intSubType == 1
		intBestRec = 17;
		strSubjectType = 'BL6';
		dblOffsetT=0;
	elseif intSubType == 2
		intBestRec = 5;
		strSubjectType = 'DBA';
		dblOffsetT=0;
	end
	indUseRecs = contains(cellSubjectType,strSubjectType);
	vecRunRecs = find(indUseRecs & ~(indRemRecs | indRemRecs2));
	%vecRunRecs = intBestRec;
	for intRecIdx=1:numel(vecRunRecs)
		intRec=vecRunRecs(intRecIdx)
		sRec = sExp(intRec);
		strName=[sRec.sJson.subject '_' sRec.sJson.date];
		cellStimType = cellfun(@(x) x.strExpType,sRec.cellBlock,'uniformoutput',false);
		vecBlocksDG = find(contains(cellStimType,'driftinggrating','IgnoreCase',true));
		for intBlockIdx=1:numel(vecBlocksDG)
			intBlock = vecBlocksDG(intBlockIdx)
			sBlock = sRec.cellBlock{intBlock};
			intPopCounter = intPopCounter + 1;
			
			if isfield(sBlock,'vecPupilStimOn')
				vecPupilStimOn = sBlock.vecPupilStimOn;
				vecPupilStimOff = sBlock.vecPupilStimOff;
			else
				vecPupilStimOn = sBlock.vecStimOnTime;
				vecPupilStimOff = sBlock.vecStimOffTime;
			end
			
			%% get pupil data
			dblSampNi = str2double(sRec.sSources.sMeta.niSampRate);
			dblFirstSamp = str2double(sRec.sSources.sMeta.firstSample);
			vecPupilTime = sRec.sPupil.vecTime;
			vecPupilLocX = sRec.sPupil.vecCenterX;
			vecPupilLocY = sRec.sPupil.vecCenterY;
			vecPupilSize = sRec.sPupil.vecRadius;
			vecPupilSync = sRec.sPupil.vecSyncLum;
			if isfield(sRec.sPupil,'vecBlinks') && ~all(sRec.sPupil.vecBlinks==0)
				vecBlinks = sRec.sPupil.vecBlinks;
			else
				%filter absvidlum
				dblLowPass = 0.01/(1/median(diff(vecPupilTime)));
				[fb,fa] = butter(2,dblLowPass,'high');
				vecAbsVidLum = zscore(filtfilt(fb,fa, sRec.sPupil.sRaw.vecPupilAbsVidLum));
				vecBlinks = vecAbsVidLum > 5;
			end
			dblFs = 1/median(diff(vecPupilTime));
			
			% get blinking
			vecBinEdges = [-inf vecPupilStimOn vecPupilStimOn(end)+median(diff(vecPupilStimOn)) inf];
			[vecCounts,vecMeans,vecSDs,cellVals,cellIDs] = makeBins(vecPupilTime,double(vecBlinks),vecBinEdges);
			vecBlinkFractionPerTrial = vecMeans(2:(end-1));
			
			%remove trials with blinking
			indRemTrials = vecBlinkFractionPerTrial > 0.1;
			
			%% prep data
			%split by ori
			sTrialObjects = sBlock.sStimObject(sBlock.vecTrialStimTypes);
			vecOrientation = cell2vec({sTrialObjects.Orientation});
			vecOrientation(indRemTrials) = [];
			[vecOriIdx,vecUnique,vecCounts,cellSelect,vecRepetition] = val2idx(vecOrientation);
			if numel(vecUnique) ~= 24,continue,end
			
			%get data matrix
			cellSpikeT = {sRec.sCluster(:).SpikeTimes};
			vecStimOnTime = sBlock.vecStimOnTime(~indRemTrials);
			vecStimOffTime = sBlock.vecStimOffTime(~indRemTrials);
			dblDur = median(vecStimOffTime-vecStimOnTime);
			matData = getSpikeCounts(cellSpikeT,vecStimOnTime,dblDur);
			
			%include?
			indUseCells = arrayfun(@(x) x.Violations1ms < 0.25 & abs(x.NonStationarity) < 0.25,sRec.sCluster(:));
			
			%% split cells into areas
			%build cell vectors
			cellCellsPerArea = cell(1,numel(cellUseAreas));
			cellAreasPerCluster = {sRec.sCluster.Area};
			for intArea=1:numel(cellUseAreas)
				cellCellsPerArea{intArea} = contains(cellAreasPerCluster,cellUseAreas{intArea},'IgnoreCase',true);
			end
			vecCellsNrPerArea = cellfun(@sum,cellCellsPerArea);
			
			%% go through adjacent stimuli
			%pre-allocate
			vecPerf = nan(numel(cellUseAreas),1);
			matPerf = nan(numel(cellUseAreas),intStimNr);
			cellPref = cell(numel(cellUseAreas),1);
			
			%params
			intTypeCV = 2; %leave repetition out
			vecOriNoDir = mod(vecOrientation,180);
			vecTrialTypesNoDir = deg2rad(vecOriNoDir)*2;
			[vecOriIdx,vecUnique,vecCounts,cellSelect,vecRepetition] = val2idx(vecTrialTypesNoDir);
			dblLambda = 100;
			
			%[vecTrialTypes,vecUnique,vecCounts,cellSelect,vecRepetition] = val2idx(vecOriNoDir);
			[vecTrialTypes,vecUnique,vecCounts,cellSelect,vecRepetition] = val2idx(vecOrientation);
			vecPriorDistribution = vecCounts;
			intStimNr = numel(vecUnique);
			%figure
			for intArea = 1:numel(cellUseAreas)
				%% select cells
				strAreaGroup =  cellAreaGroupsAbbr{intArea};
				vecSelectCells = find(indUseCells(:) & cellCellsPerArea{intArea}(:));
				if isempty(vecSelectCells)
					continue;
				end
				
				%% calc tuning curves
				matUseData = matData(vecSelectCells,:);
				sOut = getTuningCurves(matUseData,vecOrientation);
				vecPrefOri = rad2deg(sOut.matFittedParams(:,1));
				indKeepCells = sOut.vecFitP<0.05;
				if intArea==3
					indKeepCells = true(size(indKeepCells));
				end
				vecPrefOri(~indKeepCells)=[];
				cellPref{intArea} = vecPrefOri;
				
				matOriResp = sOut.matMeanResp(indKeepCells,:);
				vecMeanHz = mean(matOriResp,2);
				indRem=vecMeanHz<1;
				if sum(~indRem) == 0,continue;end
				if strcmp(strSubjectType,'BL6')
					cellAggMeanRespWt{intArea} = cat(1,cellAggMeanRespWt{intArea},matOriResp(~indRem,:)./vecMeanHz(~indRem));
				else
					cellAggMeanRespAlb{intArea} = cat(1,cellAggMeanRespAlb{intArea},matOriResp(~indRem,:)./vecMeanHz(~indRem));
				end
				
				%% decode all
				matUseData = matData(vecSelectCells,:);
				%[dblPerformanceCV,vecDecodedIndexCV,matTemplateDistsCV,dblMeanErrorDegs,matConfusion] = doCrossValidatedDecodingTM(matUseData,vecOriIdx,intTypeCV,vecPriorDistribution);
				%vecPerfTM(intArea) = dblPerformanceCV;
				[dblPerformanceLR,vecDecodedIndexCV_LR,matPosteriorProbability,dblMeanErrorDegsLR,matConfusionLR,matWeights] = ...
					doCrossValidatedDecodingLR(matUseData,vecOrientation,intTypeCV,vecPriorDistribution,dblLambda);
				vecPerf(intArea) = dblPerformanceLR;
				if strcmp(strSubjectType,'BL6')
					matAggConfusionWt(intPopCounter,intArea,:,:) = matConfusionLR;
				else
					matAggConfusionAlb(intPopCounter,intArea,:,:) = matConfusionLR;
				end
				
				%subplot(2,3,intArea)
				%imagesc(matConfusionLR)
				%title([strAreaGroup '; LR: ' strName '_' num2str(intBlock)],'interpreter','none');
				
				%% adjacent stims
				for intStimCenter=1:intStimNr
					intStim1 = mod(intStimCenter-1,intStimNr);
					intStim2 = mod(intStimCenter+1,intStimNr);
					if intStim1==0,intStim1=intStimNr;end
					if intStim2==0,intStim2=intStimNr;end
					indStim1 = cellSelect{intStim1};
					indStim2 = cellSelect{intStim2};
					vecUseTrialTypes = vecTrialTypes(indStim1 | indStim2);
					
					%select data
					matUseData = matData(vecSelectCells,indStim1 | indStim2);
					%decode
					%vecUsePriorDistribution = vecCounts([intStim1 intStim2]);
					%[dblPerformanceCV,vecDecodedIndexCV,matTemplateDistsCV,dblMeanErrorDegs,matConfusion] = doCrossValidatedDecodingTM(matUseData,vecUseTrialTypes,intTypeCV,vecUsePriorDistribution);
					%matPerfTM(intArea,intStim1) = dblPerformanceCV;
					[dblPerformanceLR,vecDecodedIndexCV,matPosteriorProbability,dblMeanErrorDegs,matConfusion,matWeights] = ...
						doCrossValidatedDecodingLR(matUseData,vecUseTrialTypes,intTypeCV,vecPriorDistribution([intStim1 intStim2]),dblLambda);
					matPerf(intArea,intStimCenter) = dblPerformanceLR;
					
				end
			end
			%maxfig;drawnow;
			%export_fig(fullpath(strTargetPath,['GratingDecoding_' strName 'B' num2str(intBlock) '.jpg']));
			%export_fig(fullpath(strTargetPath,['GratingDecoding_' strName 'B' num2str(intBlock) '.pdf']));
			%close;
			if strcmp(strSubjectType,'BL6')
				vecAggPerfWt(end+1,:,:) = vecPerf;
				matAggPerfWt(end+1,:,:) = matPerf;
				cellAggPrefWt(end+1,:) = cellPref;
			else
				vecAggPerfAlb(end+1,:,:) = vecPerf;
				matAggPerfAlb(end+1,:,:) = matPerf;
				cellAggPrefAlb(end+1,:) = cellPref;
			end
			
		end
	end
	
	%% plot
	
end
%% plot

vecPlotDegs1 = vecUnique;
vecPlotDegs1(vecPlotDegs1>=180) = vecPlotDegs1(vecPlotDegs1>=180) - 360;
[vecPlotDegs1,vecReorder1]=sort(vecPlotDegs1);
vecPlotDegs = vecUnique;
vecPlotDegs(vecPlotDegs>=180) = vecPlotDegs(vecPlotDegs>=180) - 360;
[vecPlotDegs,vecReorder]=sort(vecPlotDegs);
dblEndX = vecPlotDegs(end)+median(diff(vecPlotDegs));
vecPlotDegs = [vecPlotDegs; dblEndX];

%extract data
matWtPerfCtx = squeeze(matAggPerfWt(:,1,vecReorder));
matWtPerfCtx(any(isnan(matWtPerfCtx),2),:) = [];
matWtPerfNot = squeeze(matAggPerfWt(:,2,vecReorder));
matWtPerfNot(any(isnan(matWtPerfNot),2),:) = [];

matAlbPerfCtx = squeeze(matAggPerfAlb(:,1,vecReorder));
matAlbPerfCtx(any(isnan(matAlbPerfCtx),2),:) = [];
matAlbPerfNot = squeeze(matAggPerfAlb(:,2,vecReorder));
matAlbPerfNot(any(isnan(matAlbPerfNot),2),:) = [];

vecMeanWtPerfCtx = mean(matWtPerfCtx,1);
vecMeanWtPerfNot = mean(matWtPerfNot,1);
vecMeanAlbPerfCtx = mean(matAlbPerfCtx,1);
vecMeanAlbPerfNot = mean(matAlbPerfNot,1);
vecSemWtPerfCtx = std(matWtPerfCtx,[],1)/sqrt(size(matWtPerfCtx,1));
vecSemWtPerfNot = std(matWtPerfNot,[],1)/sqrt(size(matWtPerfNot,1));
vecSemAlbPerfCtx = std(matAlbPerfCtx,[],1)/sqrt(size(matAlbPerfCtx,1));
vecSemAlbPerfNot = std(matAlbPerfNot,[],1)/sqrt(size(matAlbPerfNot,1));

% plot pair decoding
figure
subplot(2,3,1)
hold on
plot([-180 180],[0.5 0.5],'--','Color',[0.5 0.5 0.5]);
errorbar(vecPlotDegs,[vecMeanWtPerfCtx(:); vecMeanWtPerfCtx(1)],[vecSemWtPerfCtx(:); vecSemWtPerfCtx(1)],'x-','Color',vecColBl6)
errorbar(vecPlotDegs,[vecMeanAlbPerfCtx(:); vecMeanAlbPerfCtx(1)],[vecSemAlbPerfCtx(:); vecSemAlbPerfCtx(1)],'x-','Color',vecColAlb)
xlabel('Stim ori (degs)')
ylabel('Grating pair dec. perf.')
legend({'Chance','BL6','DBA'},'Location','best')
xlim([-180 180]);
set(gca,'xtick',min(vecPlotDegs):90:max(vecPlotDegs));
title(cellAreaGroups{1});
fixfig;

subplot(2,3,4)
hold on
plot([-180 180],[0.5 0.5],'--','Color',[0.5 0.5 0.5]);
errorbar(vecPlotDegs,[vecMeanWtPerfNot(:); vecMeanWtPerfNot(1)],[vecSemWtPerfNot(:); vecSemWtPerfNot(1)],'x-','Color',vecColBl6);
errorbar(vecPlotDegs,[vecMeanAlbPerfNot(:); vecMeanAlbPerfNot(1)],[vecSemAlbPerfNot(:); vecSemAlbPerfNot(1)],'x-','Color',vecColAlb);
xlabel('Stim ori (degs)')
ylabel('Grating pair dec. perf.')
legend({'Chance','BL6','DBA'},'Location','best')
xlim([-180 180]);
set(gca,'xtick',min(vecPlotDegs):90:max(vecPlotDegs));
title(cellAreaGroups{2});
fixfig;

%{
subplot(2,3,3)
hold on
plot([0 360],[0.5 0.5],'k--');
plot(vecPlot,[vecMeanWtPerfHip; vecMeanWtPerfHip(1)],'b-');
plot(vecPlot,[vecMeanAlbPerfHip; vecMeanAlbPerfHip(1)],'r-');
xlabel('Stim ori (degs)')
ylabel('Grating pair dec. perf.')
legend({'Chance','BL6','DBA'},'Location','best')
xlim([0 360]);
set(gca,'xtick',0:90:360);
title(cellAreaGroups{3});
fixfig;
%}
normaxes;

% plot pref ori
%%{
vecPrefOriCtxWt = cell2vec(cellAggPrefWt(:,1));
vecPrefOriNotWt = cell2vec(cellAggPrefWt(:,2));
%vecPrefOriHipWt = cell2vec(cellAggPrefWt(:,3));

vecPrefOriCtxAlb = cell2vec(cellAggPrefAlb(:,1));
vecPrefOriNotAlb = cell2vec(cellAggPrefAlb(:,2));
%vecPrefOriHipAlb = cell2vec(cellAggPrefAlb(:,3));
vecBinCenters = 0:15:359;
vecBins = (0:15:375) - 15/2;
vecX = vecBins(2:end)-median(diff(vecBins))/2;

vecPlotBins = vecX;
vecPlotBins(vecPlotBins>=180) = vecPlotBins(vecPlotBins>=180) - 360;
vecPlotBins(end) = 180;
[vecPlotBins,vecReorder2]=sort(vecPlotBins);

vecCountsCtxWt = histcounts(vecPrefOriCtxWt,vecBins);
vecCountsNotWt = histcounts(vecPrefOriNotWt,vecBins);
%vecCountsHipWt = histcounts(vecPrefOriHipWt,vecBins);
vecCountsCtxAlb = histcounts(vecPrefOriCtxAlb,vecBins);
vecCountsNotAlb = histcounts(vecPrefOriNotAlb,vecBins);
%vecCountsHipAlb = histcounts(vecPrefOriHipAlb,vecBins);

%smooth
vecCountsCtxWt = imfilt(vecCountsCtxWt(vecReorder2),[1 1 1]/3);
vecCountsNotWt = imfilt(vecCountsNotWt(vecReorder2),[1 1 1]/3);
vecCountsCtxAlb = imfilt(vecCountsCtxAlb(vecReorder2),[1 1 1]/3);
vecCountsNotAlb = imfilt(vecCountsNotAlb(vecReorder2),[1 1 1]/3);

subplot(2,3,2)
hold on;
plot(vecPlotBins,vecCountsCtxWt/sum(vecCountsCtxWt),'x-','Color',vecColBl6)
plot(vecPlotBins,vecCountsCtxAlb/sum(vecCountsCtxAlb),'x-','Color',vecColAlb)
hold off
title(cellAreaGroups{1});
ylabel('Normalized # of neurons');
xlabel('Preferred orientation (degs)');
xlim([min(vecPlotBins) max(vecPlotBins)]);
set(gca,'xtick',min(vecPlotBins):90:max(vecPlotBins));
fixfig;

subplot(2,3,5)
hold on;
plot(vecPlotBins,vecCountsNotWt/sum(vecCountsNotWt),'x-','Color',vecColBl6)
plot(vecPlotBins,vecCountsNotAlb/sum(vecCountsNotAlb),'x-','Color',vecColAlb)
title(cellAreaGroups{2});
ylabel('Normalized # of neurons');
xlabel('Preferred orientation (degs)');
xlim([min(vecPlotBins) max(vecPlotBins)]);
set(gca,'xtick',min(vecPlotBins):90:max(vecPlotBins));
fixfig;

%subplot(2,3,6)
%hold on;
%plot(vecX,vecCountsHipWt,'b-')
%plot(vecX,vecCountsHipAlb,'r-')
%%}

%plot pop resp
vecMeanRespCtxWt = mean(cellAggMeanRespWt{1},1);
vecMeanRespNotWt = mean(cellAggMeanRespWt{2},1);
vecSdRespCtxWt = std(cellAggMeanRespWt{1},[],1);
vecSdRespNotWt = std(cellAggMeanRespWt{2},[],1);

vecMeanRespCtxAlb = mean(cellAggMeanRespAlb{1},1);
vecMeanRespNotAlb = mean(cellAggMeanRespAlb{2},1);
vecSdRespCtxAlb = std(cellAggMeanRespAlb{1},[],1);
vecSdRespNotAlb = std(cellAggMeanRespAlb{2},[],1);

%smooth
vecMeanRespCtxWt = imfilt(vecMeanRespCtxWt(vecReorder),[1 1 1]/3);
vecMeanRespNotWt = imfilt(vecMeanRespNotWt(vecReorder),[1 1 1]/3);
vecSdRespCtxWt = imfilt(vecSdRespCtxWt(vecReorder),[1 1 1]/3);
vecSdRespNotWt = imfilt(vecSdRespNotWt(vecReorder),[1 1 1]/3);

vecMeanRespCtxAlb = imfilt(vecMeanRespCtxAlb(vecReorder),[1 1 1]/3);
vecMeanRespNotAlb = imfilt(vecMeanRespNotAlb(vecReorder),[1 1 1]/3);
vecSdRespCtxAlb = imfilt(vecSdRespCtxAlb(vecReorder),[1 1 1]/3);
vecSdRespNotAlb = imfilt(vecSdRespNotAlb(vecReorder),[1 1 1]/3);

subplot(2,3,3)
hold on;
errorbar(vecPlotDegs1,vecMeanRespCtxWt,vecSdRespCtxWt/sqrt(size(cellAggMeanRespWt{1},1)),'x-','Color',vecColBl6)
errorbar(vecPlotDegs1,vecMeanRespCtxAlb,vecSdRespCtxAlb/sqrt(size(cellAggMeanRespAlb{1},1)),'x-','Color',vecColAlb)
xlabel('Stim ori (degs)')
ylabel('Normalized mean tuning curve')
legend({'BL6','DBA'},'Location','best')
xlim([-180 180]);
set(gca,'xtick',-180:90:180);
title(cellAreaGroups{1});
fixfig;

subplot(2,3,6)
hold on;
errorbar(vecPlotDegs1,vecMeanRespNotWt,vecSdRespNotWt/sqrt(size(cellAggMeanRespWt{2},1)),'x-','Color',vecColBl6)
errorbar(vecPlotDegs1,vecMeanRespNotAlb,vecSdRespNotAlb/sqrt(size(cellAggMeanRespAlb{2},1)),'x-','Color',vecColAlb)
xlabel('Stim ori (degs)')
ylabel('Normalized mean tuning curve')
legend({'BL6','DBA'},'Location','best')
xlim([-180 180]);
set(gca,'xtick',-180:90:180);
title(cellAreaGroups{2});
fixfig;
%{
subplot(2,3,6)
hold on;
errorbar(vecUnique,vecMeanRespHipWt,vecSdRespHipWt/sqrt(size(cellAggMeanRespWt{3},1)),'b-')
errorbar(vecUnique,vecMeanRespHipAlb,vecSdRespHipAlb/sqrt(size(cellAggMeanRespAlb{3},1)),'r-')
xlabel('Stim ori (degs)')
ylabel('Mean tuning curve')
legend({'BL6','DBA'},'Location','best')
xlim([0 360]);
set(gca,'xtick',0:90:360);
title(cellAreaGroups{3});
fixfig;
%}
drawnow;maxfig;
%% save
drawnow;
export_fig(fullpath(strTargetPath,['OriDecodingAndPreference.tif']));
export_fig(fullpath(strTargetPath,['OriDecodingAndPreference.pdf']));

%% test horizontal vs vertical
%define horizontal and vertical
figure;maxfig;
dblSurroundDegs = 30;
indVert = mod(vecUnique,180) <= 0+dblSurroundDegs | mod(vecUnique,180) >= 180-dblSurroundDegs;
indHorz = mod(vecUnique,180) >= 90-dblSurroundDegs & mod(vecUnique,180) <= 90+dblSurroundDegs;
vecVertOris = vecUnique(indVert);
vecHorzOris = vecUnique(indHorz);

%extract data
vecWtPerfCtx_Vert = mean(matWtPerfCtx(:,indVert),2);
vecWtPerfCtx_Horz = mean(matWtPerfCtx(:,indHorz),2);
vecAlbPerfCtx_Vert = mean(matAlbPerfCtx(:,indVert),2);
vecAlbPerfCtx_Horz = mean(matAlbPerfCtx(:,indHorz),2);

vecWtPerfNot_Vert = mean(matWtPerfNot(:,indVert),2);
vecWtPerfNot_Horz = mean(matWtPerfNot(:,indHorz),2);
vecAlbPerfNot_Vert = mean(matAlbPerfNot(:,indVert),2);
vecAlbPerfNot_Horz = mean(matAlbPerfNot(:,indHorz),2);

subplot(2,3,1)
hold on
errorbar([1 2]+0.05,[mean(vecWtPerfCtx_Horz) mean(vecWtPerfCtx_Vert)],[std(vecWtPerfCtx_Horz) std(vecWtPerfCtx_Vert)]./sqrt(numel(vecWtPerfCtx_Vert)),...
	'x-','Color',vecColBl6,'CapSize',20);
errorbar([1 2]-0.05,[mean(vecAlbPerfCtx_Horz) mean(vecAlbPerfCtx_Vert)],[std(vecAlbPerfCtx_Horz) std(vecAlbPerfCtx_Vert)]./sqrt(numel(vecAlbPerfCtx_Vert)),...
	'x-','Color',vecColAlb,'CapSize',20);
plot([0.5 2.5],[1 1]*0.5,'--','Color',[0.5 0.5 0.5]);
hold off
ylabel('Decoding accuracy');
set(gca,'xtick',[1 2],'xticklabel',{'Horz','Vert'});
xlim([0.5 2.5]);
ylim([0.45 0.75]);
title(cellAreaGroupsAbbr{1});
legend(gca,{'BL6','DBA','Chance'},'location','best');
fixfig;


subplot(2,3,2)
hold on
errorbar([1 2]+0.05,[mean(vecWtPerfNot_Horz) mean(vecWtPerfNot_Vert)],[std(vecWtPerfNot_Horz) std(vecWtPerfNot_Vert)]./sqrt(numel(vecWtPerfNot_Vert)),...
	'x-','Color',vecColBl6,'CapSize',20);
errorbar([1 2]-0.05,[mean(vecAlbPerfNot_Horz) mean(vecAlbPerfNot_Vert)],[std(vecAlbPerfNot_Horz) std(vecAlbPerfNot_Vert)]./sqrt(numel(vecAlbPerfNot_Vert)),...
	'x-','Color',vecColAlb,'CapSize',20);
plot([0.5 2.5],[1 1]*0.5,'--','Color',[0.5 0.5 0.5]);
hold off
ylabel('Decoding accuracy');
set(gca,'xtick',[1 2],'xticklabel',{'Horz','Vert'});
xlim([0.5 2.5]);
ylim([0.45 0.75]);
title(cellAreaGroupsAbbr{2});
legend(gca,{'BL6','DBA','Chance'},'location','best');
fixfig;

%calc diff
vecWtCtxDiff = vecWtPerfNot_Horz - vecWtPerfNot_Vert;
vecAlbCtxDiff = vecAlbPerfNot_Horz - vecAlbPerfNot_Vert;
vecWtNotDiff = vecWtPerfCtx_Horz - vecWtPerfCtx_Vert;
vecAlbNotDiff = vecAlbPerfCtx_Horz - vecAlbPerfCtx_Vert;

[h,pWCD0]=ttest(vecWtCtxDiff);
[h,pACD0]=ttest(vecAlbCtxDiff);
[h,pWND0]=ttest(vecWtNotDiff);
[h,pAND0]=ttest(vecAlbNotDiff);

[h,pCtxWvA]=ttest2(vecWtCtxDiff,vecAlbCtxDiff);
[h,pNotWvA]=ttest2(vecWtNotDiff,vecAlbNotDiff);

%plot
subplot(2,3,4)
x1=beeswarm(ones(size(vecWtCtxDiff)),vecWtCtxDiff);
x2=beeswarm(2*ones(size(vecAlbCtxDiff)),vecAlbCtxDiff);
cla;
subplot(2,3,4)
hold on
scatter(x1,vecWtCtxDiff,[],vecColBl6,'o');
scatter(x2,vecAlbCtxDiff,[],vecColAlb,'o');
plot([0.5 2.5],[0 0],'--','color',[0.5 0.5 0.5]);
hold off
xlim([0.5 2.5]);
ylim([-1.001 1.001]*max(abs(get(gca,'ylim'))));
ylabel('Horz - Vert difference Dec. perf.');
set(gca,'xtick',[1 2],'xticklabel',{'BL6','DBA'});
title(sprintf('%s: t v 0;BL6,p=%.3f; DBA,p=%.3f; 2t,p=%.3f',cellAreaGroupsAbbr{1},pWCD0,pACD0,pCtxWvA));
fixfig;grid off


subplot(2,3,5)
x1=beeswarm(ones(size(vecWtNotDiff)),vecWtNotDiff);
x2=beeswarm(2*ones(size(vecAlbNotDiff)),vecAlbNotDiff);
cla;
subplot(2,3,5)
hold on
scatter(x1,vecWtNotDiff,[],vecColBl6,'o');
scatter(x2,vecAlbNotDiff,[],vecColAlb,'o');
plot([0.5 2.5],[0 0],'--','color',[0.5 0.5 0.5]);
hold off
xlim([0.5 2.5]);
ylim([-1.001 1.001]*max(abs(get(gca,'ylim'))));
ylabel('Horz - Vert difference Dec. perf.');
set(gca,'xtick',[1 2],'xticklabel',{'BL6','DBA'});
title(sprintf('%s: t v 0;BL6,p=%.3f; DBA,p=%.3f; 2t,p=%.3f',cellAreaGroupsAbbr{2},pWND0,pAND0,pNotWvA));
fixfig;grid off
drawnow;

%% decoding performance of anti-direction vs chance for bl6 and albino
matConfCtxWt = squeeze(matAggConfusionWt(:,1,:,:));
indRemCtxWt = sum(sum(matConfCtxWt,2),3)==0;
matConfCtxWt(indRemCtxWt,:,:)=[];
matConfNotWt = squeeze(matAggConfusionWt(:,2,:,:));
indRemNotWt = sum(sum(matConfNotWt,2),3)==0;
matConfNotWt(indRemNotWt,:,:)=[];

matConfCtxAlb = squeeze(matAggConfusionAlb(:,1,:,:));
indRemCtxAlb = sum(sum(matConfCtxAlb,2),3)==0;
matConfCtxAlb(indRemCtxAlb,:,:)=[];
matConfNotAlb = squeeze(matAggConfusionAlb(:,2,:,:));
indRemNotAlb = sum(sum(matConfNotAlb,2),3)==0;
matConfNotAlb(indRemNotAlb,:,:)=[];

matProDir = diag(diag(true(intStimNr,intStimNr)));
matAntiDir = circshift(matProDir,intStimNr/2,1);

%CtxWt
intRecCtxWt = size(matConfCtxWt,1);
vecAntiDirP_CtxWt = nan(1,intRecCtxWt);
for intRec=1:intRecCtxWt
	matThisC = squeeze(matConfCtxWt(intRec,:,:));
	intAllTrials = sum(matThisC(:));
	if intAllTrials==0,continue;end
	vecAntiDirP_CtxWt(intRec) = sum(matThisC(matAntiDir))/(intAllTrials-sum(matThisC(matProDir)));
end
%NotWt
intRecNotWt = size(matConfNotWt,1);
vecAntiDirP_NotWt = nan(1,intRecNotWt);
for intRec=1:intRecNotWt
	matThisC = squeeze(matConfNotWt(intRec,:,:));
	intAllTrials = sum(matThisC(:));
	if intAllTrials==0,continue;end
	vecAntiDirP_NotWt(intRec) = sum(matThisC(matAntiDir))/(intAllTrials-sum(matThisC(matProDir)));
end
%CtxAlb
intRecCtxAlb = size(matConfCtxAlb,1);
vecAntiDirP_CtxAlb = nan(1,intRecCtxAlb);
for intRec=1:intRecCtxAlb
	matThisC = squeeze(matConfCtxAlb(intRec,:,:));
	intAllTrials = sum(matThisC(:));
	if intAllTrials==0,continue;end
	vecAntiDirP_CtxAlb(intRec) = sum(matThisC(matAntiDir))/(intAllTrials-sum(matThisC(matProDir)));
end
%NotAlb
intRecNotAlb = size(matConfNotAlb,1);
vecAntiDirP_NotAlb = nan(1,intRecNotAlb);
for intRec=1:intRecNotAlb
	matThisC = squeeze(matConfNotAlb(intRec,:,:));
	intAllTrials = sum(matThisC(:));
	if intAllTrials==0,continue;end
	vecAntiDirP_NotAlb(intRec) = sum(matThisC(matAntiDir))/(intAllTrials-sum(matThisC(matProDir)));
end

%tests
dblChanceP = 1/intStimNr;
[h,pADCW]=ttest(vecAntiDirP_CtxWt,dblChanceP);
[h,pADCA]=ttest(vecAntiDirP_CtxAlb,dblChanceP);
[h,pADNW]=ttest(vecAntiDirP_NotWt,dblChanceP);
[h,pADNA]=ttest(vecAntiDirP_NotAlb,dblChanceP);
[h,pAD2W]=ttest2(vecAntiDirP_CtxWt,vecAntiDirP_NotWt);
[h,pAD2A]=ttest2(vecAntiDirP_CtxAlb,vecAntiDirP_NotAlb);

[pAD2W_ranksum] = ranksum(vecAntiDirP_CtxWt,vecAntiDirP_NotWt,'method','exact');
[pAD2A_ranksum] = ranksum(vecAntiDirP_CtxAlb,vecAntiDirP_NotAlb,'method','exact');

% plot
%Ctx
subplot(2,3,3)
x1=beeswarm(ones(size(vecAntiDirP_CtxWt)),vecAntiDirP_CtxWt);
x2=beeswarm(2*ones(size(vecAntiDirP_NotWt)),vecAntiDirP_NotWt);
cla;
subplot(2,3,3)
hold on
scatter(x1,vecAntiDirP_CtxWt,[],vecColBl6,'o');
scatter(x2,vecAntiDirP_NotWt,[],vecColBl6,'o');
plot([0.5 2.5],dblChanceP*[1 1],'--','color',[0.5 0.5 0.5]);
hold off
xlim([0.5 2.5]);
ylim([0 0.16]);
ylabel('P(Anti-direction|Error)');
set(gca,'xtick',[1 2],'xticklabel',cellAreaGroupsAbbr);
title(sprintf('%s: t v rand;Ctx,p=%.3f; NOT,p=%.3f; 2t,p=%.3f','BL6',pADCW,pADNW,pAD2W));
fixfig;grid off
drawnow;

%Not
subplot(2,3,6)
x1=beeswarm(ones(size(vecAntiDirP_CtxAlb)),vecAntiDirP_CtxAlb);
x2=beeswarm(2*ones(size(vecAntiDirP_NotAlb)),vecAntiDirP_NotAlb);
cla;
subplot(2,3,6)
hold on
scatter(x1,vecAntiDirP_CtxAlb,[],vecColAlb,'o');
scatter(x2,vecAntiDirP_NotAlb,[],vecColAlb,'o');
plot([0.5 2.5],dblChanceP*[1 1],'--','color',[0.5 0.5 0.5]);
hold off
xlim([0.5 2.5]);
ylim([0 0.16]);
ylabel('P(Anti-direction|Error)');
set(gca,'xtick',[1 2],'xticklabel',cellAreaGroupsAbbr);
title(sprintf('%s: t v rand;Ctx,p=%.3f; NOT,p=%.3f; 2t,p=%.3f','DBA',pADCA,pADNA,pAD2A));
fixfig;grid off
drawnow;

%% save
drawnow;
export_fig(fullpath(strTargetPath,['OriDecodingHorzVert.tif']));
export_fig(fullpath(strTargetPath,['OriDecodingHorzVert.pdf']));

%% confusion matrices
matConfCtxWt = squeeze(mean(matAggConfusionWt(:,1,:,:),1));
matConfNotWt = squeeze(mean(matAggConfusionWt(:,2,:,:),1));
%matConfHipWt = squeeze(mean(matAggConfusionWt(:,3,:,:),1));
matConfCtxAlb = squeeze(mean(matAggConfusionAlb(:,1,:,:),1));
matConfNotAlb = squeeze(mean(matAggConfusionAlb(:,2,:,:),1));
%matConfHipAlb = squeeze(mean(matAggConfusionAlb(:,3,:,:),1));

vecTick = (1:6:25)-0.5;
vecTickL = 0:90:360;
figure;maxfig
subplot(2,3,1)
imagesc(matConfCtxWt);
title(['BL6 ' cellAreaGroups{1}]);
set(gca,'xtick',vecTick,'xticklabel',vecTickL);
set(gca,'ytick',vecTick,'yticklabel',vecTickL);
axis xy
fixfig;
grid off;
xlabel('Decoded orientation');
ylabel('Stimulus orientation');

subplot(2,3,2)
imagesc(matConfCtxAlb);
title(['Alb ' cellAreaGroups{1}]);
set(gca,'xtick',vecTick,'xticklabel',vecTickL);
set(gca,'ytick',vecTick,'yticklabel',vecTickL);
axis xy
fixfig;
grid off;
xlabel('Decoded orientation');
ylabel('Stimulus orientation');

%subplot(2,3,3)
%imagesc(matConfHipWt);
%title(['BL6 ' cellAreaGroups{3}]);


subplot(2,3,4)
imagesc(matConfNotWt);
title(['BL6 ' cellAreaGroups{2}]);
set(gca,'xtick',vecTick,'xticklabel',vecTickL);
set(gca,'ytick',vecTick,'yticklabel',vecTickL);
axis xy
fixfig;
grid off;
xlabel('Decoded orientation');
ylabel('Stimulus orientation');

subplot(2,3,5)
imagesc(matConfNotAlb);
title(['Alb ' cellAreaGroups{2}]);
set(gca,'xtick',vecTick,'xticklabel',vecTickL);
set(gca,'ytick',vecTick,'yticklabel',vecTickL);
axis xy
fixfig;
grid off;
xlabel('Decoded orientation');
ylabel('Stimulus orientation');

%subplot(2,3,6)
%imagesc(matConfHipAlb);
%title(['Alb ' cellAreaGroups{3}]);

%% save
drawnow;
export_fig(fullpath(strTargetPath,['OriDecodingAveraged.tif']));
export_fig(fullpath(strTargetPath,['OriDecodingAveraged.pdf']));