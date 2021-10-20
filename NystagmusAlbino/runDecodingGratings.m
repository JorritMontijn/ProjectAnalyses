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
cellUseAreas{1} = {'Primary visual','Posteromedial visual','anteromedial visual'};
%NOT
cellUseAreas{2} = {'nucleus of the optic tract'};
%hippocampus
cellUseAreas{3} = {'Hippocampal formation','Field CA1','Field CA2','Field CA3','subiculum','dentate gyrus'};
cellAreaGroups = {'Vis. ctx','NOT','Hippocampus'};
cellAreaGroupsAbbr = {'Ctx','NOT','Hip'};

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
			vecZetaP = cellfun(@min,{sRec.sCluster.ZetaP});
			%indUseCells = arrayfun(@(x) x.KilosortGood==1 | x.Contamination < 0.1,sRec.sCluster(:));
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
			figure
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
				if strcmp(strSubjectType,'BL6')
					cellAggMeanRespWt{intArea} = cat(1,cellAggMeanRespWt{intArea},sOut.matMeanResp(indKeepCells,:));
				else
					cellAggMeanRespAlb{intArea} = cat(1,cellAggMeanRespAlb{intArea},sOut.matMeanResp(indKeepCells,:));
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
				
				subplot(2,3,intArea)
				imagesc(matConfusionLR)
				title([strAreaGroup '; LR: ' strName '_' num2str(intBlock)],'interpreter','none');
				
				%% adjacent stims
				for intStim1=1:intStimNr
					intStim2 = intStim1-1;
					if intStim2==0
						intStim2=intStimNr;
					end
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
					matPerf(intArea,intStim1) = dblPerformanceLR;
					
				end
			end
			maxfig;drawnow;
			export_fig(fullpath(strTargetPath,['GratingDecoding_' strName 'B' num2str(intBlock) '.jpg']));
			export_fig(fullpath(strTargetPath,['GratingDecoding_' strName 'B' num2str(intBlock) '.pdf']));

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
%%
vecMeanWtPerfCtx = squeeze(nanmean(matAggPerfWt(:,1,:),1));
vecMeanWtPerfNot = squeeze(nanmean(matAggPerfWt(:,2,:),1));
vecMeanWtPerfHip = squeeze(nanmean(matAggPerfWt(:,3,:),1));

subplot(2,3,1)
polar(deg2rad(vecUnique),vecMeanWtPerfCtx-0.5)


subplot(2,3,2)
polar(deg2rad(vecUnique),vecMeanWtPerfNot-0.5)

subplot(2,3,3)
polar(deg2rad(vecUnique),vecMeanWtPerfHip-0.5)

vecMeanAlbPerfCtx = squeeze(nanmean(matAggPerfAlb(:,1,:),1));
vecMeanAlbPerfNot = squeeze(nanmean(matAggPerfAlb(:,2,:),1));
vecMeanAlbPerfHip = squeeze(nanmean(matAggPerfAlb(:,3,:),1));

subplot(2,3,4)
polar(deg2rad(vecUnique),vecMeanAlbPerfCtx-0.5)


subplot(2,3,5)
polar(deg2rad(vecUnique),vecMeanAlbPerfNot-0.5)

subplot(2,3,6)
polar(deg2rad(vecUnique),vecMeanAlbPerfHip-0.5)

%%
dblEndX = vecUnique(end)+median(diff(vecUnique));
figure
subplot(2,3,1)
hold on
plot([0 360],[0.5 0.5],'k--');
plot([vecUnique; dblEndX],[vecMeanWtPerfCtx; vecMeanWtPerfCtx(1)],'b-');
plot([vecUnique; dblEndX],[vecMeanAlbPerfCtx; vecMeanAlbPerfCtx(1)],'r-');
xlabel('Stim ori (degs)')
ylabel('Grating pair dec. perf.')
legend({'Chance','BL6','DBA'},'Location','best')
xlim([0 360]);
set(gca,'xtick',0:90:360);
title(cellAreaGroups{1});
fixfig;

subplot(2,3,2)
hold on
plot([0 360],[0.5 0.5],'k--');
plot([vecUnique; dblEndX],[vecMeanWtPerfNot; vecMeanWtPerfNot(1)],'b-');
plot([vecUnique; dblEndX],[vecMeanAlbPerfNot; vecMeanAlbPerfNot(1)],'r-');
xlabel('Stim ori (degs)')
ylabel('Grating pair dec. perf.')
legend({'Chance','BL6','DBA'},'Location','best')
xlim([0 360]);
set(gca,'xtick',0:90:360);
title(cellAreaGroups{2});
fixfig;


subplot(2,3,3)
hold on
plot([0 360],[0.5 0.5],'k--');
plot([vecUnique; dblEndX],[vecMeanWtPerfHip; vecMeanWtPerfHip(1)],'b-');
plot([vecUnique; dblEndX],[vecMeanAlbPerfHip; vecMeanAlbPerfHip(1)],'r-');
xlabel('Stim ori (degs)')
ylabel('Grating pair dec. perf.')
legend({'Chance','BL6','DBA'},'Location','best')
xlim([0 360]);
set(gca,'xtick',0:90:360);
title(cellAreaGroups{3});
fixfig;

normaxes;

%% plot pref ori
%{
vecPrefOriCtxWt = cell2vec(cellAggPrefWt(:,1));
vecPrefOriNotWt = cell2vec(cellAggPrefWt(:,2));
vecPrefOriHipWt = cell2vec(cellAggPrefWt(:,3));

vecPrefOriCtxAlb = cell2vec(cellAggPrefAlb(:,1));
vecPrefOriNotAlb = cell2vec(cellAggPrefAlb(:,2));
vecPrefOriHipAlb = cell2vec(cellAggPrefAlb(:,3));
vecBinCenters = 0:15:359;
vecBins = (0:15:375) - 15/2;

vecCountsCtxWt = histcounts(vecPrefOriCtxWt,vecBins);
vecCountsNotWt = histcounts(vecPrefOriNotWt,vecBins);
vecCountsHipWt = histcounts(vecPrefOriHipWt,vecBins);
vecCountsCtxAlb = histcounts(vecPrefOriCtxAlb,vecBins);
vecCountsNotAlb = histcounts(vecPrefOriNotAlb,vecBins);
vecCountsHipAlb = histcounts(vecPrefOriHipAlb,vecBins);

vecX = vecBins(2:end)-median(diff(vecBins));
subplot(2,3,4)
hold on;
plot(vecX,vecCountsCtxWt,'b-')
plot(vecX,vecCountsCtxAlb,'r-')

subplot(2,3,5)
hold on;
plot(vecX,vecCountsNotWt,'b-')
plot(vecX,vecCountsNotAlb,'r-')

subplot(2,3,6)
hold on;
plot(vecX,vecCountsHipWt,'b-')
plot(vecX,vecCountsHipAlb,'r-')
%}
%% plot pop resp
%{
vecMeanRespCtxWt = mean(bsxfun(@rdivide,cellAggMeanRespWt{1},sum(cellAggMeanRespWt{1},2)),1);
vecMeanRespNotWt = mean(bsxfun(@rdivide,cellAggMeanRespWt{2},sum(cellAggMeanRespWt{2},2)),1);
vecMeanRespHipWt = mean(bsxfun(@rdivide,cellAggMeanRespWt{3},sum(cellAggMeanRespWt{3},2)),1);
vecSdRespCtxWt = std(bsxfun(@rdivide,cellAggMeanRespWt{1},sum(cellAggMeanRespWt{1},2)),[],1);
vecSdRespNotWt = std(bsxfun(@rdivide,cellAggMeanRespWt{2},sum(cellAggMeanRespWt{2},2)),[],1);
vecSdRespHipWt = std(bsxfun(@rdivide,cellAggMeanRespWt{3},sum(cellAggMeanRespWt{3},2)),[],1);

vecMeanRespCtxAlb = mean(bsxfun(@rdivide,cellAggMeanRespAlb{1},sum(cellAggMeanRespAlb{1},2)),1);
vecMeanRespNotAlb = mean(bsxfun(@rdivide,cellAggMeanRespAlb{2},sum(cellAggMeanRespAlb{2},2)),1);
vecMeanRespHipAlb = mean(bsxfun(@rdivide,cellAggMeanRespAlb{3},sum(cellAggMeanRespAlb{3},2)),1);
vecSdRespCtxAlb = std(bsxfun(@rdivide,cellAggMeanRespAlb{1},sum(cellAggMeanRespAlb{1},2)),[],1);
vecSdRespNotAlb = std(bsxfun(@rdivide,cellAggMeanRespAlb{2},sum(cellAggMeanRespAlb{2},2)),[],1);
vecSdRespHipAlb = std(bsxfun(@rdivide,cellAggMeanRespAlb{3},sum(cellAggMeanRespAlb{3},2)),[],1);
%}
vecMeanRespCtxWt = mean(cellAggMeanRespWt{1},1);
vecMeanRespNotWt = mean(cellAggMeanRespWt{2},1);
vecMeanRespHipWt = mean(cellAggMeanRespWt{3},1);
vecSdRespCtxWt = std(cellAggMeanRespWt{1},[],1);
vecSdRespNotWt = std(cellAggMeanRespWt{2},[],1);
vecSdRespHipWt = std(cellAggMeanRespWt{3},[],1);

vecMeanRespCtxAlb = mean(cellAggMeanRespAlb{1},1);
vecMeanRespNotAlb = mean(cellAggMeanRespAlb{2},1);
vecMeanRespHipAlb = mean(cellAggMeanRespAlb{3},1);
vecSdRespCtxAlb = std(cellAggMeanRespAlb{1},[],1);
vecSdRespNotAlb = std(cellAggMeanRespAlb{2},[],1);
vecSdRespHipAlb = std(cellAggMeanRespAlb{3},[],1);

subplot(2,3,4)
hold on;
errorbar(vecUnique,vecMeanRespCtxWt,vecSdRespCtxWt/sqrt(size(cellAggMeanRespWt{1},1)),'b-')
errorbar(vecUnique,vecMeanRespCtxAlb,vecSdRespCtxAlb/sqrt(size(cellAggMeanRespAlb{1},1)),'r-')
xlabel('Stim ori (degs)')
ylabel('Mean tuning curve')
legend({'BL6','DBA'},'Location','best')
xlim([0 360]);
set(gca,'xtick',0:90:360);
title(cellAreaGroups{1});
fixfig;

subplot(2,3,5)
hold on;
errorbar(vecUnique,vecMeanRespNotWt,vecSdRespNotWt/sqrt(size(cellAggMeanRespWt{2},1)),'b-')
errorbar(vecUnique,vecMeanRespNotAlb,vecSdRespNotAlb/sqrt(size(cellAggMeanRespAlb{2},1)),'r-')
xlabel('Stim ori (degs)')
ylabel('Mean tuning curve')
legend({'BL6','DBA'},'Location','best')
xlim([0 360]);
set(gca,'xtick',0:90:360);
title(cellAreaGroups{2});
fixfig;

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

%% confusion
matConfCtxWt = squeeze(mean(matAggConfusionWt(:,1,:,:),1));
matConfNotWt = squeeze(mean(matAggConfusionWt(:,2,:,:),1));
matConfHipWt = squeeze(mean(matAggConfusionWt(:,3,:,:),1));
matConfCtxAlb = squeeze(mean(matAggConfusionAlb(:,1,:,:),1));
matConfNotAlb = squeeze(mean(matAggConfusionAlb(:,2,:,:),1));
matConfHipAlb = squeeze(mean(matAggConfusionAlb(:,3,:,:),1));

figure
subplot(2,3,1)
imagesc(matConfCtxWt);
title(['BL6 ' cellAreaGroups{1}]);
subplot(2,3,2)
imagesc(matConfNotWt);
title(['BL6 ' cellAreaGroups{2}]);
subplot(2,3,3)
imagesc(matConfHipWt);
title(['BL6 ' cellAreaGroups{3}]);
subplot(2,3,4)
imagesc(matConfCtxAlb);
title(['Alb ' cellAreaGroups{1}]);
subplot(2,3,5)
imagesc(matConfNotAlb);
title(['Alb ' cellAreaGroups{2}]);
subplot(2,3,6)
imagesc(matConfHipAlb);
title(['Alb ' cellAreaGroups{3}]);
return
%% save
drawnow;maxfig;
export_fig(fullpath(strTargetPath,['GratingTracking' getDate '.jpg']));
export_fig(fullpath(strTargetPath,['GratingTracking' getDate '.pdf']));