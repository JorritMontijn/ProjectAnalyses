%% creates pupil plot fig over time of grating
%[done/ii) spike-triggered average of eye movement in NOT
%[done/iv) plot pupil x position as function of time of grating

%% load data
clear all;
strDataPath = 'E:\DataPreProcessed';
%strDataPath = 'F:\Data\Processed\Neuropixels';
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

%% get data
%pre-allocate
vecPupilSTA_edges = [-0.5:0.02:0.5];
cellAggX = {};
cellAggY = {};
cellAggSize = {};
cellAggSync = {};
cellAggOriIdx = {};
cellSTA = cellfill(nan(0,numel(vecPupilSTA_edges)-1),[2 2])

%run
cellNameAP = arrayfun(@(x) x.sJson.file_preproAP,sExp,'uniformoutput',false);
cellExperiment = arrayfun(@(x) x.sJson.experiment,sExp,'uniformoutput',false);
cellRemove = {};%{'RecMA5_2021-03-01R01_g0_t0'};
indRemRecs = contains(cellExperiment,cellRemove);
indRemRecs2 = ~contains(cellNameAP,cellUseForEyeTracking);
cellSubjectType = arrayfun(@(x) x.sJson.subjecttype,sExp,'uniformoutput',false);
for intSubType=1:2
	if intSubType == 1
		strSubjectType = 'BL6';
		dblOffsetT=0;
		dblAverageMouseHeadTiltInSetup = -30; %average tilt due to head-bar placement; procedures changed between BL6 and DBA experiments
		boolInvertX = 1; %other eye was recorded, so temporonasal is nasotemporal
	elseif intSubType == 2
		strSubjectType = 'DBA';
		dblOffsetT=0;
		dblAverageMouseHeadTiltInSetup = 0;
		boolInvertX = 0;
	end
	indUseRecs = contains(cellSubjectType,strSubjectType);
	vecRunRecs = find(indUseRecs & ~(indRemRecs | indRemRecs2));
	matAggTE_LocX = [];
	matAggTE_LocY = [];
	matAggTE_Size = [];
	matAggTE_Sync = [];
	vecTrialCounts = [];
	matFracMoveRight = zeros(24,0);
	matFracMoveUp = zeros(24,0);
	vecAggOriIdx = [];
	vecAggCounts = zeros(24,1);
	for intRecIdx=1:numel(vecRunRecs)
		intRec=vecRunRecs(intRecIdx)
		sRec = sExp(intRec);
		
		%% pupil loc
		vecPupilT = sRec.sPupil.vecTime;
		vecPupilLocX = sRec.sPupil.vecCenterX;
		vecPupilLocY = sRec.sPupil.vecCenterY;
		vecPupilSize = sRec.sPupil.vecRadius;
		vecPupilSync = sRec.sPupil.vecSyncLum;
		if isfield(sRec.sPupil,'vecBlinks') && ~all(sRec.sPupil.vecBlinks==0)
			vecPupilBlinks = sRec.sPupil.vecBlinks;
		else
			%filter absvidlum
			dblLowPass = 0.01/(1/median(diff(vecPupilT)));
			[fb,fa] = butter(2,dblLowPass,'high');
			vecPupilAbsVidLum = zscore(filtfilt(fb,fa, sRec.sPupil.sRaw.vecPupilAbsVidLum));
			vecPupilBlinks = vecPupilAbsVidLum > 5;
		end
		dblFs = 1/median(diff(vecPupilT));
		
		%% get DG blocks
		cellStimType = cellfun(@(x) x.strExpType,sRec.cellBlock,'uniformoutput',false);
		vecBlocksDG = find(contains(cellStimType,'driftinggrating','IgnoreCase',true));
		for intBlockIdx=1:numel(vecBlocksDG)
			intBlock = vecBlocksDG(intBlockIdx);
			sBlock = sRec.cellBlock{intBlock};
			if isfield(sBlock,'vecPupilStimOn')
				vecPupilStimOn = sBlock.vecPupilStimOn;
			else
				vecPupilStimOn = sBlock.vecStimOnTime;
			end
			
			%split by ori
			sTrialObjects = sBlock.sStimObject(sBlock.vecTrialStimTypes);
			vecOrientation = cell2vec({sTrialObjects.Orientation});
			[vecOriIdx,vecUnique,vecCounts,cellSelect,vecRepetition] = val2idx(vecOrientation);
			if numel(vecUnique) ~= 24,continue,end
			
			%% are they tracking the stimulus? i.e., is more time spent moving in the direction of the stimulus?
			%[matTE,vecWindowBinCenters] = getRespMat(vecPupilT,vecPupilLocX,vecStimOn,[-0.2 1.3]);
			vecWindowBinEdges = 0.05:0.02:1;
			[matTE_LocX] = getRespMat(vecPupilT,vecPupilLocX,vecPupilStimOn,vecWindowBinEdges);
			[matTE_LocY] = getRespMat(vecPupilT,vecPupilLocY,vecPupilStimOn,vecWindowBinEdges);
			[matTE_Blink] = getRespMat(vecPupilT,double(vecPupilBlinks),vecPupilStimOn,vecWindowBinEdges);
			matTE_Blink(isnan(matTE_Blink))=0;
			vecBlinkFrac = sum(matTE_Blink==1,2)./sum(~isnan(matTE_Blink),2);
			vecWindowBinEdgesSize = -0.2:0.02:1.3;
			[matTE_Size,vecWindowBinCenters] = getRespMat(vecPupilT,vecPupilSize,vecPupilStimOn,vecWindowBinEdgesSize);
			[matTE_Sync,vecWindowBinCenters] = getRespMat(vecPupilT,vecPupilSync,vecPupilStimOn,vecWindowBinEdgesSize);
			%remove trials with >1/5 of blinking
			indRemTrials = vecBlinkFrac>(1/5);
			matTE_LocX(indRemTrials,:) = [];
			matTE_LocY(indRemTrials,:) = [];
			matTE_Size(indRemTrials,:) = [];
			matTE_Sync(indRemTrials,:) = [];
			vecTheseOris = vecOrientation(~indRemTrials);
			vecTheseOris = mod(vecTheseOris+dblAverageMouseHeadTiltInSetup,360);
			%add data
			matAggTE_LocX = cat(1,matAggTE_LocX,matTE_LocX);
			matAggTE_LocY = cat(1,matAggTE_LocY,matTE_LocY);
			matAggTE_Size = cat(1,matAggTE_Size,matTE_Size);
			matAggTE_Sync = cat(1,matAggTE_Sync,matTE_Sync);
			[vecOriIdx,vecUnique,vecCounts,cellSelect,vecRepetition] = val2idx(vecTheseOris);
			vecAggOriIdx = cat(1,vecAggOriIdx,vecOriIdx);
			
			vecAggCounts = vecAggCounts+vecCounts;
			
		end
		
		%% spike-triggered average per area
		vecPupildT = vecPupilT(2:end)-median(diff(vecPupilT))/2;
		vecPupilMoveX_clean = diff(vecPupilLocX);
		vecPupilMoveX_clean(vecPupilBlinks>0) = nan;
		indUseCells = arrayfun(@(x) x.Violations1ms < 0.25 & abs(x.NonStationarity) < 0.25,sRec.sCluster(:));
		cellAreasPerCluster = {sRec.sCluster.Area};
		for intArea=1:numel(cellUseAreas)
			% select cells
			strAreaGroup =  cellAreaGroupsAbbr{intArea};
			vecSelectCells = find(indUseCells(:) & flat(contains(cellAreasPerCluster,cellUseAreas{intArea},'IgnoreCase',true)));
			if isempty(vecSelectCells)
				continue;
			end
			vecSpikeT = sort(cell2vec({sRec.sCluster(vecSelectCells).SpikeTimes}));
			%reduce spikes to max 100k
			intMaxSpikes = 100000; 
			if numel(vecSpikeT) > intMaxSpikes
				vecUseSpikes = unique(round(linspace(1,numel(vecSpikeT),intMaxSpikes)));
			else
				vecUseSpikes = 1:numel(vecSpikeT);
			end
			[matTE,vecPupilSTA_centers] = getRespMat(vecPupildT,vecPupilMoveX_clean,vecSpikeT(vecUseSpikes),vecPupilSTA_edges);
			vecSTA = nanmean(matTE,1);
			cellSTA{intSubType,intArea}(end+1,:) = vecSTA;
		end
	end
	%save data
	cellAggX{intSubType} = matAggTE_LocX;
	cellAggY{intSubType} = matAggTE_LocY;
	cellAggSize{intSubType} = matAggTE_Size;
	cellAggSync{intSubType} = matAggTE_Sync;
	cellAggOriIdx{intSubType} = vecAggOriIdx;
end

%% plot 1
%define data
vecOriRight = find(ismember(vecUnique,mod(-15:15:15,360)));
vecOriLeft = find(ismember(vecUnique,mod(180 + [-15:15:15],360)));
vecTime = vecWindowBinEdges(2:end) - median(diff(vecWindowBinEdges))/2;
figure;maxfig;
for intSubType=1:2
	
	matLeft = cellAggX{intSubType}(ismember(cellAggOriIdx{intSubType},vecOriLeft),:);
	matLeft = matLeft - nanmean(matLeft,2);
	matRight = cellAggX{intSubType}(ismember(cellAggOriIdx{intSubType},vecOriRight),:);
	matRight = matRight - nanmean(matRight,2);
	subplot(2,3,intSubType);
	hold on
	errorbar(vecTime,nanmean(matLeft,1),nanstd(matLeft,[],1)./sqrt(size(matLeft,1)-sum(isnan(matLeft),1)),'xb');
	errorbar(vecTime,nanmean(matRight,1),nanstd(matRight,[],1)./sqrt(size(matRight,1)-sum(isnan(matRight),1)),'xr');
	hold off
	legend({'Leftward stimulus','Rightward stimulus'},'location','best');
	ylabel('Normalized pupil X position');
	xlabel('Time after stim onset (s)');
	title(cellSubjectGroups{intSubType});
	fixfig;
end

%% save 1
drawnow;maxfig;
export_fig(fullpath(strTargetPath,['EyeMovementLeftRight.tif']));
export_fig(fullpath(strTargetPath,['EyeMovementLeftRight.pdf']));

%% plot 2
figure;maxfig;
for intSubType=1:2
	if intSubType == 2
		vecCol = vecColAlb;
	else
		vecCol = vecColBl6;
	end
	for intArea=1:2
		matSTA = zscore(cellSTA{intSubType,intArea},[],2);
		subplot(2,3,intSubType+(intArea-1)*3);
		
		errorbar(vecPupilSTA_centers,mean(matSTA,1),std(matSTA,[],1)/sqrt(size(matSTA,1)),'color',vecCol);
		title(sprintf('%s, %s',cellSubjectGroups{intSubType},cellAreaGroupsAbbr{intArea}));
		xlabel('Time after spike (s)');
		ylabel('Pupil x-movement (z-score)');
		fixfig
	end
end
%cellSTA{intSubType,intArea}

%% save 2
drawnow;maxfig;
export_fig(fullpath(strTargetPath,['EyeMovementSTA.tif']));
export_fig(fullpath(strTargetPath,['EyeMovementSTA.pdf']));
