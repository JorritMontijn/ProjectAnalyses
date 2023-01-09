%% load data and define groups
%strDataPath
%cellUseForEyeTracking
%strTargetPath
%cellUseAreas{1} = {'Primary visual','Posteromedial visual'};
%cellUseAreas{2} = {'nucleus of the optic tract'};
%cellUseAreas{3} = {'superior colliculus'};
%cellAreaGroups = {'Vis. ctx','NOT','Hippocampus'};
%cellAreaGroupsAbbr = {'Ctx','NOT','Hip'};
%cellSubjectGroups = {'BL6','DBA'};
runHeaderNOT;

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
cellSTA = cellfill(nan(0,numel(vecPupilSTA_edges)-1),[2 numel(cellAreaGroups)])

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
		dblAverageMouseHeadTiltInSetup = -15; %average tilt due to head-bar placement; procedures changed between BL6 and DBA experiments
		boolInvertX = 1; %other eye was recorded, so temporonasal is nasotemporal
	elseif intSubType == 2
		strSubjectType = 'DBA';
		dblOffsetT=0;
		dblAverageMouseHeadTiltInSetup = -15;
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
		intRec=vecRunRecs(intRecIdx);
		sRec = sExp(intRec);
		strRec = sRec.Name(1:(end-7));
		
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
		
		%% get stim absence
		cellStimType = cellfun(@(x) x.strExpType,sRec.cellBlock,'uniformoutput',false);
		vecBlocksDG = find(contains(cellStimType,'driftinggrating','IgnoreCase',true));
		vecBlocksNM = find(contains(cellStimType,'naturalmovie','IgnoreCase',true));
		%get timing for DG
		intBlock = vecBlocksDG(1);
		sBlock = sRec.cellBlock{intBlock};
		if isfield(sBlock,'vecPupilStimOn')
			vecPupilLatency = sBlock.vecPupilStimOn-sBlock.vecStimOnTime;
		else
			vecPupilLatency = 0;
		end
		dblPL = median(vecPupilLatency);
		vecAllStimOnsets = [];
		vecAllStimOffsets = [];
		vecAllStimBlock = [];
		dblRecDur = max(cell2vec({sRec.sCluster.SpikeTimes}));
		for intBlockIdx=1:numel(sRec.cellBlock)
			sBlock = sRec.cellBlock{intBlockIdx};
			vecAllStimOnsets = cat(2,vecAllStimOnsets,sBlock.vecStimOnTime);
			vecAllStimOffsets = cat(2,vecAllStimOffsets,sBlock.vecStimOffTime);
			vecAllStimBlock = [];
		end
		vecISI = [vecAllStimOnsets dblRecDur] - [0 vecAllStimOffsets];
		[dblLongestBlankPeriod,intStim] = max(vecISI);
		dblBlankStart = vecAllStimOffsets(intStim-1);
		
		%get eye movement events
		sPupil = sRec.sPupil;
		dblPupilStart = dblBlankStart+dblPL+10;
		dblPupilStop = dblBlankStart+dblLongestBlankPeriod+dblPL-10;
		indUsePupil = sPupil.vecTime>dblPupilStart & sPupil.vecTime<dblPupilStop & ~sPupil.vecBlinks;
		vecPupilT = sPupil.vecTime(indUsePupil);
		vecPupilX = sPupil.vecCenterX(indUsePupil);
		vecPupilY = sPupil.vecCenterY(indUsePupil);
		vecdT = diff(vecPupilT);
		vecdX = diff(vecPupilX);
		vecdY = diff(vecPupilY);
		%vecdM = sqrt(vecdX.^2 + vecdY.^2)./vecdT;
		vecdM = vecdX;
		%vecdM = vecdY;
		strMoveType = 'horz';
		
		%5% cut-off
		[vecSortdM,vecReorder] = sort(vecdM);
		int5P = floor(numel(vecSortdM)*0.05);
		vecUpperValIdx = vecReorder((end-int5P-1):end);
		indUpperVals = false(size(vecdM));
		indUpperVals(vecUpperValIdx) = true;
		vecUpperVals = vecdM(vecUpperValIdx);
		
		%get epochs
		vecStartMoveIdx = 1+find(diff(indUpperVals)==1);
		vecStopMoveIdx = find(diff(indUpperVals)==-1);
		if vecStartMoveIdx(1) > vecStopMoveIdx(1)
			vecStartMoveIdx = cat(2,1,vecStartMoveIdx);
		end
		if vecStartMoveIdx(end) > vecStopMoveIdx(end)
			vecStopMoveIdx = cat(2,vecStopMoveIdx,vecStartMoveIdx(end));
		end
		vecMoveDurSecs = nan(size(vecStopMoveIdx));
		for intEpoch=1:numel(vecMoveDurSecs)
			vecMoveDurSecs(intEpoch) = sum(vecdT(vecStartMoveIdx(intEpoch):vecStopMoveIdx(intEpoch)));
		end
		dblMedianFrameDur = median(vecdT);
		vecCutOffs = dblMedianFrameDur*(0.5+(1:(max(vecStopMoveIdx-vecStartMoveIdx)+1)));
		vecStartMoveT = vecPupilT(vecStartMoveIdx)-dblPL;
		
		
		%% spike-triggered average per area
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
			
			%% plot grand average over all cells
			doPEP(vecSpikeT(vecUseSpikes),[-4:0.05:4],vecStartMoveT);
			intResampNum=250;
			intPlot=0;
			[dblZetaP,sZeta] = zetatest(vecSpikeT(vecUseSpikes),vecStartMoveT-2,4,intResampNum,intPlot);
			ylabel(sprintf('Total firing rate in %s (Hz)',strAreaGroup));
			xlabel(sprintf('Time after %s eye-movement (s)',strMoveType));
			title(sprintf('%s (%s) - %s %s',strRec,strSubjectType,strAreaGroup,strMoveType),'interpreter','none');
			
			%% go through individual cells
			error do this
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
