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
vecPupilSTA_edges = [-0.5:0.025:0.5];
cellAggX = {};
cellAggY = {};
cellAggSize = {};
cellAggSync = {};
cellAggOriIdx = {};
cellSTA = cellfill(nan(0,numel(vecPupilSTA_edges)-1),[2 numel(cellAreaGroups)]);

%run
cellNameAP = arrayfun(@(x) x.sJson.file_preproAP,sExp,'uniformoutput',false);
cellExperiment = arrayfun(@(x) x.sJson.experiment,sExp,'uniformoutput',false);
cellRemove = {};%{'RecMA5_2021-03-01R01_g0_t0'};
indRemRecs = contains(cellExperiment,cellRemove);
indRemRecs2 = ~contains(cellNameAP,cellUseForEyeTracking);
cellSubjectType = arrayfun(@(x) x.sJson.subjecttype,sExp,'uniformoutput',false);
for intSubType=1%:2
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
		fprintf('%s; Rec %d/%d (%s) [%s]\n',strSubjectType,intRecIdx,numel(vecRunRecs),strRec,getTime);
		
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
		if ~isfield(sPupil,'vecBlinks')
			sPupil.vecBlinks = false(size(sPupil.vecTime));
		end
		dblPupilStart = dblBlankStart+dblPL+10;
		dblPupilStop = dblBlankStart+dblLongestBlankPeriod+dblPL-10;
		indUsePupil = sPupil.vecTime>dblPupilStart & sPupil.vecTime<dblPupilStop & ~sPupil.vecBlinks;
		vecPupilT = sPupil.vecTime(indUsePupil);
		vecPupilX = sPupil.vecCenterX(indUsePupil);
		vecPupilY = sPupil.vecCenterY(indUsePupil);
		vecPupildT = diff(vecPupilT);
		vecPupildX = diff(vecPupilX);
		vecPupildY = diff(vecPupilY);
		vecPupildD = sqrt(vecPupildX.^2 + vecPupildY.^2);
		
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
			intMaxSpikes = inf;
			if numel(vecSpikeT) > intMaxSpikes
				vecUseSpikes = unique(round(linspace(1,numel(vecSpikeT),intMaxSpikes)));
			else
				vecUseSpikes = 1:numel(vecSpikeT);
			end
			[matTE,vecPupilSTA_centers] = getRespMat(vecPupilT(2:end)-median(vecPupildT/2),vecPupildD,vecSpikeT(vecUseSpikes),vecPupilSTA_edges);
			vecSTA = nanmean(matTE,1);
			cellSTA{intSubType,intArea}(end+1,:) = vecSTA;
		end
	end
end

%save data
save(fullpath(strTargetPath,'PupilPlotsSTA'),'vecPupilSTA_centers','cellSTA');
			
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
		%remove missing recs
		matSTA(any(isnan(matSTA),2),:) = [];
		subplot(2,3,intSubType+(intArea-1)*3);
		
		errorbar(vecPupilSTA_centers,mean(matSTA,1),std(matSTA,[],1)/sqrt(size(matSTA,1)),'color',vecCol);
		title(sprintf('%s, %s',cellSubjectGroups{intSubType},cellAreaGroupsAbbr{intArea}));
		xlabel('Time after spike (s)');
		ylabel('Pupil movement (z-score)');
		fixfig
	end
end
%cellSTA{intSubType,intArea}


%% save 2
drawnow;maxfig;
export_fig(fullpath(strTargetPath,['EyeMovementBlankSTA.tif']));
export_fig(fullpath(strTargetPath,['EyeMovementBlankSTA.pdf']));
