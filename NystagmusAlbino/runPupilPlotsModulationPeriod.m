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
clear all;
runHeaderNOT;
intRunMaxArea = 2;

vecColAlb = [0.9 0 0];
vecColBl6 = lines(1);

%% get data
%pre-allocate
dblMoveEventThreshold = 0.05;
vecSTA_edges = [-4:0.1:4];cellSTA = cellfill(nan(0,numel(vecSTA_edges)-1),[3 intRunMaxArea 0]); %[move-type x V1/NOT x rec]
matSTAZetaP = nan([3 intRunMaxArea 0]); %[move-type x V1/NOT x rec]
cellZetaP = cell([3 intRunMaxArea 0]); %[move-type x V1/NOT x rec]

%run
cellNameAP = arrayfun(@(x) x.sJson.file_preproAP,sExp,'uniformoutput',false);
cellExperiment = arrayfun(@(x) x.sJson.experiment,sExp,'uniformoutput',false);
cellRemove = {};%{'RecMA5_2021-03-01R01_g0_t0'};
indRemRecs = contains(cellExperiment,cellRemove);
indRemRecs2 = ~contains(cellNameAP,cellUseForEyeTracking);
cellSubjectType = arrayfun(@(x) x.sJson.subjecttype,sExp,'uniformoutput',false);
for intSubType=1
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
	
	for intRecIdx=1:numel(vecRunRecs)
		intRec=vecRunRecs(intRecIdx);
		sRec = sExp(intRec);
		strRec = sRec.Name(1:(end-7));
		
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
		vecdT = diff(vecPupilT);
		vecdX = diff(vecPupilX);
		vecdY = diff(vecPupilY);
		hFig=figure;
		maxfig(hFig);
		for intMoveType = 1:3
			if intMoveType==1
				vecdM = vecdX./vecdT;
				strMoveType = 'horz';
			elseif intMoveType==2
				vecdM = vecdY./vecdT;
				strMoveType = 'vert';
			elseif intMoveType==3
				vecdM = sqrt(vecdX.^2 + vecdY.^2)./vecdT;
				strMoveType = 'tot';
			else
				error('not possible');
			end
			
			%5% cut-off
			[vecSortdM,vecReorder] = sort(vecdM);
			intCutOff = floor(numel(vecSortdM)*dblMoveEventThreshold);
			vecUpperValIdx = vecReorder((end-intCutOff-1):end);
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
			
			%remove overlapping starts
			dblMinDist = range(vecSTA_edges)/2;
			vecCheckEvents = find(diff(vecStartMoveT)<dblMinDist)+1;
			indRemEvent = false(size(vecStartMoveT));
			for intEvIdx=1:numel(vecCheckEvents)
				intEv = vecCheckEvents(intEvIdx);
				vecDiffT = vecStartMoveT(intEv) - vecStartMoveT(~indRemEvent);
				if any(vecDiffT > 0 & vecDiffT < dblMinDist)
					indRemEvent(intEv) = true;
				end
			end
			vecStartMoveT(indRemEvent) = [];
			
			%% spike-triggered average per area
			indUseCells = arrayfun(@(x) x.Violations1ms < 0.25 & abs(x.NonStationarity) < 0.25,sRec.sCluster(:));
			cellAreasPerCluster = {sRec.sCluster.Area};
			for intArea=1:intRunMaxArea
				% select cells
				strAreaGroup =  cellAreaGroupsAbbr{intArea};
				vecSelectCells = find(indUseCells(:) & flat(contains(cellAreasPerCluster,cellUseAreas{intArea},'IgnoreCase',true)));
				if isempty(vecSelectCells)
					continue;
				end
				cellSpikes = {sRec.sCluster(vecSelectCells).SpikeTimes};
				intCellNum = numel(vecSelectCells);
				vecSpikeT = sort(cell2vec(cellSpikes));
				%reduce spikes to max 100k
				intMaxSpikes = inf;
				if numel(vecSpikeT) > intMaxSpikes
					vecUseSpikes = unique(round(linspace(1,numel(vecSpikeT),intMaxSpikes)));
				else
					vecUseSpikes = 1:numel(vecSpikeT);
				end
				
				%% plot grand average over all cells
				%real
				figure(hFig);
				hAx = subplot(2,3,intMoveType+(intArea-1)*3);
				[vecMean,vecSEM,vecWindowBinCenters,matPET] = doPEP(vecSpikeT(vecUseSpikes),vecSTA_edges,vecStartMoveT,hAx);
				intResampNum=250;
				
				%jittered
				return
			end
		end
	end
end