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
vecSTA_edges = [-0.5:0.05:0.5];

%run
cellNameAP = arrayfun(@(x) x.sJson.file_preproAP,sExp,'uniformoutput',false);
cellExperiment = arrayfun(@(x) x.sJson.experiment,sExp,'uniformoutput',false);
cellRemove = {};%{'RecMA5_2021-03-01R01_g0_t0'};
indRemRecs = contains(cellExperiment,cellRemove);
indRemRecs2 = ~contains(cellNameAP,cellUseForEyeTracking);
cellSubjectType = arrayfun(@(x) x.sJson.subjecttype,sExp,'uniformoutput',false);
intSubType=1;
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

%pre-allocate
matNSTA = nan([3 intRunMaxArea numel(vecRunRecs) numel(vecSTA_edges)-1]); %[move-type x area x rec x time]

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
			
			%select only responsive cells
			vecZetaP = ones(1,numel(cellSpikes));
			intResampNum=100;
			for intCell=1:numel(cellSpikes)
				vecZetaP(intCell) = zetatest(cellSpikes{intCell},vecStartMoveT+vecSTA_edges(1),range(vecSTA_edges),intResampNum,0);
			end
			vecZetaP(vecZetaP>0.05)=[];
			cellSpikes(vecZetaP>0.05)=[];
			intCellNum = numel(cellSpikes);
			if intCellNum==0,continue;end
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
			[vecMean,vecSEM,vecWindowBinCenters,matPET] = doPEP(vecSpikeT(vecUseSpikes),vecSTA_edges,vecStartMoveT,-1);
			intResampNum=250;
			
			%jittered
			intTrials = numel(vecStartMoveT);
			dblJitterSize = 2;
			dblUseMaxDur = dblMinDist;
			vecJitterPerTrial = dblJitterSize*linspace(-dblUseMaxDur,dblUseMaxDur,intTrials)'; %new
			matJitterPerTrial = nan(intTrials,intResampNum);
			for intResampling=1:intResampNum
				matJitterPerTrial(:,intResampling) = vecJitterPerTrial(randperm(numel(vecJitterPerTrial)));
			end
			matRandMeans = nan(intResampNum,length(vecSTA_edges)-1);
			
			for intResampling=1:intResampNum
				%% get random subsample
				vecStimUseOnTime = vecStartMoveT' + matJitterPerTrial(:,intResampling);
				
				%get temp offset
				vecRandMean = doPEP(vecSpikeT(vecUseSpikes),vecSTA_edges,vecStimUseOnTime,-1);
				matRandMeans(intResampling,:)=vecRandMean;
			end
			
			%smooth
			dblSmoothWidth = 0.8;
			vecSmooth = 1;%normpdf(-10:10,0,dblSmoothWidth) ./ sum(normpdf(-10:10,0,dblSmoothWidth));
			vecMean = imfilt(vecMean,vecSmooth);
			matRandMeans = imfilt(matRandMeans,vecSmooth);
			vecRandMeansMu = mean(matRandMeans,1);
			vecRandMeansSd = std(matRandMeans,[],1);
			vecSTAC = vecSTA_edges(2:end)-mean(diff(vecSTA_edges))/2;
			
			% save data
			vecNormMean = (vecMean - mean(matRandMeans(:)))./std(matRandMeans(:));
			matNSTA(intMoveType,intArea,intRecIdx,:) = vecNormMean; %[move-type x V1/NOT x rec]
			
			%plot
			figure(hFig);
			hAx = subplot(2,3,intMoveType+(intArea-1)*3);
			errorfill(vecSTAC,vecRandMeansMu,2*vecRandMeansSd,[0.5 0.5 0.5]);
			hold on
			plot(vecSTAC,vecMean,'color',lines(1))
			hold off
			ylabel(sprintf('Real and jittered spiking rates (Hz)'));
			xlabel(sprintf('Time after %s eye-movement (s)',strMoveType));
			title(sprintf('%s (%s) - %s %s; Max-Z=%.1f [N=%d cells]',strRec,strSubjectType,strAreaGroup,strMoveType,max(abs(vecNormMean)),intCellNum),'interpreter','none');
			fixfig;grid off
		end
	end
	%% save figure
	drawnow;
	strFigName = ['PupilMoveSpikeModulation_' strRec];
	export_fig(fullpath(strTargetPath,[strFigName '.tif']));
	export_fig(fullpath(strTargetPath,[strFigName '.pdf']));
end

%% plot overall figure
figure;maxfig;
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
	for intArea=1:intRunMaxArea
		strAreaGroup =  cellAreaGroupsAbbr{intArea};
		hAx = subplot(2,3,intMoveType+(intArea-1)*3);
		matMeansZ = squeeze(matNSTA(intMoveType,intArea,:,:));
		matMeansZ(all(isnan(matMeansZ),2),:)=[];
		
			
			[h,p]=ttest(matMeansZ);
		[h_corr,c,p_corr] = fdr_bh(p);
		intRecNum = size(matMeansZ,1);
		hP=plot(vecSTAC,matMeansZ,'color',[0.5 0.5 0.5 0.5]);
		ylim([-3 7])
		hold on
		vecMeansZmu = mean(matMeansZ,1);
		vecMeansZsd = std(matMeansZ,[],1);
		errorbar(vecSTAC,vecMeansZmu,vecMeansZsd./sqrt(intRecNum),'color',lines(1));
		scatter(vecSTAC(h_corr==1),(max(get(gca,'ylim'))-0.5)*ones(1,sum(h_corr)));
		hold off
		ylabel(sprintf('Jitter-normalized spiking rates (z)'));
		xlabel(sprintf('Time after %s eye-movement (s)',strMoveType));
		title(sprintf('Grand average (%s) - %s %s; [N=%d recs]',strSubjectType,strAreaGroup,strMoveType,intRecNum),'interpreter','none');
		fixfig;grid off
		%set(hP,'LineWidth',1);
		
	end
end

% save figure
drawnow;
strFigName = ['PupilMoveSpikeModulationSummary'];
export_fig(fullpath(strTargetPath,[strFigName '.tif']));
export_fig(fullpath(strTargetPath,[strFigName '.pdf']));