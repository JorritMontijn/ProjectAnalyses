%[done/a) plot 3D locations of cells and show LR responses

%% load data and define groups
%strDataPath
%cellUseForEyeTracking
%strTargetPath
%cellUseAreas{1} = {'Primary visual','Posteromedial visual'};
%cellUseAreas{2} = {'nucleus of the optic tract'};
%cellUseAreas{3} = {'Field CA1','Field CA2','Field CA3'};
%cellAreaGroups = {'Vis. ctx','NOT','Hippocampus'};
%cellAreaGroupsAbbr = {'Ctx','NOT','Hip'};
%cellSubjectGroups = {'BL6','DBA'};
runHeaderNOT;

strAllenCCFPath = 'F:\Data\AllenCCF';
%strAllenCCFPath = 'E:\AllenCCF';
sAtlas = AL_PrepABA(strAllenCCFPath);
tv = sAtlas.tv;
av = sAtlas.av;
st = sAtlas.st;
if ~isfield(sExp(1).sCluster,'Waveform') || ~isfield(sExp(1).sCluster,'SelfArea')
	sExpNew = [];
	%try
	%	load(fullpath(strDataPath,'ProbeLocationPreProWorkspace'));
	%catch
		for intFile=1:numel(sExp)
			return
			%%
			sAP = sExp(intFile);
			fprintf('Loading waveforms for %d/%d: %s [%s]\n',intFile,numel(sExp),sAP.Name,getTime);
			if ~isfield(sAP,'sPupil')
				sAP.sPupil = [];
			end
			%load clustering data
			strSpikePath = sAP.sSources.sClustered.folder;
			sSpikes = loadKSdir(strSpikePath);
			strImPath = sLoad.sAP.sSources.sEphysAp.folder;
			strImFile = sLoad.sAP.sSources.sEphysAp.name;
			sMetaIM = DP_ReadMeta(strImFile,strImPath);
			if ~isfield(sAP.sCluster,'Waveform')
				[vecClustIdx,matClustWaveforms] = getWaveformPerCluster(sSpikes);
			end
			sAP.sSources.sMetaIM = sMetaIM;
			
			%calculate distance to area boundary
			vecDepth = sSpikes.ycoords;
			sProbeCoords = sAP.sSources.sProbeCoords;
			sLocCh = getBrainAreasPerChannel(sProbeCoords,sAtlas,false,numel(vecDepth));
			cellAreaPerCh = sLocCh.cellAreaPerCh;
			cellParentAreaPerCh = sLocCh.cellParentAreaPerCh;
			vecParentAreaPerCh_av = sLocCh.vecParentAreaPerCh_av;
			vecAreaBoundaries = sLocCh.vecAreaBoundaries;
			vecAreaCenters = sLocCh.vecAreaCenters;
			vecAreaLabels = sLocCh.vecAreaLabels;
			vecDistToBoundaryPerCh = sLocCh.vecDistToBoundaryPerCh;
			matCoordsPerCh = sLocCh.matCoordsPerCh;
			
			%check match
			[vecClustAreaId,cellClustAreaLabel,cellClustAreaFull] = PF_GetAreaPerCluster(sProbeCoords,vecDepth);
			indSame = strcmp(cellAreaPerCh(:),cellClustAreaFull(:));
			dblOverlap = sum(indSame)/length(indSame);
			if dblOverlap < 0.5
				fprintf('Area assignment agreement for %s is %.3f, please check!\n',strName,dblOverlap);
				error
			end
			
			%get cluster depths
			[spikeAmps, vecAllSpikeDepth] = templatePositionsAmplitudes(sSpikes.temps, sSpikes.winv, sSpikes.ycoords, sSpikes.spikeTemplates, sSpikes.tempScalingAmps);
			dblProbeLength = range(sSpikes.ycoords);
			vecAllSpikeClust = sSpikes.clu;
			
			%assign cluster data
			intClustNum = numel(sAP.sCluster);
			for intClust=1:intClustNum
				intClustIdx = sAP.sCluster(intClust).IdxClust;
				intDepth = dblProbeLength-round(median(vecAllSpikeDepth(vecAllSpikeClust==intClustIdx)));
				intDominantChannel = ceil(intDepth/10);
				
				%assign
				if ~isfield(sAP.sCluster,'Waveform')
					sAP.sCluster(intClust).Waveform = matClustWaveforms(sAP.sCluster(intClust).IdxClust == vecClustIdx,:);
				end
				sAP.sCluster(intClust).BoundDist = vecDistToBoundaryPerCh(intDominantChannel);
				sAP.sCluster(intClust).ParentArea = cellParentAreaPerCh{intDominantChannel};
				sAP.sCluster(intClust).SelfArea = cellAreaPerCh{intDominantChannel};
				sAP.sCluster(intClust).CoordsABA = matCoordsPerCh(:,intDominantChannel);
			end
			
			
			if isempty(sExpNew)
				sExpNew = sAP;
			else
				sExpNew(intFile) = sAP;
			end
			
		end
		sExp = sExpNew;
		clear sExpNew;
		
		save(fullpath(strDataPath,'ProbeLocationPreProWorkspace2'),'-v7.3');
		disp done
	%end
end


%best rec BL6: 20191216B5 (rec 17)
%best rec DBA: 20210212B2 (rec 5)

vecColAlb = [0.9 0 0];
vecColBl6 = lines(1);
intUseAreaNum = numel(cellUseAreas);

%% pre-allocate
cellAggBoundDist = cell(intUseAreaNum,2);
cellAggSpikeDur = cell(intUseAreaNum,2);
cellAggSpikePTR = cell(intUseAreaNum,2);
cellAggSpikeHz = cell(intUseAreaNum,2);
cellAggSpikeRLR = cell(intUseAreaNum,2);
cellAggCoords = cell(intUseAreaNum,2);

% run
cellNameAP = arrayfun(@(x) x.sJson.file_preproAP,sExp,'uniformoutput',false);
cellExperiment = arrayfun(@(x) x.sJson.experiment,sExp,'uniformoutput',false);
cellRemove = {};%{'RecMA5_2021-03-01R01_g0_t0'};
indRemRecs = contains(cellExperiment,cellRemove);
cellSubjectType = arrayfun(@(x) x.sJson.subjecttype,sExp,'uniformoutput',false);
for intSubType=1:2
	intPopCounter = 0;
	if intSubType == 1
		intBestRec = 17;
		strSubjectType = 'BL6';
	elseif intSubType == 2
		intBestRec = 5;
		strSubjectType = 'DBA';
	end
	indUseRecs = contains(cellSubjectType,strSubjectType);
	vecRunRecs = find(indUseRecs & ~(indRemRecs));
	%vecRunRecs = intBestRec;
	for intRecIdx=1:numel(vecRunRecs)
		intRec=vecRunRecs(intRecIdx);
		sRec = sExp(intRec);
		strName=[sRec.sJson.subject '_' sRec.sJson.date];
		
		% split cells into areas
		cellCellsPerArea = cell(1,numel(cellUseAreas));
		cellAreasPerCluster = {sRec.sCluster.Area};
		cellSelfPerCluster = {sRec.sCluster.SelfArea};
		cellParentAreasPerCluster = {sRec.sCluster.ParentArea};
		for intArea=1:numel(cellUseAreas)
			cellCellsPerArea{intArea} = contains(cellAreasPerCluster(:),cellUseAreas{intArea},'IgnoreCase',true);
		end
		vecCellsNrPerArea = cellfun(@sum,cellCellsPerArea);
		dblSampRateIM = str2double(sRec.sSources.sMetaAP.imSampRate);
		dblSampRateNI = str2double(sRec.sSources.sMetaNI.niSampRate);
		
		%%
		vecSame = false(1,numel(cellAreasPerCluster));
		for i=1:numel(cellAreasPerCluster)
			strCleanOld = strrep(strrep(cellAreasPerCluster{i},'layer ',''),'/','');
			strCleanOld(any(bsxfun(@eq,strCleanOld,arrayfun(@(x) num2str(x),0:9)'),1)) = [];
			strCleanNew = strrep(strrep(cellSelfPerCluster{i},'layer ',''),'/','');
			strCleanNew(any(bsxfun(@eq,strCleanNew,arrayfun(@(x) num2str(x),0:9)'),1)) = [];
			vecSame(i) = strcmp(strCleanOld,strCleanNew);
			%fprintf('%d: %s - %s - %s\n',i,cellParentAreasPerCluster{i},cellAreasPerCluster{i},cellSelfPerCluster{i});
		end
		dblAgreement=sum(vecSame)/numel(vecSame);
		if dblAgreement < 0.5
			fprintf('Area assignment agreement for %s is %.3f, please check!\n',strName,dblAgreement);
		end
		
		%get waveform in areas
		for intArea=1:intUseAreaNum
			%include?
			indUseCells = arrayfun(@(x) x.Violations1ms < 0.25 & abs(x.NonStationarity) < 0.25,sRec.sCluster(:));
			
			%build cell vectors
			vecSelectCells = find(indUseCells(:) & cellCellsPerArea{intArea}(:));
			if isempty(vecSelectCells) || isempty(sRec.sPupil)
				continue;
			end
			
			%collect data from all DG blocks
			cellStimType = cellfun(@(x) x.strExpType,sRec.cellBlock,'uniformoutput',false);
			vecBlocksDG = find(contains(cellStimType,'driftinggrating','IgnoreCase',true));
			matData = [];
			vecOrientation = [];
			for intBlockIdx=1:numel(vecBlocksDG)
				intBlock = vecBlocksDG(intBlockIdx);
				fprintf('Processing rec %d/%d (%s), B%d, area %s [%s]\n',intRecIdx,numel(vecRunRecs),strName,intBlockIdx,cellAreaGroupsAbbr{intArea},getTime);
				sBlock = sRec.cellBlock{intBlock};
				intPopCounter = intPopCounter + 1;
				
				
				if isfield(sBlock,'vecPupilStimOn')
					vecPupilStimOn = sBlock.vecPupilStimOn;
					vecPupilStimOff = sBlock.vecPupilStimOff;
				else
					vecPupilStimOn = sBlock.vecStimOnTime;
					vecPupilStimOff = sBlock.vecStimOffTime;
				end
				
				% get pupil data
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
				
				%split by ori
				sTrialObjects = sBlock.sStimObject(sBlock.vecTrialStimTypes);
				vecThisOrientation = cell2vec({sTrialObjects.Orientation});
				vecThisOrientation(indRemTrials) = [];
				[vecOriIdx,vecUnique,vecCounts,cellSelect,vecRepetition] = val2idx(vecThisOrientation);
				if numel(vecUnique) ~= 24,continue,end
				
				%get data matrix
				cellSpikeT = {sRec.sCluster(:).SpikeTimes};
				vecStimOnTime = sBlock.vecStimOnTime(~indRemTrials);
				vecStimOffTime = sBlock.vecStimOffTime(~indRemTrials);
				dblStimDur = median(vecStimOffTime-vecStimOnTime);
				matThisData = getSpikeCounts(cellSpikeT,vecStimOnTime,dblStimDur)./dblStimDur;
				
				%concatenate data
				matData = cat(2,matData,matThisData);
				vecOrientation = cat(2,vecOrientation,vecThisOrientation(:)');
			end
			
			%calc tuning curve
			matUseData = matData(vecSelectCells,:);
			sOut = getTuningCurves(matUseData,vecOrientation);
			vecTuningP_A = sOut.vecOriAnova;
			vecTuningP_R2 = sOut.vecFitP;
			[h, crit_p, vecTuningP_R2_corr] = fdr_bh(vecTuningP_R2);
			vecPrefOri = rad2deg(sOut.matFittedParams(:,1));
			indTunedCells = vecTuningP_R2_corr(:)'<0.05;
			
			%get RLR
			vecR = mean(sOut.matFittedResp(:,ismember(sOut.vecUniqueDegs,0)),2);
			vecL = mean(sOut.matFittedResp(:,ismember(sOut.vecUniqueDegs,180)),2);
			vecRLR = vecR ./ (vecR+vecL);
			
			%remove range 0
			vecBoundDist = cell2vec({sRec.sCluster(vecSelectCells).BoundDist})';
			vecRangeHz = range(matUseData,2)';
			indRem2=vecRangeHz==0 | vecBoundDist>20 | ~indTunedCells | isnan(vecRLR(:)');
			matUseData(indRem2,:)=[];
			vecBoundDist(indRem2)=[];
			vecRLR(indRem2)=[];
			vecSelectCells = vecSelectCells(~indRem2);
			
			%get distance to boundary and mean rates
			vecSpikeHz = mean(matUseData,2)';
			
			%get waveform props
			dblRecDur = max(cellfun(@max,{sRec.sCluster(vecSelectCells).SpikeTimes})) - min(cellfun(@min,{sRec.sCluster(vecSelectCells).SpikeTimes}));
			vecSpikeRate = cellfun(@numel,{sRec.sCluster(vecSelectCells).SpikeTimes})/dblRecDur;
			matAreaWaveforms = cell2mat({sRec.sCluster(vecSelectCells).Waveform}'); %[cell x sample]
			matCABA = cell2mat({sRec.sCluster(vecSelectCells).CoordsABA}); %[cell x sample]
			intNeurons=size(matAreaWaveforms,1);
			vecSpikeDur = nan(1,intNeurons);
			vecSpikePTR = nan(1,intNeurons);
			for intNeuron=1:intNeurons
				%find trough
				[dblTroughVal,intTrough]=min(matAreaWaveforms(intNeuron,:));
				[dblPeakVal,intTroughToPeak]=max(matAreaWaveforms(intNeuron,intTrough:end));
				intPeak = intTrough + intTroughToPeak - 1;
				
				dblTroughTime = intTrough/dblSampRateIM;
				dblTroughToPeakTime = intTroughToPeak/dblSampRateIM;
				dblPeakTime = intPeak/dblSampRateIM;
				vecSpikeDur(intNeuron) = dblTroughToPeakTime;
				vecSpikePTR(intNeuron) = abs(dblPeakVal/dblTroughVal);
			end
			cellAggBoundDist{intArea,intSubType} = cat(2,cellAggBoundDist{intArea,intSubType},vecBoundDist(:)');
			cellAggSpikeDur{intArea,intSubType} = cat(2,cellAggSpikeDur{intArea,intSubType},vecSpikeDur);
			cellAggSpikePTR{intArea,intSubType} = cat(2,cellAggSpikePTR{intArea,intSubType},vecSpikePTR);
			cellAggSpikeHz{intArea,intSubType} = cat(2,cellAggSpikeHz{intArea,intSubType},vecSpikeHz(:)');
			cellAggSpikeRLR{intArea,intSubType} = cat(2,cellAggSpikeRLR{intArea,intSubType},vecRLR(:)');
			cellAggCoords{intArea,intSubType} = cat(2,cellAggCoords{intArea,intSubType},matCABA);
		end
	end
end

%% plot 1
%av is [AP x DV x ML]
dblReduceBy = 2;
vecAP = (dblReduceBy/2):dblReduceBy:size(av,1);
vecDV = (dblReduceBy/2):dblReduceBy:size(av,2);
vecML = (dblReduceBy/2):dblReduceBy:size(av,3);

avRed = av(vecAP,vecDV,vecML);
SE = strel('sphere',1);
avCenter = avRed>1;
avCenter = imfill(avCenter,'holes');
avErode = imerode(avCenter,SE);
avEdge = avCenter - avErode;

SE1 = strel('disk',1);
intPoints = 7;
matLines = dblReduceBy*getTrace3D(avEdge,intPoints);

%h = plot3(matLines(:,3), matLines(:,1), matLines(:,2), 'Color', [0 0 0 0.3]);

%% find NOT
avNot=av==(st.index(contains(st.name,'nucleus of the optic tract','ignorecase',true))+1);
vecRangeNot1 = find(sum(sum(avNot,2),3));
vecRangeNot2 = find(sum(sum(avNot,1),3));
vecRangeNot3 = find(sum(sum(avNot,1),2));
vecNot1 = (vecRangeNot1(1)-1):(vecRangeNot1(end)+1);
vecNot2 = (vecRangeNot2(1)-1):(vecRangeNot2(end)+1);
vecNot3 = (vecRangeNot3(1)-1):(vecRangeNot3(end)+1);
avNot = avNot(vecNot1,vecNot2,vecNot3);
avCenter = imfill(avNot,'holes');
avErode = imerode(avNot,SE);
avEdge = avCenter - avErode;


intPoints = 15;
matLinesNot = getTrace3D(avEdge,intPoints,0);
%add offset
matLinesNot2=bsxfun(@plus,matLinesNot,([vecNot1(1) vecNot2(1) vecNot3(1)]));
indRem = matLinesNot2(:,3) < (size(av,3)/2);
matLinesNot2(indRem,:) = [];

%% plot all NOT cells at their respective locations and colour by RLR
figure;maxfig;
hold on
h = plot3(matLinesNot2(:,3), matLinesNot2(:,1), matLinesNot2(:,2), 'Color', [1 0 0 0.3]);
h.Annotation.LegendInformation.IconDisplayStyle = 'off';

cellMarker = {'x','o'};
for intSubType=1:2
	
	vecRLR = cellAggSpikeRLR{2,intSubType};
	matCABA = cellAggCoords{2,intSubType};
	scatter3(matCABA(3,:),matCABA(1,:),matCABA(2,:),[],vecRLR,'marker',cellMarker{intSubType});
	
	
end
hold off
legend(cellSubjectGroups)
xlabel('ML');
ylabel('AP');
zlabel('DV');
axis equal;
set(gca,'Zdir','reverse');
hc=colorbar;
clabel(hc,'R-L ratio');

fixfig;
grid off;
campos([1250 250 -100]);

title(sprintf('Recording locations in NOT'));

%save plot
drawnow;
export_fig([strTargetPath filesep sprintf('RecLocNot.tif')]);
export_fig([strTargetPath filesep sprintf('RecLocNot.pdf')]);

%% single axes
figure;maxfig;
subplot(2,3,1)%DV,ML,AP
[r,p]=corr(cellAggCoords{2,1}(1,:)',cellAggSpikeRLR{2,1}');
scatter(cellAggCoords{2,1}(1,:),cellAggSpikeRLR{2,1},'x');
xlabel('DV');
ylabel('R-L ratio');
title(sprintf('BL6, r(DV,RLR)=%.3f,p=%.3f',r,p));
fixfig;

subplot(2,3,2)%DV,ML,AP
[r,p]=corr(cellAggCoords{2,1}(2,:)',cellAggSpikeRLR{2,1}');
scatter(cellAggCoords{2,1}(2,:),cellAggSpikeRLR{2,1},'x');
xlabel('ML');
ylabel('R-L ratio');
title(sprintf('BL6, r(ML,RLR)=%.3f,p=%.3f',r,p));
fixfig;

subplot(2,3,3)%DV,ML,AP
[r,p]=corr(cellAggCoords{2,1}(3,:)',cellAggSpikeRLR{2,1}');
scatter(cellAggCoords{2,1}(3,:),cellAggSpikeRLR{2,1},'x');
xlabel('AP');
ylabel('R-L ratio');
title(sprintf('BL6, r(AP,RLR)=%.3f,p=%.3f',r,p));
fixfig;

subplot(2,3,4)%DV,ML,AP
[r,p]=corr(cellAggCoords{2,2}(1,:)',cellAggSpikeRLR{2,2}');
scatter(cellAggCoords{2,2}(1,:),cellAggSpikeRLR{2,2});
xlabel('DV');
ylabel('R-L ratio');
title(sprintf('DBA, r(DV,RLR)=%.3f,p=%.3f',r,p));
fixfig;

subplot(2,3,5)%DV,ML,AP
[r,p]=corr(cellAggCoords{2,2}(2,:)',cellAggSpikeRLR{2,2}');
scatter(cellAggCoords{2,2}(2,:),cellAggSpikeRLR{2,2});
xlabel('ML');
ylabel('R-L ratio');
title(sprintf('DBA, r(ML,RLR)=%.3f,p=%.3f',r,p));
fixfig;

subplot(2,3,6)%DV,ML,AP
[r,p]=corr(cellAggCoords{2,2}(3,:)',cellAggSpikeRLR{2,2}');
scatter(cellAggCoords{2,2}(3,:),cellAggSpikeRLR{2,2});
xlabel('AP');
ylabel('R-L ratio');
title(sprintf('DBA, r(AP,RLR)=%.3f,p=%.3f',r,p));
fixfig;

%save plot
drawnow;
export_fig([strTargetPath filesep sprintf('RecLocNotPerAx.tif')]);
export_fig([strTargetPath filesep sprintf('RecLocNotPerAx.pdf')]);
