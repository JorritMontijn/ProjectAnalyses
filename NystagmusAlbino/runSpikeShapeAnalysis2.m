%% exploratory analysis, no proper controls

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
%clear all;
runHeaderNOT;

%strAllenCCFPath = 'F:\Data\AllenCCF';
strAllenCCFPath = 'E:\AllenCCF';
sAtlas = AL_PrepABA(strAllenCCFPath);
tv = sAtlas.tv;
av = sAtlas.av;
st = sAtlas.st;
%if ~isfield(sExp(1).sCluster,'Waveform') || ~isfield(sExp(1).sCluster,'BoundDist')
	sExpNew = [];
	%try
	%	load(fullpath(strTargetPath,'ProbeLocationPreProWorkspace'));
	%catch
		for intFile=1:numel(sExp)
			%%
			sAP = sExp(intFile);
			fprintf('Loading waveforms for %d/%d: %s [%s]\n',intFile,numel(sExp),sAP.Name,getTime);
			if ~isfield(sAP,'sPupil')
				sAP.sPupil = [];
			end
			strName=sAP.Name;
			
			%load metadata
			sProbeCoords = sAP.sSources.sProbeCoords;
			cellAreaOrig = sProbeCoords.sProbeAdjusted.probe_area_full_per_depth;
			probe_n_coords=numel(cellAreaOrig);
			
			%check length
			dblProbeLengthCalc = sqrt(sum((diff(sProbeCoords.sProbeAdjusted.probe_vector_cart,1,1).*sProbeCoords.VoxelSize).^2));
			dblProbeLengthOrig = sProbeCoords.sProbeAdjusted.stereo_coordinates.ProbeLength;
			if abs(dblProbeLengthCalc - dblProbeLengthOrig) > dblProbeLengthOrig/1e6
				error([mfilename ':LengthMismatch'],sprintf('Probe lengths for %s do not agree, please check!\n',strName));
			end
			
			%calculate distance to area boundary
			vecFracDepth = linspace(0,1,probe_n_coords)';
			vecDepth = vecFracDepth*dblProbeLengthCalc;%sSpikes.ycoords;
			sLocCh = getBrainAreasPerChannel(sProbeCoords,sAtlas,false,vecDepth);
			vecAreaAv = sLocCh.vecAreaPerChAv;
			cellAreaPerCh = sLocCh.cellAreaPerCh;
			cellParentAreaPerCh = sLocCh.cellParentAreaPerCh;
			vecParentAreaPerCh_av = sLocCh.vecParentAreaPerCh_av;
			vecAreaBoundaries = sLocCh.vecAreaBoundaries;
			vecAreaCenters = sLocCh.vecAreaCenters;
			vecAreaLabels = sLocCh.vecAreaLabels;
			vecDistToBoundaryPerCh = sLocCh.vecDistToBoundaryPerCh;
			matCoordsPerCh = sLocCh.matCoordsPerCh;
			
			%check #1: AP vs PF way
			probe_area_ids = PH_GetProbeAreas(sProbeCoords.sProbeAdjusted.probe_vector_cart,sAtlas.av);
			probe_area_full = sAtlas.st.name(probe_area_ids);
			indSame = strcmp(probe_area_full(:),cellAreaOrig);
			dblChOverlapPF = sum(indSame)/length(indSame);
			if dblChOverlapPF~=1
				error([mfilename ':AreaMismatch'],sprintf('UPF and Acquipix functions return different areas for %s, please check!\n',strName));
			end
			
			%check #2
			indSame = strcmp(cellAreaPerCh(:),cellAreaOrig);
			dblChOverlap1 = sum(indSame)/length(indSame);
			if dblChOverlap1~=1
				error([mfilename ':AreaMismatch'],sprintf('Original and atlas-retrieved areas for %s do not agree, please check!\n',strName));
			end
			
			%get cluster depths
			strSpikePath = sAP.sSources.sClustered.folder;
			intClustNum = numel(sAP.sCluster);
			vecDepthOnProbe = nan(1,intClustNum);
			if isfolder(strSpikePath)
				sSpikes = loadKSdir(strSpikePath);
				if ~isfield(sAP.sCluster,'Waveform')
					[vecClustIdx,matClustWaveforms] = getWaveformPerCluster(sSpikes);
				end
				[spikeAmps, vecAllSpikeDepth, templateDepths] = templatePositionsAmplitudes(sSpikes.temps, sSpikes.winv, sSpikes.ycoords, sSpikes.spikeTemplates, sSpikes.tempScalingAmps);
				vecAllSpikeClust = sSpikes.clu;
				dblProbeLengthChanMap = max(sSpikes.ycoords(:));
				for intClust=1:intClustNum
					intClustIdx = sAP.sCluster(intClust).IdxClust;
					vecDepthOnProbe(intClust) = dblProbeLengthChanMap-round(median(vecAllSpikeDepth(vecAllSpikeClust==intClustIdx)));
				end
			else
				dblProbeLengthChanMap = 3840;
				for intClust=1:intClustNum
					vecDepthOnProbe(intClust) = (sAP.sCluster(intClust).Depth/dblProbeLengthOrig)*dblProbeLengthChanMap;
				end
			end
			
			%assign cluster data
			vecOldDepth = nan(1,intClustNum);
			vecNewDepth = nan(1,intClustNum);
			indSameAreaPerClust = false(1,intClustNum);
			for intClust=1:intClustNum
				intClustIdx = sAP.sCluster(intClust).IdxClust;
				intNewDepthRaw = vecDepthOnProbe(intClust);
				intNewDepth = intNewDepthRaw*(dblProbeLengthOrig/dblProbeLengthChanMap);
				intDominantChannel = ceil(intNewDepth/10);
				
				vecOldDepth(intClust) = sAP.sCluster(intClust).Depth;
				vecNewDepth(intClust) = intNewDepth;
				%check cluster match
				[vecClustAreaId,cellClustAreaLabel,cellClustAreaFull] = PF_GetAreaPerCluster(sProbeCoords,intNewDepthRaw);
				strOldArea = sAP.sCluster(intClust).Area;
				strNewArea = cellClustAreaFull{1};
				indSameAreaPerClust(intClust) = strcmp(strOldArea,strNewArea);
				
				%assign
				if ~isfield(sAP.sCluster,'Waveform') || isempty(sAP.sCluster(intClust).Waveform)
					sAP.sCluster(intClust).Waveform = matClustWaveforms(sAP.sCluster(intClust).IdxClust == vecClustIdx,:);
				end
				if isempty(sAP.sCluster(intClust).Waveform),error;end
				sAP.sCluster(intClust).BoundDist = vecDistToBoundaryPerCh(intDominantChannel);
				sAP.sCluster(intClust).ParentArea = cellParentAreaPerCh{intDominantChannel};
				sAP.sCluster(intClust).SelfArea = cellAreaPerCh{intDominantChannel};
				sAP.sCluster(intClust).CoordsABA = matCoordsPerCh(intDominantChannel,:);
			end
			if any(abs(vecOldDepth - vecNewDepth) > dblProbeLengthOrig/1e6)
				error([mfilename ':DepthMismatch'],sprintf('Original and recalculated depths for %s do not agree, please check!\n',strName));
			end
			dblClustOverlap = sum(indSameAreaPerClust)/length(indSameAreaPerClust);
			if dblClustOverlap<0.99
				error([mfilename ':AreaMismatch'],sprintf('Original and atlas-retrieved cluster areas for %s do not agree, please check!\n',strName));
			end
			
			if isempty(sExpNew)
				sExpNew = sAP;
			else
				sExpNew(intFile) = sAP;
			end
			
		end
		sExp = sExpNew;
		clear sExpNew;
		
		save(fullpath(strTargetPath,'ProbeLocationPreProWorkspace'),'-v7.3');
		disp done
	%end
%end


%best rec BL6: 20191216B5 (rec 17)
%best rec DBA: 20210212B2 (rec 5)

vecColAlb = [0.9 0 0];
vecColBl6 = lines(1);

%% pre-allocate
intUseAreaNum = 2;%numel(cellUseAreas);
cellAggBoundDist = cell(intUseAreaNum,2);
cellAggSpikeDur = cell(intUseAreaNum,2);
cellAggSpikePTR = cell(intUseAreaNum,2);
cellAggSpikeHz = cell(intUseAreaNum,2);
cellAggSpikeRLR = cell(intUseAreaNum,2);

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
		for intArea=1:numel(cellUseAreas)
			cellCellsPerArea{intArea} = contains(cellAreasPerCluster(:),cellUseAreas{intArea},'IgnoreCase',true);
		end
		vecCellsNrPerArea = cellfun(@sum,cellCellsPerArea);
		dblSampRateIM = str2double(sRec.sSources.sMetaAP.imSampRate);
		dblSampRateNI = str2double(sRec.sSources.sMetaNI.niSampRate);
		
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
				if isempty(dblTroughToPeakTime)
					error
				end
				vecSpikeDur(intNeuron) = dblTroughToPeakTime;
				vecSpikePTR(intNeuron) = abs(dblPeakVal/dblTroughVal);
			end
			if numel(vecSpikePTR) ~= numel(vecSpikeHz),error;end
			cellAggBoundDist{intArea,intSubType} = cat(2,cellAggBoundDist{intArea,intSubType},vecBoundDist(:)');
			cellAggSpikeDur{intArea,intSubType} = cat(2,cellAggSpikeDur{intArea,intSubType},vecSpikeDur);
			cellAggSpikePTR{intArea,intSubType} = cat(2,cellAggSpikePTR{intArea,intSubType},vecSpikePTR);
			cellAggSpikeHz{intArea,intSubType} = cat(2,cellAggSpikeHz{intArea,intSubType},vecSpikeHz(:)');
			cellAggSpikeRLR{intArea,intSubType} = cat(2,cellAggSpikeRLR{intArea,intSubType},vecRLR(:)');
		end
	end
end

%% plot 1
dblPTT = 0.5;
dblSWT = 0.5/1000;
cellMarker = {'x','o'};
figure
maxfig;
for intSubType=1:2
	for intArea=1:intUseAreaNum
		if isempty(cellAggSpikeDur{intArea,intSubType})
			continue;
		end
		subplot(2,3,intArea+(intSubType-1)*3)
		colormap('redbluepurple')
		hold on;
		h1=plot(dblSWT*[1 1]*1000,[0 1.2],'color', [0.7 0.7 0.7],'linestyle','-');
		h2=plot([0 1.5],dblPTT*[1 1],'color', [0.7 0.7 0.7],'linestyle','-');
		if intArea == 1
			vecCol = [0 0 1];
		elseif intArea == 2
			vecCol = [1 0 0];
		end
		strArea = cellAreaGroupsAbbr{intArea};
		
		strSubjectType = cellSubjectGroups{intSubType};
		%scatter(1000*(cellAggSpikeDur{intArea,intSubType}+(rand(size(cellAggSpikeDur{intArea,intSubType}))-0.5)/dblSampRateIM),cellAggSpikePTR{intArea,intSubType},20+10*cellAggSpikeHz{intArea,intSubType},vecCol,'marker','.');
		scatter(1000*(cellAggSpikeDur{intArea,intSubType}+(rand(size(cellAggSpikeDur{intArea,intSubType}))-0.5)/dblSampRateIM)',cellAggSpikePTR{intArea,intSubType}',...
			50,log(1+cellAggSpikeHz{intArea,intSubType})','marker','.');
		h=colorbar;
		clabel(h,'log(1+rate)');
		hold off;
		ylim([0 1.2]);
		ylabel('Peak-to-trough ratio');
		xlabel('Spike width (ms)');
		title(sprintf('%s, %s; waveform props, mean rate=%.1fHz',strSubjectType,cellAreaGroupsAbbr{intArea},mean(cellAggSpikeHz{intArea,intSubType})));
		fixfig;
		grid off
	end
end


%BL6
intNarrowCtxBL6 = sum(cellAggSpikeDur{1,1} < dblSWT & cellAggSpikePTR{1,1} > dblPTT); %Ctx BL6
intBroadCtxBL6 = sum(cellAggSpikeDur{1,1} > dblSWT & cellAggSpikePTR{1,1} < dblPTT); %Ctx BL6
dblRatioCtxBL6 = intNarrowCtxBL6 / (intBroadCtxBL6 + intNarrowCtxBL6);

intNarrowNotBL6 = sum(cellAggSpikeDur{2,1} < dblSWT & cellAggSpikePTR{2,1} > dblPTT); %Ctx BL6
intBroadNotBL6 = sum(cellAggSpikeDur{2,1} > dblSWT & cellAggSpikePTR{2,1} < dblPTT); %Ctx BL6
dblRatioNotBL6 = intNarrowNotBL6 / (intBroadNotBL6 + intNarrowNotBL6);


[dblMu_CtxBL6,vecCI_CtxBL6] = binofit(intNarrowCtxBL6,intBroadCtxBL6 + intNarrowCtxBL6);
[dblMu_NotBL6,vecCI_NotBL6] = binofit(intNarrowNotBL6,intBroadNotBL6 + intNarrowNotBL6);
[p_DiffBL6,z]=bino2test(intNarrowCtxBL6,intBroadCtxBL6 + intNarrowCtxBL6,intNarrowNotBL6,intBroadNotBL6 + intNarrowNotBL6);

subplot(2,3,3)
hold on
errorbar(1,dblMu_CtxBL6,vecCI_CtxBL6(1) - dblMu_CtxBL6,vecCI_CtxBL6(2) - dblMu_CtxBL6,'x','color',[0 0 1],'CapSize',20);
errorbar(2,dblMu_NotBL6,vecCI_NotBL6(1) - dblMu_NotBL6,vecCI_NotBL6(2) - dblMu_NotBL6,'x','color',[1 0 0],'CapSize',20);
hold off
ylabel(sprintf('Fraction of narrow spiking (%s +/- 95%% CI)',getGreek('mu')));
set(gca,'xtick',[1 2],'xticklabel',cellAreaGroupsAbbr(1:2));
xlim([0.5 2.5]);
ylim([0 1]);
title(sprintf('BL6:Ctx,mu=%.4f,NOT,mu=%.4f; 2bino, p=%.1e',dblMu_CtxBL6,dblMu_NotBL6,p_DiffBL6));
fixfig;grid off;drawnow;

%DBA
intNarrowCtxDBA = sum(cellAggSpikeDur{1,2} < dblSWT & cellAggSpikePTR{1,2} > dblPTT); %Ctx BL6
intBroadCtxDBA = sum(cellAggSpikeDur{1,2} > dblSWT & cellAggSpikePTR{1,2} < dblPTT); %Ctx BL6
dblRatioCtxDBA = intNarrowCtxDBA / (intBroadCtxDBA + intNarrowCtxDBA);


intNarrowNotDBA = sum(cellAggSpikeDur{2,2} < dblSWT & cellAggSpikePTR{2,2} > dblPTT); %Ctx BL6
intBroadNotDBA = sum(cellAggSpikeDur{2,2} > dblSWT & cellAggSpikePTR{2,2} < dblPTT); %Ctx BL6
dblRatioNotDBA = intNarrowNotDBA / (intBroadNotDBA + intNarrowNotDBA);

[dblMu_CtxDBA,vecCI_CtxDBA] = binofit(intNarrowCtxDBA,intBroadCtxDBA + intNarrowCtxDBA);
[dblMu_NotDBA,vecCI_NotDBA] = binofit(intNarrowNotDBA,intBroadNotDBA + intNarrowNotDBA);
[p_DiffDBA,z]=bino2test(intNarrowCtxDBA,intBroadCtxDBA + intNarrowCtxDBA,intNarrowNotDBA,intBroadNotDBA + intNarrowNotDBA);


subplot(2,3,6)
hold on
errorbar(1,dblMu_CtxDBA,vecCI_CtxDBA(1) - dblMu_CtxDBA,vecCI_CtxDBA(2) - dblMu_CtxDBA,'x','color',[0 0 1],'CapSize',20);
errorbar(2,dblMu_NotDBA,vecCI_NotDBA(1) - dblMu_NotDBA,vecCI_NotDBA(2) - dblMu_NotDBA,'x','color',[1 0 0],'CapSize',20);
hold off
ylabel(sprintf('Fraction of narrow spiking (%s +/- 95%% CI)',getGreek('mu')));
set(gca,'xtick',[1 2],'xticklabel',cellAreaGroupsAbbr(1:2));
xlim([0.5 2.5]);
ylim([0 1]);
title(sprintf('DBA:Ctx,mu=%.4f,NOT,mu=%.4f; 2bino,p=%.1e',dblMu_CtxDBA,dblMu_NotDBA,p_DiffDBA));
fixfig;grid off;drawnow;

%save plot
drawnow;
export_fig([strTargetPath filesep sprintf('SpikeShapes2.tif')]);
export_fig([strTargetPath filesep sprintf('SpikeShapes2.pdf')]);

%% plot 2
try
%data
%cellAggBoundDist{intArea,intSubType} = cat(2,cellAggBoundDist{intArea,intSubType},vecBoundDist(:)');
%cellAggSpikeHz{intArea,intSubType} = cat(2,cellAggSpikeHz{intArea,intSubType},vecSpikeHz(:)');
%cellAggSpikeRLR{intArea,intSubType} = cat(2,cellAggSpikeRLR{intArea,intSubType},vecRLR(:)');
figure
maxfig;
for intSubType=1:2
	for intArea=1:intUseAreaNum
		subplot(2,3,intArea+(intSubType-1)*3)
		colormap('redbluepurple')
		hold on;
		if intArea == 1
			vecCol = [0 0 1];
		elseif intArea == 2
			vecCol = [1 0 0];
		end
		strArea = cellAreaGroupsAbbr{intArea};
		
		strSubjectType = cellSubjectGroups{intSubType};
		scatter(log(1+cellAggSpikeHz{intArea,intSubType}),cellAggSpikeRLR{intArea,intSubType},...
			50,log(1+cellAggSpikeHz{intArea,intSubType}),'marker','.');
		%h=colorbar;
		%clabel(h,'log(1+Spike rate)');
		hold off;
		ylim([0 1]);
		ylabel('Right-to-left ratio');
		xlabel('Spiking rate (log(1+Hz))');
		title(sprintf('%s, %s;',strSubjectType,cellAreaGroupsAbbr{intArea}));
		fixfig;
		grid off
	end
end

% test
[p_DiffRLR,h,stats] = ranksum(cellAggSpikeRLR{2,1},cellAggSpikeRLR{2,2});
subplot(2,3,3)
hold on
errorbar(1,mean(cellAggSpikeRLR{2,1}),std(cellAggSpikeRLR{2,1})./sqrt(numel(cellAggSpikeRLR{2,1})),'x','color',vecColBl6,'CapSize',20);
errorbar(2,mean(cellAggSpikeRLR{2,2}),std(cellAggSpikeRLR{2,2})./sqrt(numel(cellAggSpikeRLR{2,2})),'x','color',vecColAlb,'CapSize',20);
hold off
ylabel(sprintf('Right-to-left response ratio'));
set(gca,'xtick',[1 2],'xticklabel',cellSubjectGroups(1:2));
xlim([0.5 2.5]);
title(sprintf('BL6,mu=%.4f,DBA,mu=%.4f; wilcoxon, p=%.3f',mean(cellAggSpikeRLR{2,1}),mean(cellAggSpikeRLR{2,2}),p_DiffRLR));
fixfig;grid off;drawnow;


%save plot
drawnow;
export_fig([strTargetPath filesep sprintf('SpikeShapes3.tif')]);
export_fig([strTargetPath filesep sprintf('SpikeShapes3.pdf')]);
catch
end