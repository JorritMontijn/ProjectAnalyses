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
clearvars -except sExp;
runHeaderNOT;


sAtlas = AL_PrepABA(strAllenCCFPath);
tv = sAtlas.tv;
av = sAtlas.av;
st = sAtlas.st;
if ~isfield(sExp(1).sCluster,'Waveform') || ~isfield(sExp(1).sCluster,'SelfArea')
	sExpNew = [];
	try
		disp('loading workspace')
		load(fullpath(strDataPath,'ProbeLocationPreProWorkspace2'));
	catch
		%%
		for intFile=1:numel(sExp)
			return
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
				if ~isfield(sAP.sCluster,'Waveform')
					sAP.sCluster(intClust).Waveform = matClustWaveforms(sAP.sCluster(intClust).IdxClust == vecClustIdx,:);
				end
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
			
			
			%{
			%% to do: plot npx through brain
			%get locations along probe
			[probe_area_ids,probe_area_boundaries,probe_area_centers] = PH_GetProbeAreas(sProbeCoords.sProbeAdjusted.probe_vector_cart,sAtlas.av);
			probe_area_idx = probe_area_ids(round(probe_area_centers));
			probe_area_labels = sAtlas.st.acronym(probe_area_idx);
			probe_area_full = sAtlas.st.name(probe_area_idx);
			
			%add locations to GUI data
			probe_area_ids_per_depth = probe_area_ids;
			probe_area_labels_per_depth = sAtlas.st.acronym(probe_area_ids);
			probe_area_full_per_depth = sAtlas.st.name(probe_area_ids);
			dblOverlap = sum(sProbeCoords.sProbeAdjusted.probe_area_ids_per_depth==probe_area_ids_per_depth)./numel(probe_area_ids_per_depth)
			
			%PF depth calc
			vecAllSpikeTimes = sSpikes.st;
			vecAllSpikeClust = sSpikes.clu;
			[vecTemplateIdx,dummy,spike_templates_reidx] = unique(vecAllSpikeClust);
			vecClustIdx = cell2vec({sAP.sCluster.IdxClust});
			dblOverlapC = sum(vecTemplateIdx==vecClustIdx)./numel(vecTemplateIdx)
			
			vecTemplateDepths = nan(1,numel(vecTemplateIdx));
			for intCluster=1:numel(vecTemplateIdx)
				intClustIdx = vecTemplateIdx(intCluster);
				vecTemplateDepths(intCluster) = dblProbeLength - templateDepths(intClustIdx+1);
			end
			
			dblOverlapX = sum(vecTemplateDepths'==cell2vec({sAP.sCluster.Depth}))./numel(vecTemplateDepths)
			%}
		end
		sExp = sExpNew;
		clear sExpNew;
		
		save(fullpath(strDataPath,'ProbeLocationPreProWorkspace2'),'-v7.3');
		disp done
	end
end


%best rec BL6: 20191216B5 (rec 17)
%best rec DBA: 20210212B2 (rec 5)

vecColAlb = [0.9 0 0];
vecColBl6 = lines(1);
cellUseAreas(3) = [];
intUseAreaNum = numel(cellUseAreas);

%% pre-allocate
cellAggBoundDist = cell(intUseAreaNum,2);
cellAggSpikeDur = cell(intUseAreaNum,2);
cellAggSpikePTR = cell(intUseAreaNum,2);
cellAggSpikeHz = cell(intUseAreaNum,2);
cellAggSpikeRLR = cell(intUseAreaNum,2);
cellAggCoords = cell(intUseAreaNum,2);
cellAggCoords2 = cell(intUseAreaNum,2);
cellAggArea = cell(intUseAreaNum,2);
cellAggSelfArea = cell(intUseAreaNum,2);
cellAggSourceRec = cell(intUseAreaNum,2);

% run
cellNameAP = arrayfun(@(x) x.sJson.file_preproAP,sExp,'uniformoutput',false);
cellExperiment = arrayfun(@(x) x.sJson.experiment,sExp,'uniformoutput',false);
cellRemove = {};%{'RecMA5_2021-03-01R01_g0_t0'};
indRemRecs = contains(cellExperiment,cellRemove);
cellSubjectType = arrayfun(@(x) x.sJson.subjecttype,sExp,'uniformoutput',false);
for intSubType=[2 1]
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
		
		%% compare old and new areas
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
		if dblAgreement~=1
			fprintf('Area assignment agreement for %s is %.3f, please check!\n',strName,dblAgreement);
		end
		
		%% retrieve using alternative 2
		%get locations along probe
		sProbeCoords = sRec.sSources.sProbeCoords;
		[probe_area_ids,probe_area_boundaries,probe_area_centers,matLocCh] = PH_GetProbeAreas(sProbeCoords.sProbeAdjusted.probe_vector_cart,sAtlas.av);
		probe_area_idx = probe_area_ids(round(probe_area_centers));
		probe_area_labels = sAtlas.st.acronym(probe_area_idx);
		probe_area_full = sAtlas.st.name(probe_area_idx);
		
		cellAreasPerCluster3 = sProbeCoords.sProbeAdjusted.probe_area_full_per_cluster;
		
		sProbeAdjusted2 = sProbeCoords.sProbeAdjusted;
		sProbeAdjusted2.probe_area_ids = probe_area_idx;
		
		%add locations to GUI data
		sProbeAdjusted2.probe_area_ids_per_depth = probe_area_ids;
		sProbeAdjusted2.probe_area_labels_per_depth = sAtlas.st.acronym(probe_area_ids);
		sProbeAdjusted2.probe_area_full_per_depth = sAtlas.st.name(probe_area_ids);
		sLocCh = getBrainAreasPerChannel(sProbeAdjusted2,sAtlas,false,numel(sProbeAdjusted2.probe_area_full_per_depth));
		cellAreasPerChannel = sProbeAdjusted2.probe_area_full_per_depth;
		
		%get clusters
		dblProbeLengthProbe = sRec.sSources.sProbeCoords.ProbeLengthMicrons;
		dblProbeLengthPax = sRec.stereo_coordinates.ProbeLength;
		dblProbeLengthSph = sRec.sSources.sProbeCoords.sProbeAdjusted.probe_vector_sph(6)*mean(sRec.sSources.sProbeCoords.VoxelSize);
		dblProbeLengthCart = sqrt(sum((diff(sRec.sSources.sProbeCoords.sProbeAdjusted.probe_vector_cart,1,1).*sRec.sSources.sProbeCoords.VoxelSize).^2));
		dblProbeLengthBregma = sRec.sSources.sProbeCoords.sProbeAdjusted.probe_vector_bregma(6);
		
		vecDepth = cell2vec({sRec.sCluster.Depth});
		vecDepthOnProbe = (vecDepth/dblProbeLengthSph)*dblProbeLengthProbe;
		
		vecDepthAP = cell2vec({sRec.sCluster.Depth});
		vecDepthUPF = mean(sProbeCoords.VoxelSize)*sProbeCoords.sProbeAdjusted.depth_per_cluster;
		
		sProbeCoords2.sProbeAdjusted = sProbeAdjusted2;
		sProbeCoords2.ProbeLengthMicrons = sProbeCoords.ProbeLengthMicrons;
		[vecClustAreaId,cellClustAreaLabel,cellClustAreaFull,vecVoxelDepth] = PF_GetAreaPerCluster(sProbeCoords2,vecDepthOnProbe);
		vecUsedDepth = (min(max(round(vecVoxelDepth),1),floor(sProbeCoords2.sProbeAdjusted.probe_vector_sph(end)))/floor(sProbeCoords2.sProbeAdjusted.probe_vector_sph(end)))*sProbeCoords.ProbeLengthMicrons;
		vecUsedLocsPerCluster = min(max(round(vecVoxelDepth),1),floor(sProbeCoords2.sProbeAdjusted.probe_vector_sph(end)));
		matLocPerCluster = matLocCh(:,vecUsedLocsPerCluster);
		
		sLocCh = getBrainAreasPerChannel(sProbeAdjusted2,sAtlas,false,vecDepthUPF);
		cellAreasPerClusterNew = sLocCh.cellAreaPerCh;
		dblSimilarity = mean(strcmp(cellAreasPerClusterNew,cellClustAreaFull));
		
		dblSimilarityNewAPToOldUPF = mean(strcmp(cellAreasPerClusterNew,cellAreasPerCluster'));
		
		dblSimilarityNewUPFdToOldUPF = mean(strcmp(cellClustAreaFull,cellAreasPerCluster'));
		
		dblSimilarityNewAPToOldAP = mean(strcmp(cellAreasPerClusterNew,cellSelfPerCluster'));
		
		dblSimilarityNewUPFToOldAP = mean(strcmp(cellClustAreaFull,cellSelfPerCluster'));
		
		dblSimilarityOldUPFToOldAP = mean(strcmp(cellAreasPerCluster,cellSelfPerCluster));
		
		%% waveforms
		%get waveform in areas
		for intArea=1:intUseAreaNum
			%include?
			indUseCells = arrayfun(@(x) x.Violations1ms < 0.25 & abs(x.NonStationarity) < 0.25,sRec.sCluster(:));
			
			%build cell vectors
			vecSelectCells = find(indUseCells(:) & cellCellsPerArea{intArea}(:));
			if isempty(vecSelectCells)
				fprintf('No cells in %s for %s: skipping...\n',cellAreaGroupsAbbr{intArea},strName);
				continue;
			end
			if 0%isempty(sRec.sPupil)
				fprintf('No pupil data for %s: skipping...\n',strName);
				continue;
			end
			
			%% collect data from all DG blocks
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
				if isempty(sRec.sPupil)
					indRemTrials = false(size(vecPupilStimOn));
				else
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
				end
				
				%split by ori
				sTrialObjects = sBlock.sStimObject(sBlock.vecTrialStimTypes(1:numel(sBlock.vecStimOnTime)));
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
			
			%get source rec
			vecSourceRec = ones(size(vecSpikeHz))*intRec;
			
			%get waveform props
			dblRecDur = max(cellfun(@max,{sRec.sCluster(vecSelectCells).SpikeTimes})) - min(cellfun(@min,{sRec.sCluster(vecSelectCells).SpikeTimes}));
			vecSpikeRate = cellfun(@numel,{sRec.sCluster(vecSelectCells).SpikeTimes})/dblRecDur;
			matAreaWaveforms = cell2mat({sRec.sCluster(vecSelectCells).Waveform}'); %[cell x sample]
			matCABA = cell2mat({sRec.sCluster(vecSelectCells).CoordsABA}')'; %[cell x sample]
			matCUPF = matLocPerCluster(:,vecSelectCells);
			cellArea = {sRec.sCluster(vecSelectCells).Area};
			cellSelfArea = {sRec.sCluster(vecSelectCells).SelfArea};
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
			cellAggCoords2{intArea,intSubType} = cat(2,cellAggCoords2{intArea,intSubType},matCUPF);
			cellAggArea{intArea,intSubType} = cat(2,cellAggArea{intArea,intSubType},cellArea);
			cellAggSelfArea{intArea,intSubType} = cat(2,cellAggSelfArea{intArea,intSubType},cellSelfArea);
			cellAggSourceRec{intArea,intSubType} = cat(2,cellAggSourceRec{intArea,intSubType},vecSourceRec);
		end
	end
end

%% find NOT
%find bregma
vecBregma = sAtlas.Bregma;

avNot=av==(st.index(contains(st.name,'nucleus of the optic tract','ignorecase',true))+1);
vecRangeNot1 = find(sum(sum(avNot,2),3));
vecRangeNot2 = find(sum(sum(avNot,1),3));
vecRangeNot3 = find(sum(sum(avNot,1),2));
vecNot1 = (vecRangeNot1(1)-1):(vecRangeNot1(end)+1);
vecNot1(vecNot1>vecBregma(1))=[];
vecNot2 = (vecRangeNot2(1)-1):(vecRangeNot2(end)+1);
vecNot3 = (vecRangeNot3(1)-1):(vecRangeNot3(end)+1);
avNot = avNot(vecNot1,vecNot2,vecNot3);
avCenter = imfill(avNot,'holes');
SE = strel('sphere',1);
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
h = plot3(matLinesNot2(:,1), matLinesNot2(:,2), matLinesNot2(:,3), 'Color', [1 0 0 0.3]);
h.Annotation.LegendInformation.IconDisplayStyle = 'off';

cellMarker = {'x','o'};
intAreaType = 2;
%for intAreaType=1:2
	for intSubType=1:2
		
		vecRLR = cellAggSpikeRLR{intAreaType,intSubType};
		matCABA = cellAggCoords{intAreaType,intSubType};
		matCABA(1,matCABA(1,:)>vecBregma(1)) = 2*vecBregma(1) - matCABA(1,matCABA(1,:)>vecBregma(1));
		matCUPF = cellAggCoords2{intAreaType,intSubType};
		matCUPF(1,matCUPF(1,:)>vecBregma(1)) = 2*vecBregma(1) - matCUPF(1,matCUPF(1,:)>vecBregma(1));
		
		%line([matCUPF(1,:)' matCABA(1,:)']',[matCUPF(2,:)' matCABA(2,:)']',[matCUPF(3,:)' matCABA(3,:)']')
		
		
		%h= scatter3(matCABA(1,:),matCABA(2,:),matCABA(3,:),[],vecRLR,'marker',cellMarker{1});
		%cellText = cellAggSelfArea{intAreaType,intSubType};
		%text(matCABA(1,:),matCABA(2,:),matCABA(3,:),cellText);
		
		h2= scatter3(matCUPF(1,:),matCUPF(2,:),matCUPF(3,:),[],vecRLR,'marker',cellMarker{intSubType});
		%cellText = cellAggArea{intAreaType,intSubType};
		%text(matCUPF(1,:),matCUPF(2,:),matCUPF(3,:),cellText);
	end
%end
hold off
legend(cellSubjectGroups)
xlabel('ML');
ylabel('AP');
zlabel('DV');
axis equal;
hc=colorbar;
clabel(hc,'R-L ratio');
axis equal
fixfig;
grid off;
campos([1000 -50   100])
title(sprintf('Recording locations in NOT'));

%save plot
drawnow;
export_fig([strTargetPath filesep sprintf('RecLocNot.tif')]);
export_fig([strTargetPath filesep sprintf('RecLocNot.pdf')]);

%% plot shift
figure;maxfig;

for intAreaType=1:2
	for intSubType=1:2
		
		vecRLR = cellAggSpikeRLR{intAreaType,intSubType};
		matCABA = cellAggCoords{intAreaType,intSubType};
		matCABA(1,matCABA(1,:)>vecBregma(1)) = 2*vecBregma(1) - matCABA(1,matCABA(1,:)>vecBregma(1));
		matCUPF = cellAggCoords2{intAreaType,intSubType};
		matCUPF(1,matCUPF(1,:)>vecBregma(1)) = 2*vecBregma(1) - matCUPF(1,matCUPF(1,:)>vecBregma(1));
		
		vecDist = sqrt(sum((matCUPF - matCABA).^2));
		vecDepth = matCUPF(3,:);
		vecProbeLength = arrayfun(@(x) x.stereo_coordinates.ProbeLength,sExp(cellAggSourceRec{intAreaType,intSubType}));
		subplot(2,3,(intAreaType-1)*3+intSubType)%ML,AP,DV
		
		scatter3(vecDist,vecDepth,vecProbeLength);
		xlabel('ABA-UPF shift')
		ylabel('DV-Depth')
		zlabel('Probe Length')
		
	end
end

%% single axes
%swap to same hemisphere
vecML1 = cellAggCoords{2,1}(1,:);
vecML1(vecML1>vecBregma(1)) = 2*vecBregma(1) - vecML1(1,vecML1>vecBregma(1));
vecAP1 = cellAggCoords{2,1}(2,:);
vecDV1 = cellAggCoords{2,1}(3,:);
vecML2 = cellAggCoords{2,2}(1,:);
vecML2(vecML2>vecBregma(1)) = 2*vecBregma(1) - vecML2(1,vecML2>vecBregma(1));
vecAP2 = cellAggCoords{2,2}(2,:);
vecDV2 = cellAggCoords{2,2}(3,:);

figure;maxfig;
subplot(2,3,1)%DV,ML,AP
[r,p]=corr(vecML1',cellAggSpikeRLR{2,1}');
scatter(vecML1,cellAggSpikeRLR{2,1},'x');
xlabel('ML');
ylabel('R-L ratio');
title(sprintf('BL6, r(DV,RLR)=%.3f,p=%.3f',r,p));
ylim([0 1]);
fixfig;

subplot(2,3,2)%DV,ML,AP
[r,p]=corr(vecAP1',cellAggSpikeRLR{2,1}');
scatter(vecAP1,cellAggSpikeRLR{2,1},'x');
xlabel('AP');
ylabel('R-L ratio');
title(sprintf('BL6, r(ML,RLR)=%.3f,p=%.3f',r,p));
ylim([0 1]);
fixfig;

subplot(2,3,3)%DV,ML,AP
[r,p]=corr(vecDV1',cellAggSpikeRLR{2,1}');
scatter(vecDV1,cellAggSpikeRLR{2,1},'x');
xlabel('DV');
ylabel('R-L ratio');
title(sprintf('BL6, r(AP,RLR)=%.3f,p=%.3f',r,p));
ylim([0 1]);
fixfig;

subplot(2,3,4)%DV,ML,AP
[r,p]=corr(vecML2',cellAggSpikeRLR{2,2}');
scatter(vecML2,cellAggSpikeRLR{2,2});
xlabel('DV');
ylabel('R-L ratio');
title(sprintf('DBA, r(DV,RLR)=%.3f,p=%.3f',r,p));
ylim([0 1]);
fixfig;

subplot(2,3,5)%DV,ML,AP
[r,p]=corr(vecAP2',cellAggSpikeRLR{2,2}');
scatter(vecAP2,cellAggSpikeRLR{2,2});
xlabel('ML');
ylabel('R-L ratio');
title(sprintf('DBA, r(ML,RLR)=%.3f,p=%.3f',r,p));
ylim([0 1]);
fixfig;

subplot(2,3,6)%DV,ML,AP
[r,p]=corr(vecDV2',cellAggSpikeRLR{2,2}');
scatter(vecDV2,cellAggSpikeRLR{2,2});
xlabel('AP');
ylabel('R-L ratio');
title(sprintf('DBA, r(AP,RLR)=%.3f,p=%.3f',r,p));
ylim([0 1]);
fixfig;

%save plot
drawnow;
export_fig([strTargetPath filesep sprintf('RecLocNotPerAx.tif')]);
export_fig([strTargetPath filesep sprintf('RecLocNotPerAx.pdf')]);
