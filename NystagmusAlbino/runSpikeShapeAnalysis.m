%% exploratory analysis, no proper controls

%% load data
strDataPath = 'E:\DataPreProcessed';
sFiles = dir(fullpath(strDataPath,'*_AP.mat'));
if ~exist('sExp','var') || isempty(sExp) || ~isfield(sExp(1).sCluster,'Waveform')
	sExp = [];
	for intFile=1:numel(sFiles)
		fprintf('Loading %d/%d: %s [%s]\n',intFile,numel(sFiles),sFiles(intFile).name,getTime);
		sLoad = load(fullpath(sFiles(intFile).folder,sFiles(intFile).name));
		sAP = sLoad.sAP;
		if ~isfield(sAP,'sPupil')
			sAP.sPupil = [];
		end
		%load clustering data
		strSpikePath = sLoad.sAP.sSources.sClustered.folder;
		sSpikes = loadKSdir(strSpikePath);
		strImPath = sLoad.sAP.sSources.sEphysAp.folder;
		strImFile = sLoad.sAP.sSources.sEphysAp.name;
		sMetaIM = DP_ReadMeta(strImFile,strImPath);
		[vecClustIdx,matClustWaveforms] = getWaveformPerCluster(sSpikes);
		sAP.sSources.sMetaIM = sMetaIM;
		for intNeuron=1:numel(sAP.sCluster)
			sAP.sCluster(intNeuron).Waveform = matClustWaveforms(sAP.sCluster(intNeuron).IdxClust == vecClustIdx,:);
		end
		if isempty(sExp)
			sExp = sAP;
		else
			sExp(end+1) = sAP;
		end
		
	end
end

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
intUseAreaNum = numel(cellUseAreas);
vecColAlb = [0.9 0 0];
vecColBl6 = lines(1);

%% pre-allocate
cellAggSpikeDur = cell(2,2);
cellAggSpikePTR = cell(2,2);

%% run
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
		intRec=vecRunRecs(intRecIdx)
		sRec = sExp(intRec);
		strName=[sRec.sJson.subject '_' sRec.sJson.date];
		
		% split cells into areas
		cellCellsPerArea = cell(1,numel(cellUseAreas));
		cellAreasPerCluster = {sRec.sCluster.Area};
		for intArea=1:numel(cellUseAreas)
			cellCellsPerArea{intArea} = contains(cellAreasPerCluster(:),cellUseAreas{intArea},'IgnoreCase',true);
		end
		vecCellsNrPerArea = cellfun(@sum,cellCellsPerArea);
		dblSampRateIM = str2double(sAP.sSources.sMetaIM.imSampRate);
		dblSampRateNI = str2double(sAP.sSources.sMetaNI.niSampRate);
		
		%get waveform in areas
		for intArea=1:intUseAreaNum
			%include?
			indUseCells = arrayfun(@(x) x.Violations1ms < 0.25 & abs(x.NonStationarity) < 0.25,sRec.sCluster(:));
			
			%build cell vectors
			vecSelectCells = find(indUseCells(:) & cellCellsPerArea{intArea}(:));
			if isempty(vecSelectCells)
				continue;
			end
			
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
				vecSpikeDur(intNeuron) = dblTroughToPeakTime;
				vecSpikePTR(intNeuron) = abs(dblPeakVal/dblTroughVal);
			end
			cellAggSpikeDur{intArea,intSubType} = cat(2,cellAggSpikeDur{intArea,intSubType},vecSpikeDur);
			cellAggSpikePTR{intArea,intSubType} = cat(2,cellAggSpikePTR{intArea,intSubType},vecSpikePTR);
		end
	end
end

%% plot
cellMarker = {'x','o'};
figure
maxfig;
for intSubType=1:2
	subplot(2,3,intSubType)
	hold on;
	for intArea=1:intUseAreaNum
		if intArea == 1
			vecCol = [0 0 1];
		elseif intArea == 2
			vecCol = [1 0 0];
		end
		strArea = cellAreaGroupsAbbr{intArea};
		
		strSubjectType = cellSubjectGroups{intSubType};
		scatter(1000*(cellAggSpikeDur{intArea,intSubType}+(rand(size(cellAggSpikeDur{intArea,intSubType}))-0.5)/dblSampRateIM),cellAggSpikePTR{intArea,intSubType},[],vecCol,'marker','.');
	end
	hold off;
	ylim([0 1.2]);
	ylabel('Peak-to-trough ratio');
	xlabel('Spike width (ms)');
	legend(cellAreaGroupsAbbr,'location','best');
	title(sprintf('%s spike waveform properties',cellSubjectGroups{intSubType}));
	fixfig;
end

%save plot
drawnow;
export_fig([strTargetPath filesep sprintf('SpikeShapes.tif')]);
export_fig([strTargetPath filesep sprintf('SpikeShapes.pdf')]);
	