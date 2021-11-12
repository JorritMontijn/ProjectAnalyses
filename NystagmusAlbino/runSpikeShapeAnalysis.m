%% exploratory analysis, no proper controls

%% load data
strAllenCCFPath = 'F:\Data\AllenCCF';
[tv,av,st] = RP_LoadABA(strAllenCCFPath);
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
		
		%calculate distance to area boundary
		[cellAreas,probe_areas,probe_area_boundaries,probe_area_centers,probe_area_labels,vecChDistToAreaBoundary] = getBrainAreasPerChannel(sAP,tv,av,st);
		%get cluster depths
		[spikeAmps, vecAllSpikeDepth] = templatePositionsAmplitudes(sSpikes.temps, sSpikes.winv, sSpikes.ycoords, sSpikes.spikeTemplates, sSpikes.tempScalingAmps);
		dblProbeLength = max(sSpikes.ycoords);
		vecAllSpikeClust = sSpikes.clu;
	
		%assign cluster data
		intClustNum = numel(sAP.sCluster);
		for intClust=1:intClustNum
			intClustIdx = sAP.sCluster(intClust).IdxClust;
			intDepth = dblProbeLength-round(median(vecAllSpikeDepth(vecAllSpikeClust==intClustIdx)));
			intDominantChannel = ceil(intDepth/10);
			
			%assign
			sAP.sCluster(intClust).Waveform = matClustWaveforms(sAP.sCluster(intClust).IdxClust == vecClustIdx,:);
			sAP.sCluster(intClust).BoundDist = vecChDistToAreaBoundary(intDominantChannel);
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
dblPTT = 0.5;
dblSWT = 0.5/1000;
cellMarker = {'x','o'};
figure
maxfig;
for intSubType=1:2
	subplot(2,3,intSubType)
	hold on;
	h1=plot(dblSWT*[1 1]*1000,[0 1.2],'color', [0.7 0.7 0.7],'linestyle','-');
	h2=plot([0 1.5],dblPTT*[1 1],'color', [0.7 0.7 0.7],'linestyle','-');
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
	legend([{'';''}; cellAreaGroupsAbbr(:)],'location','best');
	title(sprintf('%s spike waveform properties',cellSubjectGroups{intSubType}));
	fixfig;
	grid off
	set(get(get(h1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
	set(get(get(h2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
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

subplot(2,3,4)
hold on
errorbar(1,dblMu_CtxBL6,vecCI_CtxBL6(1) - dblMu_CtxBL6,vecCI_CtxBL6(2) - dblMu_CtxBL6,'x','color',[0 0 1],'CapSize',20);
errorbar(2,dblMu_NotBL6,vecCI_NotBL6(1) - dblMu_NotBL6,vecCI_NotBL6(2) - dblMu_NotBL6,'x','color',[1 0 0],'CapSize',20);
hold off
ylabel(sprintf('Fraction of narrow spiking (%s +/- 95%% CI)',getGreek('mu')));
set(gca,'xtick',[1 2],'xticklabel',cellAreaGroupsAbbr(1:2));
xlim([0.5 2.5]);
ylim([0 1]);
title(sprintf('BL6,mu=%.4f,NOT,mu=%.4f; 2bino, p=%.1e',dblMu_CtxBL6,dblMu_NotBL6,p_DiffBL6));
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


subplot(2,3,5)
hold on
errorbar(1,dblMu_CtxDBA,vecCI_CtxDBA(1) - dblMu_CtxDBA,vecCI_CtxDBA(2) - dblMu_CtxDBA,'x','color',[0 0 1],'CapSize',20);
errorbar(2,dblMu_NotDBA,vecCI_NotDBA(1) - dblMu_NotDBA,vecCI_NotDBA(2) - dblMu_NotDBA,'x','color',[1 0 0],'CapSize',20);
hold off
ylabel(sprintf('Fraction of narrow spiking (%s +/- 95%% CI)',getGreek('mu')));
set(gca,'xtick',[1 2],'xticklabel',cellAreaGroupsAbbr(1:2));
xlim([0.5 2.5]);
ylim([0 1]);
title(sprintf('DBA,mu=%.4f,NOT,mu=%.4f; 2bino,p=%.1e',dblMu_CtxDBA,dblMu_NotDBA,p_DiffDBA));
fixfig;grid off;drawnow;

%save plot
drawnow;
export_fig([strTargetPath filesep sprintf('SpikeShapes.tif')]);
export_fig([strTargetPath filesep sprintf('SpikeShapes.pdf')]);
	