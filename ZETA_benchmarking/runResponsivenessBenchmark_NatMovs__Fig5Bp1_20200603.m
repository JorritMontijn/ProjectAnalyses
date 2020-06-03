%% load spiking data & plot tuning curves

%% set recording
%close all;
clear all;
cellUniqueAreas = {...
	'V1',...Area 1
	'SC',...Area 2
	'Poisson',...Area 3
	'Retina',...Area 4
	'CaNM',...Area 5
	'CaDG',...Area 6
	'lateral geniculate',...Area 7
	'Primary visual',...Area 8
	'Lateral posterior nucleus',...Area 9
	'Anterior pretectal nucleus',...Area 10
	'Nucleus of the optic tract',...Area 11
	'Superior colliculus',...Area 12
	'Anteromedial visual',...Area 13
	'posteromedial visual',...Area 14
	'Anterolateral visual',...Area 15
	'Lateral visual',...Area 16
	'Rostrolateral area',...Area 17
	'Anterior area',...Area 18
	'Subiculum',...Area 19
	'Field CA1',...Area 20
	'Field CA2',...Area 21
	'Field CA3',...Area 22
	'Dentate gyrus',...Area 23
	'Retrosplenial'...Area 24
	};

strDisk = 'D:';
strDataMasterPath = [strDisk '\Data\Processed\ePhys\'];
strDataTargetPath = [strDisk '\Data\Processed\ZETA\NatMovs\'];
strFigPath = [strDisk '\Data\Results\ZETA\NatMovs\'];
intMakePlots =0; %0=none, 1=normal plot, 2=including raster
vecRandTypes = [1 2];%1=normal,2=rand
vecRestrictRange = [0 inf];
boolSave = true;
intResampleNum = 100;%10:10:90;%[10:10:100];
vecRunAreas = [7:16];%[7:24];%[1:4];%1:6;%1:5;
cellRunStim = {'','RunDriftingGratings','RunNaturalMovie'};
vecRunStim = 3;%2:3;
cellRepStr = {...
	'RunDriftingGratings','-DG';...
	'RunNaturalMovie','-NM';...
	'lateral geniculate','LGN';...
	'Primary visual','V1';...
	'Lateral posterior nucleus','LP';...
	'Anterior pretectal nucleus','APN';...
	'Nucleus of the optic tract','NOT';...
	'Superior colliculus','SC';...
	'Anteromedial visual','AM';...
	'posteromedial visual','PM';...
	'Anterolateral visual','AL';...
	'Lateral visual','L';...
	'Rostrolateral area','RL';...
	'Anterior area','A';...
	'Subiculum','H-Sub';...
	'Field CA1','H-CA1';...
	'Field CA2','H-CA2';...
	'Field CA3','H-CA3';...
	'Dentate gyrus','H-DG';...
	'Retrosplenial','RetSpl';...
};

%set var
for intArea=vecRunAreas
	if intArea < 7
		vecUseRunStim = 1;
	else
		vecUseRunStim = vecRunStim;
	end
for intRunStim=vecUseRunStim
for intRandType=vecRandTypes
	%reset vars
	clearvars -except vecRestrictRange intResampleNum cellRepStr intRandType vecRandTypes intRunStim vecRunStim cellRunStim intArea vecRunAreas cellUniqueAreas boolSave strDataMasterPath strDataTargetPath strFigPath intMakePlots vecRunTypes
	strArea = cellUniqueAreas{intArea};
	strRunStim = cellRunStim{intRunStim};

if intRandType == 1
	strRunType = strArea;
elseif intRandType ==2
	strRunType = [strArea '-Rand'];
end

%% load data
if contains(strRunType,cellUniqueAreas(7:end),'IgnoreCase',true)
	%% find data
	strDataSourcePath = 'D:\Data\Processed\Neuropixels\';
	sFiles = dir([strDataSourcePath '*.mat']);
	cellFiles = {sFiles(:).name}';
	strName = replace([lower(strArea) strRunStim],lower(cellRepStr(:,1)),cellRepStr(:,2));
	
	%% go through files
	clear sAggStim;
	clear sAggNeuron;
	intNeurons = 0;
	for intFile=1:numel(cellFiles)
		%% load
		fprintf('Loading %s [%s]\n',cellFiles{intFile},getTime);
		sLoad = load([strDataSourcePath cellFiles{intFile}]);
		sAP = sLoad.sAP;
		intNewFile = 1;
		%check if neuron is in target area
		for intClust=1:numel(sAP.sCluster)
			strClustArea = sAP.sCluster(intClust).Area;
			if ~isempty(strClustArea) && contains(strClustArea,strArea,'IgnoreCase',true) && (sAP.sCluster(intClust).KilosortGood || sAP.sCluster(intClust).Contamination < 0.1)
				%% aggregate data
				%check if stim type is present
				indUseStims = ismember(cellfun(@(x) x.structEP.strFile,sAP.cellStim,'uniformoutput',false),strRunStim);
				if isempty(indUseStims) || ~any(indUseStims)
					continue;
				end
				%add data
				if intNeurons == 0
					intNewFile = 0;
					sAggNeuron(1) = sAP.sCluster(intClust);
					sAggStim(1).cellStim = sAP.cellStim(indUseStims);
					sAggStim(1).Rec = sAggNeuron(end).Rec;
				elseif ~isempty(indUseStims) && any(indUseStims)
					sAggNeuron(end+1) = sAP.sCluster(intClust);
				end
				if intNewFile
					sAggStim(end+1).cellStim = sAP.cellStim(indUseStims);
					sAggStim(end).Rec = sAggNeuron(end).Rec;
					intNewFile = 0;
				end
				intNeurons = intNeurons + 1;
			end
		end
	end
	if ~exist('sAggStim','var')
		continue;
	end
	
	cellRecIdx = {sAggStim.Rec};
	fprintf('Found %d cells from %d recordings in "%s" [%s]\n',intNeurons,numel(cellRecIdx),strRunType,getTime);
	
elseif contains(strRunType,'CaNM')
	%% find data
	strDataSourcePath = 'D:\Data\Processed\imagingGCaMP\';
	sFiles = dir([strDataSourcePath '20150511xyt02_ses.mat']);
	cellFiles = {sFiles(:).name}';
	sLoad = load([strDataSourcePath cellFiles{1}]);
	ses = sLoad.ses;
	intNeurons = numel(ses.neuron);
	vecOn = ses.structStim.FrameOn;
	cellData = {ses.neuron(:).dFoF};
end

	%% pre-allocate output variables
	vecNumSpikes = nan(1,intNeurons);
	vecZetaP = nan(1,intNeurons);
	vecP = nan(1,intNeurons);
	vecHzD = nan(1,intNeurons);
	vecAnovaP = nan(1,intNeurons);
	vecPeakT = nan(1,intNeurons);
	vecPeakW = nan(1,intNeurons);
	cellD = cell(1,intNeurons);
	cellArea = cell(1,intNeurons);
	cellDeriv = cell(1,intNeurons);
	cellInterpT = cell(1,intNeurons);
	cellPeakT = cell(1,intNeurons);
	cellPeakI = cell(1,intNeurons);

	%% analyze
	hTic=tic;
	for intNeuron=[1:intNeurons]%31
		%% message
		if toc(hTic) > 5
			fprintf('Processing neuron %d/%d [%s]\n',intNeuron,intNeurons,getTime);
			hTic=tic;
		end
		clear vecTrialStarts;
		%% load or generate data
		if contains(strRunType,cellUniqueAreas(7:end),'IgnoreCase',true)
			%% get neuronal data
			sThisNeuron = sAggNeuron(intNeuron);
			vecSpikeTimes = sThisNeuron.SpikeTimes;
			strRecIdx = sThisNeuron.Rec;
			strMouse = sThisNeuron.Mouse;
			strBlock = '';
			strArea = strName;
			strDate = sThisNeuron.Exp;
			intSU = sThisNeuron.Cluster;
			intClust = sThisNeuron.IdxClust;

			%% get matching recording data
			sThisRec = sAggStim(strcmpi(strRecIdx,cellRecIdx));
			vecStimOnTime = [];
			vecStimOffTime = [];
			vecFrameRate = [];
			for intRec=1:numel(sThisRec.cellStim)
				vecStimOnTime = cat(2,vecStimOnTime,sThisRec.cellStim{intRec}.structEP.vecStimOnTime);
				vecStimOffTime = cat(2,vecStimOffTime,sThisRec.cellStim{intRec}.structEP.vecStimOffTime);
				vecFrameRate = cat(2,vecFrameRate,sThisRec.cellStim{intRec}.structEP.FrameRate);
			end
			
			vecTrialStarts = [];
			vecTrialStarts(:,1) = vecStimOnTime;
			vecTrialStarts(:,2) = vecStimOffTime;
		elseif contains(strRunType,'CaNM')
			%% generate
			strRecIdx = ses.session;
			strMouse = getFlankedBy(ses.experiment,'_','.mat','last');
			strBlock = sprintf('xyt%02d',ses.recording);
			strDate = ses.session;
			intSU = intNeuron;
			intClust = intNeuron;

			%get data
			vecOn = ses.structStim.FrameOn;
			vecOff = ses.structStim.FrameOff;
			vecdFoF = cellData{intNeuron};
			vecTraceT = (1:numel(vecdFoF)) ./ ses.samplingFreq;
			vecStimOnTime = vecOn(:) ./ ses.samplingFreq;
			if vecOff(1) == vecOn(2)
				vecStimOffTime = vecStimOnTime + median(diff(vecStimOnTime))/2;
			else
				vecStimOffTime = vecOff(:) ./ ses.samplingFreq;
			end
			%put in output
			vecTraceAct = vecdFoF;
			vecTrialStarts(:,1) = vecStimOnTime;
			vecTrialStarts(:,2) = vecStimOffTime;
		end

		%% get visual responsiveness
		%get trial dur
		dblUseMaxDur = round(median(diff(vecTrialStarts(:,1)))*2)/2;
		%set derivative params
		if contains(strRunType,'Rand')
			dblDur = dblUseMaxDur;
			vecJitter = 2*dblDur*rand([numel(vecStimOnTime) 1])-dblDur;
			matEventTimes = bsxfun(@plus,vecTrialStarts,vecJitter);
		else
			matEventTimes = vecTrialStarts;
		end
		close;close;
		if exist('vecTraceAct','var') && ~isempty(vecTraceAct) %calcium imaging
			%continue;
			[dblZetaP,vecLatencies,sZETA] = getTraceZeta(vecTraceT,vecTraceAct,matEventTimes,dblUseMaxDur,intResampleNum,intMakePlots,0,vecRestrictRange);
			intSpikeNum = mean(vecTraceAct) + std(vecTraceAct);
			%unpack
			vecT = sZETA.vecRefT;
		else
			[dblZetaP,vecLatencies,sZETA] = getZeta(vecSpikeTimes,matEventTimes,dblUseMaxDur,intResampleNum,intMakePlots,0,vecRestrictRange);
			
			intSpikeNum = numel(vecSpikeTimes);
			%unpack
			vecT = sZETA.vecSpikeT;
		end
		
		%% get bin-wise approach
		%get data
		dblFrameDur = 1/mean(vecFrameRate);
		dblStimDur = median(diff(vecStimOnTime));
		vecBinEdges = 0:dblFrameDur:20;
		vecBinCenters = vecBinEdges(2:end)-dblFrameDur/2;
		intBins = numel(vecBinCenters);
		intTrials = numel(vecStimOnTime);
		matResp = nan(intTrials,intBins);
		for intTrial=1:intTrials
			matResp(intTrial,:) = histcounts(vecSpikeTimes,vecBinEdges+vecStimOnTime(intTrial));
		end
		
		%test
		[dblAnovaP,tbl,stats] = anova1(matResp,[],'off');
		
		%% save output
		% assign data
		vecNumSpikes(intNeuron) = intSpikeNum;
		vecZetaP(intNeuron) = dblZetaP;
		vecAnovaP(intNeuron) = dblAnovaP;
		cellArea{intNeuron} = strArea;
		
		%continue;
		%% build vector for cells to plot
		%% save plot
		if intMakePlots% && ~(exist('vecTraceAct','var') && ~isempty(vecTraceAct)))
			%plot
			vecHandles = get(gcf,'children');
			ptrFirstSubplot = vecHandles(find(contains(get(vecHandles,'type'),'axes'),1,'last'));
			axes(ptrFirstSubplot);
			title(sprintf('%s-N%d, %s-%s,U%d/C%d,p=%.3e',strArea,intNeuron,strDate,strBlock,intSU,intClust,dblAnovaP));
			drawnow;

			strFileName = sprintf('NatMov%s%s-N%d%sSU%dC%d',strArea,strDate,intNeuron,strBlock,intSU,intClust);
			export_fig([strFigPath strFileName 'Zeta' strArea 'Resamp' num2str(intResampleNum) '.tif']);
			%print(gcf, '-dpdf', [strFigPath strFileName 'Zeta' strRunType 'Resamp' num2str(intResampleNum) '.pdf']);
			export_fig([strFigPath strFileName 'Zeta_EF_' strArea 'Resamp' num2str(intResampleNum) '.pdf']);
			%pause
		end
	end
	%save
	%vecNumSpikes = nan(1,intNeurons);
	%vecZeta = nan(1,intNeurons);
	%vecP = nan(1,intNeurons);
	%vecHzD = nan(1,intNeurons);
	%vecHzP = nan(1,intNeurons);
	%cellZeta = cell(1,intNeurons);
	%cellArea = cell(1,intNeurons);
	if boolSave
		save([strDataTargetPath 'ZetaDataNatMov' strRunType strRunStim '.mat' ],...
			'vecNumSpikes','vecZetaP','vecAnovaP','cellArea');
	end
end
end
end