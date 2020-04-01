%% load spiking data & plot tuning curves

%% set recording
%close all;
clear all;
cellUniqueAreas = {'V1','SC','Poisson','Retina','CaNM','CaDG','lateral geniculate','Primary visual'};
strDataMasterPath = 'D:\Data\Processed\ePhys\';
strDataTargetPath = 'D:\Data\Results\OriMetric\Data\';
strFigPath = 'D:\Data\Results\OriMetric\TuningCurves\';
intMakePlots =0; %0=none, 1=normal plot, 2=including raster
vecRunTypes = [1 2];%[1 2];
boolSave = true;
vecResamples = [10:10:100];
vecRunAreas = 7;%[1:4];%1:6;%1:5;
cellRunStim = {'','RunDriftingGratings','RunNaturalMovie'};
vecRunStim = 2:3;

%set var
for intRunStim=vecRunStim
for intArea=vecRunAreas
for intRunType=vecRunTypes
	%reset vars
	clearvars -except intRunType vecRunTypes intRunStim vecRunStim cellRunStim intArea vecRunAreas cellUniqueAreas boolSave vecResamples strDataMasterPath strDataTargetPath strFigPath intMakePlots vecRunTypes
	strArea = cellUniqueAreas{intArea};
	strRunStim = cellRunStim{intRunStim};

if intRunType == 1
	strRunType = strArea;
elseif intRunType ==2
	strRunType = [strArea '-Rand'];
end

%% load data
if contains(strRunType,'lateral geniculate') || contains(strRunType,'Primary visual')
	%% find data
	strDataSourcePath = 'D:\Data\Processed\Neuropixels\';
	sFiles = dir([strDataSourcePath '*.mat']);
	cellFiles = {sFiles(:).name}';

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
			if ~isempty(strClustArea)
				%strClustArea
			end
			if ~isempty(strClustArea) && contains(strClustArea,strArea) && (sAP.sCluster(intClust).KilosortGood || sAP.sCluster(intClust).Contamination < 0.1)
				%% aggregate data
				if intNeurons == 0
					intNewFile = 0;
					sAggNeuron(1) = sAP.sCluster(intClust);
					indUseStims = ismember(cellfun(@(x) x.structEP.strFile,sAP.cellStim,'uniformoutput',false),strRunStim);
					sAggStim(1).cellStim = sAP.cellStim(indUseStims);
					sAggStim(1).Rec = sAggNeuron(end).Rec;
				else
					sAggNeuron(end+1) = sAP.sCluster(intClust);
				end
				if intNewFile
					indUseStims = ismember(cellfun(@(x) x.structEP.strFile,sAP.cellStim,'uniformoutput',false),strRunStim);
					sAggStim(end+1).cellStim = sAP.cellStim(indUseStims);
					sAggStim(end).Rec = sAggNeuron(end).Rec;
					intNewFile = 0;
				end
				intNeurons = intNeurons + 1;
			end
		end
	end
	cellRecIdx = {sAggStim.Rec};
	fprintf('Found %d cells from %d recordings in "%s" [%s]\n',intNeurons,numel(cellRecIdx),strRunType,getTime);
	
elseif contains(strRunType,'V1') || contains(strRunType,'SC')
	%% find data
	strDataSourcePath = [strDataMasterPath 'DriftingGratings\'];
	sFiles = dir([strDataSourcePath '*' strArea '_*.mat']);
	cellFiles = {sFiles(:).name}';

	%% go through files
	sAggRec = struct;
	sAggNeuron = struct;
	for intFile=1:numel(cellFiles)
		%% load
		sLoad = load([strDataSourcePath cellFiles{intFile}]);
		sNeuron = sLoad.sNeuron;
		sRecording = sLoad.sRecording;

		%% aggregate data
		if intFile == 1
			sAggRec= sRecording;
			sAggNeuron = sNeuron;
		else
			sAggRec(end+1) = sRecording;
			sAggNeuron((end+1):(end+numel(sNeuron))) = sNeuron;
		end
	end
	intNeurons = numel(sAggNeuron);
	cellRecIdx = {sAggRec.RecIdx};
elseif contains(strRunType,'Retina')
	%% find data
	strDataSourcePath = strDataMasterPath;
	sFiles = dir([strDataSourcePath '*' strArea '_*.mat']);
	cellFiles = {sFiles(:).name}';
	sLoad = load([strDataSourcePath cellFiles{1}]);
	cellRetData = sLoad.LightFlash_For_J;
	intNeurons = size(cellRetData,1)-2;
elseif contains(strRunType,'Poisson')
	intNeurons = 100;
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
elseif contains(strRunType,'CaDG')
	strDataSourcePath = 'D:\Data\Processed\imagingGCaMP\';
	sFiles = dir([strDataSourcePath '20150511xyt01_ses.mat']);
	cellFiles = {sFiles(:).name}';
	sLoad = load([strDataSourcePath cellFiles{1}]);
	ses = sLoad.ses;
	intNeurons = numel(ses.neuron);
	vecOn = ses.structStim.FrameOn;
	cellData = {ses.neuron(:).dFoF};
end

for intResampleIdx = 1:numel(vecResamples)
	intResampleNum = vecResamples(intResampleIdx);
	%% message
	fprintf('Processing %s, resampling %d (%d/%d) [%s]\n',strRunType,intResampleNum,intResampleIdx,numel(vecResamples),getTime);
	hTic=tic;

	%% pre-allocate output variables
	vecNumSpikes = nan(1,intNeurons);
	vecZeta = nan(1,intNeurons);
	vecP = nan(1,intNeurons);
	vecHzD = nan(1,intNeurons);
	vecHzP = nan(1,intNeurons);
	vecPeakT = nan(1,intNeurons);
	vecPeakW = nan(1,intNeurons);
	cellD = cell(1,intNeurons);
	cellArea = cell(1,intNeurons);
	cellDeriv = cell(1,intNeurons);
	cellInterpT = cell(1,intNeurons);
	cellPeakT = cell(1,intNeurons);
	cellPeakI = cell(1,intNeurons);

	%% analyze
	for intNeuron=[1:intNeurons]%31
		%% message
		if toc(hTic) > 5
			fprintf('Processing neuron %d/%d [%s]\n',intNeuron,intNeurons,getTime);
			hTic=tic;
		end
		clear vecTrialStarts;
		%% load or generate data
		if contains(strRunType,'lateral geniculate') || contains(strRunType,'Primary visual')
			%% get neuronal data
			sThisNeuron = sAggNeuron(intNeuron);
			vecSpikeTimes = sThisNeuron.SpikeTimes;
			strRecIdx = sThisNeuron.Rec;
			strMouse = sThisNeuron.Mouse;
			strArea = sThisNeuron.Area;
			strDate = sThisNeuron.Date;
			intSU = sThisNeuron.Cluster;
			intClust = sThisNeuron.IdxClust;

			%% get matching recording data
			sThisRec = sAggStim(strcmpi(strRecIdx,cellRecIdx));
			vecStimOnTime = [];
			vecStimOffTime = [];
			for intRec=1:numel(sThisRec.cellStim)
				vecStimOnTime = cat(2,vecStimOnTime,sThisRec.cellStim{intRec}.structEP.vecStimOnTime);
				vecStimOffTime = cat(2,vecStimOffTime,sThisRec.cellStim{intRec}.structEP.vecStimOffTime);
			end
			
			vecTrialStarts = [];
			vecTrialStarts(:,1) = vecStimOnTime;
			vecTrialStarts(:,2) = vecStimOffTime;
		elseif contains(strRunType,'V1') || contains(strRunType,'SC')
			%% get neuronal data
			sThisNeuron = sAggNeuron(intNeuron);
			vecSpikeTimes = sThisNeuron.SpikeTimes;
			strRecIdx = sThisNeuron.RecIdx;
			strMouse = sThisNeuron.Mouse;
			strBlock = sThisNeuron.Block;
			strArea = sThisNeuron.Area;
			strDate = sThisNeuron.Date;
			intSU = sThisNeuron.IdxSU;
			intClust = sThisNeuron.IdxClust;

			%% get matching recording data
			sThisRec = sAggRec(strcmpi(strRecIdx,cellRecIdx));
			vecStimOnTime = sThisRec.vecStimOnTime;
			vecStimOffTime = sThisRec.vecStimOffTime;
			vecStimOriDegrees = sThisRec.vecStimOriDegrees;
			vecStimOriRads = deg2rad(vecStimOriDegrees);
			vecEyeTimestamps = sThisRec.vecEyeTimestamps;
			matEyeData = sThisRec.matEyeData;
			vecTrialStarts = [];
			vecTrialStarts(:,1) = vecStimOnTime;
			vecTrialStarts(:,2) = vecStimOffTime;
		elseif contains(strRunType,'Retina')
			%% generate
			strCellID = cellRetData{intNeuron+2,1};
			cellStr = strsplit(strCellID,'_');
			cellStr2 = strsplit(cellStr{2},'-');
			strRecIdx = cellStr2{1};
			strMouse = 'Kamermans';
			strBlock = cellStr2{1};
			strDate = cellStr{1};
			intSU = intNeuron;
			intClust = str2double(cellStr2{2}(2));

			%get data
			cellData = cellRetData(3:end,2);
			matN = cellData{intNeuron};

			%get trial timing
			intOn = find(cellRetData{3,3}(:,2)==1,1);
			dblStimOn = cellRetData{3,3}(intOn,1);
			intOff = find(cellRetData{3,3}(intOn:end,2)==0,1) + intOn - 1;
			dblStimOff = cellRetData{3,3}(intOff,1);
			dblTrialDur=cellRetData{3,3}(end-1,1);
			%build vectors
			vecAllTrialStarts = (matN(:,2)-1)*dblTrialDur;
			vecTrialStartTime = ((1:max(matN(:,2)))-1)'*dblTrialDur;
			vecStimOnTime = vecTrialStartTime + dblStimOn;
			vecStimOffTime = vecTrialStartTime + dblStimOff;

			%put in output
			vecSpikeTimes = matN(:,1) + vecAllTrialStarts;
			vecTrialStarts(:,1) = vecStimOnTime;
			vecTrialStarts(:,2) = vecStimOffTime;
		elseif contains(strRunType,'CaNM') || contains(strRunType,'CaDG')
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
		elseif contains(strRunType,'Poisson')
			%% generate
			strRecIdx = 'x';
			strMouse = 'Artificial';
			strBlock = '1';
			strDate = getDate();
			intSU = intNeuron;
			intClust = intNeuron;

			%set parameters
			dblBaseRate = exprnd(5);
			dblPrefRate = dblBaseRate+exprnd(20);
			dblKappa = rand(1)*5+5;
			vecTrialAngles=repmat([0:45:360],[1 20]);
			dblTrialDur=2;
			vecStimOnTime = dblTrialDur*(1:numel(vecTrialAngles))';
			vecStimOffTime = vecStimOnTime + 1;

			vecTrialStarts(:,1) = vecStimOnTime;
			vecTrialStarts(:,2) = vecStimOffTime;
			[vecSpikeTimes,dblPrefOri] = getGeneratedSpikingData(vecTrialAngles,vecTrialStarts,dblBaseRate,dblPrefRate,dblKappa,true);
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
			if intResampleNum == vecResamples(end)% && vecPlotCells(intNeuron)
				[dblZeta,vecLatencies,sZETA,sMSD] = getTraceZeta(vecTraceT,vecTraceAct,matEventTimes,dblUseMaxDur,intResampleNum,intMakePlots);
			else
				[dblZeta,vecLatencies,sZETA,sMSD] = getTraceZeta(vecTraceT,vecTraceAct,matEventTimes,dblUseMaxDur,intResampleNum,0);
			end
			intSpikeNum = mean(vecTraceAct) + std(vecTraceAct);
			%unpack
			vecT = sZETA.vecRefT;
		else
			if intResampleNum == vecResamples(end)% && vecPlotCells(intNeuron)
				[dblZeta,vecLatencies,sZETA,sMSD] = getZeta(vecSpikeTimes,matEventTimes,dblUseMaxDur,intResampleNum,intMakePlots);
			else
				[dblZeta,vecLatencies,sZETA,sMSD] = getZeta(vecSpikeTimes,matEventTimes,dblUseMaxDur,intResampleNum,0);
			end
			
			intSpikeNum = numel(vecSpikeTimes);
			%unpack
			vecT = sZETA.vecSpikeT;
		end
		
		%unpack
		if isempty(sMSD)
			vecMSD = [];
			vecTime = [];
			vecIdx = [];
			dblPeakTime = nan;
			dblPeakWidth = nan;
		else
			vecMSD = sMSD.vecRate;
			vecTime = sMSD.vecT;
			vecIdx = sMSD.vecPeakStartStopIdx;
			dblOnset = sMSD.dblOnset;
			dblPeakTime = sMSD.dblPeakTime;
			dblPeakWidth = sMSD.dblPeakWidth;
		end 
		
		% assign data
		dblMeanD = sZETA.dblMeanD;
		vecNumSpikes(intNeuron) = intSpikeNum;
		vecZeta(intNeuron) = dblZeta;
		vecP(intNeuron) = sZETA.dblP;
		vecHzD(intNeuron) = dblMeanD;
		vecHzP(intNeuron) = sZETA.dblMeanP;
		vecPeakT(intNeuron) = dblPeakTime;
		vecPeakW(intNeuron) = dblPeakWidth;
		cellD{intNeuron} = sZETA.vecD;
		cellInterpT{intNeuron} = vecT;
		cellArea{intNeuron} = strArea;
		cellDeriv{intNeuron} = vecMSD;
		cellPeakT{intNeuron} = vecTime;%vecPeakTimes;
		cellPeakI{intNeuron} = vecIdx;%vecPeakValues;

		%continue;
		%% build vector for cells to plot
		%% save plot
		if intMakePlots && (intResampleNum == vecResamples(end))% && ~(exist('vecTraceAct','var') && ~isempty(vecTraceAct)))
			%plot
			vecHandles = get(gcf,'children');
			ptrFirstSubplot = vecHandles(find(contains(get(vecHandles,'type'),'axes'),1,'last'));
			axes(ptrFirstSubplot);
			title(sprintf('%s-N%d, %s-%s,U%d/C%d',strArea,intNeuron,strDate,strBlock,intSU,intClust));
			drawnow;

			strFileName = sprintf('%s%s-%sSU%dC%d',strArea,strDate,strBlock,intSU,intClust);
			export_fig([strFigPath strFileName 'Zeta' strRunType 'Resamp' num2str(intResampleNum) '.tif']);
			%print(gcf, '-dpdf', [strFigPath strFileName 'Zeta' strRunType 'Resamp' num2str(intResampleNum) '.pdf']);
			export_fig([strFigPath strFileName 'Zeta_EF_' strRunType 'Resamp' num2str(intResampleNum) '.pdf']);
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
		save([strDataTargetPath 'ZetaDataMSD' strRunType strRunStim 'Resamp' num2str(intResampleNum) '.mat' ],...
			'vecPeakT','vecPeakW','vecNumSpikes','vecZeta','vecP','vecHzD','vecHzP','cellD','cellInterpT','cellArea','cellDeriv','cellPeakT','cellPeakI');
	end
end
end
end
end