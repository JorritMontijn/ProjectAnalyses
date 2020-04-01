%% load spiking data & plot tuning curves

%% set recording
%close all;
clear all;
cellUniqueAreas = {'V1','SC','Poisson','Retina','GCaMP'};
strArea = 'SC'; %V1, SC, NOT
strDataMasterPath = 'D:\Data\Processed\ePhys\';
strDataTargetPath = 'D:\Data\ResultsOriMetric\Data\';
strFigPath = 'D:\Data\ResultsOriMetric\TuningCurves\';
intMakePlots = 2; %0=none, 1=normal plot, 2=including raster
vecRunTypes = [1 2];
intResampleNum = 100;
boolSave = false;

%set var
for intArea=4%3:numel(cellUniqueAreas)
	strArea = cellUniqueAreas{intArea};
	%reset vars
	clearvars -except boolSave intArea strArea cellUniqueAreas strDataMasterPath strDataTargetPath strFigPath intMakePlots vecRunTypes
	
for intRunType=vecRunTypes
if intRunType == 1
	strRunType = strArea;
elseif intRunType ==2
	strRunType = [strArea '-Rand'];
end
%% load data
if contains(strRunType,'V1') || contains(strRunType,'SC')
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
	intNeurons = 133;
elseif contains(strRunType,'GCaMP')
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
vecResamples = 50;%[2 5:5:50];%100];
vecPlotCells=false(1,intNeurons);
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
	cellZ = cell(1,intNeurons);
	cellArea = cell(1,intNeurons);
	cellDeriv = cell(1,intNeurons);
	cellInterpT = cell(1,intNeurons);
	
	%% analyze
	for intNeuron=[1:intNeurons]%31
		%% message
		if toc(hTic) > 5
			fprintf('Processing neuron %d/%d [%s]\n',intNeuron,intNeurons,getTime);
			hTic=tic;
		end
		clear vecTrialStarts;
		%% load or generate data
		if contains(strRunType,'V1') || contains(strRunType,'SC')
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
		elseif contains(strRunType,'GCaMP')
			%% generate
			strRecIdx = ses.session;
			strMouse = getFlankedBy(ses.experiment,'_','.mat','last');
			strBlock = sprintf('xyt%02d',ses.recording);
			strDate = ses.session;
			intSU = intNeuron;
			intClust = intNeuron;
			
			%get data
			vecOn = ses.structStim.FrameOn;
			vecdFoF = cellData{intNeuron};
			vecTraceT = (1:numel(vecdFoF)) ./ ses.samplingFreq;
			vecStimOnTime = vecOn(:) ./ ses.samplingFreq;
			vecStimOffTime = vecStimOnTime + median(diff(vecStimOnTime))/2;
			
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
		if contains(strRunType,'Rand')
			vecJitter = 2*median(diff(vecStimOnTime))*rand(size(vecStimOnTime));
			matStimOnOff = bsxfun(@plus,vecTrialStarts,vecJitter);
		else
			matStimOnOff = vecTrialStarts;
		end
		close;
		if exist('vecTraceAct','var') && ~isempty(vecTraceAct) %calcium imaging
			if intResampleNum == vecResamples(end) && vecPlotCells(intNeuron)
				[dblZeta,sOptOut] = getTraceZeta(vecTraceT,vecTraceAct,matStimOnOff,2,[],intResampleNum,0);
			else
				[dblZeta,sOptOut] = getTraceZeta(vecTraceT,vecTraceAct,matStimOnOff,intMakePlots,[],intResampleNum,0);
			end
			intSpikeNum = mean(vecTraceAct) + std(vecTraceAct);
		else
			if intResampleNum == vecResamples(end) && vecPlotCells(intNeuron)
				[dblZeta,sOptOut] = getZeta(vecSpikeTimes,matStimOnOff,2,[],intResampleNum,0);
			else
				[dblZeta,sOptOut] = getZeta(vecSpikeTimes,matStimOnOff,intMakePlots,[],intResampleNum,0);
			end
			intSpikeNum = numel(vecSpikeTimes);
		end
		
		% assign data
		dblHzD = sOptOut.dblHzD;
		vecNumSpikes(intNeuron) = intSpikeNum;
		vecZeta(intNeuron) = dblZeta;
		vecP(intNeuron) = sOptOut.dblP;
		vecHzD(intNeuron) = dblHzD;
		vecHzP(intNeuron) = sOptOut.dblHzP;
		cellZ{intNeuron} = sOptOut.vecZ;
		cellArea{intNeuron} = strArea;
		
		%pause
		%continue;
		%% build vector for cells to plot
		if 0%intMakePlots || (intResampleNum == vecResamples(end-1) && sOptOut.dblP < 0.05 && sOptOut.dblHzP > 0.05)
			vecPlotCells(intNeuron) = true;
		end
		%% save plot
		if intMakePlots || (intResampleNum == vecResamples(end) && vecPlotCells(intNeuron))% && ~(exist('vecTraceAct','var') && ~isempty(vecTraceAct)))
			%plot
			subplot(2,3,6)
			
			%strRecIdx = sThisNeuron.RecIdx;
			%strMouse = sThisNeuron.Mouse;
			%strBlock = sThisNeuron.Block;
			%strArea = sThisNeuron.Area;
			%strDate = sThisNeuron.Date;
			%intSU = sThisNeuron.IdxSU;
			%intClust = sThisNeuron.IdxClust;
			
			strFileName = sprintf('%s%sB%sSU%dC%d',strArea,strDate,strBlock,intSU,intClust);
			export_fig([strFigPath strFileName 'Zeta' strRunType 'Resamp' num2str(intResampleNum) '.tif']);
			print(gcf, '-dpdf', [strFigPath strFileName 'Zeta' strRunType 'Resamp' num2str(intResampleNum) '.pdf']);
		end
		%pause
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
	save([strDataTargetPath 'ZetaData' strRunType 'Resamp' num2str(intResampleNum) '.mat' ],...
		'vecNumSpikes','vecZeta','vecP','vecHzD','vecHzP','cellZ','cellInterpT','cellArea','cellDeriv');	
	end
end
end
end