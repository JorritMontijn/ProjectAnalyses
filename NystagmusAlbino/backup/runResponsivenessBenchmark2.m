%% load spiking data & plot tuning curves

%% set recording
close all;
clear all;
cellUniqueAreas = {'V1','SC','Poisson','Retina','CaNM','CaDG'};
strArea = 'V1'; %V1, SC, NOT
strDataMasterPath = 'D:\Data\Processed\ePhys\';
strDataTargetPath = 'D:\Data\ResultsOriMetric\Data\';
strFigPath = 'D:\Data\ResultsOriMetric\TuningCurves\';
intMakePlots = 0; %0=none, 1=normal plot, 2=including raster
vecRunTypes = [1 2];
intResampleNum = 100;
boolSave = true;

%set var
for intArea=5:numel(cellUniqueAreas)
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
vecResamples = [2 5:5:50];%100];
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
	cellPeakT = cell(1,intNeurons);
	cellPeakV = cell(1,intNeurons);

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
		%set derivative params
		intPeaks = 3;
		intSmoothSd = 3;
		dblBinSize = (1/1000);
		if contains(strRunType,'Rand')
			vecJitter = 2*median(diff(vecStimOnTime))*rand(size(vecStimOnTime));
			matStimOnOff = bsxfun(@plus,vecTrialStarts,vecJitter);
		else
			matStimOnOff = vecTrialStarts;
		end
		close;close;
		if exist('vecTraceAct','var') && ~isempty(vecTraceAct) %calcium imaging
			if intResampleNum == vecResamples(end) && vecPlotCells(intNeuron)
				[dblZeta,sOptOut] = getTraceZeta(vecTraceT,vecTraceAct,matStimOnOff,2,[],intResampleNum,0);
			else
				[dblZeta,sOptOut] = getTraceZeta(vecTraceT,vecTraceAct,matStimOnOff,intMakePlots,[],intResampleNum,0);
			end
			intSpikeNum = mean(vecTraceAct) + std(vecTraceAct);
			%get deriv
			vecSmoothDeriv = diff(sOptOut.vecRealDiff);
			%unpack
			vecT = sOptOut.vecRefT;
			%find peaks
			[vecPeakValues,vecPeakLocs,vecPeakWidths,vecPeakProminences] = findpeaks(vecSmoothDeriv);
			if numel(vecPeakValues) < intPeaks
				intPeaks = numel(vecPeakValues);
			end
			%sort by prominence
			[vecPeakProminences,vecReorder] = sort(vecPeakProminences,'descend');
			vecPeakTimes = vecT(vecPeakLocs);
			vecPeakTimes = vecPeakTimes(vecReorder);
			vecPeakWidths = vecPeakWidths(vecReorder);
			vecPeakValues = vecPeakValues(vecReorder);
			%select top3
			vecPeakTimes = vecPeakTimes(1:intPeaks);
			vecPeakWidths = vecPeakWidths(1:intPeaks);
			vecPeakValues = vecPeakValues(1:intPeaks);
			vecPeakProminences = vecPeakProminences(1:intPeaks);
			%get time vectors
			vecDerivT = vecT(1:(end-1)) + diff(vecT)/2;
			vecBinT = vecDerivT;
		else
			if intResampleNum == vecResamples(end) && vecPlotCells(intNeuron)
				[dblZeta,sOptOut] = getZeta(vecSpikeTimes,matStimOnOff,2,[],intResampleNum,0);
			else
				[dblZeta,sOptOut] = getZeta(vecSpikeTimes,matStimOnOff,intMakePlots,[],intResampleNum,0);
			end
			intSpikeNum = numel(vecSpikeTimes);
			%get smooth deriv
			[vecSmoothDeriv,sDerivative] = getSmoothDeriv(sOptOut.vecSpikeT,intPeaks,intSmoothSd,dblBinSize);
			%unpack
			vecPeakTimes = sDerivative.vecPeakTimes;
			vecPeakValues = sDerivative.vecPeakValues;
			vecBinT = sDerivative.vecBinT;
			vecLDD = sDerivative.vecLDD;
			vecDerivT = sDerivative.vecDerivT;
			vecT = sOptOut.vecSpikeT;
		end

		% assign data
		dblHzD = sOptOut.dblHzD;
		vecNumSpikes(intNeuron) = intSpikeNum;
		vecZeta(intNeuron) = dblZeta;
		vecP(intNeuron) = sOptOut.dblP;
		vecHzD(intNeuron) = dblHzD;
		vecHzP(intNeuron) = sOptOut.dblHzP;
		cellZ{intNeuron} = sOptOut.vecZ;
		cellInterpT{intNeuron} = vecT;
		cellArea{intNeuron} = strArea;
		cellDeriv{intNeuron} = vecSmoothDeriv;
		cellPeakT{intNeuron} = vecPeakTimes;
		cellPeakV{intNeuron} = vecPeakValues;

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
			cla;
			plot(vecBinT,vecSmoothDeriv);
			hold on
			scatter(vecPeakTimes,vecPeakValues,'kx');
			for intPeak=1:intPeaks
				text(vecPeakTimes(intPeak),vecPeakValues(intPeak),sprintf('T=%.0fms',vecPeakTimes(intPeak)/dblBinSize),'FontSize',12);
			end
			hold off
			xlabel('Time from event (s)');
			ylabel('Derivative');
			title(sprintf('%s %s-%s,U%d/C%d',strArea,strDate,strBlock,intSU,intClust));
			fixfig
			drawnow;

			strFileName = sprintf('%s%s-%sSU%dC%d',strArea,strDate,strBlock,intSU,intClust);
			export_fig([strFigPath strFileName 'Zeta' strRunType 'Resamp' num2str(intResampleNum) '.tif']);
			print(gcf, '-dpdf', [strFigPath strFileName 'Zeta' strRunType 'Resamp' num2str(intResampleNum) '.pdf']);

			%% make derivative figure
			if ~exist('vecTraceAct','var') || isempty(vecTraceAct) %calcium imaging

				%make maximized figure
				figure
				drawnow;
				jFig = get(handle(gcf), 'JavaFrame');
				jFig.setMaximized(true);
				figure(gcf);
				drawnow;

				cellDeriv={'Locally Dynamic Derivative','Temporal Offset Derivative','1/(Inter-spike interval)'};
				cellDerivAbr={'LDD','TOD','1/ISI'};
				cellTempDeriv = [];
				cellTempDeriv{1} = vecLDD;
				cellTempDeriv{2} = diff(sOptOut.vecRealDiff)./diff(sOptOut.vecSpikeT);
				cellTempDeriv{3} = 1./diff(sOptOut.vecSpikeT);
				for intDeriv=1:3
					dblMaxT = roundi(max(vecDerivT),3,'floor');
					[vecCounts,vecMeans,vecSDs,cellVals,cellIDs] = makeBins(vecDerivT,cellTempDeriv{intDeriv},vecBinT);
					vecUseVals = find(vecCounts);
					vecInterpVals = interp1([-dblBinSize vecBinT(vecUseVals) dblMaxT+dblBinSize],[0 vecMeans(vecUseVals)' 0],vecBinT);
					vecFilt = normpdf(-2*(intSmoothSd):2*intSmoothSd,0,intSmoothSd);
					vecFilt = vecFilt./sum(vecFilt);
					vecSmoothed = conv(vecInterpVals,vecFilt,'same');
					%find peaks
					[pksSm,locsSm,wSm,pSm] = findpeaks(vecSmoothed,'SortStr','descend','Npeaks',intPeaks);
					[pks,locs,w,p] = findpeaks(vecInterpVals,'SortStr','descend','Npeaks',intPeaks);

					%plot
					subplot(2,3,intDeriv)
					plot(vecBinT,vecInterpVals)
					hold on
					scatter(vecBinT(locs),pks,'kx');
					for intPeak=1:intPeaks
						text(vecBinT(locs(intPeak)),pks(intPeak),sprintf('T=%.0fms',vecBinT(locs(intPeak))/dblBinSize),'FontSize',12);
					end
					hold off
					title([cellDeriv{intDeriv}]);
					xlabel('Time from event (s)');
					ylabel('Derivative');
					fixfig

					subplot(2,3,intDeriv+3)
					plot(vecBinT,vecSmoothed)
					hold on
					scatter(vecBinT(locsSm),pksSm,'kx');
					for intPeak=1:intPeaks
						text(vecBinT(locsSm(intPeak)),pksSm(intPeak),sprintf('T=%.0fms',vecBinT(locsSm(intPeak))/dblBinSize),'FontSize',12);
					end
					hold off
					if intDeriv==1
						title([cellDerivAbr{intDeriv} ', ' getGreek('sigma','lower') sprintf('=%dms; %s %s-%s, SU %d (C%d)',intSmoothSd,strArea,strDate,strBlock,intSU,intClust)]);
					else
						title([cellDerivAbr{intDeriv} sprintf(', %s=%dms;',getGreek('sigma','lower'),intSmoothSd)]);
					end
					xlabel('Time from event (s)');
					ylabel('Derivative');
					fixfig

					drawnow;
					strFileName = sprintf('%s%s-%sSU%dC%d',strArea,strDate,strBlock,intSU,intClust);
					export_fig([strFigPath strFileName 'Deriv' strRunType 'Resamp' num2str(intResampleNum) '.tif']);
					print(gcf, '-dpdf', [strFigPath strFileName 'Deriv' strRunType 'Resamp' num2str(intResampleNum) '.pdf']);

				end
			end
			%strRecIdx = sThisNeuron.RecIdx;
			%strMouse = sThisNeuron.Mouse;
			%strBlock = sThisNeuron.Block;
			%strArea = sThisNeuron.Area;
			%strDate = sThisNeuron.Date;
			%intSU = sThisNeuron.IdxSU;
			%intClust = sThisNeuron.IdxClust;
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