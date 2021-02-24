%% load spiking data & plot tuning curves

%% set recording
close all;
clear all;
strDisk = 'F:';
strDataTargetPath = [strDisk '\Data\Processed\ZETA\Latencies\'];
strFigPath = [strDisk '\Data\Results\ZETA\Latencies\'];
intMakePlots = 2; %0=none, 1=normal plot, 2=including raster
vecRunTypes = 1;%[1 2];
boolSave = true;
vecResamples = 100;%[3 5:5:50];%100];
vecBaseRate = 2.^(-1:0.5:5);
vecJitters = [1:0.5:10];

%set var
for intJitterIdx=1:numel(vecJitters)
	dblJitter = roundi(vecJitters(intJitterIdx),1);
	for intBaseRateIdx=1:numel(vecBaseRate)
	dblBaseRate = roundi(vecBaseRate(intBaseRateIdx),1);
	strArea = sprintf('PoissonPeakRate%02.1fJitter%02.1f',dblBaseRate,dblJitter);
	%reset vars
	clearvars -except dblJitter intJitterIdx vecJitters dblBaseRate strArea intBaseRateIdx vecBaseRate vecResamples boolSave strDataTargetPath strFigPath intMakePlots vecRunTypes
	
for intRunType=vecRunTypes
if intRunType == 1
	strRunType = strArea;
elseif intRunType ==2
	strRunType = [strArea '-Rand'];
end
%% load data
intNeurons = 10;
for intResampleIdx = 1:numel(vecResamples)
	intResampleNum = vecResamples(intResampleIdx);
	%% message
	fprintf('Processing %s, resampling %d (%d/%d) [%s]\n',strRunType,intResampleNum,intResampleIdx,numel(vecResamples),getTime);
	hTic=tic;

	%% pre-allocate output variables
	vecNumSpikes = nan(1,intNeurons);
	vecZetaP = nan(1,intNeurons);
	vecP = nan(1,intNeurons);
	vecHzD = nan(1,intNeurons);
	vecHzP = nan(1,intNeurons);
	cellZ = cell(1,intNeurons);
	cellArea = cell(1,intNeurons);
	cellDeriv = cell(1,intNeurons);
	cellInterpT = cell(1,intNeurons);
	cellPeakT = cell(1,intNeurons);
	cellPeakI = cell(1,intNeurons);
	
	vecMIMIP = nan(1,intNeurons);
	vecComputTimeZETA = nan(1,intNeurons);
	vecComputTimeMIMI = nan(1,intNeurons);
	cellCoeffs = cell(1,intNeurons);
	vecLatenciesZ = nan(1,intNeurons);
	vecLatenciesM = nan(1,intNeurons);
	
	%% analyze
	for intNeuron=[1:intNeurons]%31
		%% message
		if toc(hTic) > 5
			fprintf('Processing neuron %d/%d [%s]\n',intNeuron,intNeurons,getTime);
			hTic=tic;
		end
		clear vecTrialStarts;
		
		%% generate data
		strRecIdx = 'x';
		strMouse = 'Artificial';
		strBlock = '1';
		strDate = getDate();
		intSU = intNeuron;
		intClust = intNeuron;
		
		%set parameters
		dblPrefRate = dblBaseRate;
		dblKappa = rand(1)*5+5;
		vecTrialAngles=deg2rad(repmat([0:45:359],[1 20]));
		dblTrialDur=2;
		vecStimOnTime = dblTrialDur*(1:numel(vecTrialAngles))';
		vecStimOffTime = vecStimOnTime + 1;
		
		vecTrialStarts(:,1) = vecStimOnTime;
		vecTrialStarts(:,2) = vecStimOffTime;
		[vecSpikeTimes,dblPrefOri] = getGeneratedSpikingDataWithPeak(vecTrialAngles,vecTrialStarts,dblBaseRate,dblPrefRate,dblJitter,dblKappa,true);
		
		%% get visual responsiveness
		%set derivative params
		intPeaks = 3;
		intSmoothSd = 3;
		dblBinSize = (1/1000);
		dblDur = roundi(median(diff(vecStimOnTime)),1);
		if contains(strRunType,'Rand')
			vecJitter = 2*dblDur*rand(size(vecStimOnTime))-dblDur;
			matStimOnOff = bsxfun(@plus,vecTrialStarts,vecJitter);
		else
			matStimOnOff = vecTrialStarts;
		end
		close;close;
		
		hTic = tic;
		if intResampleNum == vecResamples(end) && intNeuron == 1
			[dblZetaP,vecLatZ,sZETA,sMSD] = getZeta(vecSpikeTimes,matStimOnOff,dblDur,intResampleNum,intMakePlots,4);
		else
			[dblZetaP,vecLatZ,sZETA,sMSD] = getZeta(vecSpikeTimes,matStimOnOff,dblDur,intResampleNum,0,4);
		end
		dblComputTimeZ = toc(hTic);
		%MIMI
		hTic = tic;
		[dblMIMI_P,vecLatM,sMIMI,sRate] = getMIMI(vecSpikeTimes,matStimOnOff,dblDur);
		dblComputTimeM = toc(hTic);
		
		
		intSpikeNum = numel(vecSpikeTimes);
		%unpack
		vecT = sZETA.vecSpikeT;
		if isempty(sMSD)
			vecRate = [];
			vecTime = [];
			dblPeakWidth = [];
		else
			vecRate = sMSD.vecRate;
			vecTime = sMSD.vecT;
			dblPeakWidth = sMSD.dblPeakWidth;
		end

		% assign data
		vecNumSpikes(intNeuron) = intSpikeNum;
		vecZetaP(intNeuron) = dblZetaP;
		vecP(intNeuron) = sZETA.dblP;
		vecHzD(intNeuron) = sZETA.dblMeanD;
		vecHzP(intNeuron) = sZETA.dblMeanP;
		cellZ{intNeuron} = sZETA.vecD;
		cellInterpT{intNeuron} = vecT;
		cellArea{intNeuron} = strArea;
		cellDeriv{intNeuron} = vecRate;
		cellPeakT{intNeuron} = vecTime;%vecPeakTimes;
		%cellPeakI{intNeuron} = vecIdx;%vecPeakValues;
		
		vecMIMIP(intNeuron) = dblMIMI_P;
		vecComputTimeZETA(intNeuron) = dblComputTimeZ;
		vecComputTimeMIMI(intNeuron) = dblComputTimeM;
		cellCoeffs{intNeuron} = sMIMI.FitCoeffs;
		vecX = sMIMI.vecX;
		vecLatenciesZ(intNeuron) = vecLatZ(3);
		vecLatenciesM(intNeuron) = vecLatM(1);
		
		%continue;
		%% build vector for cells to plot
		%% save plot
		if intMakePlots && intNeuron == 1 &&(intResampleNum == vecResamples(end))% && ~(exist('vecTraceAct','var') && ~isempty(vecTraceAct)))
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
			close;
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
		save([strDataTargetPath 'ZetaDataPPL3' strRunType 'Resamp' num2str(intResampleNum) '.mat' ],...
			'dblBaseRate','vecNumSpikes','vecZetaP','vecP','vecHzD','vecHzP','cellZ','cellInterpT','cellArea','cellDeriv','cellPeakT','cellPeakI',...
			'vecMIMIP','vecComputTimeZETA','vecComputTimeMIMI','cellCoeffs','vecX','vecLatenciesZ','vecLatenciesM');
end
end
end
	end
end