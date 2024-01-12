%% load spiking data & plot tuning curves

%% set recording
%close all;
clear all;

if isfolder('F:\Drive\MontijnHeimel_TimeseriesZeta')
	strPath = 'F:\Drive\MontijnHeimel_TimeseriesZeta';
	strDataSourcePath = 'F:\Data\Processed\Neuropixels\';
else
	strPath = 'C:\Drive\MontijnHeimel_TimeseriesZeta';
	strDataSourcePath = 'E:\DataPreProcessed\';
end
strDataTargetPath = fullfile(strPath,'\Data\');
strFigPath = fullfile(strPath,'\Figs\');
intMakePlots =0; %0=none, 1=normal plot, 2=including raster
vecRandTypes = [1 2];%[1 2];%1=normal,2=rand
vecRestrictRange = [0 inf];
boolSave = true;
boolDirectQuantile = false;
intResampNum = 500;%250;%10:10:90;%[10:10:100];
optLow = 2;
optHigh = 1000;
strRunStim = 'RunNaturalMovie';
strArea = 'primary visual';

%% rep
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
vecTrialNum = [];
%set var
if boolDirectQuantile
	strQ = 'Q1';
else
	strQ = 'Q0';
end

%% load data
strName = replace([lower(strArea) strRunStim],lower(cellRepStr(:,1)),cellRepStr(:,2));
[sAggStim,sAggNeuron]=loadDataNpx(strArea,strRunStim,strDataSourcePath);
if isempty(sAggStim),return;end
cellRecIdx = {sAggStim.Exp};
intNeurons = numel(sAggNeuron);

%% pre-allocate output variables
cellNeuron = cell(intNeurons,2);
matNumSpikes = zeros(intNeurons,2);
matZetaP_Old = ones(intNeurons,2);
matZetaP_Stitch = ones(intNeurons,2);
matZetaP_NoStitch = ones(intNeurons,2);
matAnovaP = ones(intNeurons,2);
matAnovaP_optimal = ones(intNeurons,2);
matTtestP = ones(intNeurons,2);
%load([strDataTargetPath 'ZetaDataAnova' strRunType strRunStim '.mat']);

%% analyze
hTic = tic;
dblFixedBinWidth = 50/1000;
for intIdx=1:intNeurons%26
	%% message
	if toc(hTic) > 5 || intIdx==1
		fprintf('Processing neuron %d/%d [%s]\n',intIdx,intNeurons,getTime);
		hTic=tic;
	end
	clear vecTrialStarts;
	
	%% get neuronal data
	sThisNeuron = sAggNeuron(intIdx);
	vecSpikeTimes = sThisNeuron.SpikeTimes;
	strRecIdx = sThisNeuron.Exp;
	strMouse = sThisNeuron.Subject;
	strBlock = '';
	strArea = strName;
	strDate = sThisNeuron.Date;
	intSU = sThisNeuron.Cluster;
	intClust = sThisNeuron.IdxClust;
	
	%% find transitions
	boolFindTransitionsInMovie = false;
	if boolFindTransitionsInMovie
		vecDiff = nan(1,500);
		for intFrame1=1:500
			intFrame2 = modx(intFrame1+1,500);
			vecDiff(intFrame2) = sum(flat(Earthflight_WingedPlanet__CondorFlightSchool_NarratedByDavidTen(intFrame1).cdata...
				- Earthflight_WingedPlanet__CondorFlightSchool_NarratedByDavidTen(intFrame2).cdata));
		end
		vecFrameTransitions = find(vecDiff>6e7); %1    91   251   340
	else
		vecFrameTransitions = [1 91 251 340];
	end
	
	%% get matching recording data
	sThisRec = sAggStim(strcmpi(strRecIdx,cellRecIdx));
	vecStimOnTime = [];
	vecStimOffTime = [];
	for intRec=1:numel(sThisRec.cellBlock)
		vecStimOnTime = cat(2,vecStimOnTime,sThisRec.cellBlock{intRec}.vecStimOnTime);
		vecStimOffTime = cat(2,vecStimOffTime,sThisRec.cellBlock{intRec}.vecStimOffTime);
	end
	vecStimDur = vecStimOffTime-vecStimOnTime;
	vecTransition1 = vecStimOnTime+((vecFrameTransitions(1)-1)/500)*vecStimDur;
	vecTransition2 = vecStimOnTime+((vecFrameTransitions(2)-1)/500)*vecStimDur;
	vecTransition3 = vecStimOnTime+((vecFrameTransitions(3)-1)/500)*vecStimDur;
	vecTransition4 = vecStimOnTime+((vecFrameTransitions(4)-1)/500)*vecStimDur;
	dblTotDur = 500/60;
	dblFrameDur = dblTotDur/500;
	dblDur1 = median(vecTransition2-vecTransition1);
	dblDur2 = median(vecTransition3-vecTransition2);
	dblDur3 = median(vecTransition4-vecTransition3);
	dblDur4 = median(vecStimOffTime-vecTransition4);
	dblShiftBy = -0.5;
	dblUseMaxDur = 3;
	
	%% get visual responsiveness
	%get trial dur
	%set derivative params
	%get trial dur
	%set derivative params
	
	for intRandType=vecRandTypes
		matEventTimes = cat(2,vecTransition3(:),vecStimOffTime(:));
		if intRandType == 2
			dblDur = dblUseMaxDur;
			vecJitter = (2*dblDur*rand([numel(vecStimOnTime) 1])-dblDur);
			matEventTimes = bsxfun(@plus,matEventTimes,vecJitter);
		else
			matEventTimes = matEventTimes;
		end
		vecTrialNum = unique([vecTrialNum size(matEventTimes,1)]);
		intTrials = size(matEventTimes,1);
		
		%check for sufficient spikes
		intSpikeNum = numel(vecSpikeTimes);
		if intSpikeNum<10,continue;end
		
		dblBinW = dblFrameDur;
		[vecMean,vecSEM,vecWindowBinCenters,matPET] = ...
			doPEP(vecSpikeTimes,[dblShiftBy:dblBinW:(dblShiftBy+3*dblUseMaxDur)],matEventTimes(:,1));
		intSpikesInWindow = sum(matPET(:)*dblBinW);
		
		if intSpikesInWindow < 10,continue;end
		
		%% run zeta
		intGetLatencies = 0;
		intPlot = 0;
		dblZetaP_old = getZeta(vecSpikeTimes,matEventTimes+dblShiftBy,dblUseMaxDur,intResampNum,intPlot);
		dblZetaP_Stitch=zetatest(vecSpikeTimes,matEventTimes+dblShiftBy,dblUseMaxDur,intResampNum,intPlot,[],[],[],true);
		dblZetaP_NoStitch=zetatest(vecSpikeTimes,matEventTimes+dblShiftBy,dblUseMaxDur,intResampNum,intPlot,[],[],[],false);
		if 0
			if dblZetaP_NoStitch < 0.01
				doPEP(vecSpikeTimes,-0.5:0.1:dblUseMaxDur,matEventTimes(:,1));
				pause;
			end
			continue;
		end
		
		%% ANOVA optimal
		hTic2 = tic;
		[vecTrialPerSpike,vecTimePerSpike] = getSpikesInTrial(vecSpikeTimes,matEventTimes(:,1),dblUseMaxDur);
		if numel(vecTimePerSpike) < 3,continue;end
		[optN, dblC, allN, allC] = opthist(vecTimePerSpike,[],optHigh);
		if optN<optLow,optN=optLow;end %at least 2 bins
		if optN>optHigh,optN=optHigh;end %at least 2 bins
		
		dblBinWidth = dblUseMaxDur/optN;
		vecBins = 0:dblBinWidth:dblUseMaxDur;
		matPSTH = nan(intTrials,numel(vecBins)-1);
		for intTrial=1:intTrials
			matPSTH(intTrial,:) = histcounts(vecTimePerSpike(vecTrialPerSpike==intTrial),vecBins);
		end
		dblAnovaP_optimal=anova1(matPSTH,[],'off');
		dblAnovaDur = toc(hTic2);
		
		%% ANOVA fixed
		vecBins = 0:dblFixedBinWidth:dblUseMaxDur;
		matPSTH = nan(intTrials,numel(vecBins)-1);
		for intTrial=1:intTrials
			matPSTH(intTrial,:) = histcounts(vecTimePerSpike(vecTrialPerSpike==intTrial),vecBins);
		end
		dblAnovaP=anova1(matPSTH,[],'off');
		
		%% t-test
		%'vecTtestP','vecTtestTime'
		hTic3 = tic;
		vecRespBinsDur = sort(flat([matEventTimes(:,1) matEventTimes(:,2)]));
		vecR = histcounts(vecSpikeTimes,vecRespBinsDur);
		vecD = diff(vecRespBinsDur)';
		vecMu_Dur = vecR(1:2:end)./vecD(1:2:end);
		dblStart1 = min(vecRespBinsDur);
		dblFirstPreDur = dblStart1 - max([0 dblStart1 - median(vecD(2:2:end))]);
		dblR1 = sum(vecSpikeTimes > (dblStart1 - dblFirstPreDur) & vecSpikeTimes < dblStart1);
		vecMu_Pre = [dblR1 vecR(2:2:end)]./[dblFirstPreDur vecD(2:2:end)];
		
		%get metrics
		dblMeanD = mean(vecMu_Dur - vecMu_Pre) / ( (std(vecMu_Dur) + std(vecMu_Pre))/2);
		[h,dblTtestP]=ttest(vecMu_Dur,vecMu_Pre);
		dblTtestDur = toc(hTic3);
		
		%%
		% assign data
		cellNeuron{intIdx,intRandType} = [strArea strDate 'N' num2str(intSU)];
		matNumSpikes(intIdx,intRandType) = intSpikeNum;
		matZetaP_Old(intIdx,intRandType) = dblZetaP_old;
		matZetaP_Stitch(intIdx,intRandType) = dblZetaP_Stitch;
		matZetaP_NoStitch(intIdx,intRandType) = dblZetaP_NoStitch;
		matAnovaP(intIdx,intRandType) = dblAnovaP;
		matAnovaP_optimal(intIdx,intRandType) = dblAnovaP_optimal;
		matTtestP(intIdx,intRandType) = dblTtestP;
	end
end
if boolSave
	save([strDataTargetPath 'ZetaData' strArea 'Resamp' num2str(intResampNum)  '.mat' ],...
		'cellNeuron','matNumSpikes','dblShiftBy','dblUseMaxDur','matZetaP_Old',...
		'matZetaP_Stitch','matZetaP_NoStitch',...
		'matAnovaP','matAnovaP_optimal','matTtestP');
end
