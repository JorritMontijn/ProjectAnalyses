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
matTtest2 = ones(intNeurons,2);
matZeta2 = ones(intNeurons,2);
matAnova2 = ones(intNeurons,2);
matAnova2_unbalanced = ones(intNeurons,2);
matAnova2_optimal = ones(intNeurons,2);

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
	dblDur = min(dblDur2,dblDur4);
	dblUseMaxDur = 1.5;%dblDur-dblShiftBy;
	
	%stim 1
	matTrialTS1 = [];
	matTrialTS1(:,1) = vecTransition2+dblShiftBy;
	matTrialTS1(:,2) = vecTransition2+dblUseMaxDur+dblShiftBy;
	intTrialsS1 = numel(vecTransition2);
	
	%stim 2
	matTrialTS2 = [];
	matTrialTS2(:,1) = vecTransition4+dblShiftBy;
	matTrialTS2(:,2) = vecTransition4+dblUseMaxDur+dblShiftBy;
	intTrialsS2 = numel(vecTransition4);
	
	%% get visual responsiveness
	for intRandType=vecRandTypes
		%% randomize
		if intRandType == 2
			%if random, use 50% of trials from s1 and 50 from s2 for both sets
			vecUseS1for1 = randperm(intTrialsS1,round(intTrialsS1/2));
			vecUseS2for1 = randperm(intTrialsS2,round(intTrialsS2/2));
			vecUseS1for2 = find(~ismember(1:intTrialsS1,vecUseS1for1));
			vecUseS2for2 = find(~ismember(1:intTrialsS2,vecUseS2for1));
			
			matTrialT1 = cat(1,matTrialTS1(vecUseS1for1,:),matTrialTS2(vecUseS2for1,:));
			matTrialT2 = cat(1,matTrialTS1(vecUseS1for2,:),matTrialTS2(vecUseS2for2,:));
		else
			%if not random, use normal data
			matTrialT1 = matTrialTS1;
			matTrialT2 = matTrialTS2;
		end
		
		%check for sufficient spikes
		intSpikeNum = numel(vecSpikeTimes);
		if intSpikeNum<10,continue;end
		
		dblBinW = dblFrameDur;
		[vecMean,vecSEM,vecWindowBinCenters,matPET] = ...
			doPEP(vecSpikeTimes,[0:dblBinW:dblUseMaxDur],matTrialT1(:,1));
		intSpikesInWindow = sum(matPET(:)*dblBinW);
		
		if intSpikesInWindow < 10,continue;end
		
		%% run tests
		intPlot = 0;
		dblZeta2P = zetatest2(vecSpikeTimes,matTrialT1,vecSpikeTimes,matTrialT2,dblUseMaxDur,intResampNum,intPlot,boolDirectQuantile);
		
		%% ANOVA
		%if balanced
		hTic2 = tic;
		[vecTrialPerSpike1,vecTimePerSpike1] = getSpikesInTrial(vecSpikeTimes,matTrialT1(:,1),dblUseMaxDur);
		[vecTrialPerSpike2,vecTimePerSpike2] = getSpikesInTrial(vecSpikeTimes,matTrialT2(:,1),dblUseMaxDur);
		if numel(vecTimePerSpike1) < 3 && numel(vecTimePerSpike2) < 3,continue;end
		xComb = sort(cat(1,vecTimePerSpike1,vecTimePerSpike2));
		[optN, dblC, allN, allC] = opthist(xComb);
		if optN<optLow,optN=optLow;end %at least 2 bins
		if optN>optHigh,optN=optHigh;end %at least 2 bins
		
		intTrials1 = size(matTrialT1,1);
		intTrials2 = size(matTrialT2,1);
		dblBinWidth = dblUseMaxDur/optN;
		vecBins = 0:dblBinWidth:dblUseMaxDur;
		matPSTH1 = nan(intTrials1,optN);
		matLabelN1 = ones(size(matPSTH1));
		matLabelBin1 = repmat(1:optN,[intTrials1 1]);
		for intTrial=1:intTrials1
			matPSTH1(intTrial,:) = histcounts(vecTimePerSpike1(vecTrialPerSpike1==intTrial),vecBins);
		end
		matPSTH2 = nan(intTrials2,numel(vecBins)-1);
		matLabelN2 = 2*ones(size(matPSTH2));
		matLabelBin2 = repmat(1:optN,[intTrials2 1]);
		for intTrial=1:intTrials2
			matPSTH2(intTrial,:) = histcounts(vecTimePerSpike2(vecTrialPerSpike2==intTrial),vecBins);
		end
		
		%if not balanced
		y = cat(1,matPSTH1(:),matPSTH2(:));
		g1 = cat(1,matLabelN1(:),matLabelN2(:));
		g2 = cat(1,matLabelBin1(:),matLabelBin2(:));
		[vecP,tbl,stats] = anovan(y,{g1 g2},'model','interaction','display','off');
		[h crit_p adj_p]=fdr_bh(vecP([1 3]));
		dblAnova2P_optimal = min(adj_p);
		
		%fixed bin
		vecBins = 0:dblFixedBinWidth:dblUseMaxDur;
		optN = numel(vecBins)-1;
		matPSTH1 = nan(intTrials1,optN);
		matLabelN1 = ones(size(matPSTH1));
		matLabelBin1 = repmat(1:optN,[intTrials1 1]);
		for intTrial=1:intTrials1
			matPSTH1(intTrial,:) = histcounts(vecTimePerSpike1(vecTrialPerSpike1==intTrial),vecBins);
		end
		matPSTH2 = nan(intTrials2,numel(vecBins)-1);
		matLabelN2 = 2*ones(size(matPSTH2));
		matLabelBin2 = repmat(1:optN,[intTrials2 1]);
		for intTrial=1:intTrials2
			matPSTH2(intTrial,:) = histcounts(vecTimePerSpike2(vecTrialPerSpike2==intTrial),vecBins);
		end
		%if balanced
		if intTrials1==intTrials2
			matPSTH = matPSTH1 - matPSTH2;
			dblAnova2P=anova1(matPSTH,[],'off');
		else
			dblAnova2P=1;
		end
		
		y = cat(1,matPSTH1(:),matPSTH2(:));
		g1 = cat(1,matLabelN1(:),matLabelN2(:));
		g2 = cat(1,matLabelBin1(:),matLabelBin2(:));
		[vecP,tbl,stats] = anovan(y,{g1 g2},'model','interaction','display','off');
		[h crit_p adj_p]=fdr_bh(vecP([1 3]));
		dblAnova2P_unbalanced= min(adj_p);
		
		%% t-test
		%'vecTtestP','vecTtestTime'
		hTic3 = tic;
		vecRespBinsDur1 = sort(flat([matTrialT1(:,1) matTrialT1(:,2)]));
		vecR1 = histcounts(vecSpikeTimes,vecRespBinsDur1);
		vecD1 = diff(vecRespBinsDur1)';
		vecMu1 = vecR1(1:2:end)./vecD1(1:2:end);
		
		vecRespBinsDur2 = sort(flat([matTrialT2(:,1) matTrialT2(:,2)]));
		vecR2 = histcounts(vecSpikeTimes,vecRespBinsDur2);
		vecD2 = diff(vecRespBinsDur2)';
		vecMu2 = vecR2(1:2:end)./vecD2(1:2:end);
		
		%get metrics
		[h,dblTtest2P]=ttest2(vecMu1,vecMu2);
		
		%% save
		% assign data
		cellNeuron{intIdx,intRandType} = [strArea strDate 'N' num2str(intSU)];
		matTtest2(intIdx,intRandType) = dblTtest2P;
		matZeta2(intIdx,intRandType) = dblZeta2P;
		matAnova2(intIdx,intRandType) = dblAnova2P;
		matAnova2_unbalanced(intIdx,intRandType) = dblAnova2P_unbalanced;
		matAnova2_optimal(intIdx,intRandType) = dblAnova2P_optimal;
	end
end

%% save
if boolSave
	save([strDataTargetPath 'Zeta2Data' strArea 'Resamp' num2str(intResampNum) strQ '.mat' ],...
		'cellNeuron','matTtest2','matZeta2','matAnova2','matAnova2_optimal','matAnova2_unbalanced','dblShiftBy','dblUseMaxDur');
end
