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
strDataPath = fullfile(strPath,'\Data\');
strFigPath = fullfile(strPath,'\Figs\');

vecRandTypes = [1 2];
intResampNum = 500;
intRunNum = inf;
boolSave = true;%true;
boolDirectQuantile = false;

%% load data
%reset vars
strArea = 'Primary visual';
strRunStim = 'RunDriftingGratings';%'RunNaturalMovie'
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

%% load data
strName = replace([lower(strArea) strRunStim],lower(cellRepStr(:,1)),cellRepStr(:,2));
[sAggStim,sAggNeuron]=loadDataNpx(strArea,strRunStim,strDataSourcePath);
cellRecIdx = {sAggStim.Exp};
intNeurons = numel(sAggNeuron);

%% pre-allocate output variables
intRunNum = min(intNeurons,intRunNum);
cellNeuron = cell(intRunNum,2);
matTtest2 = nan(intRunNum,2);
matZeta2 = nan(intRunNum,2);
matZeta2_old = nan(intRunNum,2);
matAnova2 = nan(intRunNum,2);
matAnova2_unbalanced = nan(intRunNum,2);
matAnova2_optimal = nan(intRunNum,2);
optLow = 2;
optHigh = 1e3;

%% get neuronal data
hTicN = tic;
vecRunNeurons = sort(randperm(intNeurons,intRunNum));
for intIdx = 430:intRunNum %285
	%% message
	if toc(hTicN) > 5
		fprintf('Processing neuron %d/%d [%s]\n',intIdx,intRunNum,getTime);
		hTicN=tic;
	end
	
	%% get neuron
	intNeuron = vecRunNeurons(intIdx);
	sThisNeuron = sAggNeuron(intNeuron);
	vecSpikeTimes = sThisNeuron.SpikeTimes;
	strRecIdx = sThisNeuron.Exp;
	strMouse = sThisNeuron.Subject;
	strBlock = '';
	strArea = strName;
	strDate = sThisNeuron.Date;
	intSU = sThisNeuron.Cluster;
	intClust = sThisNeuron.IdxClust;
	
	% get matching recording data
	sThisRec = sAggStim(strcmpi(strRecIdx,cellRecIdx));
	vecStimOnTime = [];
	vecStimOffTime = [];
	vecStimTypes = [];
	for intRec=1:numel(sThisRec.cellBlock)
		vecStimOnTime = cat(2,vecStimOnTime,sThisRec.cellBlock{intRec}.vecStimOnTime);
		vecStimOffTime = cat(2,vecStimOffTime,sThisRec.cellBlock{intRec}.vecStimOffTime);
		vecStimTypes = cat(2,vecStimTypes,sThisRec.cellBlock{intRec}.vecTrialStimTypes);
	end
	intStimNum = numel(unique(vecStimTypes));
	
	% 		%check shift
	% 		matAllEvents = cat(2,vecStimOnTime',vecStimOffTime');
	% 		matT1 = matAllEvents(1:480,:);
	% 		matT2 = matAllEvents(481:960,:);
	% 		dblShift = -0.25;
	%
	% 		dblZeta2P1 = zetatest2(vecSpikeTimes,matT1,vecSpikeTimes,matT2,dblUseMaxDur,intResampNum,4,boolDirectQuantile);
	% 		dblZeta2P2 = zetatest2(vecSpikeTimes,matT1,vecSpikeTimes,matT2,dblUseMaxDur,intResampNum,intPlot,boolDirectQuantile);
	% 		dblZeta2P3 = zetatest2(vecSpikeTimes,matT1,vecSpikeTimes,matT2,dblUseMaxDur,intResampNum,intPlot,boolDirectQuantile);
	%
	% 		dblZeta2P1b = zetatest2(vecSpikeTimes,matT1+dblShift,vecSpikeTimes,matT2+dblShift,dblUseMaxDur,intResampNum,4,boolDirectQuantile);
	% 		dblZeta2P2b = zetatest2(vecSpikeTimes,matT1+dblShift,vecSpikeTimes,matT2+dblShift,dblUseMaxDur,intResampNum,intPlot,boolDirectQuantile);
	% 		dblZeta2P3b = zetatest2(vecSpikeTimes,matT1+dblShift,vecSpikeTimes,matT2+dblShift,dblUseMaxDur,intResampNum,intPlot,boolDirectQuantile);
	%
	% 		dblZeta2P1c = zetatest2(vecSpikeTimes,matT1-dblShift,vecSpikeTimes,matT2-dblShift,dblUseMaxDur,intResampNum,4,boolDirectQuantile);
	% 		dblZeta2P2c = zetatest2(vecSpikeTimes,matT1-dblShift,vecSpikeTimes,matT2-dblShift,dblUseMaxDur,intResampNum,intPlot,boolDirectQuantile);
	% 		dblZeta2P3c = zetatest2(vecSpikeTimes,matT1-dblShift,vecSpikeTimes,matT2-dblShift,dblUseMaxDur,intResampNum,intPlot,boolDirectQuantile);
	%
	%check which stim to use
	%vecDur = vecStimOffTime-vecStimOnTime;
	%vecSpikeCounts = getSpikeCounts(vecSpikeTimes,vecStimOnTime,vecStimOffTime)./vecDur;
	%[matRespNSR,vecStimTypes,vecUniqueDegs] = getStimulusResponses(vecSpikeCounts,vecStimTypes);
	%vecMuPerS = mean(matRespNSR,3);
	%[dummy,intStim1] = max(vecMuPerS);
	intStim1 = 1;
	intStim2 = modx(round(intStim1 + intStimNum/4),intStimNum);
	
	%stim 1
	indUseTrials1 = vecStimTypes==intStim1;
	matTrialTS1 = [];
	matTrialTS1(:,1) = vecStimOnTime(indUseTrials1);
	matTrialTS1(:,2) = vecStimOffTime(indUseTrials1);
	dblUseMaxDur1 = round(median(diff(vecStimOnTime))*2)/2;
	intTrialsS1 = sum(indUseTrials1);
	
	%stim 2
	indUseTrials2 = vecStimTypes==intStim2;
	matTrialTS2 = [];
	matTrialTS2(:,1) = vecStimOnTime(indUseTrials2);
	matTrialTS2(:,2) = vecStimOffTime(indUseTrials2);
	dblUseMaxDur2 = round(median(diff(vecStimOnTime))*2)/2;
	dblUseMaxDur = min(dblUseMaxDur1,dblUseMaxDur2);
	intTrialsS2 = sum(indUseTrials1);
	dblFixedBinWidth = 50/1000;
	
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
		
		%% zeta
		intPlot = 0;
		[dblZeta2P,sZETA2] = zetatest2(vecSpikeTimes,matTrialT1,vecSpikeTimes,matTrialT2,dblUseMaxDur,intResampNum,intPlot,boolDirectQuantile);
		
		%% plot?
		if 1
			if dblZeta2P < 0.01 && sZETA2.dblMeanP > 0.05 %190? 267? 279?
				dblZeta2P = zetatest2(vecSpikeTimes,matTrialT1,vecSpikeTimes,matTrialT2,dblUseMaxDur,intResampNum,4,boolDirectQuantile);
				
				
				
				subplot(2,3,5);
				boxplot(cat(2,sZETA2.vecMu1,sZETA2.vecMu2),cat(2,0*sZETA2.vecMu1,0*sZETA2.vecMu2+1));
				
				
				subplot(2,3,6);
				cla;bplot(sZETA2.vecMu1,1,'outliers');
				hold on;
				bplot(sZETA2.vecMu2,2,'outliers');
				pause;close
			end
			continue;
		end
		
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
		cellNeuron{intIdx,intRandType} = [strArea strDate 'N' num2str(intSU) 'S' num2str(intStim1)];
		matTtest2(intIdx,intRandType) = dblTtest2P;
		matZeta2(intIdx,intRandType) = dblZeta2P;
		matAnova2(intIdx,intRandType) = dblAnova2P;
		matAnova2_unbalanced(intIdx,intRandType) = dblAnova2P_unbalanced;
		matAnova2_optimal(intIdx,intRandType) = dblAnova2P_optimal;
	end
end

%% save
if boolSave
	save([strDataPath 'Zeta2DataStimDiff' strArea 'Resamp' num2str(intResampNum) '.mat' ],...
		'cellNeuron','matTtest2','matZeta2','matAnova2','matAnova2_optimal','matAnova2_unbalanced');
end
