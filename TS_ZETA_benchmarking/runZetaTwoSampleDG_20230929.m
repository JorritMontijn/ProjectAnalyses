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
intRunPairNum = 1000;
boolSave = true;%true;
boolDirectQuantile = false;
global boolWithReplacement;

%% load data
%reset vars
strArea1 = 'Primary visual';
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
strName = replace([lower(strArea1) strRunStim],lower(cellRepStr(:,1)),cellRepStr(:,2));
[sAggStim,sAggNeuron]=loadDataNpx(strArea1,strRunStim,strDataSourcePath);
cellRecIdx = {sAggStim.Exp};
intNeurons = numel(sAggNeuron);

%% pre-allocate output variables
cellNeuron = cell(intRunPairNum,2);
matTtest2 = nan(intRunPairNum,2);
matZeta2 = nan(intRunPairNum,2);
matZeta2_old = nan(intRunPairNum,2);
matAnova2 = nan(intRunPairNum,2);
matAnova2_unbalanced = nan(intRunPairNum,2);
optLow = 2;
optHigh = 1e6;

%% get neuronal data
hTicN = tic;
for intPair = 1:intRunPairNum
	vecPair = randperm(intNeurons,2);
	%% message
	if toc(hTicN) > 5
		fprintf('Processing pair %d/%d [%s]\n',intPair,intRunPairNum,getTime);
		hTicN=tic;
	end
	
	%% get neuron1
	intNeuron1 = vecPair(1);
	sThisNeuron1 = sAggNeuron(intNeuron1);
	vecSpikeTimesN1 = sThisNeuron1.SpikeTimes;
	strRecIdx1 = sThisNeuron1.Exp;
	strMouse1 = sThisNeuron1.Subject;
	strBlock1 = '';
	strArea1 = strName;
	strDate1 = sThisNeuron1.Date;
	intSU1 = sThisNeuron1.Cluster;
	intClust1 = sThisNeuron1.IdxClust;
	
	% get matching recording data
	sThisRec = sAggStim(strcmpi(strRecIdx1,cellRecIdx));
	vecStimOnTime = [];
	vecStimOffTime = [];
	for intRec=1%:numel(sThisRec.cellStim)
		vecStimOnTime = cat(2,vecStimOnTime,sThisRec.cellBlock{intRec}.vecStimOnTime);
		vecStimOffTime = cat(2,vecStimOffTime,sThisRec.cellBlock{intRec}.vecStimOffTime);
	end
	
	matTrialTN1 = [];
	matTrialTN1(:,1) = vecStimOnTime;
	matTrialTN1(:,2) = vecStimOffTime;
	dblUseMaxDur1 = round(median(matTrialTN1(:,2) - matTrialTN1(:,1))*2)/2;
	
	%% neuron 2
	intNeuron2 = vecPair(2);
	sThisNeuron2 = sAggNeuron(intNeuron2);
	vecSpikeTimesN2 = sThisNeuron2.SpikeTimes;
	strRecIdx2 = sThisNeuron2.Exp;
	intSU2 = sThisNeuron2.Cluster;
	
	% get matching recording data
	sThisRec = sAggStim(strcmpi(strRecIdx2,cellRecIdx));
	vecStimOnTime = [];
	vecStimOffTime = [];
	for intRec=1%:numel(sThisRec.cellStim)
		vecStimOnTime = cat(2,vecStimOnTime,sThisRec.cellBlock{intRec}.vecStimOnTime);
		vecStimOffTime = cat(2,vecStimOffTime,sThisRec.cellBlock{intRec}.vecStimOffTime);
	end
	
	matTrialTN2 = [];
	matTrialTN2(:,1) = vecStimOnTime;
	matTrialTN2(:,2) = vecStimOffTime;
	dblUseMaxDur2 = round(median(matTrialTN2(:,2) - matTrialTN2(:,1))*2)/2;
	dblUseMaxDur = min(dblUseMaxDur1,dblUseMaxDur2);
	
	for intRandType=vecRandTypes
		%% randomize
		if intRandType == 2
			%if random, split data of neuron into two times 50% of spikes
			matTrialT1 = matTrialTN1;
			matTrialT2 = matTrialTN1;
			intSpikesN1 = numel(vecSpikeTimesN1);
			intTakeSpikesN1 = round(intSpikesN1/2);
			vecSpikesN1 = randperm(intSpikesN1,intTakeSpikesN1);
			vecSpikesN2 = find(~ismember(1:intSpikesN1,vecSpikesN1));
			vecSpikeTimes1 = sort(vecSpikeTimesN1(vecSpikesN1));
			vecSpikeTimes2 = sort(vecSpikeTimesN1(vecSpikesN2));
		else
			%if not random, take 50% of spikes from neuron 1 and 50% of spikes from neuron 2
			matTrialT1 = matTrialTN1;
			matTrialT2 = matTrialTN2;
			intTakeSpikesN1 = round(numel(vecSpikeTimesN1)/2);
			intTakeSpikesN2 = round(numel(vecSpikeTimesN2)/2);
			vecSpikeTimes1 = sort(vecSpikeTimesN1(randperm(numel(vecSpikeTimesN1),intTakeSpikesN1)));
			vecSpikeTimes2 = sort(vecSpikeTimesN2(randperm(numel(vecSpikeTimesN2),intTakeSpikesN2)));
		end
	
		%% run tests
		intPlot = 0;
		boolWithReplacement = false; %randperm
		dblZeta2P = zetatest2b(vecSpikeTimes1,matTrialT1,vecSpikeTimes2,matTrialT2,dblUseMaxDur,intResampNum,intPlot);
		boolWithReplacement = true; %randi
		dblZeta2P_withrep = zetatest2b(vecSpikeTimes1,matTrialT1,vecSpikeTimes2,matTrialT2,dblUseMaxDur,intResampNum,intPlot);
		
		if 0%sZETA.dblZETA > 2.5 %sZETA.dblMeanZ < 2 && sZETA.dblZETA > 2
			intPlot = 4;
			[dblZeta2P,sZETA] = zetatest2b(vecSpikeTimes1,matTrialT1,vecSpikeTimes2,matTrialT2,dblUseMaxDur,intResampNum,intPlot);
			subplot(2,3,5)
			boxplot([sZETA.vecMu1' sZETA.vecMu2'],'notch','on');
			ylabel('Spiking rate per trial (Hz)');
			set(gca,'xtick',[1 2],'xticklabel',{'Neuron 1','Neuron 2'});
			fixfig;
			return
		end
		
		%% ANOVA
		%if balanced
		hTic2 = tic;
		[vecTrialPerSpike1,vecTimePerSpike1] = getSpikesInTrial(vecSpikeTimes1,matTrialT1(:,1),dblUseMaxDur);
		[vecTrialPerSpike2,vecTimePerSpike2] = getSpikesInTrial(vecSpikeTimes2,matTrialT2(:,1),dblUseMaxDur);
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
		%if balanced
		if intTrials1==intTrials2
			matPSTH = matPSTH1 - matPSTH2;
			dblAnova2P=anova1(matPSTH,[],'off');
		else
			dblAnova2P=1;
		end
		
		%if not balanced
		y = cat(1,matPSTH1(:),matPSTH2(:));
		g1 = cat(1,matLabelN1(:),matLabelN2(:));
		g2 = cat(1,matLabelBin1(:),matLabelBin2(:));
		[vecP,tbl,stats] = anovan(y,{g1 g2},'continuous',[2],'model','interaction','display','off');
		[h crit_p adj_p]=fdr_bh(vecP([1 3]));
		dblAnova2P_unbalanced = min(adj_p);
		
		%% t-test
		%'vecTtestP','vecTtestTime'
		hTic3 = tic;
		vecRespBinsDur1 = sort(flat([matTrialT1(:,1) matTrialT1(:,2)]));
		vecR1 = histcounts(vecSpikeTimes1,vecRespBinsDur1);
		vecD1 = diff(vecRespBinsDur1)';
		vecMu1 = vecR1(1:2:end)./vecD1(1:2:end);
		
		vecRespBinsDur2 = sort(flat([matTrialT2(:,1) matTrialT2(:,2)]));
		vecR2 = histcounts(vecSpikeTimes2,vecRespBinsDur2);
		vecD2 = diff(vecRespBinsDur2)';
		vecMu2 = vecR2(1:2:end)./vecD2(1:2:end);
		
		%get metrics
		[h,dblTtest2P]=ttest2(vecMu1,vecMu2);
		
		%% save
		% assign data
		cellNeuron{intPair,intRandType} = [strArea1 strDate1 'N' num2str(intSU1) 'N' num2str(intSU2)];
		matTtest2(intPair,intRandType) = dblTtest2P;
		matZeta2(intPair,intRandType) = dblZeta2P;
		matZeta2_old(intPair,intRandType) = dblZeta2P_withrep;
		matAnova2(intPair,intRandType) = dblAnova2P;
		matAnova2_unbalanced(intPair,intRandType) = dblAnova2P_unbalanced;
	end
end

%% save
if boolSave
	save([strDataPath 'Zeta2DataAnova' strArea1 'Resamp' num2str(intResampNum) '.mat' ],...
		'cellNeuron','matTtest2','matZeta2','matZeta2_old','matAnova2','matAnova2_unbalanced');
end
