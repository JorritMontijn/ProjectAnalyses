%% load spiking data & plot tuning curves

%% set recording
%close all;
clear all;
if isfolder('F:\Drive\MontijnHeimel_TimeseriesZeta')
	strPath = 'F:\Drive\MontijnHeimel_TimeseriesZeta';
	strDataSourcePath = 'D:\Data\Processed\Steinmetz\';
else
	strPath = 'C:\Drive\MontijnHeimel_TimeseriesZeta';
	strDataSourcePath = 'E:\Steinmetz\';
end
strDataTargetPath = fullfile(strPath,'\Data\');
strFigPath = fullfile(strPath,'\Figs\');
intMakePlots =0; %0=none, 1=normal plot, 2=including raster
vecRandTypes = [1 2];%[1 2];%1=normal,2=rand
vecRestrictRange = [0 inf];
boolSave = true;
intResampNum = 250;%250;%10:10:90;%[10:10:100];
optLow = 2;
optHigh = 100;
boolCombineAll = true;

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

%% load data
sDir = dir(fullpath(strDataSourcePath,'*.mat'));
cellFiles = {sDir.name};

for intRec=1%2:numel(sDir)
	sLoad = load(fullpath(sDir(intRec).folder,sDir(intRec).name));
	sAP = sLoad.sAP;
	strRec = sAP.Rec;
	
	%% prep neural data
	sBehaviour = sAP.sBehaviour;
	sStim = sAP.sStim;
	sCluster = sAP.sCluster;
	sCluster = sCluster(cellfun(@(x) ~isempty(x),{sCluster.KilosortGood}));
	indUseCells = cell2vec({sCluster.KilosortGood})==1;
	cellSpikeT = {sCluster(indUseCells).SpikeTimes};
	cellArea = {sCluster(indUseCells).Area};
	intNeurons = numel(cellArea);
	
	%% prep behavioural
	%remove excluded trials
	sCombos = getStimulusCombos(sStim,{'vecContrastLeft','vecContrastRight'});
	indInclude0 = sStim.indIncluded;
	indCorrectResp = sStim.indCorrect;
	intTrialNum = numel(indInclude0);
	
	%select only unambiguous trials
	indIncludeCombos = (sCombos.matComboVal(:,1) ~= sCombos.matComboVal(:,2));
	indIncludeComboTrials = sum(sCombos.matComboTrials(indIncludeCombos,:));
	
	matComboIdx = sCombos.matComboIdx(indIncludeCombos,:);
	matComboVal = sCombos.matComboVal(indIncludeCombos,:);
	matComboTrials = sCombos.matComboTrials(indIncludeCombos,:);
	
	%use only combos with >2 repetitions per conditions
	intUseCombos = size(matComboTrials,1);
	vecCorrectCounts = nan(1,intUseCombos);
	vecIncorrectCounts = nan(1,intUseCombos);
	matComboTrialsPerResp = false(intUseCombos,intTrialNum,2);
	for intCombo=1:intUseCombos
		indComboT = matComboTrials(intCombo,:);
		%correct
		indCorrect = indComboT(:) & indInclude0 & indCorrectResp;
		vecCorrectCounts(intCombo) = sum(indCorrect);
		matComboTrialsPerResp(intCombo,:,1) = indCorrect;
		
		%incorrect
		indIncorrect = indComboT(:) & indInclude0 & ~indCorrectResp;
		vecIncorrectCounts(intCombo) = sum(indIncorrect);
		matComboTrialsPerResp(intCombo,:,2) = indIncorrect;
	end
	indUseCombos = (vecIncorrectCounts > 2) & (vecCorrectCounts > 2);
	
	matComboIdx = matComboIdx(indUseCombos,:);
	matComboVal = matComboVal(indUseCombos,:);
	matComboTrials = matComboTrials(indUseCombos,:);
	matComboTrialsPerResp = matComboTrialsPerResp(indUseCombos,:,:);
	vecCorrectCounts = vecCorrectCounts(indUseCombos);
	vecIncorrectCounts = vecIncorrectCounts(indUseCombos);
	intUseCombos = sum(indUseCombos);
	indUseTrials = sum(sum(matComboTrialsPerResp,1),3)==1;
	
	dblOffsetT = 0;
	vecStimOnTime = sStim.vecStimOnTime(indUseTrials);
	vecCueDur = sStim.vecGoCueTime(indUseTrials) - sStim.vecStimOnTime(indUseTrials);
	vecStimDur = sStim.vecResponseTime(indUseTrials) - sStim.vecStimOnTime(indUseTrials);
	dblUseMaxDur = min(vecCueDur)-dblOffsetT;
	
	%% run
	if boolCombineAll,intUseCombos=1;end
	matZeta2P = nan(intNeurons,intUseCombos,2);
	matTtest2P = nan(intNeurons,intUseCombos,2);
	matAnova2P = nan(intNeurons,intUseCombos,2);
	matTrialNum = nan(2,intUseCombos,2);
	for intUseCombo=1:intUseCombos
		%% message
		fprintf('Processing %s (%d/%d), combination %d/%d [%s]\n',strRec,intRec,numel(sDir),intUseCombo,intUseCombos,getTime);
		hTic=tic;
		
		%% analyze
		for intNeuron=1:intNeurons%31
			%% message
			if toc(hTic) > 5
				fprintf('Processing neuron %d/%d [%s]\n',intNeuron,intNeurons,getTime);
				hTic=tic;
			end
			clear vecTrialStarts;
			
			for intRand=[1 2]
				%% prep data
				if boolCombineAll
					indCorrectTrials = sum(matComboTrialsPerResp(:,:,1),1)==1;
					indIncorrectTrials = sum(matComboTrialsPerResp(:,:,2),1)==1;
				else
					indCorrectTrials = matComboTrialsPerResp(intUseCombo,:,1);
					indIncorrectTrials = matComboTrialsPerResp(intUseCombo,:,2);
				end
				vecEventOnCorrect = sStim.vecResponseTime(indCorrectTrials) - dblOffsetT;
				matTrialT1 = cat(2,vecEventOnCorrect,vecEventOnCorrect+dblUseMaxDur);
				
				vecEventOnIncorrect = sStim.vecResponseTime(indIncorrectTrials) - dblOffsetT;
				matTrialT2 = cat(2,vecEventOnIncorrect,vecEventOnIncorrect+dblUseMaxDur);
				
				if intRand==2
					matTrialT12 = cat(1,matTrialT1,matTrialT2);
					intT1 = size(matTrialT1,1);
					intT2 = size(matTrialT2,1);
					intTotT = intT1+intT2;
					
					vecUseRand1 = randi(intTotT,[1,intT1]);
					vecUseRand2 = randi(intTotT,[1,intT2]);
					
					vecT1 = sort(matTrialT12(vecUseRand1,1));
					matTrialT1 = cat(2,vecT1,vecT1+dblUseMaxDur);
					
					vecT2 = sort(matTrialT12(vecUseRand2,1));
					matTrialT2 = cat(2,vecT2,vecT2+dblUseMaxDur);
				end
				
				%% get data
				vecSpikeTimes = cellSpikeT{intNeuron};
				
				%% run tests
				intPlot = 0;
				[dblZeta2P,sZETA] = zetatest2(vecSpikeTimes,matTrialT1,vecSpikeTimes,matTrialT2,dblUseMaxDur,intResampNum,intPlot);
				
				%% ANOVA
				[vecTrialPerSpike1,vecTimePerSpike1] = getSpikesInTrial(vecSpikeTimes,matTrialT1(:,1),dblUseMaxDur);
				[vecTrialPerSpike2,vecTimePerSpike2] = getSpikesInTrial(vecSpikeTimes,matTrialT2(:,1),dblUseMaxDur);
				if numel(vecTimePerSpike1) < 3 && numel(vecTimePerSpike2) < 3,continue;end
				vecAllSpikes = cat(1,vecTimePerSpike1,vecTimePerSpike2);
				[optN, dblC, allN, allC] = opthist(vecAllSpikes);
				if optN<optLow,optN=optLow;end %at least 2 bins
				if optN>optHigh,optN=optHigh;end %at least 2 bins
				%optN = dblUseMaxDur/0.01;
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
				dblAnova2P = min(adj_p);
				%dblAnova2P_unbalanced=1;
				
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
				matTrialNum(1,intUseCombo,intRand) = sum(indCorrectTrials);
				matTrialNum(2,intUseCombo,intRand) = sum(indIncorrectTrials);
				matZeta2P(intNeuron,intUseCombo,intRand) = dblZeta2P;
				matTtest2P(intNeuron,intUseCombo,intRand) = dblTtest2P;
				matAnova2P(intNeuron,intUseCombo,intRand) = dblAnova2P;
			end
		end
		%% create summary
		matZeta2P_summary = nan(intNeurons,2);
		matTtest2P_summary = nan(intNeurons,2);
		matAnova2P_summary = nan(intNeurons,2);
		for intRand=[1 2]
			for intNeuron=1:intNeurons
				dblP = bonf_holm(matZeta2P(intNeuron,~isnan(matZeta2P(intNeuron,intUseCombo,intRand)),intRand));
				if isempty(dblP),dblP=1;end
				matZeta2P_summary(intNeuron,intRand) = min(dblP);
				
				dblP = bonf_holm(matTtest2P(intNeuron,~isnan(matTtest2P(intNeuron,intUseCombo,intRand)),intRand));
				if isempty(dblP),dblP=1;end
				matTtest2P_summary(intNeuron,intRand) = min(dblP);
				
				dblP = bonf_holm(matAnova2P(intNeuron,~isnan(matAnova2P(intNeuron,intUseCombo,intRand)),intRand));
				if isempty(dblP),dblP=1;end
				matAnova2P_summary(intNeuron,intRand) = min(dblP);
			end
			if intRand == 1
				strRand = 'Inclusion';
			else
				strRand = 'FPR';
			end
			fprintf('%s %s, combo %d/%d at alpha=0.05, Z=%d%%,T=%d%%,A=%d%%\n',...
				strRec,strRand,intUseCombo,intUseCombos,...
				round(100*(sum(matZeta2P_summary(:,intRand)<0.05)/intNeurons)),...
				round(100*(sum(matTtest2P_summary(:,intRand)<0.05)/intNeurons)),...
				round(100*(sum(matAnova2P_summary(:,intRand)<0.05)/intNeurons)));
		end
		
	end
	%% create summary
	matZeta2P_summary = nan(intNeurons,2);
	matTtest2P_summary = nan(intNeurons,2);
	matAnova2P_summary = nan(intNeurons,2);
	for intRand=[1 2]
		for intNeuron=1:intNeurons
			dblP = bonf_holm(matZeta2P(intNeuron,~isnan(matZeta2P(intNeuron,:,intRand)),intRand));
			%dblP = squeeze(min(matZeta2P(intNeuron,:,intRand),[],2));
			if isempty(dblP),dblP=nan;end
			matZeta2P_summary(intNeuron,intRand) = nanmin(dblP);
			
			dblP = bonf_holm(matTtest2P(intNeuron,~isnan(matTtest2P(intNeuron,:,intRand)),intRand));
			%dblP = squeeze(min(matTtest2P(intNeuron,:,intRand),[],2));
			if isempty(dblP),dblP=nan;end
			matTtest2P_summary(intNeuron,intRand) = nanmin(dblP);
			
			dblP = bonf_holm(matAnova2P(intNeuron,~isnan(matAnova2P(intNeuron,:,intRand)),intRand));
			%dblP = squeeze(min(matAnova2P(intNeuron,:,intRand),[],2));
			if isempty(dblP),dblP=nan;end
			matAnova2P_summary(intNeuron,intRand) = nanmin(dblP);
		end
		if intRand == 1
			strRand = 'Inclusion';
		else
			strRand = 'FPR';
		end
		fprintf('%s %s at alpha=0.05, Z=%d%%,T=%d%%,A=%d%%\n',...
			strRec, strRand,...
			round(100*(sum(matZeta2P_summary(:,intRand)<0.05)/intNeurons)),...
			round(100*(sum(matTtest2P_summary(:,intRand)<0.05)/intNeurons)),...
			round(100*(sum(matAnova2P_summary(:,intRand)<0.05)/intNeurons)));
	end
	%% save data
	if boolSave
		if boolCombineAll
			strCombine = 'Lumped';
		else
			strCombine = 'Split';
		end
		save([strDataTargetPath 'Zeta2Steinmetz' strCombine strRec '.mat' ],...
			'cellArea','matTrialNum','matZeta2P','matTtest2P','matAnova2P');
	end
end