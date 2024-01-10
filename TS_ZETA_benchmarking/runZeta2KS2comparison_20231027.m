%% load spiking data & plot tuning curves

%% set recording
%close all;
clear all;
cellUniqueAreas = {...
	'HeteroPoissonPeak',...Area 1
	'TriPhasic',...Area 2
	'QuadriPhasic',...Area 3
	'iidGaussian',...Area 4
	'TriHomoPhasic',...Area 5
	'',...Area 6
	'lateral geniculate',...Area 7
	'Primary visual',...Area 8
	};

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
boolSave = true;
intResampleNum = 1000;%250;%10:10:90;%[10:10:100];
intNeurons = 10000;%10000
vecRunTrialNums = [2 10:10:100 200:100:1000];
intRunTrialNums = numel(vecRunTrialNums);
			
cellNeuron = cell(intRunTrialNums,intNeurons,2);
matNumSpikes = nan(intRunTrialNums,intNeurons,2);
matZetaP = nan(intRunTrialNums,intNeurons,2);
matTtestP = nan(intRunTrialNums,intNeurons,2);
matKsP = nan(intRunTrialNums,intNeurons,2);
%set var
for intTrialNumIdx=2%19:intRunTrialNums
	intTrialNum = vecRunTrialNums(intTrialNumIdx);
	for intRandType=vecRandTypes
		%reset vars
		
		if intRandType == 1
			strRunType = 'KS2Z';
			fprintf('Prepping normal... [%s]\n',getTime);
		elseif intRandType ==2
			strRunType = 'KS2Z-Rand';
			fprintf('Prepping random... [%s]\n',getTime);
		end
		
		%% message
		fprintf('Processing %s, trial # %d (%d/%d) [%s]\n',strRunType,intTrialNum,intTrialNumIdx,intRunTrialNums,getTime);
		hTic=tic;
		
		%% analyze
		dblStimDur = 1;%5
		dblBinW = 0.5;
		for intNeuron=1:intNeurons%26
			%% message
			if toc(hTic) > 5
				fprintf('Processing neuron %d/%d [%s]\n',intNeuron,intNeurons,getTime);
				hTic=tic;
			end
			clear vecTrialStarts;
			
			%% prep data
			strDate = getDate();
			intSU = intNeuron;
			
			%neuron 1
			dblBaseRate = 0.1;
			dblFirstRate1 = exprnd(3);
			dblFirstRate2 = exprnd(6);
			dblUseMaxDur = 2;
			dblTrialDur = 4;
			vecDurs = [1 1];
			vecRates1 = [dblFirstRate1 0];
			vecRates2 = [dblFirstRate2 0];
			intNumT = intTrialNum;
			vecTrialDur=dblUseMaxDur:dblUseMaxDur:(intNumT*dblUseMaxDur);
			vecRepStarts = 5+cumsum(vecTrialDur);
			dblEndT = vecRepStarts(end)+5;
			vecSpikeTimes1 = getGeneratedMultiPhasicR(dblBaseRate,vecRates1,vecDurs,vecRepStarts,dblEndT);
			if contains(strRunType,'Rand')
				%generate duplicate
				vecSpikeTimes2 = getGeneratedMultiPhasicR(dblBaseRate,vecRates1,vecDurs,vecRepStarts,dblEndT);
			else
				%generate difference
				vecSpikeTimes2 = getGeneratedMultiPhasicR(dblBaseRate,vecRates2,vecDurs,vecRepStarts,dblEndT);
			end
			matEventTimes = vecRepStarts(:);
			matEventTimes(:,2) = matEventTimes(:,1) + dblUseMaxDur;
			
			%% get visual responsiveness
			%get trial dur
			%set derivative params
			intTrials = size(matEventTimes,1);
			intSpikeNum = numel(vecSpikeTimes1)+numel(vecSpikeTimes2);
			[vecTrialPerSpike1,vecTimePerSpike1] = getSpikesInTrial(vecSpikeTimes1,matEventTimes(:,1),dblUseMaxDur);
			[vecTrialPerSpike2,vecTimePerSpike2] = getSpikesInTrial(vecSpikeTimes2,matEventTimes(:,1),dblUseMaxDur);
			if intSpikeNum<3 || numel(vecTimePerSpike1) < 1|| numel(vecTimePerSpike2)< 1,continue;end
			
			%[vecMean,vecSEM,vecWindowBinCenters] = ...
			%	doPEP(vecSpikeTimes,[-0.5:dblBinW:(2.5)],matEventTimes(:,1));
			%pause
			%continue;
			
			%if size(matEventTimes,1) > 0,continue;end
			%%{
			intPlot = 4;
			boolDirectQuantile = false;
			[dblZetaP,sZeta2]=zetatest2(vecSpikeTimes1,matEventTimes,vecSpikeTimes2,matEventTimes,dblUseMaxDur,intResampleNum,intPlot,boolDirectQuantile);
			dblTtestP = sZeta2.dblMeanP;
			
			%% KS
			hTic2 = tic;
			
			vecKsSpikes1 = sort(vecTimePerSpike1);
			vecKsSpikes2 = sort(vecTimePerSpike2);
			
			[h,pKS2] = kstest2(vecKsSpikes1,vecKsSpikes2);
			
			if intPlot > 0
				subplot(2,3,5);cla
				title(sprintf('KS2, p=%.3f',pKS2));
				hold on
				stairs([0 vecKsSpikes1'],[linspace(0,1,1+numel(vecKsSpikes1))]);
				stairs([0 vecKsSpikes2'],[linspace(0,1,1+numel(vecKsSpikes2))]);
				hold off
				ylabel('Cumulative Probability');
				xlabel('Time of spike (s)');
				fixfig;
				
				%save
				drawnow;
				export_fig(fullpath(strFigPath,sprintf('ZETA2_KS2_T%dExample.tif',intTrialNum)));
				export_fig(fullpath(strFigPath,sprintf('ZETA2_KS2_T%dExample.pdf',intTrialNum)));
				return
			end
			%%
			% assign data
			cellNeuron{intTrialNumIdx,intNeuron,intRandType} = [strDate 'N' num2str(intSU)];
			matNumSpikes(intTrialNumIdx,intNeuron,intRandType) = intSpikeNum;
			matZetaP(intTrialNumIdx,intNeuron,intRandType) = dblZetaP;
			matTtestP(intTrialNumIdx,intNeuron,intRandType) = dblTtestP;
			matKsP(intTrialNumIdx,intNeuron,intRandType) = pKS2;
		end
	end
end

if boolSave
	save([strDataTargetPath 'Zeta2DataKsResamp' num2str(intResampleNum) '.mat' ],...
		'cellNeuron','matNumSpikes','matTtestP','matZetaP','matKsP','vecRunTrialNums');
end