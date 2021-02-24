%% load spiking data & plot tuning curves

%% set recording
close all;
clear all;

strDisk = 'F:';
strDataTargetPath = [strDisk '\Data\Processed\ZETA\Latencies\'];
strFigPath = [strDisk '\Data\Results\ZETA\Latencies\'];
vecBinDurs = (1.5.^(0:10))/1000;
intMakePlots = 0; %0=none, 1=normal plot, 2=including raster
vecRunTypes = 1;%[1 2];
boolSave = true;
vecBaseRate = 2^5;%2.^(-1:1:5);
vecJitters = [1:1:10];
intNeurons = 10;
intTrialReps = 20;
intBinNum = numel(vecBinDurs);

%set var
for intJitterIdx=1:numel(vecJitters)
	dblJitter = roundi(vecJitters(intJitterIdx),1);
	for intBaseRateIdx=1:numel(vecBaseRate)
		dblBaseRate = roundi(vecBaseRate(intBaseRateIdx),1);
		strArea = sprintf('PoissonPeakRate%02.1fJitter%02.1f',dblBaseRate,dblJitter);
		fprintf('Processing %s [%s]\n',strArea,getTime);
		%reset vars
		clearvars -except intTrialReps vecBinDurs intNeurons intBinNum dblJitter intJitterIdx vecJitters dblBaseRate strArea intBaseRateIdx vecBaseRate vecResamples boolSave strDataTargetPath strFigPath intMakePlots vecRunTypes
		hTic = tic;
		for intRunType=vecRunTypes
			if intRunType == 1
				strRunType = strArea;
			elseif intRunType ==2
				strRunType = [strArea '-Rand'];
			end
			
			%% pre-allocate output variables
			vecLatencies = nan(1,intNeurons);
			matBinLatencies = nan(intBinNum,intNeurons);
			vecRealLatencies = nan(1,intNeurons);
			vecComputIFR = nan(1,intNeurons);
			vecLatenciesM = nan(1,intNeurons);
			vecComputM = nan(1,intNeurons);
				
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
				vecTrialAngles=deg2rad(repmat([0:45:359],[1 intTrialReps]));
				dblTrialDur=2;
				vecStimOnTime = dblTrialDur*(1:numel(vecTrialAngles))';
				vecStimOffTime = vecStimOnTime + 1;
				
				matEventTimes(:,1) = vecStimOnTime;
				matEventTimes(:,2) = vecStimOffTime;
				[vecSpikeTimes,dblPrefOri] = getGeneratedSpikingDataWithPeak(vecTrialAngles,matEventTimes,dblBaseRate,dblPrefRate,dblJitter,dblKappa,true);
				dblRandOnset = 0.01*((rand(1)*2) - 1);
				vecSpikeTimes = vecSpikeTimes + dblRandOnset;
				
				%% get visual responsiveness
				if intNeuron == 1 && intMakePlots > 0
					[dblZeta,vecLatencies,sZETA,sMSD] = getZeta(vecSpikeTimes,vecStimOnTime,[],[],intMakePlots);
				end
				
				%get IFR
				hTic=tic;
				[vecIFR,sIFR] = getIFR(vecSpikeTimes,matEventTimes,[],3,[],1.2);
				[dblOnsetVal,dblOnset] = getPeak(vecIFR,sIFR.vecT,[0 1]);
				dblComputIFR=toc(hTic);
				vecLatencies(intNeuron) = dblOnset;
				vecRealLatencies(intNeuron) = 0.1+dblRandOnset;
				vecComputIFR(intNeuron) = dblComputIFR;
				%get mimi
				hTic=tic;
				[dblMIMI_P,vecLatM,sMIMI,sRate] = getMIMI(vecSpikeTimes,matEventTimes);
				dblComputMIMI=toc(hTic);
				vecLatenciesM(intNeuron) = vecLatM(1);
				vecComputM(intNeuron) = dblComputMIMI;
				
				%% get bin-wise approach
				%get data
				for intBinIdx=1:intBinNum
					dblFrameDur = vecBinDurs(intBinIdx);
					dblStimDur = median(diff(vecStimOnTime));
					vecBinEdges = 0:dblFrameDur:2;
					vecBinCenters = vecBinEdges(2:end)-dblFrameDur/2;
					intBins = numel(vecBinCenters);
					intTrials = numel(vecStimOnTime);
					matResp = nan(intTrials,intBins);
					for intTrial=1:intTrials
						matResp(intTrial,:) = histcounts(vecSpikeTimes,vecBinEdges+vecStimOnTime(intTrial,1));
					end
					
					%test
					vecR = mean(matResp,1);
					[dblOnsetVal,dblOnset] = getPeak(vecR,vecBinCenters,[0 1]);
					matBinLatencies(intBinIdx,intNeuron) = dblOnset;
				end
			end
			%save
			if boolSave
				save([strDataTargetPath 'ZetaDataPPL3b' strRunType '.mat' ],...
					'dblBaseRate','dblJitter','vecBinDurs','vecLatencies','matBinLatencies','vecRealLatencies',...
					'vecComputIFR','vecLatenciesM','vecComputM');
			end
		end
	end
end