function [matMeanRate,indTuned,cellSpikeTimes,sOut,cellSpikeTimesPerCellPerTrial,...
		vecStimOnStitched,vecNonStat,dblBC,dblMaxDevFrac,indResp] = ...
		NpxPrepData(cellSpikeTimesIn,vecStimOnTime,vecStimOffTime,vecStimTypes)
	%NpxPrepData Summary of this function goes here
	%   [matMeanRate,indTuned,cellSpikeTimes,sOut,cellSpikeTimesPerCellPerTrial,vecStimOnStitched,vecNonStat,dblBC,dblMaxDevFrac] = ...
	%		NpxPrepData(cellSpikeTimesIn,vecStimOnTime,vecStimOffTime,vecStimTypes)
	
	%% set cut-offs
	%get dur
	cellSpikeTimes = cellSpikeTimesIn;
	dblDur = median(vecStimOffTime-vecStimOnTime);
	dblMinRate = 0.1;%2/(numel(vecStimOnTime)*dblDur);
	matRawData = getSpikeCounts(cellSpikeTimes,vecStimOnTime,dblDur)./dblDur;
	indResp = sum(matRawData,2)'>(size(matRawData,2)/dblDur)*dblMinRate;
	matMeanRate = matRawData(indResp,:);
	
	%% remove non-responsive cells
	[vecStimIdx,vecUnique,vecPriorDistribution,cellSelect,vecRepetition] = val2idx(vecStimTypes);
	sOut = getTuningCurves(matMeanRate,vecStimTypes,0);
	cellSpikeTimes(~indResp)=[];
	indTuned = sOut.vecOriAnova<0.05;
	intTrialNum = numel(vecStimTypes);
	intRespN = size(matMeanRate,1);
	intStimNum = numel(vecUnique);
	intRepNum = min(vecPriorDistribution);
	
	%% check non-stationarity
	%get spikes per trial per neuron
	cellSpikeTimesPerCellPerTrial = cell(intRespN,intTrialNum);
	vecNonStat = nan(1,intRespN);
	boolDiscardEdges = true;
	vecPseudoEventT = [];
	for intN=1:intRespN
		% build pseudo data, stitching stimulus periods
		[vecPseudoSpikeTimes,vecPseudoEventT] = getPseudoSpikeVectors(cellSpikeTimes{intN},vecStimOnTime,dblDur,boolDiscardEdges);
		
		%real
		[vecTrialPerSpike,vecTimePerSpike] = getSpikesInTrial(vecPseudoSpikeTimes,vecPseudoEventT,dblDur);
		for intTrial=1:intTrialNum
			vecSpikeT = vecTimePerSpike(vecTrialPerSpike==intTrial);
			cellSpikeTimesPerCellPerTrial{intN,intTrial} = vecSpikeT;
		end
		
		%calc non-stationarity
		vecSortedSpikeTimes = sort(vecPseudoSpikeTimes,'ascend') - min(vecPseudoSpikeTimes);
		dblAUC = sum(vecSortedSpikeTimes);
		dblLinAUC = (max(vecSortedSpikeTimes) * numel(vecSortedSpikeTimes) ) / 2;
		vecNonStat(intN) = (dblAUC - dblLinAUC) / dblLinAUC;
	end
	vecStimOnStitched = vecPseudoEventT;
	matDataZ = zscore(log(1+matMeanRate),[],2);
	vecMeanZ = mean(matDataZ,1);
	vecFilt = normpdf(-4:4,0,1)/sum(normpdf(-4:4,0,1));
	vecFiltM = imfilt(vecMeanZ,vecFilt);
	
	if intRespN > 0
	%calc metrics
		[h,pKS,ksstat,cv] = kstest(vecFiltM);
		[BF, dblBC] = bimodalitycoeff(vecFiltM);
		dblMaxDevFrac = max(abs(vecFiltM));
	else
		dblBC = nan;
		dblMaxDevFrac = nan;
	end
end
