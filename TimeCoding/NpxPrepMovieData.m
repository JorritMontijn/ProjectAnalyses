function [matMeanRate,cellSpikeTimes,cellSpikeTimesPerCellPerTrial,...
		vecStimOnStitched,vecNonStat,dblBC,dblMaxDevFrac] = ...
		NpxPrepMovieData(cellSpikeTimesIn,vecStimOnTime,vecStimOffTime,vecStimTypes)
	%NpxPrepMovieData Summary of this function goes here
	%   [matMeanRate,cellSpikeTimes,cellSpikeTimesPerCellPerTrial,vecStimOnStitched,vecNonStat,dblBC,dblMaxDevFrac] = ...
	%		NpxPrepMovieData(cellSpikeTimesIn,vecStimOnTime,vecStimOffTime,vecStimTypes)
	
	%% set cut-offs
	dblMinRate = 0.1;
	
	%get dur
	cellSpikeTimes = cellSpikeTimesIn;
	dblDur = median(vecStimOffTime-vecStimOnTime);
	matRawData = getSpikeCounts(cellSpikeTimes,vecStimOnTime,dblDur)./dblDur;
	indResp = sum(matRawData,2)'>(size(matRawData,2)/dblDur)*dblMinRate;
	
	%% remove non-responsive cells
	[vecStimIdx,vecUnique,vecPriorDistribution,cellSelect,vecRepetition] = val2idx(vecStimTypes);
	matMeanRate = matRawData(indResp,:);
	cellSpikeTimes(~indResp)=[];
	intTrialNum = numel(vecStimTypes);
	intRespN = size(matMeanRate,1);
	intStimNum = numel(vecUnique);
	intRepNum = min(vecPriorDistribution);
	
	%% check non-stationarity
	if nargout < 3
		return
	end
	%get spikes per trial per neuron
	cellSpikeTimesPerCellPerTrial = cell(intRespN,intTrialNum);
	vecNonStat = nan(1,intRespN);
	boolDiscardEdges = true;
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
	
	%calc metrics
	[h,pKS,ksstat,cv] = kstest(vecFiltM);
	[BF, dblBC] = bimodalitycoeff(vecFiltM);
	dblMaxDevFrac = max(abs(vecFiltM));
	
