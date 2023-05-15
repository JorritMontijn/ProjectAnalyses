function [sUseNeuron,vecStimOnTime,vecStimOffTime,vecStimIdx,structStim] = NpxPrepMovies(sAggNeuron,sThisRec,cellUseAreas)
	%NpxPrepGrating Summary of this function goes here
	%   [sUseNeuron,vecStimOnTime,vecStimOffTime,vecStimIdx,structStim] = NpxPrepMovies(sAggNeuron,sThisRec,cellUseAreas);
	
	%remove stimulus sets that are not 100 reps
	sThisRec.cellBlock(cellfun(@(x) x.intNumRepeats,sThisRec.cellBlock) < 100) = [];
	
	% concatenate stimulus structures
	structStim = catstim(sThisRec.cellBlock(1));
	vecOrigStimOnTime = structStim.vecStimOnTime;
	dblDur = median(diff(vecOrigStimOnTime));
	vecOrigStimOffTime = vecOrigStimOnTime+dblDur;
	intFramesInMovie = 500;
	dblBinAvg = 5;
	intUseBinsInMovie = intFramesInMovie/dblBinAvg;
	dblBinRate = round(intUseBinsInMovie/dblDur);
	dblBinDur = 1/dblBinRate;
	vecBinEdges = linspace(0,dblBinDur*intUseBinsInMovie,intUseBinsInMovie+1);
	
	%generate fake stimulus vectors
	vecUniqueStims = 1:intUseBinsInMovie;
	vecFrameIdx = repmat(vecUniqueStims,[1 numel(vecOrigStimOnTime)]);
	vecStimOnTime = flat(vecBinEdges(1:(end-1))' + vecOrigStimOnTime)';
	vecStimOffTime = vecStimOnTime + dblBinDur;
	
	[vecStimIdx,vecUniqueStims,vecRepNum,cellSelect,vecTrialRepetition] = val2idx(vecFrameIdx);
	
	%remove neurons from other recs
	strRec = sThisRec.Exp;
	indQualifyingNeurons = contains({sAggNeuron.Exp},strRec);
	
	%remove neurons in incorrect areas
	indConsiderNeurons = contains({sAggNeuron.Area},cellUseAreas,'IgnoreCase',true);
	
	%remove bad clusters
	indGoodNeurons = (cell2vec({sAggNeuron.KilosortGood}) == 1) | (cell2vec({sAggNeuron.Contamination}) < 0.1);
	
	%subselect from total
	indUseNeurons = indQualifyingNeurons(:) & indConsiderNeurons(:) & indGoodNeurons(:);
	sUseNeuron = sAggNeuron(indUseNeurons);
	
end

