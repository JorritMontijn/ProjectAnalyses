function [sUseNeuron,vecStimOnTime,vecStimOffTime,vecOrientation,structStim] = NpxPrepGrating(sAggNeuron,sThisRec,cellUseAreas)
	%NpxPrepGrating Summary of this function goes here
	%   [sUseNeuron,vecStimOnTime,vecStimOffTime,vecOrientation,structStim] = NpxPrepGrating(sAggNeuron,sThisRec,cellUseAreas);
	
	%remove stimulus sets that are not 24 stim types
	sThisRec.cellBlock(cellfun(@(x) x.intTrialNum/x.intNumRepeats,sThisRec.cellBlock) ~= 24) = [];
	sThisRec.cellBlock(cellfun(@(x) ~strcmp(x.strExpType,'RunDriftingGrating'),sThisRec.cellBlock)) = [];
	
	% concatenate stimulus structures
	structStim = catstim(sThisRec.cellBlock(1:end));
	vecStimOnTime = structStim.vecStimOnTime;
	vecStimOffTime = structStim.vecStimOffTime;
	vecOrientation = cell2vec({structStim.sStimObject(structStim.vecTrialStimTypes).Orientation})';
	vecTempFreq = cell2vec({structStim.sStimObject(structStim.vecTrialStimTypes).TemporalFrequency})';
	vecPhase = structStim.Phase;
	vecDelayTimeBy = vecPhase./vecTempFreq;
	[vecOriIdx,vecUniqueOris,vecRepNum,cellSelect,vecTrialRepetition] = val2idx(vecOrientation);
	indRem=vecTrialRepetition>min(vecRepNum);
	vecOrientation(indRem) = [];
	vecDelayTimeBy(indRem) = [];
	vecStimOnTime(indRem) = [];
	vecStimOffTime(indRem) = [];
	
	%remove neurons from other recs
	strRec = sThisRec.Exp;
	indQualifyingNeurons = contains({sAggNeuron.Exp},strRec);
	
	%remove neurons in incorrect areas
	indConsiderNeurons = contains({sAggNeuron.Area},cellUseAreas,'IgnoreCase',true);
	
	%remove bad clusters
	indGoodNeurons = ((cell2vec({sAggNeuron.KilosortGood}) == 1) | (cell2vec({sAggNeuron.Contamination}) < 0.1))...
		& (abs(cell2vec({sAggNeuron.NonStationarity})) < 0.25);
	
	%subselect from total
	indUseNeurons = indQualifyingNeurons(:) & indConsiderNeurons(:) & indGoodNeurons(:);
	sUseNeuron = sAggNeuron(indUseNeurons);
	
end

