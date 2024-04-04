function [sUseNeuron,vecStimOnTime,vecStimOffTime,vecOrientation,structStim] = NpxPrepWhisking(sAggNeuron,sThisRec,cellUseAreas)
	%NpxPrepGrating Summary of this function goes here
	%   [sUseNeuron,vecStimOnTime,vecStimOffTime,vecOrientation,structStim] = NpxPrepWhisking(sAggNeuron,sThisRec,cellUseAreas);
	
	%remove stimulus sets that are not 24 stim types
	sThisRec.cellBlock(cellfun(@(x) ~strcmp(x.strExpType,'RunOptoWhiskerStim'),sThisRec.cellBlock)) = [];
	
	%check if any sets exist
	sUseNeuron=[];
	vecStimOnTime=[];
	vecStimOffTime=[];
	vecOrientation=[];
	structStim=[];
	if isempty(sThisRec.cellBlock),return;end
	for i=1:numel(sThisRec.cellBlock)
		sBlock = sThisRec.cellBlock{i};
		intStimNumber = numel(sBlock.vecStimOnTime);
		sBlock.intStimNumber = intStimNumber;
		sBlock.cellPulseData = sBlock.cellPulseData(1:intStimNumber);
		sBlock.cellPulseDelay = sBlock.cellPulseDelay(1:intStimNumber);
		sBlock.cellPulseITI = sBlock.cellPulseITI(1:intStimNumber);
		sBlock.cellPulseDur = sBlock.cellPulseDur(1:intStimNumber);
		sThisRec.cellBlock{i} = sBlock;
	end
	
	% concatenate stimulus structures
	structStim = catstim(sThisRec.cellBlock);
	if strcmp(structStim.sStimParams.strPortOut_whisker,'ao1')
		intOptoStim = 1;
		intWhiskerStim = 2;
	else
		intOptoStim = 2;
		intWhiskerStim = 1;
	end
	
	%select whisker stim trials and create fake trials during spontaneous periods
	vecStimOnTime = structStim.vecStimOnTime;
	vecStimOffTime = structStim.vecStimOffTime;
	intRepNum = numel(vecStimOnTime);
	cellOptoWhiskStim = structStim.cellPulseData;
	cellDelays = structStim.cellPulseDelay;
	dblSamplingRate = structStim.dblSamplingRate;
	vecWhiskT = nan(1,intRepNum);
	for intRep=1:intRepNum
		vecWhiskStim = [false;cellOptoWhiskStim{intRep}(:,intWhiskerStim);false];
		vecOptoStim = [false;cellOptoWhiskStim{intRep}(:,intOptoStim);false];
		vecDelays = cellDelays{intRep};
		
		[dummy,vecPeaksWhisk]=findpeaks(double(vecWhiskStim));
		[dummy,vecPeaksOpto]=findpeaks(double(vecOptoStim));
		vecWhiskTimes = (vecPeaksWhisk-1)./dblSamplingRate;
		vecOptoTimes = (vecPeaksOpto-1)./dblSamplingRate;
		vecRealDelays = vecWhiskTimes - vecOptoTimes;
		if ~all((vecDelays(:) - vecRealDelays(:)) < 1e-3)
			error uh-oh
		end
		intUseStim = find(vecDelays==1);
		vecWhiskT(intRep) = vecWhiskTimes(intUseStim) + vecStimOnTime(intRep);
	end
	vecNoWhiskT = vecWhiskT - 1;
	vecStimOnTime = cat(2,vecNoWhiskT,vecWhiskT);
	vecOrientation = cat(2,zeros(size(vecNoWhiskT)),ones(size(vecWhiskT)));
	[vecStimOnTime,vecReorder]=sort(vecStimOnTime);
	vecOrientation = vecOrientation(vecReorder);
	vecStimOffTime = vecStimOnTime+1;
	
	%remove neurons from other recs
	strRec = sThisRec.Exp;
	indQualifyingNeurons = contains({sAggNeuron.Exp},strRec);
	
	%remove neurons in incorrect areas
	indConsiderNeurons = contains({sAggNeuron.Area},cellUseAreas,'IgnoreCase',true);
	
	%remove bad clusters
	%indGoodNeurons = ((cell2vec({sAggNeuron.KilosortGood}) == 1) | (cell2vec({sAggNeuron.Contamination}) < 0.1))...
	%	& (abs(cell2vec({sAggNeuron.NonStationarity})) < 0.25);
	indGoodNeurons = ~strcmp({sAggNeuron.bc_unitType},'NOISE');
	
	%subselect from total
	indUseNeurons = indQualifyingNeurons(:) & indConsiderNeurons(:) & indGoodNeurons(:);
	sUseNeuron = sAggNeuron(indUseNeurons);
end

