%get neuronal tuning
if intMouse==8
	%remove last trial
	vecRem = true(size(cellMultiSes{1}.structStim.Orientation));
	vecRem((end-47):end) = false;
	cellMultiSes{1}.structStim = remel(cellMultiSes{1}.structStim,vecRem);
	%recalc dfof
	%cellMultiSes{1} = doRecalcdFoF(cellMultiSes{1},3);
end

%[indTuned,vecNeuronPrefStim,vecOSI,cellAggregate] = getTunedStimDetectionNeurons(cellSes(vecBlock==intPopulation));
[indPresent,vecNeuronPrefStim,vecOSI,cellAggregate] = getTunedStimDetectionNeurons(cellMultiSes{intPopulation},[],boolOnlyPresence);
structStim = cellMultiSes{intPopulation}.structStim;

%remove low-fluorescent astrocytes
%vecMeanF = arrayfun(@(x) mean(x.F),sObject);
%[dummy,vecReorder]=sort(vecMeanF,'ascend');
%sObject(vecReorder(1:round(numel(vecReorder)/3))) = [];

% remove trials with reaction time <100ms
dblRemSecs = 0.1;
indTooFast = (structStim.vecTrialRespSecs-structStim.SecsOn)<dblRemSecs & structStim.vecTrialResponse == 1;
fprintf('Removed %d trials with responses <%dms\n',sum(indTooFast),round(dblRemSecs*1000));
%structStim.vecTrialResponse(indTooFast) = 0;
structStim = remel(structStim,~indTooFast);
%structStim.FrameOff = structStim.FrameOn+1;

%remove trials that were too slow
dblRemSecs = 3;
indTooSlow = (structStim.vecTrialRespSecs-structStim.SecsOn)>dblRemSecs;
fprintf('Removed %d trials with responses >%dms\n',sum(indTooSlow),round(dblRemSecs*1000));
structStim.vecTrialResponse(indTooSlow) = 0;
structStim.vecTrialRespSecs(indTooSlow) = nan;
%structStim.FrameOff = structStim.FrameOn+1;

%remove locomotion trials
if boolExcludeLocomotor
	indLoco = false(1,numel(structStim.Orientation));
	for intTrial=1:numel(structStim.Orientation)
		vecMoveTimes = structStim.cellMoveSecs{intTrial};
		indLoco(intTrial) = sum(vecMoveTimes > structStim.SecsOn(intTrial) & vecMoveTimes < structStim.SecsOff(intTrial)) > 0;
	end
	structStim = remel(structStim,~indLoco);
	fprintf('Removed %d trials with locomotion\n',sum(indLoco));
end

%take opposite directions as the same
structStim.Orientation = mod(structStim.Orientation,180);
vecOrientations = unique(structStim.Orientation);
cellMultiSes{intPopulation}.structStim = structStim;

%get orientation-based trial selection vectors
sTypesOri = getStimulusTypes(structStim,{'Orientation'});
cellSelectOri = getSelectionVectors(structStim,sTypesOri);

%remove non-tuned neurons
cellMultiSes{intPopulation}.neuron(~indPresent) = [];
vecNeuronPrefStim(~indPresent) = [];
intNeurons = numel(cellMultiSes{intPopulation}.neuron);
intTrials = length(structStim.Orientation);
intOris = length(vecOrientations);
%cellMultiSes{intPopulation}.structStim = structStim;
