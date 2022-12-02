%set params
%stim params
vecOris = deg2rad(0:15:359);
intReps = 10;
dblTrialT = 1.5;
dblStimDur = 1;

%tuning params
dblBaseRate = 1;
dblPrefRate = 2;
dblKappa = 5;
boolDoublePeaked = true;
dblPrefOri = 2*pi*rand(1);

%peak params
dblJitter = 0.05;
intAddSpikes = 0;
dblStartDelay = 0.1;

%derived variables
vecTrialAngles = repmat(vecOris,[1 intReps]);
intTrialNum = numel(vecTrialAngles);
vecStimStartT = 1 + 0:dblTrialT:(dblTrialT*intTrialNum);
vecStimStopT = vecStimStartT + dblStimDur;
matTrialT = cat(1,vecStimStartT,vecStimStopT)';




%generate data
[vecSpikeTimes,dblPrefOri] = getGeneratedSpikingDataWithPeak(...
	vecTrialAngles,matTrialT,dblBaseRate,dblPrefRate,dblJitter,dblKappa,boolDoublePeaked,...
	dblPrefOri,intAddSpikes,dblStartDelay);

[vecSpikeTimes2,dblPrefOri] = getGeneratedSpikingDataWithPeak(...
	vecTrialAngles,matTrialT,dblBaseRate,dblBaseRate,dblJitter,dblKappa,boolDoublePeaked,...
	dblPrefOri,intAddSpikes,dblStartDelay);

%calc zeta
[vecSpikeT,vecRealDiff,vecRealFrac,vecRealFracLinear,cellRandT,cellRandDiff,dblZetaP,dblZETA,intZETALoc] = ...
	calcZetaOne(vecSpikeTimes,matTrialT,dblTrialT,...
	250,false,2,true);

[vecSpikeT2,vecRealDiff2,vecRealFrac2,vecRealFracLinear2,cellRandT2,cellRandDiff2,dblZetaP2,dblZETA2,intZETALoc2] = ...
	calcZetaOne(vecSpikeTimes2,matTrialT,dblTrialT,...
	250,false,2,true);

%2-sample zeta test
[dblZetaP,sZETA] = zetatest2(vecSpikeTimes,matTrialT,vecSpikeTimes2,matTrialT,false);


%
[vecSpikeTC,vecRealDiffC,vecRealFracC,vecRealFracLinearC,cellRandTC,cellRandDiffC,dblZetaPC,dblZETAC,intZETALocC] = ...
	calcZetaOne(sort(cat(1,vecSpikeTimes,vecSpikeTimes2)),matTrialT,dblTrialT,...
	250,false,2,true);
