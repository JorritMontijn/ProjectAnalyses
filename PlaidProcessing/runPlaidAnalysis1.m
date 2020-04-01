
%which set
intType = 1;

%save figs and data?
boolSaveFigs = true;
boolSaveData = true;
	
if intType == 1
	vecRunSets = [1];%[-27:-21 -14 -13 -12 -10 -7:-1];
end

for intLoadSet=vecRunSets
	clearvars -except intSubSample* vecDoShuff vecRunSets intLoadSet bool*
	boolLoad = true;
	
	%% get simulation name [strSimulation] from [intLoadSim]
	loadPlaid;
	
	%% RUN: #header
	strBlockNr = getFlankedBy(mfilename,'Analysis','');
	strBlockNr = strBlockNr(1);
	strFigDir = ['D:\Data\ResultsPlaids\Figures' strBlockNr '\'];
	strDataDir = ['D:\Data\ResultsPlaids\Data' strBlockNr '\'];
	if isempty(strBlockNr),error;end
	
	if boolLoad
		runPlaidHeader;
	end
	
	%% get response matrix
	matResp = getRespMat(sSesAggregate);
	matRespZ = zscore(matResp,[],1); %z-scored at each trial to remove mean population activity
	 
	%% get tuning curves
	%select only gratings
	indGratings = sSesAggregate.structStim.Contrast==100;
	matRespZG = matRespZ(:,indGratings);
	matRespG = matResp(:,indGratings);
	vecStimOriDegrees = sSesAggregate.structStim.Orientation(indGratings);
	
	vecResp = matRespZG(1,:);
	vecAngles = deg2rad(vecStimOriDegrees);
	dblDeltaPrimeZ = getDeltaPrime(matRespZG,vecAngles);
	dblDeltaPrime = getDeltaPrime(matRespG,vecAngles);
	%[dblDeltaPrime,dblSignificance,dblShuffledVariance] = testDeltaPrime(matRespG,vecAngles);
	
	return
end