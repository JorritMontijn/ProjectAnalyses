
%which set
intType = 1;

%save figs and data?
if intType == 1
	vecRunSets = [1:8];%[-27:-21 -14 -13 -12 -10 -7:-1];
end
figure
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
	
	vecAngles = deg2rad(vecStimOriDegrees);
	sOut = getTuningCurves(matRespZG,vecStimOriDegrees);
	
	%% plot
	matC = lines(2);
	vecPrefOri = sOut.matFittedParams(:,1)/2;
	vecMeanR = mean(matRespG,2);
	vecSignificant = sOut.vecOriTtest < 0.05;
	subplot(3,3,intLoadSet)
	hold on
	scatter(rad2deg(vecPrefOri(vecSignificant)),vecMeanR(vecSignificant),[],matC(1,:),'x');
	scatter(rad2deg(vecPrefOri(~vecSignificant)),vecMeanR(~vecSignificant),[],matC(2,:),'x');
	hold off
	xlabel('Pref ori (deg)');
	ylabel('Mean resp (dF/F0)');
	fixfig;
end
legend(subplot(3,3,8),{'p<0.05','p>0.05'},'Location','Best');
fixfig;
maxfig;