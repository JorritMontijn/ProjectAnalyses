function  epochMatrix = MOL_PupilBeh_epoch(Data, sParams)
%% MOL_PupilBeh_epoch this function combine pupil and trial data from different session and extract all epoches
% of specific length

%% COMBINE PUPIL DATA FROM DIFFERENT SESSIONS
pupil = [];
pupilTs = [];
sessionID = 0;

for sessionN = 1: length(Data.pupilData.area)
    pupil = [pupil, Data.pupilData.area{sessionN}'];
    pupilTs = [pupilTs, Data.pupilData.ts{sessionN}'];
    sessionID(length(sessionID):(length(Data.pupilData.area{sessionN})+length(sessionID)-1)) = Data.pupilData.session_ID(sessionN); 
end

%% take a z-score of pupil size data
zPupil = (pupil - nanmean(pupil))/nanstd(pupil);

%% additionall data
secBef = sParams.secBeforeStim;
secAft = sParams.secAfterStim;
stimChangeTs =  Data.trialData.stimChange;
pupilFs = Data.pupilData.fs(1);

pupilSectionL = pupilFs*(secBef+secAft); %how many frames do you consider?

epochMatrix = zeros(length(stimChangeTs), pupilSectionL);

sessionsN = Data.pupilData.session_ID;

%% GET TIME-LOCKED PUPIL EPOCHS
trialCounter = 0;
for numTrial = 1:length(sessionsN)
    sessionTL = length(stimChangeTs(Data.trialData.session_ID==sessionsN(numTrial)));
    for sessionT = 1:sessionTL
        trialCounter = trialCounter+1;
        pupStart = stimChangeTs(trialCounter)-secBef*1e+06;
        pupEnd = stimChangeTs(trialCounter)+secAft*1e+06;
        pupilLogical = pupilTs(sessionID==sessionsN(numTrial))>=pupStart & pupilTs(sessionID==sessionsN(numTrial))<=pupEnd;
        pupilSection = zPupil(pupilLogical);
        
        % pupilSection may not have always have the same omount of frames,
        % correct for this
        if sum(pupilLogical)==0
            epochMatrix(trialCounter, :) = zeros(1,pupilSectionL);
        elseif    sum(pupilLogical) < pupilSectionL
            pupilSection(pupilSectionL) = pupilSection(end);
            epochMatrix(trialCounter, :) = pupilSection;
        elseif sum(pupilLogical) > pupilSectionL
            pupilSection = pupilSection(1:pupilSectionL);
        else
            epochMatrix(trialCounter, :) = pupilSection;
        end
    end
end
end