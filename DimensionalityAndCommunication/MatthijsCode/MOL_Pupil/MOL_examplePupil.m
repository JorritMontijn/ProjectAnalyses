%% MOL_Opto_V1PPC_Behavior
%% Get input arguments:

[Data] = MOL_GetData('E:','CHDET',{'ChangeDetectionConflict' 'BehaviorConflict'},{'2012'},{'2018-08-14_14-30-15'},{'sessionData' 'trialData' 'pupilData'});
sessionData     = Data.sessionData;
trialData       = Data.trialData;
pupilData       = Data.pupilData;

%% Remove last 20 trials:
trialData = MOL_RemoveLastnTrials(trialData,20);


%%
figure;
exampletrial = 35;
trialData.stimChange(exampletrial)
idx = pupilData.ts{1}>trialData.stimChange(exampletrial)-1e6 & pupilData.ts{1}<trialData.stimChange(exampletrial)+5e6;

plot(pupilData.ts{1}(idx),pupilData.area{1}(idx),'k','LineWidth',2);
xlabel('Time (in us)')
ylabel('Pupil area (in pixels)')



