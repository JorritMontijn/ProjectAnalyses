%{
1) Meerdink_20200207_001_Split1_normcorr_SPSIG?
2) Meerdink_20200207_001_Split1_normcorr_SPSIG_Res

In 1) staan 'sigCorrected' (df/f x neurons) en 'deconCorrected' (spikes x neurons). Spikes zijn #spikes per frame, niet spiketimes. DeconCorrected moet ik dus nog even beter naar kijken.
In 2) staan 'Res' en 'info'.

'Res' zijn matrices van verschillende signalen (frames x trials x cells) en 'ax', een as in tijd per trial (t = -2 tot aan t = 7, t = 0 is stimonset).

In 'info' staan relevante dingen voor jou:

info.StimTimes
Tijd van stim onset voor imaging data

info.Stim.log.motionlog
(:,1) = rensnelheid in cm/s
(:,2) = tijd 

info.Stim.log.licklog
(:,1) = detector identifier (irrelevant voor jou)
(:,2) = lickTimes

info.Stim.log.stimlog
(:,1) = stimOnset (is dus in tijd van stimulus pc voor gedragsdata, komt dus niet overeen met info.StimTimes, maar kun je gebruiken om ze te alignen met elkaar)
(:,2) = stimOffset (zou 2 seconden later moeten zijn)
(:,3) = image type (1=rewarded, 0=not rewarded, als in go versus no-go, hoeft dus niet per se rewarded te zijn, dit staat puur voor of het een go of no-go image is)
(:,4) = type of trial (1=active, 2=passive, 0=NA: not rewarded image, in deze recording zijn er als het goed is bijna alleen maar active trials en NA trials)
(:,5) = hit or miss (1 = hit, 0=miss/passive trial)
(:,6) = reward time
(:,7) = stimulus type (1:16, er zijn vier plaatjes, met elk vier condities (alleen hier zijn al die condities FullScreen, dus er is geen verschil), 1:4 = plaatje 1, 5:8 = plaatje 2, etc.)

Timeline voor de taak is:
t = -2:0 > muis mag niet likken
t = 0:2 > stimulus
t = 2:3 > wait (muis mag likken, maar gebeurt niks)
t = 3 > als passive trial: reward. Als active trial: start lickwindow
t = 5 > als active trial: einde lickwindow
%}

strPath = 'D:\Data\Processed\imagingKoen\';
strFile1 = 'Meerdink_20200207_001_Split1_normcorr_SPSIG.mat';
strFile2 = 'Meerdink_20200207_001_Split1_normcorr_SPSIG_Res.mat';

sLoad1 = load(fullfile(strPath,strFile1));
sLoad2 = load(fullfile(strPath,strFile2));

matAct = sLoad1.sigCorrected-1;
sInfo = sLoad2.info;

clear sLoad1;
clear sLoad2;

vecFrameTimes = sInfo.Frametimes;
vecStimTimes = sInfo.StimTimes;
vecLickTimes = sInfo.Stim.log.licklog(:,2);

vecRunTimes = sInfo.Stim.log.motionlog(:,2);
vecRunSpeed = sInfo.Stim.log.motionlog(:,1);

%get data
matStimLog = sInfo.Stim.log.stimlog(1:numel(vecStimTimes),:);
vecStimOn = matStimLog(:,1);
vecStimOff = matStimLog(:,2);
vecImageType = matStimLog(:,3); %1=rewarded,0=non-rewarded
vecTrialType = matStimLog(:,4); %1=active, 2=passive, 0=NA: not rewarded image, in deze recording zijn er als het goed is bijna alleen maar active trials en NA trials
vecHitMiss = matStimLog(:,5);% hit or miss (1 = hit, 0=miss/passive trial)
vecRewardTime = matStimLog(:,6);% reward time
vecStimTypeOrig= matStimLog(:,7);% stimulus type (1:16, er zijn vier plaatjes, met elk vier condities (alleen hier zijn al die condities FullScreen, dus er is geen verschil), 1:4 = plaatje 1, 5:8 = plaatje 2, etc.)
vecStimType = ceil(vecStimTypeOrig./4);
vecStimTypes = unique(vecStimType); %transform

%correct stim times
dblPreTime = 0;
%
vecDiff = vecStimTimes - vecStimOn;
vecStimOn = vecStimOn + vecDiff;
vecStimOff = vecStimOff + vecDiff;
matStimOnOff = [vecStimOn vecStimOff] - dblPreTime;

%%
vecZeta = nan(1,intNeurons);
vecZetaP = nan(1,intNeurons);
vecHzP = nan(1,intNeurons);
matLatencies = nan(3,intNeurons);
intNeurons = size(matAct,2);
dblTrialDur = 7;%min(diff(vecStimOn));
for intNeuron=1:intNeurons
	intNeuron
	vecTrace=matAct(:,intNeuron);
	[dblZETA,vecLatencies,sZETA,sMSD] = getTraceZeta(vecFrameTimes,vecTrace,matStimOnOff,dblTrialDur,50,0,3);
	vecZeta(intNeuron) = dblZETA;
	vecZetaP(intNeuron) = sZETA.dblP;
	vecHzP(intNeuron) = sZETA.dblHzP;
	matLatencies(:,intNeuron) = vecLatencies - dblPreTime;
end
	
%% neuron 42, 81
[dblZETA,vecLatencies,sZETA,sMSD] = getTraceZeta(vecFrameTimes,matAct(:,42),vecStimOn,dblTrialDur,50,2,3);