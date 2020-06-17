%{
 Je vind er twee .mat files:
1) Delier_20191017_002_Split1_normcorr_SPSIG
2) Delier_20191017_002_Split1_normcorr_SPSIG_Res

In 1) staan 'sigCorrected' (df/f x neurons) en 'deconCorrected' (spikes x neurons). Spikes zijn #spikes per frame, niet spiketimes. DeconCorrected moet ik dus nog even beter naar kijken.
In 2) staan 'Res' en 'info'.

'Res' zijn matrices van verschillende signalen (frames x trials x cells) en 'ax', een as in tijd per trial (t = -2 tot aan t = 7, t = 0 is stimonset).

In 'info' staan relevante dingen voor jou:

info.StimTimes
Tijd van stim onset voor imaging data

info.Stim.log.licklog
(:,1) = detector identifier (irrelevant voor jou)
(:,2) = lickTimes

info.Stim.log.motionlog
(:,1) = rensnelheid in cm/s
(:,2) = tijd
(:,3) = afstand in de virtuele tunnel

info.Stim.log.stimlog
(:,1) = stimOnset/trialOnset (is dus in tijd van stimulus pc voor gedragsdata, komt dus niet overeen met info.StimTimes, maar kun je gebruiken om ze te alignen met elkaar, dit is de begin van de tunnel, distance = 0).
(:,2) = rewardOnset (moment van beloningsafgifte, is 3 seconde na grayOnset)
(:,3) = grayOnset (tijd waarop muis door einde van de tunnel (distance = 0.5) loopt. Is dus variabel per trial)
(:,4) = mismatchOnset (tijd waarop er een mismatch plaatsvindt, zie hieronder)
(:,5) = mismatchLocation (distance waarop er een mismatch plaatsvindt, zie hieronder)

In diezelfde folder vind je een .ppt met 1 slide die wat info over de tunnel geeft. Tunnel is 200 cm lang, maar muis ziet maar 100 cm, dat is 0.5 in info.Stimlog.motionlog(:,3), daarna wordt het grijs voor de muis. Zodra de muis door 0.5 distance loopt (100 cm in .ppt) gaat de beloningsfase in. 1 seconde na de finish hoort de muis een toon, 2 seconde daarna krijgt hij een beloning. Hij mag gewoon doorlopen (gebeurt niks in de tunnel want het is grijs, ze zien dus geen beweging), maar dat doen ze meestal niet. 6 sec na beloning wordt hij terug getransporteerd naar de start.
In een aantal trials druk ik tijdens de taak op de mismatch knop, en dan stop ik 0.5 sec met renderen van de tunnel, dan staat de tunnel dus stil ongeacht of de muis rent -> mismatch. Cellen reageren dan vaak zeer sterk. In deze dataset drukte ik pas in trial 100 voor het eerst op de mismatch knop.

Hoop dat je er wat aan hebt. Ik hoor het wel als je vragen hebt.


 Wat betreft de locaties van de stimuli ga ik er weer even vanuit dat de corridor 200 cm lang is, maar dat de muizen maar 100 cm lang stimuli zien (t/m distance = 0.5), daarna is het dus grijs. Het midden van de stimuli bevindt zich dan op locaties:
[22.2200   33.3300   44.4440   55.5500   66.6600   77.7770] (dus bv de eerste stimulus bevindt zich op 22.22% van 0.5)

De grootte vd stimuli is maar 0.027 op dezelfde schaal (0.027% van 0.5), dus dat is zeer klein en ik weet niet of je daar rekening mee wilt houden.
%}

%% load data
clear all;
strDataPath = 'F:\Data\Processed\VirtualTunnel\Delier\';
strDataFile = 'DelierPreProSpikes.mat';
load([strDataPath strDataFile]);


%% transform times
licktimes = sInfo.Stim.log.licklog(:,2); %in frame or stim time?

vecStartT = sInfo.StimTimes(:,1);
intTrials = numel(vecStartT);
vecLocalT = sInfo.Stim.log.stimlog(1:intTrials,1);
vecRewardT = sInfo.Stim.log.stimlog(1:intTrials,2);
vecRewardT(sInfo.Stim.log.stimlog(1:intTrials,2)==0) = [];
vecGreyT = sInfo.Stim.log.stimlog(1:intTrials,3);
vecGreyT(sInfo.Stim.log.stimlog(1:intTrials,3)==0) = [];
vecMismatchTrials = find(sInfo.Stim.log.stimlog(1:intTrials,4)~=0);
vecNormalTrials = find(sInfo.Stim.log.stimlog(1:intTrials,4)==0);
vecMismatchT = sInfo.Stim.log.stimlog(1:intTrials,4);
vecMismatchT(vecNormalTrials) = [];
vecMismatchLoc = sInfo.Stim.log.stimlog(1:intTrials,5);
vecMismatchLoc(sInfo.Stim.log.stimlog(1:intTrials,5)==0) = [];
intMismatchNum = numel(vecMismatchLoc);

%% select control locations for mismatch at same position
cellControlT = cell(1,intMismatchNum);
vecAllLocs = sInfo.Stim.log.motionlog(:,3);
vecAllT = sInfo.Stim.log.motionlog(:,2);
dblTrialLocLength = 2;
vecTrialNr = sum(bsxfun(@gt,vecAllT',vecLocalT),1);%option 1
%vecTrialNr = sum(bsxfun(@gt,vecAllT',vecStartT),1);%option 2
vecTrialLocStarts = vecTrialNr(:).*dblTrialLocLength;
vecEventOn = unique(vecTrialLocStarts);
vecAllTrialLocs = vecAllLocs + vecTrialLocStarts;
vecMismatchTrialLocs = vecMismatchLoc + vecMismatchTrials.*dblTrialLocLength;
vecRealMismatchT = vecMismatchT + vecStartT(vecMismatchTrials) - vecLocalT(vecMismatchTrials);
vecRealNormalT = vecStartT(vecNormalTrials);

%% calculate visually responsive cells
intNeurons = numel(cellSpikeTimes);
intResampNum = 250;
intLatencyPeaks = 4;
hTic=tic;
vecStimLocs = ([22.2200   33.3300   44.4440   55.5500   66.6600   77.7770]./100)*2;

intIters=100;
vecVisLocZetaP = nan(intIters,intNeurons);
vecVisTimZetaP = nan(intIters,intNeurons);
vecMisMatZetaP = nan(1,intNeurons);
intPlot = 0;

for intNeuron=1:intNeurons
		if toc(hTic) > 5
			hTic = tic;
			fprintf('Neuron %d/%d [%s]\n',intNeuron,intNeurons,getTime);
		end
	%%
	for intIter=1:intIters
		vecRandNormalTrials = sort(vecNormalTrials(randperm(numel(vecNormalTrials),numel(vecRealMismatchT))));
	vecSpikeLocations = interp1(vecAllT,vecAllTrialLocs,cellSpikeTimes{intNeuron});
	%location normal
	vecNormalOn = vecEventOn(vecRandNormalTrials);
	[dblZetaP,vecLatencies,sZETA,sMSD] = getZeta(vecSpikeLocations,vecNormalOn,0.5,intResampNum,intPlot,intLatencyPeaks);
	vecVisLocZetaP(intIter,intNeuron) = sZETA.dblZETA;
	
	%time normal
	[dblZetaP,vecLatencies,sZETA,sMSD] = getZeta(cellSpikeTimes{intNeuron},vecStartT(vecRandNormalTrials),2,intResampNum,intPlot,intLatencyPeaks);
	vecVisTimZetaP(intIter,intNeuron) = sZETA.dblZETA;
	end
	
	%mismatch
	[dblZetaP,vecLatencies,sZETA,sMSD] = getZeta(cellSpikeTimes{intNeuron},vecRealMismatchT,2,intResampNum,intPlot,intLatencyPeaks);
	vecMisMatZetaP(intNeuron) = sZETA.dblZETA;
	
end
vecVisLocZetaP = mean(vecVisLocZetaP,1);
vecVisTimZetaP = mean(vecVisTimZetaP,1);

%{
		Time-locked		Location-locked		Mismatch-mod
Time-L	
Loc-L	r=0.24,p=0.03
MM-Mod	r=0.005,p=0.96	r=-0.34,p=0.002		

Time-L, N=57
Loc-L, N=64
MM-Mod, N=31
%}
%% plot
figure;
subplot(2,3,1)
scatter(vecVisLocZetaP,vecVisTimZetaP,[],vecMisMatZetaP)
%scatter(vecVisLocZetaP,vecVisTimZetaP,[],log10(vecMisMatZetaP))
colormap('bluepurplered');
xlabel('Location-modulation (ZETA)');
ylabel('Time-modulation (ZETA)');
h=colorbar;
h.Label.String = 'Mismatch-modulation (ZETA)';
fixfig;
title(sprintf('%d neurons,%d mismatch trials',numel(vecVisLocZetaP),numel(vecRealMismatchT)));
%set(gca,'xscale','log','yscale','log')

%[rLT,pLT]=corr(log10(vecVisLocZetaP)',log10(vecVisTimZetaP'));
[rLT,pLT]=corr(vecVisLocZetaP',vecVisTimZetaP');
subplot(2,3,2)
scatter(vecVisLocZetaP,vecVisTimZetaP)
xlabel('Location-modulation (ZETA)');
ylabel('Time-modulation (ZETA)');
%set(gca,'xscale','log','yscale','log')
title(sprintf('r(Loc,Time)=%.3f,p=%.3f',rLT,pLT));
fixfig;

%[rLM,pLM]=corr(log10(vecVisLocZetaP)',log10(vecMisMatZetaP'));
[rLM,pLM]=corr(vecVisLocZetaP',vecMisMatZetaP');
subplot(2,3,4)
scatter(vecVisLocZetaP,vecMisMatZetaP)
xlabel('Location-modulation (ZETA)');
ylabel('Mismatch-modulation (ZETA)');
%set(gca,'xscale','log','yscale','log')
title(sprintf('r(Loc,MisM)=%.3f,p=%.3f',rLM,pLM));
fixfig;

%[rTM,pTM]=corr(log10(vecVisTimZetaP)',log10(vecMisMatZetaP'));
[rTM,pTM]=corr(vecVisTimZetaP',vecMisMatZetaP');
subplot(2,3,5)
scatter(vecVisTimZetaP,vecMisMatZetaP)
xlabel('Time-modulation (ZETA)');
ylabel('Mismatch-modulation (ZETA)');
%set(gca,'xscale','log','yscale','log')
title(sprintf('r(Time,MisM)=%.3f,p=%.3f',rTM,pTM));
fixfig;
maxfig;
normaxes('xy');

%% calculate mismatch-modulated cells
vecMismatchTrials;
vecMismatchT;
vecMismatchLoc;
intMismatchNum;


%% compare mismatch and visual responsiveness