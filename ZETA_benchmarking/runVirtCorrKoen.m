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
%}

%% load data
strDataPath = 'F:\Data\Processed\VirtualTunnel\Delier\';
strDataFile = 'DelierPreProSpikes.mat';
load([strDataPath strDataFile]);


%% transform times
licktimes = sInfo.Stim.log.licklog(:,2); %in frame or stim time?

vecStartT = sInfo.StimTimes(:,1);
intTrials = numel(vecStartT);
vecLocalT = sInfo.Stim.log.stimlog(1:intTrials,1);
vecRewardT = sInfo.Stim.log.stimlog(1:intTrials,2) - vecLocalT + vecStartT;
vecRewardT(sInfo.Stim.log.stimlog(1:intTrials,2)==0) = [];
vecGreyT = sInfo.Stim.log.stimlog(1:intTrials,3) - vecLocalT + vecStartT;
vecGreyT(sInfo.Stim.log.stimlog(1:intTrials,3)==0) = [];
vecMismatchT = sInfo.Stim.log.stimlog(1:intTrials,4) - vecLocalT + vecStartT;
vecMismatchT(sInfo.Stim.log.stimlog(1:intTrials,4)==0) = [];

%% prep analysis
dblUseMaxDur = round(min(diff(vecStartT)));
intResampNum = 100;
intPlot = 0;
intLatencyPeaks = 4;
boolVerbose = false;

%% run analysis
cellAlignType = {'Start','Grey','Reward','Mismatch'};
intNeurons = numel(cellSpikeTimes);
matRespD = nan(4,intNeurons);
vecFilt = normpdf(-10:10,0,5)./sum(normpdf(-10:10,0,5));
parfor intNeuron=1:intNeurons
	intNeuron
	clf;
	for intAlignType=1:4
		if intAlignType == 1
			vecEventT = vecStartT;
			
		elseif intAlignType == 2
			vecEventT = vecGreyT;
			
		elseif intAlignType == 3
			vecEventT = vecRewardT;
			
		elseif intAlignType == 4
			vecEventT = vecMismatchT;
			
		end
		strAlignType = cellAlignType{intAlignType};
		
		%trial start
		[dblZETA,vecLatencies,sZETA,sMSD] = getZeta(cellSpikeTimes{intNeuron},vecEventT,dblUseMaxDur,intResampNum,0,intLatencyPeaks);
		matRespD(intAlignType,intNeuron) = max(imfilt(sMSD.vecRate',vecFilt));
		
		%plot
		if intPlot
		subplot(4,3,1+((intAlignType-1)*3))
		plotRaster(cellSpikeTimes{intNeuron},vecEventT,dblUseMaxDur,10000);
		xlabel('Time from event (s)');
		ylabel('Event #');
		title(sprintf('Raster; neuron %d, %s',intNeuron, strAlignType));
		fixfig;
		grid off;
		
		subplot(4,3,2+((intAlignType-1)*3))
		hold on
		line(repmat(sZETA.vecSpikeT,[1 50]),sZETA.matRandDiff(:,randperm(intResampNum,50)),'color',[0.5 0.5 0.5]);
		plot(sZETA.vecSpikeT,sZETA.vecD,'Color',lines(1));
		hold off
		xlabel('Time (s)');
		ylabel('Spiking anomaly (s)');
		title(sprintf('ZETA=%.3f (p=%.3f)',dblZETA,sZETA.dblP));
		fixfig
		
		subplot(4,3,3+((intAlignType-1)*3))
		hold on
		stairs(sMSD.vecT,sMSD.vecRate)
		vecLatencyVals = sZETA.vecLatencyVals;
		scatter(vecLatencies(1),vecLatencyVals(1),'bx');
		scatter(vecLatencies(2),vecLatencyVals(2),'b*');
		scatter(vecLatencies(3),vecLatencyVals(3),'gx');
		scatter(vecLatencies(4),vecLatencyVals(4),'rx');
		hold off
		xlabel('Time (s)');
		ylabel('Spiking rate (Hz)');
		title(sprintf('ZETA=%.0fms,-ZETA=%.0fms,Pk=%.0fms,On=%.2fms',vecLatencies(1)*1000,vecLatencies(2)*1000,vecLatencies(3)*1000,vecLatencies(4)*1000));
		fixfig
		end
	end
	if intPlot,pause;end
end
%%
subplot(2,2,1);
matCorr = corr(matRespD');
imagesc(matCorr,[-1 1]);colormap(redblue);colorbar

subplot(2,2,2)
scatter(matRespD(1,:),matRespD(4,:))