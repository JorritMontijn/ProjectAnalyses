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
strDataFile1 = 'Delier_20191017_002_Split1_normcorr_SPSIG.mat';
strDataFile2 = 'Delier_20191017_002_Split1_normcorr_SPSIG_Res.mat';
fprintf('Loading raw data from "%s" [%s]\n',strDataPath,getTime);
sLoad1 = load([strDataPath strDataFile1]);
sLoad2 = load([strDataPath strDataFile2]);

matCa = sLoad1.sigCorrected;

dblSamplingFreq = sLoad1.freq;
matSomaF = sLoad1.sig;
matNpF = sLoad1.sigBack;

%pre-alloc
[intFrames,intNeurons] = size(matSomaF);
matSoma_dFoF = nan(size(matSomaF));
matNp_dFoF = nan(size(matSomaF));
%% recalc dfof
fprintf('Recalculating dF/F0 for %d cells [%s]\n',intNeurons,getTime);
parfor intNeuron=1:intNeurons
	%soma
	matSoma_dFoF(:,intNeuron) = calcdFoF(matSomaF(:,intNeuron), dblSamplingFreq);
	
	%np
	matNp_dFoF(:,intNeuron) = calcdFoF(matNpF(:,intNeuron), dblSamplingFreq);
end

% calculate dFoF
dblNpFactor = 0.7;
mat_dFoF = matSoma_dFoF - dblNpFactor*matNp_dFoF;
%mat_dFoF = mat_dFoF - (mean(mat_dFoF,1) - min(mat_dFoF,[],1))/2;
mat_dFoF = mat_dFoF - sum(mat_dFoF.*(mat_dFoF<0),1)./intFrames;

%% detect spikes
fprintf('Running spike detection on %d cells [%s]\n',intNeurons,getTime)
%parameters
dblSpikeTau = 0.4;
dblThresholdFactor = 0.4;
intBlockSize = 4133;

%pre-allocate
cellFramesAP = cell(1,intNeurons);
cellNumberAP = cell(1,intNeurons);
matSpikes = nan(size(mat_dFoF));
matExpFit = nan(size(mat_dFoF));
cellSpikeTimes = cell(1,intNeurons);

%recalc dfof
vecT = (1:intFrames)/dblSamplingFreq;
dblTotDur = vecT(end);
for intNeuron=1:intNeurons
	%soma
	vec_dFoF = mat_dFoF(:,intNeuron)';
	[vecFramesAP, vecNumberAP, vecSpikes, vecExpFit, vecSpikeTimes] = doDetectSpikes(vec_dFoF,dblSamplingFreq,dblSpikeTau,intBlockSize,dblThresholdFactor);
	cellFramesAP{intNeuron} = vecFramesAP;
	cellNumberAP{intNeuron} = vecNumberAP;
	matSpikes(:,intNeuron) = vecSpikes;
	matExpFit(:,intNeuron) = vecExpFit;
	cellSpikeTimes{intNeuron} = vecSpikeTimes;
	
	%%{
	clf;
	hold on
	plot(vecT,vec_dFoF)
	plot(vecT,vecExpFit)
	ylabel('dF/F0');
	xlabel('Time (s)');
	hold off;
	title(sprintf('Neuron %d/%d; mean rate is %.1fHz',intNeuron,intNeurons,numel(vecSpikeTimes)/dblTotDur))
	fixfig;
	pause;
	%%}
end

%% save data
sInfo = sLoad2.info;
strFileOut = [strDataPath 'DelierPreProSpikes.mat'];
fprintf('Saving data to "%s" [%s]   ...   \n',strFileOut,getTime);
save(strFileOut,...
	'sInfo','cellSpikeTimes','mat_dFoF','matSoma_dFoF','matNp_dFoF','matExpFit','matSpikes','cellNumberAP','cellFramesAP',...
	'dblSamplingFreq','dblNpFactor','intBlockSize','dblThresholdFactor','dblSpikeTau')
fprintf('\bDone!\n')