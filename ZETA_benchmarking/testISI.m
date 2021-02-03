%{
ZETA alternative
generate bootstrapped PSTHs by taking random samples from real inter-spike
intervals; compare with true PSTH
- ZETA takes into account entire history dependence, MIMI only accounts for
1-spike back
=> counter-example: bursty cell

IFR alternative
MIMI

%}





%% select all neurons in LP and drifting grating stimuli
clearvars;
strArea = 'primary visual';
[sAggStim,sAggNeuron]=loadDataNpx(strArea,'driftinggrating');

%% get data for lateral geniculate neuron #16 (original neuron #30)
intNeuronToAnalyze = 3;
sNeuron = sAggNeuron(intNeuronToAnalyze);
vecSpikeTimes = sNeuron.SpikeTimes;
vecStimOnTime = sAggStim(1).cellStim{1}.structEP.vecStimOnTime;
vecStimOffTime = sAggStim(1).cellStim{1}.structEP.vecStimOffTime;
dblTrialDur = min(vecStimOffTime-vecStimOnTime);
vecOrientation = sAggStim(1).cellStim{1}.structEP.Orientation;
boolJitter = false;
strRand = '';
strFigPath = 'F:\Data\Results\ZETA\Inclusion\';

%jitter
if boolJitter
vecStimOnTime = vecStimOnTime + 2*(rand(size(vecStimOnTime))-0.5)*dblTrialDur*2;
strRand = 'Jittered';
end
vecStimOnTime = sort(vecStimOnTime);

matEventTimes = cat(1,vecStimOnTime,vecStimOnTime+dblTrialDur);
return
%% fit
[dblMIMI_P,vecLatencies,sMIMI,sRate] = getMIMI(vecSpikeTimes,vecStimOnTime,dblTrialDur,4,2);

strFigFile = sprintf('TMZ_Example_%sN%d%s',strArea,intNeuronToAnalyze,strRand);
drawnow;
export_fig([strFigPath strFigFile '.tif']);
export_fig([strFigPath strFigFile '.pdf']);
