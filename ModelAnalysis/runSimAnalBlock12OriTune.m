%% analyze input strength dependency for dimensionality of pop responses

%% initialize
intType = 6;
%close all
clearvars -except intType
vecRunAreas = [1];

%settings
boolLoadSync = false;
boolLoadLGN = false;
boolLoadSpikeTimes = false;
boolLoadTrialSplit = false;
boolDoSplitAnalysis = false;
boolUseZscore = true;
intSamples = 100;

%save figs and data?
boolSaveFigs = true;
boolSaveData = true;

intLoadSim=999;
boolLoad = true;

%% get simulation name [strSimulation] from [intLoadSim]
loadSim;

%% RUN: #header
strBlockNr = getFlankedBy(mfilename,'Block','');
strBlockNr = strBlockNr(1);
strFigDir = ['F:\Data\Results\SimFigs\Block' strBlockNr '\'];
strDataDir = ['F:\Data\Results\SimFigs\Data' strBlockNr '\'];
if isempty(strBlockNr),error;end

if boolLoad
	runAnalysisHeader;
end
intRepetitions = min(sum(bsxfun(@eq,vecTrialStimType,vecUniqueStimTypes'),2));
intNeurons = size(matData,1);
strType;

%save(['DataForRex_' strSimulation '.mat'],'matModelResp','vecCellArea','vecTrialStimType','vecStimTypeOris','vecStimTypeOriNoise');
%return

cellIn = strsplit(strSimulation,'_');
strFiller = cellIn{1};
strType = cell2mat(cellIn(2:(end-1)));
strDate = cellIn{end};
strTag = [strType '_' strDate];

%% start
vecUseStimTypes = unique(matCompareTypes);
return
%% dependence of dimensionality of stimulus intensity
for intNeuron=638:size(matData,2)
[sOut] = getTuningCurves(matData(intNeuron,:),vecTrialOris,true);
pause
end
export_fig(['F:\Data\Results\SimFigs\Neuron' num2str(intNeuron) '_' strTag '.tif']);
export_fig(['F:\Data\Results\SimFigs\Neuron' num2str(intNeuron) '_' strTag '.pdf']);
%638
%%
[sOut] = getTuningCurves(matData,vecTrialOris,false);
%%
dblStep = 0.1;
vecBins = 0:dblStep:1;
vecBinC = (dblStep/2):dblStep:(1-dblStep/2);
vecCounts = histcounts(sOut.vecFitR2,vecBins)
figure
bar(vecBinC,vecCounts)
xlabel('R^2 von Mises fit');
ylabel('# of neurons');
drawnow;pause(0.1);
export_fig(['F:\Data\Results\SimFigs\R2_fits' strTag '.tif']);
export_fig(['F:\Data\Results\SimFigs\R2_fits' strTag '.pdf']);

%%
figure
scatter(vecPrefSF,sOut.vecFitR2)
xlabel('Preferred spatial frequency')
ylabel('R^2 tuning curve')
drawnow;pause(0.1);
export_fig(['F:\Data\Results\SimFigs\SF_R2_' strTag '.tif']);
export_fig(['F:\Data\Results\SimFigs\SF_R2_' strTag '.pdf']);

%%
scatter(sOut.vecFitR2,mean(sOut.matFittedParams(:,sOut.indKappa),2))


