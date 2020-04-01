%% analyze input strength dependency for dimensionality of pop responses

%% initialize
%settings
boolLoadSync = false;
boolLoadLGN = false;
boolLoadSpikeTimes = true;
boolLoadTrialSplit = false;
boolDoSplitAnalysis = false;

intLoadSim = 100;
clearvars -except intSubSample* vecDoShuff vecRunSims intLoadSim bool* vecRunAreas intSamples
boolLoad = true;

%% get simulation name [strSimulation] from [intLoadSim]
loadSim;

%% RUN: #header
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
if boolDoSplitAnalysis
	strTag = ['TS_' strTag];
	strType = ['TS_' strType];
end

%% plot
intNeurons = numel(cellSpikeTimesCortex);
intTrial = 4;
dblBinStep = 0.001;
dblStartT = vecStimStartSecs(intTrial);
dblStopT = vecStimStopSecs(intTrial);
vecPlotWindow = [dblStartT - 0.2 dblStopT + 0.2];
vecPlotWindow = round(vecPlotWindow / dblBinStep) * dblBinStep;
vecBins = vecPlotWindow(1):dblBinStep:vecPlotWindow(end);
vecPlotBins = vecBins(2:end) - dblBinStep/2;

%get high-res spikes
matSpikes = zeros(intNeurons,length(vecPlotBins));
for intN=1:intNeurons
	vecCounts = histcounts(cellSpikeTimesCortex{intN},vecBins);
	matSpikes(intN,:) = vecCounts;
end
matSpikes = matSpikes>0;

%get low-res spikes
dblMeanBinStep = 0.025;
vecMeanBins = vecPlotWindow(1):dblMeanBinStep:vecPlotWindow(end);
vecMeanPlotBins = vecMeanBins(2:end) - dblMeanBinStep/2;

matMeanSpikes = zeros(intNeurons,length(vecMeanPlotBins));
for intN=1:intNeurons
	vecCounts = histcounts(cellSpikeTimesCortex{intN},vecMeanBins);
	matMeanSpikes(intN,:) = vecCounts;
end

%plot
figure
subplot(2,1,1)
errorbar(vecMeanPlotBins,mean(matMeanSpikes,1)/dblMeanBinStep,(std(matMeanSpikes,[],1)/dblMeanBinStep)/sqrt(intNeurons))
subplot(2,1,2)
implot(matSpikes);
vecTicks = get(gca,'xtick');
set(gca,'xticklabel',vecBins(vecTicks+1))
hold on
plot([find(vecBins>dblStartT,1,'first') find(vecBins>dblStopT,1,'first')],intNeurons*[1 1]/2,'b')
hold off


%% make example of picture to LGN
matRetLGN = eye(32*32,32*32);
matPlusFilt = normpdf(-10:10,0,6);
matPlusFilt = matPlusFilt ./ sum(matPlusFilt(:));
matMinFilt = normpdf(-50:50,0,15);
matMinFilt = matMinFilt ./ sum(matMinFilt(:));
matRetLGN = conv2(matRetLGN,matPlusFilt,'same') - conv2(matRetLGN,matMinFilt,'same');

dblMax = max(abs(matRetLGN(:)));
figure
imagesc(matRetLGN,dblMax*[-1 1]);colormap(redblue)

