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
	boolLoadSpikeTimes = true;
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
for intNeuron=638:size(matData,1)
[sOutTemp] = getTuningCurves(matData(intNeuron,:)./median(vecStimStopSecs-vecStimStartSecs),vecTrialOris,true);
fixfig
ylabel('Spiking rate (Hz)');
xlabel('Stimulus orientation (degs)');
set(gca,'xtick',0:90:360,'xticklabel',0:45:180)
title([strTag '; Neuron ' num2str(intNeuron)],'interpreter','none')
pause%(0.1)
end
return
%%
export_fig(['F:\Data\Results\SimFigs\Neuron' num2str(intNeuron) '_' strTag '.tif']);
export_fig(['F:\Data\Results\SimFigs\Neuron' num2str(intNeuron) '_' strTag '.pdf']);
%638
%%
[sOut] = getTuningCurves(matData./median(vecStimStopSecs-vecStimStartSecs),vecTrialOris,false);
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
%%
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


%% bandwidth
vecGoodNeurons = sOut.vecFitR2 > 0.1;
histx(rad2deg(real(sOut.matBandwidth(vecGoodNeurons)))/4); %/4 because the output doubles the angle to create full circle and gives FWHM not HWHM
xlabel('Bandwidth (HWHM) in degs')
ylabel('Neuron count')
fixfig

return
%%
export_fig(['F:\Data\Results\SimFigs\Bandwidth_' strTag '.tif']);
export_fig(['F:\Data\Results\SimFigs\Bandwidth_' strTag '.pdf']);

%% outside stim
hTic = tic;
matModelRespITI = nan(size(matModelResp));
vecBinsTimeITI = vecBinsTime;
vecBinsTimeITI(1:2:end) = vecBinsTimeITI(1:2:end)+0.05;
dblITI_Dur = mean(vecBinsTimeITI(4:2:end) - vecBinsTimeITI(3:2:end));

vecBinsTimeStim = vecBinsTime;
vecBinsTimeStim(2:2:end) = vecBinsTimeStim(2:2:end)+0.05;
dblStim_Dur = mean(vecBinsTimeStim(3:2:end) - vecBinsTimeStim(2:2:(end-1)));
matModelRespStim = nan(size(matModelResp));
for intNeuron=1:intNeurons
		[vecCounts,edges,bin] = histcounts(cellSpikeTimesCortex{intNeuron},vecBinsTimeITI);
		matModelRespITI(intNeuron,:) = uint16(vecCounts(3:2:end))/dblITI_Dur; %Counts, not Hz!
		[vecCounts,edges,bin] = histcounts(cellSpikeTimesCortex{intNeuron},vecBinsTimeITI);
		matModelRespStim(intNeuron,:) = uint16(vecCounts(2:2:end))/dblStim_Dur; %Counts, not Hz!
		
		if toc(hTic) > 5 || intNeuron==1
			hTic = tic;
			fprintf('Now at neuron %d/%d [%s]\n',intNeuron,intNeurons,getTime);
		end
	end
	%%
	[sOutITI] = getTuningCurves(matModelRespITI,vecTrialOris,false);

	subplot(2,2,1)
	histx(max(sOutITI.matMeanResp,[],2))
	
	subplot(2,2,2)
	histx(min(sOutITI.matMeanResp,[],2))
	
	[sOutStim] = getTuningCurves(matModelRespStim,vecTrialOris,false);

	subplot(2,2,3)
	histx(max(sOutStim.matMeanResp,[],2))
	
	subplot(2,2,4)
	histx(min(sOutStim.matMeanResp,[],2))
	
	%% overall FR
	(sum(cellfun(@numel,cellSpikeTimesCortex))/intNeurons)/vecStimStopSecs(end)