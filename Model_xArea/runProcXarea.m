

%% get simulation name [strSimulation] from [intLoadSim]
clear all;close all;
strSimulation = 'Simulation_xAreaDistributed_SG18_2019-07-04';
strSimDataPath = 'D:\Data\SimAggregates\';
load([strSimDataPath strSimulation '.mat']);
strFigPath = 'D:\Data\ResultsModelXarea\';

%% RUN: #header
cellIn = strsplit(strSimulation,'_');
strFiller = cellIn{1};
strType = cell2mat(cellIn(2:(end-1)));
strDate = cellIn{end};
strTag = [strType '_' strDate];

%% run analysis
boolMakePlots = true;
vecFitOriDegs = 0:359;
cellArea = {'V1','V2'};
intNeurons = numel(cellSpikeTimesCortex);
dblStepT = mean(diff(vecOverallT));
for intNeuron=60:intNeurons
	%% rename variables
	strArea = cellArea{vecCellArea(intNeuron)};
	vecSpikeTimes = dblStepT*double(cellSpikeTimesCortex{intNeuron}');
	vecStimOnTime = vecStimStartSecs;
	vecStimOffTime = vecStimStopSecs;
	vecStimOriDegrees = vecTrialOris;
	vecStimOriRads = deg2rad(vecStimOriDegrees);
	
	%% get trial data
	vecStimCounts = getSpikeCounts(vecSpikeTimes,vecStimOnTime,vecStimOffTime);
	vecStimResp = bsxfun(@rdivide,vecStimCounts,(vecStimOffTime-vecStimOnTime)); %transform to Hz
	dblBaseDur = 0.5;
	vecBaseCounts = getSpikeCounts(vecSpikeTimes,vecStimOnTime-dblBaseDur,vecStimOnTime);
	vecBaseResp = bsxfun(@rdivide,vecBaseCounts,dblBaseDur); %transform to Hz
	
	%% get tuning curves & parameters
	dblRho = getTuningRho(vecStimResp,vecStimOriRads);
	dblDeltaPrime = getDeltaPrime(vecStimResp,vecStimOriRads);
	sTuning = getTuningCurves(vecStimResp,vecStimOriDegrees);
	[matRespNSR,vecStimTypes,vecUniqueDegs] = getStimulusResponses(vecStimResp,vecStimOriDegrees);
		
	%% save data
	cellArea{intNeuron} = strArea;
	vecParams = sTuning.matFittedParams;
	matFitParams(:,intNeuron) = vecParams;
	vecPrefDir(intNeuron) = vecParams(1);
	matFitResp(:,intNeuron) = feval(sTuning.funcFit,vecParams,deg2rad(vecFitOriDegs));
	%vecParams(1) = pi/2;
	%rho
	vecRho(intNeuron) = dblRho;
	%d'
	vecDeltaPrime(intNeuron) = dblDeltaPrime;
	%raw ori, mean per ori, and sd per ori
	cellRawOri{intNeuron} = sTuning.vecUniqueRads;
	cellRawMean{intNeuron} = sTuning.matMeanResp;
	cellRawSD{intNeuron} = sTuning.matSDResp;
	%peak CV + BW
	vecPeaksCV(:,intNeuron) = sTuning.matVariance;
	vecPeaksBW(:,intNeuron) = sTuning.matBandwidth;
	
	%%% plot raw tuning curves
	%vecPlotTheta = [cellRawOri{intNeuron} cellRawOri{intNeuron}(1)];
	%vecPlotRho = cellRawMean{intNeuron}/max(cellRawMean{intNeuron});
	%vecPlotRho(end+1) = vecPlotRho(1);
	%intPlot = find(ismember(cellUniqueAreas,cellArea{intNeuron}));
	%axes(vecH(intPlot));
	%polarplot(vecPlotTheta,vecPlotRho);
	
	
	%% get cluster quality
	boolMakePlotsCQ = true;
	close;close;
	sClustQual = getClusterQuality(vecSpikeTimes,vecStimOnTime,boolMakePlotsCQ);
	vecNonStatIdx(intNeuron) = sClustQual.dblNonstationarityIndex;
	vecViolIdx(intNeuron) = sClustQual.dblViolIdx1ms;
	
	%% plot neuron
	if boolMakePlotsCQ
	%build figure name
	strFileName = sprintf('%s%s_N%04d',strArea,strTag,intNeuron);
	
	%add orientation tuning curve
	subplot(2,2,4)
	vecMean = nanmean(matRespNSR,3);
	vecSD = nanstd(matRespNSR,[],3);
	[dblMaxResp,intMaxIdx] = max(vecMean);
	vecShiftDegs = circshift(1:numel(vecUniqueDegs),[0 (round(numel(vecUniqueDegs)/2))-intMaxIdx]);
	errorbar(vecUniqueDegs-180,vecMean(vecShiftDegs),vecSD(vecShiftDegs))
	xlim([-180 180]);
	set(gca,'xtick',-180:90:180);
	ylim([0 max(get(gca,'ylim'))]);
	xlabel('Angle from strongest response (deg)');
	ylabel('Firing rate (Hz)');
	title(sprintf('%s %s, N%d, %s''=%.3f, %s=%.3f',strArea,strTag,intNeuron,getGreek('delta','lower'),dblDeltaPrime,getGreek('rho','lower'),dblRho),'interpreter','none');
	fixfig
	
	%save plot
	drawnow;
	export_fig([strFigPath strFileName 'ClustQual.tif']);
	export_fig([strFigPath strFileName 'ClustQual.pdf']);
	end
	
	%% get visual responsiveness
	matStimOnOff = cat(2,vecStimOnTime',vecStimOffTime');
	[dblZ,vecInterpT,vecZ,matDiffTest,dblHzD,dblP] = getVisualResponsiveness(vecSpikeTimes,matStimOnOff,boolMakePlots);
	hAx=get(gcf,'Children');
	axes(hAx(end));
	title(sprintf('%s %s, N%d, d=%.3f, d(Hz)=%.3f',strArea,strTag,intNeuron,dblZ,dblHzD),'interpreter','none');
	% assign data
	vecVisZ(intNeuron) = dblZ;
	vecVisHzD(intNeuron) = dblHzD;
	
	%% save plot
	if boolMakePlots
	drawnow;
	export_fig([strFigPath strFileName 'VisResp.tif']);
	export_fig([strFigPath strFileName 'VisResp.pdf']);
	
	pause(0.1);
	close
	close
	end
	
end