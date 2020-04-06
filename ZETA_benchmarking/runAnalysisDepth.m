clear all;%close all;
%% recordings
%sites
cellRec{1}{1} = 'P:\Montijn\DataNeuropixels\Exp2019-11-20\20191120_MP2_RunDriftingGratingsR01_g0'; %115/298, SUBCORT (+some CORT), nice responses
cellRec{1}{2} = 'P:\Montijn\DataNeuropixels\Exp2019-11-21\20191121_MP2_RunDriftingGratingsR01_g0'; %72/285, nice responses
cellRec{1}{3} = 'P:\Montijn\DataNeuropixels\Exp2019-11-22\20191122_MP2_RunDriftingGratingsR01_g0'; %54/571, a few very nice responses, B2 only good up to T=2500
cellRec{1}{4} = 'P:\Montijn\DataNeuropixels\Exp2019-11-22\20191122_MP2_R02_RunDriftingGratingsR01_g0'; %check original log files; assignment problem?
cellRec{2}{1} = 'P:\Montijn\DataNeuropixels\Exp2019-12-10\20191210_MP3_RunDriftingGratingsR01_g0'; %CORT{PM} (+some subcort)53/283, some nice responses
cellRec{2}{2} = 'P:\Montijn\DataNeuropixels\Exp2019-12-11\20191211_MP3_RunDriftingGratingsR01_g0'; %CORT{V1}+SUBCORT{LP}, 182/417, many cells...
cellRec{2}{3} = 'P:\Montijn\DataNeuropixels\Exp2019-12-12\20191212_MP3_RunNaturalMovieR01_g0'; %CORT{}+SUBCORT{NOT?}, 120/388, B1 bad, B2 good, B3 trials 1050-1150 bad, B5 excellent
cellRec{2}{4} = 'P:\Montijn\DataNeuropixels\Exp2019-12-13\20191213_MP3_RunDriftingGratingsR01_g0'; %196/512, 98 CORT{}+ 145 SUBCORT{NOT?}, beautiful
cellRec{2}{5} = 'P:\Montijn\DataNeuropixels\Exp2019-12-16\20191216_MP3_RunNaturalMovieR01_g0'; %232/621, CORT{V1}+SUBCORT{LGN}, amazing...
cellRec{2}{6} = 'P:\Montijn\DataNeuropixels\Exp2019-12-17\20191217_MP3_RunDriftingGratingsR01_g0'; %72/407, CORT{}, few SUBCORT
cellRec{3}{1} = 'P:\Montijn\DataNeuropixels\Exp2020-01-15\20200115_MP4_RunDriftingGratingsR01_g0'; %133/398, SUBCORT (+some CORT) very nice responses
cellRec{3}{2} = 'P:\Montijn\DataNeuropixels\Exp2020-01-16\20200116_MP4_RunDriftingGratingsR01_g0'; %47/325, CORT (+some SUBCORT) quite poor
cellRec{3}{3} = 'P:\Montijn\DataNeuropixels\Exp2020-01-16\20200116_MP4_RunDriftingGratingsR02_g0'; %51/216, SUBCORT, but very nice cells

vecDepthCorrection{1}{1} = 0;%01
vecDepthCorrection{1}{2} = 0;%02
vecDepthCorrection{1}{3} = 0;%03
vecDepthCorrection{1}{4} = nan;%04
vecDepthCorrection{2}{1} = 0;%05
vecDepthCorrection{2}{2} = 0;%06
vecDepthCorrection{2}{3} = 0;%07
vecDepthCorrection{2}{4} = 0;%08
vecDepthCorrection{2}{5} = 0;%09
vecDepthCorrection{2}{6} = 0;%10
vecDepthCorrection{3}{1} = 0;%11
vecDepthCorrection{3}{2} = 0;%12
vecDepthCorrection{3}{3} = 0;%13

vecRunAnal = [3 2];

%% load data
cellRecPath = strsplit(cellRec{vecRunAnal(1)}{vecRunAnal(2)},filesep);
strProcPath = 'D:\Data\Processed\Neuropixels\';
strExp = cellRecPath{end-1};
strMouse = strcat('M',getFlankedBy(cellRecPath{end},'_M','_'));
strRecIdx = strcat('S',num2str(vecRunAnal(1)),'L',num2str(vecRunAnal(2))); %subject / location
strFileAP = strjoin({strExp,strMouse,strRecIdx,'AP.mat'},'_');
strFilePathAP = fullfile(strProcPath,strFileAP); %sAP.sCluster(i).Contamination
%strFileDG = strjoin({strExp,strMouse,'DG.mat'},'_');
%strFilePathDG = fullfile(strProcPath,strFileDG); %sDG.sB(j).vecZeta

%load
load(strFilePathAP);
%load(strFilePathDG);

%% output
strFigDir = ['D:\Data\Results\OriMetric\TuningCurves\' strExp strRecIdx filesep];
if ~exist(strFigDir,'dir'),mkdir(strFigDir);end

%% depth
matZetaP =cell2mat({sAP.sCluster(:).ZetaP}');
matMeanP =cell2mat({sAP.sCluster(:).MeanP}');
vecContam=cell2vec({sAP.sCluster(:).Contamination});
vecKSG = cell2vec({sAP.sCluster(:).KilosortGood});
vecDepth =cell2vec({sAP.sCluster(:).Depth});
intNeurons = numel(vecContam);
vecOnsets = nan(1,intNeurons);
vecZETA = nan(1,intNeurons);
vecPrefDeg = nan(1,intNeurons);
vecDprime = nan(1,intNeurons);
vecOPI = nan(1,intNeurons);
vecDSI = nan(1,intNeurons);
cellArea = cell(1,intNeurons);
matAggResp = nan(24,intNeurons);
hTic = tic;

%% include
%indIncludeZ = any(matMeanP(:,2) > 0.05 & matZetaP(:,2)<0.05,2)
indIncludeZ = 1;%any(matZetaP<(0.05/size(matZetaP,2)),2);
indIncludeC = (vecContam < 0.1) | vecKSG;
indInclude = indIncludeZ(:) & indIncludeC(:);

%%
intMaxBlock = numel(sAP.cellStim);
intResampNum = 50;
intPlot = 0;%3;
intLatencyPeaks =4;
vecRestrictRange = [];%[0.025 0.1];
boolVerbose = false;
hTic=tic;
for intBlock=1:intMaxBlock
	strStimType = getDataAP(sAP,'strFile','stimblock',intBlock);
	if contains(strStimType,'NaturalMovie') && intPlot == 0
		continue;
	end
			
	[dummy,vecReorder]=sort(matZetaP(:,intBlock),'ascend');
	vecStimOriDegrees = getDataAP(sAP,'Orientation','stimblock',intBlock);
	vecN = 1:intNeurons;
	vecN = vecN(vecReorder);
	vecN(~indInclude(vecReorder)) = [];
	vecN = vecN(:)';
	if intPlot > 0
		vecN = vecN(1:10);
	end
	for intNeuronIdx=1:numel(vecN)%1:2
		intNeuron = vecN(intNeuronIdx);
		if toc(hTic) > 5 || intPlot > 0
			fprintf('Now at block %d, neuron %d/%d (%d) [%s]\n',intBlock,intNeuronIdx,numel(vecN),intNeuron,getTime);
			hTic = tic;
		end
		[strArea] = getDataAP(sAP,intNeuron,'Area');
		[vecSpikeTimes,vecStimOn,vecStimOff,vecStimType,sStimObject] = getDataAP(sAP,intNeuron,'stimblock',intBlock);
		strStimType = sStimObject(1).StimType;
		
		dblUseMaxDur = round(median(diff(vecStimOn))*2)/2;
		[dummy,vecUnique,vecRepsPerType,cellSelect,vecRepetition] = label2idx(vecStimType);
		if intPlot > 0
			close;
		end
		%varEventTimes = [vecStimOn(1:1600);vecStimOff(1:1600)]';
		varEventTimes = [vecStimOn;vecStimOff]';
		[dblZETA,vecLatencies,sZETA,sRate] = getZeta(vecSpikeTimes,varEventTimes,dblUseMaxDur,intResampNum,intPlot,intLatencyPeaks,vecRestrictRange,boolVerbose);
		if intLatencyPeaks > 0
			vecOnsets(intNeuron) = vecLatencies(4);
		end
		vecZETA(intNeuron) = dblZETA;
		if intPlot,
			title(subplot(2,3,1),strFileAP,'interpreter','none');
			title(subplot(2,3,2),sprintf('B%d:%s, N%d; depth %d',intBlock,strStimType,intNeuron,round(vecDepth(intNeuron))),'interpreter','none');
			title(subplot(2,3,3),sprintf('Total cells: %d, include: %d',numel(indInclude),sum(indInclude)),'interpreter','none');
			pause;
		end
		
		%get tuning curves
		[vecStimResp,vecBaseResp] = getTrialResp(vecSpikeTimes,vecStimOn,vecStimOff);
		[sOut] = getTuningCurves(vecStimResp,vecStimOriDegrees);
		intOriNum = numel(sOut.vecUniqueDegs);
		intReps = numel(vecStimOriDegrees) / intOriNum;
		%errorbar(vecPlotDegs,sOut.matMeanResp,sOut.matSDResp./sqrt(intReps));
		%hold on;
		%scatter(vecStimOriDegrees+5*rand(size(vecStimOriDegrees)),vecStimResp+0.3*rand(size(vecStimResp)))
		%plot(vecPlotDegs,sOut.matFittedResp);
		%hold off;
		x=feval(sOut.funcFit,sOut.matFittedParams,[sOut.matFittedParams(:,1) sOut.matFittedParams(:,1)+pi]);
		dblDSI = 1-(min(x)/max(x));
		if sOut.matFittedParams(:,3) > 1
			dblPrefDeg = rad2deg(mod(sOut.matFittedParams(:,1)+pi,2*pi));
		else
			dblPrefDeg = rad2deg(sOut.matFittedParams(:,1));
		end
		dblDeltaPrime = getDeltaPrime(vecStimResp,deg2rad(vecStimOriDegrees));
		dblOPI = getOPI(vecStimResp,deg2rad(vecStimOriDegrees));
		matAggResp(:,intNeuron) = sOut.matMeanResp;
		vecDSI(intNeuron) = dblDSI;
		vecPrefDeg(intNeuron) = dblPrefDeg;
		vecDprime(intNeuron) = dblDeltaPrime;
		vecOPI(intNeuron) = dblOPI;
	end
	if ~contains(strStimType,'NaturalMovie') && intPlot == 0
		break;
	end
end

if 0
	%%
	export_fig([strFigDir sprintf('%sB%02d_N%03d_%s.tif',strRecIdx,intBlock,intNeuron,strStimType)]);
	export_fig([strFigDir sprintf('%sB%02d_N%03d_%s.pdf',strRecIdx,intBlock,intNeuron,strStimType)]);
end
%%
indInclude2 = indInclude;% & (vecZETA(:) > 2);
vecPlotDegs = mod(sOut.vecUniqueDegs+180,360)
cellArea = {sAP.sCluster(:).Area};
indNeuronsNOT = contains(cellArea(:),'optic tract') | contains(cellArea(:),'Lateral posterior nucleus');
vecCorrectedDepth = -(vecDepth + vecDepthCorrection{vecRunAnal(1)}{vecRunAnal(2)});
%close all
figure
subplot(2,3,1)
%scatter3(vecDepth(indInclude),abs(vecZETA(indInclude)),vecContam(indInclude));
scatter(vecOnsets(indInclude2&indNeuronsNOT),vecCorrectedDepth(indInclude2&indNeuronsNOT),[],abs(vecZETA(indInclude2&indNeuronsNOT)),'x');
hold on
scatter(vecOnsets(indInclude2&~indNeuronsNOT),vecCorrectedDepth(indInclude2&~indNeuronsNOT),[],abs(vecZETA(indInclude2&~indNeuronsNOT)),'o');
hold off

ylim([-3000 0]);
xlim([0 0.1]);
xlabel('Onset latency (s)');
ylabel('Depth from pia (\mum)');
h=colorbar;
ylabel(h,'ZETA');
fixfig;
title(sprintf('Cortex: %d; subcortex: %d',sum(vecCorrectedDepth(indInclude2) > -1000),sum(vecCorrectedDepth(indInclude2) < -1500)))

h=subplot(2,3,2)
scatter(abs(vecZETA(indInclude2&indNeuronsNOT)),vecCorrectedDepth(indInclude2&indNeuronsNOT),[],vecPrefDeg(indInclude2&indNeuronsNOT),'x');
hold on
scatter(abs(vecZETA(indInclude2&~indNeuronsNOT)),vecCorrectedDepth(indInclude2&~indNeuronsNOT),[],vecPrefDeg(indInclude2&~indNeuronsNOT),'o');
hold off
colormap(h,hsv);
ylim([-3000 0]);
%xlim([0 360]);
xlabel('ZETA');
ylabel('Depth from pia (\mum)');
h=colorbar;
ylabel(h,'pref-deg');
fixfig;
% subplot(2,2,2)
% vecBinEdges = -3000:100:0;
% [vecCounts,vecMeans,vecSDs,cellVals,cellIDs] = makeBins(vecCorrectedDepth(indInclude),vecOnsets(indInclude),vecBinEdges);
% vecBins = vecBinEdges(2:end)-diff(vecBinEdges);
% herrorbar(vecMeans,vecBins,vecSDs,vecSDs);
% xlim([0 0.1]);
% xlabel('Onset latency (s)');
% ylabel('Depth from pia (\mum)');
% fixfig;

subplot(2,3,3)
%vecPrefDeg = mod(vecPrefDeg+180,360);
scatter(vecPrefDeg(indInclude2&indNeuronsNOT),vecCorrectedDepth(indInclude2&indNeuronsNOT),[],vecDSI(indInclude2&indNeuronsNOT),'x');
hold on
scatter(vecPrefDeg(indInclude2&~indNeuronsNOT),vecCorrectedDepth(indInclude2&~indNeuronsNOT),[],vecDSI(indInclude2&~indNeuronsNOT),'o');
hold off
ylim([-3000 0]);
xlim([0 360]);
xlabel('Pref ori (deg)');
ylabel('Depth from pia (\mum)');
h=colorbar;
ylabel(h,'DSI');
fixfig;

%%
cellSpikes = {sAP.sCluster(:).SpikeTimes};
matStimResp = getTrialResp(cellSpikes,vecStimOn,vecStimOff);

vecRpO = nanmean(matAggResp(:,indNeuronsNOT & indInclude2),2);
		
subplot(2,3,4)
plot(vecPlotDegs,vecRpO)
xlim([0 360]);
xlabel('Stim ori (deg)');
ylabel('Mean rate');
title(sprintf('Resp LP'));
h=colorbar;
fixfig;


dblTopNOT = -1000;
dblBotNOT = -1800;
vecRpO = nanmean(matAggResp(:,vecCorrectedDepth < dblTopNOT & vecCorrectedDepth > dblBotNOT),2);
		
subplot(2,3,5)
plot(vecPlotDegs,vecRpO)
xlim([0 360]);
xlabel('Stim ori (deg)');
ylabel('Mean rate');
title(sprintf('Resp %d - %d micron',abs(dblTopNOT),abs(dblBotNOT)));
h=colorbar;
fixfig;


dblTopNOT = -3000;
dblBotNOT = -3500;
vecRpO = nanmean(matAggResp(:,vecCorrectedDepth < dblTopNOT & vecCorrectedDepth > dblBotNOT),2);
		
subplot(2,3,6)
plot(vecPlotDegs,vecRpO)
xlim([0 360]);
xlabel('Stim ori (deg)');
ylabel('Mean rate');
title(sprintf('Resp %d - %d micron',abs(dblTopNOT),abs(dblBotNOT)));
h=colorbar;
fixfig;
