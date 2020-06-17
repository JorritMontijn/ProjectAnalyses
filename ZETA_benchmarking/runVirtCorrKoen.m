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
vecMismatchT = sInfo.Stim.log.stimlog(1:intTrials,4);
vecMismatchT(sInfo.Stim.log.stimlog(1:intTrials,4)==0) = [];
vecMismatchLoc = sInfo.Stim.log.stimlog(1:intTrials,5);
vecMismatchLoc(sInfo.Stim.log.stimlog(1:intTrials,5)==0) = [];
intMismatchNum = numel(vecMismatchLoc);

%% select control locations for mismatch at same position
cellControlT = cell(1,intMismatchNum);
vecAllLocs = sInfo.Stim.log.motionlog(:,3);
vecAllT = sInfo.Stim.log.motionlog(:,2);
for intMM=1:intMismatchNum
	%get potential locations
	dblMismatchLoc = vecMismatchLoc(intMM);
	vecPotIdx = find(diff(vecAllLocs > dblMismatchLoc)==1);
	vecPotLocFloor = vecAllLocs(vecPotIdx);
	vecPotLocCeil = vecAllLocs(vecPotIdx+1);
	matLoc = [vecAllLocs(vecPotIdx) vecAllLocs(vecPotIdx+1)];
	matT = [vecAllT(vecPotIdx) vecAllT(vecPotIdx+1)];
	vecPotT = nan(size(vecPotIdx));
	for intPot=1:numel(vecPotIdx)
		vecPotT(intPot) = interp1(matLoc(intPot,:),matT(intPot,:),dblMismatchLoc);
	end
	vecPotTrials = sum(bsxfun(@gt,vecPotT,vecLocalT'),2);
	
	%remove trials
	[a,b,c]=unique(vecPotTrials);
	indFirstOccurrence = false(size(vecPotTrials));
	indFirstOccurrence(b) = true;
	cellControlT{intMM} = vecPotT(indFirstOccurrence & vecPotTrials > 0 & vecPotTrials <= intTrials & ~ismember(vecPotTrials,vecMismatchTrials));
end

%% prep analysis
dblUseMaxDur = 3;%round(min(diff(vecStartT)));
intResampNum = 100;
intPlot = 0;
intLatencyPeaks = 4;
boolVerbose = false;

%% align time
cellControlAlignT = cell(size(cellControlT));
vecAlignMismatchT = nan(size(vecMismatchT));
for intMM=1:intMismatchNum
	%control
	cellControlAlignT{intMM} = zeros(size(cellControlT));
	for intEntry=1:numel(cellControlT{intMM})
		dblOldT = cellControlT{intMM}(intEntry);
		intTargetTrial = find(dblOldT>vecLocalT,1,'last');
		dblOffsetT = -vecLocalT(intTargetTrial) + vecStartT(intTargetTrial);
		cellControlAlignT{intMM}(intEntry) = dblOldT + dblOffsetT;
	end
	
	%mismatch
	dblOldT = vecMismatchT(intMM);
	intTargetTrial = find(dblOldT>vecLocalT,1,'last');
	dblOffsetT = -vecLocalT(intTargetTrial) + vecStartT(intTargetTrial);
	
	vecAlignMismatchT(intMM) = dblOldT + dblOffsetT;
end

%% run analysis
vecBinSizes = [(2.^[-4:-1])./dblSamplingFreq (1/dblSamplingFreq)*[1:floor(dblSamplingFreq)]];
intBinNum = numel(vecBinSizes);
cellAlignType = {'Start','Grey','Reward','Mismatch'};
intNeurons = numel(cellSpikeTimes);

intIters = 3;
matRespD = nan(2,intNeurons,intIters);
matRespD2 = nan(2,intNeurons,intIters);
matRespBins = nan(2,intNeurons,intBinNum,intIters);
intPlot=4
%figure;
hTic = tic;
for intIter=1:intIters
	for intNeuron=6%1:intNeurons
		if toc(hTic) > 5
			hTic = tic;
			fprintf('Iter %d/%d, neuron %d/%d [%s]\n',intIter,intIters,intNeuron,intNeurons,getTime);
		end
	%%
	clf;
	sOptions = struct;
	sOptions.handleFig = -1;
	
	for intAlignType=1:2
		if intAlignType == 1
			% control
			vecEventT = sort(cellfun(@(x) x(randi(numel(x))),cellControlAlignT));
			strAlignType = 'Control';
		elseif intAlignType == 2
			% mismatch
			vecEventT = sort(vecAlignMismatchT);
			strAlignType = 'Mismatch';
		end
		%zeta
		[dblZetaP,vecLatencies,sZETA,sMSD] = getZeta(cellSpikeTimes{intNeuron},vecEventT,dblUseMaxDur,intResampNum,0,intLatencyPeaks);
		matRespD(intAlignType,intNeuron,intIter) = abs(sZETA.dblZETA);
		if isempty(sMSD) || isempty(sMSD.vecRate)
			dblRate = 0;
		else
			dblRate = max(sMSD.vecRate);
		end
		matRespD2(intAlignType,intNeuron,intIter) = dblRate;
		
		%bins
		%matRespD(intAlignType,intNeuron) = max(imfilt(sMSD.vecRate',vecFilt));
		for intBinIdx=1:intBinNum
			dblBinSize = vecBinSizes(intBinIdx);
			[vecMean,vecSEM,vecWindowBinCenters] = doPEP(cellSpikeTimes{intNeuron},0:dblBinSize:(dblUseMaxDur+dblBinSize/2),vecEventT,sOptions);
			matRespBins(intAlignType,intNeuron,intBinIdx,intIter) = max(vecMean);
			
			if dblBinSize==(1/dblSamplingFreq) && intPlot
				subplot(4,2,3+(intAlignType-1))
				errorbar(vecWindowBinCenters,vecMean,vecSEM)
				xlabel('Time from event (s)');
				ylabel('Mean spiking rate (HZ)');
				title(sprintf('PETH; neuron %d, %s',intNeuron, strAlignType));
				fixfig;
				grid off;

			end
		end
		
		%plot
		if intPlot
			%%
			subplot(4,2,1+(intAlignType-1))
			plotRaster(cellSpikeTimes{intNeuron},vecEventT,dblUseMaxDur,10000);
			xlabel('Time from event (s)');
			ylabel('Event #');
			title(sprintf('Raster; neuron %d, %s',intNeuron, strAlignType));
			fixfig;
			grid off;
			
			subplot(4,2,5+(intAlignType-1))
			hold on
			line(repmat(sZETA.vecSpikeT,[1 50]),sZETA.matRandD(:,randperm(intResampNum,50)),'color',[0.5 0.5 0.5]);
			plot(sZETA.vecSpikeT,sZETA.vecD,'Color',lines(1));
			hold off
			xlabel('Time (s)');
			ylabel('Spiking anomaly (s)');
			title(sprintf('ZETA=%.3f (p=%.3f)',dblZetaP,sZETA.dblP));
			fixfig
			
			subplot(4,2,7+(intAlignType-1))
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
	
	if intPlot,
		%normalize y axes
	normaxes(gcf,'y',[1 5]);
	normaxes(gcf,'y',[2 6]);
	normaxes(gcf,'y',[3 7]);
	normaxes(gcf,'y',[4 8]);
		pause;
	end
end
if 0
	%% save
	strTargetDir = 'F:\Data\Results\ZETA\';
	strFigName = sprintf('ExampleNeuronv2VirtCorr_N%d',intNeuron);
	fprintf('Saving figures [%s] ... \n',getTime)
	export_fig([strTargetDir strFigName '.tif']);
	export_fig([strTargetDir strFigName '.pdf']);
	fprintf('\bDone! [%s]\n',getTime);
end
end
%% get effect sizes per iters
vecD = nan(intIter,1);
vecD2 = nan(intIter,1);
matD_Bins = nan(intBinNum,intIter);
intSingleFrameBin = find(vecBinSizes==(1/dblSamplingFreq));
for intIter=1:intIters
	[h2,p2]=ttest(matRespD(1,:,intIter),matRespD(2,:,intIter));
	vecD(intIter) = getCohensD(matRespD(2,:,intIter),matRespD(1,:,intIter));
	vecD2(intIter) = getCohensD(matRespD2(2,:,intIter),matRespD2(1,:,intIter));
	
	
	for intBinIdx=1:intBinNum
		matD_Bins(intBinIdx,intIter) = getCohensD(matRespBins(2,:,intBinIdx,intIter),matRespBins(1,:,intBinIdx,intIter));
	end
end
matAvRespBins = nanmean(matRespBins(:,:,:,1),4);
matAvRespD = nanmean(matRespD(:,:,1),3);
matAvRespD2 = nanmean(matRespD2(:,:,1),3);

%% plot
figure
subplot(2,3,1)
hold on
scatter(matAvRespBins(1,:,intSingleFrameBin),matAvRespBins(2,:,intSingleFrameBin))
vecLim = [0 max([max(get(gca,'xlim')) max(get(gca,'ylim'))])];
plot(vecLim,vecLim,'k--');
hold off
xlim(vecLim);
ylim(vecLim);
title(sprintf('%.1fms bins (1 frame); D=%.3f',vecBinSizes(intSingleFrameBin)*1000,mean(matD_Bins(intSingleFrameBin,:))));
ylabel('Peak rate mismatches (Hz)');
xlabel('Peak rate controls (Hz)');
fixfig;

subplot(2,3,2)
hold on
scatter(matAvRespD(1,:),matAvRespD(2,:))
vecLim = [0 max([max(get(gca,'xlim')) max(get(gca,'ylim'))])];
plot(vecLim,vecLim,'k--');
hold off
xlim(vecLim);
ylim(vecLim);
title(sprintf('ZETA; D=%.3f',mean(vecD)));
ylabel('ZETA (\sigma Hz) mismatches');
xlabel('ZETA (\sigma Hz) controls');
fixfig;

subplot(2,3,3)
hold on
scatter(matAvRespD2(1,:),matAvRespD2(2,:))
vecLim = [0 max([max(get(gca,'xlim')) max(get(gca,'ylim'))])];
plot(vecLim,vecLim,'k--');
hold off
xlim(vecLim);
ylim(vecLim);
title(sprintf('ZETA MSD instantaneous rate; D=%.3f',mean(vecD2)));
ylabel('Peak rate mismatches (Hz)');
xlabel('Peak rate controls (Hz)');
fixfig;

%t-test on best bin
[dummy,intBest] = max(mean(matD_Bins,2));
[h,p]=ttest(matD_Bins(intBest,:),vecD');

subplot(2,2,3)
hold on
errorbar(vecBinSizes,mean(matD_Bins,2),std(matD_Bins,[],2)./sqrt(intIters),'kx-');
errorbar([min(vecBinSizes) max(vecBinSizes)],mean(vecD)*[1 1],(std(vecD)./sqrt(intIters))*[1 1],'color',[0.1 0.1 0.8]);
errorbar([min(vecBinSizes) max(vecBinSizes)],mean(vecD2)*[1 1],(std(vecD2)./sqrt(intIters))*[1 1],'color',[0.8 0.1 0.1]);
scatter(vecBinSizes(intSingleFrameBin),mean(matD_Bins(intSingleFrameBin,:)),[],[0.1 0.8 0.1]);
hold off
ylabel('Mismatch/normal effect size (Cohen''s d)');
xlabel('Bin size (s)');
title(sprintf('Blue: ZETA, red: MSD, black: binning, green: single frame; best bin vs zeta,p=%.3e',p));
fixfig;
maxfig;
%set(gca,'xscale','log')

%test mismatch activation per cell
vecNeuronP_Hz = nan(1,intNeurons);
vecNeuronP_Z = nan(1,intNeurons);
for intNeuron=1:intNeurons
	[h,pZ]=ttest(matRespD(1,intNeuron,:),matRespD(2,intNeuron,:));
	vecNeuronP_Z(intNeuron) = pZ;
	
	[h,pHz]=ttest(matRespBins(1,intNeuron,intSingleFrameBin,:),matRespBins(2,intNeuron,intSingleFrameBin,:));
	vecNeuronP_Hz(intNeuron) = pHz;
end

subplot(2,2,4)
scatter(vecNeuronP_Hz,vecNeuronP_Z)
ylabel('p-value ZETA');
xlabel('p-value Hz');
set(gca,'xscale','log','yscale','log')
vecLims = [min([get(gca,'xlim') get(gca,'ylim')]) max([get(gca,'xlim') get(gca,'ylim')])];
xlim(vecLims);
ylim(vecLims);


if 0
	%% save
	strTargetDir = 'F:\Data\Results\ZETA\';
	strFigName = sprintf('SummaryEffectSizeVirtCorr');
	fprintf('Saving figures [%s] ... \n',getTime)
	export_fig([strTargetDir strFigName '.tif']);
	export_fig([strTargetDir strFigName '.pdf']);
	fprintf('\bDone! [%s]\n',getTime);
end