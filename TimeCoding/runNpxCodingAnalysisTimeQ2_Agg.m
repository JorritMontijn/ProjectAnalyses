
%% set parameters
cellDataTypes = {'Npx','Sim','ABI','SWN'};%topo, model, allen, nora
intRunDataType = 1;
strRunStim = 'DG';%DG or NM? => superseded to WS by SWN
cellTypes = {'Real','Poiss','ShuffTid','Shuff','PoissGain','Uniform'};
dblRemOnset = 0; %remove onset period in seconds; 0.125 for sim, 0.25 for npx
runHeaderPopTimeCoding;

%% go through recs
tic
for intRec=1:intRecNum %19 || weird: 11
	close all;
	%% prep ABI or Npx data
	if strcmp(strRunType,'ABI')
		error to be updated
		runRecPrepABI;
		strThisRec = strRec;
	elseif strcmp(strRunType,'Sim')
		%load
		runRecPrepSim;
		
		%edit vars
		strThisRec = strRec;
		strDataPathT0=strDataPathSimT0;
		vecOri180 = mod(vecOrientation,180)*2;
		vecStimIdx = vecOri180;
		
		%% move onset
		%remove first x ms
		vecStimOnTime = vecStimOnTime + dblRemOnset;
		
		%get data matrix
		[matData,indTuned,cellSpikeTimes,sOut,cellSpikeTimesPerCellPerTrial,indResp] = ...
			SimPrepData(cellSpikeTimesRaw,vecStimOnTime,vecStimOffTime,vecOrientation);
		
		%get cell props
		%vecNeuronPrefOri = pi-[sLoad.sNeuron.PrefOri];
		vecNeuronType = [sLoad.sNeuron.Types]; %1=pyr,2=interneuron
		
	elseif strcmp(strRunType,'Npx') || strcmp(strRunType,'SWN')
		%prep
		runRecPrepNpx;
		strThisRec = strRec;
		strDataPathT0 = strTargetDataPath;
		
		%% move onset
		%remove first x ms
		vecStimOnTime = vecStimOnTime + dblRemOnset;
		
		%get data matrix
		[matData,indTuned,cellSpikeTimes,sTuning24,cellSpikeTimesPerCellPerTrial] = ...
			NpxPrepData(cellSpikeTimesRaw,vecStimOnTime,vecStimOffTime,vecOrientation);
		
		%get cell props
		vecNeuronType = ones(size(indTuned)); %1=pyr,2=interneuron
		%narrow vs broad not done
	else
		error impossible
	end
	
	%get ori vars
	intTunedN = sum(indTuned);
	intRespN = size(matData,1);
	intNumN = numel(cellSpikeTimes);
	intTrialNum = numel(vecStimOnTime);
	vecOri180 = mod(vecOrientation,180);
	[vecOriIdx,vecUnique,vecPriorDistribution,cellSelect,vecRepetition] = val2idx(vecOri180);
	sTuning = getTuningCurves(matData,vecOri180,0);
	vecNeuronPrefOri = sTuning.matFittedParams(:,1);
	vecNeuronBandwidth = real(sTuning.matBandwidth);
	
	intOriNum = numel(vecUnique);
	intRepNum = min(vecPriorDistribution);
	dblStimDur = median(vecStimOffTime - vecStimOnTime);
	if mean(sum(matData)) < 90%90 / 50
		fprintf('Avg # of spikes per trial was %.1f for %s; skipping...\n',mean(sum(matData)),strThisRec);
		continue;
	end
	
	
	%types: Real, UniformTrial, ShuffTid, PoissGain
	clear sAggData;
	for intType=1:numel(cellTypes)
		strType = cellTypes{intType};
		
		%% load prepped data and ifr
		sSource = load(fullpath(strDataPathT0,sprintf('T0Data_%s%s%s%s%s',strThisRec,strRunType,strRunStim,strType)));
		vecTime = sSource.vecTime;
		vecIFR = sSource.vecIFR;
		vecAllSpikeTime = sSource.vecAllSpikeTime;
		vecAllSpikeNeuron = sSource.vecAllSpikeNeuron;
		cellSpikeTimes = sSource.cellSpikeTimes;
		if isempty(vecTime),continue;end
		
		%% build trial-neuron cell matrix
		cellSpikeTimesPerCellPerTrial = cell(intNumN,intOrigTrialNum);
		for intN=1:intNumN
			%real
			[vecTrialPerSpike,vecTimePerSpike] = getSpikesInTrial(cellSpikeTimes{intN},vecStimOnTime,dblStimDur);
			for intTrial=1:intOrigTrialNum
				vecSpikeT = vecTimePerSpike(vecTrialPerSpike==intTrial);
				cellSpikeTimesPerCellPerTrial{intN,intTrial} = vecSpikeT;
			end
		end
		
		%% calc IFRs per trial
		cellIFR_perTrial = cell(intTrialNum,1);
		cellTimeIFR_perTrial = cell(intTrialNum,1);
		cellISI_perTrial = cell(intTrialNum,1);
		for intTrial=1:intTrialNum
			%real
			vecAllSpikes = sort(cell2vec(cellSpikeTimesPerCellPerTrial(:,intTrial)));
			vecISI0 = [vecAllSpikes(2:end) - vecAllSpikes(1:(end-1)); inf];
			vecAllSpikes(vecISI0==0)=vecAllSpikes(vecISI0==0)-(10^-5)*rand();
			vecAllSpikes = uniquetol(vecAllSpikes,1e-7);
			[vecTimeIFR,vecIFR] = getIFR(vecAllSpikes,0,dblMaxDur,[],[],[],0);
			cellIFR_perTrial{intTrial} = vecIFR;
			cellTimeIFR_perTrial{intTrial} = vecTimeIFR;
			cellISI_perTrial{intTrial} = diff(vecAllSpikes);
		end
		
		%% linear cv = constant SperTrial/MperTrial
		boolPlot = true;
		vecRperTrial = zeros(intTrialNum,1);
		vecSperTrial = zeros(intTrialNum,1);
		vecHperTrial = zeros(intTrialNum,1);
		vecLperTrial = zeros(intTrialNum,1);
		vecMperTrial = zeros(intTrialNum,1);
		intQuantiles = 10;
		
		for intTrial=1:intTrialNum
			%real
			vecAllSpikes = sort(cell2vec(cellSpikeTimesPerCellPerTrial(:,intTrial)));
			if numel(vecAllSpikes) < 5,continue;end
			vecISI0 = [vecAllSpikes(2:end) - vecAllSpikes(1:(end-1)); inf];
			vecAllSpikes(vecISI0==0)=vecAllSpikes(vecISI0==0)-(10^-5)*rand();
			vecAllSpikes = uniquetol(vecAllSpikes,1e-7);
			intSpikes = numel(vecAllSpikes);
			vecD = diff(vecAllSpikes);
			vecD1 = vecD(1:(end-1));
			vecD2 = vecD(2:end);
			[r,p]=corr(vecD1,vecD2);
			vecRperTrial(intTrial) = r;
			vecIFR = cellIFR_perTrial{intTrial};
			vecTimeIFR = cellTimeIFR_perTrial{intTrial};
			vecISI = cellISI_perTrial{intTrial};
			
			vecR_sorted = sort(vecIFR);
			intHighQ = max(1,round(numel(vecR_sorted)/intQuantiles));
			intLowQ = max(1,round(numel(vecR_sorted)/intQuantiles));
			vecHperTrial(intTrial) = nanmean(vecR_sorted((1+end-1*intHighQ):(end-0*intHighQ)));
			vecLperTrial(intTrial) = nanmean(vecR_sorted((1+0*intLowQ):(1*intLowQ)));
			vecMperTrial(intTrial) = nanmean(vecR_sorted);
			vecSperTrial(intTrial) = nanstd(vecR_sorted);
		end
		
		%% calc sparseness
		matResp = cellfun(@numel,cellSpikeTimesPerCellPerTrial);
		intNeurons = size(matResp,1);
		vecPopSparseness = 1-mean(matResp,1).^2 ./ sum((matResp.^2)./intNeurons,1);
		dblActBinW = 50;
		vecActBins = 0:dblActBinW:700;
		vecActBinsC = vecActBins(2:end)-dblActBinW/2;
		
		%%
		figure;maxfig;
		subplot(2,3,1)
		histx(vecMperTrial)
		ylabel('Count (trials)')
		xlabel('Population activity, mean over time (Hz)')
		title(strThisRec,'interpreter','none');
		
		subplot(2,3,2)
		histx(vecSperTrial)
		ylabel('Count (trials)')
		xlabel('Sd of population activity, mean over time (Hz)')
		
		subplot(2,3,3)
		[vecCounts,vecMeansV,vecSDsV] = makeBins(vecMperTrial,vecSperTrial,vecActBins);
		indPlotBins = vecCounts>10;
		hold on
		scatter(vecMperTrial,vecSperTrial(:),[],[0.5 0.5 1],'.');
		errorbar(vecActBinsC(indPlotBins),vecMeansV(indPlotBins),vecSDsV(indPlotBins)./sqrt(vecCounts(indPlotBins)),'color',[0 0 1])
		xlabel('Mean population activity (Hz)')
		ylabel('Sd of population activity (Hz)')
		title([strType ';sample=mu/sd of IFR over all spikes in trial']);
		
		subplot(2,3,4)
		vecCVperTrial = vecSperTrial./vecMperTrial;
		[vecCounts,vecMeansCV,vecSDsCV] = makeBins(vecMperTrial,vecCVperTrial,vecActBins);
		hold on
		scatter(vecMperTrial,vecCVperTrial(:),[],[0.5 0.5 1],'.');
		errorbar(vecActBinsC(indPlotBins),vecMeansCV(indPlotBins),vecSDsCV(indPlotBins)./sqrt(vecCounts(indPlotBins)),'color',[0 0 1])
		xlabel('Mean population activity (Hz)')
		ylabel('CV of population activity (Hz)')
		fixfig;
		
		%[r2b,p2b]=corr(vecMperTrial,vecCVperTrial);
		%title(sprintf('Corr(M,Sd)=%.3f, p=%.1e',r2b,p2b));
		
		%%
		export_fig(fullpath(strFigurePathSR,sprintf('Q2A_PopActStatistics_%s_%s_%s.tif',strThisRec,strType,strOnset)));
		export_fig(fullpath(strFigurePathSR,sprintf('Q2A_PopActStatistics_%s_%s_%s.pdf',strThisRec,strType,strOnset)));
		
		
		%% calculate activity during low/high epochs
		%are only cells tuned to the orientation active during low phases? Is this different from high 5%?
		
		%pre-alloc
		intQuantileNum = 3; %5=20%
		
		%pre-allocate
		matAggR = zeros(3,intTrialNum,intTunedN);%save middle/lower/upper
		for intTrial=1:intTrialNum
			%% define data
			%get IFRs
			vecTrialIFR = cellIFR_perTrial{intTrial}(2:(end-1));
			vecTrialTimeIFR = cellTimeIFR_perTrial{intTrial}(2:(end-1));
			cellSpikes = cellSpikeTimesPerCellPerTrial(:,intTrial);
			
			
			if numel(vecTrialIFR)<intQuantileNum,continue;end
			%% divide quantiles
			[vecIFR_sorted,vecReorder] = sort(vecTrialIFR);
			vecTimeIFR_sorted = vecTrialTimeIFR(vecReorder);
			intSamples = numel(vecIFR_sorted);
			intEndLow = max(1,round(intSamples/intQuantileNum));
			dblLow_UpperBound = vecIFR_sorted(intEndLow);
			
			intStartHigh = intSamples-intEndLow+1;
			dblHigh_LowerBound = vecIFR_sorted(intStartHigh);
			vecHighLowIdx = ones(1,intSamples);
			vecHighLowIdx(vecTrialIFR<dblLow_UpperBound) = 2; %low
			vecHighLowIdx(vecTrialIFR>dblHigh_LowerBound) = 3; %high
			
			%plot
			if 0
				%% plot
				plot(vecTrialTimeIFR,vecTrialIFR,'k');
				hold on
				scatter(vecTrialTimeIFR(vecHighLowIdx==2),vecTrialIFR(vecHighLowIdx==2),'b.');
				scatter(vecTrialTimeIFR(vecHighLowIdx==3),vecTrialIFR(vecHighLowIdx==3),'r.');
				hold off
			end
			
			%% determine epochs
			vecChanges = find(diff(vecHighLowIdx)~=0);
			intEpochs = numel(vecChanges)+1;
			vecEpochType = nan(1,intEpochs);
			vecEpochStarts = nan(1,intEpochs);
			vecEpochStops = nan(1,intEpochs);
			dblLastEpoch = 0;
			for intEpoch=1:(intEpochs-1)
				intChangeIdx = vecChanges(intEpoch);
				dblNewEpoch = (vecTrialTimeIFR(intChangeIdx) + vecTrialTimeIFR(intChangeIdx+1))/2;
				vecEpochType(intEpoch) = vecHighLowIdx(intChangeIdx);
				vecEpochStarts(intEpoch) = dblLastEpoch;
				vecEpochStops(intEpoch) = dblNewEpoch;
				dblLastEpoch = dblNewEpoch;
			end
			vecEpochType(intEpochs) = vecHighLowIdx(end);
			vecEpochStarts(intEpochs) = dblNewEpoch;
			vecEpochStops(intEpochs) = dblMaxDur;
			
			%% assign epochs
			for intNeuron=1:intTunedN
				vecSpikes = cellSpikes{intNeuron};
				
				%do stuff here
				vecSpikeQ = zeros(size(vecSpikes));
				vecCountsPerType = zeros(1,3);
				for intSpikeIdx=1:numel(vecSpikes)
					intEpochType = vecEpochType(vecEpochStarts < vecSpikes(intSpikeIdx) & vecEpochStops > vecSpikes(intSpikeIdx));
					vecCountsPerType(intEpochType) = vecCountsPerType(intEpochType) + 1;
				end
				
				%save
				matAggR(:,intTrial,intNeuron) = vecCountsPerType; %low
			end
		end
		
		% mean per trial per quantile
		matLowR = squeeze(matAggR(2,:,:));
		matHighR = squeeze(matAggR(3,:,:));
		vecOri360 = vecOri180*2;
		vecOriLow = vecOri360;
		vecOriHigh = vecOri360;
		
		% ori tuning is stable across quantiles
		sOutLow = getTuningCurves(matLowR',vecOri360,0);
		vecPrefRadLow = sOutLow.matFittedParams(:,1);
		sOutHigh = getTuningCurves(matHighR',vecOri360,0);
		vecPrefRadHigh = sOutHigh.matFittedParams(:,1);
		%ori tuning is stable
		
		% ori tuning diff within low
		intHalfTrials = floor(intTrialNum/2);
		sOutLow1 = getTuningCurves(matLowR(1:intHalfTrials,:)',vecOri360(1:intHalfTrials),0);
		vecPrefRadLow1 = sOutLow1.matFittedParams(:,1);
		sOutLow2 = getTuningCurves(matLowR((intHalfTrials+1):end,:)',vecOri360((intHalfTrials+1):end),0);
		vecPrefRadLow2 = sOutLow2.matFittedParams(:,1);
		% ori tuning diff within high
		intHalfTrials = floor(intTrialNum/2);
		sOutHigh1 = getTuningCurves(matHighR(1:intHalfTrials,:)',vecOri360(1:intHalfTrials),0);
		vecPrefRadHigh1 = sOutHigh1.matFittedParams(:,1);
		sOutHigh2 = getTuningCurves(matHighR((intHalfTrials+1):end,:)',vecOri360((intHalfTrials+1):end),0);
		vecPrefRadHigh2 = sOutHigh2.matFittedParams(:,1);
		
		vecDiffHL1 = abs(circ_dist(vecPrefRadHigh1,vecPrefRadLow2));
		vecDiffHL2 = abs(circ_dist(vecPrefRadHigh2,vecPrefRadLow1));
		vecDiffHL = rad2deg(vecDiffHL1);%+vecDiffHL2)/2;
		vecDiffLL = rad2deg(abs(circ_dist(vecPrefRadLow1,vecPrefRadLow2)));
		vecDiffHH = rad2deg(abs(circ_dist(vecPrefRadHigh1,vecPrefRadHigh2)));
		
		% plot
		figure;maxfig;
		subplot(2,3,1);
		vecEdges = 0:90:360;
		[matCounts,matValMeans,matValSDs,cellVals,cellIDs] = ...
			makeBins2(vecPrefRadLow,vecPrefRadHigh,ones(size(vecPrefRadHigh)),vecEdges,vecEdges);
		%imagesc(matCounts);axis xy
		scatter(rad2deg(vecPrefRadLow),rad2deg(vecPrefRadHigh),'x');
		xlim([-10 370]);
		set(gca,'xtick',0:90:360);
		set(gca,'ytick',0:90:360);
		ylim([-10 370]);
		xlabel('Preferred ori low q')
		ylabel('Preferred ori high q')
		title([strThisRec,strType,strOnset],'interpreter','none');
		
		% ori tuning diff within low
		subplot(2,3,2);
		scatter(rad2deg(vecPrefRadLow1),rad2deg(vecPrefRadLow2),[],vecDiffLL,'x');
		xlim([-10 370]);
		set(gca,'xtick',0:90:360);
		set(gca,'ytick',0:90:360);
		ylim([-10 370]);
		xlabel('Preferred ori low q, 1st half')
		ylabel('Preferred ori low q, 2nd half')
		fixfig;
		
		
		% ori tuning diff within high
		subplot(2,3,3);
		scatter(rad2deg(vecPrefRadHigh1),rad2deg(vecPrefRadHigh2),[],vecDiffHH,'x');
		xlim([-10 370]);
		set(gca,'xtick',0:90:360);
		set(gca,'ytick',0:90:360);
		ylim([-10 370]);
		xlabel('Preferred ori high q, 1st half')
		ylabel('Preferred ori high q, 2nd half')
		fixfig;
		
		% ori tuning diff within high
		subplot(2,3,4);
		scatter(rad2deg(vecPrefRadHigh1),rad2deg(vecPrefRadLow2),[],vecDiffHL,'x');
		xlim([-10 370]);
		set(gca,'xtick',0:90:360);
		set(gca,'ytick',0:90:360);
		ylim([-10 370]);
		xlabel('Preferred ori high q, 1st half')
		ylabel('Preferred ori low q, 2nd half')
		fixfig;
		
		dblBinS = 22.5;
		vecBinE = 0:dblBinS:180;
		vecBinC = vecBinE(2:end)-dblBinS/2;
		vecCHL = histcounts(vecDiffHL,vecBinE);
		vecCLL = histcounts(vecDiffLL,vecBinE);
		vecCHH = histcounts(vecDiffHH,vecBinE);
		
		[h,pHL_LL]=ttest(vecDiffHL,vecDiffLL);
		[h,pHL_HH]=ttest(vecDiffHL,vecDiffHH);
		[h,pHH_LL]=ttest(vecDiffHH,vecDiffLL);
		subplot(2,3,5)
		hold on
		plot([0.5 3.5],[90 90],'--','color',[0.5 0.5 0.5]);
		errorbar(1,mean(vecDiffHL),std(vecDiffHL)./sqrt(numel(vecDiffHL)),'xk');
		errorbar(2,mean(vecDiffLL),std(vecDiffLL)./sqrt(numel(vecDiffLL)),'xb');
		errorbar(3,mean(vecDiffHH),std(vecDiffHH)./sqrt(numel(vecDiffHH)),'xr');
		hold off
		ylabel('Angular diff. pref. ori.');
		set(gca,'xtick',[1 2 3],'xticklabel',{'High-low','Low-low','High-high'});
		title(sprintf('T-tests: HL-LL,p=%.3f, HL-HH,p=%.3f, HH-LL,p=%.3f',pHL_LL,pHL_HH,pHH_LL));
		fixfig;
		
		%%
		export_fig(fullpath(strFigurePathSR,sprintf('Q2B_Qsplit_OriCoding_%s_%s_%s.tif',strThisRec,strType,strOnset)));
		export_fig(fullpath(strFigurePathSR,sprintf('Q2B_Qsplit_OriCoding_%s_%s_%s.pdf',strThisRec,strType,strOnset)));
		
		
		%%
		sData = struct;
		sData.strType = strType;
		sData.dblRemOnset = dblRemOnset;
		sData.strThisRec = strThisRec;
		sData.cellIFR_perTrial = cellIFR_perTrial;
		sData.cellTimeIFR_perTrial = cellTimeIFR_perTrial;
		sData.cellISI_perTrial = cellISI_perTrial;
		sData.vecRperTrial = vecRperTrial;
		sData.vecSperTrial = vecSperTrial;
		sData.vecHperTrial = vecHperTrial;
		sData.vecLperTrial = vecLperTrial;
		sData.vecMperTrial = vecMperTrial;
		sData.vecPopSparseness = vecPopSparseness;
		%add to superstructure
		sAggData(intType) = sData;
	end
	
	%% save agg data
	save(fullpath(strTargetDataPath,sprintf('Q2Data%s_%s.mat',strThisRec,strOnset)),...
		'sAggData');
	
	%% plot real vs shuffled
	cellMarkers = {'x','o'};
	intReal = 1;
	intShuff = 3;%shufftid
	vecMperTrial = sAggData(intReal).vecMperTrial;
	vecSperTrial = sAggData(intReal).vecSperTrial;
	vecLperTrial = sAggData(intReal).vecLperTrial;
	vecHperTrial = sAggData(intReal).vecHperTrial;
	vecPopSparseness = sAggData(intReal).vecPopSparseness;
	
	strShuff = sAggData(intShuff).strType;
	vecMperTrial_S = sAggData(intShuff).vecMperTrial;
	vecSperTrial_S = sAggData(intShuff).vecSperTrial;
	vecLperTrial_S = sAggData(intShuff).vecLperTrial;
	vecHperTrial_S = sAggData(intShuff).vecHperTrial;
	vecPopSparseness_S = sAggData(intShuff).vecPopSparseness;
	
	figure;maxfig;
	[rML,pML]=corr(vecMperTrial(:),vecLperTrial(:));
	[rMH,pMH]=corr(vecMperTrial,vecHperTrial(:));
	[rHL,pHL]=corr(vecLperTrial(:),vecHperTrial(:));
	
	subplot(2,3,1)
	[vecCounts,vecMeansL,vecSDsL] = makeBins(vecMperTrial,vecLperTrial,vecActBins);
	[vecCounts,vecMeansH,vecSDsH] = makeBins(vecMperTrial,vecHperTrial,vecActBins);
	[vecCounts_S,vecMeansL_S,vecSDsL_S] = makeBins(vecMperTrial_S,vecLperTrial_S,vecActBins);
	[vecCounts_S,vecMeansH_S,vecSDsH_S] = makeBins(vecMperTrial_S,vecHperTrial_S,vecActBins);hold on
	indPlotBins = vecCounts>10;
	
	scatter(cat(1,vecMperTrial_S,vecMperTrial_S),cat(1,vecLperTrial_S(:),vecHperTrial_S(:)),[],[0.7 0.7 0.7],'.');
	scatter(vecMperTrial,vecLperTrial(:),[],[0.5 0.5 1],'.');
	scatter(vecMperTrial,vecHperTrial(:),[],[1 0.5 0.5],'.');
	
	errorbar(vecActBinsC(indPlotBins),vecMeansL_S(indPlotBins),vecSDsL_S(indPlotBins)./sqrt(vecCounts_S(indPlotBins)),'color',[0.5 0.5 0.5])
	errorbar(vecActBinsC(indPlotBins),vecMeansH_S(indPlotBins),vecSDsH_S(indPlotBins)./sqrt(vecCounts_S(indPlotBins)),'color',[0.5 0.5 0.5])
	errorbar(vecActBinsC(indPlotBins),vecMeansL(indPlotBins),vecSDsL(indPlotBins)./sqrt(vecCounts(indPlotBins)),'color',[0 0 1])
	errorbar(vecActBinsC(indPlotBins),vecMeansH(indPlotBins),vecSDsH(indPlotBins)./sqrt(vecCounts(indPlotBins)),'color',[1 0 0])
	
	hold off
	fixfig;
	xlabel('Mean firing rate during trial (Hz)');
	ylabel('FR in upper/lower 5% (Hz)');
	title(['Real rate fluctuations vs ' strShuff])
	fixfig;
	
	subplot(2,3,2)
	dblStep = 0.05;
	vecBinE = -1:dblStep:1.5;
	vecBinC = vecBinE(2:end)-dblStep/2;
	vecCountsLow = histcounts((vecLperTrial-vecLperTrial_S)./vecLperTrial_S,vecBinE);
	vecCountsHigh = histcounts((vecHperTrial-vecHperTrial_S)./vecHperTrial_S,vecBinE);
	plot(vecBinC*100,vecCountsLow,'color',[0 0 1]);
	hold on
	plot(vecBinC*100,vecCountsHigh,'color',[1 0 0]);
	hold off
	xlabel(['% change in FR over ' strShuff]);
	legend({'Lowest 5%','Highest 5%'});
	ylabel('Number of trials (count)');
	vecL = (vecLperTrial-vecLperTrial_S)./vecLperTrial_S;
	vecH = (vecHperTrial-vecHperTrial_S)./vecHperTrial_S;
	[h,pL]=ttest(vecL);
	[h,pH]=ttest(vecH);
	title(sprintf('Low, avg=%.1f%%, p=%.1e; high, avg=+%.1f%%, p=%.1e',mean(vecL)*100,pL,mean(vecH)*100,pH));
	fixfig;
	
	%real
	r1=corr(vecLperTrial(:),vecPopSparseness(:));
	r2=corr(vecMperTrial(:),vecPopSparseness(:));
	[rSpH,pSpH]=corr(vecMperTrial(:),vecPopSparseness(:));
	vecActBinsH = 0:dblActBinW:1700;
	vecActBinsHC = vecActBinsH(2:end)-dblActBinW/2;
	[vecCounts2,vecMeans2,vecSDs2] = makeBins(vecMperTrial,vecPopSparseness,vecActBinsH);
	%shuff
	r1_S=corr(vecLperTrial_S(:),vecPopSparseness_S(:));
	r2_S=corr(vecMperTrial_S(:),vecPopSparseness_S(:));
	[rSpH_S,pSpH_S]=corr(vecMperTrial_S(:),vecPopSparseness_S(:));
	vecActBinsH_S = 0:dblActBinW:1700;
	vecActBinsHC_S = vecActBinsH_S(2:end)-dblActBinW/2;
	[vecCounts2_S,vecMeans2_S,vecSDs2_S] = makeBins(vecMperTrial_S,vecPopSparseness_S,vecActBinsH_S);
	
	%fit
	indPlotBins2 = vecCounts2>10 | vecCounts2_S>10;
	
	mdl = fitlm(vecMperTrial,vecPopSparseness);
	ci = coefCI(mdl);
	[ypred,yci] = predict(mdl,vecActBinsHC(indPlotBins2)');
	
	mdl_S = fitlm(vecMperTrial_S,vecPopSparseness_S);
	ci_S = coefCI(mdl_S);
	[ypred_S,yci_S] = predict(mdl_S,vecActBinsHC_S(indPlotBins2)');
	
	subplot(2,3,3)
	hold on;
	scatter(vecMperTrial,vecPopSparseness,[],1-(1-lines(1))*(2/3),'.');
	scatter(vecMperTrial_S,vecPopSparseness_S,[],[0.7 0.7 0.7],'.');
	plot(vecActBinsHC(indPlotBins2),yci(:,1),'--','color',lines(1));
	plot(vecActBinsHC(indPlotBins2),ypred,'color',lines(1));
	plot(vecActBinsHC(indPlotBins2),yci(:,2),'--','color',lines(1));
	plot(vecActBinsHC_S(indPlotBins2),yci_S(:,1),'--','color',[0.5 0.5 0.5]);
	plot(vecActBinsHC_S(indPlotBins2),ypred_S,'color',[0.5 0.5 0.5]);
	plot(vecActBinsHC_S(indPlotBins2),yci_S(:,2),'--','color',[0.5 0.5 0.5]);
	hold off;
	hold off;
	ylabel('Pop. sparseness per trial');
	xlabel('Mean pop. firing rate per trial (Hz) ');
	title(sprintf('Real=%.2f, p=%.1e; %s=%.2f, p=%.1e',rSpH,pSpH,strType,rSpH_S,pSpH_S));
	fixfig;
	
	%split population in highest 50/lowest 50
	intUseUpperCells = min(sum(matResp>0,1));
	vecHighAct = nan(intTrialNum,1);
	vecLowAct =  nan(intTrialNum,1);
	vecQuantiles = [1/3 1/2 2/3];
	vecQuantileIdx = round(vecQuantiles*intNeurons);
	matQuantileAct = nan(intTrialNum,numel(vecQuantileIdx));
	vecMeanOfActiveCells = nan(intTrialNum,1);
	for intTrial=1:intTrialNum
		vecR = sort(matResp(:,intTrial));
		matQuantileAct(intTrial,:) = vecR(vecQuantileIdx);
		vecMeanOfActiveCells(intTrial) = mean(vecR((end-intUseUpperCells+1):end));
	end
	[vecHsorted,vecReorder]=sort(vecHperTrial);
	
	%save fig
	export_fig(fullpath(strFigurePathSR,sprintf('Q2C_QuantileDeviations%s_%s_%s.tif',strThisRec,strOnset)));
	export_fig(fullpath(strFigurePathSR,sprintf('Q2C_QuantileDeviations%s_%s_%s.pdf',strThisRec,strOnset)))
	
end
toc
