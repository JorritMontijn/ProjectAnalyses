%% overview
%jitter spikes within certain window, then compare real with jittered sd/mean as function of timescale

%% set parameters
cellDataTypes = {'Npx','Sim','ABI','SWN'};%topo, model, allen, nora
intRunDataType = 1;
strRunStim = 'DG';%DG or NM? => superseded to WS by SWN
dblRemOnset = 0; %remove onset period in seconds; 0.125 for sim, 0.25 for npx
cellTypes = {'Real','Poiss','ShuffTid','ShuffTxClass','Uniform'};%,'Poiss','ShuffTid','Shuff','PoissGain','Uniform','ShuffTxClass'};
runHeaderPopTimeCoding;
boolMakeFigs = true;
%vecTimescales = 0.01:0.01:10;%10;1.5;
vecTimescales = linspace(1e-4,1e1,10000);%linspace(1e-2,1e1,1000)
vecJitter = [0 1e3];%[0 (2.^(-9:10))];
intPopSize = inf; %24 (smallest pop of all recs) or inf (uses full pop for each rec)

%% go through recs
tic
for intRec=1:intRecNum %19 || weird: 11
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
		vecOrigStimOnTime = vecStimOnTime + dblRemOnset;
		
		%get data matrix
		[matData,indTuned,cellSpikeTimes,sOut,cellSpikeTimesPerCellPerTrial,indResp] = ...
			SimPrepData(cellSpikeTimesRaw,vecOrigStimOnTime,vecStimOffTime,vecOrientation);
		
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
		vecOrigStimOnTime = vecStimOnTime + dblRemOnset;
		
		%get data matrix
		[matData,indTuned,cellSpikeTimes,sTuning24,cellSpikeTimesPerCellPerTrial] = ...
			NpxPrepData(cellSpikeTimesRaw,vecOrigStimOnTime,vecStimOffTime,vecOrientation);
		
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
	intTrialNum = numel(vecOrigStimOnTime);
	vecOri180 = mod(vecOrientation,180);
	[vecOriIdx,vecUnique,vecPriorDistribution,cellSelect,vecRepetition] = val2idx(vecOri180);
	sTuning = getTuningCurves(matData,vecOri180,0);
	vecNeuronPrefOri = sTuning.matFittedParams(:,1);
	vecNeuronBandwidth = real(sTuning.matBandwidth);
	
	intOriNum = numel(vecUnique);
	intRepNum = min(vecPriorDistribution);
	dblStimDur = median(vecStimOffTime - vecOrigStimOnTime);
	if mean(sum(matData)) < 90%90 / 50
		fprintf('Avg # of spikes per trial was %.1f for %s; skipping...\n',mean(sum(matData)),strThisRec);
		continue;
	else
		fprintf('Running %s... [%s]\n',strThisRec,getTime);
	end
	%% go through types
	clear sAggData;
	vecRunTypes = 1:numel(cellTypes);
	if 0
		%load old data
		load(fullpath(strTargetDataPath,sprintf('Q2cData%s_%s.mat',strThisRec,strOnset)));
		vecRunTypes=(numel(sAggData)+1):numel(cellTypes);
	end
	for intType=vecRunTypes
		strType = cellTypes{intType};
		
		%% load prepped data and ifr
		sSource = load(fullpath(strDataPathT0,sprintf('T0Data_%s%s%s%s%s',strThisRec,strRunType,strRunStim,strType)));
		cellOrigSpikeTimes = sSource.cellSpikeTimes;
		vecPseudoStartT = vecOrigStimOnTime;
		
		%% take only period during stimuli
		cellSpikeTimes = cell(size(cellOrigSpikeTimes));
		for i=1:numel(cellOrigSpikeTimes)
			[vecPseudoSpikeTimes,vecPseudoStartT] = getPseudoSpikeVectors(cellOrigSpikeTimes{i},vecOrigStimOnTime,dblStimDur,true);
			
			dblStart = vecPseudoStartT(1);
			dblEnd = vecPseudoStartT(end)+dblStimDur;
			dblTotDur = dblEnd-dblStart;
			vecPseudoSpikeTimes(vecPseudoSpikeTimes<dblStart | vecPseudoSpikeTimes>dblEnd) = [];
			vecPseudoSpikeTimes = vecPseudoSpikeTimes-dblStart;
			dblLastSpike = max(vecPseudoSpikeTimes);
			
			cellSpikeTimes{i} = vecPseudoSpikeTimes;
		end
		
		%% run for whole pop and for single neurons
		intTimescaleNum = numel(vecTimescales);
		intJitterNum = numel(vecJitter);
		%whole pop
		matMean = nan(intJitterNum,intTimescaleNum);
		matSd = nan(intJitterNum,intTimescaleNum);
		matCV = nan(intJitterNum,intTimescaleNum);
		
		vecR2_Exp = nan(intJitterNum,1);
		vecHalfLife_Exp = nan(intJitterNum,1);
		vecAsymptote_Exp = nan(intJitterNum,1);
		vecScale_Exp = nan(intJitterNum,1);
		matCVFit_Exp = nan(intJitterNum,intTimescaleNum);
		
		vecR2_Root = nan(intJitterNum,1);
		vecAsymptote_Root = nan(intJitterNum,1);
		vecScale_Root = nan(intJitterNum,1);
		vecExponent_Root = nan(intJitterNum,1);
		matCVFit_Root = nan(intJitterNum,intTimescaleNum);
		
		%single neurons
		mat1CV = nan(intTimescaleNum,intJitterNum,intNumN);
		mat1R2_Exp = nan(intJitterNum,intNumN);
		mat1HalfLife_Exp = nan(intJitterNum,intNumN);
		mat1Asymptote_Exp = nan(intJitterNum,intNumN);
		mat1Scale_Exp = nan(intJitterNum,intNumN);
		
		mat1R2_Root = nan(intJitterNum,intNumN);
		mat1Asymptote_Root = nan(intJitterNum,intNumN);
		mat1Scale_Root = nan(intJitterNum,intNumN);
		mat1Exponent_Root = nan(intJitterNum,intNumN);
		
		
		for intJitterIdx=1:numel(vecJitter)
			%% jitter
			dblJitter = vecJitter(intJitterIdx);
			fprintf('Running %s (%d/%d) - %s, jitter %.2e (%d/%d) [%s]\n',strRec,intRec,intRecNum,strType,dblJitter,intJitterIdx,numel(vecJitter),getTime);
			
			%% whole pop
			vecAllSpikeTime = cell2vec(cellSpikeTimes);
			vecApplyJitter = randn(size(vecAllSpikeTime))*dblJitter;
			vecAllSpikeTime = sort(mod(vecAllSpikeTime + vecApplyJitter,dblTotDur));
			
			%get fits
			boolFixedAsymptote = false;
			[sFits,vecCV,vecMean,vecSd] = getFitsCV(vecAllSpikeTime,dblTotDur,vecTimescales,boolFixedAsymptote);
			
			matMean(intJitterIdx,:) = vecMean;
			matSd(intJitterIdx,:) = vecSd;
			matCV(intJitterIdx,:) = vecCV;
			vecR2_Exp(intJitterIdx) = sFits.Exp.R2;
			vecHalfLife_Exp(intJitterIdx) = sFits.Exp.HalfLife;
			vecAsymptote_Exp(intJitterIdx) = sFits.Exp.Asymptote;
			vecScale_Exp(intJitterIdx) = sFits.Exp.Scale;
			matCVFit_Exp(intJitterIdx,:) = sFits.Exp.FitY;
			
			vecR2_Root(intJitterIdx) = sFits.Root.R2;
			vecAsymptote_Root(intJitterIdx) = sFits.Root.Asymptote;
			vecScale_Root(intJitterIdx) = sFits.Root.Scale;
			vecExponent_Root(intJitterIdx) = sFits.Root.Exponent;
			matCVFit_Root(intJitterIdx,:) = sFits.Root.FitY;
			
			%% single neurons
			hTic=tic;
			for intN=1:intNumN
				if toc(hTic) > 5
					fprintf('  Neuron %d/%d [%s]\n',intN,intNumN,getTime);
					hTic=tic;
				end
				vecAllSpikeTime = cellSpikeTimes{intN};
				vecApplyJitter = randn(size(vecAllSpikeTime))*dblJitter;
				vecAllSpikeTime = sort(mod(vecAllSpikeTime + vecApplyJitter,dblTotDur));
				
				%get fits
				boolFixedAsymptote = false;
				[sFitsSingle,vecCV] = getFitsCV(vecAllSpikeTime,dblTotDur,vecTimescales,boolFixedAsymptote);
				mat1CV(:,intJitterIdx,intN) = vecCV;
				mat1R2_Exp(intJitterIdx,intN) = sFitsSingle.Exp.R2;
				mat1HalfLife_Exp(intJitterIdx,intN) = sFitsSingle.Exp.HalfLife;
				mat1Asymptote_Exp(intJitterIdx,intN) = sFitsSingle.Exp.Asymptote;
				mat1Scale_Exp(intJitterIdx,intN) = sFitsSingle.Exp.Scale;
				
				mat1R2_Root(intJitterIdx,intN) = sFitsSingle.Root.R2;
				mat1Asymptote_Root(intJitterIdx,intN) = sFitsSingle.Root.Asymptote;
				mat1Scale_Root(intJitterIdx,intN) = sFitsSingle.Root.Scale;
				mat1Exponent_Root(intJitterIdx,intN) = sFitsSingle.Root.Exponent;
			end
		end
		
		
		% plot
		if boolMakeFigs
			figure;maxfig;
			colormap(parula)
			vecPlotJitterIdx = 1:intJitterNum;%[1 2 ceil(intJitterNum/2) intJitterNum];
			vecPlotJitter = vecJitter(vecPlotJitterIdx);
			vecPlotJitter(1) = 1e-5; %to show on log axis
			subplot(2,4,1)
			matCV = matSd./matMean;
			plot(vecTimescales,matCV(vecPlotJitterIdx,:)')
			cellL=(cellfun(@(x) sprintf('%.2f',x),vec2cell(vecJitter),'UniformOutput',false));
			h=colorbar;
			vecShowJitterTicks = unique(round(linspace(1,numel(vecPlotJitter),5)));
			%vecShowJitterTicks = [1 6 11 15 20];
			h.Ticks = imnorm(vecShowJitterTicks);
			h.TickLabels = cellL(vecShowJitterTicks);
			title(sprintf('Color=jitter; %s - %s',strRec,strType),'interpreter','none');
			set(gca, 'ColorOrder', parula(numel(vecPlotJitterIdx)))
			xlabel('Timescale (bin size in s)');
			ylabel('CV (sd/mu)');
			
			subplot(2,4,2);
			plot(vecTimescales,(matCV(vecPlotJitterIdx,:) - matCVFit_Exp(vecPlotJitterIdx,:))')
			title('exp residuals')
			set(gca, 'ColorOrder', parula(numel(vecPlotJitterIdx)))
			xlabel('Timescale (bin size in s)');
			ylabel('Residuals exp decay fit');
			
			subplot(2,4,3);
			plot(vecPlotJitter,vecR2_Exp)
			title('exp decay r2')
			xlabel('Spike jitter (s)');
			ylabel('R^2 exp decay fit');
			set(gca,'xscale','log');
			
			subplot(2,4,4);
			plot(vecPlotJitter,vecHalfLife_Exp)
			title('exp decay half-life')
			xlabel('Spike jitter (s)');
			ylabel('Half-life decay fit');
			set(gca,'xscale','log');
			
			subplot(2,4,5);
			plot(matMean',matSd')
			title('Sd/mean')
			xlabel('Mean spike count');
			ylabel('Sd spike count');
			set(gca, 'ColorOrder', parula(numel(vecPlotJitterIdx)))
			
			subplot(2,4,6);
			plot(vecTimescales,(matCV(vecPlotJitterIdx,:) - matCVFit_Root(vecPlotJitterIdx,:))')
			title('root residuals')
			set(gca, 'ColorOrder', parula(numel(vecPlotJitterIdx)))
			xlabel('Timescale (bin size in s)');
			ylabel('Residuals root fit');
			
			subplot(2,4,7);
			plot(vecPlotJitter,vecR2_Root)
			title('root r2')
			xlabel('Spike jitter (s)');
			ylabel('R^2 root fit');
			set(gca,'xscale','log');
			
			subplot(2,4,8);
			plot(vecPlotJitter,vecExponent_Root)
			title('root exponent')
			xlabel('Spike jitter (s)');
			ylabel('Exponent of root fit');
			set(gca,'xscale','log');
			fixfig;
			
			%% save figure
			export_fig(fullpath(strFigurePathSR,sprintf('Q2cA_TimescaleJitterCV_%s_%s_%s.tif',strThisRec,strType,strOnset)));
			export_fig(fullpath(strFigurePathSR,sprintf('Q2cA_TimescaleJitterCV_%s_%s_%s.pdf',strThisRec,strType,strOnset)));
		end
		
		%% save data
		sData = struct;
		sData.strRec = strRec;
		sData.strType = strType;
		sData.dblRemOnset = dblRemOnset;
		sData.strThisRec = strThisRec;
		sData.vecJitter = vecJitter;
		sData.vecTimescales = vecTimescales;
		
		%whole pop
		sData.matMean = matMean;
		sData.matSd = matSd;
		sData.matCV = matCV;
		
		sData.vecR2_Exp = vecR2_Exp;
		sData.vecHalfLife_Exp = vecHalfLife_Exp;
		sData.vecAsymptote_Exp = vecAsymptote_Exp;
		sData.vecScale_Exp = vecScale_Exp;
		sData.matCVFit_Exp = matCVFit_Exp;
		
		sData.vecR2_Root = vecR2_Root;
		sData.vecAsymptote_Root = vecAsymptote_Root;
		sData.vecScale_Root = vecScale_Root;
		sData.vecExponent_Root = vecExponent_Root;
		sData.matCVFit_Root = matCVFit_Root;
		
		%single neurons
		sData.mat1CV = mat1CV;
		sData.mat1R2_Exp = mat1R2_Exp;
		sData.mat1HalfLife_Exp = mat1HalfLife_Exp;
		sData.mat1Asymptote_Exp = mat1Asymptote_Exp;
		sData.mat1Scale_Exp = mat1Scale_Exp;
		
		sData.mat1R2_Root = mat1R2_Root;
		sData.mat1Asymptote_Root = mat1Asymptote_Root;
		sData.mat1Scale_Root = mat1Scale_Root;
		sData.mat1Exponent_Root = mat1Exponent_Root;

		sAggData(intType) = sData;
	end
	
	%% save agg data
	save(fullpath(strTargetDataPath,sprintf('Q2cData%s_%s.mat',strThisRec,strOnset)),...
		'sAggData');
	close all;
end
toc
