%% overview
%In the limit of window=0, the poissgain becomes the real distribution while for window=inf it
%becomes the baseline Poisson. With a continuum of window sizes one can therefore check the dominant
%timescales. Within a window, use n=original spikes, and distribute uniform randomly. Make window by
%timescale plot of CV. Do for per-neuron, neuron-id shuffled per spike (no coordination), and for
%whole pop/random neuron (all same avg FR)

%% set parameters
cellDataTypes = {'Npx','Sim','ABI','SWN'};%topo, model, allen, nora
intRunDataType = 1;
strRunStim = 'DG';%DG or NM? => superseded to WS by SWN
dblRemOnset = 0; %remove onset period in seconds; 0.125 for sim, 0.25 for npx
runHeaderPopTimeCoding;
boolMakeFigs = true;
vecTimescales = linspace(1e-2,1e2,1000);%logspace(-2,3,1000);%0.01:0.01:10;%10;1.5;
vecRandomizationWindow = [1.9.^(-9:9) inf]-0.001;
cellRandTypes = {'PerNeuron','ShuffNid','HomoPop'};
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
		vecStimOnTime = vecStimOnTime + dblRemOnset;
		
		%get data matrix
		[matData,indTuned,cellUseSpikeTimes,sOut,cellSpikeTimesPerCellPerTrial,indResp] = ...
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
		[matData,indTuned,cellUseSpikeTimes,sTuning24,cellSpikeTimesPerCellPerTrial] = ...
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
	intNumN = numel(cellUseSpikeTimes);
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
	else
		fprintf('Running %s... [%s]\n',strThisRec,getTime);
	end
	
	%% load prepped data and ifr
	sSource = load(fullpath(strDataPathT0,sprintf('T0Data_%s%s%s%s%s',strThisRec,strRunType,strRunStim,'Real')));
	cellRealSpikeTimes = sSource.cellSpikeTimes;
	%if isempty(vecTime),continue;end
	
	%In the limit of window=0, the poissgain becomes the real distribution while for window=inf it
	%becomes the baseline Poisson. With a continuum of window sizes one can therefore check the dominant
	%timescales. Within a window, use n=original spikes, and distribute uniform randomly. Make window by
	%timescale plot of CV. Do for per-neuron, neuron-id shuffled per spike (no coordination), and for
	%whole pop/random neuron (all same avg FR)
	
	%% take only period during stimuli
	for i=1:numel(cellRealSpikeTimes)
		[vecPseudoSpikeTimes,vecPseudoStartT] = getPseudoSpikeVectors(cellRealSpikeTimes{i},vecStimOnTime,dblStimDur,true);
		cellRealSpikeTimes{i} = vecPseudoSpikeTimes;
	end
	
	%% run
	clear sAggData;
	for intType=1%:numel(cellRandTypes)
		strType = cellRandTypes{intType};
		
		%pre-allocate
		intTimescaleNum = numel(vecTimescales);
		intRandWidthNum = numel(vecRandomizationWindow);
		matMean = nan(intRandWidthNum,intTimescaleNum);
		matSd = nan(intRandWidthNum,intTimescaleNum);
		matCV = nan(intRandWidthNum,intTimescaleNum);
		matCVFit_Exp = nan(intRandWidthNum,intTimescaleNum);
		matCVFit_Root = nan(intRandWidthNum,intTimescaleNum);
		
		vecSlope_Lin = nan(intRandWidthNum,1);
		vecR2_Lin = nan(intRandWidthNum,1);
		
		vecR2_Exp = nan(intRandWidthNum,1);
		vecHalfLife_Exp = nan(intRandWidthNum,1);
		vecAsymptote_Exp = nan(intRandWidthNum,1);
		vecScale_Exp = nan(intRandWidthNum,1);
		
		vecR2_Root = nan(intRandWidthNum,1);
		vecAsymptote_Root = nan(intRandWidthNum,1);
		vecScale_Root = nan(intRandWidthNum,1);
		vecExponent_Root = nan(intRandWidthNum,1);
		
		%% randomize neuron id of each spike
		if strcmp(strType,'PerNeuron')
			%keep original assignment
			cellUseSpikeTimes = cellRealSpikeTimes;
			vecAllSpikeTime = cell2vec(cellUseSpikeTimes);
		elseif strcmp(strType,'ShuffNid')
			%neuron-id shuffled per spike (no coordination), but keep FR per neuron
			intTotS = sum(cellfun(@numel,cellRealSpikeTimes));
			vecAllSpikeTime = nan(1,intTotS);
			vecAllSpikeNeuron = zeros(1,intTotS,'int16');
			intS = 1;
			for intN=1:intNumN
				%add spikes
				intThisS = numel(cellRealSpikeTimes{intN});
				vecSpikeT = cellRealSpikeTimes{intN};
				vecAllSpikeTime(intS:(intS+intThisS-1)) = vecSpikeT;
				vecAllSpikeNeuron(intS:(intS+intThisS-1)) = intN;
				intS = intS + intThisS;
			end
			
			%sort
			[vecAllSpikeTime,vecReorder] = sort(vecAllSpikeTime);
			vecAllSpikeNeuron = vecAllSpikeNeuron(vecReorder);
			
			%randomize ids
			vecAllSpikeNeuron = vecAllSpikeNeuron(randperm(intTotS));
			
			%rebuild spike cell array
			cellUseSpikeTimes = cell(size(cellRealSpikeTimes));
			for intN=1:intNumN
				%add spikes
				cellUseSpikeTimes{intN} = vecAllSpikeTime(vecAllSpikeNeuron==intN);
			end
		elseif strcmp(strType,'HomoPop')
			%randomly assign spikes to neurons (same FR), then randomize times
			vecAllSpikeTime = cell2vec(cellRealSpikeTimes);
			vecRandN = randi(intNumN,size(vecAllSpikeTime));
			
			%rebuild spike cell array
			cellUseSpikeTimes = cell(size(cellRealSpikeTimes));
			for intN=1:intNumN
				%add spikes
				cellUseSpikeTimes{intN} = vecAllSpikeTime(vecRandN==intN);
			end
		else
			not possible
		end
		dblStartT = min(vecAllSpikeTime);
		dblStopT = max(vecAllSpikeTime);
		
		%% randomize spike times per neuron
		for intRandIdx=1:numel(vecRandomizationWindow)
			%msg
			dblRandWindow = vecRandomizationWindow(intRandIdx);
			fprintf('Running %s - %s, %.2e s (%d/%d) [%s]\n',strRec,strType,dblRandWindow,intRandIdx,numel(vecRandomizationWindow),getTime);
		
			%% randomize spike times
			dblInitRand = min(dblRandWindow,1);
			vecWindowEdges = [(dblStartT-dblInitRand*rand(1)):dblRandWindow:dblStopT dblStopT];
			cellRandSpikeTimes = cellUseSpikeTimes;
			for intN=1:intNumN
				vecThisSpikeT = cellUseSpikeTimes{intN};
				vecSpikesPerBin = histcounts(vecThisSpikeT,vecWindowEdges);
				vecNonZeroBins = find(vecSpikesPerBin>0);
				intSpikeC = 1;
				for intBin=1:numel(vecNonZeroBins)
					intBinIdx = vecNonZeroBins(intBin);
					dblStartBin = vecWindowEdges(intBinIdx);
					if intBinIdx == numel(vecSpikesPerBin)
						dblWindowWidth = dblStopT-dblStartBin;
					else
						dblWindowWidth = dblRandWindow;
					end
					intSpikeCount = vecSpikesPerBin(intBinIdx);
					vecNewT = dblStartBin + dblWindowWidth*rand(intSpikeCount,1);
					
					vecAssign = intSpikeC:(intSpikeC+intSpikeCount-1);
					vecThisSpikeT(vecAssign) = vecNewT;
					intSpikeC = intSpikeC + intSpikeCount;
				end
				cellRandSpikeTimes{intN} = vecThisSpikeT;
			end
			%compile and sort
			vecAllSpikeTime = sort(cell2vec(cellRandSpikeTimes));
			
			%% get mean and sd
			vecMean = nan(size(vecTimescales));
			vecSd = nan(size(vecTimescales));
			intK=1;
			for intScale=1:intTimescaleNum
				vecBins = dblStartT:vecTimescales(intScale):dblStopT;
				vecCounts = histcounts( vecAllSpikeTime,vecBins);
				vecMean(intScale) = mean(vecCounts);
				vecSd(intScale) = std(vecCounts);
			end
			
			%% get lin fit for mean/sd
			%fit linear model
			mdl = fitlm(vecMean,vecSd);
			dblSlope_SdMean = mdl.Coefficients.Estimate(2);
			vecFitY =  mdl.Fitted;
			[dblR2_SdMean,dblSS_tot,dblSS_res,dblT,dblP,dblR2_adjusted,dblR2_SE] = getR2(vecSd,vecFitY,intK);
			
			%% fit cv/timescale pop
			fExp = fittype('a+b*exp(-x/c)',...
				'dependent',{'y'},'independent',{'x'},...
				'coefficients',{'a','b','c'});
			
			boolFixedAsymptote = true;
			vecX = vecTimescales';
			vecY = (vecSd./vecMean)';
			vecStartCoeffs = [vecY(end) vecY(1)/vecTimescales(1) 1/2];
			vecUpper = [1e16 1e16 1];
			vecLower = [0 0 0];
			if boolFixedAsymptote
				fRoot = fittype('(1/((b*x)^c))',...
					'dependent',{'y'},'independent',{'x'},...
					'coefficients',{'b','c'});
				vecStartCoeffs = vecStartCoeffs(2:end);
				vecUpper = vecUpper(2:end);
				vecLower = vecLower(2:end);
			else
				fRoot = fittype('a+(1/((b*x)^c))',...
					'dependent',{'y'},'independent',{'x'},...
					'coefficients',{'a','b','c'});
			end
			
			intK=1;
			%fit linear model
			mdl = fitlm(vecX,vecY);
			
			dblSlope_Lin = mdl.Coefficients.Estimate(2);
			vecFitY =  mdl.Fitted;
			[dblR2_Lin,dblSS_tot,dblSS_res,dblT,dblP,dblR2_adjusted,dblR2_SE] = getR2(vecY,vecFitY,intK);
			
			%fit exp model
			vecStartCoeffsExp = [vecY(end) vecY(1)-vecY(end) 0.1];
			[fitobject,gof,output] = fit(vecX,vecY,fExp,'lower',[0 0 1e-6],'upper',[1e16 1e16 1e16],'startpoint',vecStartCoeffsExp);
			vecFitExp = fitobject(vecX);
			[dblR2_Exp,dblSS_tot,dblSS_res,dblT,dblP,dblR2_adjusted,dblR2_SE] = getR2(vecY,vecFitExp,intK);
			dblHalfLife_Exp = fitobject.c*log(2);
			dblAsymptote_Exp = fitobject.a;
			dblScale_Exp = fitobject.b;
			
			%fit root model
			[fitobject,gof,output] = fit(vecX,vecY,fRoot,'lower',vecLower,'upper',vecUpper,'startpoint',vecStartCoeffs);
			vecFitRoot = fitobject(vecX);
			[dblR2_Root,dblSS_tot,dblSS_res,dblT,dblP,dblR2_adjusted,dblR2_SE] = getR2(vecY,vecFitRoot,intK);
			if boolFixedAsymptote
				dblAsymptote_Root=0;
			else
				dblAsymptote_Root = fitobject.a;
			end
			dblScale_Root = fitobject.b;
			dblExponent_Root = fitobject.c;
			
			%% save
			matMean(intRandIdx,:) = vecMean;
			matSd(intRandIdx,:) = vecSd;
			vecSlope_Lin(intRandIdx) = dblSlope_Lin;
			vecR2_Lin(intRandIdx) = dblR2_Lin;
			
			matCV(intRandIdx,:) = vecSd./vecMean;
			matCVFit_Exp(intRandIdx,:) = vecFitExp;
			matCVFit_Root(intRandIdx,:) = vecFitRoot;
			
			vecR2_Exp(intRandIdx) = dblR2_Exp;
			vecHalfLife_Exp(intRandIdx) = dblHalfLife_Exp;
			vecAsymptote_Exp(intRandIdx) = dblAsymptote_Exp;
			vecScale_Exp(intRandIdx) = dblScale_Exp;
			
			vecR2_Root(intRandIdx) = dblR2_Root;
			vecAsymptote_Root(intRandIdx) = dblAsymptote_Root;
			vecScale_Root(intRandIdx) = dblScale_Root;
			vecExponent_Root(intRandIdx) = dblExponent_Root;
		end
		
		
		% plot
		if boolMakeFigs
			figure;maxfig;
			vecPlotRandWin = vecRandomizationWindow;
			vecPlotRandWin(end)=vecPlotRandWin(end-1)*2;
			colormap(parula)
			vecPlotJitters = 1:intRandWidthNum;%[1 2 ceil(intJitterNum/2) intJitterNum];
			subplot(2,4,1)
			matCV = matSd./matMean;
			plot(vecTimescales,matCV(vecPlotJitters,:)')
			cellL=(cellfun(@(x) sprintf('%.2f',x),vec2cell(vecRandomizationWindow),'UniformOutput',false));
			h=colorbar;
			vecShowJitterTicks = round(linspace(1,numel(vecRandomizationWindow),5));
			h.Ticks = imnorm(vecShowJitterTicks);
			h.TickLabels = cellL(vecShowJitterTicks);
			title(sprintf('Color=jitter; %s - %s',strRec,strType),'interpreter','none');
			set(gca, 'ColorOrder', parula(numel(vecPlotJitters)))
			xlabel('Timescale (bin size in s)');
			ylabel('CV (sd/mu)');
			
			subplot(2,4,2);
			plot(vecTimescales,(matCV(vecPlotJitters,:) - matCVFit_Exp(vecPlotJitters,:))')
			title('exp residuals')
			set(gca, 'ColorOrder', parula(numel(vecPlotJitters)))
			xlabel('Timescale (bin size in s)');
			ylabel('Residuals exp decay fit');
			
			subplot(2,4,3);
			plot(vecPlotRandWin,vecR2_Exp)
			title('exp decay r2')
			xlabel('Spike jitter (s)');
			ylabel('R^2 exp decay fit');
			set(gca,'xscale','log');
			
			subplot(2,4,4);
			plot(vecPlotRandWin,vecHalfLife_Exp)
			title('exp decay half-life')
			xlabel('Spike jitter (s)');
			ylabel('Half-life decay fit');
			set(gca,'xscale','log');
			
			subplot(2,4,5);
			plot(matMean',matSd')
			title('Sd/mean')
			xlabel('Mean spike count');
			ylabel('Sd spike count');
			set(gca, 'ColorOrder', parula(numel(vecPlotJitters)))
			
			subplot(2,4,6);
			plot(vecTimescales,(matCV(vecPlotJitters,:) - matCVFit_Root(vecPlotJitters,:))')
			title('root residuals')
			set(gca, 'ColorOrder', parula(numel(vecPlotJitters)))
			xlabel('Timescale (bin size in s)');
			ylabel('Residuals root fit');
			
			subplot(2,4,7);
			plot(vecPlotRandWin,vecR2_Root)
			title('root r2')
			xlabel('Spike jitter (s)');
			ylabel('R^2 root fit');
			set(gca,'xscale','log');
			
			subplot(2,4,8);
			plot(vecPlotRandWin,vecExponent_Root)
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
		sData.vecRandomizationWindow = vecRandomizationWindow;
		sData.vecTimescales = vecTimescales;
		
		sData.matMean = matMean;
		sData.matSd = matSd;
		sData.vecSlope_Lin = vecSlope_Lin;
		sData.matCV = matCV;
		sData.matCVFit_Exp = matCVFit_Exp;
		sData.matCVFit_Root = matCVFit_Root;
		sData.vecR2_Exp = vecR2_Exp;
		
		sData.vecHalfLife_Exp = vecHalfLife_Exp;
		sData.vecAsymptote_Exp = vecAsymptote_Exp;
		sData.vecScale_Exp = vecScale_Exp;
		
		sData.vecR2_Root = vecR2_Root;
		sData.vecAsymptote_Root = vecAsymptote_Root;
		sData.vecScale_Root = vecScale_Root;
		sData.vecExponent_Root = vecExponent_Root;
		
		sAggData(intType) = sData;
	end
	return
	%% save agg data
	save(fullpath(strTargetDataPath,sprintf('Q2cData%s_%s.mat',strThisRec,strOnset)),...
		'sAggData');
	close all;
end
toc
