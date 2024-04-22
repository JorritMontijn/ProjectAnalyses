%% overview
%jitter spikes within certain window, then compare real with jittered sd/mean as function of timescale

%% set parameters
cellDataTypes = {'Npx','Sim','ABI','SWN'};%topo, model, allen, nora
intRunDataType = 1;
strRunStim = 'DG';%DG or NM? => superseded to WS by SWN
dblRemOnset = 0; %remove onset period in seconds; 0.125 for sim, 0.25 for npx
runHeaderPopTimeCoding;
boolMakeFigs = true;
vecTimescales = 0.01:0.01:10;%10;1.5;
vecJitter = [0 0.01:0.01:0.1 0.2:0.1:1];
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
	else
		fprintf('Running %s... [%s]\n',strThisRec,getTime);
	end
	%% go through types
	clear sAggData;
	figure;maxfig;hold on
	matCol=lines(6);
	vecHandles = [];
	for intType=1:numel(cellTypes)
		strType = cellTypes{intType};
		fprintf('Running %s (%d/%d) - %s [%s]\n',strRec,intRec,intRecNum,strType,getTime);
		
		%% load prepped data and ifr
		sSource = load(fullpath(strDataPathT0,sprintf('T0Data_%s%s%s%s%s',strThisRec,strRunType,strRunStim,strType)));
		cellSpikeTimes = sSource.cellSpikeTimes;
		%if isempty(vecTime),continue;end
		
		%% take only period during stimuli
		for i=1:numel(cellSpikeTimes)
			[vecPseudoSpikeTimes,vecPseudoStartT] = getPseudoSpikeVectors(cellSpikeTimes{i},vecStimOnTime,dblStimDur,true);
			cellSpikeTimes{i} = vecPseudoSpikeTimes;
		end
		vecAllSpikeTime = sort(cell2vec(cellSpikeTimes));
		
		%true uniform
		%vecAllSpikeTime = sort(rand(size(vecAllSpikeTime))*range(vecAllSpikeTime)+min(vecAllSpikeTime));
		vecISI = diff(sort(vecAllSpikeTime));
		
		vecBinEdges = 0:0.001:0.02;
		vecBinCenters = vecBinEdges(2:end)-diff(vecBinEdges)/2;
		vecCounts = histcounts(vecISI,vecBinEdges);
		
		fExp = @(b,c,x) b*exp(-x/c);
		fitExp = fittype(fExp,...
			'dependent',{'y'},'independent',{'x'},...
			'coefficients',{'b','c'});
		
		vecX = 1000*vecBinCenters';
		vecY = vecCounts';
		
		subplot(2,3,intType);hold on
		vecHandles(intType) = plot(vecX,vecY,'color',matCol(intType,:));
		options = fitoptions;
		vecStartCoeffsExp = [vecY(1) 0.1];
		%fExp = @(p,x) p(1)*exp(-x/p(2));
		%[pFit,Resnorm,FVAL,EXITFLAG,OUTPUT,LAMBDA,JACOB] = curvefitfun(fExp,vecStartCoeffsExp,vecX,vecY,[0 0],[1e10 1e6]);
		%vecFitY3 = fExp(pFit,vecX);
		
		[fitobject,gof,output] = fit(vecX,vecY,fitExp,'lower',[0 0],'upper',[1e10 1e10],'startpoint',vecStartCoeffsExp);
		vecFitExp = fitobject(vecX);
		[dblR2_Exp,dblSS_tot,dblSS_res,dblT,dblP,dblR2_adjusted,dblR2_SE] = getR2(vecY,vecFitExp,1);
		dblHalfLife_Exp = fitobject.c*log(2);
		dblScale_Exp = fitobject.b;
		
		plot(vecX,vecFitExp,'--','color',matCol(intType,:));
		if intType==1
			vecFitReal = vecFitExp;
		else
			plot(vecX,vecFitReal,'--','color',matCol(1,:));
		end
		title(strType);
		set(gca,'yscale','log');
		xlabel('Inter-spike interval (ms)');
		ylabel('# of spikes');
	end
	%set(gca,'yscale','log');
	%legend(vecHandles,cellTypes,'location','best')
	fixfig;
	return
	%%
	for i=1
		%%
		continue
		%pre-allocate
		intTimescaleNum = numel(vecTimescales);
		intJitterNum = numel(vecJitter);
		matMean = nan(intJitterNum,intTimescaleNum);
		matSd = nan(intJitterNum,intTimescaleNum);
		matCV = nan(intJitterNum,intTimescaleNum);
		matCVFit_Exp = nan(intJitterNum,intTimescaleNum);
		matCVFit_Root = nan(intJitterNum,intTimescaleNum);
		
		vecSlope_Lin = nan(intJitterNum,1);
		vecR2_Lin = nan(intJitterNum,1);
		
		vecR2_Exp = nan(intJitterNum,1);
		vecHalfLife_Exp = nan(intJitterNum,1);
		vecAsymptote_Exp = nan(intJitterNum,1);
		vecScale_Exp = nan(intJitterNum,1);
		
		vecR2_Root = nan(intJitterNum,1);
		vecAsymptote_Root = nan(intJitterNum,1);
		vecScale_Root = nan(intJitterNum,1);
		vecExponent_Root = nan(intJitterNum,1);
		
		clear sAggData;
		for intJitterIdx=1:numel(vecJitter)
			%% jitter
			dblJitter = vecJitter(intJitterIdx);
			vecAllSpikeTime = cell2vec(cellSpikeTimes);
			vecRand = (rand(size(vecAllSpikeTime))-0.5)*2; %[-1 +1]
			vecAllSpikeTime = sort(vecAllSpikeTime + dblJitter*vecRand);
			
			%% get lin fit for mean/sd
			vecSlopes = nan(size(vecTimescales));
			vecR2 = nan(size(vecTimescales));
			vecMean = nan(size(vecTimescales));
			vecSd = nan(size(vecTimescales));
			intK=1;
			for intScale=1:intTimescaleNum
				vecBins = sSource.dblStartEpoch:vecTimescales(intScale):(sSource.dblStartEpoch+sSource.dblEpochDur);
				vecCounts = histcounts( vecAllSpikeTime,vecBins);
				vecMean(intScale) = mean(vecCounts);
				vecSd(intScale) = std(vecCounts);
			end
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
			vecUpper = [1e6 1e6 1];
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
			[fitobject,gof] = fit(vecX,vecY,fExp,'lower',[0 0 1e-6],'upper',[1e6 1e6 1e6],'startpoint',vecStartCoeffsExp);
			vecFitExp = fitobject(vecX);
			[dblR2_Exp,dblSS_tot,dblSS_res,dblT,dblP,dblR2_adjusted,dblR2_SE] = getR2(vecY,vecFitExp,intK);
			dblHalfLife_Exp = fitobject.c*log(2);
			dblAsymptote_Exp = fitobject.a;
			dblScale_Exp = fitobject.b;
			
			%fit root model
			[fitobject,gof] = fit(vecX,vecY,fRoot,'lower',vecLower,'upper',vecUpper,'startpoint',vecStartCoeffs);
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
			matMean(intJitterIdx,:) = vecMean;
			matSd(intJitterIdx,:) = vecSd;
			vecSlope_Lin(intJitterIdx) = dblSlope_Lin;
			vecR2_Lin(intJitterIdx) = dblR2_Lin;
			
			matCV(intJitterIdx,:) = vecSd./vecMean;
			matCVFit_Exp(intJitterIdx,:) = vecFitExp;
			matCVFit_Root(intJitterIdx,:) = vecFitRoot;
			
			vecR2_Exp(intJitterIdx) = dblR2_Exp;
			vecHalfLife_Exp(intJitterIdx) = dblHalfLife_Exp;
			vecAsymptote_Exp(intJitterIdx) = dblAsymptote_Exp;
			vecScale_Exp(intJitterIdx) = dblScale_Exp;
			
			vecR2_Root(intJitterIdx) = dblR2_Root;
			vecAsymptote_Root(intJitterIdx) = dblAsymptote_Root;
			vecScale_Root(intJitterIdx) = dblScale_Root;
			vecExponent_Root(intJitterIdx) = dblExponent_Root;
		end
		
		
		% plot
		if boolMakeFigs
			figure;maxfig;
			colormap(parula)
			vecPlotJitters = 1:intJitterNum;%[1 2 ceil(intJitterNum/2) intJitterNum];
			subplot(2,4,1)
			matCV = matSd./matMean;
			plot(vecTimescales,matCV(vecPlotJitters,:)')
			cellL=(cellfun(@(x) sprintf('%.2f',x),vec2cell(vecJitter),'UniformOutput',false));
			h=colorbar;
			vecShowJitterTicks = [1 6 11 15 20];
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
			plot(vecJitter,vecR2_Exp)
			title('exp decay r2')
			xlabel('Spike jitter (s)');
			ylabel('R^2 exp decay fit');
			
			subplot(2,4,4);
			plot(vecJitter,vecHalfLife_Exp)
			title('exp decay half-life')
			xlabel('Spike jitter (s)');
			ylabel('Half-life decay fit');
			
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
			plot(vecJitter,vecR2_Root)
			title('root r2')
			xlabel('Spike jitter (s)');
			ylabel('R^2 root fit');
			
			subplot(2,4,8);
			plot(vecJitter,vecExponent_Root)
			title('root exponent')
			xlabel('Spike jitter (s)');
			ylabel('Exponent of root fit');
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
	
	%% save agg data
	save(fullpath(strTargetDataPath,sprintf('Q2cData%s_%s.mat',strThisRec,strOnset)),...
		'sAggData');
	close all;
end
toc
