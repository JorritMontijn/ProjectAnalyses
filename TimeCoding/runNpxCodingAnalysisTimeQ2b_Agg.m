
%% set parameters
cellDataTypes = {'Npx','Sim','ABI','SWN'};%topo, model, allen, nora
intRunDataType = 1;
strRunStim = 'DG';%DG or NM? => superseded to WS by SWN
cellTypes = {'Real','Poiss','ShuffTid','Shuff','PoissGain','Uniform'};
dblRemOnset = 0; %remove onset period in seconds; 0.125 for sim, 0.25 for npx
runHeaderPopTimeCoding;
boolMakeFigs = false;

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
	else
		fprintf('Running %s... [%s]\n',strThisRec,getTime);
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
		
		%%
		vecAllSpikeTime = sort(cell2vec(cellSpikeTimes));
		vecTimescales = 0.01:0.01:1;
		vecMean = nan(size(vecTimescales));
		vecSd = nan(size(vecTimescales));
		for intScale=1:numel(vecTimescales)
			vecBins = sSource.dblStartEpoch:vecTimescales(intScale):(sSource.dblStartEpoch+sSource.dblEpochDur);
			vecCounts = histcounts( vecAllSpikeTime,vecBins);
			vecMean(intScale) = mean(vecCounts);
			vecSd(intScale) = std(vecCounts);
		end
		
		%subselect neurons
		vecN = 1:intNumN;
		matMean = nan(intNumN,numel(vecTimescales));
		matSd = nan(intNumN,numel(vecTimescales));
		matTime = nan(intNumN,numel(vecTimescales));
		for intN=1:intNumN
			intN
			vecUseN = randperm(intNumN,intN);
			vecAllSpikeTime = sort(cell2vec(cellSpikeTimes(vecUseN)));
			for intScale=1:numel(vecTimescales)
				vecBins = sSource.dblStartEpoch:vecTimescales(intScale):(sSource.dblStartEpoch+sSource.dblEpochDur);
				vecCounts = histcounts( vecAllSpikeTime,vecBins);
				matMean(intN,intScale) = mean(vecCounts);
				matSd(intN,intScale) = std(vecCounts);
				matTime(intN,intScale) = vecTimescales(intScale);
			end
		end
		
		%single neurons
		vecN = 1:intNumN;
		matMeanSingle = nan(intNumN,numel(vecTimescales));
		matSdSingle = nan(intNumN,numel(vecTimescales));
		matTimeSingle = nan(intNumN,numel(vecTimescales));
		for intN=1:intNumN
			intN
			vecAllSpikeTime = sort(cell2vec(cellSpikeTimes(intN)));
			for intScale=1:numel(vecTimescales)
				vecBins = sSource.dblStartEpoch:vecTimescales(intScale):(sSource.dblStartEpoch+sSource.dblEpochDur);
				vecCounts = histcounts( vecAllSpikeTime,vecBins);
				matMeanSingle(intN,intScale) = mean(vecCounts);
				matSdSingle(intN,intScale) = std(vecCounts);
				matTimeSingle(intN,intScale) = vecTimescales(intScale);
			end
		end
		
		%% plot
		figure
		subplot(3,4,1)
		cline(vecMean,vecSd,vecTimescales)
		xlabel('Mean of spike counts over bins');
		ylabel('Sd of spike counts over bins');
		title(sprintf('Color=timescale; %s - %s',strRec,strType),'interpreter','none');
		colorbar('Ticks',[0 0.5 1],'Ticklabels',(vecTimescales([1 ceil(numel(vecTimescales)/2) end])));
		
		subplot(3,4,2);
		hold on
		scatter(matMean(:),matSd(:),[],matTime(:),'.');
		hold off
		xlabel('Mean of spike counts over bins');
		ylabel('Sd of spike counts over bins');
		title(sprintf('Neuron subselection for pop size 1-%d',intNumN));
		colorbar;
		
		vecEdgesMu = 0:2:700;
		vecEdgesSd = 0:2:180;
		matDensity = makeBins2(matMean(:),matSd(:),ones(size(matTime(:))),vecEdgesMu,vecEdgesSd);
		matNormDens = matDensity./sum(matDensity,1);
		subplot(3,4,3)
		imagesc(vecEdgesMu,vecEdgesSd,matNormDens);
		%set(gca,'xscale','log');
		xlabel('Mean of spike counts over bins');
		ylabel('Sd of spike counts over bins');
		axis xy;
		colormap(gca,redwhite);
		colorbar;
		title('Density normalized per mu-bin');
		
		%fit sd/mean
		vecSlopes = nan(1,intNumN);
		vecR2 = nan(1,intNumN);
		intK=1;
		for intN=1:intNumN
			vecX = matMean(intN,:)';
			vecY = matSd(intN,:)';
			dblSlope = ((vecX' * vecX) \ vecX') * vecY;
			vecFitY = vecX*dblSlope;
			[dblR2,dblSS_tot,dblSS_res,dblT,dblP,dblR2_adjusted,dblR2_SE] = getR2(vecY,vecFitY,intK);
			
			vecSlopes(intN) = dblSlope;
			vecR2(intN) = dblR2;
		end
		
		%plot slope of sd/mean as function of # of neurons
		subplot(6,4,4);
		plot(1:intNumN,vecSlopes);
		ylabel('Slope of sd/mu (CV)');
		
		subplot(6,4,8);
		plot(1:intNumN,vecR2);
		xlabel('Population size (# of neurons)');
		ylabel('Linearity (R^2 of sd/mu)');
		
		%% single neurons
		subplot(3,4,5)
		cline(matMeanSingle(1,:),matSdSingle(1,:),matMeanSingle(1,:))
		xlabel('Mean of spike counts over bins');
		ylabel('Sd of spike counts over bins');
		title('Color=timescale');
		colorbar('Ticks',[0 0.5 1],'Ticklabels',(vecTimescales([1 ceil(numel(vecTimescales)/2) end])));
		
		subplot(3,4,6);
		hold on
		scatter(matMeanSingle(:),matSdSingle(:),[],matTimeSingle(:),'.');
		hold off
		xlabel('Mean of spike counts over bins');
		ylabel('Sd of spike counts over bins');
		title(sprintf('%d single neurons',intNumN));
		colorbar;
		
		vecEdgesMu = 0:1:50;
		vecEdgesSd = 0:1:20;
		matDensity = makeBins2(matMeanSingle(:),matSdSingle(:),ones(size(matTimeSingle(:))),vecEdgesMu,vecEdgesSd);
		matNormDens = matDensity./sum(matDensity,1);
		subplot(3,4,7)
		imagesc(vecEdgesMu,vecEdgesSd,matNormDens);
		%set(gca,'xscale','log');
		xlabel('Mean of spike counts over bins');
		ylabel('Sd of spike counts over bins');
		axis xy;
		colormap(gca,redwhite);
		colorbar;
		title('Density normalized per mu-bin');
		
		%fit sd/mean
		vecSlopes_Single = nan(1,intNumN);
		vecR2_Single = nan(1,intNumN);
		intK=1;
		for intN=1:intNumN
			vecX = matMeanSingle(intN,:)';
			vecY = matSdSingle(intN,:)';
			dblSlope = ((vecX' * vecX) \ vecX') * vecY;
			vecFitY = vecX*dblSlope;
			[dblR2,dblSS_tot,dblSS_res,dblT,dblP,dblR2_adjusted,dblR2_SE] = getR2(vecY,vecFitY,intK);
			
			vecSlopes_Single(intN) = dblSlope;
			vecR2_Single(intN) = dblR2;
		end
		
		%plot slope of sd/mean as function of # of neurons
		subplot(6,4,12);
		plot(1:intNumN,vecSlopes_Single);
		ylabel('Slope of sd/mu (CV)');
		
		subplot(6,4,16);
		plot(1:intNumN,vecR2_Single);
		xlabel('Single neuron ID');
		ylabel('Linearity (R^2 of sd/mu)');
		
		
		%% timescales
		subplot(3,4,9)
		matCV = matSd./matMean;
		cline(matTime(end,:),matCV(end,:),matMean(end,:))
		xlabel('Timescale (bin width in s)');
		ylabel('CV (sd/mu) of spike counts over bins');
		title('Color=mean spike count');
		colorbar('Ticks',[0 0.5 1],'Ticklabels',roundi(matMean(end,[1 ceil(numel(vecTimescales)/2) end]),1));
		
		subplot(3,4,10);
		hold on
		scatter(matTime(:),matCV(:),[],matMean(:),'.');
		hold off
		xlabel('Timescale (bin width in s)');
		ylabel('CV (sd/mu) of spike counts over bins');
		title('CV/timescale for all pop sizes');
		colorbar;
		
		vecEdgesTime = [vecTimescales(1)-median(diff(vecTimescales)/2) vecTimescales(2:end)-diff(vecTimescales)/2 vecTimescales(end)+median(diff(vecTimescales)/2)];
		vecEdgesCV = 0:0.1:2;
		matDensity = makeBins2(matTime(:),matCV(:),ones(size(matTime(:))),vecEdgesTime,vecEdgesCV);
		matNormDens = matDensity./sum(matDensity,1);
		subplot(3,4,11)
		imagesc(vecEdgesTime,vecEdgesCV,matNormDens);
		%set(gca,'xscale','log');
		xlabel('Timescale (bin width in s)');
		ylabel('Sd of spike counts over bins');
		axis xy;
		colormap(gca,redwhite);
		colorbar;
		title('Density normalized per t-bin');
		
		
		g = fittype('a+b*exp(-x/c)',...
			'dependent',{'y'},'independent',{'x'},...
			'coefficients',{'a','b','c'});
		warning('off','curvefit:fit:noStartPoint');
		%fit sd/mean
		vecSlopes_Time = nan(1,intNumN);
		vecR2_Time = nan(1,intNumN);
		vecLambdaExp_Time = nan(1,intNumN);
		vecR2Exp_Time = nan(1,intNumN);
		intK=1;
		for intN=1:intNumN
			%fit linear model
			vecX = vecTimescales';
			vecY = matCV(intN,:)';
			mdl = fitlm(vecX,vecY);
			
			dblSlope = mdl.Coefficients.Estimate(2);
			vecFitY =  mdl.Fitted;
			[dblR2,dblSS_tot,dblSS_res,dblT,dblP,dblR2_adjusted,dblR2_SE] = getR2(vecY,vecFitY,intK);
			vecSlopes_Time(intN) = dblSlope;
			vecR2_Time(intN) = dblR2;
			
			%fit exp model
			[fitobject,gof] = fit(vecX,vecY,g,'lower',[0 0 1e-6],'upper',[1e6 1e6 1e6]);
			vecFitExp = fitobject(vecX);
			[dblR2,dblSS_tot,dblSS_res,dblT,dblP,dblR2_adjusted,dblR2_SE] = getR2(vecY,vecFitExp,intK);
			vecLambdaExp_Time(intN) = fitobject.c*log(2);
			vecR2Exp_Time(intN) = dblR2;
			
		end
		
		%plot slope of sd/mean as function of # of neurons
		subplot(6,4,20);
		plot(1:intNumN,vecR2_Time);
		ylabel('CV-linearity with timescale (R^2 of mu/sd)');
		
		subplot(6,4,24);
		plot(1:intNumN,vecR2Exp_Time);
		xlabel('Population size (# of neurons)');
		ylabel('CV-exponentiality with timescale (R^2 of mu/sd)');
		
		%% save figure
		export_fig(fullpath(strFigurePathSR,sprintf('Q2bA_TimescaleCV_%s_%s_%s.tif',strThisRec,strType,strOnset)));
		export_fig(fullpath(strFigurePathSR,sprintf('Q2bA_TimescaleCV_%s_%s_%s.pdf',strThisRec,strType,strOnset)));
		
		%% save data
		sData = struct;
		sData.strRec = strRec;
		sData.strType = strType;
		sData.dblRemOnset = dblRemOnset;
		sData.strThisRec = strThisRec;
		sData.vecTimescales = vecTimescales;
		sData.vecSlopes = vecSlopes;
		sData.vecR2 = vecR2;
		sData.vecSlopes_Single = vecSlopes_Single;
		sData.vecR2_Single = vecR2_Single;
		sData.vecSlopes_Time = vecSlopes_Time;
		sData.vecR2_Time = vecR2_Time;
		sData.vecLambdaExp_Time = vecLambdaExp_Time;
		sData.vecR2Exp_Time = vecR2Exp_Time;
		
		%add to superstructure
		sAggData(intType) = sData;
		
	end
	close all;
	%% save agg data
	save(fullpath(strTargetDataPath,sprintf('Q2bData%s_%s.mat',strThisRec,strOnset)),...
		'sAggData');
end
toc
