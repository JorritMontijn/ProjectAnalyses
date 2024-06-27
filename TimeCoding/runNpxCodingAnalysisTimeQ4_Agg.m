%% overview
%what is better model for pop activity? multivariate gaussian noise or fixed sd/mu over a
%gain-scaling axis?

%% set parameters
cellDataTypes = {'Npx','Sim','ABI','SWN'};%topo, model, allen, nora
intRunDataType = 1;
strRunStim = 'DG';%DG or NM? => superseded to WS by SWN
dblRemOnset = 0; %remove onset period in seconds; 0.125 for sim, 0.25 for npx
cellTypes = {'Real','Poiss','ShuffTid','RandTid','RandTxClass','Uniform'};%{'Real','Poiss','ShuffTid','RandTid','RandTxClass','Uniform'};
runHeaderPopTimeCoding;
boolMakeFigs = false;
vecTimescales = [1e-2 5e-2 1e-1 5e-1 1e-0]; %timescales for rapid flucts, noise corrs, tuning
vecJitter = [0 1e3];%[0 (2.^(-9:10))];
intPopSize = inf; %24 (smallest pop of all recs) or inf (uses full pop for each rec)

%% go through recs
tic
matAggTuning = [];
vecSourceRec = [];
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
	if mean(sum(matData)) < 90%90 / 50
		fprintf('Avg # of spikes per trial was %.1f for %s; skipping...\n',mean(sum(matData)),strThisRec);
		continue;
	else
		fprintf('Running %s... [%s]\n',strThisRec,getTime);
	end
	
	%get ori vars
	intTunedN = sum(indTuned);
	intRespN = size(matData,1);
	vecOri180 = mod(vecOrientation,180);
	[matRespNSR,vecStimTypes,vecUniqueDegs] = getStimulusResponses(matData,vecOri180);
	matMedianR = mean(matRespNSR,3);
	matAggTuning = cat(1,matAggTuning,matMedianR);
	vecSourceRec = cat(1,vecSourceRec,intRec*ones(intRespN,1));
	continue;
	
	%% go through types
	clear sAggData;
	vecRunTypes = 1:numel(cellTypes);
	for intType=vecRunTypes
		strType = cellTypes{intType};
		
		%% load prepped data and ifr
		sSource = load(fullpath(strDataPathT0,sprintf('T0Data_%s%s%s%s%s',strThisRec,strRunType,strRunStim,strType)));
		cellOrigSpikeTimes = sSource.cellSpikeTimes;
		vecPseudoStartT = vecOrigStimOnTime;
		
		vecTrialT = sSource.vecStimOnTime;
		vecTrialT_end = vecTrialT+1;
		vecEdges = sort(cat(2,vecTrialT,vecTrialT_end));
		vecSpikeCounts = histcounts(sSource.vecAllSpikeTime,vecEdges);
		vecSpikesPerTrial = vecSpikeCounts(1:2:end);
		
		[a,b,c,d,e,cellValsX]=getQuantiles(vecSpikesPerTrial,vecSpikesPerTrial,7);
		vecMu = cellfun(@mean,cellValsX);
		vecSd = cellfun(@std,cellValsX);
		plot(vecMu,vecSd)
		return
	end
	%% save agg data
	save(fullpath(strTargetDataPath,sprintf('Q3Data%s_%s.mat',strThisRec,strOnset)),...
		'sAggData');
end
toc

%%
vecSourceRec = val2idx(vecSourceRec);
vecS=sum(matAggTuning);
dblEquidistanceness = 1-std(vecS)/mean(vecS)
vecRecs = unique(vecSourceRec);
matStimMeans = [];
vecNperR = [];
vecComboN = [];
vecComboEquiD = [];



for i=1:numel(vecRecs)
	intRec=vecRecs(i)
	vecN = vecSourceRec==intRec;
	matThisR = matAggTuning(vecN,:);
	matStimMeans(i,:) = mean(matThisR);
	vecNperR(i) = sum(vecN);
	
	matAllCombs = nchoosek(vecRecs,i);
	for j = 1:size(matAllCombs,2)
		vecChoose = matAllCombs(:,j);
		vecN = ismember(vecSourceRec,vecChoose);
		matThisR = matAggTuning(vecN,:);
		vecStimMeans = mean(matThisR);
		vecComboN(end+1) = sum(vecN);
		vecComboEquiD(end+1) = 1 - std(vecStimMeans) ./ mean(vecStimMeans);
	end
end

vecEquiD = 1 - std(matStimMeans,[],2) ./ mean(matStimMeans,2);

%
figure
vecPlotR = matStimMeans(1,:);
vecPlotR = sum(matAggTuning);
%vecPlotR = sum(matStimMeans);
dblLim = 4000;%700;
dblExampleEquiD = 1 - std(vecPlotR) ./ mean(vecPlotR);
subplot(2,3,4)
polarplot(deg2rad(vecUniqueDegs)*2,vecPlotR);

subplot(2,3,1)
[x,y]=pol2cart(deg2rad(vecUniqueDegs)*2,vecPlotR);
matCol = circcol(numel(x));
colormap(matCol);
cline([x x(1)],[y y(1)],[1:numel(x) 1]);
hold on
[x2,y2]=pol2cart(deg2rad(vecUniqueDegs)*2,ones(size(vecUniqueDegs))*mean(vecPlotR));
plot([x2 x2(1)],[y2 y2(1)],'--','color',[0.5 0.5 0.5]);
colorbar
axis equal
ylim(dblLim*[-1 1]);
xlim(dblLim*[-1 1]);
title(sprintf('%.1f%% equidistance symmetry',100*dblExampleEquiD))


subplot(2,3,2)
scatter(vecComboN,vecComboEquiD);
[r,p]=corr(vecComboN',vecComboEquiD');
vecBinsE = 0:50:800;
vecBinsC = vecBinsE(2:end)-diff(vecBinsE(1:2))/2;
[vecCounts,vecMeans,vecSDs,cellVals,cellIDs] = makeBins(vecComboN',vecComboEquiD',vecBinsE);
indRem = isnan(vecMeans);
vecMeans(indRem) = [];
vecBinsC(indRem) = [];
hold on
plot(vecBinsC,vecMeans);

%fit logistic growth curve
	%	beta(1) = L (asymptote) [default 1]
	%	beta(2) = k (slope) [default 1]
	%	beta(3) = x0 (x-offset) [default 0]
	%	beta(4) = y0 (y-offset) [default 0]

	sOptions = optimset('Display','iter','TolX',1e-20,'TolFun',1e-20,'MaxIter',1e4,'MaxFunEvals',1e4);
vecP0 = [0.05 0.01 0 0.9];
xFit = vecBinsC;
yFit = vecMeans';
%xFit = vecComboN;
%yFit = vecComboEquiD;

vecP = lsqcurvefit(@logisticfitx0, vecP0, xFit, yFit,[0 0 0 0],[1 1e3 0 1e3],sOptions);
xPlot = 10:10:1000;
y = logisticfitx0(vecP,xPlot);
dblCeiling = vecP(1) + vecP(4);
plot(xPlot,y);
y2 = logisticfitx0(vecP0,xFit);
scatter(xFit,y2);

%ylim([0 1]);
title(sprintf('Across each rec: %.1f%% equidistance symmetry',100*mean(vecEquiD)))

subplot(2,3,3);
scatter(ones(size(vecEquiD)),vecEquiD);
ylim([0 1]);
title(sprintf('Across each rec: %.1f%% equidistance symmetry',100*mean(vecEquiD)))
hold on
errorbar(1,mean(vecEquiD),std(vecEquiD)./sqrt(numel(vecEquiD)));