%get real data
dblBaseDur = 0.4;
dblEp1Dur = 0.2;
dblEp2Dur = 0.2;
dblMeanRateEp1 = 4;
dblMeanRateEp2 = 3.5;

dblShape = 2;
dblScale = dblMeanRateEp2/dblShape;
vecGam = random('gamma',dblShape,dblScale,[1 1000]);
dblMean = dblShape * dblScale;

vecMeanRateHz = cellfun(@numel,cellSpikes)./max(cellfun(@max,cellSpikes) - cellfun(@min,cellSpikes));
vecBaseRateHz = mean(matRespBase,2)./dblBaseDur;
vecEp1RateHz = mean(matRespEp1,2)./dblEp1Dur;
vecEp2RateHz = mean(matRespEp2,2)./dblEp2Dur;


%get neuronal activity from t=0 - t=0.4
dblBinSizeSecs = 200/1000;%25/1000;
vecBinEdges = 0:dblBinSizeSecs:0.4;
vecBinCenters = vecBinEdges(2:end)-dblBinSizeSecs/2;
intBins = numel(vecBinCenters);

%set params
intN = size(matRespBase,1);
intT = size(matRespBase,2);
vecTrialEps=zeros(1,intT);
dblTrialDur = 0.4;
dblStimDur = dblTrialDur/2;
vecTrialStarts = (dblTrialDur:dblTrialDur:dblTrialDur*intT)';
vecTrialStops = vecTrialStarts + dblStimDur;
matTrialT = cat(2,vecTrialStarts,vecTrialStops);
cellGenSpikes = cell(1,intN);
for intNeuron=1:intN
	%build baseline activity for X and Y
	dblBaseRate = vecEp1RateHz(intNeuron);
	dblPrefRate = vecEp2RateHz(intNeuron);
	dblPrefOri = 0;
	if dblBaseRate > dblPrefRate
		dblPrefOri = 1;
		dblBaseRate = vecEp2RateHz(intNeuron);
		dblPrefRate = vecEp1RateHz(intNeuron);
	end
	[vecSpikeTimes,dblPrefOri] = getGeneratedSpikingData(vecTrialEps,matTrialT,dblBaseRate,dblPrefRate,1,false,dblPrefOri);
	cellGenSpikes{intNeuron} = vecSpikeTimes;
end
%check rates
vecExpRate = (vecEp1RateHz + vecEp2RateHz)/2;
vecGenRate = (cellfun(@numel,cellGenSpikes)./intT)./dblTrialDur;

%% transform data
vecGenEp1Edges = sort(cat(1,vecTrialStarts,vecTrialStarts+dblStimDur));
vecGenEp2Edges = sort(cat(1,vecTrialStops,vecTrialStops+dblStimDur));
%calculate response matrix
matGenRespEp1 = nan(intN,intT);
matGenRespEp2 = nan(intN,intT);
%bin spikes
matGenActBin = zeros(intBins,intN,intT);
for intNeuron=1:intN
	vecEp1R = histcounts(cellGenSpikes{intNeuron},vecGenEp1Edges);
	matGenRespEp1(intNeuron,:) = vecEp1R(1:2:end);
	vecEp2R = histcounts(cellGenSpikes{intNeuron},vecGenEp2Edges);
	matGenRespEp2(intNeuron,:) = vecEp2R(1:2:end);
	[vecTrialPerSpikeEp1,vecTimePerSpike] = getSpikesInTrial(cellGenSpikes{intNeuron},vecTrialStarts);
	%bin spikes
	for intTrial=1:intTrials
		vecSpikeTimes = vecTimePerSpike(vecTrialPerSpikeEp1==intTrial);
		matGenActBin(:,intNeuron,intTrial) = histcounts(vecSpikeTimes,vecBinEdges);
	end
end

%plot
figure
subplot(2,2,1)
imagesc(matRespEp1)

subplot(2,2,2)
imagesc(matRespEp2)

subplot(2,2,3)
imagesc(matGenRespEp1)

subplot(2,2,4)
imagesc(matGenRespEp2)

	%% prep data
	matGenData = cat(1,matGenRespEp1,matGenRespEp2);
	vecCellArea = cat(1,ones(intN,1),ones(intN,1)*2);
	vecTrialStimType = ones(1,intT);
	intStimTypes = numel(unique(vecTrialStimType));

		%% get data splits
		%set parameters
		dblCutOff = 0.9;
		sParams=struct;
		intUseSize = intN;%16%floor(intCells/2);
		sParams.intSizeX = intUseSize;%floor(intCells/2);
		sParams.intSizeY = intUseSize;%floor(intCells/2);
		sParams.intResamplings = 10;
		sParams.vecCellArea = vecCellArea;
		sParams.intWithinArea = 3;%0, across; 1, within 1; 2, within 2
		sParams.vecUseStimTypes = 1;
		cellStr = {'B1','B2','B1-B2'};
		strC = cellStr{sParams.intWithinArea};
		
		%get splits
		[cellGenMatX1,cellGenMatX2,cellGenNeuronsX,cellGenMatY1,cellGenMatY2,cellGenNeuronsY,cellGenTrials12,cellGenTrials1,cellGenTrials2] = ...
			doDimDataSplitsCV(matGenData,vecTrialStimType,sParams);
		
		%% analyze
		intMaxDim = sParams.intSizeY;
		[matAvS,matAvPX,matAvPY,matAggS,matAggPX,matAggPY,matNorms] = doDimSharedPrivateAnalysisCV(matActBin,cellGenNeuronsX,cellGenNeuronsY,cellGenMatX1,cellGenMatY1,cellGenMatX2,cellGenMatY2,dblLambda,boolShuffle);
		
%% plot
figure
subplot(4,6,1)
matCovX = cov(matX2);
imagesc(matCovX,max(abs(matCovX(:)))*[-1 1]);
title('A1) Cov(X)')
fixfig;grid off;

subplot(4,6,2)
imagesc(matTransferXY,max(abs(matTransferXY(:)))*[-1 1]);
title('A2) Transfer(X->Y)')
fixfig;grid off;

subplot(4,6,7)
matCovY = cov(matY2);
imagesc(matCovY,max(abs(matCovY(:)))*[-1 1]);
title('A3) Cov(Y)')
fixfig;grid off;

subplot(4,6,8)
imagesc(matTransferYX,max(abs(matTransferYX(:)))*[-1 1]);
title('A4) Transfer(Y->X)')
fixfig;grid off;

%check magnitudes
vecXS = mean(abs(matX2 * matSubspace_Shared),1)
vecXPX = mean(abs(matX2 * matSubspace_PrivateX),1)
vecXPY = mean(abs(matX2 * matSubspace_PrivateY),1)

vecYS = mean(abs(matY2 * matSubspace_Shared),1)
vecYPX = mean(abs(matY2 * matSubspace_PrivateX),1)
vecYPY = mean(abs(matY2 * matSubspace_PrivateY),1)

subplot(2,3,2)
plot(vecXS,'m')
hold on
plot(vecXPX,'r')
plot(vecXPY,'b')
hold off
xlabel('Dimension (neuron #)');
title('B) X proj in subsp')
legend({'Shared','X-priv','Y-priv',},'Location','Best')
fixfig;

subplot(2,3,3)
plot(vecYS,'m')
hold on
plot(vecYPX,'r')
plot(vecYPY,'b')
hold off
xlabel('Dimension (neuron #)');
title('C) Y proj in subsp')
legend({'Shared','X-priv','Y-priv',},'Location','Best')
fixfig;

subplot(2,3,4)
imagesc(matSubspace_Shared,max(abs(matSubspace_Shared(:)))*[-1 1])
title('D) Shared subspace');
xlabel('Dimension (neuron #)');
ylabel('Dimension (neuron #)');
colormap(redblue)
fixfig;grid off;

subplot(2,3,5)
imagesc(matSubspace_PrivateX,max(abs(matSubspace_PrivateX(:)))*[-1 1])
title('E) Private X subspace');
fixfig;grid off;
xlabel('Dimension (neuron #)');
ylabel('Dimension (neuron #)');

subplot(2,3,6)
imagesc(matSubspace_PrivateY,max(abs(matSubspace_PrivateY(:)))*[-1 1])
title('F) Private Y subspace');
fixfig;grid off;
xlabel('Dimension (neuron #)');
ylabel('Dimension (neuron #)');

maxfig;
%export_fig('Example_SharedPrivate.tif')
%export_fig('Example_SharedPrivate.pdf')

%check magnitudes
dblXS = mean(vecXS(:))
dblXPX = mean(vecXPX(:))
dblXPY = mean(vecXPY(:))

dblYS = mean(vecYS(:))
dblYPX = mean(vecYPX(:))
dblYPY = mean(vecYPY(:))
