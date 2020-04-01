%% initialize
%close all;
clearvars;
boolLoad = false;
boolSaveFigs = true;
if ~exist('matModelResp','var')
	clearvars -except bool*;%close all;
	intLoadSim = 13; %% SET WHICH SIM TO LOAD HERE; 41 vs 42; 2 vs. 5
	boolLoad = true;
end
if intLoadSim == 11 && boolLoad
	strSimulation = 'xAreaDistributed_OriDrift18Att0_2017-11-16';
elseif intLoadSim == 12 && boolLoad
	strSimulation = 'xAreaDistributed_OriDrift18Att1_2017-11-28';
elseif intLoadSim == 13 && boolLoad
		strSimulation = 'xAreaDistributed_OriDrift18Att0Simple_2017-11-30';
		
elseif intLoadSim == 21 && boolLoad
	strSimulation = 'xAreaDistributed_OriDrift2Att0_2017-11-16'; %new connectivity
elseif intLoadSim == 22 && boolLoad
	strSimulation = 'xAreaDistributed_OriDrift2Att1_2017-11-16'; %18 oris
	
end

%% RUN: #header
strPath = 'A:\SimAggregates\';
strFigDir = 'D:\Data\Results\Block0\';
if boolLoad
	runModelHeader;
end

%% get data and define binning
intAnalType = 2;
dblTrialDur = unique(roundi(vecStimStopSecs - vecStimStartSecs,3));
matModelRespP = double(matModelResp)/dblTrialDur;
intNeurons = size(matModelRespP,1);
intRawSizeT = numel(vecOverallT);
dblEndT = max(vecOverallT);
dblBinSize = 10/1000;%0.5/3; %50 ms bins
intBinsT = floor(dblEndT/dblBinSize);

vecBinsTime = linspace(0,dblEndT,intBinsT+1);
vecBinsStimStart = round(vecStimStartSecs/dblBinSize);
vecBinsStimStop = round(vecStimStopSecs/dblBinSize);
hTic = tic;
if ~exist('matSpikeCounts2D','var') %&& intAnalType == 1
	if exist([strPath 'Simulation_' strSimulation '_prepro3.mat'],'file');
		load([strPath 'Simulation_' strSimulation '_prepro3.mat']);
	else
		matSpikeCounts2D = zeros(intNeurons,intBinsT,'uint8');
		for intNeuron=1:intNeurons
			matSpikeCounts2D(intNeuron,:) = histcounts(cellSpikeTimesCortex{intNeuron},vecBinsTime);
			%if toc(hTic) > 5 || intNeuron==1
			%	hTic = tic;
			if mod(intNeuron,100) == 0
				fprintf('Now at neuron %d/%d [%s]\n',intNeuron,intNeurons,getTime);
			end
		end

		%save
		if intBinsT > 10000
			fprintf('Processing completed; saving prepro data [%s]\n',getTime);
			save([strPath 'Simulation_' strSimulation '_prepro3.mat'],'matSpikeCounts2D');
		end
	end
end
if ~isa(matSpikeCounts2D,'double'),matSpikeCounts2D=double(matSpikeCounts2D);end
vecBinsTime = vecBinsTime(1:(end-1)); %left edge bins

%% build selection vectors and get stimulus responses
if numel(unique(vecTrialStimType)) > 1
	if ~exist('vecPrefOri','var')
		intColumns = 252; %252
		vecColumnOris = 0:pi/intColumns:(pi-pi/intColumns);
		vecPrefOri = sort(reshape(vecColumnOris' * ones(5,1)',[intColumns*5 1]))';
	end
	intStim2 = numel(unique(vecTrialStimType));
	intStim1 = intStim2-1;
	intTrials = length(vecTrialStimType);
	intStimTypes = numel(unique(vecTrialStimType));
	indSelectPyramidsV1 = vecCellTypes(:)==1 & vecCellArea(:)==1;
	indSelectInterneuronsV1 = vecCellTypes(:)==2 & vecCellArea(:)==1;
	
	indSelectPyramidsV2 = vecCellTypes(:)==1 & vecCellArea(:)==2;
	indSelectInterneuronsV2 = vecCellTypes(:)==2 & vecCellArea(:)==2;
	
	if intAnalType == 1
		cellStimTimePoints = cellfill(false(1,intBinsT),[1 intStimTypes]);
		for intTrial=1:intTrials
			cellStimTimePoints{vecTrialStimType(intTrial)}((vecBinsStimStart(intTrial)+1):vecBinsStimStop(intTrial)) = true;
		end
		matAct = matSpikeCounts2D/dblBinSize;
		matActS1 = matSpikeCounts2D(:,cellStimTimePoints{intStim1})/dblBinSize;
		matActS2 = matSpikeCounts2D(:,cellStimTimePoints{intStim2})/dblBinSize;
	else
		matAct = matModelRespP;
		matActS1 = matAct(:,vecTrialStimType==1);
		matActS2 = matAct(:,vecTrialStimType==2);
	end
	matActPyrS1 = matActS1(indSelectPyramidsV1,:);
	matActPyrS2 = matActS2(indSelectPyramidsV1,:);
	vecPrefAngles=-vecPrefOri+pi;
	dblOriS1 = rad2ang(circ_mean(2*vecPrefAngles(indSelectPyramidsV1)',mean(matActPyrS1,2)))/2;
	dblOriS2 = rad2ang(circ_mean(2*vecPrefAngles(indSelectPyramidsV1)',mean(matActPyrS2,2)))/2;
	
	%%
	figure
	subplot(2,3,1)
	scatter(vecPrefAngles(indSelectPyramidsV1)-0.005,xmean(matActPyrS1,2)','b.');
	hold on
	scatter(vecPrefAngles(indSelectPyramidsV1)+0.005,xmean(matActPyrS2,2)','r.');
	hold off
	ylim([0 max(get(gca,'ylim'))]);
	xlim([0 pi+0.01]);
	vecAngles = 0:(0.125*pi):pi;
	set(gca,'xtick',vecAngles,'xticklabel',rad2ang(vecAngles));
	xlabel('Preferred orientation of V1 pyramid (degrees)');
	ylabel('Mean single neuron spiking rate (Hz)');
	title(sprintf('Activity of V1 Pyramids; circular mean (deg): S1 (blue),%.3f; S2 (red),%.3f',dblOriS1,dblOriS2));
	fixfig
	
	%
	matActPyrV1S1 = matActS1(indSelectPyramidsV1,:);
	matActPyrV1S2 = matActS2(indSelectPyramidsV1,:);
	vecActV1S1 = mean(matActPyrV1S1,2);
	vecActV1S2 = mean(matActPyrV1S2,2);
	vecDiffV1 = vecActV1S1-vecActV1S2;
	[a,vecReorder]=sort(vecDiffV1,'descend');
	indPrefS2 = vecDiffV1<0;
	vecActV1S1R = vecActV1S1(vecReorder);
	vecActV1S2R = vecActV1S2(vecReorder);
	indPrefS2R = indPrefS2(vecReorder);
	
	vecPlotX1 = 1:sum(~indPrefS2R);
	vecPlotX2 = (1:sum(indPrefS2R))+sum(~indPrefS2R);
	
	subplot(2,3,2)
	line([vecPlotX1;vecPlotX1],[vecActV1S1R(~indPrefS2R) vecActV1S2R(~indPrefS2R)]','Color',[0 0 1])
	hold on
	line([vecPlotX2;vecPlotX2],[vecActV1S1R(indPrefS2R) vecActV1S2R(indPrefS2R)]','Color',[1 0 0])
	scatter(1:sum(indSelectPyramidsV1),vecActV1S1R,'b.');
	scatter(1:sum(indSelectPyramidsV1),vecActV1S2R,'r.');
	
	hold off
	xlim([0 sum(indSelectPyramidsV1)]);
	xlabel('V1 Pyramid sorted by stimulus preference');
	ylabel('Mean single neuron spiking rate (Hz)');
	title(sprintf('Activity of V1 Pyramids'));
	fixfig
	
	%%
	vecBinsPlot = 0:5:70;
	vecPlotLocs = vecBinsPlot(2:end)-(vecBinsPlot(2)-vecBinsPlot(1))/2;
	subplot(2,3,3)
	vecActHzPerPyrV1 = mean(matAct(indSelectPyramidsV1,:),2);
	vecCountsPyrV1 = histcounts(vecActHzPerPyrV1,vecBinsPlot);
	vecActHzPerIntV1 = mean(matAct(indSelectInterneuronsV1,:),2);
	vecCountsIntV1 = histcounts(vecActHzPerIntV1,vecBinsPlot);
	
	
	plot(vecPlotLocs,vecCountsPyrV1,'b-');
	hold on;
	plot(vecPlotLocs,vecCountsIntV1,'b--');
	hold off
	title('Distro mean act per neuron, V1')
	legend({'Pyramids','Interneurons'});
	xlabel('Spiking rate (Hz)');
	ylabel('Number of neurons');
	fixfig;
	
	%% V2
	if sum(vecCellArea==2) > 0
		indSelectNeuronsV2 = vecCellTypes(:)==1 & vecCellArea(:)==2;
		matActPyrV2S1 = matActS1(indSelectNeuronsV2,:);
		matActPyrV2S2 = matActS2(indSelectNeuronsV2,:);
		vecActV2S1 = mean(matActPyrV2S1,2);
		vecActV2S2 = mean(matActPyrV2S2,2);
		vecDiffV2 = vecActV2S1-vecActV2S2;
		[a,vecReorder]=sort(vecDiffV2,'descend');
		indPrefS2 = vecDiffV2<0;
		vecActV2S1R = vecActV2S1(vecReorder);
		vecActV2S2R = vecActV2S2(vecReorder);
		indPrefS2R = indPrefS2(vecReorder);
		
		vecPlotX1 = 1:sum(~indPrefS2R);
		vecPlotX2 = (1:sum(indPrefS2R))+sum(~indPrefS2R);
		subplot(2,3,5)
		line([vecPlotX1;vecPlotX1],[vecActV2S1R(~indPrefS2R) vecActV2S2R(~indPrefS2R)]','Color',[0 0 1])
		hold on
		line([vecPlotX2;vecPlotX2],[vecActV2S1R(indPrefS2R) vecActV2S2R(indPrefS2R)]','Color',[1 0 0])
		scatter(1:sum(indSelectNeuronsV2),vecActV2S1R,'b.');
		scatter(1:sum(indSelectNeuronsV2),vecActV2S2R,'r.');
		
		hold off
		xlim([0 sum(indSelectNeuronsV2)]);
		xlabel('V2 Pyramid sorted by stimulus preference');
		ylabel('Mean single neuron spiking rate (Hz)');
		title(sprintf('Activity of V2 Pyramids'));
		fixfig
		
		subplot(2,3,6)
		vecActHzPerPyrV2 = mean(matAct(indSelectPyramidsV2,:),2);
		vecCountsPyrV2 = histcounts(vecActHzPerPyrV2,vecBinsPlot);
		vecActHzPerIntV2 = mean(matAct(indSelectInterneuronsV2,:),2);
		vecCountsIntV2 = histcounts(vecActHzPerIntV2,vecBinsPlot);
		
		
		plot(vecPlotLocs,vecCountsPyrV2,'b-');
		hold on;
		plot(vecPlotLocs,vecCountsIntV2,'b--');
		hold off
		title('Distro mean act per neuron, V2')
		legend({'Pyramids','Interneurons'});
		xlabel('Spiking rate (Hz)');
		ylabel('Number of neurons');
		fixfig;
	end
	%full screen
	drawnow;
	jFig = get(handle(gcf), 'JavaFrame');
	jFig.setMaximized(true);
	figure(gcf);
	drawnow;
	
	%save figure
	if boolSaveFigs
		export_fig([strFigDir 'PopulationResponses' strSimulation '.tif']);
		export_fig([strFigDir 'PopulationResponses' strSimulation '.pdf']);
	end
end

%% Fano factors
indSelectNeurons = vecCellArea==1;
matActThis = matAct(indSelectNeurons,vecTrialStimType==2);
%matActThis = matAct(indSelectNeurons,:);
dblMin  = 0;%floor(min([log(xvar(matActThis,2));log(xmean(matActThis,2))]));
dblMax  = 4;%ceil(max([log(xvar(matActThis,2));log(xmean(matActThis,2))]));

vecLimX = [dblMin dblMax];
vecLimY = [dblMin dblMax];
dblBin = 0.25;
vecBins = 0:dblBin:4;
vecBinsPlot = (dblBin/2):dblBin:4;

figure
subplot(2,3,1)
vecFanos = var(matActThis,0,2)./mean(matActThis,2);
vecY = histcounts(vecFanos,vecBins);
plot(vecBinsPlot,vecY);
xlim([0 4]);
xlabel('Fano Factor')
ylabel('Number of neurons')
title('All neurons')
fixfig;

subplot(2,3,2)
vecY = histcounts(vecFanos(vecCellTypes(indSelectNeurons)==1),vecBins);
plot(vecBinsPlot,vecY);
xlim([0 4]);
xlabel('Fano Factor')
ylabel('Number of pyramidal neurons')
title('Only excitatory')
fixfig;

subplot(2,3,3)
vecY = histcounts(vecFanos(vecCellTypes(indSelectNeurons)==2),vecBins);
plot(vecBinsPlot,vecY);
xlim([0 4]);
xlabel('Fano Factor')
ylabel('Number of inhibitory neurons')
title('Only inhibitory')
fixfig;


subplot(2,3,4)
scatter(log(xmean(matActThis,2)),log(xvar(matActThis,2)));
title('All neurons; Each point is one neuron')
xlabel('log(mean)')
ylabel('log(var)')
xlim(vecLimX);ylim(vecLimY);fixfig;


subplot(2,3,5)
scatter(log(xmean(matActThis(vecCellTypes(indSelectNeurons)==1,:),2)),log(xvar(matActThis(vecCellTypes(indSelectNeurons)==1,:),2)));
title('Only excitatory; Each point is one neuron')
xlabel('log(mean)')
ylabel('log(var)')
xlim(vecLimX);ylim(vecLimY);fixfig;


subplot(2,3,6)
scatter(log(xmean(matActThis(vecCellTypes(indSelectNeurons)==2,:),2)),log(xvar(matActThis(vecCellTypes(indSelectNeurons)==2,:),2)));
title('Only inhibitory; Each point is one neuron')
xlabel('log(mean)')
ylabel('log(var)')
xlim(vecLimX);ylim(vecLimY);fixfig;


%full screen
drawnow;
jFig = get(handle(gcf), 'JavaFrame');
jFig.setMaximized(true);
figure(gcf);
drawnow;

if boolSaveFigs
	export_fig([strFigDir 'FanoFactors' strSimulation '.tif']);
	export_fig([strFigDir 'FanoFactors' strSimulation '.pdf']);
end
%% factor analysis to check eigenvalues of first 5 modes of shared covariance
figure
%full screen
drawnow;
jFig = get(handle(gcf), 'JavaFrame');
jFig.setMaximized(true);
figure(gcf);
drawnow;

intUseNs = 50; %50 by Chengcheng
vecRandNeuronsV1 = randperm(intCellsV1,intUseNs);
vecRandNeuronsV2 = randperm(intCellsV2,intUseNs)+intCellsV1;
matDataFA = matModelRespP(vecRandNeuronsV2,vecTrialStimType==1)';
[matL, mu, vecPsi] = factorAnal(matDataFA,5);
matS_shared = matL*matL';
matPsi = diag(vecPsi);
matCov = cov(matDataFA);

subplot(2,3,1)
imagesc(matCov);
colorbar
title([sprintf('Covariance') ' \Sigma']);
fixfig;grid off;

subplot(2,3,2)
imagesc(matS_shared);
colorbar
title([sprintf('Shared cov, LL^T')])
fixfig;grid off;

subplot(2,3,3)
imagesc(matS_shared + matPsi);
colorbar
title([sprintf('Shared+Indep cov, LL^T') ' + \Psi'])
fixfig;grid off;

subplot(2,3,4)
[vecV,vecD] = eig(matCov);
vecD = diag(vecD);
vecLargestLambdas = findmax(vecD,5);
plot(vecLargestLambdas)
title('Eigen decomposition; first 5 el''s');ylabel('Eigenvalue');xlabel('Component #');
fixfig;

subplot(2,3,5)
[vecV,vecD] = eig(matS_shared);
vecD = diag(vecD);
vecLargestLambdas = findmax(vecD,5);
plot(vecLargestLambdas)
title('Eigen decomposition; first 5 el''s');ylabel('Eigenvalue');xlabel('Component #');
fixfig;

subplot(2,3,6)
[vecV,vecD] = eig(matS_shared + matPsi);
vecD = diag(vecD);
vecLargestLambdas = findmax(vecD,5);
plot(vecLargestLambdas)
title('Eigen decomposition; first 5 el''s');ylabel('Eigenvalue');xlabel('Component #');
fixfig;

if boolSaveFigs
	drawnow;
	export_fig([strFigDir 'EigenDecompCov' strSimulation '.tif']);
	export_fig([strFigDir 'EigenDecompCov' strSimulation '.pdf']);
end

return

%% plot activity over time
vecBinsTime = linspace(0,dblEndT,intBinsT+1);
vecBinsStimStart = round(vecStimStartSecs/dblBinSize);
vecBinsStimStop = round(vecStimStopSecs/dblBinSize);
vecBinsStim = zeros(size(vecBinsTime));
for intStim=1:intStimTypes
	vecTheseStims = find(vecTrialStimType==intStim);
	for intTrialOfStim=vecTheseStims
		intStart = vecBinsStimStart(intTrialOfStim);
		intStop = vecBinsStimStop(intTrialOfStim);
		vecBinsStim(intStart+1:intStop) = intStim;
	end
end

%matActV1 = matModelRespP(1:intCellsV1,:);
matActV1 = matSpikeCounts2D(1:intCellsV1,:);

vecPrefOri
vecPrefSF
vPPsi = (vecPrefPsi/max(vecPrefPsi))*nanmax(unique(diff(vecPrefOri)));
%vecPrefRF_X
%vecPrefRF_Y
vecTypesV1 = vecCellTypes(1:intCellsV1)';
vecLimX = [0 3.2];
vecLimY = [-2 2];
for i=1:size(matActV1,2)
	intStim = vecBinsStim(i);
	
	indSpiking = matActV1(:,i) & vecTypesV1==1;
	vecSpikes = matActV1(indSpiking,i);
	vecOriPs = vecPrefOri(indSpiking);
	vecPsiPs = vPPsi(indSpiking);
	
	vecSFPs = vecPrefSF(indSpiking);
	vecXPs = vecPrefRF_X(indSpiking);
	vecYPs = vecPrefRF_Y(indSpiking);
	
scatter(vecOriPs+vecPsiPs,vecXPs+vecYPs,[],vecSpikes,'*');
xlim(vecLimX);
ylim(vecLimY);
%scatter(vecXPs,vecYPs,vecSpikes,vecOriPs);
title(sprintf('t=%d/%d (%.3fs); stim %d',i,size(matActV1,2),vecBinsTime(i),intStim));
drawnow;
pause(0.05);
end

%% tuning curves
vecMeanOriRespAll = nan(1,intStimTypes);
vecSDOriRespAll =  nan(1,intStimTypes);

vecMeanOriResp = nan(intNeurons,intStimTypes);
vecSDOriResp = nan(intNeurons,intStimTypes);
for intStimType=1:intStimTypes
	vecMeanOriRespAll(intStimType) = xmean(xmean(matModelRespP(:,vecTrialStimType==intStimType),1),2);
	vecSDOriRespAll(intStimType) = xstd(xmean(matModelRespP(:,vecTrialStimType==intStimType),1),2);
	
	vecMeanOriResp(:,intStimType) = xmean(matModelRespP(:,vecTrialStimType==intStimType),2);
	vecSDOriResp(:,intStimType) = xstd(matModelRespP(:,vecTrialStimType==intStimType),2);
end
vecStimRads = ang2rad(vecStimTypeOris);
%%
for intNeuron=1:intNeurons
	clf;
	hold on
	%for  intStimType=1:intStimTypes
	%	indTrials = vecTrialStimType==intStimType;
	%scatter(3*rand(1,sum(indTrials))+ones(1,sum(indTrials))*vecStimTypeOris(intStimType),rand(1,sum(indTrials))+matModelRespP(intNeuron,indTrials));
	%end
	errorbar(vecStimTypeOris,vecMeanOriResp(intNeuron,:),vecSDOriResp(intNeuron,:));
	title(sprintf('Neuron %d',intNeuron));
	xlabel('Orientation');
	ylabel('Spiking activity (Hz)');
	drawnow;
	ylim([0 max(get(gca,'ylim'))]);
	drawnow;
	%pause(0.1);
	
end
%
%%
for t=1:1400
	imagesc(sStimInputs.cellLGN_ON{2}(:,:,t));
	title(sprintf('t=%d',t))
	drawnow;
	pause(0.01)
end

%%
%full screen
drawnow;
jFig = get(handle(gcf), 'JavaFrame');
jFig.setMaximized(true);
figure(gcf);
drawnow;

if boolSaveFigs
	%export_fig([strFigDir 'FanoFactors' strSimulation '.tif']);
	%export_fig([strFigDir 'FanoFactors' strSimulation '.pdf']);
end

return

%% analyze connectivity
fprintf('Building connectivity matrix... [%s]\n',getTime);
vecCellTypes = vecCellTypes';
vecSynWeight = ones(size(vecSynWeight));
vecSign = vecSynExcInh;%vecCortSynType;
vecSign(vecSign==2) = -1;
matConn2D = zeros(intNeurons,intNeurons);
matConn2D = getFillGrid(matConn2D,matSynFromTo(:,1),matSynFromTo(:,2),vecSynConductance.*vecSign.*vecSynWeight);

vecTypesV1 = vecCellTypes(1:intCellsV1);
[d,vecReorderV1] = sort(vecTypesV1);
vecIdxV2 = (intCellsV1+1):intNeurons;
matConn2D = matConn2D([vecReorderV1 vecIdxV2],[vecReorderV1 vecIdxV2]);


intV1Pyrs = sum(d==1);
matConn2DV2 = matConn2D(1:intV1Pyrs,vecIdxV2);
vecMeanAngleInputs = nan(1,intCellsV2);
for intV2=1:intCellsV2
	vecProjFromV1 = find(matConn2DV2(:,intV2)>0);
	vecRad = (vecProjFromV1/intV1Pyrs)*2*pi;
	vecMeanAngleInputs(intV2) = mod(circ_mean(vecRad),2*pi);
end

[d,vecReorderV2] = sort(vecMeanAngleInputs);
vecReorderV2 = vecReorderV2 + intCellsV1;
vecIdxV1 = 1:intCellsV1;
vecTypesV2 = vecCellTypes(vecReorderV2);
matConn2D = matConn2D([vecIdxV1 vecReorderV2],[vecIdxV1 vecReorderV2]);

[d,vecReReorderV2] = sort(vecTypesV2);
vecReReorderV2 = vecReReorderV2 + intCellsV1;
matConn2D = matConn2D([vecIdxV1 vecReReorderV2],[vecIdxV1 vecReReorderV2]);

%%
matConn2D = matConn2D(1:1200,1:1200);
%%
matConn2D = matConn2D(1:intV1Pyrs,1:intV1Pyrs);

%%
% pca
intRank = rank(matConn2D);
[coeff,score,latent,tsquared,explained,mu] = pca(matConn2D);
vecR2_PCA = explained/100;

%% plot
figure
%{
%activity matrix
subplot(2,3,6)
%imagesc(vecBinsTime,1:intNeurons,-(matSpikeCounts2D>0));colormap(grey);freezeColors;
hold on
plot([min(vecBinsTime) max(vecBinsTime)],[1 1]*sum(vecCellArea==1),'b--');
hold off
ylabel('Neuron #');
xlabel('Time (s)');
title('Spike raster; dotted line is V1/V2 border');
fixfig
grid off;
ylim([0 intNeurons]);
%}
subplot(2,3,[1 2 4 5])
%surf(matConn2D,'EdgeColor','none');
%imagesc(matConn2D);
imagesc(matConn2D,[-0.6 0.6+eps]);
colormap(redblue);
xlabel('Post-synaptic neuron');
ylabel('Pre-synaptic neuron');
title(['Connectivity matrix; rank is ' num2str(intRank)]);
cH = colorbar;
cH.Label.String = 'Conductance (G)';
xlim([0 intNeurons]);
ylim([0 intNeurons]);
fixfig;
grid off;

%pca
subplot(2,3,3)
plot(cumsum(vecR2_PCA));
xlabel('Principal component');
ylabel('Cumulative explained variance');
title('PCA on connectivity matrix');
xlim([0 intNeurons])
fixfig

%full screen
drawnow;
jFig = get(handle(gcf), 'JavaFrame');
jFig.setMaximized(true);
figure(gcf);
drawnow;

strSimulation= 'New'
%save figure
if boolSaveFigs
	export_fig([strFigDir 'ConnectivityMatrix' strSimulation '.tif']);
	export_fig([strFigDir 'ConnectivityMatrix' strSimulation '.pdf']);
end

%% V1 fields
cellType = {'Pyramid','Interneuron'};
figure
vecRand = randperm(size(matPrefGabors,3),6);
for i=1:6
	subplot(2,3,i)
	intRand = vecRand(i);
	dblMax = max(max(abs(matPrefGabors(:,:,intRand))));
	imagesc(matPrefGabors(:,:,intRand),dblMax*[-1 1]);
	title(sprintf('V1 neuron %d (#%d); %s',intRand,intRand,cellType{vecCellTypes(intRand)}));
end
colormap(redblue)


%% V2 fields
cellType = {'Pyramid','Interneuron'};
figure
vecRand = randperm(size(matFieldsV2,3),6);
for i=1:6
	subplot(2,3,i)
	intRand = vecRand(i);
	dblMax = max(max(abs(matFieldsV2(:,:,intRand))));
	imagesc(matFieldsV2(:,:,intRand),dblMax*[-1 1]);
	title(sprintf('V2 neuron %d (#%d); %s',intRand,intRand+intCellsV1,cellType{vecCellTypes(intRand+intCellsV1)}));
end
colormap(redblue)

%% decoding crap
matData = matModelRespP';
vecOriRads=ang2rad(vecTrialOris(vecTrialStimType))';
vecNeurons = 1:1200;
for intTrial=1:intTrials
	indTrials = true(intTrials,1);
	indTrials(intTrial) = false;
	vecThisOriRads = vecOriRads(indTrials);
	vecThisAct = matData(intTrial,vecNeurons);
	matThisData = matData(indTrials,vecNeurons);
	dblDecodedAngle = doGlobalOriEstimator(vecThisAct, matThisData, vecThisOriRads);

end
