%% initialize
%clearvars;
boolLoad = false;
boolSaveFigs = true;
if ~exist('matModelResp','var')
	clearvars -except bool*;%close all;
	intLoadSim = 11; %% SET WHICH SIM TO LOAD HERE; 41 vs 42; 2 vs. 5
	boolLoad = true;
end
if intLoadSim == 11 && boolLoad
	strSimulation = 'xAreaDistributed_OriFull_2017-06-15';
	
elseif intLoadSim == 12 && boolLoad
	strSimulation = 'OldRerunSquareGrating_736803.458055_2017-04-18';
	
	
elseif intLoadSim == 21 && boolLoad
	strSimulation = 'xAreaTransferUniform_736797.970520_2017-04-12';
	
elseif intLoadSim == 22 && boolLoad
	strSimulation = 'xAreaTransferFeedbackUniform_736797.882023_2017-04-12';
	
elseif intLoadSim == 23 && boolLoad
	strSimulation = 'PPxAreaTransferFeedbackInhUniform_736798.614076_2017-04-13';
elseif intLoadSim == 24 && boolLoad
	strSimulation = 'PPxAreaTransferFeedbackInh2Uniform_736798.670873_2017-04-13';
elseif intLoadSim == 25 && boolLoad
	strSimulation = 'PPxAreaTransferFeedback2Uniform_736798.746731_2017-04-13';
	
	
elseif intLoadSim == 31 && boolLoad
	strSimulation = 'xAreaDistributed_491429_2017-04-21';
elseif intLoadSim == 32 && boolLoad
	strSimulation = 'xAreaDistributed_498422_2017-04-28';
elseif intLoadSim == 33 && boolLoad
	strSimulation = 'xAreaDistributed_491429_2017-05-02';
	
	
elseif intLoadSim == 41 && boolLoad
	strSimulation = 'xAreaDistributed_491426_2017-04-03';
	
elseif intLoadSim == 42 && boolLoad
	strSimulation = 'xAreaDistributed_493425_2017-04-05';
	
elseif intLoadSim == 43 && boolLoad
	strSimulation = 'SquareGrating2017-04-05';
	
elseif intLoadSim == 5 && boolLoad
	strSimulation = 'SquareGrating2017-01-25';
	%Cells: 1440
	%Cort Conn: 403200
	%LGN Conn: 32256*2
	%Trials: 3000 / Types: 2 / Reps: 1500
	%Oris: [42.5 47.5]
elseif intLoadSim == 51 && boolLoad
	strSimulation = 'OldRerun150SquareGrating_736796.588630_2017-04-11';
	%original connectivity structure
elseif intLoadSim == 52 && boolLoad
	strSimulation = 'OldRerun150FixedSquareGrating_736796.589038_2017-04-11';
	
elseif intLoadSim == 53 && boolLoad
	strSimulation = 'PPOldRerunFixedSquareGrating_736796.741130_2017-04-11';
	
elseif intLoadSim == 6 && boolLoad
	strSimulation = 'SquareGrating2017-01-26';
	%Cells: 1440
	%Cort Conn: 403200
	%LGN Conn: 32256*2
	%Trials: 3600 / Types: 12 / Reps: 300
	%Oris: [0    15    30    45    60    75    90   105   120   135   150   165]
	%elseif intLoadSim == 7 && boolLoad
	%	strSimulation = 'SquareGrating2017-03-15';
	%Cells: 240
	%Cort Conn: 16800
	%LGN Conn: 5376*2
	%Trials: 20 / Types: 2 / Reps: 10
	%Oris: [0    15    30    45    60    75    90   105   120   135   150   165]
	%elseif intLoadSim == 8 && boolLoad
	%strSimulation = 'xAreaLine2017-03-21';
	%	strSimulation = 'xAreaSquareGrating2017-03-22';
	%elseif intLoadSim == 9 && boolLoad
	%strSimulation = 'xAreaLine2017-03-21';
	%	strSimulation = 'SquareGrating_736781.722426_2017-03-27';
end

%% RUN: #header
strFigDir = 'D:\Data\Results\V1_LIFmodel\figures\';
if boolLoad
	runModelHeader;
end

%% get data and define binning
matModelRespP = matModelResp;
intNeurons = size(matModelRespP,1);
intRawSizeT = numel(vecOverallT);
dblEndT = max(vecOverallT);
dblBinSize = 10/1000;%0.5/3; %50 ms bins
intBinsT = floor(dblEndT/dblBinSize);

vecBinsTime = linspace(0,dblEndT,intBinsT+1);
vecBinsStimStart = round(vecStimStartSecs/dblBinSize);
vecBinsStimStop = round(vecStimStopSecs/dblBinSize);
hTic = tic;
if ~exist('matSpikeCounts2D','var')
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
		save(['D:\Data\Processed\V1_LIFmodel\Simulation_' strSimulation '_prepro.mat'],'matModelResp','matSpikeCounts2D');
	end
end
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
	indSelectNeuronsV1 = vecCellTypes(:)==1 & vecCellArea(:)==1;
	cellStimTimePoints = cellfill(false(1,intBinsT),[1 intStimTypes]);
	for intTrial=1:intTrials
		cellStimTimePoints{vecTrialStimType(intTrial)}((vecBinsStimStart(intTrial)+1):vecBinsStimStop(intTrial)) = true;
	end
	matActPyrS1 = matSpikeCounts2D(indSelectNeuronsV1,cellStimTimePoints{intStim1})/dblBinSize;
	matActPyrS2 = matSpikeCounts2D(indSelectNeuronsV1,cellStimTimePoints{intStim2})/dblBinSize;
	vecPrefAngles=-vecPrefOri+pi;
	dblOriS1 = rad2ang(circ_mean(2*vecPrefAngles(indSelectNeuronsV1)',mean(matActPyrS1,2)))/2;
	dblOriS2 = rad2ang(circ_mean(2*vecPrefAngles(indSelectNeuronsV1)',mean(matActPyrS2,2)))/2;
	
	%%
	figure
	subplot(2,3,1)
	scatter(vecPrefAngles(indSelectNeuronsV1)-0.005,xmean(matActPyrS1,2)','b.');
	hold on
	scatter(vecPrefAngles(indSelectNeuronsV1)+0.005,xmean(matActPyrS2,2)','r.');
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
	matActPyrV1S1 = matSpikeCounts2D(indSelectNeuronsV1,cellStimTimePoints{intStim1})/dblBinSize;
	matActPyrV1S2 = matSpikeCounts2D(indSelectNeuronsV1,cellStimTimePoints{intStim2})/dblBinSize;
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
	scatter(1:sum(indSelectNeuronsV1),vecActV1S1R,'b.');
	scatter(1:sum(indSelectNeuronsV1),vecActV1S2R,'r.');
	
	hold off
	xlim([0 sum(indSelectNeuronsV1)]);
	xlabel('V1 Pyramid sorted by stimulus preference');
	ylabel('Mean single neuron spiking rate (Hz)');
	title(sprintf('Activity of V1 Pyramids'));
	fixfig
	
	
	%
	matActHz = matSpikeCounts2D/dblBinSize;
	vecFilt = normpdf(-2:2,0,1);
	vecFilt=vecFilt/sum(vecFilt);
	vecActHz = conv(xmean(matActHz(indSelectNeuronsV1,:),1),vecFilt,'same');
	
	subplot(2,3,3)
	cellColor = {[0 0 1],[1 0 0]};
	vecLimX = [0 dblEndT];
	vecLimY = [0 max(vecActHz(:))];
	hold on;
	intC = 0;
	for intType=[intStim1 intStim2]
		intC=intC+1;
		for intTrial=find(vecTrialStimType==intType)
			dblB = vecBinsStimStart(intTrial);
			dblE = vecBinsStimStop(intTrial);
			patch([dblB dblE dblE dblB]*dblBinSize,[vecLimY(1) vecLimY(1) vecLimY(end) vecLimY(end)],cellColor{intC},'EdgeColor','none','FaceAlpha',.5);
		end
	end
	plot(vecBinsTime,vecActHz,'k');
	hold off
	ylabel('Mean population spiking rate (Hz)');
	xlabel('Time (s)');
	title(sprintf('Population activity of V1 Pyramids over time; S1, blue; S2, red'));
	fixfig
	
	%% V2
	if sum(vecCellArea==2) > 0
		indSelectNeuronsV2 = vecCellTypes(:)==1 & vecCellArea(:)==2;
		matActPyrV2S1 = matSpikeCounts2D(indSelectNeuronsV2,cellStimTimePoints{intStim1})/dblBinSize;
		matActPyrV2S2 = matSpikeCounts2D(indSelectNeuronsV2,cellStimTimePoints{intStim2})/dblBinSize;
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
		
		
		vecActHzV2 = conv(xmean(matActHz(indSelectNeuronsV2,:),1),vecFilt,'same');
		subplot(2,3,6)
		cellColor = {[0 0 1],[1 0 0]};
		vecLimX = [0 dblEndT];
		vecLimY = [0 max(vecActHzV2(:))];
		hold on;
		intC = 0;
		for intType=[intStim1 intStim2]
			intC=intC+1;
			for intTrial=find(vecTrialStimType==intType)
				dblB = vecBinsStimStart(intTrial);
				dblE = vecBinsStimStop(intTrial);
				patch([dblB dblE dblE dblB]*dblBinSize,[vecLimY(1) vecLimY(1) vecLimY(end) vecLimY(end)],cellColor{intC},'EdgeColor','none','FaceAlpha',.5);
			end
		end
		plot(vecBinsTime,vecActHzV2,'k');
		hold off
		ylabel('Mean population spiking rate (Hz)');
		xlabel('Time (s)');
		title(sprintf('Population activity of V2 Pyramids over time; S1, blue; S2, red'));
		fixfig
		
		
	end
	%full screen
	drawnow;
	jFig = get(handle(gcf), 'JavaFrame');
	jFig.setMaximized(true);
	figure(gcf);
	drawnow;
	
	%save figure
	if boolSaveFigs
		if (sum(vecTrialStimType==intStim1) + sum(vecTrialStimType==intStim2)) > 1000 %remove panels
			subplot(2,3,3);cla;
			subplot(2,3,6);cla;
		end
		export_fig([strFigDir 'PopulationResponses' strSimulation '.tif']);
		export_fig([strFigDir 'PopulationResponses' strSimulation '.pdf']);
	end
end
%% analyze connectivity
fprintf('Building connectivity matrix... [%s]\n',getTime);
vecSynWeight = ones(size(vecSynWeight));
vecSign = vecSynExcInh;%vecCortSynType;
vecSign(vecSign==2) = -1;
matConn2D = zeros(intNeurons,intNeurons);
%matConn2D = getFillGrid(matConn2D,matCortConn(:,1),matCortConn(:,2),vecCortConductance.*vecSign);
matConn2D = getFillGrid(matConn2D,matSynFromTo(:,1),matSynFromTo(:,2),vecSynConductance.*vecSign.*vecSynWeight);
intRank = rank(matConn2D);

% pca
[coeff,score,latent,tsquared,explained,mu] = pca(matConn2D);
vecR2_PCA = explained/100;



%% plot
figure

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

subplot(2,3,[1 2 4 5])
%surf(matConn2D,'EdgeColor','none');
%imagesc(matConn2D,[-0.6 0.6+eps]);
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

%save figure
if boolSaveFigs
	export_fig([strFigDir 'ConnectivityMatrix' strSimulation '.tif']);
	export_fig([strFigDir 'ConnectivityMatrix' strSimulation '.pdf']);
end


%%
figure
intNeuronsV2=size(matFieldsV2,3);
for intNeuron=1:15%intNeuronsV2
	subplot(3,5,intNeuron);
	matF = matFieldsV2(:,:,intNeuron);
	imagesc(matF);%,max(abs(matF(:)))*[-1 1]);
end
%colormap(redblue)
