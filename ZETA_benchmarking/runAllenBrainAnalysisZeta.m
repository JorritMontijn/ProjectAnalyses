
%% define locations & load data
clear all;clc;

strDisk = 'F:';
strDataSource = '\Data\Processed\AllenBrainVisualEphys\';
strFigTarget = '\Data\Results\ZETA\AllenBrainVisualEphys\';

%load metadata
sCSV = loadcsv(strcat(strDisk,strDataSource,'sessions.csv'));
%select VIP
vecSessions = sCSV.id(contains(sCSV.genotype,'Vip'));

strDataSourceAgg = strcat(strDataSource,'Aggregates\');
strDataFile = 'AggSes2020-06-24.mat'; %VIP
%strDataFile = 'AggSes2020-07-01.mat'; %VIP
load(strcat(strDataSourceAgg,strDataFile));

%% classify neurons by onset latency
cellCellType = {'none','artifact','VIP','late','sustained'};
intSignZETA = 0;
intConsider = 0;
dblEpochEnd = 100/1000;
vecPeakWindow = [10/1000 30/1000];
for intSes=4:numel(sSes)
	%get peak times
	vecPeaks = sSes(intSes).matLatencies(3,:);
	vecOnsets = sSes(intSes).matLatencies(4,:);

	%classify
	dblStartT = 0.5;
	vecCellType = ones(size(vecPeaks));
	indConsider = vecPeaks > dblStartT;
	indArtifact = indConsider & (vecPeaks < (dblStartT + 0.001));
	indVIP = indConsider & (vecPeaks < (dblStartT + 0.010)) & ~indArtifact;
	indLate = indConsider & (vecPeaks < (dblStartT + vecPeakWindow(2))) & ~indArtifact & ~indVIP;
	indPossiblyVisual = indConsider & (vecPeaks > (dblStartT + vecPeakWindow(2)));
	vecCellType = vecCellType + indArtifact*1 + indVIP*2 + indLate*3 + indPossiblyVisual*4;
	fprintf('Ses %d: none=%d, artifact=%d, VIP=%d, late=%d, sustained=%d\n',intSes,sum(vecCellType==1),sum(vecCellType==2),sum(vecCellType==3),sum(vecCellType==4),sum(vecCellType==5));
	intSignZETA = intSignZETA + sum(~isnan(vecPeaks));
	intConsider = intConsider + sum(indConsider);
	
	%rate
	hTic = tic;
	vecPreRate = nan(size(vecCellType));
	vecDurRate = nan(size(vecCellType));
	vecPostRate = nan(size(vecCellType));
	vecLateRate = nan(size(vecCellType));
	vecSustRate = nan(size(vecCellType));
	
	%save data
	cellPreRate = cell(size(vecCellType));
	cellDurRate = cell(size(vecCellType));
	cellPostRate = cell(size(vecCellType));
	cellLateRate = cell(size(vecCellType));
	cellSustRate = cell(size(vecCellType));
	cellRateT = cell(size(vecCellType));
	cellRate = cell(size(vecCellType));
	cellBinRateT = cell(size(vecCellType));
	cellBinRate = cell(size(vecCellType));
	
	vecTestCells = find(indVIP | indLate | indPossiblyVisual);
	vecRealVIPs = [];
	intResampNum = 100;
	vecPreOnsets = sSes(intSes).vecOptoEventsT - dblStartT;
	intPlot = 3;
	intLatencyPeaks = 4;
	for intUseNeuron = 1:numel(vecTestCells)
		intNeuron = vecTestCells(intUseNeuron);
		if toc(hTic) > 5,hTic=tic;fprintf('Proc %d.3: Neuron %d/%d [%s]\n',intSes,intUseNeuron,numel(vecTestCells),getTime);end
		
		%check if rate is enhanced during 10-ms stim window
		vecSpikeTimes = sSes(intSes).cellSpikes{intNeuron};
		vecBinE = sort(cat(1,sSes(intSes).vecOptoEventsT-50/1000,...
			sSes(intSes).vecOptoEventsT,...
			sSes(intSes).vecOptoEventsT+10/1000,...
			sSes(intSes).vecOptoEventsT+vecPeakWindow(2),...
			sSes(intSes).vecOptoEventsT+60/1000,...
			sSes(intSes).vecOptoEventsT+100/1000));
		vecBinDur = diff(vecBinE);
		
		vecC = histcounts(vecSpikeTimes,vecBinE);
		vecR = vecC(:)./vecBinDur;
		vecThisPreRate = vecR(1:6:end);
		vecThisDurRate = vecR(2:6:end);
		vecThisPostRate = vecR(3:6:end);
		vecThisLateRate = vecR(4:6:end);
		vecThisSustRate = vecR(5:6:end);
		%save data
		vecPreRate(intNeuron) = mean(vecThisPreRate);
		vecDurRate(intNeuron) = mean(vecThisDurRate);
		vecPostRate(intNeuron) = mean(vecThisPostRate);
		vecLateRate(intNeuron) = mean(vecThisLateRate);
		vecSustRate(intNeuron) = mean(vecThisSustRate);
		%save data
		cellPreRate{intNeuron} = vecThisPreRate;
		cellDurRate{intNeuron} = vecThisDurRate;
		cellPostRate{intNeuron} = vecThisPostRate;
		cellLateRate{intNeuron} = vecThisLateRate;
		cellSustRate{intNeuron} = vecThisSustRate;
		
		%change type from VIP if reduction during stim
		if strcmpi(cellCellType{vecCellType(intNeuron)},'VIP')
			if vecPreRate(intNeuron) > vecDurRate(intNeuron)
				vecCellType(intNeuron) = vecCellType(intNeuron) + 1;
			end
		end
		
		%plot
		dblStimT = dblEpochEnd;
		dblUseMaxDur = dblStimT*2;
		matEventTimes = sSes(intSes).vecOptoEventsT-dblStimT;
		intPlot=0;
		%[dblZetaP,vecLatencies,sZETA,sRate] = getZeta(vecSpikeTimes,matEventTimes,dblUseMaxDur,2,intPlot,intLatencyPeaks,vecPeakWindow);
		[vecIFR,sRate] = getIFR(vecSpikeTimes,matEventTimes,dblUseMaxDur);
		title(subplot(2,3,1),sprintf('Ses%dN%d,Pre=%.1fHz,Dur=%.1fHz,Post=%.1fHz,Late=%.1fHz,Sust=%.1fHz',...
			intSes,intNeuron,mean(vecThisPreRate),mean(vecThisDurRate),mean(vecThisPostRate),mean(vecThisLateRate),mean(vecThisSustRate)));
		
		%save IFR
		if isempty(sRate)
			sRate.vecT = [];
			sRate.vecRate = [];
		end
		
		cellRateT{intNeuron} = sRate.vecT;
		cellRate{intNeuron} = sRate.vecRate;
		
		%calculate normal bins
		dblStep = 2.5/1000;
		vecBinEdges = 0:dblStep:(dblUseMaxDur);
		[vecMean,vecSEM,vecWindowBinCenters] = doPEP(vecSpikeTimes,vecBinEdges-dblUseMaxDur/2,sSes(intSes).vecOptoEventsT,-1);
		
		cellBinRateT{intNeuron} = vecWindowBinCenters;
		cellBinRate{intNeuron} = vecMean;
	
		if 0%indLate(intNeuron)
			%%
			[dblZetaP,vecLatencies,sZETA,sRate] = getZeta(vecSpikeTimes,matEventTimes,dblUseMaxDur,2,3,intLatencyPeaks);
			title(subplot(2,3,1),sprintf('Ses%dN%d,Pre=%.1fHz,Dur=%.1fHz,Post=%.1fHz,Late=%.1fHz,Sust=%.1fHz',...
				intSes,intNeuron,mean(vecThisPreRate),mean(vecThisDurRate),mean(vecThisPostRate),mean(vecThisLateRate),mean(vecThisSustRate)));
			pause
			close;
			%% save
			strDir = 'F:\Data\Results\ZETA\AllenBrainVisualEphys\';
			strFile = sprintf('ExampleNeuron_%dN%d',vecSessions(intSes),intNeuron);
			%export_fig([strDir strFile '.tif']);
			%export_fig([strDir strFile '.pdf']);
		end
	end
	
	%% edit real vips & save data
	sSes(intSes).vecCellType = vecCellType;
	sSes(intSes).vecPreRate = vecPreRate;
	sSes(intSes).vecDurRate = vecDurRate;
	sSes(intSes).vecPostRate = vecPostRate;
	sSes(intSes).vecLateRate = vecLateRate;
	sSes(intSes).vecSustRate = vecSustRate;
	sSes(intSes).cellPreRate = cellPreRate;
	sSes(intSes).cellDurRate = cellDurRate;
	sSes(intSes).cellPostRate = cellPostRate;
	sSes(intSes).cellLateRate = cellLateRate;
	sSes(intSes).cellSustRate = cellSustRate;
	sSes(intSes).cellRateT = cellRateT;
	sSes(intSes).cellRate = cellRate;
	sSes(intSes).cellBinRateT = cellBinRateT;
	sSes(intSes).cellBinRate = cellBinRate;
	
end
intSignZETA
intConsider
	

%% aggregate all cells
vecAggMeanRate = [];
vecAggPeaks = [];
vecAggOnsets = [];
vecAggCellType = [];
vecAggPreRate = [];
vecAggDurRate = [];
vecAggPostRate = [];
vecAggLateRate = [];
vecAggSustRate = [];
cellAggPreRate = [];
cellAggDurRate = [];
cellAggPostRate = [];
cellAggLateRate = [];
cellAggSustRate = [];
cellAggRateT = [];
cellAggRate = [];
cellAggBinRateT = [];
cellAggBinRate = [];

%build bins
dblStep = (2.5/1000);
vecBinEdges = (-dblStep/2):dblStep:(dblUseMaxDur+dblStep);
%vecPlotT = (vecBinEdges(2:end)-median(diff(vecBinEdges))/2-dblUseMaxDur/2)*1000;
vecPlotT = vecWindowBinCenters*1000;
matAct = [];
for intSes=4:numel(sSes)
	indKeep = ~isnan(sSes(intSes).vecPreRate);
	dblSesDur = (max(cell2vec(sSes(intSes).cellSpikes(indKeep)))-min(cell2vec(sSes(intSes).cellSpikes(indKeep))));
	vecAggMeanRate = cat(2,vecAggMeanRate,cellfun(@numel,sSes(intSes).cellSpikes(indKeep))./dblSesDur);
	vecAggPeaks = cat(2,vecAggPeaks,sSes(intSes).matLatencies(3,indKeep));
	vecAggCellType = cat(2,vecAggCellType,sSes(intSes).vecCellType(indKeep));
	vecAggPreRate = cat(2,vecAggPreRate,sSes(intSes).vecPreRate(indKeep));
	vecAggDurRate = cat(2,vecAggDurRate,sSes(intSes).vecDurRate(indKeep));
	vecAggPostRate = cat(2,vecAggPostRate,sSes(intSes).vecPostRate(indKeep));
	vecAggLateRate = cat(2,vecAggLateRate,sSes(intSes).vecLateRate(indKeep));
	vecAggSustRate = cat(2,vecAggSustRate,sSes(intSes).vecSustRate(indKeep));
	cellAggPreRate = cat(2,cellAggPreRate,sSes(intSes).cellPreRate(indKeep));
	cellAggDurRate = cat(2,cellAggDurRate,sSes(intSes).cellDurRate(indKeep));
	cellAggPostRate = cat(2,cellAggPostRate,sSes(intSes).cellPostRate(indKeep));
	cellAggLateRate = cat(2,cellAggLateRate,sSes(intSes).cellLateRate(indKeep));
	cellAggSustRate = cat(2,cellAggSustRate,sSes(intSes).cellSustRate(indKeep));
	cellAggRateT = cat(2,cellAggRateT,sSes(intSes).cellRateT(indKeep));
	cellAggRate = cat(2,cellAggRate,sSes(intSes).cellRate(indKeep));
	cellAggBinRateT = cat(2,cellAggBinRateT,sSes(intSes).cellBinRateT(indKeep));
	cellAggBinRate = cat(2,cellAggBinRate,sSes(intSes).cellBinRate(indKeep));

end
%%
%build mean activity vector
matActBinned = [];
intKeepN = numel(cellAggRateT);
vecPeakPosT = nan(1,intKeepN);
vecPeakNegT = nan(1,intKeepN);
for intNeuron=1:intKeepN
	%vecAct = makeBins(cellAggRateT{intNeuron},cellAggRate{intNeuron},vecBinEdges);
	vecActBinned = cellAggBinRate{intNeuron};
	matActBinned(intNeuron,:) = vecActBinned;
	
	vecT = cellAggRateT{intNeuron}-dblStimT;
	if isempty(vecT),continue;end
	vecAct = cellAggRate{intNeuron};
	indPossT = vecT>vecPeakWindow(1) & vecT<vecPeakWindow(2);
	vecPossT = vecT(indPossT);
	[a,intPosT]=max(vecAct(indPossT));
	dblPeakPosT = vecPossT(intPosT);
	if isempty(dblPeakPosT)
		vecPeakPosT(intNeuron) = nan;
	else
		vecPeakPosT(intNeuron) = dblPeakPosT;
	end
	
	[a,intNegT]=min(vecAct(indPossT));
	dblPeakNegT = vecPossT(intNegT);
	if isempty(dblPeakNegT)
		vecPeakNegT(intNeuron) = nan;
	else
		vecPeakNegT(intNeuron) = dblPeakNegT;
	end
end
matActBinnedZ = zscore(matActBinned,[],2);
vecPreMean = mean(matActBinnedZ(:,1:floor(size(matActBinnedZ,2)/2)),2);
matActBinnedZ = bsxfun(@minus,matActBinnedZ,vecPreMean);

matC = lines(4);
%VIP
indVip = vecAggCellType==3;
matR = matActBinnedZ(indVip,:);
vecMeanVip = mean(matR,1);
vecSemVip = std(matR,[],1)./sqrt(sum(indVip));
figure
subplot(2,2,2)
hold on
%errorbar(vecPlotT,vecMeanVip,vecSemVip);
plot(vecPlotT,vecMeanVip,'color',matC(1,:));
plot(vecPlotT,vecMeanVip-vecSemVip,'color',matC(1,:));
plot(vecPlotT,vecMeanVip+vecSemVip,'color',matC(1,:));

%Inh
[a,vecSigPostRateDiff]=cellfun(@ttest,cellAggPreRate,cellAggPostRate);
indInh = vecAggCellType>3 ...
	& ((vecAggPostRate < vecAggPreRate) & (vecSigPostRateDiff<0.05));

matR = matActBinnedZ(indInh,:);
vecMeanInh = mean(matR,1);
vecSemInh = std(matR,[],1)./sqrt(sum(indInh));
%errorbar(vecPlotT,vecMeanInh,vecSemInh);
plot(vecPlotT,vecMeanInh,'color',matC(2,:));
plot(vecPlotT,vecMeanInh-vecSemInh,'color',matC(2,:));
plot(vecPlotT,vecMeanInh+vecSemInh,'color',matC(2,:));


%act
[a,vecSigPostRateDiff]=cellfun(@ttest,cellAggPreRate,cellAggPostRate);
indAct = ~indInh & vecAggCellType>3 ...
	& ((vecAggPostRate > vecAggPreRate) & (vecSigPostRateDiff<0.05));

matR = matActBinnedZ(indAct,:);
vecMeanAct = mean(matR,1);
vecSemAct = std(matR,[],1)./sqrt(sum(indAct));


%errorbar(vecPlotT,vecMeanAct,vecSemAct);
plot(vecPlotT,vecMeanAct,'color',matC(3,:));
plot(vecPlotT,vecMeanAct-vecSemAct,'color',matC(3,:));
plot(vecPlotT,vecMeanAct+vecSemAct,'color',matC(3,:));


%other
indOther = vecAggCellType~=2 & ~indInh & ~indAct & ~indVip;
matR = matActBinnedZ(indOther,:);
vecMeanOther = mean(matR,1);
vecSemOther = std(matR,[],1)./sqrt(sum(indOther));
%errorbar(vecPlotT,vecMeanOther,vecSemOther,'k');
plot(vecPlotT,vecMeanOther,'color',matC(4,:));
plot(vecPlotT,vecMeanOther-vecSemOther,'color',matC(4,:));
plot(vecPlotT,vecMeanOther+vecSemOther,'color',matC(4,:));


xlabel('Time after opto start (ms)')
ylabel('Normalized firing rate');
legend({'VIP','Inh','dAct','Other'});
fixfig;

indKeep = indVip | indInh | indAct;
vecType = indVip*1 + indInh*2 + indAct*3;
matActSubBinnedZ = matActBinnedZ(indKeep,:);
[dummy,vecReorder] = sort(vecType(indKeep));
hS1=subplot(2,2,1)
imagesc(vecPlotT,1:size(matActSubBinnedZ,1),matActSubBinnedZ(vecReorder,:))
ylabel('Neuron');
xlabel('Time after opto start (ms)');
title(sprintf('Normalized activity of significant neurons; VIP n=%d',sum(indVip)));
h=colorbar;
ylabel(h,'Norm. Act');
fixfig;
grid off;

[p,h,stats] = ranksum(vecPeakNegT(indInh),vecPeakPosT(indAct));
subplot(2,2,3)
hold on;
bplot(1000*vecPeakNegT(indInh),1);
bplot(1000*vecPeakPosT(indAct),2);
set(gca,'xtick',[1 2],'xticklabel',{'Inh','Act'});
ylabel('Peak latency (ms)');
hold off
title(sprintf('Mann-Whitney,p=%.3e;inh n=%d (%.3fs),act n=%d (%.3fs)',p,sum(indInh),nanmedian(vecPeakNegT(indInh)),sum(indAct),nanmedian(vecPeakPosT(indAct))));
ylim(1000*[0 vecPeakWindow(2)])
fixfig;
maxfig;

% cluster
matC=lines(4);
hS3=subplot(2,2,4);
dblStep = 2.5/1000;
vecDistroE = (vecPeakWindow(1)-dblStep):dblStep:vecPeakWindow(2);
vecDistroC = vecDistroE(2:end) - median(diff(vecDistroE))/2;
vecDistroNeg = histcounts(vecPeakNegT(indInh),vecDistroE);
vecDistroNegStairs = [vecDistroNeg vecDistroNeg(end)];
vecDistroPos = histcounts(vecPeakPosT(indAct),vecDistroE);
vecDistroPosStairs = [vecDistroPos vecDistroPos(end)];
stairs(vecDistroE(1:end)*1000,vecDistroNegStairs/sum(vecDistroNegStairs),'color',matC(2,:))
hold on;
stairs(vecDistroE(1:end)*1000,vecDistroPosStairs/sum(vecDistroPosStairs),'color',matC(3,:))
hold off;
xlabel('Peak/trough time (ms)')
ylabel('Fraction of cells');
xlim([0 vecPeakWindow(2)]*1000);
fixfig;

%{
%% cluster
vecD = pdist(matActSubBinnedZ);
matD = squareform(vecD);
[intClustNum,vecSilhouetteD]=doFastClustering(matD,50);

Z = linkage(matActSubBinnedZ,'ward');
vecClusterID = cluster(Z,'maxclust',intClustNum);

hS5=subplot(2,3,5);
[dummy,vecReorder2] = sort(vecClusterID);

imagesc(vecPlotT,1:size(matActSubBinnedZ,1),matActSubBinnedZ(vecReorder2,:))
ylabel('Neuron');
xlabel('Time after opto start (ms)');
h=colorbar;
ylabel(h,'Norm. Act');
fixfig;

% cluster
hS6=subplot(2,3,6);
matCorr = corr(matActSubBinnedZ(vecReorder2,:)');
imagesc(matCorr,[-1 1]);
colorbar;
colormap(hS6,redblue);
%}
%%
return
%% save figure
export_fig([strFigTarget 'VIP_disinhibition.tif']);
export_fig([strFigTarget 'VIP_disinhibition.pdf']);