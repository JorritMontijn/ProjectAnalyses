clear all;
%close all;
strPath = 'F:\Data\Processed\ZETA\Inclusion\';
strFigPath = 'F:\Data\Results\ZETA\Inclusion\';
cellUniqueAreas = {...
	'lateral geniculate',...Area 7
	'Primary visual',...Area 8
	'Lateral posterior nucleus',...Area 9
	'Anterior pretectal nucleus',...Area 10
	'Nucleus of the optic tract',...Area 11
	'Superior colliculus',...Area 12
	'Anteromedial visual',...Area 13
	'posteromedial visual',...%,...Area 14
	'Anterolateral visual',...Area 15
	'Lateral visual',...Area 16
	'Rostrolateral area',...Area 17
	'Anterior area'};%,...Area 18
cellRunRand = {...
	'',...Rand 1
	'-Rand',...Rand 2
	};
cellRepStr = {...
	'RunDriftingGratings','';...
	'RunNaturalMovie','-NM';...
	'lateral geniculate','LGN';...
	'Primary visual','V1';...
	'Lateral posterior nucleus','LP';...
	'Anterior pretectal nucleus','APN';...
	'Nucleus of the optic tract','NOT';...
	'Superior colliculus','SC';...
	'Anteromedial visual','AM';...
	'posteromedial visual','PM';...
	'Anterolateral visual','AL';...
	'Lateral visual','L';...
	'Rostrolateral area','RL';...
	'Anterior area','A';...
	'Subiculum','H-Sub';...
	'Field CA1','H-CA1';...
	'Field CA2','H-CA2';...
	'Field CA3','H-CA3';...
	'Dentate gyrus','H-DG';...
	'Retrosplenial','RetSpl';...
	};

%% prep
boolSavePlots = true;
cellDatasetNames = {};
cellNumSpikes = [];

sSignif = [];

%comput time
sComput = [];

%p-value vectors
sP = [];

strArea = ''; %V1, SC, Retina, Poisson, GCaMP

%% get all files
sDir=dir([strPath 'ZetaTtest*Resamp*.mat']);
cellNames = {sDir(:).name};
indRand = contains(cellNames,'Rand');
cellResamps = cellfun(@getFlankedBy,cellNames,cellfill('Resamp',size(cellNames)),cellfill('.mat',size(cellNames)),'UniformOutput',false);
vecResamps = str2double(cellResamps);
vecUniqueResamps = sort(unique(vecResamps));
intResamps = numel(vecUniqueResamps);
for intRandType=1:2
	%set var
	vecR = [0 1];
	intRand = vecR(intRandType); %[1 0]
	strRand = cellRunRand{intRandType};
	
	for intResampIdx=1:intResamps
		intThisResamp = vecUniqueResamps(intResampIdx);
		intFile = find(intThisResamp==vecResamps & indRand == intRand);
		sLoad=load([strPath cellNames{intFile}]);
		
		%sig at alpha=0.05
		dblAlpha = 0.05;
		sSignif.matBISIP_g(intResampIdx,intRandType) = sum(sLoad.vecBISIP_g<dblAlpha)/numel(sLoad.vecBISIP_g);
		sSignif.matBISIP_ks(intResampIdx,intRandType) = sum(sLoad.vecBISIP_ks<dblAlpha)/numel(sLoad.vecBISIP_ks);
		sSignif.matHzP(intResampIdx,intRandType) = sum(sLoad.vecHzP<dblAlpha)/numel(sLoad.vecHzP);
		sSignif.matISIP_IntG(intResampIdx,intRandType) = sum(sLoad.vecISIP_IntG<dblAlpha)/numel(sLoad.vecISIP_IntG);
		sSignif.matISIP_g(intResampIdx,intRandType) = sum(sLoad.vecISIP_g<dblAlpha)/numel(sLoad.vecISIP_g);
		sSignif.matISIP_ks(intResampIdx,intRandType) = sum(sLoad.vecISIP_ks<dblAlpha)/numel(sLoad.vecISIP_ks);
		sSignif.matPoissP(intResampIdx,intRandType) = sum(sLoad.vecPoissP<dblAlpha)/numel(sLoad.vecPoissP);
		sSignif.matZetaP(intResampIdx,intRandType) = sum(sLoad.vecZetaP<dblAlpha)/numel(sLoad.vecZetaP);
		
		%comput time
		sComput.cellComputTimeBISI{intResampIdx,intRandType} = sLoad.vecComputTimeBISI;
		sComput.cellComputTimePoiss{intResampIdx,intRandType} = sLoad.vecComputTimePoiss;
		sComput.cellComputTimeISI{intResampIdx,intRandType} = sLoad.vecComputTimeISI;
		sComput.cellComputTimeZETA{intResampIdx,intRandType} = sLoad.vecComputTimeZETA;
		
		%p-value vectors
		sP.Zeta{intResampIdx,intRandType} = sLoad.vecZetaP; %light blue
		sP.Z_ISI_g{intResampIdx,intRandType} = sLoad.vecBISIP_g; %dark blue
		sP.Ttest{intResampIdx,intRandType} = sLoad.vecHzP; %black
		sP.ISI_IntG{intResampIdx,intRandType} = sLoad.vecISIP_IntG; %purple
		sP.ISI_g{intResampIdx,intRandType} = sLoad.vecISIP_g; %red
		sP.Z_ISI_ks{intResampIdx,intRandType} = sLoad.vecBISIP_ks; %red-purple burgundy
		sP.ISI_ks{intResampIdx,intRandType} = sLoad.vecISIP_ks; %orange
		sP.Poiss{intResampIdx,intRandType} = sLoad.vecPoissP; %grey
		
		%other
		cellNumSpikes{intResampIdx,intRandType} = sLoad.vecNumSpikes;
		matResamps(intResampIdx,intRandType) = intThisResamp;
	end
	
end

%% calc AUCs
%prep
cellTests = fieldnames(sP);
intTests = numel(cellTests);
matAUC = nan(intTests,intResamps);
matAUC_ci = nan(intTests,intResamps,2);
cellColor = {[0 0.4 0.8],...
	[0 0 0.5],...
	[0 0 0],...
	[0.5 0 1],...
	[1 0 0],...
	[0.8 0 0.3],...
	[1 0.5 0],...
	[0.5 0.5 0.5]};
vecAUC = nan(1,intTests);
%run
figure
hold on
for intTest=1:intTests
	%intTest = intTest + 1;
	cellData = sP.(cellTests{intTest});
	for intResampIdx=1:intResamps
		intThisResamp = vecUniqueResamps(intResampIdx);
		
		intCells = numel(cellData{intResampIdx,1});
		vecBothData = cat(2,cellData{intResampIdx,1},cellData{intResampIdx,2});
		vecBothLabels = cat(2,zeros(size(cellData{intResampIdx,1})),ones(size(cellData{intResampIdx,2})));
		vecThresholds = sort(vecBothData);
		vecRealP = cellData{intResampIdx,1};
		vecShuffP = cellData{intResampIdx,2};
		
		vecTP = sum(vecRealP<vecThresholds',2)/intCells;
		vecFP = sum(vecShuffP<vecThresholds',2)/intCells;
		[dblAUC,Aci] = auc(cat(1,vecBothLabels,vecBothData)');
		matAUC(intTest,intResampIdx) = dblAUC;
		matAUC_ci(intTest,intResampIdx,:) = Aci;
	end
	
	%plot
	plot(vecUniqueResamps,matAUC(intTest,:),'color',cellColor{intTest});
	plot(vecUniqueResamps,matAUC_ci(intTest,:,1),'--','color',cellColor{intTest});
	plot(vecUniqueResamps,matAUC_ci(intTest,:,2),'--','color',cellColor{intTest});
end
set(gca,'xscale','log')
return
%% plot ROC
if ~isempty(sP) && size(sP.Zeta,1) >= intResampIdx && ~isempty(sP.Zeta{intResampIdx,1})
	%prep data
	cellTests = fieldnames(sP);
	intTests = numel(cellTests);
	cellColor = {[0 0.4 0.8],...
		[0 0 0.5],...
		[0 0 0],...
		[0.5 0 1],...
		[1 0 0],...
		[0.8 0 0.3],...
		[1 0.5 0],...
		[0.5 0.5 0.5]};
	vecAUC = nan(1,intTests);
	cellLegend = cellTests;
	%plot
	
	figure;
	%vecH(intResampNpx) = subplot(4,3,intResampNpx);
	subplot(2,2,1)
	maxfig;
	hold on;
	
	for intTest=1:intTests
		%intTest = intTest + 1;
		cellData = sP.(cellTests{intTest});
		
		intCells = numel(cellData{intResampIdx,1});
		vecBothData = cat(2,cellData{intResampIdx,1},cellData{intResampIdx,2});
		vecBothLabels = cat(2,zeros(size(cellData{intResampIdx,1})),ones(size(cellData{intResampIdx,2})));
		vecThresholds = sort(vecBothData);
		vecRealP = cellData{intResampIdx,1};
		vecShuffP = cellData{intResampIdx,2};
		
		strTest = cellLegend{intTest};
		vecTP = sum(vecRealP<vecThresholds',2)/intCells;
		vecFP = sum(vecShuffP<vecThresholds',2)/intCells;
		
		plot(vecFP,vecTP,'Color',cellColor{intTest});
		
		[dblAUC,Aci] = auc(cat(1,vecBothLabels,vecBothData)');
		vecAUC(intTest) = dblAUC;
		cellLegend{intTest} = [strTest sprintf(', AUC=%.3f',dblAUC)];
	end
	hold off;
	xlabel('False positive fraction');
	ylabel('Inclusion fraction');
	fixfig;
	legend(cellLegend,'location','best','interpreter','none')
	
	%{
		%shuffled p MIMI
		subplot(2,2,2)
		dblStep = 0.05;
		vecBinsP = 0:dblStep:1;
		vecBinsPlot = (dblStep/2):dblStep:(1-dblStep/2);
		vecCountsHzP = histcounts(cellHzP{intResamp,2},vecBinsP);
		vecCountsISIP = histcounts(cellISIP{intResamp,2},vecBinsP);
		vecCountsZETAP = histcounts(cellZetaP{intResamp,2},vecBinsP);
		hold on
		plot(vecBinsPlot,vecCountsHzP,'k');
		plot(vecBinsPlot,vecCountsISIP,'r');
		plot(vecBinsPlot,vecCountsZETAP,'b');
		hold off
		title(sprintf('False positive p-value distribution'));
		xlabel('p-value of shuffled control');
		ylabel('Number of cells (count)');
		fixfig;
		legend({'Mean-rate t-test',[strTest ' shuffle'],'ZETA'},'location','best')
	%}
	%comput time
	subplot(2,2,2)
	vecComputTimeBISI = sComput.cellComputTimeBISI{intResampIdx,1};
	vecComputTimePoiss = sComput.cellComputTimePoiss{intResampIdx,1};
	vecComputTimeISI = sComput.cellComputTimeISI{intResampIdx,1};
	vecComputTimeZETA = sComput.cellComputTimeZETA{intResampIdx,1};
	
	vecBinsC = logspace(log10(min(vecComputTimeZETA))-1,log10(max(vecComputTimeZETA))+1,31);
	vecBinsPlotC = vecBinsC;
	vecCountsISIC = histcounts(vecComputTimeISI,vecBinsC);
	vecCountsZETAC = histcounts(vecComputTimeZETA,vecBinsC);
	vecCountsBISIC = histcounts(vecComputTimeBISI,vecBinsC);
	vecCountsPoiss = histcounts(vecComputTimePoiss,vecBinsC);
	hold on
	plot(vecBinsPlotC(2:end),vecCountsZETAC,'color',cellColor{1});
	plot(vecBinsPlotC(2:end),vecCountsBISIC,'color',cellColor{2});
	plot(vecBinsPlotC(2:end),vecCountsISIC,'color',cellColor{4});
	plot(vecBinsPlotC(2:end),vecCountsPoiss,'color',cellColor{end});
	hold off
	title(sprintf('Computation time comparison'));
	set(gca,'xscale','log');
	xlabel('Computation time per cell (s)');
	ylabel('Number of cells (count)');
	fixfig;
	legend({'ZETA','ZETA-ISI','ISI','Poisson'},'location','best')
	%}
	if boolSavePlots
		% save figure
		drawnow;
		strFigName = sprintf('%s_ROC_%s','ZetaBurst',strArea);
		export_fig([strFigPath strFigName '.tif']);
		export_fig([strFigPath strFigName '.pdf']);
	end
	
	%get quantiles
	vecQ1 = quantile(cellNumSpikes{1}/720,[0.25, 0.5, 0.75])
	vecQ2 = quantile(cellNumSpikes{2}/720,[0.25, 0.5, 0.75])
	%mean(cellNumSpikes{1}),std(cellNumSpikes{1}),mean(cellNumSpikes{2}),std(cellNumSpikes{2})
	%fprintf('Stim-mod, mean=%.0f,sd=%.0f; rand, mean=%.0f,sd=%.0f\n',a)
end

