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

intIdxNpx = 0;
intIdx = 1;
strArea = 'BurstyTest'; %V1, SC, Retina, Poisson, GCaMP

for intRandType=1:2
	%set var
	strRand = cellRunRand{intRandType};
	
	
	%% load data
	strRunType = [strArea strRand];
	sDir=dir([strPath 'ZetaBurst' strRunType 'Resamp300.mat']);
	intFiles=numel(sDir);
	for intFile=1:intFiles
		strFile = sDir(intFile).name;
		intResampNum = str2double(getFlankedBy(strFile,'Resamp','.mat'));
		sLoad=load([strPath strFile]);
		
		%sig at alpha=0.05
		dblAlpha = 0.05;
		sSignif.matBISIP_g(intIdx,intRandType) = sum(sLoad.vecBISIP_g<dblAlpha)/numel(sLoad.vecBISIP_g);
		sSignif.matBISIP_ks(intIdx,intRandType) = sum(sLoad.vecBISIP_ks<dblAlpha)/numel(sLoad.vecBISIP_ks);
		sSignif.matHzP(intIdx,intRandType) = sum(sLoad.vecHzP<dblAlpha)/numel(sLoad.vecHzP);
		sSignif.matISIP_IntG(intIdx,intRandType) = sum(sLoad.vecISIP_IntG<dblAlpha)/numel(sLoad.vecISIP_IntG);
		sSignif.matISIP_g(intIdx,intRandType) = sum(sLoad.vecISIP_g<dblAlpha)/numel(sLoad.vecISIP_g);
		sSignif.matISIP_ks(intIdx,intRandType) = sum(sLoad.vecISIP_ks<dblAlpha)/numel(sLoad.vecISIP_ks);
		sSignif.matPoissP(intIdx,intRandType) = sum(sLoad.vecPoissP<dblAlpha)/numel(sLoad.vecPoissP);
		sSignif.matZetaP(intIdx,intRandType) = sum(sLoad.vecZetaP<dblAlpha)/numel(sLoad.vecZetaP);
		
		%comput time
		sComput.cellComputTimeBISI{intIdx,intRandType} = sLoad.vecComputTimeBISI;
		sComput.cellComputTimePoiss{intIdx,intRandType} = sLoad.vecComputTimePoiss;
		sComput.cellComputTimeISI{intIdx,intRandType} = sLoad.vecComputTimeISI;
		sComput.cellComputTimeZETA{intIdx,intRandType} = sLoad.vecComputTimeZETA;
		
		%p-value vectors
		sP.Zeta{intIdx,intRandType} = sLoad.vecZetaP; %light blue
		sP.Z_ISI_g{intIdx,intRandType} = sLoad.vecBISIP_g; %dark blue
		sP.Ttest{intIdx,intRandType} = sLoad.vecHzP; %black
		sP.ISI_IntG{intIdx,intRandType} = sLoad.vecISIP_IntG; %purple
		sP.ISI_g{intIdx,intRandType} = sLoad.vecISIP_g; %red
		sP.Z_ISI_ks{intIdx,intRandType} = sLoad.vecBISIP_ks; %red-purple burgundy
		sP.ISI_ks{intIdx,intRandType} = sLoad.vecISIP_ks; %orange
		sP.Poiss{intIdx,intRandType} = sLoad.vecPoissP; %grey
		
		%other
		cellNumSpikes{intIdx,intRandType} = sLoad.vecNumSpikes;
	end
	
end

%% plot ROC
if ~isempty(sP) && size(sP.Zeta,1) >= intIdx && ~isempty(sP.Zeta{intIdx,1})
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
	%vecH(intIdxNpx) = subplot(4,3,intIdxNpx);
	subplot(2,2,1)
	maxfig;
	hold on;
	
	for intTest=1:intTests
		%intTest = intTest + 1;
		cellData = sP.(cellTests{intTest});
		
		intCells = numel(cellData{intIdx,1});
		vecBothData = cat(2,cellData{intIdx,1},cellData{intIdx,2});
		vecBothLabels = cat(2,zeros(size(cellData{intIdx,1})),ones(size(cellData{intIdx,2})));
		vecThresholds = sort(vecBothData);
		vecRealP = cellData{intIdx,1};
		vecShuffP = cellData{intIdx,2};
		
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
	
	%test AUC
	intZeta = 1;
	intTtest = 3;
	
	cellData = sP.(cellTests{intZeta});
	vecBothDataZ = cat(2,cellData{intIdx,1},cellData{intIdx,2});
	vecBothLabelsZ = cat(2,zeros(size(cellData{intIdx,1})),ones(size(cellData{intIdx,2})));
	[AUC_Z,AUC_Z_ci,AUC_Z_se] = auc(cat(1,vecBothLabelsZ,vecBothDataZ)',0.05,'mann-whitney');
	
	cellData = sP.(cellTests{intTtest});
	vecBothDataT = cat(2,cellData{intIdx,1},cellData{intIdx,2});
	vecBothLabelsT = cat(2,zeros(size(cellData{intIdx,1})),ones(size(cellData{intIdx,2})));
	[AUC_T,AUC_T_ci,AUC_T_se] = auc(cat(1,vecBothLabelsT,vecBothDataT)',0.05,'mann-whitney');
	
% Observed data
m0 = AUC_T - AUC_Z;
s0 = (AUC_T_se + AUC_Z_se)/2;
z = m0/s0;
AUC_p = 1 - abs(normcdf(z)-normcdf(-z));
return

	%{
		%shuffled p MIMI
		subplot(2,2,2)
		dblStep = 0.05;
		vecBinsP = 0:dblStep:1;
		vecBinsPlot = (dblStep/2):dblStep:(1-dblStep/2);
		vecCountsHzP = histcounts(cellHzP{intIdx,2},vecBinsP);
		vecCountsISIP = histcounts(cellISIP{intIdx,2},vecBinsP);
		vecCountsZETAP = histcounts(cellZetaP{intIdx,2},vecBinsP);
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
	vecComputTimeBISI = sComput.cellComputTimeBISI{intIdx,1};
	vecComputTimePoiss = sComput.cellComputTimePoiss{intIdx,1};
	vecComputTimeISI = sComput.cellComputTimeISI{intIdx,1};
	vecComputTimeZETA = sComput.cellComputTimeZETA{intIdx,1};
	
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

