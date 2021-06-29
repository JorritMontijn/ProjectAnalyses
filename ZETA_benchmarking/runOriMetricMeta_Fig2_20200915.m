clear all;
%close all;
strPath = 'F:\Data\Processed\ZETA\Inclusion\';
strFigPath = 'F:\Data\Results\ZETA\Inclusion\';
cellUniqueAreas = {...
	'V1',...Area 1
	'SC',...Area 2
	...'Poisson',...Area 3
	'Retina',...Area 4
	...%'CaNM',...Area 5
	'CaDG',...Area 6
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

cellRunStim = {...
	'',...Stim 1
	'RunDriftingGratings',...Stim 2
	%'RunNaturalMovie'...Stim 3
	};
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
strArea = 'Primary visual';%cellUniqueAreas{intArea}; %V1, SC, Retina, Poisson, GCaMP
strStim = 'RunDriftingGratings';
strName = replace([strArea strStim],cellRepStr(:,1),cellRepStr(:,2));
matZetaP = [];
matHzP = [];
for intRandType=1:2
	%set var
	strRand = cellRunRand{intRandType};
	
	%{
	%% load data
	strRunType = [strArea strRand strStim];
	sDir=dir([strPath 'ZetaDataMSD' strRunType 'Resamp100*']);
	intFiles=numel(sDir);
	for intFile=1:intFiles
		strFile = sDir(intFile).name;
		intResampNum = str2double(getFlankedBy(strFile,'Resamp','.mat'));
		sLoad=load([strPath strFile]);
		vecZeta = abs(sLoad.vecZeta);
		matZetaP(intRandType,:)=1-(normcdf(abs(vecZeta))-normcdf(-abs(vecZeta)));
	end
	%}
	%% load data
	strRunType = [strArea strRand strStim];
	sDir=dir([strPath 'ZetaData2MSD' strRunType 'Resamp100*']);
	intFiles=numel(sDir);
	for intFile=1:intFiles
		strFile = sDir(intFile).name;
		intResampNum = str2double(getFlankedBy(strFile,'Resamp','.mat'));
		sLoad=load([strPath strFile]);
		matHzP(intRandType,:) = sLoad.vecHzP;
		matHzZ(intRandType,:) = -norminv(sLoad.vecHzP/2);
		matZetaP(intRandType,:)= sLoad.vecP;
		matZetaZ(intRandType,:)= abs(sLoad.vecZeta);
		
	end
end


%% load data
strRunStim = 'RunDriftingGratings';
[sAggStim,sAggNeuron]=loadDataNpx(strArea,strRunStim);
cellRecIdx = {sAggStim.Rec};
intNeurons = numel(sAggNeuron);
cellDate = cellfun(@strrep,{sAggNeuron.Date},cellfill('-',size({sAggNeuron.Date})),cellfill('',size({sAggNeuron.Date})),'uniformoutput',false);
cellClust = cellfun(@num2str,{sAggNeuron.Cluster},'uniformoutput',false);
cellId = cellfun(@strcat,cellDate,cellfill('U',size({sAggNeuron.Date})),cellClust,'uniformoutput',false);
%strFind1 = '20200306U466';
%%
intChoose1 = 20;%7,11,17,20
intChoose2 = 3;

vecEx1 = find(matZetaP(1,:)<0.05 & matHzP(1,:)<0.05);
cellFind = {'',''};
cellFind{1} = cellId{vecEx1(intChoose1)};

%strFind2 = '20200306U458';
vecEx2 = find(matZetaP(1,:)<0.05 & matHzP(1,:)>0.05);
cellFind{2} = cellId{vecEx2(intChoose2)}
vecNeurons = [find(strcmp(cellId,cellFind{1})) find(strcmp(cellId,cellFind{2}))];


% get neuronal data
figure
h1=subplot(2,3,2);
hold on;
h2=subplot(2,3,5);
hold on;

intNeuronC = 0;
for intNeuron=vecNeurons
	intNeuronC = intNeuronC + 1;
	sThisNeuron = sAggNeuron(intNeuron);
	vecSpikeTimes = sThisNeuron.SpikeTimes;
	strRecIdx = sThisNeuron.Rec;
	strMouse = sThisNeuron.Mouse;
	strBlock = '';
	strDate = sThisNeuron.Date;
	intSU = sThisNeuron.Cluster;
	intClust = sThisNeuron.IdxClust;
	strNeuron = cellFind{intNeuronC};
	
	%% get matching recording data
	sThisRec = sAggStim(strcmpi(strRecIdx,cellRecIdx));
	vecStimOnTime = [];
	vecStimOffTime = [];
	for intRec=1:numel(sThisRec.cellStim)
		vecStimOnTime = cat(2,vecStimOnTime,sThisRec.cellStim{intRec}.structEP.vecStimOnTime);
		vecStimOffTime = cat(2,vecStimOffTime,sThisRec.cellStim{intRec}.structEP.vecStimOffTime);
	end
	
	vecTrialStarts = [];
	vecTrialStarts(:,1) = vecStimOnTime;
	vecTrialStarts(:,2) = vecStimOffTime;
	
	%% run zeta
	%get trial dur
	for intRand=1:2
		if intRand==2
			strRand = '-Rand';
		else
			strRand = '';
		end
		strRunType = [strName strRand strStim];
		
		dblUseMaxDur = round(median(diff(vecTrialStarts(:,1)))*2)/2;
		%set derivative params
		if contains(strRunType,'Rand')
			dblDur = dblUseMaxDur;
			vecJitter = 4*dblDur*rand([numel(vecStimOnTime) 1])-dblDur;
			matEventTimes = bsxfun(@plus,vecTrialStarts,vecJitter);
		else
			matEventTimes = vecTrialStarts;
		end
		
		%if size(matEventTimes,1) > 0,continue;end
		vecRestrictRange=[0 inf];
		intMakePlots = 0;
		intResamp = 100;
		[dblZetaP,vecLatencies,sZETA] = getZeta(vecSpikeTimes,matEventTimes,dblUseMaxDur,intResamp,intMakePlots,0,vecRestrictRange);
		
		%update matrices
		matZetaP(intRand,intNeuron) = dblZetaP;
		matHzP(intRand,intNeuron) = sZETA.dblMeanP;
		matZetaZ(intRand,intNeuron) = abs(sZETA.dblZETA);
		matHzZ(intRand,intNeuron) = -norminv(sZETA.dblMeanP/2);
		
		if intRand == 1
			axes(h1);
		else
			axes(h2);
		end
		scatter(matHzZ(intRand,intNeuron),matZetaZ(intRand,intNeuron),300,[0 0 0],'x');
		
		
		% plot
		subplot(4,6,5+6*(intNeuronC-1)+12*(intRand-1))
		%bin
		sOpt = struct;
		sOpt.handleFig =-1;
		if dblUseMaxDur < 0.5
			dblBinSize = dblUseMaxDur/20;
		else
			dblBinSize = 0.025;
		end
		vecBins = 0:dblBinSize:dblUseMaxDur;
		[vecMean,vecSEM,vecWindowBinCenters] = doPEP(vecSpikeTimes,vecBins,matEventTimes(:,1),sOpt);
		errorbar(vecWindowBinCenters,vecMean,vecSEM);
		ylim([0 max(get(gca,'ylim'))]);
		xlabel('Time from event (s)');
		ylabel('Mean spiking rate (Hz)');
		title(strNeuron);
		fixfig
		
		subplot(4,6,6+6*(intNeuronC-1)+12*(intRand-1))
		%zeta
		cla;
		hold all
		intPlotIters = 50;
		for intOffset=1:intPlotIters
			plot(sZETA.vecSpikeT,sZETA.matRandD(:,intOffset),'Color',[0.5 0.5 0.5]);
		end
		plot(sZETA.vecSpikeT,sZETA.vecD,'Color',lines(1));
		hold off
		xlabel('Time from event (s)');
		ylabel('Offset of data from linear (s)');
		if intRand==2
			title(sprintf('Z,p=%.3f, T,p=%.3f)',dblZetaP,sZETA.dblMeanP));
		else
			title(sprintf('Z,p=%.1e, T,p=%.1e)',dblZetaP,sZETA.dblMeanP));
		end
		fixfig
	end
end
%
matHzP(matHzZ(:)==0) = 1e-29;
matZetaP(matZetaP(:)==0) = 1e-29;
axes(h1);
matC = [0.5 0.5 0.5;...
	0 0.8 0;...
	0.8 0 0;...
	0 0 0.8];
vecColor1 = 1 + (matZetaP(1,:) < 0.05 & matHzP(1,:) > 0.05) + 2*(matZetaP(1,:) > 0.05 & matHzP(1,:) < 0.05) + 3*(matZetaP(1,:) < 0.05 & matHzP(1,:) < 0.05);
scatter(matHzZ(1,:),matZetaZ(1,:),500,vecColor1,'.');
colormap(h1,matC(1:max(vecColor1),:));
%xlim([0 1]);ylim([0 1]);
xlabel('Z-statistic rate-based t-test (\Phi^-^1(1-p/2))')
ylabel('ZETA (\zeta_c)')
title(sprintf('Inclusion: %s=%.3f, %s=%.3f',getGreek('zeta'),sum(matZetaP(1,:)<0.05)/numel(matZetaP(1,:)),getGreek('mu'),sum(matHzP(1,:)<0.05)/numel(matHzP(1,:))))
%set(gca,'xscale','log','yscale','log');
fixfig;

axes(h2)
vecColor2 = 1 + 1*(matZetaP(2,:) > 0.05 & matHzP(2,:) < 0.05) + 2*(matZetaP(2,:) < 0.05 & matHzP(2,:) > 0.05) + 3*(matZetaP(2,:) < 0.05 & matHzP(2,:) < 0.05);
scatter(matHzZ(2,:),matZetaZ(2,:),500,vecColor1,'.');
colormap(h2,matC(1:max(vecColor1),:));
%xlim([0 1]);ylim([0 1]);
xlabel('Z-statistic rate-based t-test (\Phi^-^1(1-p/2))')
ylabel('ZETA (\zeta_c)')
title(sprintf('False alarms: %s=%.3f, %s=%.3f',getGreek('zeta'),sum(matZetaP(2,:)<0.05)/numel(matZetaP(2,:)),getGreek('mu'),sum(matHzP(2,:)<0.05)/numel(matHzP(2,:))))
%set(gca,'xscale','log','yscale','log');
fixfig;maxfig;
return
%%
export_fig(sprintf('Fig2_%s.tif',getDate));
export_fig(sprintf('Fig2_%s.pdf',getDate));