
%% define qualifying areas
clear all;%close all;
if isfolder('F:\Drive\PopTimeCoding') && isfolder('F:\Data\Processed\Neuropixels\')
	strDataPath = 'F:\Data\Processed\Neuropixels\';
	strFigurePathSR = 'F:\Drive\PopTimeCoding\single_recs';
	strFigurePath = 'F:\Drive\PopTimeCoding\figures\';
	strTargetDataPath = 'F:\Drive\PopTimeCoding\data\';
else
	strDataPath = 'E:\DataPreProcessed\';
	strFigurePathSR = 'C:\Drive\PopTimeCoding\single_recs';
	strFigurePath = 'C:\Drive\PopTimeCoding\figures\';
	strTargetDataPath = 'C:\Drive\PopTimeCoding\data\';
end

%% define parameters
%single rec plots
cellTypes =  {'Real','Poiss','ShuffTid','RandTid','RandTxClass','Uniform'};
boolSingleRecPlots = false;
dblRemOnset = 0; %remove onset period in seconds; 0.125 for sim, 0.25 for npx

%% onset string
if dblRemOnset == 0
	strOnset = '';
else
	strOnset = sprintf('%.2f',dblRemOnset);
end
strTag = 'Q3Data';
%strTargetDataPath = 'C:\Drive\PopTimeCoding\data\q2c\';
sFiles = dir ([strTargetDataPath strTag '*' strOnset '.mat']);
intRecNum = numel(sFiles);

%% pre-allocate
for intFile=1:intRecNum
	%% load
	sLoad=load(fullpath(sFiles(intFile).folder,sFiles(intFile).name));
	sAggData = sLoad.sAggData;
	strRec = getFlankedBy(sFiles(intFile).name,strTag,'_g0_t0');
	cellFields = fieldnames(sAggData);
	cellTheseTypes = cellfun(@(x) sAggData.(x).strType,cellFields,'uniformoutput',false);
	vecTimescales = [sAggData.(cellFields{1}).dblTimescale];
	
	%% plot & save data
	for intType=1:numel(cellTypes)
		strType = cellTypes{intType};
		intUseEntry = find(strcmp(cellTheseTypes,strType));
		if isempty(intUseEntry)
			warning(sprintf('   Error on %s-%d, no %s',strRec,intType,strType));
			continue;
		end
		sData = sAggData.(cellFields{intUseEntry});
		
		%% save data
		sCombo.(strType).matR2_Gauss(intFile,:) = [sData.dblR2_Gauss];
		sCombo.(strType).matR2_Gain(intFile,:) = [sData.dblR2_Gain];
		sCombo.(strType).matR2_GainFull(intFile,:) = [sData.dblR2_GainFull];
		sCombo.(strType).matR2_Gauss_Single(intFile,:) = [sData.dblR2_Gauss_Single];
		sCombo.(strType).matR2_Gain_Single(intFile,:) = [sData.dblR2_Gain_Single];
		sCombo.(strType).matR2_GainFull_Single(intFile,:) = [sData.dblR2_GainFull_Single];
	end
end

%% plot R^2
error switch plotting: plot gauss/gain/gainfull in single plot, plot real/poiss/etc in separate plots
make one plot for pop, one for single

figure;maxfig;
cellFields = fieldnames(sCombo.(strType));
matCol=lines(6);
for i=1:numel(cellFields)
	strField = cellFields{i};
	strName = strField(7:end);
	subplot(2,3,i);hold on
	for intType=1:numel(cellTypes)
		strType = cellTypes{intType};
		matData = sCombo.(strType).(strField);
		
		errorbar(1:numel(vecTimescales),mean(matData),std(matData)./sqrt(intRecNum),'color',matCol(intType,:));
	end
	set(gca,'xtick',1:numel(vecTimescales),'xticklabel',vecTimescales);
	xlabel('Timescale (s)');
	ylabel('Spike count distro R^2');
	title(sprintf('%s',strName),'interpreter','none');
	legend(cellTypes,'location','best');
	ylim([0 1]);
end
fixfig;
