
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
		sCombo.(strType).matR2_LinNorm(intFile,:) = [sData.dblR2_LinNorm];
		sCombo.(strType).matR2_LogNorm(intFile,:) = [sData.dblR2_LogNorm];
		sCombo.(strType).matR2_Gain(intFile,:) = [sData.dblR2_Gain];
		sCombo.(strType).matR2_GainFull(intFile,:) = [sData.dblR2_GainFull];
		sCombo.(strType).matR2_LinNorm_Single(intFile,:) = [sData.dblR2_LinNorm_Single];
		sCombo.(strType).matR2_LogNorm_Single(intFile,:) = [sData.dblR2_LogNorm_Single];
		sCombo.(strType).matR2_Gain_Single(intFile,:) = [sData.dblR2_Gain_Single];
		sCombo.(strType).matR2_GainFull_Single(intFile,:) = [sData.dblR2_GainFull_Single];
		
		sCombo.(strType).matAIC_LinNorm(intFile,:) = [sData.dblAIC_LinNorm];
		sCombo.(strType).matAIC_LogNorm(intFile,:) = [sData.dblAIC_LogNorm];
		sCombo.(strType).matAIC_Gain(intFile,:) = [sData.dblAIC_Gain];
		sCombo.(strType).matAIC_GainFull(intFile,:) = [sData.dblAIC_GainFull];
		sCombo.(strType).matAIC_LinNorm_Single(intFile,:) = [sData.dblAIC_LinNorm_Single];
		sCombo.(strType).matAIC_LogNorm_Single(intFile,:) = [sData.dblAIC_LogNorm_Single];
		sCombo.(strType).matAIC_Gain_Single(intFile,:) = [sData.dblAIC_Gain_Single];
		sCombo.(strType).matAIC_GainFull_Single(intFile,:) = [sData.dblAIC_GainFull_Single];
		
		sCombo.(strType).matN(intFile,:) = [sData.intN];
		sCombo.(strType).matN_Single(intFile,:) = [sData.intN_Single];
		sCombo.(strType).vecNumN(intFile) = sData(1).intNumN;
	end
end

%% plot R^2
figure;maxfig;
cellFields = fieldnames(sCombo.(strType));
cellModels = cellFields(contains(cellFields,'matR2_'));
intModelNum = numel(cellModels)/2;
matCol=lines(intModelNum);
for intType=1:numel(cellTypes)
	subplot(2,3,intType);hold on
	strType = cellTypes{intType};
	cellLegend = {};
	for i=1:intModelNum
		strField = cellModels{i};
		strName = strField(7:end);
		cellLegend{i} = strName;
		matData = sCombo.(strType).(strField);
		matData(matData<0)=0;
		errorbar(1:numel(vecTimescales),mean(matData),std(matData)./sqrt(intRecNum),'color',matCol(i,:));
	end
	set(gca,'xtick',1:numel(vecTimescales),'xticklabel',vecTimescales);
	xlabel('Timescale (s)');
	ylabel('Spike count distro R^2');
	title(sprintf('%s - Pop',strType),'interpreter','none');
	legend(cellLegend,'location','best','interpreter','none');
	ylim([0 1]);
end
fixfig;

drawnow;
export_fig(fullpath(strFigurePath,sprintf('Q3_GainModelR2_Pop.tif')));
export_fig(fullpath(strFigurePath,sprintf('Q3_GainModelR2_Pop.pdf')));

figure;maxfig;
for intType=1:numel(cellTypes)
	subplot(2,3,intType);hold on
	strType = cellTypes{intType};
	cellLegend = {};
	for i=1:intModelNum
		strField = cellModels{i+intModelNum};
		strName = strField(7:end);
		cellLegend{i} = strName;
		matData = sCombo.(strType).(strField);
		matData(matData<0)=0;
		
		errorbar(1:numel(vecTimescales),mean(matData),std(matData)./sqrt(intRecNum),'color',matCol(i,:));
	end
	set(gca,'xtick',1:numel(vecTimescales),'xticklabel',vecTimescales);
	xlabel('Timescale (s)');
	ylabel('Spike count distro R^2');
	title(sprintf('%s - Single',strType),'interpreter','none');
	legend(cellLegend,'location','best','interpreter','none');
	ylim([0 1]);
end
fixfig;

drawnow;
export_fig(fullpath(strFigurePath,sprintf('Q3_GainModelR2_Single.tif')));
export_fig(fullpath(strFigurePath,sprintf('Q3_GainModelR2_Single.pdf')));

%% plot AIC
%compare
figure;maxfig;
for intType=1:numel(cellTypes)
	subplot(2,3,intType);hold on
	strType = cellTypes{intType};
	cellLegend = {};
	for i=1:intModelNum
		strField = strrep(cellModels{i},'R2','AIC');
		strName = strField(8:end);
		cellLegend{i} = strName;
		matData = sCombo.(strType).(strField);
		%matData(matData<0)=0;
		
		errorbar(1:numel(vecTimescales),mean(matData),std(matData)./sqrt(intRecNum),'color',matCol(i,:));
	end
	set(gca,'xtick',1:numel(vecTimescales),'xticklabel',vecTimescales);
	xlabel('Timescale (s)');
	ylabel('Spike count distro AIC');
	title(sprintf('%s - Pop',strType),'interpreter','none');
	legend(cellLegend,'location','best','interpreter','none');
	%ylim([0 1]);
end
fixfig;

drawnow;
export_fig(fullpath(strFigurePath,sprintf('Q3_GainModelAIC_Pop.tif')));
export_fig(fullpath(strFigurePath,sprintf('Q3_GainModelAIC_Pop.pdf')));

figure;maxfig;
for intType=1:numel(cellTypes)
	subplot(2,3,intType);hold on
	strType = cellTypes{intType};
	cellLegend = {};
	for i=1:intModelNum
		strField = cellModels{i+intModelNum};
		strField = strrep(strField,'R2','AIC');
		strName = strField(8:end);
		cellLegend{i} = strName;
		matData = sCombo.(strType).(strField);
		%matData(matData<0)=0;
		
		errorbar(1:numel(vecTimescales),mean(matData),std(matData)./sqrt(intRecNum),'color',matCol(i,:));
	end
	set(gca,'xtick',1:numel(vecTimescales),'xticklabel',vecTimescales);
	xlabel('Timescale (s)');
	ylabel('Spike count distro AIC');
	title(sprintf('%s - Pop',strType),'interpreter','none');
	legend(cellLegend,'location','best','interpreter','none');
	%ylim([0 1]);
end
fixfig;

drawnow;
export_fig(fullpath(strFigurePath,sprintf('Q3_GainModelAIC_Single.tif')));
export_fig(fullpath(strFigurePath,sprintf('Q3_GainModelAIC_Single.pdf')));

%% plot normalized AIC for real pop
%compare
figure;
hold on
matAIC = [];
for i=1:intModelNum
	strField = strrep(cellModels{i},'R2','AIC');
	strName = strField(8:end);
	cellLegend{i} = strName;
	matAIC(:,:,i) = sCombo.(strType).(strField);
	%matData(matData<0)=0;
end

%normalize to worst model for each timescale
%(AICmin ? AICi)/2 is the relative log-likelihood
matRelAIC = matAIC;
for intScale=1:size(matAIC,2)
	matThisAIC = matAIC(:,intScale,:);
	%matRelAIC(:,intScale,:) = (matThisAIC - min(matThisAIC(:)))./2;
	matRelAIC(:,intScale,:) = -zscore(matThisAIC,[],'all');
end
	
for i=1:intModelNum
	matThisRelAIC = matRelAIC(:,:,i);
	strField = strrep(cellModels{i},'R2','AIC');
	strName = strField(8:end);
	errorbar(1:numel(vecTimescales),mean(matThisRelAIC),std(matThisRelAIC)./sqrt(intRecNum),'color',matCol(i,:));
end
set(gca,'xtick',1:numel(vecTimescales),'xticklabel',vecTimescales);
xlabel('Timescale (s)');
ylabel(sprintf('Model evidence (z-score) \n(normalized -AIC per timescale)'));
title(sprintf('Real - Pop'),'interpreter','none');fixfig;
legend(cellLegend,'location','best','interpreter','none');
%set(gca,'yscale','log');

fixfig;

drawnow;
export_fig(fullpath(strFigurePath,sprintf('Q3_GainModelEvidence.tif')));
export_fig(fullpath(strFigurePath,sprintf('Q3_GainModelEvidence.pdf')));
