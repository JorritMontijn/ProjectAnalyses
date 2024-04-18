%% aim
%{
https://www.nature.com/articles/s41598-021-95037-z
https://www.nature.com/articles/s42003-021-02437-y

q1: Progression of orientation information over time: how does trial activity evolve? what is the function of the onset peak?
> is decoding better when matched for stimulus phase? => no.

q2: "Spike time and rate coding can be represented within a single model of spiking probability as a
function of time: rate codes are uniform over a certain period tau, while spike time codes are
temporally localized peaks"
> is this true?

q3: Rate codes do not exist; a rate code is simply a subset of spike time codes where the temporal
integration window is very large. But what about multi dim codes? Those are all rate based. Can we
formulate a multidimensional spike-time code? I.e., can we make a taxonomy of neural codes?

q4: How does information evolve over time, is initial peak indeed less tuned? Is pop activity
rhythmic? Are stimuli encoded invariant to brain state? Eg, high arousal, low arousal. Or is
stimulus manifold dynamic over time? Does manifold scale with arousal? => How does manifold depend
on binning size? What is the optimal time window?

%}
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
boolSingleRecPlots = false;
boolUseMax = true;
dblRemOnset = 0; %remove onset period in seconds; 0.125 for sim, 0.25 for npx

%% onset string
if dblRemOnset == 0
	strOnset = '';
else
	strOnset = sprintf('%.2f',dblRemOnset);
end
sFiles = dir ([strTargetDataPath 'Q2bData*' strOnset '.mat']);
intRecNum = numel(sFiles);

%% pre-allocate
figure;maxfig;
cellTypes =  {'Real','Poiss','Shuff'};
intTypeNum = numel(cellTypes);
vecH = nan([1 intTypeNum]);
vecH2 = nan([1 intTypeNum]);
for intType=1:intTypeNum
	vecH(intType) = subplot(2,intTypeNum,intType);hold on;
	vecH2(intType) = subplot(2,intTypeNum,intType+intTypeNum);hold on;
end


cellSlopes = {};
cellR2 = {};
cellR2Exp = {};
cellExpHalflife = {};
cellExpScale = {};
cellExpAsymptote = {};
matCol=lines(numel(cellTypes));
for intFile=1:intRecNum
	%% load
	sLoad=load(fullpath(sFiles(intFile).folder,sFiles(intFile).name));
	sAggData = sLoad.sAggData;
	strRec = getFlankedBy(sFiles(intFile).name,'Q2bData','_g0_t0');
	cellTheseTypes = {sAggData.strType};
	vecTimescales = sAggData(1).vecTimescales;
	
	%% plot & save data
	for intType=1:numel(cellTypes)
		strType = cellTypes{intType};
		intUseEntry = find(strcmp(cellTheseTypes,strType));
		if isempty(intUseEntry)
			warning(sprintf('   Error on %s-%d, no %s',strRec,intType,strType));
			continue;
		end
		
		%% save data
		cellSlopes{intFile,intType} = sAggData(intType).vecSlopes;
		cellR2{intFile,intType} = sAggData(intType).vecR2;
		cellR2Exp{intFile,intType} = sAggData(intType).vecR2Exp_Time;
		cellExpHalflife{intFile,intType} = sAggData(intType).vecHalfLifeExp_Time;
		cellExpScale{intFile,intType} = sAggData(intType).vecScaleExp_Time;
		cellExpAsymptote{intFile,intType} = sAggData(intType).vecAsymptoteExp_Time;
		
		%% plot sd/mean
		plot(vecH(intType),sAggData(intUseEntry).vecMean,sAggData(intUseEntry).vecSd,'color',matCol(intType,:));
		
		%% plot cv/timescale
		plot(vecH2(intType),sAggData(intUseEntry).vecTimescales,sAggData(intUseEntry).vecCV,'color',matCol(intType,:));
		
		
	end
end
for intType=1:numel(cellTypes)
	%finish plots 1; mean and labels
	xlabel(vecH(intType),'Mean of spike counts');
	ylabel(vecH(intType),'Sd of spike counts');
	title(vecH(intType),sprintf('%s',cellTypes{intType}),'interpreter','none');
	
	%finish plots 2; timescale and cv
	xlabel(vecH2(intType),'Timescale (s)');
	ylabel(vecH2(intType),'CV (sd/mu)');
	title(vecH2(intType),sprintf('%s',cellTypes{intType}),'interpreter','none');
	set(vecH2(intType),'ylim',[0 1.2]);%[0 max(get(vecH2(intType),'ylim'))]);
end
fixfig;

drawnow;
export_fig(fullpath(strFigurePath,sprintf('Q2b_CV_Timescales.tif')));
export_fig(fullpath(strFigurePath,sprintf('Q2b_CV_Timescales.pdf')));

%% plot slopes/r2 for sd/mean and half-life/r2 for cv/timescale
%use same pop size for all, or max for each?
if boolUseMax
	vecUseEntries = cellfun(@numel,cellSlopes(:,1));
else
	vecUseEntries = ones(1,intRecNum)*min(cellfun(@numel,cellSlopes(:,1)));
end
matSlopes = nan(intRecNum,intTypeNum);
matR2 = nan(intRecNum,intTypeNum);
matExpR2 = nan(intRecNum,intTypeNum);
matExpHalfLife = nan(intRecNum,intTypeNum);
matExpAsymptote = nan(intRecNum,intTypeNum);
matExpScale = nan(intRecNum,intTypeNum);
for intType=1:intTypeNum
	for intRec=1:intRecNum
		matSlopes(intRec,intType) = cellSlopes{intRec,intType}(vecUseEntries(intRec));
		matR2(intRec,intType) = cellR2{intRec,intType}(vecUseEntries(intRec));
		matExpR2(intRec,intType) = cellR2Exp{intRec,intType}(vecUseEntries(intRec));
		matExpHalfLife(intRec,intType) = cellExpHalflife{intRec,intType}(vecUseEntries(intRec));
		matExpAsymptote(intRec,intType) = cellExpAsymptote{intRec,intType}(vecUseEntries(intRec));
		matExpScale(intRec,intType) = cellExpScale{intRec,intType}(vecUseEntries(intRec));
	end
end

% plot
figure;maxfig;
vecH3 = nan([1 6]);
cellPlotType = {'R^2 linear fit sd/mu of spike counts per bin',...
	'Slope sd/mu of spike counts per bin',...
	'R^2 exp decay fit CV/timescale',...
	'Half-life CV/timescale',...
	'Asymptote CV/timescale',...
	'Scale CV/timescale'};
for intPlot=1:2
	vecH3(intPlot) = subplot(2,4,intPlot);hold on;
	title(cellPlotType{intPlot});
end
for intPlot=3:6
	vecH3(intPlot) = subplot(2,4,intPlot+2);hold on;
	title(cellPlotType{intPlot});
end

for intType=1:intTypeNum
%slope sd/mu
%swarmchart(vecH3(1),intType*ones(1,intRecNum),matR2(:,intType),'jitterwidth',0.5);
errorbar(vecH3(1),intType,mean(matR2(:,intType)),std(matR2(:,intType))./sqrt(intRecNum),'x','capsize',20);

%linearity sd/mu
%swarmchart(vecH3(2),intType*ones(1,intRecNum),matSlopes(:,intType),'jitterwidth',0.5);
errorbar(vecH3(2),intType,mean(matSlopes(:,intType)),std(matSlopes(:,intType))./sqrt(intRecNum),'x','capsize',20);

%half-life cv/time
%swarmchart(vecH3(3),intType*ones(1,intRecNum),matExpR2(:,intType),'jitterwidth',0.5);
errorbar(vecH3(3),intType,mean(matExpR2(:,intType)),std(matExpR2(:,intType))./sqrt(intRecNum),'x','capsize',20);

%R^2 cv/time
%swarmchart(vecH3(4),intType*ones(1,intRecNum),matExpLambda(:,intType),'jitterwidth',0.5);
errorbar(vecH3(4),intType,mean(matExpHalfLife(:,intType)),std(matExpHalfLife(:,intType))./sqrt(intRecNum),'x','capsize',20);

errorbar(vecH3(5),intType,mean(matExpAsymptote(:,intType)),std(matExpAsymptote(:,intType))./sqrt(intRecNum),'x','capsize',20);

errorbar(vecH3(6),intType,mean(matExpScale(:,intType)),std(matExpScale(:,intType))./sqrt(intRecNum),'x','capsize',20);

end

%tests
[h,p12]=ttest(matExpR2(:,1),matExpR2(:,2));
[h,p13]=ttest(matExpR2(:,1),matExpR2(:,3));
[h,p23]=ttest(matExpR2(:,2),matExpR2(:,3));

subplot(2,4,4)
title(sprintf('R^2 exp decay p, R-P=%.2e;R-S=%.2e;S-P=%.2e',p12,p13,p23));
axis off

%finish figs
ylabel(vecH3(1),'R^2 lin fit sd/mu');
set(vecH3(1),'xtick',1:intTypeNum,'xticklabel',cellTypes,'xlim',[0.5 numel(cellTypes)+0.5]);
ylabel(vecH3(2),'Slope lin fit (sd/mu)');
set(vecH3(2),'xtick',1:intTypeNum,'xticklabel',cellTypes,'xlim',[0.5 numel(cellTypes)+0.5]);
ylabel(vecH3(3),'R^2 exponential decay (CV/s)');
set(vecH3(3),'xtick',1:intTypeNum,'xticklabel',cellTypes,'xlim',[0.5 numel(cellTypes)+0.5]);
ylabel(vecH3(4),'Half-life of CV/timescale (\lambda_1_/_2)');
set(vecH3(4),'xtick',1:intTypeNum,'xticklabel',cellTypes,'xlim',[0.5 numel(cellTypes)+0.5]);
ylabel(vecH3(5),'Asymptote of CV/timescale');
set(vecH3(5),'xtick',1:intTypeNum,'xticklabel',cellTypes,'xlim',[0.5 numel(cellTypes)+0.5]);
ylabel(vecH3(6),'Scale of CV/timescale');
set(vecH3(6),'xtick',1:intTypeNum,'xticklabel',cellTypes,'xlim',[0.5 numel(cellTypes)+0.5]);
fixfig;

drawnow;
export_fig(fullpath(strFigurePath,sprintf('Q2b_CV_ExpDecay.tif')));
export_fig(fullpath(strFigurePath,sprintf('Q2b_CV_ExpDecay.pdf')));
