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
strTag = 'Q2cData';
sFiles = dir ([strTargetDataPath strTag '*' strOnset '.mat']);
intRecNum = numel(sFiles);

%% pre-allocate
figure;maxfig;
cellTypes =  {'Real','Poiss','Shuff','PoissGain','ShuffTid','Uniform'};
intTypeNum = numel(cellTypes);
vecH = nan([1 intTypeNum]);
vecH1_5 = nan([1 intTypeNum]);
vecH2 = nan([1 intTypeNum]);
for intType=1:intTypeNum
	vecH(intType) = subplot(3,intTypeNum,intType);hold on;
	vecH1_5(intType) = subplot(3,intTypeNum,intType+intTypeNum);hold on;
	vecH2(intType) = subplot(3,intTypeNum,intType+2*intTypeNum);hold on;
end


cellSlopes = {};
cellR2 = {};
cellExpR2 = {};
cellExpHalflife = {};
cellExpScale = {};
cellExpAsymptote = {};
matCol=lines(numel(cellTypes));
matMaxSlopes = [];
for intFile=1:intRecNum
	%% load
	sLoad=load(fullpath(sFiles(intFile).folder,sFiles(intFile).name));
	sAggData = sLoad.sAggData;
	strRec = getFlankedBy(sFiles(intFile).name,strTag,'_g0_t0');
	cellTheseTypes = {sAggData.strType};
	vecTimescales = sAggData(1).vecTimescales;
	vecJitter = sAggData(1).vecJitter;
		
	%% plot & save data
	for intType=1:numel(cellTypes)
		strType = cellTypes{intType};
		intUseEntry = find(strcmp(cellTheseTypes,strType));
		if isempty(intUseEntry)
			warning(sprintf('   Error on %s-%d, no %s',strRec,intType,strType));
			continue;
		end
		
		%% save data
		cellExpR2{intFile,intType} = sAggData(intUseEntry).vecR2_Exp;
		cellExpHalflife{intFile,intType} = sAggData(intUseEntry).vecHalfLife_Exp;
		cellExpScale{intFile,intType} = sAggData(intUseEntry).vecScale_Exp;
		cellExpAsymptote{intFile,intType} = sAggData(intUseEntry).vecAsymptote_Exp;
		
		
		cellRootR2{intFile,intType} = sAggData(intUseEntry).vecR2_Root;
		cellRootAsymptote{intFile,intType} = sAggData(intUseEntry).vecAsymptote_Root;
		cellRootScale{intFile,intType} = sAggData(intUseEntry).vecScale_Root;
		cellRootExponent{intFile,intType} = sAggData(intUseEntry).vecExponent_Root;
		
		matMean = sAggData(intUseEntry).matMean;
		matSd = sAggData(intUseEntry).matSd;
		matCV = sAggData(intUseEntry).matCV;
		
	end
end

for intType=1:numel(cellTypes)
	%finish plots 1; mean and sd
	xlabel(vecH(intType),'Mean of spike counts');
	ylabel(vecH(intType),'Sd of spike counts');
	title(vecH(intType),sprintf('%s',cellTypes{intType}),'interpreter','none');
	
	%finish plots 1.5; mean and var
	xlabel(vecH1_5(intType),'Mean of spike counts');
	ylabel(vecH1_5(intType),'Var of spike counts');
	title(vecH1_5(intType),sprintf('%s',cellTypes{intType}),'interpreter','none');
	
	%finish plots 2; timescale and cv
	xlabel(vecH2(intType),'Timescale (s)');
	ylabel(vecH2(intType),'CV (sd/mu)');
	title(vecH2(intType),sprintf('%s',cellTypes{intType}),'interpreter','none');
	set(vecH2(intType),'ylim',[0 1.2]);%[0 max(get(vecH2(intType),'ylim'))]);
end
fixfig;

drawnow;
export_fig(fullpath(strFigurePath,sprintf('Q2c_CV_Timescales.tif')));
export_fig(fullpath(strFigurePath,sprintf('Q2c_CV_Timescales.pdf')));


%% plot slopes/r2 for sd/mean and half-life/r2 for cv/timescale
%use same pop size for all, or max for each?
vecUseEntries = cellfun(@numel,cellSlopes(:,1));

matExpR2 = nan(intRecNum,intTypeNum);
matExpHalfLife = nan(intRecNum,intTypeNum);
matExpAsymptote = nan(intRecNum,intTypeNum);
matExpScale = nan(intRecNum,intTypeNum);
for intType=1:intTypeNum
	for intRec=1:intRecNum
		matExpR2(intRec,intType) = cellExpR2{intRec,intType}(vecUseEntries(intRec));
		matExpHalfLife(intRec,intType) = cellExpHalflife{intRec,intType}(vecUseEntries(intRec));
		matExpAsymptote(intRec,intType) = cellExpAsymptote{intRec,intType}(vecUseEntries(intRec));
		matExpScale(intRec,intType) = cellExpScale{intRec,intType}(vecUseEntries(intRec));
		
		matRootR2(intRec,intType) = cellRootR2{intRec,intType}(vecUseEntries(intRec));
		matRootExponent(intRec,intType) = cellRootExponent{intRec,intType}(vecUseEntries(intRec));
		matRootAsymptote(intRec,intType) = cellRootAsymptote{intRec,intType}(vecUseEntries(intRec));
		matRootScale(intRec,intType) = cellRootScale{intRec,intType}(vecUseEntries(intRec));
	end
end

% plot
figure;maxfig;
cellH3 = cell(1,10);
cellPlotType = {'R^2 linear fit sd/mu of spike counts per bin',...
	'Slope sd/mu of spike counts per bin',...
	'R^2 exp decay fit CV/timescale',...
	'Half-life exp CV/timescale',...
	'Asymptote exp CV/timescale',...
	'Scale CV/timescale',...
	'R^2 root fit CV/timescale',...
	'Exponent root CV/timescale',...
	'Asymptote root CV/timescale',...
	'Scale root CV/timescale'};
for intPlot=1:2
	cellH3{intPlot} = subplot(3,4,intPlot);hold on;
	title(cellPlotType{intPlot});
end
for intPlot=3:10
	cellH3{intPlot} = subplot(3,4,intPlot+2);hold on;
	title(cellPlotType{intPlot});
end

for intType=1:intTypeNum
%slope sd/mu
%swarmchart(vecH3(1),intType*ones(1,intRecNum),matR2(:,intType),'jitterwidth',0.5);
errorbar(cellH3{1},intType,mean(matR2(:,intType)),std(matR2(:,intType))./sqrt(intRecNum),'x','capsize',20);

%linearity sd/mu
%swarmchart(vecH3(2),intType*ones(1,intRecNum),matSlopes(:,intType),'jitterwidth',0.5);
errorbar(cellH3{2},intType,mean(matSlopes(:,intType)),std(matSlopes(:,intType))./sqrt(intRecNum),'x','capsize',20);

%half-life cv/time
%swarmchart(vecH3(3),intType*ones(1,intRecNum),matExpR2(:,intType),'jitterwidth',0.5);
errorbar(cellH3{3},intType,mean(matExpR2(:,intType)),std(matExpR2(:,intType))./sqrt(intRecNum),'x','capsize',20);

%R^2 cv/time
%swarmchart(vecH3(4),intType*ones(1,intRecNum),matExpLambda(:,intType),'jitterwidth',0.5);
errorbar(cellH3{4},intType,mean(matExpHalfLife(:,intType)),std(matExpHalfLife(:,intType))./sqrt(intRecNum),'x','capsize',20);

errorbar(cellH3{5},intType,mean(matExpAsymptote(:,intType)),std(matExpAsymptote(:,intType))./sqrt(intRecNum),'x','capsize',20);

errorbar(cellH3{6},intType,mean(matExpScale(:,intType)),std(matExpScale(:,intType))./sqrt(intRecNum),'x','capsize',20);

%% root
errorbar(cellH3{7},intType,mean(matRootR2(:,intType)),std(matRootR2(:,intType))./sqrt(intRecNum),'x','capsize',20);

errorbar(cellH3{8},intType,mean(matRootExponent(:,intType)),std(matRootExponent(:,intType))./sqrt(intRecNum),'x','capsize',20);

errorbar(cellH3{9},intType,mean(matRootAsymptote(:,intType)),std(matRootAsymptote(:,intType))./sqrt(intRecNum),'x','capsize',20);

errorbar(cellH3{10},intType,mean(matRootScale(:,intType)),std(matRootScale(:,intType))./sqrt(intRecNum),'x','capsize',20);

end

%tests
[h,p12]=ttest(matExpR2(:,1),matExpR2(:,2));
[h,p13]=ttest(matExpR2(:,1),matExpR2(:,3));
[h,p23]=ttest(matExpR2(:,2),matExpR2(:,3));

subplot(3,4,4)
title(sprintf('R^2 exp decay p, R-P=%.2e;R-S=%.2e;S-P=%.2e',p12,p13,p23));
axis off

%finish figs
ylabel(cellH3{1},'R^2 lin fit sd/mu');
set(cellH3{1},'xtick',1:intTypeNum,'xticklabel',cellTypes,'xlim',[0.5 numel(cellTypes)+0.5]);
ylabel(cellH3{2},'Slope lin fit (sd/mu)');
set(cellH3{2},'xtick',1:intTypeNum,'xticklabel',cellTypes,'xlim',[0.5 numel(cellTypes)+0.5]);
ylabel(cellH3{3},'R^2 exponential decay (CV/s)');
set(cellH3{3},'xtick',1:intTypeNum,'xticklabel',cellTypes,'xlim',[0.5 numel(cellTypes)+0.5]);
ylabel(cellH3{4},'Half-life of CV/timescale (\lambda_1_/_2)');
set(cellH3{4},'xtick',1:intTypeNum,'xticklabel',cellTypes,'xlim',[0.5 numel(cellTypes)+0.5]);
ylabel(cellH3{5},'Asymptote of CV/timescale');
set(cellH3{5},'xtick',1:intTypeNum,'xticklabel',cellTypes,'xlim',[0.5 numel(cellTypes)+0.5]);
ylabel(cellH3{6},'Scale of CV/timescale');
set(cellH3{6},'xtick',1:intTypeNum,'xticklabel',cellTypes,'xlim',[0.5 numel(cellTypes)+0.5]);

ylabel(cellH3{7},'R^2 root fit (CV/s)');
set(cellH3{7},'xtick',1:intTypeNum,'xticklabel',cellTypes,'xlim',[0.5 numel(cellTypes)+0.5]);
ylabel(cellH3{8},'Exponent of root fit');
set(cellH3{8},'xtick',1:intTypeNum,'xticklabel',cellTypes,'xlim',[0.5 numel(cellTypes)+0.5]);
ylabel(cellH3{9},'Asymptote of root fit');
set(cellH3{9},'xtick',1:intTypeNum,'xticklabel',cellTypes,'xlim',[0.5 numel(cellTypes)+0.5]);
ylabel(cellH3{10},'Scale of root fit');
set(cellH3{10},'xtick',1:intTypeNum,'xticklabel',cellTypes,'xlim',[0.5 numel(cellTypes)+0.5]);
fixfig;

drawnow;
export_fig(fullpath(strFigurePath,sprintf('Q2c_CV_ExpDecay.tif')));
export_fig(fullpath(strFigurePath,sprintf('Q2c_CV_ExpDecay.pdf')));
