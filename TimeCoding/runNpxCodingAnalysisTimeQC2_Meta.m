%% aim
%{
https://www.nature.com/articles/s41598-021-95037-z
https://www.nature.com/articles/s42003-021-02437-y

q1: can single-neuron variability be explained as population-level gain multiplication?
> estimate tuning curve from real data, then apply trial-by-trial gain to all neurons
>model predicting firing rate as combination of mean-rate multiplied by tuning curve => what do
residuals look like?
A: gain axis and mean-rate axis do not align, but ori tuning is distributed as conic manifold around
gain axis. Using stim-specific gain coupling calculated cross-validated for each neuron, allows pop
response to be predicted with R^2 of 0.72

%see articles:
https://elifesciences.org/articles/8998
https://www.nature.com/articles/nn.3711
https://www.jneurosci.org/content/39/37/7344.abstract
etc

%}
%% define parameters
clear all;
boolSaveFigs = true;
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

%% find data
strStim = 'DG';%DG/NM
dblBimoThreshold = inf;%0.5;%0.4
dblDevThreshold = 0.7;%0.7;%0.7

cellTypes = {'Real','TShuff','TPoiss','xPerm'};
intTypeNum = numel(cellTypes);
sDir = dir([strTargetDataPath 'QC2*.mat']); %or ABA if old
intFileNum = numel(sDir);

%% aggregate data
%[3 x 1] cell (type (real/poiss/shuff), content=vector with r2 per recording
cellR2_Tune = cell(1,intTypeNum); %only tuning
cellR2_Mean = cell(1,intTypeNum);%tuning * pop-mean
cellR2_Gain1 = cell(1,intTypeNum);%tuning * unitary pop-gain
cellR2_Gain = cell(1,intTypeNum);%tuning * stim-specific pop-gain

%[3 x 1] cell (type (real/poiss/shuff), content=vector with average r2 per neuron, mean over all cells
cellR2_TuneN = cell(1,intTypeNum); %only tuning
cellR2_MeanN = cell(1,intTypeNum);%tuning * pop-mean
cellR2_Gain1N = cell(1,intTypeNum);%tuning * unitary pop-gain
cellR2_GainN = cell(1,intTypeNum);%tuning * stim-specific pop-gain
cellFilenames = cell(1,intTypeNum);%names

%projection data;[3 x 1] cell (type (real/poiss/shuff)
cellReflMuDistToOrth = cell(1,intTypeNum);
cellReflMuDistToSelf = cell(1,intTypeNum);
cellReflGainDistToOrth = cell(1,intTypeNum);
cellReflGainDistToSelf = cell(1,intTypeNum);
cellReflRandDistToOrth = cell(1,intTypeNum);
cellReflRandDistToSelf = cell(1,intTypeNum);
cellGainSymmetry = cell(1,intTypeNum);
cellMeanSymmetry = cell(1,intTypeNum);
cellRandSymmetry = cell(1,intTypeNum);
cellAngleMeanAndGain = cell(1,intTypeNum);
cellAngleMeanAndGainRand = cell(1,intTypeNum);
cellSdProjMean = cell(1,intTypeNum);
cellSdProjGain = cell(1,intTypeNum);
cellSdProjGain1 = cell(1,intTypeNum);
cellSdProjRand = cell(1,intTypeNum);
cellSdProjOrth = cell(1,intTypeNum);
cellSdProjAdja = cell(1,intTypeNum);
cellDistFromOrigin = cell(1,intTypeNum);
cellDistFromOrthOri = cell(1,intTypeNum);
cellDistFromAdjaOri = cell(1,intTypeNum);
cellMatDistFromOrthOri = cell(1,intTypeNum);
cellMatDistFromOrigin = cell(1,intTypeNum);
cellMatDistFromAdjaOri = cell(1,intTypeNum);

matPupilPredG = nan(1,intTypeNum);
matPupilPredM = nan(1,intTypeNum);

hTic=tic;
vecCounter = zeros(1,intTypeNum);
for intFile=1:intFileNum
	%% load data
	strFolder = sDir(intFile).folder;
	strFile = sDir(intFile).name;
	strType = strrep(strrep(getFlankedBy(strFile,'QC2Data','_','first'),'ABI',''),['_' strStim],'');
	intType = find(ismember(cellTypes,strType));
	if toc(hTic) > 5
		hTic = tic;
		fprintf('Aggregating %s... Now at %d/%d [%s]',strRunArea,intFile,numel(sDir),getTime);
	end
	%get data
	sData = load(fullpath(strFolder,strFile));
	vecCounter(intType) = vecCounter(intType) + 1;
	
	%save pop pred r2s
	sPrediction = sData.sPrediction;
	%dblBC = sData.dblBC;
	%dblMaxDevFrac = sData.dblMaxDevFrac;
	%if size(sPrediction.matMeanRate,2) < 600,continue;end
	%if size(sData.cellLRActPerQ,2) < 8 || any(flat(cellfun(@(x) any(isnan(x(:))),sData.cellLRActPerQ)))
	%	continue;
	%end
	
	%if sData.dblBC > dblBimoThreshold || sData.dblMaxDevFrac > dblDevThreshold
	%	continue;
	%end
	
	cellR2_Tune{intType}(end+1) = sPrediction.dblR2_Tune; %only tuning
	cellR2_Mean{intType}(end+1) = sPrediction.dblR2_Mean;%tuning * pop-mean
	cellR2_Gain{intType}(end+1) = sPrediction.dblR2_Gain;%tuning * pop-gain
	cellR2_Gain1{intType}(end+1) = sPrediction.dblR2_Gain1;%tuning * pop-gain
	
	%save neuron pred r2s
	cellR2_TuneN{intType}(end+1) = sPrediction.dblPredPerNeuronTune; %only tuning
	cellR2_MeanN{intType}(end+1) = sPrediction.dblPredPerNeuronMean;%tuning * pop-mean
	cellR2_GainN{intType}(end+1) = sPrediction.dblPredPerNeuronGain;%tuning * pop-gain
	cellR2_Gain1N{intType}(end+1) = sPrediction.dblPredPerNeuronGain1;%tuning * pop-gain
	%save name
	cellFilenames{intType}{end+1} = sData.strRec;
	
	%save projection data
	sProjection = sData.sProjection;
	cellFields = fieldnames(sProjection);
	for intField=1:numel(cellFields)
		strField = cellFields{intField};
		strVarName = strcat('cell',strField(4:end));
		eval([strVarName '{intType}(end+1) = mean(sProjection.' strField ');']);
	end
	cellMatDistFromOrthOri{intType}(end+1,:) = sProjection.vecDistFromOrthOri;
	cellMatDistFromOrigin{intType}(end+1,:) = sProjection.vecDistFromOrigin;
	cellMatDistFromAdjaOri{intType}(end+1,:) = sProjection.vecDistFromAdjaOri;
	
	%pupil
	matPupilPredG(vecCounter(intType),intType) = sData.dblRG;
	matPupilPredM(vecCounter(intType),intType) = sData.dblRM;
end
intRecNum = mean(vecCounter);
%% make plot 1
%real pop mean+sd
figure;maxfig
%real
matCol = lines(intTypeNum);
vecOffset = [0.1 -0.1 -0.1 0.1]*0.5;
subplot(2,intTypeNum+1,intTypeNum+1);
hold on;
subplot(2,intTypeNum+1,intTypeNum*2+2);
hold on;
vecX = 1:intTypeNum;
for intType=1:intTypeNum
	strType = cellTypes{intType};
	
	%pop predictions
	subplot(2,intTypeNum+1,intType)
	[h,pMeanGain]=ttest(cellR2_Mean{intType},cellR2_Gain{intType});
	cellLabels = {'Tune','Tune*Mu','Tune*G1','Tune*G'};
	vecMeans = [mean(cellR2_Tune{intType})...
		mean(cellR2_Mean{intType})...
		mean(cellR2_Gain1{intType})...
		mean(cellR2_Gain{intType})...
		];
	vecSds = [std(cellR2_Tune{intType})...
		std(cellR2_Mean{intType})...
		std(cellR2_Gain1{intType})...
		std(cellR2_Gain{intType})...
		];
	intN = numel(cellR2_Tune{intType});
	vecSems = vecSds./sqrt(intN);
	errorbar(vecX,vecMeans,vecSems,'x','color',matCol(intType,:));
	set(gca,'xtick',vecX,'xticklabel',cellLabels);
	ylabel('Cross-val prediction (R^2)');
	title(sprintf('%s; pop, Mu/G,p=%.1e',strType,pMeanGain));
	xlim([vecX(1)-0.5 vecX(end)+0.5]);
	fixfig;
	
	subplot(2,intTypeNum+1,intTypeNum+1);
	errorbar(vecX+vecOffset(intType),vecMeans,vecSems,'x','color',matCol(intType,:));
	
	%neuron predictions
	subplot(2,intTypeNum+1,intTypeNum+1+intType)
	[h,pMeanGain]=ttest(cellR2_MeanN{intType},cellR2_GainN{intType});
	cellLabels = {'Tune','Tune*Mu','Tune*G1','Tune*G'};
	vecMeans = [mean(cellR2_TuneN{intType})...
		mean(cellR2_MeanN{intType})...
		mean(cellR2_Gain1N{intType})...
		mean(cellR2_GainN{intType})...
		];
	vecSds = [std(cellR2_TuneN{intType})...
		std(cellR2_MeanN{intType})...
		std(cellR2_Gain1N{intType})...
		std(cellR2_GainN{intType})...
		];
	intN = numel(cellR2_TuneN{intType});
	vecSems = vecSds./sqrt(intN);
	errorbar(vecX,vecMeans,vecSems,'x','color',matCol(intType,:));
	set(gca,'xtick',vecX,'xticklabel',cellLabels);
	ylabel('Cross-val prediction (R^2)');
	title(sprintf('%s; neuron, Mu/G,p=%.1e',strType,pMeanGain));
	xlim([vecX(1)-0.5 vecX(end)+0.5]);
	fixfig;
	
	subplot(2,intTypeNum+1,intTypeNum*2+2);
	errorbar(vecX+vecOffset(intType),vecMeans,vecSems,'x','color',matCol(intType,:));
	
end

%finish multi-graphs
subplot(2,intTypeNum+1,intTypeNum+1);
hold off;
set(gca,'xtick',vecX,'xticklabel',cellLabels);
ylabel('Cross-val prediction (R^2)');
title(sprintf('Comparison; pop R2, n=%d',intN));
xlim([vecX(1)-0.5 vecX(end)+0.5]);
fixfig;

subplot(2,intTypeNum+1,2*intTypeNum+2);
hold off;
set(gca,'xtick',vecX,'xticklabel',cellLabels);
ylabel('Cross-val prediction (R^2)');
title(sprintf('Comparison; neuron R2, n=%d',intN));
xlim([vecX(1)-0.5 vecX(end)+0.5]);
fixfig;

if boolSaveFigs
	%% save fig
	export_fig(fullpath(strFigurePath,sprintf('QC2Meta_Agg1%s_ActivityPrediction.tif',strStim)));
	export_fig(fullpath(strFigurePath,sprintf('QC2Meta_Agg1%s_ActivityPrediction.pdf',strStim)));
end

%% plot projections
matCol = [0 0 0;...mean
	0 0 0.8;...gain
	0 0.8 0;...orth
	0.8 0 0;...self
	0.5 0 0.5;...(rand)
	];
for intType=1:intTypeNum
	strType = cellTypes{intType};
	%get data
	vecReflMuDistToOrth = cellReflMuDistToOrth{intType};
	vecReflMuDistToSelf = cellReflMuDistToSelf{intType};
	vecReflGainDistToOrth = cellReflGainDistToOrth{intType};
	vecReflGainDistToSelf = cellReflGainDistToSelf{intType};
	vecReflRandDistToOrth = cellReflRandDistToOrth{intType};
	vecReflRandDistToSelf = cellReflRandDistToSelf{intType};
	
	vecGainSymmetry = cellGainSymmetry{intType};
	vecMeanSymmetry = cellMeanSymmetry{intType};
	vecRandSymmetry = cellRandSymmetry{intType};
	
	vecAngleMeanAndGain = cellAngleMeanAndGain{intType};
	vecAngleMeanAndGainRand = cellAngleMeanAndGainRand{intType};
	vecSdProjMean = cellSdProjMean{intType};
	vecSdProjGain1 = cellSdProjGain1{intType};
	vecSdProjGain = cellSdProjGain{intType};
	%vecSdProjRand = cellSdProjRand{intType};
	vecSdProjOrth = cellSdProjOrth{intType};
	vecSdProjAdja = cellSdProjAdja{intType};
	vecDistFromOrigin = cellDistFromOrigin{intType};
	vecDistFromOrthOri = cellDistFromOrthOri{intType};
	vecDistFromAdjaOri = cellDistFromAdjaOri{intType};
	intNumRecs = numel(vecDistFromAdjaOri);
	
	% plot
	[h,pMeanGain]=ttest(vecSdProjMean,vecSdProjGain);
	[h,pMeanOrth]=ttest(vecSdProjMean,vecSdProjOrth);
	[h,pGainOrth]=ttest(vecSdProjGain,vecSdProjOrth);
	figure;maxfig;
	subplot(2,3,1)
	hold on
	vecLoc = [0.2 0.5 0.8 1.1 1.4];
	dblJit = 0.05;
	scatter(vecLoc(1)+(rand(1,intNumRecs)-0.5)*2*dblJit,vecSdProjMean,[200],[0.5 0.5 0.5],'.')
	scatter(vecLoc(2)+(rand(1,intNumRecs)-0.5)*2*dblJit,vecSdProjGain1,[200],[0.5 0.5 0.5],'.')
	scatter(vecLoc(3)+(rand(1,intNumRecs)-0.5)*2*dblJit,vecSdProjGain,[200],[0.5 0.5 0.5],'.')
	scatter(vecLoc(4)+(rand(1,intNumRecs)-0.5)*2*dblJit,vecSdProjOrth,[200],[0.5 0.5 0.5],'.')
	scatter(vecLoc(5)+(rand(1,intNumRecs)-0.5)*2*dblJit,vecSdProjAdja,[200],[0.5 0.5 0.5],'.')
	%scatter(vecLoc(5)+(rand(1,intStimNr)-0.5)*2*dblJit,vecSdProjRand,[200],[0.5 0.5 0.5],'.')
	errorbar(vecLoc(1),mean(vecSdProjMean),std(vecSdProjMean)/sqrt(intNumRecs),'x','color',matCol(1,:),'CapSize',20)
	errorbar(vecLoc(2),mean(vecSdProjGain1),std(vecSdProjGain1)/sqrt(intNumRecs),'x','color',matCol(2,:),'CapSize',20)
	errorbar(vecLoc(3),mean(vecSdProjGain),std(vecSdProjGain)/sqrt(intNumRecs),'x','color',matCol(2,:),'CapSize',20)
	errorbar(vecLoc(4),mean(vecSdProjOrth),std(vecSdProjOrth)/sqrt(intNumRecs),'x','color',matCol(3,:),'CapSize',20)
	errorbar(vecLoc(5),mean(vecSdProjAdja),std(vecSdProjAdja)/sqrt(intNumRecs),'x','color',matCol(4,:),'CapSize',20)
	%errorbar(vecLoc(5),mean(vecSdProjRand),std(vecSdProjRand)/sqrt(intStimNr),'x','color',matCol(1,:),'CapSize',20)
	hold off;
	set(gca,'xtick',vecLoc,'xticklabel',{'Pop mean','Pop gain1','Pop gain','S-Orth','S-Adja'});
	xlabel('Projection axis');
	ylabel('Range of spiking rates (\sigmaHz)');
	title(sprintf('%s; MG,p=%.2e; MO,p=%.3f; GO,p=%.2e',strType,pMeanGain,pMeanOrth,pGainOrth));
	fixfig;
	ylim(gca,[0 max(get(gca,'ylim'))]);
	xlim([vecLoc(1)-0.2 vecLoc(end)+0.2]);
	
	%what is distance from origin to mean R compared with distance between ori Rs?
	subplot(2,3,2);
	hold on
	vecLoc = [0.2 0.5 0.8];
	dblJit = 0.05;
	scatter(vecLoc(1)+(rand(1,intNumRecs)-0.5)*2*dblJit,vecDistFromOrigin,[200],[0.5 0.5 0.5],'.')
	errorbar(vecLoc(1),mean(vecDistFromOrigin),std(vecDistFromOrigin)/sqrt(intNumRecs),'x','color',matCol(2,:),'CapSize',20)
	scatter(vecLoc(2)+(rand(1,intNumRecs)-0.5)*2*dblJit,vecDistFromOrthOri,[200],[0.5 0.5 0.5],'.')
	errorbar(vecLoc(2),mean(vecDistFromOrthOri),std(vecDistFromOrthOri)/sqrt(intNumRecs),'x','color',matCol(3,:),'CapSize',20)
	scatter(vecLoc(3)+(rand(1,intNumRecs)-0.5)*2*dblJit,vecDistFromAdjaOri,[200],[0.5 0.5 0.5],'.')
	errorbar(vecLoc(3),mean(vecDistFromAdjaOri),std(vecDistFromAdjaOri)/sqrt(intNumRecs),'x','color',matCol(4,:),'CapSize',20)
	hold off
	set(gca,'xtick',vecLoc,'xticklabel',{'Origin','Orth. ori','Adja. ori'});
	%ylim([0 60]);
	xlim([vecLoc(1)-0.2 vecLoc(end)+0.2]);
	xlabel('Distance of pop. stim resp to');
	ylabel('Pop. spiking rate distance (\DeltaHz)');
	title(sprintf('Mu +/- sem, n=%d recs',intNumRecs));
	fixfig;
	
	%cos sim of gain with mean vs rand with mean
	%note that random vectors do not lead to cossim of 0 because the space is limited to
	%positive values only (Hz > 0)
	[h,pCosSim] = ttest2(vecAngleMeanAndGain,vecAngleMeanAndGainRand);
	subplot(2,3,3);
	hold on
	vecLoc = [0.2 0.8];
	dblJit = 0.05;
	scatter(vecLoc(1)+(rand(1,intNumRecs)-0.5)*2*dblJit,vecAngleMeanAndGain,[200],[0.5 0.5 0.5],'.')
	errorbar(vecLoc(1),mean(vecAngleMeanAndGain),std(vecAngleMeanAndGain)/sqrt(intNumRecs),'x','color',matCol(2,:),'CapSize',20)
	scatter(vecLoc(2)+(rand(1,intNumRecs)-0.5)*2*dblJit,vecAngleMeanAndGainRand,[200],[0.5 0.5 0.5],'.')
	errorbar(vecLoc(2),mean(vecAngleMeanAndGainRand),std(vecAngleMeanAndGainRand)/sqrt(intNumRecs),'x','color',matCol(1,:),'CapSize',20)
	hold off
	set(gca,'xtick',vecLoc,'xticklabel',{'Gain','Rand'});
	xlim([vecLoc(1)-0.2 vecLoc(end)+0.2]);
	xlabel('Neural axis');
	ylabel('Alignment with mean-axis (cos sim)');
	title(sprintf('T-test, cos-sim mean with gain vs rand: p=%.2e',pCosSim));
	fixfig;
	
	%symmetry around mean or gain axis
	vecLoc = [0.2 0.8 2.2 2.8 4.2 4.8];
	dblJit = 0.05;
	subplot(2,3,4);
	hold on
	errorbar(vecLoc(1),mean(vecReflRandDistToOrth),std(vecReflRandDistToOrth)/sqrt(intNumRecs),'x','color',matCol(3,:),'CapSize',20)
	errorbar(vecLoc(2),mean(vecReflRandDistToSelf),std(vecReflRandDistToSelf)/sqrt(intNumRecs),'x','color',matCol(4,:),'CapSize',20)
	errorbar(vecLoc(3),mean(vecReflMuDistToOrth),std(vecReflMuDistToOrth)/sqrt(intNumRecs),'x','color',matCol(3,:),'CapSize',20)
	errorbar(vecLoc(4),mean(vecReflMuDistToSelf),std(vecReflMuDistToSelf)/sqrt(intNumRecs),'x','color',matCol(4,:),'CapSize',20)
	errorbar(vecLoc(5),mean(vecReflGainDistToOrth),std(vecReflGainDistToOrth)/sqrt(intNumRecs),'x','color',matCol(3,:),'CapSize',20)
	errorbar(vecLoc(6),mean(vecReflGainDistToSelf),std(vecReflGainDistToSelf)/sqrt(intNumRecs),'x','color',matCol(4,:),'CapSize',20)
	hold off
	set(gca,'xtick',vecLoc,'xticklabel',{'Rand/Orth','Self','Mu/Orth','Self','Gain/Orth','Self'});
	xlim([vecLoc(1)-0.2 vecLoc(end)+0.2]);
	xlabel('Reflection over/stim ori');
	ylabel('Prediction error after reflection (\DeltaHz)');
	fixfig;
	
	% symmetry around gain axis
	vecLoc = [1:3];
	subplot(2,3,5);
	[h,pSym]=ttest2(vecMeanSymmetry,vecGainSymmetry);
	hold on
	errorbar(1,mean(vecRandSymmetry),std(vecRandSymmetry)/sqrt(intNumRecs),'d','color',matCol(5,:),'CapSize',20)
	errorbar(2,mean(vecMeanSymmetry),std(vecMeanSymmetry)/sqrt(intNumRecs),'o','color',matCol(1,:),'CapSize',20)
	errorbar(3,mean(vecGainSymmetry),std(vecGainSymmetry)/sqrt(intNumRecs),'o','color',matCol(2,:),'CapSize',20)
	hold off
	ylabel('Manifold symmetry (%)');
	set(gca,'xtick',vecLoc,'xticklabel',{'Random axis','Pop. mean','Pop. gain'});
	xlabel('Axis of symmetry');
	xlim([vecLoc(1)-0.2 vecLoc(end)+0.2]);
	title(sprintf('T-test, p=%.2e',pSym));
	fixfig;
	
	%make example figures of the above
	vecS = 10:12;%randperm(intStimNr,2);
	intS1 = 7;
	intS1Adja = 8;
	intS2 = 12;
	
	vecAvgDistFromOrthOri = mean(cellMatDistFromOrthOri{intType},1);
	vecAvgDistFromOrigin = mean(cellMatDistFromOrigin{intType},1);
	vecAvgDistFromAdjaOri = mean(cellMatDistFromAdjaOri{intType},1);
	
	dblDistToGainCenter = mean(vecAvgDistFromOrigin);
	vecGainAx = [0.1 1];
	dblAlignment=cossim(vecGainAx',ones(size(vecGainAx))');
	dblAng = acos(dblAlignment)/2+deg2rad(45);
	vecGainCenter = [cos(dblAng) sin(dblAng)]*dblDistToGainCenter;
	
	dblAng1 = dblAng + atan((mean(vecAvgDistFromOrthOri(intS1)))/dblDistToGainCenter);
	dblAng1A = dblAng + atan((mean(vecAvgDistFromOrthOri(intS1)) - mean(vecAvgDistFromAdjaOri(intS1))/4)/dblDistToGainCenter);
	dblAng2 = dblAng - atan((mean(vecAvgDistFromOrthOri(intS2)))/dblDistToGainCenter);
	
	vecCenter1 = [cos(dblAng1) sin(dblAng1)]*vecAvgDistFromOrigin(intS1);
	vecCenter1A = [cos(dblAng1A) sin(dblAng1A)]*vecAvgDistFromOrigin(intS1Adja);
	vecCenter2 = [cos(dblAng2) sin(dblAng2)]*vecAvgDistFromOrigin(intS2);
	
	%reflect
	vecX = vecCenter1';
	[dummy,vecXprimeG]=getProjOnLine(vecX,vecGainCenter');
	vecXreflG = 2*vecXprimeG - vecX;
	[dummy,vecXprime]=getProjOnLine(vecX,ones(size(vecX)));
	vecXrefl = 2*vecXprime - vecX;
	[dummy,vecXprimeR]=getProjOnLine(vecX,rand(size(vecX)));
	vecXreflR = 2*vecXprimeR - vecX;
	
	subplot(2,3,6)%[3 4 7 8]);
	cla;
	hold on
	%mean
	dblMaxLim = ceil((max([vecCenter1 vecCenter2]) + max(vecSdProjGain))/5)*5;
	intOffset = 1;
	plot([0 dblMaxLim],[0 dblMaxLim],'-','color',matCol(1,:))
	text(dblMaxLim*0.2+intOffset,dblMaxLim*0.2-intOffset,'Pop. mean axis','FontSize',16,'color',matCol(1,:),'Rotation',45)
	%gain
	plot([0 vecGainCenter(1)*1.4],[0 vecGainCenter(2)*1.4],'-','color',matCol(2,:))
	text(vecGainCenter(1)*0.3+intOffset,vecGainCenter(2)*0.3-intOffset,'Pop. gain axis','FontSize',16,'color',matCol(2,:),'Rotation',rad2deg(dblAng))
	
	%centers
	scatter(vecCenter1(1),vecCenter1(2),[],matCol(4,:),'x');
	scatter(vecCenter2(1),vecCenter2(2),[],matCol(3,:),'x');
	scatter(vecCenter1A(1),vecCenter1A(2),[],[0.5 0 0.5],'x');
	
	
	%covars
	ellipse(vecCenter1(1),vecCenter1(2),mean(vecSdProjGain),mean(vecSdProjOrth),dblAng1,'color',matCol(4,:),'LineStyle','-');
	ellipse(vecCenter2(1),vecCenter2(2),mean(vecSdProjGain),mean(vecSdProjOrth),dblAng2,'color',matCol(3,:),'LineStyle','-');
	
	%draw reflections
	plot([vecCenter1(1) vecXreflG(1)],[vecCenter1(2) vecXreflG(2)],'k--');
	plot([vecCenter1(1) vecXrefl(1)],[vecCenter1(2) vecXrefl(2)],'--','color',[0.5 0.5 0.5]);
	scatter(vecXreflG(1),vecXreflG(2),[],matCol(2,:),'o');
	scatter(vecXrefl(1),vecXrefl(2),[],matCol(1,:),'d');
	
	%text
	text(vecCenter1(1)+intOffset,vecCenter1(2)+intOffset,'Real','FontSize',16,'color',matCol(4,:))
	text(vecCenter2(1)+intOffset,vecCenter2(2)+intOffset,'Real Orth','FontSize',16,'color',matCol(3,:))
	text(vecXreflG(1)+intOffset,vecXreflG(2)-intOffset,'Gain-reflect','FontSize',16,'color',matCol(2,:))
	text(vecXrefl(1)+intOffset,vecXrefl(2)+intOffset,'Mean-reflect','FontSize',16,'color',matCol(1,:))
	text(vecXreflR(1)+intOffset,vecXreflR(2)+intOffset,'Rand-reflect','FontSize',16,'color',matCol(5,:))
	
	%finish fig
	hold off
	axis equal;
	%xlim([0 dblMaxLim*(2/3)]);
	%ylim([0 dblMaxLim]);
	xlabel('Spiking rate axis X (Hz)')
	ylabel('Spiking rate axis Y (Hz)')
	fixfig;grid off;
	
	if boolSaveFigs
		%% save fig
		export_fig(fullpath(strFigurePath,sprintf('QC2Meta_Agg2%s_Projections_%s.tif',strStim,strType)));
		export_fig(fullpath(strFigurePath,sprintf('QC2Meta_Agg2%s_Projections_%s.pdf',strStim,strType)));
	end
end

%% var over gain vs adja
vecColTune = [0.8 0 0];
vecColMean = [0 0 0];
vecColGain1 = [0 0 0.8];
vecColGain = lines(1);

figure;maxfig;
subplot(2,3,1);
vecRatioR2_Tune = (cellR2_Tune{1}./cellR2_Tune{2}-1)*100;
vecRatioR2_Mean = (cellR2_Mean{1}./cellR2_Mean{2}-1)*100;
vecRatioR2_Gain1 = (cellR2_Gain1{1}./cellR2_Gain1{2}-1)*100;
vecRatioR2_Gain = (cellR2_Gain{1}./cellR2_Gain{2}-1)*100;

vecRatioR2_Tune = cellR2_Tune{1}-cellR2_Tune{2};
vecRatioR2_Mean = cellR2_Mean{1}-cellR2_Mean{2};
vecRatioR2_Gain1 = cellR2_Gain1{1}-cellR2_Gain1{2};
vecRatioR2_Gain = cellR2_Gain{1}-cellR2_Gain{2};

[h,pMeanGain] = ttest(vecRatioR2_Mean,vecRatioR2_Gain);
hold on
plot(cat(1,vecRatioR2_Tune,vecRatioR2_Mean,vecRatioR2_Gain1,vecRatioR2_Gain),'color',[0.7 0.7 0.7]);
errorbar(1,mean(vecRatioR2_Tune),std(vecRatioR2_Tune)./sqrt(numel(vecRatioR2_Tune)),'x','color',vecColTune);
errorbar(2,mean(vecRatioR2_Mean),std(vecRatioR2_Mean)./sqrt(numel(vecRatioR2_Mean)),'x','color',vecColMean);
errorbar(3,mean(vecRatioR2_Gain1),std(vecRatioR2_Gain1)./sqrt(numel(vecRatioR2_Gain1)),'x','color',vecColGain1);
errorbar(4,mean(vecRatioR2_Gain),std(vecRatioR2_Gain)./sqrt(numel(vecRatioR2_Gain)),'x','color',vecColGain);
hold off
ylabel('Pop resp R^2 improvement over shuffled');
set(gca,'xtick',1:4,'xticklabel',{'Tune','Tune*Mu','Tune*G1','Tune*G'});
xlim([0.5 4.5]);
title(sprintf('T-test mean-gain, p=%.2e',pMeanGain));
fixfig;grid off;


subplot(2,3,2)
hold on
for intType=1:intTypeNum
	matGains(intType,:) = cellSdProjAdja{intType}./(cellSdProjGain{intType}+cellSdProjAdja{intType});
end
plot(matGains,'color',[0.7 0.7 0.7]);
for intType=1:intTypeNum
	errorbar(intType,mean(matGains(intType,:)),std(matGains(intType,:))./sqrt(numel(matGains(intType,:))),'x','color',matCol(intType,:));
end
[h,pRealShuff] = ttest(matGains(1,:),matGains(2,:));


hold off
ylabel('Fraction of info-limiting variability');
set(gca,'xtick',1:intTypeNum,'xticklabel',cellTypes);
xlim([0.5 intTypeNum+0.5]);
title(sprintf('T-test real-shuff, p=%.2e',pRealShuff));
fixfig;grid off;

subplot(2,3,3)
vecGain1Real = cellSdProjGain1{1};
vecGain1Shuff = cellSdProjGain1{2};

vecGainReal = cellSdProjGain{1};
vecGainShuff = cellSdProjGain{2};

vecGain1RealShuffDiff = (vecGain1Real - vecGain1Shuff);
vecGainRealShuffDiff = (vecGainReal - vecGainShuff);


[h,pRealShuff2] = ttest(vecGain1RealShuffDiff,vecGainRealShuffDiff);
hold on
plot(cat(1,vecGain1RealShuffDiff,vecGainRealShuffDiff),'color',[0.7 0.7 0.7]);
errorbar(1,mean(vecGain1RealShuffDiff),std(vecGain1RealShuffDiff)./sqrt(numel(vecGain1RealShuffDiff)),'x','color',vecColGain1);
errorbar(2,mean(vecGainRealShuffDiff),std(vecGainRealShuffDiff)./sqrt(numel(vecGainRealShuffDiff)),'x','color',vecColGain);
hold off
ylabel('Real-Shuff diff sd');
set(gca,'xtick',1:2,'xticklabel',{'Gain1','Gain'});
xlim([0.5 2.5]);
title(sprintf('T-test real-shuff, p=%.2e',pRealShuff2));
fixfig;grid off;

subplot(2,3,4)
vecRatioR2_Mean_MeanNorm = 100*(vecRatioR2_Mean - vecRatioR2_Mean)./vecRatioR2_Mean;
vecRatioR2_Gain_MeanNorm = 100*(vecRatioR2_Gain - vecRatioR2_Mean)./vecRatioR2_Mean;

[h,pMeanGain] = ttest(vecRatioR2_Mean_MeanNorm,vecRatioR2_Gain_MeanNorm);
hold on
plot(cat(1,vecRatioR2_Mean_MeanNorm,vecRatioR2_Gain_MeanNorm),'color',[0.7 0.7 0.7]);
errorbar(1,mean(vecRatioR2_Mean_MeanNorm),std(vecRatioR2_Mean_MeanNorm)./sqrt(numel(vecRatioR2_Mean_MeanNorm)),'x','color',vecColMean);
errorbar(2,mean(vecRatioR2_Gain_MeanNorm),std(vecRatioR2_Gain_MeanNorm)./sqrt(numel(vecRatioR2_Gain_MeanNorm)),'x','color',vecColGain);
hold off
ylabel('Pred. improv. norm. to Tune*Mu (%)');
set(gca,'xtick',1:2,'xticklabel',{'Tune*Mu','Tune*G'});
xlim([0.5 2.5]);
title(sprintf('T-test mean-gain, p=%.2e',pMeanGain));
fixfig;grid off;

subplot(2,3,5)
vecRatioR2_Gain1_Gain1Norm = 100*(vecRatioR2_Gain1 - vecRatioR2_Gain1)./vecRatioR2_Gain1;
vecRatioR2_Gain_Gain1Norm = 100*(vecRatioR2_Gain - vecRatioR2_Gain1)./vecRatioR2_Gain1;

[h,pGain1Gain] = ttest(vecRatioR2_Gain1_Gain1Norm,vecRatioR2_Gain_Gain1Norm);
hold on
plot(cat(1,vecRatioR2_Gain1_Gain1Norm,vecRatioR2_Gain_Gain1Norm),'color',[0.7 0.7 0.7]);
errorbar(1,mean(vecRatioR2_Gain1_Gain1Norm),std(vecRatioR2_Gain1_Gain1Norm)./sqrt(numel(vecRatioR2_Gain1_Gain1Norm)),'x','color',vecColGain1);
errorbar(2,mean(vecRatioR2_Gain_Gain1Norm),std(vecRatioR2_Gain_Gain1Norm)./sqrt(numel(vecRatioR2_Gain_Gain1Norm)),'x','color',vecColGain);
hold off
ylabel('Pred. improv. norm. to Gain1 (%)');
set(gca,'xtick',1:2,'xticklabel',{'Tune*G1','Tune*G'});
xlim([0.5 2.5]);
title(sprintf('T-test gain1-gain, p=%.2e',pGain1Gain));
fixfig;grid off;

subplot(2,3,6)
vecAvgGainSymmetry = cellfun(@mean,cellGainSymmetry);
vecSdGainSymmetry = cellfun(@std,cellGainSymmetry);
vecAvgMeanSymmetry = cellfun(@mean,cellMeanSymmetry);
vecSdMeanSymmetry = cellfun(@std,cellMeanSymmetry);
vecAvgRandSymmetry = cellfun(@mean,cellRandSymmetry);
vecSdRandSymmetry = cellfun(@std,cellRandSymmetry);
errorbar(1:intTypeNum,vecAvgGainSymmetry,vecSdGainSymmetry./sqrt(intRecNum));hold on
errorbar(1:intTypeNum,vecAvgMeanSymmetry,vecSdMeanSymmetry./sqrt(intRecNum),'color','k');
errorbar(1:intTypeNum,vecAvgRandSymmetry,vecSdRandSymmetry./sqrt(intRecNum),'color','r');
ylim([0 max(get(gca,'ylim'))]);
set(gca,'xtick',1:intTypeNum,'xticklabel',cellTypes);
xlim([0.5 intTypeNum+0.5]);
ylabel('Manifold symmetry across axis (%)');
legend(gca,{'Gain','Mean','Rand'});
fixfig;



if boolSaveFigs
	%% save fig
	export_fig(fullpath(strFigurePath,sprintf('QC2Meta_Agg3%s_PredAndInfoLimVar.tif',strStim)));
	export_fig(fullpath(strFigurePath,sprintf('QC2Meta_Agg3%s_PredAndInfoLimVar.pdf',strStim)));
end

%%
figure;hold on
matPupilPredGRem = matPupilPredG;
matPupilPredMRem = matPupilPredM;
indRem = any(matPupilPredGRem < -0.3 | isnan(matPupilPredGRem) | isnan(matPupilPredMRem),2);
matPupilPredGRem(indRem,:) = [];
matPupilPredMRem(indRem,:) = [];
vecMuG = mean(matPupilPredGRem,1);
vecSdG = std(matPupilPredGRem,[],1);

vecMuM = mean(matPupilPredMRem,1);
vecSdM = std(matPupilPredMRem,[],1);
intGoodRecs = sum(~indRem);

matCol=lines(intTypeNum);
vecH = [];
for intType=1:intTypeNum
	vecH(intType) = errorbar(0.75+0.1*intType,vecMuG(intType),vecSdG(intType)./sqrt(intGoodRecs),'x','color',matCol(intType,:));
	errorbar(1.75+0.1*intType,vecMuM(intType),vecSdM(intType)./sqrt(intGoodRecs),'x','color',matCol(intType,:));
end
set(gca,'xtick',[1 2],'xticklabel',{'Gain','Mean'});
xlim([0.5 2.5])
xlabel('Pop activity metric');
ylabel('Pearson corr with pupil size (r)');
legend(vecH,cellTypes,'location','best');
fixfig;

if boolSaveFigs
	%% save fig
	export_fig(fullpath(strFigurePath,sprintf('QC2Meta_Pupil%s.tif',strStim)));
	export_fig(fullpath(strFigurePath,sprintf('QC2Meta_Pupil%s.pdf',strStim)));
end
