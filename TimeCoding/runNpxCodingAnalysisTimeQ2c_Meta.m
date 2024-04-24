
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
cellTypes =  {'Real','Poiss'};
intPlotNum = 2*numel(cellTypes);
vecH = nan([1 intPlotNum]);
vecH1_5 = nan([1 intPlotNum]);
vecH2 = nan([1 intPlotNum]);
for intType=1:(intPlotNum)
	vecH(intType) = subplot(3,intPlotNum,intType);hold on;
	vecH1_5(intType) = subplot(3,intPlotNum,intType+intPlotNum);hold on;
	vecH2(intType) = subplot(3,intPlotNum,intType+2*intPlotNum);hold on;
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
		
		cellCV{intFile,intType} = sAggData(intUseEntry).matCV;
		
		matMean = sAggData(intUseEntry).matMean;
		matSd = sAggData(intUseEntry).matSd;
		matCV = sAggData(intUseEntry).matCV;
		
		vecMeanNoJitter = matMean(1,:);
		vecMeanFullJitter = matMean(end,:);
		vecSdNoJitter = matSd(1,:);
		vecSdFullJitter = matSd(end,:);
		vecCVNoJitter = matCV(1,:);
		vecCVFullJitter = matCV(end,:);
		
		%singles
		cellCV_Single{intFile,intType} = sAggData(intUseEntry).mat1CV;
		cellRootR2_Single{intFile,intType} = sAggData(intUseEntry).mat1R2_Root;
		cellRootExponent_Single{intFile,intType} = sAggData(intUseEntry).mat1Exponent_Root;
		
		%% plot sd/mean
		intPlot = (intType-1)*2+1;
		plot(vecH(intPlot),vecMeanNoJitter,vecSdNoJitter,'color',matCol(intType,:));
		plot(vecH(intPlot+1),vecMeanFullJitter,vecSdFullJitter,'color',matCol(intType,:));
		
		plot(vecH1_5(intPlot),vecMeanNoJitter,vecSdNoJitter.^2,'color',matCol(intType,:));
		plot(vecH1_5(intPlot+1),vecMeanFullJitter,vecSdFullJitter.^2,'color',matCol(intType,:));
		
		%% plot cv/timescale
		plot(vecH2(intPlot),vecTimescales,vecCVNoJitter,'color',matCol(intType,:));
		plot(vecH2(intPlot+1),vecTimescales,vecCVFullJitter,'color',matCol(intType,:));
		
	end
end

for intType=1:numel(cellTypes)
	intPlot = (intType-1)*2+1;
	%finish plots 1; mean and sd
	xlabel(vecH(intPlot),'Mean of spike counts');
	ylabel(vecH(intPlot),'Sd of spike counts');
	title(vecH(intPlot),sprintf('No jitter, %s',cellTypes{intType}),'interpreter','none');
	%finish plots 1; mean and sd
	xlabel(vecH(intPlot+1),'Mean of spike counts');
	ylabel(vecH(intPlot+1),'Sd of spike counts');
	title(vecH(intPlot+1),sprintf('Full jitter, %s',cellTypes{intType}),'interpreter','none');
	
	%finish plots 1.5; mean and var
	xlabel(vecH1_5(intPlot),'Mean of spike counts');
	ylabel(vecH1_5(intPlot),'Var of spike counts');
	title(vecH1_5(intPlot),sprintf('No jitter, %s',cellTypes{intType}),'interpreter','none');
	%finish plots 1.5; mean and var
	xlabel(vecH1_5(intPlot+1),'Mean of spike counts');
	ylabel(vecH1_5(intPlot+1),'Var of spike counts');
	title(vecH1_5(intPlot+1),sprintf('Full jitter, %s',cellTypes{intType}),'interpreter','none');
	
	%finish plots 2; timescale and cv
	xlabel(vecH2(intPlot),'Timescale (s)');
	ylabel(vecH2(intPlot),'CV (sd/mu)');
	title(vecH2(intPlot),sprintf('No jitter,%s',cellTypes{intType}),'interpreter','none');
	set(vecH2(intPlot),'ylim',[0 1.2]);
	%finish plots 2; timescale and cv
	xlabel(vecH2(intPlot+1),'Timescale (s)');
	ylabel(vecH2(intPlot+1),'CV (sd/mu)');
	title(vecH2(intPlot+1),sprintf('Full jitter,%s',cellTypes{intType}),'interpreter','none');
	set(vecH2(intPlot+1),'ylim',[0 1.2]);
end
fixfig;

drawnow;
export_fig(fullpath(strFigurePath,sprintf('Q2c_CV_Timescales.tif')));
export_fig(fullpath(strFigurePath,sprintf('Q2c_CV_Timescales.pdf')));


%% plot cv/timescale for real and poiss with/without jitter and for full pop and single neurons [2 x 2 x 2]
vecTimescales = sAggData(1).vecTimescales;
vecJitter = sAggData(1).vecJitter;

%whole pops
vecR2_RealNoJitter = cellfun(@(x) x(1), cellRootR2(:,1));
vecR2_PoissNoJitter = cellfun(@(x) x(1), cellRootR2(:,2));
vecR2_RealFullJitter = cellfun(@(x) x(end), cellRootR2(:,1));
vecR2_PoissFullJitter = cellfun(@(x) x(end), cellRootR2(:,2));
clear matR2;
matR2(:,1,1) = vecR2_RealNoJitter;
matR2(:,2,1) = vecR2_PoissNoJitter;
matR2(:,1,2) = vecR2_RealFullJitter;
matR2(:,2,2) = vecR2_PoissFullJitter;

vecExponent_RealNoJitter = cellfun(@(x) x(1), cellRootExponent(:,1));
vecExponent_PoissNoJitter = cellfun(@(x) x(1), cellRootExponent(:,2));
vecExponent_RealFullJitter = cellfun(@(x) x(end), cellRootExponent(:,1));
vecExponent_PoissFullJitter = cellfun(@(x) x(end), cellRootExponent(:,2));
clear matExponent;
matExponent(:,1,1) = vecExponent_RealNoJitter;
matExponent(:,2,1) = vecExponent_PoissNoJitter;
matExponent(:,1,2) = vecExponent_RealFullJitter;
matExponent(:,2,2) = vecExponent_PoissFullJitter;

matCV_RealNoJitter = cell2mat(cellfun(@(x) x(1,:)', cellCV(:,1),'UniformOutput',false)'); %[timescale x neuron]
matCV_PoissNoJitter = cell2mat(cellfun(@(x) x(1,:)', cellCV(:,2),'UniformOutput',false)'); %[timescale x neuron]
matCV_RealFullJitter = cell2mat(cellfun(@(x) x(end,:)', cellCV(:,1),'UniformOutput',false)'); %[timescale x neuron]
matCV_PoissFullJitter = cell2mat(cellfun(@(x) x(end,:)', cellCV(:,2),'UniformOutput',false)'); %[timescale x neuron]
clear matCV
matCV(:,:,1,1) = matCV_RealNoJitter;
matCV(:,:,2,1) = matCV_PoissNoJitter;
matCV(:,:,1,2) = matCV_RealFullJitter;
matCV(:,:,2,2) = matCV_PoissFullJitter;

%single neurons [real/poiss with/without-jitter neuron-id]
vecSingleR2_RealNoJitter = cell2vec(cellfun(@(x) x(1,:),cellRootR2_Single(:,1),'UniformOutput',false));
vecSingleR2_PoissNoJitter = cell2vec(cellfun(@(x) x(1,:),cellRootR2_Single(:,2),'UniformOutput',false));
vecSingleR2_RealFullJitter = cell2vec(cellfun(@(x) x(end,:),cellRootR2_Single(:,1),'UniformOutput',false));
vecSingleR2_PoissFullJitter = cell2vec(cellfun(@(x) x(end,:),cellRootR2_Single(:,2),'UniformOutput',false));
clear matSingleR2;
matSingleR2(:,1,1) = vecSingleR2_RealNoJitter;
matSingleR2(:,2,1) = vecSingleR2_PoissNoJitter;
matSingleR2(:,1,2) = vecSingleR2_RealFullJitter;
matSingleR2(:,2,2) = vecSingleR2_PoissFullJitter;

vecSingleExponent_RealNoJitter = cell2vec(cellfun(@(x) x(1,:),cellRootExponent_Single(:,1),'UniformOutput',false));
vecSingleExponent_PoissNoJitter = cell2vec(cellfun(@(x) x(1,:),cellRootExponent_Single(:,2),'UniformOutput',false));
vecSingleExponent_RealFullJitter = cell2vec(cellfun(@(x) x(end,:),cellRootExponent_Single(:,1),'UniformOutput',false));
vecSingleExponent_PoissFullJitter = cell2vec(cellfun(@(x) x(end,:),cellRootExponent_Single(:,2),'UniformOutput',false));
clear matSingleExponent;
matSingleExponent(:,1,1) = vecSingleExponent_RealNoJitter;
matSingleExponent(:,2,1) = vecSingleExponent_PoissNoJitter;
matSingleExponent(:,1,2) = vecSingleExponent_RealFullJitter;
matSingleExponent(:,2,2) = vecSingleExponent_PoissFullJitter;

%concatenate cvs
matSingleCV_RealNoJitter = cell2mat(cellfun(@(x) squeeze(x(:,1,:)),cellCV_Single(:,1),'UniformOutput',false)'); %[timescale x neuron]
matSingleCV_PoissNoJitter = cell2mat(cellfun(@(x) squeeze(x(:,1,:)),cellCV_Single(:,2),'UniformOutput',false)'); %[timescale x neuron]
matSingleCV_RealFullJitter = cell2mat(cellfun(@(x) squeeze(x(:,end,:)),cellCV_Single(:,1),'UniformOutput',false)'); %[timescale x neuron]
matSingleCV_PoissFullJitter = cell2mat(cellfun(@(x) squeeze(x(:,end,:)),cellCV_Single(:,2),'UniformOutput',false)'); %[timescale x neuron]
clear matSingleCV;
matSingleCV(:,:,1,1) = matSingleCV_RealNoJitter;
matSingleCV(:,:,2,1) = matSingleCV_PoissNoJitter;
matSingleCV(:,:,1,2) = matSingleCV_RealFullJitter;
matSingleCV(:,:,2,2) = matSingleCV_PoissFullJitter;
%matSingleCV: [timescale x neuron x real/poiss x (no-)/jitter]
%matCV: [timescale x pop x real/poiss x (no-)/jitter]

%% plot
%plot cv/timescale for real and poiss with/without jitter and for full pop and single neurons [2 x 2 x 2]
cellPop = {'Pop','Single'};
cellJit = {'NoJit','FullJit'};
hSummaryFig = figure;maxfig;
hPlotCVs(1) = subplot(2,3,1);hold on;
hPlotCVs(2) = subplot(2,3,4);hold on;
matColorPopSingle = [lines(1);[0.8 0 0]];
cellLineType = {'-','-';'--','--'};
cellSubLegend = {};
hAllFig=figure;maxfig;
cellLegend = {};

for intRealPoiss=1:2
	%real/poiss
	for intPopSingle=1:2
		if intPopSingle==1
			%pop
			matPlotCV = matCV;
		else
			%single
			matPlotCV = matSingleCV;
		end
		for intJit=1:2
			%(no-)jitter
			matData = matPlotCV(:,:,intRealPoiss,intJit);
			vecMean = mean(matData,2);
			vecSd = std(matData,[],2);
			
			intPlotNr = (intRealPoiss-1)*4+(intPopSingle-1)*2+intJit;
			subplot(2,4,intPlotNr);
			plot(vecTimescales,vecMean,'color',lines(1));
			hold on
			plot(vecTimescales,vecMean-vecSd,'--','color',lines(1));
			plot(vecTimescales,vecMean+vecSd,'--','color',lines(1));
			strTag = sprintf('%s - %s - %s',cellPop{intPopSingle},cellTypes{intRealPoiss},cellJit{intJit});
			title(strTag);
			cellLegend{intPlotNr} = strTag;
			xlabel('Timescale (s)');
			ylabel('CV (sd/mu)');
			
			%plot mean
			intPlotNr2 = (intPopSingle-1)*2+intJit;
			plot(hPlotCVs(intRealPoiss),vecTimescales,vecMean,cellLineType{intJit,intRealPoiss},'color',matColorPopSingle(intPopSingle,:));
			cellSubLegend{intPlotNr2} = sprintf('%s - %s',cellPop{intPopSingle},cellJit{intJit});
		end
	end
	legend(hPlotCVs(intRealPoiss),cellSubLegend);
	xlabel(hPlotCVs(intRealPoiss),'Timescale (s)');
	ylabel(hPlotCVs(intRealPoiss),'CV (sd/mu)');
	title(hPlotCVs(intRealPoiss),cellTypes{intRealPoiss});
	set(hPlotCVs(intRealPoiss),'xscale','log');
	set(hPlotCVs(intRealPoiss),'yscale','log');
end

%plot exponents
Xmean=0;
figure(hSummaryFig);
vecColTheory = [0.5 0.5 0.5];
for intRealPoiss=1:2
	subplot(2,3,2+(intRealPoiss-1)*3);hold on;
	plot([0.5 4.5],[0.5 0.5],'--','color',vecColTheory);
	text(0.1,0.525,'Theory','fontsize',14,'color',vecColTheory);
	for intPopSingle=1:2
		if intPopSingle==1
			%pop
			matPlotExp = matExponent;
		else
			%single
			matPlotExp = matSingleExponent;
		end
		
		%real/poiss
		for intJit=1:2
			%(no-)jitter
			vecData = matPlotExp(:,intRealPoiss,intJit);
			Xprev=Xmean;
			Xmean=mean(vecData);
			intPlotNr = (intPopSingle-1)*2+intJit;
			errorbar(intPlotNr,Xmean,std(vecData)./sqrt(intRecNum),'x','color',matColorPopSingle(intPopSingle,:));
			if intJit==2
				plot([intPlotNr-1 intPlotNr],[Xprev Xmean],'color',matColorPopSingle(intPopSingle,:));
			end
		end
	end
	set(gca,'xtick',1:4,'xticklabel',cellSubLegend);
	xtickangle(gca,45);
	xlim([0 5]);
	ylim([0 0.6]);
	ylabel('Root fit exponent');
	title([cellTypes{intRealPoiss} ' ;Mu +/- sem of beta fit; a+(1/((b*x)^c))']);
end

%plot R2
for intRealPoiss=1:2
	subplot(2,3,3+(intRealPoiss-1)*3);hold on;
	for intPopSingle=1:2
		if intPopSingle==1
			%pop
			matPlotR2 = matR2;
		else
			%single
			matPlotR2 = matSingleR2;
		end
		%real/poiss
		for intJit=1:2
			%(no-)jitter
			vecData = matPlotR2(:,intRealPoiss,intJit);
			
			dblAlpha = 2-2*normcdf(2);
			[phat,pci] = betafit(vecData,dblAlpha);
			Xprev=Xmean;
			Xupper = betainv(dblAlpha,pci(1,1),pci(1,2));
			Xmean = betainv(dblAlpha,phat(1),phat(2));
			Xlower = betainv(dblAlpha,pci(2,1),pci(2,2));
			
			intPlotNr = (intPopSingle-1)*2+intJit;
			errorbar(intPlotNr,Xmean,Xlower-Xmean,Xupper-Xmean,'x','color',matColorPopSingle(intPopSingle,:));
			if intJit==2
				plot([intPlotNr-1 intPlotNr],[Xprev Xmean],'color',matColorPopSingle(intPopSingle,:));
			end
		end
	end
	set(gca,'xtick',1:4,'xticklabel',cellSubLegend);
	xtickangle(gca,45);
	xlim([0 5]);
	ylim([0 1]);
	ylabel('Root fit R^2');
	title([cellTypes{intRealPoiss} ' ;Mu +/- 95-CI of beta fit; a+(1/((b*x)^c))']);
end


fixfig;