
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
%strTargetDataPath = 'C:\Drive\PopTimeCoding\data\q2c\';
sFiles = dir ([strTargetDataPath strTag '*' strOnset '.mat']);
intRecNum = numel(sFiles);

%% pre-allocate
figure;maxfig;
cellTypes =  {'Real','Poiss','ShuffTid','RandTid','RandTxClass','Uniform'};
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
		
		%linfit
		vecUseEntries = 1:numel(vecTimescales);%round(logspace(0,3,15));
		intK=1;
		fLin = fittype('a*x',...
			'dependent',{'y'},'independent',{'x'},...
			'coefficients',{'a'});
		intJitNum = numel(vecJitter);
		vecLinSlopesMuSd = nan(1,intJitNum);
		vecLinR2MuSd = nan(1,intJitNum);
		for intJit=1:intJitNum
			vecX = matMean(intJit,vecUseEntries);
			vecY = matSd(intJit,vecUseEntries);
			[fitobject,gof] = fit(vecX',vecY',fLin,'lower',0,'upper',1e16,'startpoint',mean(vecY)/mean(vecX));
			dblSlope_SdMean = fitobject.a;
			vecFitY_SdMean = fitobject(vecX);
			[dblR2_Lin,dblSS_tot,dblSS_res,dblT,dblP,dblR2_adjusted,dblR2_SE] = getR2(vecY,vecFitY_SdMean,intK);
			
			vecLinSlopesMuSd(intJit) = dblSlope_SdMean;
			vecLinR2MuSd(intJit) = dblR2_Lin;
		end
		
		cellLinR2{intFile,intType} = vecLinR2MuSd;
		cellLinSlope{intFile,intType} = vecLinSlopesMuSd;
		
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
export_fig(fullpath(strFigurePath,sprintf('Q2c_CV_TimescalesAllPlots.tif')));
export_fig(fullpath(strFigurePath,sprintf('Q2c_CV_TimescalesAllPlots.pdf')));


%% plot cv/timescale for real and poiss with/without jitter and for full pop and single neurons [2 x 2 x 2]
vecTimescales = sAggData(1).vecTimescales;
vecJitter = sAggData(1).vecJitter;
intTypeNum = numel(cellTypes);

%whole pops
%concatenate r2s
clear matR2;
for intType=1:intTypeNum
	matR2(:,intType,1) = cell2mat(cellfun(@(x) squeeze(x(1,:)),cellRootR2(:,intType),'UniformOutput',false)');
	matR2(:,intType,2) = cell2mat(cellfun(@(x) squeeze(x(end,:)),cellRootR2(:,intType),'UniformOutput',false)');
end

%concatenate exps
clear matExponent;
for intType=1:intTypeNum
	matExponent(:,intType,1) = cell2mat(cellfun(@(x) squeeze(x(1,:)),cellRootExponent(:,intType),'UniformOutput',false)');
	matExponent(:,intType,2) = cell2mat(cellfun(@(x) squeeze(x(end,:)),cellRootExponent(:,intType),'UniformOutput',false)');
end

%concatenate cvs
clear matCV;
for intType=1:intTypeNum
	matCV(:,:,intType,1) = cell2mat(cellfun(@(x) x(1,:)',cellCV(:,intType),'UniformOutput',false)');
	matCV(:,:,intType,2) = cell2mat(cellfun(@(x) x(end,:)',cellCV(:,intType),'UniformOutput',false)');
end

%single neurons [real/poiss with/without-jitter neuron-id]
%concatenate r2s
clear matSingleR2;
for intType=1:intTypeNum
	matSingleR2(:,intType,1) = cell2mat(cellfun(@(x) x(1,:),cellRootR2_Single(:,intType),'UniformOutput',false)');
	matSingleR2(:,intType,2) = cell2mat(cellfun(@(x) x(end,:),cellRootR2_Single(:,intType),'UniformOutput',false)');
end

%concatenate exps
clear matSingleExponent;
for intType=1:intTypeNum
	matSingleExponent(:,intType,1) = cell2mat(cellfun(@(x) x(1,:),cellRootExponent_Single(:,intType),'UniformOutput',false)');
	matSingleExponent(:,intType,2) = cell2mat(cellfun(@(x) x(end,:),cellRootExponent_Single(:,intType),'UniformOutput',false)');
end

%concatenate cvs
clear matSingleCV;
for intType=1:intTypeNum
	matSingleCV(:,:,intType,1) = cell2mat(cellfun(@(x) squeeze(x(:,1,:)),cellCV_Single(:,intType),'UniformOutput',false)');
	matSingleCV(:,:,intType,2) = cell2mat(cellfun(@(x) squeeze(x(:,end,:)),cellCV_Single(:,intType),'UniformOutput',false)');
end

%matSingleCV: [timescale x neuron x real/poiss x (no-)/jitter]
%matCV: [timescale x pop x real/poiss x (no-)/jitter]

%% plot
%plot cv/timescale for real and poiss with/without jitter and for full pop and single neurons [2 x 2 x 2]
cellPop = {'Pop','Single'};
cellJit = {'NoJit','FullJit'};
hSummaryFig = figure;maxfig;
hPlotCVS = [];
for intType=1:intTypeNum
	hPlotCVs(intType) = subplot(2,intTypeNum,intType);hold on;
end
matColorPopSingle = lines(2);
cellLineType = repmat({'-';'--'},[1 intTypeNum]);
cellSubLegend = {};
cellLegend = {};
boolPlotAllFigs = false;
if boolPlotAllFigs,hAllFig=figure;maxfig;end
for intType=1:intTypeNum
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
			matData = matPlotCV(:,:,intType,intJit);
			vecMean = mean(matData,2);
			vecSd = std(matData,[],2);
			
			intPlotNr = (intType-1)*4+(intPopSingle-1)*2+intJit;
			strTag = sprintf('%s - %s - %s',cellPop{intPopSingle},cellTypes{intType},cellJit{intJit});
			cellLegend{intPlotNr} = strTag;
			if boolPlotAllFigs
				subplot(intTypeNum,4,intPlotNr);
				plot(vecTimescales,vecMean,'color',lines(1));
				hold on
				plot(vecTimescales,vecMean-vecSd,'--','color',lines(1));
				plot(vecTimescales,vecMean+vecSd,'--','color',lines(1));
				title(strTag);
				xlabel('Timescale (s)');
				ylabel('CV (sd/mu)');
			end
			%plot mean
			intPlotNr2 = (intPopSingle-1)*2+intJit;
			%plot(hPlotCVs(intType),vecTimescales,matData',cellLineType{intJit,intType},'color',matColorPopSingle(intPopSingle,:));
			plot(hPlotCVs(intType),vecTimescales,vecMean,cellLineType{intJit,intType},'color',matColorPopSingle(intPopSingle,:));
			cellSubLegend{intPlotNr2} = sprintf('%s - %s',cellPop{intPopSingle},cellJit{intJit});
		end
	end
	legend(hPlotCVs(intType),cellSubLegend,'location','best');
	xlabel(hPlotCVs(intType),'Timescale (s)');
	ylabel(hPlotCVs(intType),'CV (sd/mu)');
	title(hPlotCVs(intType),cellTypes{intType});
	set(hPlotCVs(intType),'xscale','log');
	set(hPlotCVs(intType),'yscale','log');
end

%plot exponents
Xmean=0;
figure(hSummaryFig);
vecColTheory = [0.5 0.5 0.5];
for intType=1:intTypeNum
	subplot(2,intTypeNum,intTypeNum+intType);hold on;
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
			vecData = matPlotExp(:,intType,intJit);
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
	title([cellTypes{intType} ' ;Mu +/- sem of beta fit; a+(1/((b*x)^c))']);
end
fixfig;

drawnow;
export_fig(fullpath(strFigurePath,sprintf('Q2c_CV_TimescaleFits.tif')));
export_fig(fullpath(strFigurePath,sprintf('Q2c_CV_TimescaleFits.pdf')));

%% plot d primes and lin slopes
%run tests
%matSingleCV: [timescale x neuron x real/poiss x (no-)/jitter]
%matCV: [timescale x pop x real/poiss x (no-)/jitter]

% plot
cellTypeLine = {'-','--',':','-.','--',':'};
figure;maxfig;
subplot(2,3,1);hold on
cellLegendSub3 = {};
for intType=1:intTypeNum
	%get data
	matPopCV_NoJit = matCV(:,:,intType,1);
	matPopCV_Jit = matCV(:,:,intType,2);
	[dummy,vecP]=ttest(log10(matPopCV_NoJit'),log10(matPopCV_Jit'));
	vecP_Pop = vecP*numel(vecP);
	vecDprime_Pop = getdprime2(matPopCV_NoJit,matPopCV_Jit,2);
	
	matSingleCV_NoJit = matSingleCV(:,:,intType,1);
	matSingleCV_Jit = matSingleCV(:,:,intType,2);
	[dummy,vecP]=ttest(log10(matSingleCV_NoJit'),log10(matSingleCV_Jit'));
	vecP_SinglePoiss = vecP*numel(vecP);
	vecDprime_Single = getdprime2(matSingleCV_NoJit,matSingleCV_Jit,2);
	
	%plot
	plot(vecTimescales,vecDprime_Pop,cellTypeLine{intType},'color',matColorPopSingle(1,:));
	plot(vecTimescales,vecDprime_Single,cellTypeLine{intType},'color',matColorPopSingle(2,:));
	cellLegendSub3{intType*2-1} = ['Pop-' cellTypes{intType}];
	cellLegendSub3{intType*2} = ['Single-' cellTypes{intType}];
end

set(gca,'xscale','log');
ylabel('d'' NoJit/FullJit');
legend(cellLegendSub3,'location','northwest');
xlabel('Timescale (s)');

subplot(2,3,2);hold on
Xmean=0;
cellSubLegend2={};
matColJit = [0 0 0; 0.8 0 0];
for intType=1:intTypeNum
	%real/poiss
	vecJitEntry = [1 numel(cellLinR2{1,intType})];
	strLine = cellTypeLine{intType};
	for intJit=1:2
		%(no-)jitter
		vecData = cellfun(@(x) x(vecJitEntry(intJit)),cellLinR2(:,intType));
		Xprev=Xmean;
		Xmean=mean(vecData);
		intPlotNr = (intType-1)*2+intJit;
		errorbar(intPlotNr,Xmean,std(vecData)./sqrt(intRecNum),'x','color',matColJit(intJit,:));
		if intJit==2
			plot([intPlotNr-1 intPlotNr],[Xprev Xmean],strLine,'color',matColJit(intJit,:));
		end
		cellSubLegend2{intPlotNr} = sprintf('%s - %s',cellTypes{intType},cellJit{intJit});
	end
end
set(gca,'xtick',1:numel(cellSubLegend2),'xticklabel',cellSubLegend2);
xtickangle(gca,45);
xlim([0 numel(cellSubLegend2)+1]);
ylabel('Linearity of sd/mu (R^2)');
title('Real data is gain/timescale-invariant');

fixfig;

drawnow;
export_fig(fullpath(strFigurePath,sprintf('Q2c_CV_Invariance.tif')));
export_fig(fullpath(strFigurePath,sprintf('Q2c_CV_Invariance.pdf')));

return

%%
%{
%plot R2
for intType=1:2
	subplot(2,3,3+(intType-1)*3);hold on;
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
			vecData = matPlotR2(:,intType,intJit);
			
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
	title([cellTypes{intType} ' ;Mu +/- 95-CI of beta fit; a+(1/((b*x)^c))']);
end
fixfig;
%}