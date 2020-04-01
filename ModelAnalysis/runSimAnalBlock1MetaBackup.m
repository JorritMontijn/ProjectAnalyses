%% perform meta plots

%% initialize
clearvars;clc;

%% set vars
strSourceDir = 'D:\Data\Processed\V1_LIFmodel\';
strFigDir = 'D:\Data\Results\Block1\';
boolSaveFigs = true;

%% load directory contents
sFiles = dir(strSourceDir);
intFiles = numel(sFiles);
strProcessingID = '';
vecProcessed = false(1,intFiles);

%% prep vars
intRunType = 2;
if intRunType == 1
	strRunType = 'FullConn';
	cellRuns = {'ContFull','LumFull'}; %with recurrent connectivity
elseif intRunType == 2
	strRunType = 'NoRecur';
	cellRuns = {'ContNone','LumNone'}; %no recurrent connectivity
end

%% loop through files
for intFile = [intFiles 1:(intFiles-1)]
	strFile = sFiles(intFile).name;
	if length(strFile) > 4 && strcmpi(strFile(end-3:end),'.mat') && ~vecProcessed(intFile)
		strSimRun = getFlankedBy(strFile,'xAreaDistributed_','_201');
		strAnalBlock = getFlankedBy(strFile,'_AB','_Area');
		strArea = getFlankedBy(strFile,'_Area','.mat');
		intArea = str2double(strArea);
		if strcmp(strAnalBlock,'1') && intArea >= 0 && ismember(strSimRun,cellRuns)
			%% load
			sLoad = load([strSourceDir strFile]);
			intAreaIdx = intArea;
			
			strParam = sLoad.strParam;
			if strcmp(strParam,'Contrast')
				intType = 1;
			elseif strcmp(strParam,'Luminance')
				intType = 2;
			end
			fprintf('Loaded file %s [%s] (%d/%d)\n',strFile,sLoad.strParam,intFile,intFiles);
			
			%% check prepro for activity matrix
			%sLoadPrePro = load([strSourceDir getFlankedBy(strFile,'','_AB1_Area') '_prepro.mat']);
			%matModelResp = sLoadPrePro.matModelResp;
			
			
			%% prep
			cellPredDimDepR2 = sLoad.cellPredDimDepR2;
			intIntensities = numel(cellPredDimDepR2);
			mapC = redbluepurple(intIntensities);
			vecIntensities = linspace(0,100,intIntensities);
			
			
			
			
			%set plot
			%subplot(2,3,(intAreaIdx)+(intType-1)*3)
			%hold on
			
			matParams = [];
			vecR2 = [];
			for intIntensity = 1:intIntensities
				%% calculate Fisher information things
				matInfo = sLoad.cellI{intIntensity};
				if isempty(matInfo)
					dblSlope=0;
					dblSlopeSD=0;
				else
					vecGroupSize = 2.^(0:(size(matInfo,2)-1));
					matGroupSize=repmat(vecGroupSize,[size(matInfo,1) 1]);
					if all(isnan(matInfo(:))),matInfo(:)=0;end
					vecX = matGroupSize(:);
					vecY = matInfo(:);
					vecNans = isnan(vecY);
					vecX(vecNans) = [];
					vecY(vecNans) = [];
					[fitobject,gof,output] = fit(vecX,vecY,'poly1');
					dblSlope = fitobject.p1;
					matCI = confint(fitobject,normcdf(1)-normcdf(-1)); %1 SD
					dblSlopeSD = range(matCI(:,1)')/2;
					
					dblR2 = gof.rsquare;
				end
				matInfoPerNeuronMean(intIntensity,intType,intAreaIdx) = dblSlope;
				matInfoPerNeuronSD(intIntensity,intType,intAreaIdx) = dblSlopeSD;
				
				%% calculate prediction and dimensionality things
				options = optimset('Algorithm','trust-region-reflective',...
					'TolFun',10^-6,'GradObj','on','MaxFunEvals',2500,'Display','off');
				matPredDimDepR2 =cellPredDimDepR2{intIntensity};
				vecDims = 1:size(matPredDimDepR2,3);
				matPredR2 = squeeze(xmean(matPredDimDepR2,2));
				if any(size(matPredR2)==1),
					intFolds=1;
					vecMeanPredR2 = matPredR2;
				else
					intFolds = size(matPredR2,1);
					vecMeanPredR2 = xmean(matPredR2,1);
				end
				matDims = meshgrid(vecDims,1:intFolds);
				
				vecFitPointsY = matPredR2(:);
				vecFitPointsX = matDims(:);
				vecNans=isnan(vecFitPointsY);
				vecFitPointsY(vecNans)=[];
				vecFitPointsX(vecNans)=[];
				if isempty(vecFitPointsY)
					vecFitP = [nan nan nan];
					dblR2 = nan;
				else
					lb = [0,0,-inf];
					ub = [1,inf,inf];
					fnExp = @(vecP,vecX) vecP(1)-(vecP(1)./exp((vecX*vecP(2)-vecP(3))));
					%[x,resnorm,residual,exitflag] = lsqcurvefit(@logisticfit00,[max(vecMeanPredR2) 1],vecFitPointsX,vecFitPointsY,lb,ub,options);
					%[x,resnorm,residual,exitflag] = lsqcurvefit(@getCumGauss,[0 1 0 1],vecFitPointsX,vecFitPointsY,lb,ub,options);
					[vecFitP,resnorm,residual,exitflag] = lsqcurvefit(fnExp,[0.3 0.1 0],vecFitPointsX,vecFitPointsY,lb,ub,options);
					
					dblSSTot = sum((vecFitPointsY - mean(vecFitPointsY(:))).^2);
					vecFit = fnExp(vecFitP,vecFitPointsX);
					vecRes = vecFitPointsY - vecFit;
					dblSSRes = sum(vecRes.^2);
					dblR2 = 1 - (dblSSRes / dblSSTot);
				end
				%save
				matParams(intIntensity,:,intAreaIdx) = vecFitP;
				vecR2(intIntensity,intAreaIdx) = dblR2;
				%plot
				%scatter(vecDims,vecMeanPredR2,[],mapC(intIntensity,:));
				%plot(vecDims,fnExp(vecFitP,vecDims),'Color',mapC(intIntensity,:))
				
				%predictability dimensionality
				vecPD = sLoad.cellPredictiveDimensions{intIntensity}(:);
				matPredDimSizeMean(intIntensity,intType,intAreaIdx) = mean(vecPD);
				matPredDimSizeSD(intIntensity,intType,intAreaIdx) = max([std(vecPD)  mean(vecPD)/100]);
				
				vecPP = sLoad.cellPredictionsR2{intIntensity}(:);
				matFullPredPerfMean(intIntensity,intType,intAreaIdx) = mean(vecPP);
				matFullPredPerfSD(intIntensity,intType,intAreaIdx) = std(vecPP);
			end
			%hold off
			%save fits
			cellParams{intType} = matParams;
			%matDimGrowth(:,intType,(intAreaIdx)) = matParams(:,intType,(intAreaIdx));
			
		end
	end
end

%% plot summaries
figure
%full screen
drawnow;
jFig = get(handle(gcf), 'JavaFrame');
jFig.setMaximized(true);
figure(gcf);
drawnow;
cellStrArea = {'V1','V2','V1-V2'};
for intArea=1:3
	%% plot pred dim size vs intensity
	subplot(2,3,intArea)
	errorbar(vecIntensities,matPredDimSizeMean(:,1,intArea),matPredDimSizeSD(:,1,intArea),'k-'); %contrast
	hold on;
	errorbar(vecIntensities,matPredDimSizeMean(:,2,intArea),matPredDimSizeSD(:,2,intArea),'k--'); %luminance
	hold off
	title(sprintf('Communication subspace %s; %s',cellStrArea{intArea},strRunType));
	xlabel('Stimulation intensity');
	ylabel('Dimensionality of predictive subspace');
	fixfig;
	legend({'Contrast','Luminance'},'Location','Best')
	xlim([-2 102]);
	ylim([0 max(get(gca,'ylim'))]);
	
	%% plot pred dim size vs predictability
	subplot(2,3,intArea+3)
	herrorbar(matFullPredPerfMean(:,1,intArea),matPredDimSizeMean(:,1,intArea),matFullPredPerfSD(:,1,intArea),'k-'); %contrast
	hold on;
	herrorbar(matFullPredPerfMean(:,2,intArea),matPredDimSizeMean(:,2,intArea),matFullPredPerfSD(:,2,intArea),'k--'); %luminance
	errorbar(matFullPredPerfMean(:,1,intArea),matPredDimSizeMean(:,1,intArea),matPredDimSizeSD(:,1,intArea),'k','LineStyle','none'); %contrast
	errorbar(matFullPredPerfMean(:,2,intArea),matPredDimSizeMean(:,2,intArea),matPredDimSizeSD(:,2,intArea),'k','LineStyle','none'); %luminance
	hold off
	title(sprintf('Communication subspace %s; %s',cellStrArea{intArea},strRunType));
	xlabel('Predictive performance full model');
	ylabel('Dimensionality of predictive subspace');
	fixfig;
	ylim([0 max(get(gca,'ylim'))]);
	xlim([0 1]);
	%legend({'Contrast','Luminance'},'Location','Best')
end
%save figure
if boolSaveFigs
	drawnow;
	export_fig([strFigDir 'AnalBlock1Meta1_' strRunType '_Area_' getDate '.tif']);
	export_fig([strFigDir 'AnalBlock1Meta1_' strRunType '_Area_' getDate '.pdf']);
end

%% plot prediction/dimensionality as function of mean activation level
figure
%full screen
drawnow;
jFig = get(handle(gcf), 'JavaFrame');
jFig.setMaximized(true);
figure(gcf);
drawnow;
for intAreaIdx=1:2
	intPlot=(intAreaIdx-1)*3;
	intArea=intAreaIdx;
	%% plot info slope vs intensity
	subplot(2,3,intPlot+1)
	errorbar(vecIntensities,matInfoPerNeuronMean(:,1,intArea),matInfoPerNeuronSD(:,1,intArea),'k-'); %contrast
	hold on;
	errorbar(vecIntensities,matInfoPerNeuronMean(:,2,intArea),matInfoPerNeuronSD(:,2,intArea),'k--'); %luminance
	hold off
	title(sprintf('Stimulus information in %s; %s',cellStrArea{intArea},strRunType));
	xlabel('Stimulation intensity');
	ylabel('Fisher information per neuron');
	fixfig;
	legend({'Contrast','Luminance'},'Location','Best')
	xlim([-2 102]);
	ylim([0 max(get(gca,'ylim'))]);
	
	%% plot info slope vs mean pred perf
	subplot(2,3,intPlot+2)
	hold on;
	herrorbar(matFullPredPerfMean(:,1,intArea),matInfoPerNeuronMean(:,1,intArea),matFullPredPerfSD(:,1,intArea),'k-'); %contrast
	herrorbar(matFullPredPerfMean(:,2,intArea),matInfoPerNeuronMean(:,2,intArea),matFullPredPerfSD(:,2,intArea),'k--'); %luminance
	errorbar(matFullPredPerfMean(:,1,intArea),matInfoPerNeuronMean(:,1,intArea),matInfoPerNeuronSD(:,1,intArea),'k','LineStyle','none'); %contrast
	errorbar(matFullPredPerfMean(:,2,intArea),matInfoPerNeuronMean(:,2,intArea),matInfoPerNeuronSD(:,2,intArea),'k','LineStyle','none'); %luminance
	hold off
	title(sprintf('Stimulus information in %s; %s',cellStrArea{intArea},strRunType));
	xlabel('Predictive performance full model');
	ylabel('Fisher information per neuron');
	fixfig;
	ylim([0 max(get(gca,'ylim'))]);
	xlim([0 1]);
	%legend({'Contrast','Luminance'},'Location','Best')
	
	%% info slope vs dimensionality
	subplot(2,3,intPlot+3)
	hold on;
	herrorbar(matPredDimSizeMean(:,1,intArea),matInfoPerNeuronMean(:,1,intArea),matPredDimSizeSD(:,1,intArea),'k-'); %contrast
	herrorbar(matPredDimSizeMean(:,2,intArea),matInfoPerNeuronMean(:,2,intArea),matPredDimSizeSD(:,2,intArea),'k--'); %luminance
	errorbar(matPredDimSizeMean(:,1,intArea),matInfoPerNeuronMean(:,1,intArea),matInfoPerNeuronSD(:,1,intArea),'k','LineStyle','none'); %contrast
	errorbar(matPredDimSizeMean(:,2,intArea),matInfoPerNeuronMean(:,2,intArea),matInfoPerNeuronSD(:,2,intArea),'k','LineStyle','none'); %luminance
	hold off
	title(sprintf('Stimulus information in %s; %s',cellStrArea{intArea},strRunType));
	xlabel('Dimensionality of predictive subspace');
	ylabel('Fisher information per neuron');
	fixfig;
	ylim([0 max(get(gca,'ylim'))]);
	xlim([0 max(get(gca,'xlim'))]);
	
	%legend({'Contrast','Luminance'},'Location','Best')
	
	
end
%save figure
if boolSaveFigs
	drawnow;
	export_fig([strFigDir 'AnalBlock1Meta2_' strRunType '_Area_' getDate '.tif']);
	export_fig([strFigDir 'AnalBlock1Meta2_' strRunType '_Area_' getDate '.pdf']);
end

