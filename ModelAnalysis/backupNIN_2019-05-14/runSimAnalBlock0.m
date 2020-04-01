%% analyze input strength dependency for dimensionality of pop responses

%% initialize
%close all;
clearvars;
intType = 1;
%for intType =[2 4]
	close all
	clearvars -except intType
vecRunAreas = [1];

boolLoadSync = false;
boolLoadLGN = false;
boolLoadSpikeTimes = false;
boolLoadTrialSplit = false;
boolDoRegularization = true;
boolDoSplitAnalysis = false;
%vecRunSims = [-1 -2 11:14];

	
if intType == 1
	vecRunSims = 20;%[21:26 31];
elseif intType == 2
	vecRunSims = [-10 -11 -21 -1:-1:-7];
end
for intLoadSim=vecRunSims
	clearvars -except vecRunSims boolDoRegularization intLoadSim bool* intUseRepetitionMax intUseIndepTrials vecRunAreas intSizeX intSizeY
	boolLoad = true;
	%% get simulation name [strSimulation] from [intLoadSim]
	loadSim;
	
	%% RUN: #header
	strFigDir = 'D:\Data\Results\Block0\';
	if boolLoad
		runAnalysisHeader;
	end
	
	%save(['DataForRex_' strSimulation '.mat'],'matModelResp','vecCellArea','vecTrialStimType');
	%return
	
	cellIn = strsplit(strSimulation,'_');
	strFiller = cellIn{1};
	strType = cell2mat(cellIn(2:(end-1)));
	strDate = cellIn{end};
	strTag = [strType '_' strDate];
	if boolDoSplitAnalysis
		strTag = ['TS_' strTag];
		strType = ['TS_' strType];
	end
	
	
	%% set parameters
	boolShowThree = true; %show only pred, rand, full CV
	boolIncludeShuff = false;
	boolBC = false;
	if ~exist('boolDoRegularization','var'),boolDoRegularization = true;end
	intFoldK = 10;
	if ~exist('intUseRepetitionMax','var'),intUseRepetitionMax = 4000;end
	if ~exist('intSizeX','var'),intSizeX = 80;end
	if ~exist('intSizeY','var'),intSizeY = 30;end
	boolSaveFigs = true;
	boolSaveData = false;
	dblLambdaInfo = 0;
	dblLambdaPred = 0;
	dblNoiseLevel = 0;
	intMaxDimAnal = min(30,intSizeY);
	%intWithinArea = 1; %1, within V1; 2, within V2; 3, across V1-V2;
	intIters = 50;
	
	if ~boolDoRegularization
		strTag = ['NoR_' strTag];
		strType = ['NoR_' strType];
	end
	
	%% calculate which parameter to use
	strSizeXY = sprintf('X%dY%d',intSizeX,intSizeY);
	strSizeT = strcat('T',num2str(intUseRepetitionMax));
	vecRanges = range(matStimTypeCombos(2:end,:),2);
	[d,intUseParam]=max(vecRanges); %exclude orientation
	intUseParam = intUseParam + 1;
	if all(vecRanges==0), intUseParam = 0;vecParamStimType = [0];strParam='None';
		vecParamIdx = ones(1,size(matStimTypeCombos,2));
		vecParamVals = 1;
		intParamNum=1;
		mapC = redbluepurple(intParamNum);
	else
		if intUseParam == 1,vecParamStimType = vecStimTypeOris;strParam='Ori';
		elseif intUseParam == 2,vecParamStimType = vecStimTypeSFs;strParam='SF';
		elseif intUseParam == 3,vecParamStimType = vecStimTypeTFs;strParam='TF';
		elseif intUseParam == 4,vecParamStimType = vecStimTypeContrasts;strParam='Contrast';
		elseif intUseParam == 5,vecParamStimType = vecStimTypeLuminance;strParam='Luminance';
		end
		vecParamIdx = matStimTypeCombos(intUseParam,:);
		vecParamVals = unique(vecParamStimType);
		intParamNum=numel(vecParamVals);
		mapC = redbluepurple(intParamNum);
	end
	cellStrArea = {'V1 from V1','V2 from V2','V2 from V1'};
	
	%% build comparison matrix for stim types
	vecUniqueStimTypes = unique(vecTrialStimType);
	intMultiComp = false;
	intStimTypes = numel(vecUniqueStimTypes);
	if ~exist('matCompareTypes','var')
		if sum(range(matStimTypeCombos,2)>0)>1
			intMultiComp = true;
			matCompareTypes = [ones(1,size(matStimTypeCombos,2)-1)' (2:size(matStimTypeCombos,2))'];
			matCompareTypes = [1 2; 1 3; 1 4];
		else
			
			matCompareTypes = [(1:intStimTypes)' circshift((1:intStimTypes)',-1)];
		end
	end
	
	%% get data
	if boolLoadTrialSplit
		fprintf('Transforming data [%s]\n',getTime);
		[vecTrialStimTypeUnsplit,matRespUnsplit,vecTrialStimTypeSplit,matRespSplitNorm,matRespSplit,matSplitTrialIdx] = ...
			prepSplitData(matModelRespTS3,vecTrialStimType);
		
		if boolDoSplitAnalysis
			matData = matRespSplitNorm;
			matDataFull = matRespSplit;
			vecTrialStimType = vecTrialStimTypeSplit;
		else
			matData = matRespUnsplit;
			vecTrialStimType = vecTrialStimTypeUnsplit;
		end
	else
		matData = double(matModelResp);
		vecTrialStimType = vecTrialStimType;
	end
	%[~,~,dblInfoCheck]=getSeparation(matData,vecTrialStimType,false,range(vecStimTypeOris([1 2])))
	
	%% run analyses on different areas
	for intWithinArea=vecRunAreas
		%msg
		strPredArea = cellStrArea{intWithinArea};
		fprintf('Starting area %d; predicting %s [%s]\n',intWithinArea,strPredArea,getTime);
		
		
		%% start
		vecUseStimTypes = unique(matCompareTypes);
		indUseTrials = ismember(vecTrialStimType,vecUseStimTypes);
		
		%% set default parameters
		dblDiffTheta = range(vecStimTypeOris(matCompareTypes(1,:)));
		
		%check if no spikes
		if (sum(sum(matData(1:intCellsV1,vecTrialStimType==vecUseStimTypes(1)))) + sum(sum(matData(1:intCellsV1,vecTrialStimType==vecUseStimTypes(end))))) == 0
			boolAllZero = true;else boolAllZero = false;end
		indRemUnits = false(intCellsV1,1);
		for intStimType=(vecUseStimTypes(:)')
			indRemUnits(range(matData(:,vecTrialStimType==intStimType),2)==0)=true;
		end
		matData(indRemUnits,:) = [];
		if exist('matDataFull','var')
			matDataFull(indRemUnits,:) = [];
		end
		vecCellArea(indRemUnits) = [];
		
		%% do analysis
		%compare signal stationarity with stationary Poisson and Gaussian responses
		vecMean=xmean(matData,2);
		vecVar=xvar(matData,2);
		[intCells,intTrials] = size(matData);
		matPoisson = poissrnd(repmat(vecMean,[1 intTrials]));
		matGauss = normrnd(repmat(vecMean,[1 intTrials]),repmat(sqrt(vecVar),[1 intTrials]));
		
		matFFT_Data = nan(intCells,intTrials);
		matFFT_Gauss = nan(intCells,intTrials);
		matFFT_Poisson = nan(intCells,intTrials);
		for intCell=1:intCells
			matFFT_Data(intCell,:) = abs(fft(matData(intCell,:)));
			matFFT_Gauss(intCell,:) = abs(fft(matGauss(intCell,:)));
			matFFT_Poisson(intCell,:) = abs(fft(matPoisson(intCell,:)));
		end
		matFFT_Data = matFFT_Data/sum(matFFT_Data(:));
		matFFT_Gauss = matFFT_Gauss/sum(matFFT_Gauss(:));
		matFFT_Poisson = matFFT_Poisson/sum(matFFT_Poisson(:));
		
		figure
		vecPlot = (1:intTrials)/intTrials;
		errorbar(vecPlot,xmean(matFFT_Data,1),xstd(matFFT_Data,1)/sqrt(intCells));
		hold on
		errorbar(vecPlot,xmean(matFFT_Gauss,1),xstd(matFFT_Gauss,1)/sqrt(intCells));
		errorbar(vecPlot,xmean(matFFT_Poisson,1),xstd(matFFT_Poisson,1)/sqrt(intCells));
		legend({'Data','Gauss','Poisson'})
		
		set(gca,'YScale','log')
		set(gca,'XScale','log')
		title(strSimulation,'Interpreter','none');
		xlabel('Frequency (1/trials)');
		ylabel('Normalized power (integral=1)');
		fixfig
		grid off;
		
		%save
		if boolSaveFigs
			%figure(hFigDD);
			drawnow;
			export_fig([strFigDir  'Block0_Fourier_Area' num2str(intWithinArea) '_' strTag '.tif']);
			export_fig([strFigDir  'Block0_Fourier_Area' num2str(intWithinArea) '_' strTag '.pdf']);
		end
	end
end
%end
