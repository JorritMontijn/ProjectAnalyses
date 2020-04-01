%% analyze input strength dependency for dimensionality of pop responses

%% initialize
close all;
clearvars;
intType = 1;
vecRunAreas = [1];

boolLoadSync = false;
boolLoadLGN = false;
boolLoadSpikeTimes = true;
boolDoSplitAnalysis = false;
%vecRunSims = [-1 -2 11:14];

if intType == 1
	boolLoadTrialSplit = true;
	vecRunSims = 21;%[20 21 1 23];
	intUseRepetitionMax = inf;
elseif intType == 2
	boolLoadTrialSplit = true;
	vecRunSims = [21];
	intUseRepetitionMax = 400;
end
for intLoadSim=vecRunSims
	clearvars -except vecRunSims intLoadSim bool* intUseRepetitionMax intUseIndepTrials vecRunAreas
	boolLoad = true;
		
	%% get simulation name [strSimulation] from [intLoadSim]
	loadSim;
	
	%% RUN: #header
	strBlockNr = getFlankedBy(mfilename,'Block','');
	strBlockNr = strBlockNr(1);
	strFigDir = ['D:\Data\Results\Block' strBlockNr '\'];
	strDataDir = ['D:\Data\Results\Data' strBlockNr '\'];
	if isempty(strBlockNr),error;end
	
	if boolLoad
		runAnalysisHeader;
	end
	
	%save(['DataForRex_' strSimulation '.mat'],'matModelResp','vecCellArea','vecTrialStimType','vecStimTypeOris','vecStimTypeOriNoise');
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
	
	%% calculate which parameter to use
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
	cellStrArea = {'1','2'};
	
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
	
	
	%define data to be used
	matData3D = double(matModelRespTS3);
	intNeurons = size(matModelRespTS3,1);
	intTrials = size(matModelRespTS3,3);
	intBins = size(matModelRespTS3,2);
	intComparisons = size(matCompareTypes,1);
	dblDiffTheta = range(vecStimTypeOris([1 2]));
	intFoldK = 10;
	dblLambda = 0;
	
	% set parameters for Fisher information
	sParamsInfo=struct;
	sParamsInfo.intFoldK = intFoldK;
	sParamsInfo.dblLambda = dblLambda;
	sParamsInfo.dblDiffTheta = dblDiffTheta;
	sParamsInfo.boolBiasCorrection = true;
	sParamsInfo.boolDirectI = true;
	sParamsInfo.boolLinear = true;
	
	%save(['DataForRex_' strSimulation '.mat'],'matModelResp','vecCellArea','vecCellTypes','vecTrialStimType');
	%% run analyses on different areas
	for intWithinArea=vecRunAreas
		%msg
		strPredArea = cellStrArea{intWithinArea};
		fprintf('Starting area %d; predicting %s [%s]\n',intWithinArea,strPredArea,getTime);
		
		%% select data
		vecUseCells = find(vecCellArea==intWithinArea);
		
		%% start
		vecUseStimTypes = unique(round(vecTrialStimType));
		indUseTrials = ismember(vecTrialStimType,vecUseStimTypes);
		
		%%{
		% create figures
		close all;
		hFigCorr = figure;
		hFigTC = figure;
		hFigPP = figure;
		
		%% PSTH
		vecPlotEpoch = [14.1 15.9];
		vecStimEpoch = [vecStimStartSecs(10)  vecStimStopSecs(10)];
		dblBinSize = 0.01;
		vecBins = vecPlotEpoch(1):dblBinSize:vecPlotEpoch(2);
		vecBinPlot = vecBins(2:end)-dblBinSize/2;
		vecNeurons = [2500:3000]-1800;
		intNumN = numel(vecNeurons);
		cellSubSpikes = cellSpikeTimesCortex(vecNeurons);
		vecStarts = cellfun(@sum,cellfun(@lt,cellSubSpikes,cellfill(vecPlotEpoch(1),size(cellSubSpikes)),'uniformoutput',false));
		vecStops = cellfun(@sum,cellfun(@lt,cellSubSpikes,cellfill(vecPlotEpoch(2),size(cellSubSpikes)),'uniformoutput',false));
		matCounts = nan(intNumN,numel(vecBins)-1);
		for intN=1:intNumN
			cellSubSpikes{intN} = cellSubSpikes{intN}(vecStarts(intN):vecStops(intN));
			matCounts(intN,:) = histcounts(cellSubSpikes{intN},vecBins);
		end
		figure;
		subplot(4,1,[1 2])
		vecMean = xmean(matCounts,1);
		vecFilt = normpdf(-4:1:4,0,2);
		vecMeanFilt = conv(vecMean,vecFilt/sum(vecFilt(:)),'same');
		vecSDFilt = conv(xstd(matCounts,1),vecFilt/sum(vecFilt(:)),'same');
		errorbar(vecBinPlot,vecMeanFilt/dblBinSize,(vecSDFilt/dblBinSize)/sqrt(intNumN));
		hold on
		plot(vecStimEpoch(1)*[1 1],[1 5],'k--');
		plot(vecStimEpoch(2)*[1 1],[1 5],'k--');
		hold off
		xlim([vecBinPlot(1) vecBinPlot(end)]);
		
		subplot(4,1,[3 4])
		implot(matCounts)
		%colormap(hot)
		
		%% correlations
		figure(hFigCorr)
		%full screen
		drawnow;
		jFig = get(handle(gcf), 'JavaFrame');
		jFig.setMaximized(true);
		figure(gcf);
		drawnow;
		
		subplot(2,3,[2 3 5 6])
		mD1=matData(vecUseCells,vecTrialStimType==1);
		mC=corr(mD1');
		imagesc(mC,[-1 1])
		xlabel('V1 Neuron')
		ylabel('V1 Neuron')
		colormap(redblue)
		colorbar
		title(sprintf('Noise Corrs, %s',strType))
		fixfig
		grid off
		
		subplot(2,3,1)
		mS = tril(true(size(mC)),-1);
		vC = mC(mS);
		histx(vC)
		xlim([-1 1]);
		xlabel('Noise correlations')
		ylabel('Number of pairs')
		title(sprintf('Noise Corrs, mean=%.3f, sd=%.3f',mean(vC),std(vC)));
		fixfig
		
		subplot(2,3,4)
		vM = mean(mD1,2)/dblStimDur;
		histx(vM)
		xlabel('Spiking rate (Hz)')
		ylabel('Number of neurons')
		title(sprintf('Mean act, mean=%.3f, sd=%.3f',mean(vM),std(vM)));
		fixfig
		
		
		%% save
		drawnow;
		export_fig([strFigDir  'Block' strBlockNr '_MeanCorr_Area' num2str(intWithinArea) '_' strTag '.tif']);
		export_fig([strFigDir  'Block' strBlockNr '_MeanCorr_Area' num2str(intWithinArea) '_' strTag '.pdf']);
		return
		
		%% example RFs
		figure;
		vecUseCells = find(vecCellArea==intWithinArea);
		cellType = {'Pyramid','Interneuron'};
		if intWithinArea == 1,matFields = matPrefGabors;intOffset=0;
		elseif intWithinArea == 2,matFields = matFieldsV2;intOffset=size(matPrefGabors,3);end
		
		vecN = sort(vecUseCells(randi(numel(vecUseCells),[1 6])));
		vecMaxVal = nan(size(vecN));
		
		for intN=1:6
			subplot(2,3,intN)
			intNeuron = vecN(intN);
			matRF = matFields(:,:,intNeuron-intOffset);
			vecMaxVal(intN) = max(abs(matRF(:)));
			imagesc(matRF,vecMaxVal(intN)*[-1 1]);
			title(sprintf('Neuron %d',intNeuron));
			colormap(redblue)
		end
		
		
		%% save
		drawnow;
		export_fig([strFigDir  'Block' strBlockNr '_ExampleFields_Area' num2str(intWithinArea) '_' strTag '.tif']);
		export_fig([strFigDir  'Block' strBlockNr '_ExampleFields_Area' num2str(intWithinArea) '_' strTag '.pdf']);
		return
		
		%% example tuning curves
		cellType = {'Pyramid','Interneuron'};
		if intWithinArea == 1
			vecUseNeurons = randperm(intCellsV1,6);
		else
			vecUseNeurons = randperm(intCellsV2,6)+intCellsV1;
		end
		vecStimOriDegrees = vecStimTypeOris(vecTrialStimType);
		[sOut] = getTuningCurves(matModelResp(vecUseNeurons,:)/dblStimDur,vecStimOriDegrees);
		
		
		%% plot
		
		for intIdx=1:size(sOut.matMeanResp,1)
			intNeuron = vecUseNeurons(intIdx);
			
			%% plot tuning curve
			figure(hFigTC);
			subplot(2,3,intIdx)
			[a,intRotate] =max(sOut.matFittedResp(intIdx,:));
			vecShift = [0 -intRotate+5];
			vecMean = circshift(sOut.matMeanResp(intIdx,:),vecShift);
			vecSD = circshift(sOut.matSDResp(intIdx,:),vecShift);
			vecFit =  circshift(sOut.matFittedResp(intIdx,:),vecShift);
			vecOriPlot = sOut.vecUniqueDegs/2;
			dblOriStep = (unique(diff(vecOriPlot)));
			vecOriPlot = [vecOriPlot vecOriPlot(end)+dblOriStep];
			
			%get type
			if ~exist('vecCellTypes','var')
				strCell = 'Unknown';
			else
				strCell = cellType{vecCellTypes(intNeuron)};
			end
			%plot
			errorbar(vecOriPlot,[vecMean vecMean(1)],[vecSD vecSD(1)]);
			hold on
			plot(vecOriPlot,[vecFit vecFit(1)]);
			hold off
			ylim([0 max(get(gca,'ylim'))]);
			set(gca,'xtick',vecOriPlot);
			set(gca,'xticklabel',vecOriPlot-90);
			xlim([-dblOriStep/2+vecOriPlot(1) dblOriStep/2+vecOriPlot(end)]);
			ylabel('Spiking rate (Hz)');
			xlabel('\theta from peak response (degs)');
			title(sprintf('Neuron %d (%s)',intNeuron,strCell));
			legend({'Mean +/- SD','von Mises'},'Location','Best')
			fixfig
			drawnow
			%pause(2)
			
			%% plot PSTH
			figure(hFigPP);
			subplot(2,3,intIdx)
			matR = squeeze(matData3D(intNeuron,:,vecTrialStimType==intRotate));
			vecMean = xmean(matR,2);
			vecSD = xstd(matR,2);
			vecBins=(1:intBins)*dblBinSizeSecs;
			
			errorbar(vecBins,vecMean,vecSD);
			ylim([0 max(matR(:))]);
			xlim([0 1.1]);
			ylabel('Spikes per bin (mean +/-)');
			xlabel('Trial bin (s)');
			title(sprintf('Neuron %d (%s)',intNeuron,strCell));
			fixfig
			
		end
		return
		%%
		export_fig([strFigDir  'AnalBlock3_TuningCurves_Area' num2str(intWithinArea) '_' strSizeXY strSizeT strTag '.tif']);
		export_fig([strFigDir  'AnalBlock3_TuningCurves_Area' num2str(intWithinArea) '_' strSizeXY strSizeT strTag '.pdf']);
		return
		%%
		export_fig([strFigDir  'AnalBlock3_PSTH_Area' num2str(intWithinArea) '_' strSizeXY strSizeT strTag '.tif']);
		export_fig([strFigDir  'AnalBlock3_PSTH_Area' num2str(intWithinArea) '_' strSizeXY strSizeT strTag '.pdf']);
		
		return
		%% example PSTHs
		vecBins=(1:intBins)*dblBinSizeSecs;
		% plot
		for intIdx=1:size(sOut.matMeanResp,1)
			
			intNeuron = vecUseNeurons(intIdx);
			matR = squeeze(matData3D(intNeuron,:,:));
			for intStim=1:8
				intStimIdx = intStim+double(intStim>4);
				subplot(3,3,intStimIdx)
				
				vecMean = xmean(matR(:,vecTrialStimType==intStim),2);
				vecSD = xstd(matR(:,vecTrialStimType==intStim),2);
				
				errorbar(vecBins,vecMean,vecSD);
				ylim([0 max(matR(:))]);
				xlim([0 1.1]);
				ylabel('Spikes per bin (mean +/-)');
			end
			subplot(3,3,5)
			if ~exist('vecCellTypes','var')
				strCell = 'Unknown';
			else
				strCell = cellType{vecCellTypes(intNeuron)};
			end
			title(sprintf('Neuron %d (%s)',intNeuron,strCell));
			pause
		end
		return
		%%}
		%% dependence of information on time bins
		intIndivOrCumul = 2;
		if intIndivOrCumul==1
			strAnalType = 'Individual';
		else
			strAnalType = 'Cumulative';
		end
		for intUseBinNr=1:intBins
			%% build data matrix
			if intIndivOrCumul == 1
				vecUseTimeBins=intUseBinNr;
			else
				vecUseTimeBins=1:intUseBinNr;
			end
			matData = squeeze(sum(matData3D(:,vecUseTimeBins,:),2));
			
			%% get data splits
			%set parameters
			sPA=struct;
			sPA.intSizeX = intSizeX;
			sPA.intSizeY = intSizeY;
			sPA.intResamplings = intIters;
			sPA.vecCellArea = vecCellArea;
			sPA.intWithinArea = intWithinArea;
			sPA.vecUseStimTypes = vecUseStimTypes;
			sPA.intMaxReps = intUseRepetitionMax;
			%get splits
			[cellMatX,cellNeuronsX,cellMatY,cellNeuronsY,cellTrials] = doDimDataSplits(matData,vecTrialStimType,sPA);
			fprintf(' Bin %d/%d: created %dx%d data splits;\n',intUseBinNr,intBins,size(cellMatX));
			
			%% perform analysis for different bins
			[matI_LogReg_CV,sOut] = doFisherFull(cellMatX,[],sParamsInfo);
			matBinInfo(:,:,intUseBinNr) = xmean(matI_LogReg_CV,3);
			fprintf('\b  Mean Fisher: %.3f [%s]\n',mean(matI_LogReg_CV(:)),getTime);
			
		end
		
		%% plot
		vecMeanInfo = squeeze(xmean(matBinInfo(:,1,:),1));
		vecSDInfo = squeeze(xstd(matBinInfo(:,1,:),1));
		errorbar((1:intBins)*dblBinSizeSecs,vecMeanInfo,vecSDInfo);
		ylim([0 max(get(gca,'ylim'))]);
		xlabel('Trial time used (s)');
		xlim([0 intBins*dblBinSizeSecs+0.01]);
		ylabel('Fisher information');
		title(sprintf('%s;Area=%d,pop=%d; %s',strAnalType,intWithinArea,intSizeX,strType),'Interpreter','none');
		fixfig
		
		drawnow;
		export_fig([strFigDir  'Block3_Information' strAnalType '_Area' num2str(intWithinArea) '_' strSizeXY strSizeT strTag '.tif']);
		export_fig([strFigDir  'Block3_Information' strAnalType '_Area' num2str(intWithinArea) '_' strSizeXY strSizeT strTag '.pdf']);
		
		
		
	end
end
