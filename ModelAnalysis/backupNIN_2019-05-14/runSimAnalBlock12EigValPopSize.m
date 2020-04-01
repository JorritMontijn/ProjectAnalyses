%% analyze input strength dependency for dimensionality of pop responses

%% initialize
intType = 1;
	close all
	clearvars -except intType
vecRunAreas = [1];

%settings
boolLoadSync = false;
boolLoadLGN = false;
boolLoadSpikeTimes = false;
boolLoadTrialSplit = false;
boolDoSplitAnalysis = false;

vecDoShuff = [0 1];%0=none, 1=reps per neuron, 2=neurons per trial, 3=xTrials per neuron

%save figs and data?
boolSaveFigs = true;
boolSaveData = true;
		
if intType == 1
	%vecPopSize = [10 20:20:100];
	vecPopSize = [25 50 100 200 400 800 1600 3200];
	%vecRunSims = [-3 -4 -5];%[1 21 22 23];
	vecRunSims = [1 2];%[1 21 22 23];
	intLimX = vecPopSize(1); %1 / end
	dblLimFracX = 0.1; %0.1 / 1
end
for intLoadSim=vecRunSims
	clearvars -except intLimX dblLimFracX vecPopSize vecDoShuff vecRunSims intLoadSim bool* vecRunAreas
	boolLoad = true;
	
	%% get simulation name [strSimulation] from [intLoadSim]
	loadSim;
	
	%% RUN: #header
	strBlockNr = getFlankedBy(mfilename,'Block','');
	strBlockNr = strBlockNr(~isnan(str2double(vec2cell(strBlockNr))) & ~eq(strBlockNr,'i') & ~eq(strBlockNr,'j'));
	strFigDir = ['D:\Data\Results\Block' strBlockNr '\'];
	strDataDir = ['D:\Data\Results\Data' strBlockNr '\'];
	if isempty(strBlockNr),error;end
	
	if boolLoad
		runAnalysisHeader;
	end
	
	%save(['DataForRex_' strSimulation '.mat'],'matModelResp','vecCellArea','vecTrialStimType','vecStimTypeOris','vecStimTypeOriNoise');
	%return
	
	%% run analyses on different areas
	cellStrArea = {'1','2'};
	for intWithinArea=vecRunAreas
		if sum(vecCellArea==intWithinArea) == 0,continue;end
		dblDiffTheta = range(vecStimTypeOris(matCompareTypes(1,:)));
		
		%msg
		strPredArea = cellStrArea{intWithinArea};
		strSizeT = strcat('T',num2str(sum(vecTrialStimType==1)));
		strLimX = strcat('X',num2str(intLimX));
		fprintf('Starting area %d; predicting %s [%s]\n',intWithinArea,strPredArea,getTime);
		
		%% start
		vecUseStimTypes = unique(matCompareTypes);
		
		%% set default parameters
		%check if no spikes
		if (sum(sum(matData(vecCellArea==1,vecTrialStimType==vecUseStimTypes(1)))) + sum(sum(matData(vecCellArea==1,vecTrialStimType==vecUseStimTypes(end))))) == 0
			boolAllZero = true;else boolAllZero = false;end
		indRemUnits = false(intCellsV1,1);
		for intStimType=(vecUseStimTypes(:)')
			indRemUnits(range(matData(:,vecTrialStimType==intStimType),2)==0)=true;
		end
		matData(indRemUnits,:) = []; %#ok<*SAGROW>
		vecCellArea(indRemUnits) = [];
		
		%% dependence of dimensionality of stimulus intensity
		if ~boolAllZero
			%% run analysis
			matX1 = matData(vecCellArea==intWithinArea,vecTrialStimType==1);
			matL = nan(numel(vecPopSize),max(vecPopSize));
			for intPopSizeIdx=1:numel(vecPopSize)
				intPopSize = vecPopSize(intPopSizeIdx);
				matCov = cov(matX1(randperm(size(matX1,1),intPopSize),:)');
				[matV,vecL] = eig(matCov,'vector'); %eigenvectors V
				[vecLsorted,vecReorder] = sort(vecL,'descend');
				matL(intPopSizeIdx,1:intPopSize) = vecLsorted;
			end
			
			
			%% get vector loadings per neuron
			if intLoadSim > 0
			matCov = cov(matX1');
			[matV,vecL] = eig(matCov,'vector'); %eigenvectors V
			[vecLsorted,vecReorder] = sort(vecL,'descend');
			matVsorted = matV(vecReorder,:);
			
			vecIdxE = vecCellTypes(vecCellArea==intWithinArea)==1;
			intE = sum(vecIdxE);
			vecIdxI = vecCellTypes(vecCellArea==intWithinArea)==2;
			intI = sum(vecIdxI);
			
			vecE1 = matVsorted(1,vecIdxE);
			vecI1 = matVsorted(1,vecIdxI);
			
			vecE2 = matVsorted(2,vecIdxE);
			vecI2 = matVsorted(2,vecIdxI);
			
			dblStep = 0.01;
			dblMax = 0.1;
			vecBins = -dblMax:dblStep:dblMax;
			vecBinsPlot = (-dblMax+dblStep/2):dblStep:(dblMax-dblStep/2);
			vecCountsE1 = histcounts(vecE1,vecBins);
			vecCountsI1 = histcounts(vecI1,vecBins);
			vecCountsE2 = histcounts(vecE2,vecBins);
			vecCountsI2 = histcounts(vecI2,vecBins);
			
			figure;
			hold on;
			plot(vecBinsPlot,vecCountsE1/intE,'r-');
			plot(vecBinsPlot,vecCountsI1/intI,'b-');
			plot(vecBinsPlot,vecCountsE2/intE,'r--');
			plot(vecBinsPlot,vecCountsI2/intI,'b--');
			hold off
			
			title(sprintf('Noise %.2f degs',mean(vecStimTypeOriNoise)),'Interpreter','none');
			legend({'Exc, PC1','Inh, PC1','Exc, PC2','Inh, PC2'});
			xlabel('PC loadings')
			ylabel('Fraction of cells')
			fixfig;
			
			%% save figure
			if boolSaveFigs
				%figure(hFigDD);
				drawnow;
				export_fig([strFigDir  'Block' strBlockNr 'LoadingsPCs_Area' num2str(intWithinArea) '_' strLimX strSizeT strTag '.tif']);
				export_fig([strFigDir  'Block' strBlockNr 'LoadingsPCs_Area' num2str(intWithinArea) '_' strLimX strSizeT strTag '.pdf']);
			end
			end
			
			%% plot
			if intLoadSim < 0
				strTitle = strType;
			else
				strTitle = sprintf('Noise=%.2f',mean(vecStimTypeOriNoise));
			end
			
			figure;
			plot(matL');
			%set(gca, 'YScale', 'log');
			xlim([0 intLimX]);
			xlabel('Component of \Sigma, sorted by Eigen value');
			%ylabel('log(\lambda)');
			ylabel('\lambda');
			h=legend(cellfun(@num2str,vec2cell(vecPopSize),'uniformoutput',false));
			set(h,'Location','Best');
			title(sprintf('%s; colours=pop sizes',strTitle),'Interpreter','none');
			fixfig;
			
			
			%% save figure
			if boolSaveFigs
				%figure(hFigDD);
				drawnow;
				export_fig([strFigDir  'Block' strBlockNr '_Area' num2str(intWithinArea) '_' strLimX strSizeT strTag '.tif']);
				export_fig([strFigDir  'Block' strBlockNr '_Area' num2str(intWithinArea) '_' strLimX strSizeT strTag '.pdf']);
			end
			
			%% plot #2, log(lambda) as function of pop size at d=25
			figure;
			intComponents = 10;%vecPopSize(4);
			matC = redbluepurple(intComponents);
			matC = matC(end:-1:1,:);
			hold on
			for intComponent=1:intComponents
				plot(vecPopSize,matL(:,intComponent),'Color',matC(intComponent,:));
			end
			colormap(matC);
			caxis([1 intComponents]);
			colorbar;
			set(gca, 'YScale', 'log', 'XScale', 'log');
			%ylabel('\lambda at i''th component of \Sigma');
			ylabel('log(\lambda) at i''th component of \Sigma');
			xlabel('Population size');
			xlim([0 max(vecPopSize)]);
			title(sprintf('%s; color=prin comp i',strTitle),'Interpreter','none');
			fixfig;
			
			
			%% save figure
			if boolSaveFigs
				%figure(hFigDD);
				drawnow;
				export_fig([strFigDir  'Block' strBlockNr 'LofP25_Area' num2str(intWithinArea) '_' strLimX strSizeT strTag '.tif']);
				export_fig([strFigDir  'Block' strBlockNr 'LofP25_Area' num2str(intWithinArea) '_' strLimX strSizeT strTag '.pdf']);
			end
			
			
			%% plot #3
			figure;
			hold all
			matC2 = redbluepurple(numel(vecPopSize));
			for intPopSizeIdx=1:numel(vecPopSize)
				intPopSize = vecPopSize(intPopSizeIdx);
				plot((1:intPopSize)/intPopSize,matL(intPopSizeIdx,1:intPopSize),'Color',matC2(intPopSizeIdx,:));
			end
			%set(gca, 'YScale', 'log');
			%ylabel('log(\lambda)');
			ylabel('\lambda');
			xlim([0 dblLimFracX]);
			xlabel('Fractional position of Eigen component of \Sigma');
			h=legend(cellfun(@num2str,vec2cell(vecPopSize),'uniformoutput',false));
			set(h,'Location','Best');
			title([strTitle '; \mu(\lambda)=' sprintf('%.2f;',nanmean(matL,2))]);
			fixfig;
			
			%% save figure
			if boolSaveFigs
				%figure(hFigDD);
				drawnow;
				export_fig([strFigDir  'Block' strBlockNr 'X' num2str(dblLimFracX) '_Area' num2str(intWithinArea) '_' strLimX strSizeT strTag '.tif']);
				export_fig([strFigDir  'Block' strBlockNr 'X' num2str(dblLimFracX) '_Area' num2str(intWithinArea) '_' strLimX strSizeT strTag '.pdf']);
			end
		end
	end
end
%end
