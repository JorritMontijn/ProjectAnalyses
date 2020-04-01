%% analyze input strength dependency for dimensionality of pop responses

%% initialize
intType = 2;
	close all
	clearvars -except intType
vecRunAreas = [1];

%settings
boolLoadSync = false;
boolLoadLGN = false;
boolLoadSpikeTimes = false;
boolLoadTrialSplit = false;
boolDoSplitAnalysis = false;
intSubSampleN = 101; %mean N=101
intSubSampleT = 802; %mean T=802
vecDoShuff = [0 1];%0=none, 1=reps per neuron, 2=neurons per trial, 3=xTrials per neuron

%save figs and data?
boolSaveFigs = true;
boolSaveData = true;
		
if intType == 1
	vecRunSims = [-27:-21 -14 -13 -12 -10 -7:-1 1];
end
if intType == 2
	vecRunSims = [11 20 21 22 23 34];%[21 11 22 23];
end
for intLoadSim=vecRunSims
	clearvars -except intSubSampleN intSubSampleT vecDoShuff vecRunSims intLoadSim bool* vecRunAreas
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
	
	
	%% run analyses on different areas
	for intWithinArea=vecRunAreas
		if sum(vecCellArea==intWithinArea) == 0,continue;end
		dblDiffTheta = range(vecStimTypeOris(matCompareTypes(1,:)));
		
		%msg
		strPredArea = cellStrArea{intWithinArea};
		strSizeT = strcat('T',num2str(sum(vecTrialStimType==1)));
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
		matData(indRemUnits,:) = [];
		vecCellArea(indRemUnits) = [];
		
		
		%% params
		intSamples = 100;
		intRepetitons=sum(vecTrialStimType==1);
		intTotN = min(sum(vecCellArea==intWithinArea),intRepetitons);
		if intLoadSim > 0
			vecDim = 1:intSubSampleN;
			intUseT = intSubSampleT;
			strSizeT = strcat('T',num2str(intUseT));
		else
			vecDim = 1:intTotN;
			intUseT = [];
		end
		
		%% dependence of dimensionality of stimulus intensity
		if ~boolAllZero
			%run analysis
			matX = matData(vecCellArea==intWithinArea,:);
			matFisher = doPopFisher(matX,vecTrialStimType,matCompareTypes,vecDoShuff,dblDiffTheta,vecDim,intSamples,intUseT);
			
			matRawFisher = squeeze(mean(matFisher(:,:,:,1),2));
			vecMaxRawI = matRawFisher(end,:);
			matRawFisher = bsxfun(@rdivide,matRawFisher,vecMaxRawI);
			matShuffFisher = squeeze(mean(matFisher(:,:,:,2),2));
			vecMaxShuffI = matShuffFisher(end,:);
			matShuffFisher = bsxfun(@rdivide,matShuffFisher,vecMaxShuffI);
			
			intN = size(matRawFisher,1);
			strSizeN = sprintf('N%d',intN);
			
			%run merged
			matFisherMerged = doPopFisherMerge(matX,vecTrialStimType,matCompareTypes,vecDoShuff,dblDiffTheta,vecDim,intSamples,intUseT);
			matRawFisherMerged = squeeze(mean(matFisherMerged(:,:,:,1),2));
			vecMaxRawIMerged = matRawFisherMerged(end,:);
			matRawFisherMerged = bsxfun(@rdivide,matRawFisherMerged,vecMaxRawIMerged);
			matShuffFisherMerged = squeeze(mean(matFisherMerged(:,:,:,2),2));
			vecMaxShuffIMerged = matShuffFisherMerged(end,:);
			matShuffFisherMerged = bsxfun(@rdivide,matShuffFisherMerged,vecMaxShuffIMerged);
			
			%% plot
			figure;
			dblStep = intN/10;
			vecX = round([dblStep:dblStep:intN]);
			hold on
			plot([0:intN],[0 xmean(matRawFisher,2)'],'b');
			plot([0:intN],[0 xmean(matShuffFisher,2)'],'r');
			errorbar([0 vecX],[0 xmean(matRawFisher(vecX,:),2)'],[0 xstd(matRawFisher(vecX,:),2)'],'xb');
			errorbar([0 vecX],[0 xmean(matShuffFisher(vecX,:),2)'],[0 xstd(matShuffFisher(vecX,:),2)'],'xr');
			hold off
			xlim([0 intN]);
			xlabel('Number of neurons in sample');
			ylabel('Fisher I');
			title(sprintf('Area %s; %s,%s',strPredArea,strType,strSizeT),'Interpreter','none');
			fixfig;
			h=legend({sprintf('Unshuffled,I_t=%.3f',mean(vecMaxRawI)),sprintf('Shuffled,I_t=%.3f',mean(vecMaxShuffI))});
			set(h,'Location','Best');
			
			%% save figure
			if boolSaveFigs
				%figure(hFigDD);
				drawnow;
				export_fig([strFigDir  'Block' strBlockNr '_Area' num2str(intWithinArea) '_' strSizeN strSizeT strTag '.tif']);
				export_fig([strFigDir  'Block' strBlockNr '_Area' num2str(intWithinArea) '_' strSizeN strSizeT strTag '.pdf']);
			end
			
			%% plot merged
			figure;
			dblStep = intN/10;
			vecX = round([dblStep:dblStep:intN]);
			hold on
			plot([0:intN],[0 xmean(matRawFisherMerged,2)'],'b');
			plot([0:intN],[0 xmean(matShuffFisherMerged,2)'],'r');
			errorbar([0 vecX],[0 xmean(matRawFisherMerged(vecX,:),2)'],[0 xstd(matRawFisherMerged(vecX,:),2)'],'xb');
			errorbar([0 vecX],[0 xmean(matShuffFisherMerged(vecX,:),2)'],[0 xstd(matShuffFisherMerged(vecX,:),2)'],'xr');
			hold off
			xlim([0 intN]);
			xlabel('Number of neurons in sample');
			ylabel('Fisher I');
			title(sprintf('Merged,Area %s; %s,%s',strPredArea,strType,strSizeT),'Interpreter','none');
			fixfig;
			h=legend({sprintf('Unshuffled,I_t=%.3f',mean(vecMaxRawIMerged)),sprintf('Shuffled,I_t=%.3f',mean(vecMaxShuffIMerged))});
			set(h,'Location','Best');
			
			%% save figure
			if boolSaveFigs
				%figure(hFigDD);
				drawnow;
				export_fig([strFigDir  'Block' strBlockNr 'Merged_Area' num2str(intWithinArea) '_' strSizeN strSizeT strTag '.tif']);
				export_fig([strFigDir  'Block' strBlockNr 'Merged_Area' num2str(intWithinArea) '_' strSizeN strSizeT strTag '.pdf']);
			end
			
			%% save data
			if boolSaveData
				sSave = struct;
				%ID, params
				%if ~isempty(strOriNoise)
				%	strName = ['Ori' strOriNoise strSizeN strSizeT];
				%else
					strName = strType;
				%end
				sSave.strName = strType;
				sSave.vecRank = 1:intN;
				sSave.dTheta = dblDiffTheta;
				sSave.intNeurons = size(matX,1);
				sSave.intRepetitions = sum(vecTrialStimType==matCompareTypes(1,1));
				
				%cos sim f' projection
				sSave.matRawFisher = matRawFisher;
				sSave.matShuffFisher = matShuffFisher;
				sSave.matRawFisherMerged = matRawFisherMerged;
				sSave.matShuffFisherMerged = matShuffFisherMerged;
				
				%save
				save([strDataDir 'Block' strBlockNr 'Subsp_Area' num2str(intWithinArea) '_' strSizeN strSizeT strTag],'-struct','sSave');
			end
		end
	end
end
%end
