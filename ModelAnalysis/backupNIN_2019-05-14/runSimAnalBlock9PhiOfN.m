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
	vecRunSims = [-3 -4 -5];
end
if intType == 2
	vecRunSims = [11 20 21 22 23 34];%[21 11 22 23];
end
if intType == 3
	vecRunSims = 21;%[20 21 22 23];%[21 11 22 23];
end
for intLoadSim=vecRunSims
	clearvars -except intSubSample* vecNeuronN vecDoShuff vecRunSims intLoadSim bool* vecRunAreas
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
	
	%check if subdir exists, otherwise create
	sDir=dir([strFigDir strType '\']);
	if isempty(sDir)
		mkdir(strFigDir,strType);
	end
	strFigDir = [strFigDir strType '\'];
	
	sDir=dir([strDataDir strType '\']);
	if isempty(sDir)
		mkdir(strDataDir,strType);
	end
	strDataDir = [strDataDir strType '\'];
	
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
		
		%set which sizes to run
		vecNeuronN = 20:20:min(500, sum(vecCellArea==intWithinArea));
		
		%% params
		for intNeuronIdx = 1:numel(vecNeuronN)
			intUseN = vecNeuronN(intNeuronIdx);
			
			intSamples = 100;
			intRepetitons=sum(vecTrialStimType==1);
			intComparisons = size(matCompareTypes,1);
			intUseT = inf;
			fprintf('Running sample size %d (%d/%d) [%s]\n',intUseN,intNeuronIdx,numel(vecNeuronN),getTime);
			
			%% dependence of dimensionality of stimulus intensity
			if ~boolAllZero
				%cellRawCos = cell
				matRawCos3 = nan(intUseN,intComparisons,intSamples);
				matShuffCos3 = nan(intUseN,intComparisons,intSamples);
				matRawCosMerged3 = nan(intUseN,intComparisons,intSamples);
				matShuffCosMerged3 = nan(intUseN,intComparisons,intSamples);
				
				for intSample=1:intSamples
					
					
					%subsample?
					matX = matData(vecCellArea==intWithinArea,:);
					vecUseTrials = 1:size(matX,2);
					intMaxN = sum(vecCellArea==intWithinArea);
					
					if intUseN < intMaxN
						matX = matX(randperm(intMaxN,intUseN),:);
					end
					intMaxT = sum(vecTrialStimType==1);
					
					if intUseT < intMaxT
						vecUseTrials = [];
						for intStimType=1:intStimTypes
							vecT = find(vecTrialStimType==intStimType);
							vecUseTrials = [vecUseTrials vecT(randperm(numel(vecT),intSubSampleT))];
						end
						matX = matX(:,vecUseTrials);
					end
					
					%run analysis
					boolVerbose = false;
					matCosSim = doDimProj(matX,vecTrialStimType(vecUseTrials),matCompareTypes,vecDoShuff,boolVerbose);
					matRawCos3(:,:,intSample) = matCosSim(:,:,1);
					matShuffCos3(:,:,intSample) = matCosSim(:,:,2);
					
					%run merged
					matCosSim = doDimProjMerge(matX,vecTrialStimType(vecUseTrials),matCompareTypes,vecDoShuff,boolVerbose);
					matRawCosMerged3(:,:,intSample) = matCosSim(:,:,1);
					matShuffCosMerged3(:,:,intSample) = matCosSim(:,:,2);
				end
				%calc means
				matRawCos = mean(matRawCos3,3);
				matShuffCos = mean(matShuffCos3,3);
				matRawCosMerged = mean(matRawCosMerged3,3);
				matShuffCosMerged = mean(matShuffCosMerged3,3);
				
				intN = size(matRawCos,1);
				strSizeN = sprintf('N%d',intN);
				strSizeT = strcat('T',num2str(sum(vecTrialStimType(vecUseTrials)==1)));
				
				%% plot
				figure;
				dblStep = intN/10;
				vecX = round([dblStep:dblStep:intN]);
				hold on
				plot([0:intN],[0 xmean(matRawCos,2)'],'b');
				plot([0:intN],[0 xmean(matShuffCos,2)'],'r');
				errorbar([0 vecX],[0 xmean(matRawCos(vecX,:),2)'],[0 xstd(matRawCos(vecX,:),2)'],'xb');
				errorbar([0 vecX],[0 xmean(matShuffCos(vecX,:),2)'],[0 xstd(matShuffCos(vecX,:),2)'],'xr');
				hold off
				ylim([0 min(max(get(gca,'ylim')),1)]);
				xlim([0 intN]);
				xlabel('Eigen-subspace of \Sigma, ordered by \lambda');
				ylabel('cos(f'',\Sigma)');
				title(sprintf('Area %s; %s,%s, %s',strPredArea,strSizeN,strSizeT,strType),'Interpreter','none');
				fixfig;
				h=legend({'Unshuffled','Shuffled'});
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
				plot([0:intN],[0 xmean(matRawCosMerged,2)'],'b');
				plot([0:intN],[0 xmean(matShuffCosMerged,2)'],'r');
				errorbar([0 vecX],[0 xmean(matRawCosMerged(vecX,:),2)'],[0 xstd(matRawCosMerged(vecX,:),2)'],'xb');
				errorbar([0 vecX],[0 xmean(matShuffCosMerged(vecX,:),2)'],[0 xstd(matShuffCosMerged(vecX,:),2)'],'xr');
				hold off
				ylim([0 min(max(get(gca,'ylim')),1)]);
				xlim([0 intN]);
				xlabel('Eigen-subspace of \Sigma, ordered by \lambda');
				ylabel('cos(f'',\Sigma)');
				title(sprintf('Merged,Area %s; %s,%s, %s',strPredArea,strSizeN,strSizeT,strType),'Interpreter','none');
				fixfig;
				h=legend({'Unshuffled','Shuffled'});
				set(h,'Location','Best');
				
				%% save figure
				if boolSaveFigs
					%figure(hFigDD);
					drawnow;
					export_fig([strFigDir  'Block' strBlockNr 'Merged_Area' num2str(intWithinArea) '_' strSizeN strSizeT strTag '.tif']);
					export_fig([strFigDir  'Block' strBlockNr 'Merged_Area' num2str(intWithinArea) '_' strSizeN strSizeT strTag '.pdf']);
				end
				close all;
				
				%% save data
				if boolSaveData
					sSave = struct;
					strOriNoise = getFlankedBy(strType,'Ori','');
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
					sSave.matRawCos = matRawCos;
					sSave.matShuffCos = matShuffCos;
					sSave.matRawCosMerged = matRawCosMerged;
					sSave.matShuffCosMerged = matShuffCosMerged;
					
					%save
					save([strDataDir 'Block' strBlockNr '_Area' num2str(intWithinArea) '_' strSizeN strSizeT strTag],'-struct','sSave');
				end
			end
		end
	end
	
	%% run meta analysis
	runMetaBlock9PhiOfN;
end
%end
