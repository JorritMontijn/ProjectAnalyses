%% analyze input strength dependency for dimensionality of pop responses

%% initialize
intType = 5;
	%close all
	clearvars -except intType
vecRunAreas = [1];

%experimental N&T
%name,neurons,trials
cellMeta = ...
	{'monp31','121','1305';...
	'monp33','103','640';...
	'cadp13','113','978';...
	'cadp185','88','680'
	'cadp428','84','548';...
	'cadp429','73','698';...
	'alcp42','119','555',;...
	};

vecNeuronN = cellfun(@str2double,cellMeta(:,2))';
vecTrialT = cellfun(@str2double,cellMeta(:,3))';

%settings
boolLoadSync = false;
boolLoadLGN = false;
boolLoadSpikeTimes = false;
boolLoadTrialSplit = false;
boolDoSplitAnalysis = false;
intSubSampleN = 73;%100;%round(mean(vecNeuronN)); %mean N=100
intSubSampleT = inf;%805;%805;%round(mean(vecTrialT)); %mean T=805
vecDoShuff = [0 1];%0=none, 1=reps per neuron, 2=neurons per trial, 3=xTrials per neuron
boolUseZscore = true;
intSamples = 100;

%save figs and data?
boolSaveFigs = true;
boolSaveData = true;
	
if intType == 1
	vecRunSims = -14;%[-3 -4 -8 -11 -12 -13 -14];%[-27:-21 -14 -13 -12 -10 -7:-1];
end
if intType == 2
	vecRunSims = [100 110 150 1 2 102 104 106 108];%[21 11 22 23];
end
if intType == 3
	vecRunSims = -13;%[20 21 22 23];%[21 11 22 23];
end
if intType == 4
	vecRunSims =[ 102:2:110 150];%[112:2:120];%200:2:220;
end
if intType == 5
	vecRunSims =[100];%[112:2:120];%200:2:220;
end
for intLoadSim=vecRunSims
	clearvars -except intSubSample* vecDoShuff vecRunSims intLoadSim bool* vecRunAreas intSamples
	boolLoad = true;
	
	%% get simulation name [strSimulation] from [intLoadSim]
	loadSim;
	
	%% RUN: #header
	strBlockNr = getFlankedBy(mfilename,'Block','');
	strBlockNr = strBlockNr(1);
	strFigDir = ['F:\Data\Results\SimFigs\Block' strBlockNr '\'];
	strDataDir = ['F:\Data\Results\SimFigs\Data' strBlockNr '\'];
	if isempty(strBlockNr),error;end
	
	if boolLoad
		runAnalysisHeader;
	end
	intRepetitions = min(sum(bsxfun(@eq,vecTrialStimType,vecUniqueStimTypes'),2));
	intNeurons = size(matData,1);
	strType;

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
	%if ~exist('matCompareTypes','var')
		%if sum(range(matStimTypeCombos,2)>0)>1
		%	intMultiComp = true;
		%	matCompareTypes = [ones(1,size(matStimTypeCombos,2)-1)' (2:size(matStimTypeCombos,2))'];
		%	matCompareTypes = [1 2; 1 3; 1 4];
		%else
			%matStimTypeCombos
			%matCompareTypes = [(1:intStimTypes)' circshift((1:intStimTypes)',-1)];
		%end
	%end
	
	%% check if monkey is cadet and remove 0 degrees
	if strfind(strType,'cadet')
		vecRemStim = find(matStimTypeCombos(1,:)==0);
	else
		vecRemStim = [];
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
		vecDiffTheta = range(vecStimTypeOris(matCompareTypes),2);
		indRem = vecDiffTheta>22.5;
		if ~isempty(vecRemStim)
			indRem0 = any(ismember(matCompareTypes,vecRemStim),2);
		else
			indRem0 = false(size(indRem));
		end
		matCompareTypes(indRem | indRem0,:) = [];
		dblDiffTheta = mean(vecDiffTheta(~indRem));
		%add reverse comparisons
		vecRevIdx = [2 1];
		for intComp=1:size(matCompareTypes,1)
			vecComp = matCompareTypes(intComp,:);
			vecRev = vecComp(vecRevIdx);
			if ~any(all(bsxfun(@eq,matCompareTypes,vecRev),2))
				%comparison does not exist yet
				matCompareTypes(end+1,:) = vecRev;
			end
		end
		matCompareTypes
		
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
		
		%% params
		intRepetitons=sum(vecTrialStimType==1);
		intComparisons = size(matCompareTypes,1);
		%if intLoadSim > 0
			intUseN = min(sum(vecCellArea==intWithinArea),intSubSampleN);
			intUseT = intSubSampleT;
		%else
		%	intUseN = min(sum(vecCellArea==intWithinArea),intRepetitons);
		%	intUseT = inf;
		%end
		
		%% dependence of dimensionality of stimulus intensity
		if ~boolAllZero
			%cellRawCos = cell
			matRawCos3 = nan(intUseN,intComparisons,intSamples);
			matShuffCos3 = nan(intUseN,intComparisons,intSamples);
			matRawCosMerged3 = nan(intUseN,intComparisons,intSamples);
			matShuffCosMerged3 = nan(intUseN,intComparisons,intSamples);
			matRawCosRandom3 = nan(intUseN,intComparisons,intSamples);
			matShuffCosRandom3 = nan(intUseN,intComparisons,intSamples);
			
			for intSample=1:intSamples
				fprintf('Running sample %d/%d [%s]\n',intSample,intSamples,getTime);
				
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
				
				%% run analysis
				boolVerbose = false;
				matCosSim = doDimProj(matX,vecTrialStimType(vecUseTrials),matCompareTypes,vecDoShuff,boolVerbose,boolUseZscore);
				matRawCos3(:,:,intSample) = matCosSim(:,:,1);
				matShuffCos3(:,:,intSample) = matCosSim(:,:,2);
				
				%run merged
				matCosSimMerged = doDimProjMerge(matX,vecTrialStimType(vecUseTrials),matCompareTypes,vecDoShuff,boolVerbose,boolUseZscore);
				matRawCosMerged3(:,:,intSample) = matCosSimMerged(:,:,1);
				matShuffCosMerged3(:,:,intSample) = matCosSimMerged(:,:,2);
				
				%run random
				matCosSimRand = doDimProjRand(matX,vecTrialStimType(vecUseTrials),matCompareTypes,vecDoShuff,boolVerbose,boolUseZscore);
				matRawCosRandom3(:,:,intSample) = matCosSimRand(:,:,1);
				matShuffCosRandom3(:,:,intSample) = matCosSimRand(:,:,2);
				
			end
			%calc means
			matRawCos = nanmean(matRawCos3,3);
			matShuffCos = nanmean(matShuffCos3,3);
			matRawCosMerged = nanmean(matRawCosMerged3,3);
			matShuffCosMerged = nanmean(matShuffCosMerged3,3);
			matRawCosRandom = nanmean(matRawCosRandom3,3);
			matShuffCosRandom = nanmean(matShuffCosRandom3,3);
			
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
				export_fig([strFigDir  'Block' strBlockNr '_Area' num2str(intWithinArea) '_' strSizeN strSizeT strTag '.jpg']);
				export_fig([strFigDir  'Block' strBlockNr '_Area' num2str(intWithinArea) '_' strSizeN strSizeT strTag '.pdf']);
			end
			
			%% plot distribution
			matRawCos = nanmean(matRawCos3,3);
			matShuffCos = nanmean(matShuffCos3,3);
			matDiff = matRawCos3 - matShuffCos3;
			vecPhi = flat(sum(matDiff,1) ./ size(matDiff,1));
			
			
			
			%{
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
			
			%% plot random
			figure;
			dblStep = intN/10;
			vecX = round([dblStep:dblStep:intN]);
			hold on
			plot([0:intN],[0 xmean(matRawCosRandom,2)'],'b');
			plot([0:intN],[0 xmean(matShuffCosRandom,2)'],'r');
			errorbar([0 vecX],[0 xmean(matRawCosRandom(vecX,:),2)'],[0 xstd(matRawCosRandom(vecX,:),2)'],'xb');
			errorbar([0 vecX],[0 xmean(matShuffCosRandom(vecX,:),2)'],[0 xstd(matShuffCosRandom(vecX,:),2)'],'xr');
			hold off
			ylim([0 min(max(get(gca,'ylim')),1)]);
			xlim([0 intN]);
			xlabel('Eigen-subspace of \Sigma, ordered by \lambda');
			ylabel('cos(f'',\Sigma)');
			title(sprintf('Random,Area %s; %s,%s, %s',strPredArea,strSizeN,strSizeT,strType),'Interpreter','none');
			fixfig;
			h=legend({'Unshuffled','Shuffled'});
			set(h,'Location','Best');
			
			%% save figure
			if boolSaveFigs
				%figure(hFigDD);
				drawnow;
				export_fig([strFigDir  'Block' strBlockNr 'Random_Area' num2str(intWithinArea) '_' strSizeN strSizeT strTag '.tif']);
				export_fig([strFigDir  'Block' strBlockNr 'Random_Area' num2str(intWithinArea) '_' strSizeN strSizeT strTag '.pdf']);
			end
			%}
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
				if exist('vecStimTypeOriNoise','var')
					sSave.dblNoise = vecStimTypeOriNoise(1);
				else
					sSave.dblNoise = nan;
				end
				
				%cos sim f' projection
				sSave.matRawCos = matRawCos;
				sSave.matShuffCos = matShuffCos;
				sSave.matRawCosMerged = matRawCosMerged;
				sSave.matShuffCosMerged = matShuffCosMerged;
				sSave.matRawCosRandom = matRawCosRandom;
				sSave.matShuffCosRandom = matShuffCosRandom;
				
				%save
				save([strDataDir 'Block' strBlockNr '_Area' num2str(intWithinArea) '_' strSizeN strSizeT strTag],'-struct','sSave');
			end
		end
	end
end
%end
