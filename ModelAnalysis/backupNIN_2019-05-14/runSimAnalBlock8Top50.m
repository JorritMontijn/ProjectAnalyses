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
	vecRunSims = [-14 -13 -12 -10 -7:-1 11];
end
for intLoadSim=vecRunSims
	clearvars -except vecDoShuff vecRunSims intLoadSim bool* vecRunAreas
	boolLoad = true;
	
	%% get simulation name [strSimulation] from [intLoadSim]
	loadSim;
	
	%% RUN: #header
	strBlockNr = getFlankedBy(mfilename,'Block','');
	strBlockNr = [strBlockNr(1) 'Top50'];
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
		
		%% dependence of dimensionality of stimulus intensity
		if ~boolAllZero
			%run analysis
			matX = matData(vecCellArea==intWithinArea,:);
			matCosSim = doDimProjTop50(matX,vecTrialStimType,matCompareTypes,vecDoShuff);
			matRawCos = matCosSim(:,:,1);
			matShuffCos = matCosSim(:,:,2);
			
			intN = size(matRawCos,1);
			strSizeN = sprintf('N%d',intN);
		
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
			title(sprintf('Area %s; Block%s,%s,%s',strPredArea,strBlockNr,strType,strSizeT),'Interpreter','none');
			fixfig;
			h=legend({'Unshuffled','Shuffled'});
			set(h,'Location','Best');
			
			%% save figure
			if boolSaveFigs
				%figure(hFigDD);
				drawnow;
				export_fig([strFigDir  'Block' strBlockNr 'Subsp_Area' num2str(intWithinArea) '_' strSizeN strSizeT strTag '.tif']);
				export_fig([strFigDir  'Block' strBlockNr 'Subsp_Area' num2str(intWithinArea) '_' strSizeN strSizeT strTag '.pdf']);
			end
			
			%% save data
			if boolSaveData
				sSave = struct;
				%ID, params
				sSave.strName = strType;
				sSave.vecRank = 1:intN;
				sSave.dTheta = range(vecStimTypeOris(matCompareTypes(1,:)));
				
				%cos sim f' projection
				sSave.matRawCos = matRawCos;
				sSave.matShuffCos = matShuffCos;
				
				%save
				save([strDataDir 'Block' strBlockNr 'Subsp_Area' num2str(intWithinArea) '_' strSizeN strSizeT strTag],'-struct','sSave');
			end
		end
	end
end
%end
