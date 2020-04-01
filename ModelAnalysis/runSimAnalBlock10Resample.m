%% analyze input strength dependency for dimensionality of pop responses

%% initialize
intType = 1;
	close all
	clearvars -except intType
vecRunAreas = [1];

%experimental N&T
%name,neurons,trials
cellMeta = ...
	{'106r01','113','400';...
	'107l03','112','400';...
	'p113al','128','200';...
	'p115al','142','200';...
	%'JGP011','98','120';...
	%'JGP012','72','72';...
	%'JGP014','96','96';...
	%'JGP016','65','108';...
	%'JGS006','75','88';...
	'R128hs','92','2000';...
	'R129hs','75','800';...
	'cadp13','113','978';...
	'monp31','121','1305';...
	'monp33','103','640'};
vecNeuronN = cellfun(@str2double,cellMeta(:,2))';
vecTrialT = cellfun(@str2double,cellMeta(:,3))';
intResamples = 1000;

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
	vecRunSims = [20];%[-14 -13 -12 -10 -7:-1 11];
end
for intLoadSim=vecRunSims
	clearvars -except intResamples vecNeuronN vecTrialT vecDoShuff vecRunSims intLoadSim bool* vecRunAreas
	boolLoad = true;
	
	%% get simulation name [strSimulation] from [intLoadSim]
	loadSim;
	
	%% RUN: #header
	strBlockNr = getFlankedBy(mfilename,'Block','');
	strBlockNr = strBlockNr(~isnan(str2double(vec2cell(strBlockNr))));
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
			%% loop through samples for normal run & merged
			intComps = size(matCompareTypes,1);
			matPhi = nan(intComps,intResamples);
			matPhiMerged = nan(intComps,intResamples);
			vecSizeN = nan(intComps,intResamples);
			vecSizeT = nan(intComps,intResamples);
			for intSample=1:intResamples
				%msg
				fprintf('Running sample %d/%d [%s]\n',intSample,intResamples,getTime);
				
				%vars
				intSubSampleN = vecNeuronN(randi(numel(vecNeuronN)));
				intSubSampleT = vecTrialT(randi(numel(vecTrialT)));
				vecSizeN(intSample) = intSubSampleN;
				vecSizeT(intSample) =intSubSampleT;
				
				
				%subsample?
				intMaxN = sum(vecCellArea==intWithinArea);
				matX = matData(vecCellArea==intWithinArea,:);
				if intSubSampleN < intMaxN
					matX = matX(randperm(intMaxN,intSubSampleN),:);
				end
				intMaxT = sum(vecTrialStimType==1);
				vecUseTrials = 1:size(matX,2);
				if intSubSampleT < intMaxT
					vecUseTrials = [];
					for intStimType=1:intStimTypes
						vecT = find(vecTrialStimType==intStimType);
						vecUseTrials = [vecUseTrials vecT(randperm(numel(vecT),intSubSampleT))];
					end
					matX = matX(:,vecUseTrials);
				end
				
				%run analysis
				matCosSim = doDimProj(matX,vecTrialStimType(vecUseTrials),matCompareTypes,vecDoShuff,false);
				matRawCos = matCosSim(:,:,1);
				matShuffCos = matCosSim(:,:,2);
				
				%calc phi
				matPhi(:,intSample) = sum((matRawCos - matShuffCos),1)/intSubSampleN;% ./ (1 - matRelShuff + 1/1000);
				
				%run merged
				matCosSimMerged = doDimProjMerge(matX,vecTrialStimType(vecUseTrials),matCompareTypes,vecDoShuff,false);
				matRawCosMerged = matCosSimMerged(:,:,1);
				matShuffCosMerged = matCosSimMerged(:,:,2);
				
				%calc phi
				matPhiMerged(:,intSample) = sum((matRawCosMerged - matShuffCosMerged),1)/intSubSampleN;% ./ (1 - matRelShuff + 1/1000);
				
			end
			
			%% plot
			figure;
			dblStep = 0.02;
			vecX = -0.2:dblStep:0.5;
			vecBins = [vecX-dblStep/2 vecX(end)+dblStep/2];
			[vecCounts] = histcounts(matPhi(:),vecBins);
			[vecCountsMerged] = histcounts(matPhiMerged(:),vecBins);
			
			hold on
			plot(vecX,vecCounts,'b');
			plot(vecX,vecCountsMerged,'r');
			hold off
			xlim([min(vecX) max(vecX)]);
			xlabel('Diff. corr. strength, \Phi (AUC) for cos(f'',\Sigma)');
			ylabel('Count');
			title(sprintf('Area %s; %d samples, %s',strPredArea,intResamples,strType),'Interpreter','none');
			h=legend({'Normal','Merged'});
			set(h,'Location','Best');
			fixfig;
			
			%% save figure
			if boolSaveFigs
				%figure(hFigDD);
				drawnow;
				export_fig([strFigDir  'Block' strBlockNr 'Subsp_Area' num2str(intWithinArea) '_' strTag '.tif']);
				export_fig([strFigDir  'Block' strBlockNr 'Subsp_Area' num2str(intWithinArea) '_' strTag '.pdf']);
			end
			
			%% save data
			if boolSaveData
				sSave = struct;
				strOriNoise = getFlankedBy(strType,'Ori','');
				if ~isempty(strOriNoise)
					strName = ['Ori' strOriNoise];
				else
					strName = strType;
				end
				%ID, params
				sSave.strName = strName;
				sSave.dTheta = range(vecStimTypeOris(matCompareTypes(1,:)));
				
				%cos sim f' projection
				sSave.matPhiMerged = matPhiMerged;
				sSave.matPhi = matPhi;
				sSave.vecSizeN = vecSizeN;
				sSave.vecSizeT = vecSizeT;
				
				%save
				save([strDataDir 'Block' strBlockNr 'Subsp_Area' num2str(intWithinArea) '_' strTag],'-struct','sSave');
			end
		end
	end
end
%end
