%% ccg analysis
%v) make ccgs in NOT to check for local connections

%% load data
strDataPath = 'E:\DataPreProcessed';
sFiles = dir(fullpath(strDataPath,'*_AP.mat'));
if ~exist('sExp','var') || isempty(sExp)
	sExp = [];
	for intFile=1:numel(sFiles)
		fprintf('Loading %d/%d: %s [%s]\n',intFile,numel(sFiles),sFiles(intFile).name,getTime);
		sLoad = load(fullpath(sFiles(intFile).folder,sFiles(intFile).name));
		sAP = sLoad.sAP;
		if ~isfield(sAP,'sPupil')
			sAP.sPupil = [];
		end
		
		if isempty(sExp)
			sExp = sAP;
		else
			sExp(end+1) = sAP;
		end
		
	end
end

strTargetPath = 'D:\Data\Results\AlbinoProject';

%best rec BL6: 20191216B5 (rec 17)
%best rec DBA: 20210212B2 (rec 5)

%% define area categories
%cortex
cellUseAreas = [];
cellUseAreas{1} = {'Primary visual','Posteromedial visual','anteromedial visual'};
%NOT
cellUseAreas{2} = {'nucleus of the optic tract'};
cellAreaGroups = {'Vis. ctx','NOT'};
cellAreaGroupsAbbr = {'Ctx','NOT'};
cellSubjectGroups = {'BL6','DBA'};
intUseAreaNum = numel(cellUseAreas);
vecColAlb = [0.9 0 0];
vecColBl6 = lines(1);

%% pre-allocate
cellCCG = cell(2,2);

% run
cellNameAP = arrayfun(@(x) x.sJson.file_preproAP,sExp,'uniformoutput',false);
cellExperiment = arrayfun(@(x) x.sJson.experiment,sExp,'uniformoutput',false);
cellRemove = {};%{'RecMA5_2021-03-01R01_g0_t0'};
indRemRecs = contains(cellExperiment,cellRemove);
cellSubjectType = arrayfun(@(x) x.sJson.subjecttype,sExp,'uniformoutput',false);
for intSubType=1:2
	intPopCounter = 0;
	if intSubType == 1
		intBestRec = 17;
		strSubjectType = 'BL6';
	elseif intSubType == 2
		intBestRec = 5;
		strSubjectType = 'DBA';
	end
	indUseRecs = contains(cellSubjectType,strSubjectType);
	vecRunRecs = find(indUseRecs & ~(indRemRecs));
	%vecRunRecs = intBestRec;
	for intRecIdx=1:numel(vecRunRecs)
		intRec=vecRunRecs(intRecIdx)
		sRec = sExp(intRec);
		strName=[sRec.sJson.subject '_' sRec.sJson.date];
		
		% split cells into areas
		cellCellsPerArea = cell(1,numel(cellUseAreas));
		cellAreasPerCluster = {sRec.sCluster.Area};
		for intArea=1:numel(cellUseAreas)
			cellCellsPerArea{intArea} = contains(cellAreasPerCluster(:),cellUseAreas{intArea},'IgnoreCase',true);
		end
		vecCellsNrPerArea = cellfun(@sum,cellCellsPerArea);
		
		%get waveform in areas
		for intArea=1:intUseAreaNum
			%include?
			indUseCells = arrayfun(@(x) x.Violations1ms < 0.25 & abs(x.NonStationarity) < 0.25,sRec.sCluster(:));
			
			%build cell vectors
			vecSelectCells = find(indUseCells(:) & cellCellsPerArea{intArea}(:));
			if isempty(vecSelectCells) || isempty(sRec.sPupil)
				continue;
			end
			
			%calc ccgs
			vecWindow = [-10 10]/1000;
			intResampNum = 50;
			dblUseMaxDur = range(vecWindow);
			dblJitterFactor = 1/max(abs(vecWindow));
			matCCG_P = nan(numel(vecSelectCells),numel(vecSelectCells));
			for intIdxN1=1:numel(vecSelectCells)
				intN1 = vecSelectCells(intIdxN1)
				vecSpikeTimesN1 = sRec.sCluster(intN1).SpikeTimes;
				for intIdxN2=1:numel(vecSelectCells)
					intN2 = vecSelectCells(intIdxN2);
					vecSpikeTimesN2 = sRec.sCluster(intN2).SpikeTimes;

					dblZetaP = getZeta(vecSpikeTimesN1,vecSpikeTimesN2+vecWindow(1),dblUseMaxDur,intResampNum,0,0,[],[],dblJitterFactor);
					matCCG_P(intIdxN1,intIdxN2) = dblZetaP;
				end
			end
			cellCCG{intArea,intSubType}{end+1} = matCCG_P;
		end
	end
end
%save data
save(fullpath(strTargetPath,'CCG_workspace'));
% plot 1
