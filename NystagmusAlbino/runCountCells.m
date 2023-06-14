%{
ori vs dir tuning
%}
%% load data and define groups
%strDataPath
%cellUseForEyeTracking
%strTargetPath
%cellUseAreas{1} = {'Primary visual','Posteromedial visual'};
%cellUseAreas{2} = {'nucleus of the optic tract'};
%cellUseAreas{3} = {'superior colliculus'};
%cellAreaGroups = {'Vis. ctx','NOT','Hippocampus'};
%cellAreaGroupsAbbr = {'Ctx','NOT','Hip'};
%cellSubjectGroups = {'BL6','DBA'};

clear all;
runHeaderNOT;

%% pre-allocate
vecTotCells_AlbCtx = [];
vecKeepCells_AlbCtx = [];


vecTotCells_Bl6Ctx = [];
vecKeepCells_Bl6Ctx = [];

vecTotCells_Bl6NOT = [];
vecKeepCells_Bl6NOT = [];


vecTotCells_AlbNOT = [];
vecKeepCells_AlbNOT = [];
					
%% run
cellNameAP = arrayfun(@(x) x.sJson.file_preproAP,sExp,'uniformoutput',false);
cellExperiment = arrayfun(@(x) x.sJson.experiment,sExp,'uniformoutput',false);
cellRemove = {};%{'RecMA5_2021-03-01R01_g0_t0'};
indRemRecs = contains(cellExperiment,cellRemove);
indRemRecs2 = ~contains(cellNameAP,cellUseForEyeTracking);
cellSubjectType = arrayfun(@(x) x.sJson.subjecttype,sExp,'uniformoutput',false);
dblAverageMouseHeadTiltInSetup = -15;
for intSubType=1:2
	intPopCounter = 0;
	if intSubType == 1
		intBestRec = 17;
		strSubjectType = 'BL6';
		dblOffsetT=0;
	elseif intSubType == 2
		intBestRec = 5;
		strSubjectType = 'DBA';
		dblOffsetT=0;
	end
	indUseRecs = contains(cellSubjectType,strSubjectType);
	vecRunRecs = find(indUseRecs & ~(indRemRecs));% | indRemRecs2));
	%vecRunRecs = intBestRec;
	for intRecIdx=1:numel(vecRunRecs)
		intRec=vecRunRecs(intRecIdx);
		sRec = sExp(intRec);
		
		%include?
		indUseCells = arrayfun(@(x) x.Violations1ms < 0.25 & abs(x.NonStationarity) < 0.25,sRec.sCluster(:));
		
		%% split cells into areas
		%build cell vectors
		cellCellsPerArea = cell(1,numel(cellUseAreas));
		cellAreasPerCluster = {sRec.sCluster.Area};
		for intArea=1:numel(cellUseAreas)
			cellCellsPerArea{intArea} = contains(cellAreasPerCluster,cellUseAreas{intArea},'IgnoreCase',true);
		end
		vecCellsNrPerArea = cellfun(@sum,cellCellsPerArea);
		
		for intArea = 1:2%numel(cellUseAreas)
			%% select cells
			strAreaGroup =  cellAreaGroupsAbbr{intArea};
			vecSelectCells = find(indUseCells(:) & cellCellsPerArea{intArea}(:));
			intTotCells = sum(cellCellsPerArea{intArea}(:));
			intKeepCells = numel(vecSelectCells);
			
			%% save data
			if intArea == 1
				if strcmpi(strSubjectType,'DBA')
					vecTotCells_AlbCtx(end+1) = intTotCells;
					vecKeepCells_AlbCtx(end+1) = intKeepCells;
				elseif strcmpi(strSubjectType,'Bl6')
					vecTotCells_Bl6Ctx(end+1) = intTotCells;
					vecKeepCells_Bl6Ctx(end+1) = intKeepCells;
				end
			elseif intArea == 2
				if strcmpi(strSubjectType,'DBA')
					vecTotCells_AlbNOT(end+1) = intTotCells;
					vecKeepCells_AlbNOT(end+1) = intKeepCells;
				elseif strcmpi(strSubjectType,'Bl6')
					vecTotCells_Bl6NOT(end+1) = intTotCells;
					vecKeepCells_Bl6NOT(end+1) = intKeepCells;
				end
			end
		end
	end
end
fprintf([...
	'Cells,   All / Incl / Frac indcluded\n'...
	'CtxBl6: %04d / %04d / %.3f\n'...
	'NotBl6: %04d / %04d / %.3f\n'...
	'CtxAlb: %04d / %04d / %.3f\n'...
	'NotAlb: %04d / %04d / %.3f\n'],...
	sum(vecTotCells_Bl6Ctx),sum(vecKeepCells_Bl6Ctx),sum(vecKeepCells_Bl6Ctx)/sum(vecTotCells_Bl6Ctx),...
	sum(vecTotCells_Bl6NOT),sum(vecKeepCells_Bl6NOT),sum(vecKeepCells_Bl6NOT)/sum(vecTotCells_Bl6NOT),...
	sum(vecTotCells_AlbCtx),sum(vecKeepCells_AlbCtx),sum(vecKeepCells_AlbCtx)/sum(vecTotCells_AlbCtx),...
	sum(vecTotCells_AlbNOT),sum(vecKeepCells_AlbNOT),sum(vecKeepCells_AlbNOT)/sum(vecTotCells_AlbNOT)...
	);
