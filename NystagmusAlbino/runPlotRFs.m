
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
runHeaderNOT;

%% load RF data
strDataSource = 'D:\Data\Results\AlbinoProject\RF_data';
sFiles = dir([strDataSource filesep '*Topo*.mat']);
cellNames = {sExp.Name};

for intFile=1:numel(sFiles)
	strFile = sFiles(intFile).name;
	sLoad=load(fullpath(strDataSource,strFile));
	
	strRec = getFlankedBy(strFile,'RF_','B');
	intRec = find(contains(cellNames,strRec));
	
	%select cells
	sRec = sExp(intRec);
	indUseCells = arrayfun(@(x) x.Violations1ms < 0.25 & abs(x.NonStationarity) < 0.25,sRec.sCluster(:));
			
	%build cell vectors
	cellCellsPerArea = cell(1,numel(cellUseAreas));
	cellAreasPerCluster = {sRec.sCluster.Area};
	for intArea=1:numel(cellUseAreas)
		cellCellsPerArea{intArea} = contains(cellAreasPerCluster,cellUseAreas{intArea},'IgnoreCase',true);
	end
	vecCellsNrPerArea = cellfun(@sum,cellCellsPerArea);
	
	%ctx
	vecSelectCellsCtx = find(indUseCells(:) & cellCellsPerArea{1}(:));
	
	matZetaOn = sLoad.matZetaOn(:,:,vecSelectCellsCtx);
	matZetaOff = sLoad.matZetaOff(:,:,vecSelectCellsCtx);
	
% 	for intCell=1:size(matZetaOn,3)
% 		figure
% 		subplot(2,3,1)
% 		matOnZ = -norminv(matZetaOn(:,:,intCell)/2);
% 		matOnZ = matOnZ-1;
% 		matOnZ(matOnZ<0)=0;
% 		matOffZ = -norminv(matZetaOff(:,:,intCell)/2);
% 		matOffZ = matOffZ-1;
% 		matOffZ(matOffZ<0)=0;
% 		imagesc(matOnZ)
% 		
% 		subplot(2,3,2)
% 		imagesc(matOffZ)
% 		
% 		subplot(2,3,3)
% 		imagesc(matOnZ-matOffZ)
% 		
% 	end
	
	%overall
	figure
		subplot(2,3,1)
		matOnZ = -norminv(matZetaOn(:,:,:)/2);
		matOnZ = matOnZ-1;
		matOnZ(matOnZ<0)=0;
		matOffZ = -norminv(matZetaOff(:,:,:)/2);
		matOffZ = matOffZ-1;
		matOffZ(matOffZ<0)=0;
		
		matOnZ=mean(matOnZ,3);
		matOffZ=mean(matOffZ,3);
		
		imagesc(matOnZ)
		colorbar
		
		subplot(2,3,2)
		imagesc(matOffZ)
		colorbar
		
		subplot(2,3,3)
		imagesc(matOnZ+matOffZ)
		colorbar
		
	%mean counts
	matMeanOn = sLoad.matMeanCountsOn(:,:,vecSelectCellsCtx);
	matMeanOff = sLoad.matMeanCountsOff(:,:,vecSelectCellsCtx);
% 	
% 	for intCell=1:size(matZetaOn,3)
% 		figure
% 		subplot(2,3,1)
% 		matOnZ = -norminv(matZetaOn(:,:,intCell)/2);
% 		matOnZ = matOnZ-1;
% 		matOnZ(matOnZ<0)=0;
% 		matOffZ = -norminv(matZetaOff(:,:,intCell)/2);
% 		matOffZ = matOffZ-1;
% 		matOffZ(matOffZ<0)=0;
% 		imagesc(matOnZ)
% 		
% 		subplot(2,3,2)
% 		imagesc(matOffZ)
% 		
% 		subplot(2,3,3)
% 		imagesc(matOnZ-matOffZ)
% 		
% 	end
% 	
	%overall
		subplot(2,3,4)
		
		matOnZ=mean(matMeanOn,3);
		matOffZ=mean(matMeanOff,3);
		
		imagesc(matOnZ)
		colorbar
		
		subplot(2,3,5)
		imagesc(matOffZ)
		colorbar
		
		subplot(2,3,6)
		imagesc(matOnZ+matOffZ)
		colorbar
		
		title(sprintf('%s',strRec),'interpreter','none');
end
