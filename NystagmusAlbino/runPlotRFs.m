
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

sAtlas = AL_PrepABA(strAllenCCFPath);
tv = sAtlas.tv;
av = sAtlas.av;
st = sAtlas.st;

%% load RF data
strDataSource = 'D:\Data\Results\AlbinoProject\RF_data';
sFiles = dir([strDataSource filesep '*.mat']);
cellNames = {sExp.Name};
cellFilesRF = {sFiles.name};

%pre-allocate
intUseAreaNum = 2;
cellAggMapRF = cell(intUseAreaNum,2);
cellAggCenterRF = cell(intUseAreaNum,2);
cellAggCoords = cell(intUseAreaNum,2);
cellAggArea = cell(intUseAreaNum,2);
cellAggSourceRec = cell(intUseAreaNum,2);

% run
cellNameAP = arrayfun(@(x) x.sJson.file_preproAP,sExp,'uniformoutput',false);
cellExperiment = arrayfun(@(x) x.sJson.experiment,sExp,'uniformoutput',false);
cellRemove = {};%{'RecMA5_2021-03-01R01_g0_t0'};
indRemRecs = contains(cellExperiment,cellRemove);
cellSubjectType = arrayfun(@(x) x.sJson.subjecttype,sExp,'uniformoutput',false);
for intSubType=[2 1]
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
		intRec=vecRunRecs(intRecIdx);
		sRec = sExp(intRec);
		strRec=[sRec.sJson.subject '_' sRec.sJson.date];
		
		%load RF file
		intFileRF = find(contains(cellFilesRF,strRec));
		if isempty(intFileRF)
			fprintf('No RF file for %s: skipping...\n',strRec);
			continue;
		else
			strFile = sFiles(intFileRF).name;
			fprintf('Running %s - %s\n',strRec,strFile);
		end
		sLoad=load(fullpath(strDataSource,strFile));
		
		%select cells
		sRec = sExp(intRec);
		indUseCells = true;%arrayfun(@(x) x.Violations1ms < 0.25 & abs(x.NonStationarity) < 0.25,sRec.sCluster(:));
		
		%build cell vectors
		cellCellsPerArea = cell(1,numel(cellUseAreas));
		cellAreasPerCluster = {sRec.sCluster.Area};
		for intArea=2%numel(cellUseAreas)
			cellCellsPerArea{intArea} = contains(cellAreasPerCluster,cellUseAreas{intArea},'IgnoreCase',true);
			vecCellsNrPerArea = cellfun(@sum,cellCellsPerArea);
			
			indNanZeta = squeeze(any(any(isnan(sLoad.matZetaOn) |  isnan(sLoad.matZetaOff),1),2));
			matZetaOn = sLoad.matZetaOn;
			matZetaOff = sLoad.matZetaOff;
			dblCritVal = 0.01;%/(2*prod(size(sLoad.matZetaOn,[1 2])));
			indSignificant = squeeze(any(any(matZetaOn<dblCritVal | matZetaOff<dblCritVal,1),2));
			vecSelectCells = find(indSignificant(:) & ~indNanZeta(:) & indUseCells(:) & cellCellsPerArea{intArea}(:));
			strArea = cellAreaGroupsAbbr{intArea};
			
			if isempty(vecSelectCells)
				fprintf('No (significant) cells in %s for %s: skipping...\n',strArea,strRec);
				continue;
			end
			
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
			dblLowerBound = 0;
			matOnZ = -norminv(matZetaOn(:,:,vecSelectCells)/2);
			%matOnZ = matZetaOn(:,:,vecSelectCells)<dblCritVal;
			matOnZ = matOnZ-dblLowerBound;
			matOnZ(matOnZ<0)=0;
			matOffZ = -norminv(matZetaOff(:,:,vecSelectCells)/2);
			%matOffZ = matZetaOff(:,:,vecSelectCells)<dblCritVal;
			matOffZ = matOffZ-dblLowerBound;
			matOffZ(matOffZ<0)=0;
			
			matOnZAvg=mean(matOnZ+dblLowerBound,3);
			matOffZAvg=mean(matOffZ+dblLowerBound,3);
			
			%filter
			matFilt = normpdf(-2:2,0,0.5)'*normpdf(-2:2,0,0.5);
			matFilt = matFilt./sum(matFilt(:));
			matOnZAvg = imfilt(matOnZAvg,matFilt);
			matOffZAvg = imfilt(matOffZAvg,matFilt);
			matOnOffAvg = (matOnZAvg+matOffZAvg)./2;
			if size(matOnOffAvg,1)==12
				matOnOffAvg = (matOnOffAvg(1:2:12,1:2:20) + matOnOffAvg(2:2:12,2:2:20))/2;
			end
			[intMaxRow,intMaxCol]=find(matOnOffAvg==max(matOnOffAvg(:)));
			
			imagesc(matOnZAvg)
			colorbar
			title(sprintf('%s %s\nZETA On',strRec,strArea),'interpreter','none');
			
			subplot(2,3,2)
			imagesc(matOffZAvg)
			colorbar
			title('ZETA Off','interpreter','none');
			
			subplot(2,3,3)
			imagesc(matOnOffAvg)
			hold on
			scatter(intMaxCol,intMaxRow,'rx');
			hold off
			colorbar
			title('ZETA On+Off','interpreter','none');
			
			%mean counts
			matMeanOn = sLoad.matMeanCountsOn(:,:,vecSelectCells);
			matMeanOff = sLoad.matMeanCountsOff(:,:,vecSelectCells);
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
			
			%overall
			subplot(2,3,4)
			
			matOnMAvg=mean(matMeanOn,3);
			matOffMAvg=mean(matMeanOff,3);
			matOnMAvg = imfilt(matOnMAvg,matFilt);
			matOffMAvg = imfilt(matOffMAvg,matFilt);
			
			imagesc(matOnMAvg)
			colorbar
			title('Mean On','interpreter','none');
			
			subplot(2,3,5)
			imagesc(matOffMAvg)
			colorbar
			title('Mean Off','interpreter','none');
			
			subplot(2,3,6)
			imagesc(matOnMAvg+matOffMAvg)
			colorbar
			title('Mean On+Off','interpreter','none');
			
			%% get probe locations
			%get locations along probe
			sProbeCoords = sRec.sSources.sProbeCoords;
			[probe_area_ids,probe_area_boundaries,probe_area_centers,matLocCh] = PH_GetProbeAreas(sProbeCoords.sProbeAdjusted.probe_vector_cart,sAtlas.av);
			
			%get clusters
			dblProbeLengthProbe = sRec.sSources.sProbeCoords.ProbeLengthMicrons;
			dblProbeLengthSph = sRec.sSources.sProbeCoords.sProbeAdjusted.probe_vector_sph(6)*mean(sRec.sSources.sProbeCoords.VoxelSize);
			
			vecDepth = cell2vec({sRec.sCluster.Depth});
			vecDepthOnProbe = (vecDepth/dblProbeLengthSph)*dblProbeLengthProbe;
			
			[vecClustAreaId,cellClustAreaLabel,cellClustAreaFull,vecVoxelDepth] = PF_GetAreaPerCluster(sProbeCoords,vecDepthOnProbe);
			vecUsedDepth = (min(max(round(vecVoxelDepth),1),floor(sProbeCoords.sProbeAdjusted.probe_vector_sph(end)))/floor(sProbeCoords.sProbeAdjusted.probe_vector_sph(end)))*sProbeCoords.ProbeLengthMicrons;
			vecUsedLocsPerCluster = min(max(round(vecVoxelDepth),1),floor(sProbeCoords.sProbeAdjusted.probe_vector_sph(end)));
			matLocPerCluster = matLocCh(:,vecUsedLocsPerCluster);
			
			%get source rec
			vecSourceRec = ones(1,numel(vecSelectCells))*intRec;
			
			%% locations per cluster
			matCUPF = matLocPerCluster(:,vecSelectCells);
			cellAggMapRF{intArea,intSubType} = cat(3,cellAggMapRF{intArea,intSubType},matOnOffAvg);
			cellAggCenterRF{intArea,intSubType} = cat(2,cellAggCenterRF{intArea,intSubType},[intMaxCol;intMaxRow]);%x,y
			cellAggCoords{intArea,intSubType} = cat(2,cellAggCoords{intArea,intSubType},mean(matCUPF,2));
			cellAggSourceRec{intArea,intSubType} = cat(2,cellAggSourceRec{intArea,intSubType},mean(vecSourceRec,2));
			cellAggArea{intArea,intSubType} = cat(2,cellAggArea{intArea,intSubType},cellClustAreaFull(vecSelectCells)');
		end
	end
end

%% find NOT
%find bregma
vecBregma = sAtlas.Bregma;

avNot=av==(st.index(contains(st.name,'nucleus of the optic tract','ignorecase',true))+1);
vecRangeNot1 = find(sum(sum(avNot,2),3));
vecRangeNot2 = find(sum(sum(avNot,1),3));
vecRangeNot3 = find(sum(sum(avNot,1),2));
vecNot1 = (vecRangeNot1(1)-1):(vecRangeNot1(end)+1);
vecNot1(vecNot1>vecBregma(1))=[];
vecNot2 = (vecRangeNot2(1)-1):(vecRangeNot2(end)+1);
vecNot3 = (vecRangeNot3(1)-1):(vecRangeNot3(end)+1);
avNot = avNot(vecNot1,vecNot2,vecNot3);
avCenter = imfill(avNot,'holes');
SE = strel('sphere',1);
avErode = imerode(avNot,SE);
avEdge = avCenter - avErode;

intPoints = 15;

matNot2D = sum(avEdge,3)>0;
matNotPoly=mask2poly(matNot2D','outer','MINDIST')+[vecNot1(1) vecNot2(1)];

%% plot
figure;maxfig;
for intXY=1:2
	subplot(2,2,intXY)
	h = plot(matNotPoly([1:(end-1) 1],1)', matNotPoly([1:(end-1) 1],2)', 'Color', [1 0 0 0.3]);
	axis equal;
	%h.Annotation.LegendInformation.IconDisplayStyle = 'off';
	hold on
	cellMarker = {'x','*'};
	intAreaType = 2;
	for intSubType=2%1%:2
		
		matCenterRF = cellAggCenterRF{intAreaType,intSubType};
		matCoords = cellAggCoords{intAreaType,intSubType};
		matCoords(1,matCoords(1,:)>vecBregma(1)) = 2*vecBregma(1) - matCoords(1,matCoords(1,:)>vecBregma(1));
		
		%line([matCUPF(1,:)' matCABA(1,:)']',[matCUPF(2,:)' matCABA(2,:)']',[matCUPF(3,:)' matCABA(3,:)']')
		
		
		%h= scatter3(matCABA(1,:),matCABA(2,:),matCABA(3,:),[],vecRLR,'marker',cellMarker{1});
		%cellText = cellAggSelfArea{intAreaType,intSubType};
		%text(matCABA(1,:),matCABA(2,:),matCABA(3,:),cellText);
		
		h2= scatter(matCoords(1,:)',matCoords(2,:)',[],matCenterRF(intXY,:)','marker',cellMarker{intSubType});
		h2.SizeData=100;
		%cellText = cellAggArea{intAreaType,intSubType};
		%text(matCUPF(1,:),matCUPF(2,:),matCUPF(3,:),cellText);
		
		%find most predictive angle
		fMinFunc = @(x) -getCorrAtAngle(x,matCoords(1,:)',matCoords(2,:)',matCenterRF(intXY,:)');
		[dblOptAngle,fval,exitflag,output] = fminbnd(fMinFunc,-2*pi,2*pi);
		dblOptCorr = getCorrAtAngle(dblOptAngle,matCoords(1,:)',matCoords(2,:)',matCenterRF(intXY,:)');
		
		%compare with random (shuffled) rf locations
		
		%runRandIters = min(100,factorial(n))
		
	end
	fixfig;grid off;
	if intXY == 1
		title(sprintf('Horizontal RF location'));
	elseif intXY == 2
		title(sprintf('Vertical RF location'));
		
	end
end