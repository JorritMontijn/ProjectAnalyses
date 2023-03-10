
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
if ~exist('sExp','var') || ~exist('sAtlas','var')
	clear all;
	runHeaderNOT;
	
	sAtlas = AL_PrepABA(strAllenCCFPath);
	tv = sAtlas.tv;
	av = sAtlas.av;
	st = sAtlas.st;
end

%% load RF data
%set RF parameters
dblDistMouseToScreenCenterX = 8.5;%8.5
dblDistMouseToScreenCenterY = 24;%24
dblDistMouseToScreenCenterZ = -4;%-4
dblScreenWidth_cm = 51;
dblScreenHeight_cm = 29;

%data location
%strDataSource = 'D:\Data\Results\AlbinoProject\RF_data_new';
strDataSource = 'C:\Drive\VisMotorNOT\RF_data_new';
sFiles = dir([strDataSource filesep '*.mat']);
cellNames = {sExp.Name};
cellFilesRF = {sFiles.name};
boolSaveFigs = true;

%pre-allocate
intUseOnOrOff = 3;%1=on, 2=off, other=both
intUseAreaNum = 2;
cellAggMapRF = cell(intUseAreaNum,2);
cellAggCenterRF = cell(intUseAreaNum,2);
cellAggCoords = cell(intUseAreaNum,2);
cellAggArea = cell(intUseAreaNum,2);
cellAggSourceRec = cell(intUseAreaNum,2);
cellAggCoordsAll = cell(1,2);

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
		intRec=vecRunRecs(intRecIdx);
		sRec = sExp(intRec);
		strRec=[sRec.sJson.subject '_' sRec.sJson.date];
		% get probe locations
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
		cellAreasPerCluster = {sRec.sCluster.Area};
		vecSelectCells = contains(cellAreasPerCluster,cellUseAreas{2},'IgnoreCase',true);
		matAllCUPF = matLocPerCluster(:,vecSelectCells);
		cellAggCoordsAll{intSubType} = cat(2,cellAggCoordsAll{intSubType},matAllCUPF);
		
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
		indUseCells = arrayfun(@(x) x.Violations1ms < 0.25 & abs(x.NonStationarity) < 0.25,sRec.sCluster(:));
		
		%build cell vectors
		cellCellsPerArea = cell(1,numel(cellUseAreas));
		for intArea=2%numel(cellUseAreas)
			cellCellsPerArea{intArea} = contains(cellAreasPerCluster,cellUseAreas{intArea},'IgnoreCase',true);
			vecCellsNrPerArea = cellfun(@sum,cellCellsPerArea);
			
			indNanZeta = squeeze(any(any(isnan(sLoad.matZetaOn) |  isnan(sLoad.matZetaOff),1),2));
			matZetaOn = sLoad.matZetaOn;
			matZetaOff = sLoad.matZetaOff;
			%dblCritVal = 0.05/(2*prod(size(sLoad.matZetaOn,[1 2])));
			dblCritVal = 0.01;
			indSignificant = squeeze(any(any(matZetaOn<dblCritVal | matZetaOff<dblCritVal,1),2));
			vecSelectCells = find(indSignificant(:) & ~indNanZeta(:) & indUseCells(:) & cellCellsPerArea{intArea}(:));
			strArea = cellAreaGroupsAbbr{intArea};
			
			if numel(vecSelectCells) < 1
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
			dblS = 0.5;
			matFilt = normpdf(-4:4,0,dblS)'*normpdf(-4:4,0,dblS);
			matFilt = matFilt./sum(matFilt(:));
			matOnZAvg = imfilt(matOnZAvg,matFilt);
			matOffZAvg = imfilt(matOffZAvg,matFilt);
			if intUseOnOrOff == 1
			matOnOffAvg = matOnZAvg;
			elseif intUseOnOrOff == 2
				matOnOffAvg = matOffZAvg;
			else
				matOnOffAvg = (matOnZAvg.*matOffZAvg)./2;;
			end
			
			if size(matOnOffAvg,1)==12
				matOnOffAvg = (matOnOffAvg(1:2:12,1:2:20) + matOnOffAvg(2:2:12,2:2:20))/2;
			end
			vecMaxVals = findmax(matOnOffAvg(:),1);
			[intMaxRow,intMaxCol]=find(matOnOffAvg==vecMaxVals(end));
			
			%mean counts
			matMeanOn = sLoad.matMeanCountsOn(:,:,vecSelectCells);
			matMeanOff = sLoad.matMeanCountsOff(:,:,vecSelectCells);
			%matMeanOn = sLoad.matMeanCountsOn(:,:,vecSelectCells)./sLoad.matSdCountsOn(:,:,vecSelectCells);
			%matMeanOff = sLoad.matMeanCountsOff(:,:,vecSelectCells)./sLoad.matSdCountsOff(:,:,vecSelectCells);
			matOnMAvg=mean(matMeanOn,3);
			matOffMAvg=mean(matMeanOff,3);
			matOnMAvg = imfilt(matOnMAvg,matFilt);
			matOffMAvg = imfilt(matOffMAvg,matFilt);
			matOnOffMAvg = matOffMAvg+matOnMAvg;
			vecMaxVals = findmax(matOnOffMAvg(:),1);
			%[intMaxRow,intMaxCol]=find(matOnOffMAvg==vecMaxVals(end));
			
			if 1
			figure
			subplot(2,3,1)
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
			
			imagesc(matOnMAvg)
			colorbar
			title('Mean On','interpreter','none');
			
			subplot(2,3,5)
			imagesc(matOffMAvg)
			colorbar
			title('Mean Off','interpreter','none');
			
			subplot(2,3,6)
			imagesc(matOnOffMAvg)
			colorbar
			title('Mean On+Off','interpreter','none');
			
			%save plot
			maxfig;drawnow;
			if boolSaveFigs
			export_fig([strTargetPath filesep 'single_recs' filesep sprintf('RF_maps_%s.tif',strRec)]);
			export_fig([strTargetPath filesep 'single_recs' filesep sprintf('RF_maps_%s.pdf',strRec)]);
			end
			end
			
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

%{
%% plot NOT 3D v2
avNot=av==(st.index(contains(st.name,'nucleus of the optic tract','ignorecase',true))+1);
avNot=avNot(1:size(avNot,1)/2,:,:);
%switch dims 1&2 for isosurface/patch
fv = isosurface(avNot);
fv.vertices = fv.vertices(:,[2 1 3]);
fv.faces = fv.faces(:,[2 1 3]);
p = patch(fv);
p.FaceColor = [0.2 0 0];
p.EdgeColor = 'none';
p.FaceAlpha = 0.2;
daspect([1 1 1])
view(3); 
axis tight
camlight 
lighting gouraud
%}
%% plot 3D NOT
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


intPoints = 15;%15
matLinesNot = getTrace3D(avEdge,intPoints,0);
%add offset
matLinesNot2=bsxfun(@plus,matLinesNot,([vecNot1(1) vecNot2(1) vecNot3(1)]));
%transform to microns
matLinesNot2 = sAtlas.VoxelSize.*(matLinesNot2 - sAtlas.Bregma);
matLinesNot2(:,1) = abs(matLinesNot2(:,1));

%indRem = matLinesNot2(:,3) < (size(av,3)/2);
%matLinesNot2(indRem,:) = [];

% plot all NOT cells at their respective locations
figure;maxfig;
hold on
h = plot3(matLinesNot2(:,1), matLinesNot2(:,2), matLinesNot2(:,3), 'Color', [0.8 0 0 0.3]);
h.Annotation.LegendInformation.IconDisplayStyle = 'off';

%get flat
matNot2D = sum(avEdge,3)>0;
matNotPoly=mask2poly(matNot2D','outer','MINDIST')+[vecNot1(1) vecNot2(1)];
%transform to microns
matNotPoly = sAtlas.VoxelSize(1:2).*(matNotPoly - sAtlas.Bregma(1:2));
matNotPoly(:,1) = abs(matNotPoly(:,1));
dblBaseDB = -2900;

% plot cells
cellCol = {};
cellCol{1} = lines(1);
cellCol{2} = [1 0 0];
cellMarker = {'x','o'};
hold on
for intSubType=1%1:2
	
	matCABA = cellAggCoordsAll{intSubType};
	matCABA(1,matCABA(1,:)>vecBregma(1)) = 2*vecBregma(1) - matCABA(1,matCABA(1,:)>vecBregma(1));
	%transform to microns
	matCABA = sAtlas.VoxelSize'.*(matCABA - sAtlas.Bregma');
	matCABA(1,:) = abs(matCABA(1,:));
	%line([matCUPF(1,:)' matCABA(1,:)']',[matCUPF(2,:)' matCABA(2,:)']',[matCUPF(3,:)' matCABA(3,:)']')
	
	
	%h= scatter3(matCABA(1,:),matCABA(2,:),matCABA(3,:),[],vecRLR,'marker',cellMarker{1});
	%cellText = cellAggSelfArea{intAreaType,intSubType};
	%text(matCABA(1,:),matCABA(2,:),matCABA(3,:),cellText);
	
	h2= scatter3(matCABA(1,:),matCABA(2,:),matCABA(3,:),[],cellCol{intSubType},'filled');
	h3= scatter3(matCABA(1,:),matCABA(2,:),dblBaseDB*ones(size(matCABA(3,:))),[],cellCol{intSubType},'filled');
	%cellText = cellAggArea{intAreaType,intSubType};
	%text(matCUPF(1,:),matCUPF(2,:),matCUPF(3,:),cellText);
end
%plot flat projection at DV=dblBaseDB
h = plot3(matNotPoly([1:(end-1) 1],1)', matNotPoly([1:(end-1) 1],2)', dblBaseDB*ones(size(matNotPoly([1:(end-1) 1],1)')),'Color', [0.8 0 0 1]);
	

hold off
xlabel('ML');
ylabel('AP');
zlabel('DV');
axis equal;
axis xy;
fixfig;
grid off;
%campos([857 -1500 750])%1.4941    0.8571    0.7502%[1000 -50   100]
%campos([848.0622 -182.0584  825.3987])
campos([-4000    -9500   350]);
title(sprintf('Recording locations in NOT'));

%save plot
drawnow;
if boolSaveFigs
export_fig([strTargetPath filesep sprintf('RecLocsNOT3D.tif')]);
export_fig([strTargetPath filesep sprintf('RecLocsNOT3D.pdf')]);
end

%% find NOT
%find bregma
vecBregma = sAtlas.Bregma;

avNot=av==(st.index(contains(st.name,'nucleus of the optic tract','ignorecase',true))+1);
vecRangeNot1 = find(sum(sum(avNot,2),3));
vecRangeNot2 = find(sum(sum(avNot,1),3));
vecRangeNot3 = find(sum(sum(avNot,1),2));
vecNot1 = (vecRangeNot1(1)-1):(vecRangeNot1(end)+1);
vecNot1(vecNot1<vecBregma(1))=[];
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

%transform atlas to microns relative to bregma
intAreaType = 2;
matCenterRF = cellAggCenterRF{intAreaType,1};
%matCenterRF = cat(2,matCenterRF,cellAggCenterRF{intAreaType,2});
matCenterAlbRF = cellAggCenterRF{intAreaType,2};
matCoords = cellAggCoords{intAreaType,1};
%matCoords = cat(2,matCoords,cellAggCoords{intAreaType,2});
matCoords(1,matCoords(1,:)<vecBregma(1)) = 2*vecBregma(1) - matCoords(1,matCoords(1,:)<vecBregma(1));
matCoordsMu = sAtlas.VoxelSize'.*(matCoords - sAtlas.Bregma');
matCoordsMu(1,:) = abs(matCoordsMu(1,:));
matNotPolyMu= sAtlas.VoxelSize(1:2).*(matNotPoly - sAtlas.Bregma(1:2));
matNotPolyMu(:,1) = abs(matNotPolyMu(:,1));

%% test retinotopy with 90-deg constraint between x and y; readout is mean corr of x & y
%find most predictive angle
fMinFunc = @(x) -getMeanRFcorrWithAngleConstraint(x,matCoordsMu(1,:)',matCoordsMu(2,:)',matCenterRF(1,:)',matCenterRF(2,:)');
[dblRealOptAngle,fval,exitflag,output] = fminbnd(fMinFunc,-2*pi,2*pi);
[dblRealOptCorr,dblCorr1,vecProjectedLocation1,matProjectedPoints1,dblCorr2,vecProjectedLocation2,matProjectedPoints2] = getMeanRFcorrWithAngleConstraint(dblRealOptAngle,matCoordsMu(1,:)',matCoordsMu(2,:)',matCenterRF(1,:)',matCenterRF(2,:)');

%compare with random (shuffled) rf locations
intRunNum = 1000;
intP = numel(matCenterRF(1,:)');
intRandIters = intRunNum;%min(intRunNum,factorial(intP));
vecRandCorr = nan(1,intRandIters);
vecRandAngle = nan(1,intRandIters);
hTic=tic;
for intIter=1:intRandIters
	%randomize RF location
	vecRandRFs1 = matCenterRF(1,randperm(intP));
	vecRandRFs2 = matCenterRF(2,randperm(intP));
	
	%find most predictive angle
	fMinFunc = @(x) -getMeanRFcorrWithAngleConstraint(x,matCoordsMu(1,:)',matCoordsMu(2,:)',vecRandRFs1',vecRandRFs2');
	[dblOptAngle,fval,exitflag,output] = fminbnd(fMinFunc,-2*pi,2*pi);
	dblOptCorr = getMeanRFcorrWithAngleConstraint(dblOptAngle,matCoordsMu(1,:)',matCoordsMu(2,:)',vecRandRFs1',vecRandRFs2');
	
	%save
	vecRandCorr(intIter) = dblOptCorr;
	vecRandAngle(intIter) = dblOptAngle;
	
	%msg
	if toc(hTic) > 5
		fprintf('%d/%d...\n',intIter,intRandIters);
		hTic=tic;
	end
end

%plot rand
dblCorrP = sum(vecRandCorr>dblRealOptCorr)./intRandIters;

%% plot
%plot
figure;maxfig;

%transform RF block location to retinal degrees
dblSubjectPosX_cm = dblDistMouseToScreenCenterX; %left/right component to center of screen; negative is left of screen center
dblSubjectPosY_cm = dblDistMouseToScreenCenterZ; %up/down component to center of screen; negative is lower than screen center
dblScreenDistance_cm = dblDistMouseToScreenCenterY; %forward component to center of screen

%assume screen is a spherical patch
dblScreenCenterX = -dblSubjectPosX_cm;
dblScreenCenterZ = -dblSubjectPosY_cm;
dblScreenCenterD = dblScreenDistance_cm;
[azimuth,elevation,distance] = cart2sph(dblScreenCenterX,dblScreenCenterD,dblScreenCenterZ);
azimuth = azimuth - 0.5*pi; %make forward azimuth=0
dblScreenWidth_rad = pi - 2*atan2(distance,dblScreenWidth_cm/2);
dblScreenHeight_rad = pi - 2*atan2(distance,dblScreenHeight_cm/2);
az_Left = rad2deg(azimuth+dblScreenWidth_rad/2); %positive azimuth is left
az_Right = rad2deg(azimuth-dblScreenWidth_rad/2); %negative azimuth is right
el_Up = rad2deg(elevation+dblScreenHeight_rad/2); %positive elevation is up
el_Down = rad2deg(elevation-dblScreenHeight_rad/2); %negative elevation is down
az_Tot = az_Left-az_Right;
el_Tot = el_Up-el_Down;
az_Center = rad2deg(azimuth);
el_Center = rad2deg(elevation);

%get RF locations in degrees
[intBlocksVert,intBlocksHorz] = size(matOnOffAvg);
dblHorzBlockCenter = (intBlocksHorz+1)/2;
dblVertBlockCenter = (intBlocksVert+1)/2;
dblHorzDegsPerBlock = az_Tot/intBlocksHorz;
dblVertDegsPerBlock = el_Tot/intBlocksVert;
vecHorzRF_deg = az_Center+dblHorzDegsPerBlock*(matCenterRF(1,:)'-dblHorzBlockCenter);
vecVertRF_deg = el_Center+dblVertDegsPerBlock*(dblVertBlockCenter-matCenterRF(2,:)'); %RF blocks are from top to bottom, not bottom to top

if isempty(matCenterAlbRF)
	vecHorzAlbRF_deg = [];
	vecVertAlbRF_deg = [];
	matCoordsAlb = [];
	matCoordsMuAlb = [];
else
	%prep albinos
	vecHorzAlbRF_deg = az_Center+dblHorzDegsPerBlock*(matCenterAlbRF(1,:)'-dblHorzBlockCenter);
	vecVertAlbRF_deg = el_Center+dblVertDegsPerBlock*(dblVertBlockCenter-matCenterAlbRF(2,:)'); %RF blocks are from top to bottom, not bottom to top
matCoordsAlb = cellAggCoords{2,2};
matCoordsAlb(1,matCoordsAlb(1,:)<vecBregma(1)) = 2*vecBregma(1) - matCoordsAlb(1,matCoordsAlb(1,:)<vecBregma(1));
matCoordsMuAlb = sAtlas.VoxelSize'.*(matCoordsAlb - sAtlas.Bregma');
matCoordsMuAlb(1,:) = abs(matCoordsMuAlb(1,:));
end


%get mid point
dblCenterX = mean(matCoordsMu(1,:));
dblCenterY = mean(matCoordsMu(2,:));
dblArrowLength = 10; %degrees

%fit RFs to estimate retinotopic map
%center x/y
vecXY0=mean(matCoordsMu(1:2,:),2);
matXY0=matCoordsMu(1:2,:)-vecXY0;
[dblRealOptCorr,dblCorr1,vecProjectedLocation1,matProjectedPoints1,dblCorr2,vecProjectedLocation2,matProjectedPoints2] = ...
	getMeanRFcorrWithAngleConstraint(dblRealOptAngle,matCoordsMu(1,:)',matCoordsMu(2,:)',vecHorzRF_deg,vecVertRF_deg);

if intUseOnOrOff == 1
	strOnOff = 'On';
elseif intUseOnOrOff == 2
	strOnOff = 'Off';
else
	strOnOff = 'OnOff';
end

%prep albinos
for intXY=1:2
	%rotate reference vector
	if intXY == 1
		vecRF_deg = vecHorzRF_deg;
		vecAlbRF_deg = vecHorzAlbRF_deg;
		dblRot = deg2rad(0);
		strAzEl = 'azimuth';
		strHorzVert = 'Horizontal';
	else
		vecRF_deg = vecVertRF_deg;
		vecAlbRF_deg = vecVertAlbRF_deg;
		dblRot = deg2rad(90);
		strAzEl = 'elevation';
		strHorzVert = 'Vertical';
	end
	matRot = [cos(dblRealOptAngle+dblRot) sin(dblRealOptAngle+dblRot);...
		-sin(dblRealOptAngle+dblRot) cos(dblRealOptAngle+dblRot)];
	vecRefVector=[1;0];
	vecRotRef = matRot * vecRefVector;
	
	%calc corr
	[vecProjectedLocation,matProjectedPoints] = getProjOnLine(matXY0,vecRotRef);
	dblCorr1 = corr(vecProjectedLocation,vecRF_deg);
	
	%lin fit
	p1 = polyfit(vecProjectedLocation,vecRF_deg,1);
	vecRF0degs = vecRotRef*(-p1(2)/p1(1))+vecXY0;
	vecRF10degs = vecRotRef*((10-p1(2))/p1(1))+vecXY0;
	
	%get min/max extent of RF degs in NOT
	[vecProjectedLocationNot,matProjectedPoints] = getProjOnLine(matNotPolyMu'-vecXY0,vecRotRef);
	dblMinRF_degs = min(vecProjectedLocationNot)*p1(1)+p1(2);
	dblMaxRF_degs = max(vecProjectedLocationNot)*p1(1)+p1(2);
	if dblMaxRF_degs < dblMinRF_degs
		[dblMaxRF_degs,dblMinRF_degs] = swap(dblMinRF_degs,dblMaxRF_degs);
	end
	vecMinRF_loc = vecRotRef*((dblMinRF_degs-p1(2))/p1(1))+vecXY0;
	vecMaxRF_loc = vecRotRef*((dblMaxRF_degs-p1(2))/p1(1))+vecXY0;
	
	%RF x
	subplot(2,3,intXY)
	h = plot(matNotPolyMu([1:(end-1) 1],1)', matNotPolyMu([1:(end-1) 1],2)', 'Color', [1 0 0 0.3]);
	axis equal;
	hold on
	
	%plot bl6
	h2= scatter(matCoordsMu(1,:)',matCoordsMu(2,:)',[],vecRF_deg,'filled');
	h2.SizeData=100;
	h2.MarkerEdgeColor = 'k';
	
	%plot dba
	if ~isempty(matCoordsMuAlb)
		h2= scatter(matCoordsMuAlb(1,:)',matCoordsMuAlb(2,:)',[],vecAlbRF_deg,'square','filled');
		h2.SizeData=100;
		h2.MarkerEdgeColor = 'k';
	end
	
	hold on
	plot([0 vecRotRef(1)*100]+vecXY0(1),[0 vecRotRef(2)*100]+vecXY0(2))
	scatter(vecRF0degs(1),vecRF0degs(2),'kx');
	plot([vecRF0degs(1) vecRF10degs(1)],[vecRF0degs(2) vecRF10degs(2)],'k');
	hold off
	xlabel('Anatomical ML location (microns)');
	ylabel('Anatomical AP location (microns)');
	hC=colorbar;
	ylabel(hC,sprintf('RF %s (degs)',strAzEl));
	title(sprintf('%s RF center',strHorzVert));
	fixfig;grid off;
	
	%2nd plot
	subplot(2,3,intXY+3)
	h = plot(matNotPolyMu([1:(end-1) 1],1)', matNotPolyMu([1:(end-1) 1],2)', 'Color', [1 0 0 0.3]);
	axis equal;
	hold on
	set(gca,'clim',[dblMinRF_degs dblMaxRF_degs]);
	plot([vecMinRF_loc(1) vecMaxRF_loc(1)],[vecMinRF_loc(2) vecMaxRF_loc(2)],'b');
	
	%plot bl6
	h2= scatter(matCoordsMu(1,:)',matCoordsMu(2,:)',[],vecRF_deg,'filled');
	h2.SizeData=100;
	h2.MarkerEdgeColor = 'k';
	
	for intRec=1:numel(vecRF_deg)
		intExp = cellAggSourceRec{intArea,1}(intRec);
		strRec=[sExp(intExp).sJson.subject '_' sExp(intExp).sJson.date];
		text(matCoordsMu(1,intRec)',matCoordsMu(2,intRec)',strRec,'interpreter','none');
	end
	
	%plot dba
	if ~isempty(matCoordsMuAlb)
		h2= scatter(matCoordsMuAlb(1,:)',matCoordsMuAlb(2,:)',[],vecAlbRF_deg,'square','filled');
		h2.SizeData=100;
		h2.MarkerEdgeColor = 'k';
	end
	
	for intRec=1:numel(vecAlbRF_deg)
		intExp = cellAggSourceRec{intArea,2}(intRec);
		strRec=[sExp(intExp).sJson.subject '_' sExp(intExp).sJson.date];
		text(matCoordsMuAlb(1,intRec)',matCoordsMuAlb(2,intRec)',strRec,'interpreter','none');
	end
	
	scatter(vecRF0degs(1),vecRF0degs(2),'kx');
	xlabel('Anatomical ML location (microns)');
	ylabel('Anatomical AP location (microns)');
	hC=colorbar;
	ylabel(hC,sprintf('RF %s (degs)',strAzEl));
	title(sprintf('%s %s RF center, min = %.1f degs, max=%.1f degs',strHorzVert,strOnOff,dblMinRF_degs,dblMaxRF_degs));
	fixfig;grid off;
	
end

%significance
subplot(2,3,3)

vecBinsE = 0:0.02:1;
vecBinsC = vecBinsE(2:end) - median(diff(vecBinsE))/2;
vecCounts = histcounts(vecRandCorr,vecBinsE);
plot(vecBinsC,vecCounts);
hold on
plot([dblRealOptCorr dblRealOptCorr],[0 0.9*max(get(gca,'ylim'))],'r--')
%text(dblRealOptCorr,0.85*max(get(gca,'ylim')),'Real','color',[1 0 0],'VerticalAlignment','bottom','HorizontalAlignment','center','fontsize',14);
hold off;
xlabel('r(anatomical location,RF center)');
ylabel('Number of shuffled RF centers (count)')
legend({'Shuffled','Real'});
title(sprintf('Retinotopy permutation test, p=%.4f',dblCorrP));
fixfig;grid off;

%save plot
drawnow;
if boolSaveFigs
export_fig([strTargetPath filesep sprintf('%sRetinotopyNOT.tif',strOnOff)]);
export_fig([strTargetPath filesep sprintf('%sRetinotopyNOT.pdf',strOnOff)]);
end

%% for single direction
%{
%find most predictive angle
fMinFunc = @(x) -getCorrAtAngle(x,matCoords(1,:)',matCoords(2,:)',matCenterRF(intXY,:)');
[dblRealOptAngle,fval,exitflag,output] = fminbnd(fMinFunc,-2*pi,2*pi);
dblRealOptCorr = getCorrAtAngle(dblRealOptAngle,matCoords(1,:)',matCoords(2,:)',matCenterRF(intXY,:)');

%compare with random (shuffled) rf locations
intP = numel(matCenterRF(intXY,:)');
intRandIters = min(10000,factorial(intP));
vecRandCorr = nan(1,intRandIters);
vecRandAngle = nan(1,intRandIters);
hTic=tic;
for intIter=1:intRandIters
	%randomize RF location
	%vecRandRFs = matCenterRF(intXY,randperm(intP));
	vecRandRFs = randi(intRange,[1 intP]);
	
	%find most predictive angle
	fMinFunc = @(x) -getCorrAtAngle(x,matCoords(1,:)',matCoords(2,:)',vecRandRFs');
	[dblOptAngle,fval,exitflag,output] = fminbnd(fMinFunc,-2*pi,2*pi);
	dblOptCorr = getCorrAtAngle(dblOptAngle,matCoords(1,:)',matCoords(2,:)',vecRandRFs');
	
	%save
	vecRandCorr(intIter) = dblOptCorr;
	vecRandAngle(intIter) = dblOptAngle;
	
	%msg
	if toc(hTic) > 5
		fprintf('%d/%d...\n',intIter,intRandIters);
		hTic=tic;
	end
end

%plot rand
dblCorrP = sum(vecRandCorr>dblRealOptCorr)./intRandIters
%}