%run data with NOT
%S1L3
%S2L1
%S2L4 => eye tracking wrong?
%S3L1 => probably not NOT
%S3L3
clearvars;
%cellRun = {'S1L3','S2L1',};
%cellRun(end+1:end+2) = {'S3L3','S2L4','S3L1'};
cellRun = {'S1L1','S1L2','S1L3',...
	'S2L1','S2L2','S2L3','S2L4','S2L5','S2L6',...
	'S3L1'};%,'S3L2','S3L3'};
boolAllPlots = false;
boolPlotChange = true;

%get data
strSearchFormat = '\d{4}[-_]?\d{2}[-_]?\d{2}';
strOutputPath = 'F:\Data\Results\EyeTracking\Decoding\';
strPupilSourcePath = 'F:\Drive\EyeTrackingProcessed\';
strDataSourcePath = 'F:\Data\Processed\Neuropixels\';
clear sFiles;
cellPupilFiles = [];
for intFile=1:numel(cellRun)
	sFile = dir([strDataSourcePath '*' cellRun{intFile} '*.mat']);
	sFiles(intFile) = sFile(1);
	%get date
	[intB,intE] = regexp(sFiles(intFile).name,strSearchFormat);
	strDate = sFiles(intFile).name(intB:intE);
	%find matching eye-tracking data
	sPupilFiles1 = dir([strPupilSourcePath '*' strDate '*Processed*.mat']);
	sPupilFiles2 = dir([strPupilSourcePath '*Processed*' strDate '*.mat']);
	vecDist1 =strdist(sFiles(intFile).name, {sPupilFiles1.name});
	vecDist2 =strdist(sFiles(intFile).name, {sPupilFiles2.name});
	if min(vecDist1) < min(vecDist2)
		[v,i]=min(vecDist1);
		strPupilFile = sPupilFiles1(i).name;
	else
		[v,i]=min(vecDist2);
		strPupilFile = sPupilFiles2(i).name;
	end
	
	cellPupilFiles{intFile} = strPupilFile;
end
cellFiles = {sFiles(:).name}';
cellPupilFiles = cellPupilFiles';
%close all;


%% pre-allocate
cellRunAreas = {...
	'lateral geniculate',...Area 1
	'Primary visual',...Area
	'Lateral posterior nucleus',...Area
	'Anterior pretectal nucleus',...Area
	'Nucleus of the optic tract',...Area
	'Superior colliculus',...Area
	'Anteromedial visual',...Area
	'posteromedial visual',...Area
	'Anterolateral visual',...Area
	'Lateral visual',...Area
	'Rostrolateral area',...Area
	'Anterior area',...Area
	'Subiculum',...Area
	'Field CA1',...Area
	'Field CA2',...Area
	'Field CA3',...Area
	'Dentate gyrus',...Area
	'Retrosplenial'...Area
	};

dblFrameDur = 0.0099;
vecPlotWindow = [-1 1.5];
vecBinsT = vecPlotWindow(1):dblFrameDur:vecPlotWindow(2);
matMeanMoveX = nan(numel(vecBinsT),1);
matMeanMoveY = nan(numel(vecBinsT),1);
matMeanRadius = nan(numel(vecBinsT),1);
cellAggRawData = cell(3,1);
cellAggNormData = cell(3,1);
vecAggLabels = [];
vecAggRec = [];
matR2 = [];
matP = [];

%% go through files
for intFile=1:numel(cellFiles)
	strTarget = cellFiles{intFile};
	strRec = strTarget((end-10):(end-7));
	strSource = strTarget(1:(end-7));
	sLoad = load(fullfile(strDataSourcePath,strTarget));
	sAP = sLoad.sAP;
	sLoad = load(fullfile(strPupilSourcePath,cellPupilFiles{intFile}));
	sPupil = sLoad.sPupil;
	
	%% run pupil prep header
	fprintf('Running %d/%d: %s [%s]\n',intFile,numel(cellFiles),strSource,getTime());
	runPupilHeaderHybrid;
	%outputs:
	matMeanAct;
	vecOrientation;
	vecPupilStimOnTime;
	vecGlitchesPerTrials;
	vecStimOnTime;
	vecPupilSize;
	vecPupilLocX;
	vecPupilLocY;
	vecPupilSize_Change;
	vecPupilLocX_Change;
	vecPupilLocY_Change;
	
	%select data
	if boolPlotChange
		vecR = vecPupilSize_Change;
		vecX = vecPupilLocX_Change;
		vecY = vecPupilLocY_Change;
		cellVarName = {'pupil change','horz-movement','vert-movement'};
		strType = 'Change';
	else
		vecR = vecPupilSize;
		vecX = vecPupilLocX;
		vecY = vecPupilLocY;
		cellVarName = {'pupil size','horz-location','vert-location'};
		strType = 'Loc';
	end
	
	%% decode pupil size,horizontal pupil location,vertical pupil location
	% regroup
	vecRunAreaIdx = groupInto(cellClustAreas,cellRunAreas);
	[varDataOut1,vecUnique,vecCounts,cellSelect] = label2idx(vecRunAreaIdx);
	
	%assign neurons that belong to an area with fewer than 5 neurons to "other"
	indRemGroups = vecUnique == 0 | vecCounts < 5;
	indRemove = any(cell2mat(cellSelect(indRemGroups)),2);
	vecUnique(indRemGroups) = [];
	
	%select & z-score
	indKeepNeurons = (range(matMeanAct,2)>0) & (~indRemove);
	intNeurons = sum(indKeepNeurons);
	varTypeCV = 1;
	dblLambda = 100;
	matMeanAct = matMeanAct(indKeepNeurons,:);
	vecRunAreaIdx = vecRunAreaIdx(indKeepNeurons);
	matMeanZ = zscore(matMeanAct,[],1)';
	if boolAllPlots
		figure
		maxfig;
	end
	for intVar=1:3
		if intVar==1
			matY = (vecR)';
		elseif intVar==2
			matY = (vecX)';
		else
			matY = (vecY)';
		end
		[matY_hat,dblR2_CV,matB] = doCrossValidatedDecodingRR(matMeanZ,matY,varTypeCV,dblLambda);
		strVar = cellVarName{intVar};
		
		%get R^2
		[dblR2,dblSS_tot,dblSS_res,dblT,dblP] = getR2(matY,matY_hat,intNeurons);
		if boolAllPlots
			% plot
			vecLimAx = [min(cat(1,matY_hat,matY)) max(cat(1,matY_hat,matY))];
			subplot(2,3,intVar);
			h2=plot(vecLimAx,vecLimAx,'--','Color',[0.5 0.5 0.5]);
			hold on
			h1=scatter(matY,matY_hat,[],lines(1),'.');
			hold off
			xlabel(sprintf('Real %s (z-score)',strVar));
			ylabel(sprintf('CV-decoded %s (z-score)',strVar));
			title(sprintf('%s; %s; CV R^2=%.3f, p=%.3f',strRec,strVar,dblR2,dblP));
			fixfig;
			h1.SizeData=72;
			h2.LineWidth=1;
			
			%plot
			intPlotAreas = numel(vecUnique);
			subplot(2,3,intVar+3);hold on;
			dblMaxVal = max(abs(matB));
			dblStep = 0.01;
			dblEnd = roundi(dblMaxVal+dblStep,2);
			vecBins=0:dblStep:dblEnd;
			vecBinCenters = (dblStep/2):dblStep:(dblEnd);
			cellLabel = cell(1,intPlotAreas);
			
			for intPlotArea=1:intPlotAreas
				intIdx = vecUnique(intPlotArea);
				if intIdx == 0
					strArea = 'other';
				else
					strArea = cellRunAreas{intIdx};
				end
				vecThisB = abs(matB(vecRunAreaIdx==intIdx,1));
				cellLabel{intPlotArea} = strArea;
				%vecC = histcounts(vecThisB,vecBins);
				%plot(vecBinCenters,vecC/sum(vecC));
				bplot(vecThisB,intPlotArea,'points');
			end
			ylabel('Absolute regression coeffs per neuron');
			hold off
			%legend(cellLabelX);
			set(gca,'xtick',1:intPlotAreas,'xticklabel',cellLabel);
			xtickangle(gca,20);fixfig;
			fixfig;
		end
		%% concatenate data across recordings
		cellAggRawData{intVar} = cat(1,cellAggRawData{intVar},abs(matB));
		cellAggNormData{intVar} = cat(1,cellAggNormData{intVar},abs(matB)/std(abs(matB)));
		matR2(intVar,intFile) = dblR2;
		matP(intVar,intFile) = dblP;
	end
	%% save figure
	if boolAllPlots
		strFig = [strSource 'Pupil' strType];
		drawnow;
		export_fig([strOutputPath strFig '.tif']);
		export_fig([strOutputPath strFig '.pdf']);
	end
	
	%% concatenate ids
	vecAggLabels = cat(1,vecAggLabels,vecRunAreaIdx);
	vecAggRec= cat(1,vecAggRec,intFile*ones(size(vecRunAreaIdx)));
end

%% get settings
%boolPlotChange = false;
if boolPlotChange
	cellVarName = {'pupil change','horz-movement','vert-movement'};
	strType = 'Change';
else
	cellVarName = {'pupil size','horz-location','vert-location'};
	strType = 'Loc';
end

boolNorm = false;
if boolNorm
	strRegType = 'Normalized';
	cellAggData = cellAggNormData;
else
	strRegType = 'Absolute';
	cellAggData = cellAggRawData;
end
[varDataOut1,vecUnique,vecCounts,cellSelect] = label2idx(vecAggLabels);
[varSplitRec,vecUniqueRec,vecCountsRec,cellSelectRec] = label2idx(vecAggRec);

%remove data
indRemData = all(matP>0.05,1);

%% plot decoding performance averaged over recordings
h=figure;
intPlotVars = 2;
hold on;
for intVar=1:intPlotVars
	%reorder
	strVar = cellVarName{intVar};
	dblMeanR2 = mean(matR2(intVar,:));
	dblSdR2 = std(matR2(intVar,:));
	
	%plot
	bplot(matR2(intVar,~indRemData),intVar,'points');
end
hold off
ylim([0 1]);
xlabel('Decoded variable');
ylabel('Decoding performance (CV R^2)');
%legend(cellLabelX);
set(gca,'xtick',1:intPlotVars,'xticklabel',cellVarName(1:intPlotVars));
title(['Average decoding performance']);
fixfig;
%% save figure
strFig = ['DecodingPerfPupil' strType '_' getDate()];
drawnow;
export_fig([strOutputPath strFig '.tif']);
export_fig([strOutputPath strFig '.pdf']);

%% grand average decoding weights
h=figure;
maxfig;
intPlotVars = 2;
for intVar=1:intPlotVars
	%reorder
	strVar = cellVarName{intVar};
	vecMeans = splitapply(@median,cellAggData{intVar},varDataOut1);
	[vecSortedMeans,vecReorder]=sort(vecMeans,'ascend');
	vecAreaSorted = vecUnique(vecReorder);
	
	%plot
	intPlotAreas = numel(vecAreaSorted);
	subplot(1,intPlotVars,intVar);hold on;
	cellLabel = cell(1,intPlotAreas);
	for intPlotArea=1:intPlotAreas
		intIdx = vecAreaSorted(intPlotArea);
		if intIdx == 0
			strArea = 'other';
		else
			strArea = cellRunAreas{intIdx};
		end
		vecThisB = cellAggData{intVar}(vecAggLabels == intIdx);
		cellLabel{intPlotArea} = strArea;
		%vecC = histcounts(vecThisB,vecBins);
		%plot(vecBinCenters,vecC/sum(vecC));
		bplot(vecThisB,intPlotArea,'horiz');
	end
	xlabel([strRegType ' regression coeffs per neuron']);
	hold off
	%legend(cellLabelX);
	set(gca,'ytick',1:intPlotAreas,'yticklabel',cellLabel);
	%xtickangle(gca,20);fixfig;
	title(['Decoding contribution ' strVar]);
	fixfig;
end
jFig = maxfig(h,1,0.75);

%% save figure
strFig = ['RegressionCoeffsPupil' strType '_' getDate()];
drawnow;
export_fig([strOutputPath strFig '.tif']);
export_fig([strOutputPath strFig '.pdf']);

%% save data
strDataOut = ['DataGrandAveragePupil' strRegType strType getDate() '.mat'];
save([strOutputPath strDataOut],'boolPlotChange','boolNorm','vecAggLabels','vecAggRec','matR2','cellAggRawData','cellAggNormData','cellRunAreas');