%run data with NOT
%S1L3
%S2L1
%S2L4 => eye tracking wrong?
%S3L1 => probably not NOT
%S3L3
cellRun = {'S1L3','S2L1','S3L3'};
cellRun(end+1:end+2) = {'S2L4','S3L1'};
	
%get data
strDataSourcePath = 'F:\Data\Processed\Neuropixels\';
clear sFiles;
for intFile=1:numel(cellRun)
	sFile = dir([strDataSourcePath '*' cellRun{numel(cellRun)} '*.mat']);
	sFiles(intFile) = sFile(1);
end
cellFiles = {sFiles(:).name}';
close all;

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

%% go through files
for intFile=1:numel(cellFiles)
	strTarget = cellFiles{intFile};
	strRec = strTarget((end-10):(end-7));
	sLoad = load(fullfile(strDataSourcePath,strTarget));
	sAP = sLoad.sAP;
	
	%% run pupil prep header
	runPupilHeader;
	%outputs:
	matMeanAct;
	vecOrientation;
	vecPupilStimOnTime;
	vecGlitchesPerTrials;
	vecStimOnTime;
	vecPupilSize;
	vecPupilLocX;
	vecPupilLocY;

	%% decode pupil size,horizontal pupil location,vertical pupil location
	%z-score
	matMeanZ = zscore(matMeanAct,[],1)';
	varTypeCV = 1;
	dblLambda = 200;
	% regroup
	vecRunAreaIdx = groupInto(cellClustAreas,cellRunAreas);
	[varDataOut1,vecUnique,vecCounts,cellSelect] = label2idx(vecRunAreaIdx);
	
	%assign neurons that belong to an area with fewer than 5 neurons to "other"
	indRemGroups = vecCounts < 2;
	indRemove = any(cell2mat(cellSelect(indRemGroups)),2);
	vecRunAreaIdx(indRemove) = 0;
	vecUnique(indRemGroups) = [];
	figure
	maxfig;
	for intVar=1:3
		if intVar==1
			matY = zscore(vecPupilSize)';
			strVar = 'pupil size';
		elseif intVar==2
			matY = zscore(vecPupilLocX)';
			vecHorzLocZ = zscore(vecPupilLocX)';
			
			strVar = 'horz-location';
		else
			matY = zscore(vecPupilLocY)';
			strVar = 'vert-location';
		end
		[matY_hat,dblR2_CV,matB] = doCrossValidatedDecodingRR(matMeanZ,matY,varTypeCV,dblLambda);
		
		% plot
		dblLimAx = max(abs(cat(1,matY_hat,matY)));
		subplot(2,3,intVar);
		h2=plot(dblLimAx.*[-1 1],dblLimAx.*[-1 1],'--','Color',[0.5 0.5 0.5]);
		hold on
		h1=scatter(matY,matY_hat,[],lines(1),'.');
		hold off
		xlabel(sprintf('Real %s (z-score)',strVar));
		ylabel(sprintf('CV-decoded %s (z-score)',strVar));
		title(sprintf('%s; Trial-average %s; CV R^2=%.3f',strRec,strVar,dblR2_CV));
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
		cellLabelX = cell(1,intPlotAreas);
		for intPlotArea=1:intPlotAreas
			intIdx = vecUnique(intPlotArea);
			if intIdx == 0
				strArea = 'other';
			else
				strArea = cellRunAreas{intIdx};
			end
			vecThisB = abs(matB(vecRunAreaIdx==intIdx,1));
			cellLabelX{intPlotArea} = strArea;
			%vecC = histcounts(vecThisB,vecBins);
			%plot(vecBinCenters,vecC/sum(vecC));
			bplot(vecThisB,intPlotArea);
		end
		ylabel('Absolute regression coeffs per neuron');
		hold off
		%legend(cellLabelX);
		set(gca,'xtick',1:intPlotAreas,'xticklabel',cellLabelX);
		xtickangle(gca,20);fixfig;
		fixfig;
		
		%get R^2
		vecMu = mean(matY,1);
		dblSSRes_ridge = sum(sum((matY - matY_hat).^2));
		dblSSTot = sum(sum(bsxfun(@minus,matY,vecMu).^2));
		dblR2_CV = 1 - dblSSRes_ridge / dblSSTot;
		
	end
end
