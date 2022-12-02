
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

vecColAlb = [0.9 0 0];
vecColBl6 = lines(1);

%% pre-allocate
cellAreasPerExp = cell(1,numel(sExp));
cellStimsPerExp = cell(1,numel(sExp));
for intExp=1:numel(sExp)
	sRec = sExp(intExp);
	strName=[sRec.sJson.subject '_' sRec.sJson.date];
	cellAreasPerExp{intExp} = unique({sExp(intExp).sCluster.Area})';
	cellStimsPerExp{intExp} = [strName ':' char(cell2vec(cellfun(@(x) x.strExpType,sExp(intExp).cellBlock,'UniformOutput',false)))'];
end

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
		intBestRec = 17;%topo3_20220201_3
		strSubjectType = 'BL6';
		dblOffsetT=0;
	elseif intSubType == 2
		intBestRec = 5;
		strSubjectType = 'DBA';
		dblOffsetT=0;
	end
	indUseRecs = contains(cellSubjectType,strSubjectType);
	vecRunRecs = find(indUseRecs & ~(indRemRecs | indRemRecs2));
	%vecRunRecs = intBestRec;
	for intRecIdx=1:numel(vecRunRecs)
		intRec=vecRunRecs(intRecIdx)
		sRec = sExp(intRec);
		strName=[sRec.sJson.subject '_' sRec.sJson.date];
		cellStimType = cellfun(@(x) x.strExpType,sRec.cellBlock,'uniformoutput',false);
		vecBlocksDG = find(contains(cellStimType,'driftinggrating','IgnoreCase',true));
		vecBlocksNM = find(contains(cellStimType,'naturalmovie','IgnoreCase',true));
		%get timing for DG
		intBlock = vecBlocksDG(1);
		sBlock = sRec.cellBlock{intBlock};
		if isfield(sBlock,'vecPupilStimOn')
			vecPupilLatency = sBlock.vecPupilStimOn-sBlock.vecStimOnTime;
		else
			vecPupilLatency = 0;
		end
		
		for intBlockIdx=1:numel(vecBlocksNM)
			intBlock = vecBlocksNM(intBlockIdx);
			sBlock = sRec.cellBlock{intBlock};
			
			vecPupilStimOn = sBlock.vecStimOnTime+median(vecPupilLatency);
			vecPupilStimOff = sBlock.vecStimOffTime+median(vecPupilLatency);
			
			%% get pupil data
			dblSampNi = str2double(sRec.sSources.sMeta.niSampRate);
			dblFirstSamp = str2double(sRec.sSources.sMeta.firstSample);
			vecPupilTime = sRec.sPupil.vecTime;
			vecPupilLocX = sRec.sPupil.vecCenterX;
			vecPupilLocY = sRec.sPupil.vecCenterY;
			vecPupilSize = sRec.sPupil.vecRadius;
			vecPupilSync = sRec.sPupil.vecSyncLum;
			if isfield(sRec.sPupil,'vecBlinks') && ~all(sRec.sPupil.vecBlinks==0)
				vecBlinks = sRec.sPupil.vecBlinks;
			else
				%filter absvidlum
				dblLowPass = 0.01/(1/median(diff(vecPupilTime)));
				[fb,fa] = butter(2,dblLowPass,'high');
				vecAbsVidLum = zscore(filtfilt(fb,fa, sRec.sPupil.sRaw.vecPupilAbsVidLum));
				vecBlinks = vecAbsVidLum > 5;
			end
			dblFs = 1/median(diff(vecPupilTime));
			
			% get blinking
			vecBinEdges = [-inf vecPupilStimOn vecPupilStimOn(end)+median(diff(vecPupilStimOn)) inf];
			[vecCounts,vecMeans,vecSDs,cellVals,cellIDs] = makeBins(vecPupilTime,double(vecBlinks),vecBinEdges);
			vecBlinkFractionPerTrial = vecMeans(2:(end-1));
			
			%remove trials with blinking
			indRemTrials = vecBlinkFractionPerTrial > 0;
			
			%% prep data
			cellCellsPerArea = cell(1,numel(cellUseAreas));
			cellAreasPerCluster = {sRec.sCluster.Area};
			for intArea=1:numel(cellUseAreas)
				cellCellsPerArea{intArea} = contains(cellAreasPerCluster,cellUseAreas{intArea},'IgnoreCase',true);
				vecCellsNrPerArea = cellfun(@sum,cellCellsPerArea);
				
				%% select cells
				strAreaGroup =  cellAreaGroupsAbbr{intArea};
				%get data matrix
				for intType=1:2
					if intType == 1
						figure;
						vecStimOnTime = sRec.sSources.cellBlock{intBlock}.structEP.ActOnNI(~indRemTrials) - dblFirstSamp/dblSampNi;
						vecStimOffTime = sRec.sSources.cellBlock{intBlock}.structEP.ActOffNI(~indRemTrials) - dblFirstSamp/dblSampNi;
						strTit = 'NI-time';
					else
						vecStimOnTime = sBlock.vecStimOnTime(~indRemTrials);
						vecStimOffTime = sBlock.vecStimOffTime(~indRemTrials);
						strTit = 'Stim-time';
					end
					if numel(vecStimOnTime) <= 10,close;continue;end
					intPopCounter = intPopCounter + 1;
					cellSpikeT = {sRec.sCluster(:).SpikeTimes};
					
					%include?
					vecZetaP = cellfun(@min,{sRec.sCluster.ZetaP});
					%indUseCells = vecZetaP(:)<0.05 & arrayfun(@(x) x.KilosortGood==1 | x.Contamination < 0.1,sRec.sCluster(:));
					indUseCells = arrayfun(@(x) x.Violations1ms < 0.25 & abs(x.NonStationarity) < 0.25,sRec.sCluster(:));
			
					%% split cells into areas
					%build cell vectors
					vecSelectCells = find(indUseCells(:) & cellCellsPerArea{intArea}(:));
					if isempty(vecSelectCells)
						continue;
					end
					
					
					%% split data into bins
					intTypeCV = 2;
					dblLambda = 100;
					dblStimDur = median(vecStimOffTime-vecStimOnTime);
					if dblStimDur < 9
						dblMovieDur = 8; %new; actually 8.33
						intBinNr = 8;
					else
						dblMovieDur = 20; %old
						intBinNr = 20;
					end
					dblBinDur = dblMovieDur/intBinNr;
					vecBinOnset = linspace(0,dblMovieDur-dblBinDur,intBinNr);
					
					%build artificial stim vector
					intRepNum = numel(vecStimOnTime);
					vecBinOnT = flat(vecStimOnTime(:)' + vecBinOnset(:));
					vecBinIdx = flat(repmat((1:intBinNr)',[1 intRepNum]));
					matData = getSpikeCounts(cellSpikeT,vecBinOnT,dblBinDur);
					
					
					%% plot
					%[dblZetaP,vecLatencies,sZETA,sRate] = getZeta(cellSpikeT{vecSelectCells(5)},vecStimOnTime,20,[],3)
					
					%% decode movie
					matUseData = matData(vecSelectCells,:);
					vecPriorDistribution = intRepNum*ones(1,intBinNr);
					[dblPerformanceCV,vecDecodedIndexCV,matPostProbability,dblMeanErrorDegs,matConfusionML] = doCrossValidatedDecodingML(matUseData,vecBinIdx,intTypeCV,vecPriorDistribution);
					[dblPerformanceTM,vecDecodedIndexCV,matTemplateDistsCV,dblMeanErrorDegs,matConfusionTM] = doCrossValidatedDecodingTM(matUseData,vecBinIdx,intTypeCV,vecPriorDistribution);
					[dblPerformanceLR,vecDecodedIndexCV,matPostProbability,dblMeanErrorDegs,matConfusionLR] = ...
						doCrossValidatedDecodingLR(matUseData,vecBinIdx,intTypeCV,[],dblLambda);
					[dblPerformanceLR2,vecDecodedIndexCV,matPostProbability,dblMeanErrorDegs,matConfusionLR2] = ...
						doCrossValidatedDecodingLR(matUseData,vecBinIdx,intTypeCV,vecPriorDistribution,dblLambda);
					
					
					subplot(2,4,1+(intType-1)*4)
					imagesc(matConfusionML)
					title([strTit '; ' strAreaGroup '; ML: ' strName '_' num2str(intBlock)],'interpreter','none');
					
					subplot(2,4,2+(intType-1)*4)
					imagesc(matConfusionTM)
					title(['TM: ' strName '_' num2str(intBlock)],'interpreter','none');
					
					subplot(2,4,3+(intType-1)*4)
					imagesc(matConfusionLR)
					title(['LR: ' strName '_' num2str(intBlock)],'interpreter','none');
					
					subplot(2,4,4+(intType-1)*4)
					imagesc(matConfusionLR2)
					title(['LR2: ' strName '_' num2str(intBlock)],'interpreter','none');
					
					if intType == 2
						maxfig;drawnow;
						export_fig(fullpath(strTargetPath,['NatMovDecoding_' strAreaGroup '_' strName 'B' num2str(intBlock) '.jpg']));
						export_fig(fullpath(strTargetPath,['NatMovDecoding_' strAreaGroup '_' strName 'B' num2str(intBlock) '.pdf']));
					end
				end
			end
		end
	end
	%% plot
	
end
%% save
