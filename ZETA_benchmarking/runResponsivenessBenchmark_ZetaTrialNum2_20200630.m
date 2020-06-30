%% load spiking data & plot tuning curves

%% set recording
%close all;
clear all;
cellUniqueAreas = {...
	'V1',...Area 1
	'SC',...Area 2
	'Poisson',...Area 3
	'Retina',...Area 4
	'CaNM',...Area 5
	'CaDG',...Area 6
	'lateral geniculate',...Area 7
	'Primary visual',...Area 8
	'Lateral posterior nucleus',...Area 9
	'Anterior pretectal nucleus',...Area 10
	'Nucleus of the optic tract',...Area 11
	'Superior colliculus',...Area 12
	'Anteromedial visual',...Area 13
	'posteromedial visual',...Area 14
	'Anterolateral visual',...Area 15
	'Lateral visual',...Area 16
	'Rostrolateral area',...Area 17
	'Anterior area',...Area 18
	'Subiculum',...Area 19
	'Field CA1',...Area 20
	'Field CA2',...Area 21
	'Field CA3',...Area 22
	'Dentate gyrus',...Area 23
	'Retrosplenial'...Area 24
	};

strDisk = 'F:';
strDataSourcePath = [strDisk '\Data\Processed\Neuropixels\'];
strDataTargetPath = [strDisk '\Data\Processed\ZETA\TrialNum\'];
strFigPath = [strDisk '\Data\Results\ZETA\Examples\'];
intMakePlots =0; %0=none, 1=normal plot, 2=including raster
intLatencyPeaks = 0;
vecRandTypes = [1 2];%1=normal,2=rand
vecRestrictRange = [0 inf];
boolSave = true;
intResampleNum = 100;%10:10:90;%[10:10:100];
vecSubsampleTrials = 100:100:4000;
intSubsampleNum = numel(vecSubsampleTrials);
vecRunAreas = 8%7:16%[8];%[1 8];%[7:24];%[1:4];%1:6;%1:5;
cellRunStim = {'','RunDriftingGratings','RunNaturalMovie'};
vecRunStim = 2;%2:3;
cellRepStr = {...
	'RunDriftingGratings','-DG';...
	'RunNaturalMovie','-NM';...
	'lateral geniculate','LGN';...
	'Primary visual','V1';...
	'Lateral posterior nucleus','LP';...
	'Anterior pretectal nucleus','APN';...
	'Nucleus of the optic tract','NOT';...
	'Superior colliculus','SC';...
	'Anteromedial visual','AM';...
	'posteromedial visual','PM';...
	'Anterolateral visual','AL';...
	'Lateral visual','L';...
	'Rostrolateral area','RL';...
	'Anterior area','A';...
	'Subiculum','H-Sub';...
	'Field CA1','H-CA1';...
	'Field CA2','H-CA2';...
	'Field CA3','H-CA3';...
	'Dentate gyrus','H-DG';...
	'Retrosplenial','RetSpl';...
	};
vecTrialNum = [];
%set var
for intArea=vecRunAreas
	if intArea < 7
		intStim = 1;
	else
		intStim = 2;
	end
	for intRandType=vecRandTypes
		%reset vars
		clearvars -except intStim intSubsampleNum vecSubsampleTrials vecTrialNum intResampleNum cellRepStr intRandType vecRandTypes intRunStim vecRunStim cellRunStim intArea vecRunAreas cellUniqueAreas boolSave vecResamples strDataSourcePath strDataTargetPath strFigPath intMakePlots vecRunTypes
		strArea = cellUniqueAreas{intArea};
		strRunStim = cellRunStim{intStim};
		
		if intRandType == 1
			strRunType = strArea;
		elseif intRandType ==2
			strRunType = [strArea '-Rand'];
		end
		
		%% load data
		if contains(strRunType,cellUniqueAreas(7:end),'IgnoreCase',true)
			%% find data
			sFiles = dir([strDataSourcePath '*.mat']);
			cellFiles = {sFiles(:).name}';
			strName = replace([lower(strArea) strRunStim],lower(cellRepStr(:,1)),cellRepStr(:,2));
			
			%% go through files
			clear sAggStim;
			clear sAggNeuron;
			intNeurons = 0;
			for intFile=1:numel(cellFiles)
				%% load
				fprintf('Loading %s [%s]\n',cellFiles{intFile},getTime);
				sLoad = load([strDataSourcePath cellFiles{intFile}]);
				sAP = sLoad.sAP;
				intNewFile = 1;
				%check if neuron is in target area
				for intClust=1:numel(sAP.sCluster)
					strClustArea = sAP.sCluster(intClust).Area;
					if ~isempty(strClustArea) && contains(strClustArea,strArea,'IgnoreCase',true) && (sAP.sCluster(intClust).KilosortGood || sAP.sCluster(intClust).Contamination < 0.1)
						%% aggregate data
						%check if stim type is present
						indUseStims = ismember(cellfun(@(x) x.structEP.strFile,sAP.cellStim,'uniformoutput',false),strRunStim);
						if isempty(indUseStims) || ~any(indUseStims)
							continue;
						end
						%add data
						if intNeurons == 0
							intNewFile = 0;
							sAggNeuron(1) = sAP.sCluster(intClust);
							sAggStim(1).cellStim = sAP.cellStim(indUseStims);
							sAggStim(1).Rec = sAggNeuron(end).Rec;
						elseif ~isempty(indUseStims) && any(indUseStims)
							sAggNeuron(end+1) = sAP.sCluster(intClust);
						end
						if intNewFile
							sAggStim(end+1).cellStim = sAP.cellStim(indUseStims);
							sAggStim(end).Rec = sAggNeuron(end).Rec;
							intNewFile = 0;
						end
						intNeurons = intNeurons + 1;
					end
				end
			end
			if ~exist('sAggStim','var')
				continue;
			end
			
			cellRecIdx = {sAggStim.Rec};
			fprintf('Found %d cells from %d recordings in "%s" [%s]\n',intNeurons,numel(cellRecIdx),strRunType,getTime);
			
		end
		
		%% message
		hTic=tic;
		
		%% pre-allocate output variables
		matUseP = nan(intNeurons,intSubsampleNum,24);
		matAllP = nan(intNeurons,intSubsampleNum,24);
		matCorrP = nan(intNeurons,intSubsampleNum);
		
		%% analyze
		for intNeuron=[1:intNeurons]%31
			%% get neuronal data
			sThisNeuron = sAggNeuron(intNeuron);
			vecSpikeTimes = sThisNeuron.SpikeTimes;
			strRecIdx = sThisNeuron.Rec;
			strMouse = sThisNeuron.Mouse;
			strBlock = '';
			strArea = strName;
			strDate = sThisNeuron.Date;
			intSU = sThisNeuron.Cluster;
			intClust = sThisNeuron.IdxClust;
			
			%% get matching recording data
			sThisRec = sAggStim(strcmpi(strRecIdx,cellRecIdx));
			vecStimOnTime = [];
			vecStimOffTime = [];
			vecOrientation = [];
			for intRec=1:numel(sThisRec.cellStim)
				vecOrientation = cat(2,vecOrientation,sThisRec.cellStim{intRec}.structEP.Orientation);
				vecStimOnTime = cat(2,vecStimOnTime,sThisRec.cellStim{intRec}.structEP.vecStimOnTime);
				vecStimOffTime = cat(2,vecStimOffTime,sThisRec.cellStim{intRec}.structEP.vecStimOffTime);
			end
			
			vecTrialStarts = [];
			vecTrialStarts(:,1) = vecStimOnTime;
			vecTrialStarts(:,2) = vecStimOffTime;
			
			%% jitter
			%get trial dur
			dblUseMaxDur = round(median(diff(vecTrialStarts(:,1)))*2)/2;
			%set derivative params
			if contains(strRunType,'Rand')
				dblDur = dblUseMaxDur;
				vecJitter = 2*dblDur*rand([numel(vecStimOnTime) 1])-dblDur;
				matEventTimes = bsxfun(@plus,vecTrialStarts,vecJitter);
			else
				matEventTimes = vecTrialStarts;
			end
			
			%% build response vectors per trial
			vecRespPre = [];
			vecRespDur = [];
			vecRespBinsDur = sort(flat([matEventTimes(:,1) matEventTimes(:,2)]));
			vecR = histcounts(vecSpikeTimes,vecRespBinsDur);
			vecD = diff(vecRespBinsDur)';
			vecHz_Dur = vecR(1:2:end)./vecD(1:2:end);
			
			vecRespBinsPre = sort(flat([matEventTimes(:,1)-0.3 matEventTimes(:,1)]));
			vecR = histcounts(vecSpikeTimes,vecRespBinsPre);
			vecD = diff(vecRespBinsPre)';
			vecHz_Pre = vecR(1:2:end)./vecD(1:2:end);
			
			for intSubsampleTrialIdx=1:numel(vecSubsampleTrials)
				
				%% message
				if toc(hTic) > 5
					fprintf('Processing neuron %d/%d, %d trials [%s]\n',intNeuron,intNeurons,vecSubsampleTrials(intSubsampleTrialIdx),getTime);
					hTic=tic;
				end
				
				
				%% subsample
				vecTrialNum = unique([vecTrialNum size(matEventTimes,1)]);
				% subsample trials
				intTrialNum = vecSubsampleTrials(intSubsampleTrialIdx);
				if intTrialNum > numel(vecHz_Dur),break;end
				vecUseTrials = sort(randperm(numel(vecHz_Dur),intTrialNum));
				
				%% split by orientation
				dblStep = 360/24;
				vecBinOriEdges = (-dblStep/2):dblStep:(360-dblStep/2);
				[vecReps,dummy,dummy,dummy,cellIDs] = makeBins(vecOrientation(vecUseTrials),vecOrientation(vecUseTrials),vecBinOriEdges);
				vecSubPre = vecHz_Dur(vecUseTrials);
				vecSubDur = vecHz_Pre(vecUseTrials);
				intUseComps = sum(vecReps>1);
				vecP = ones(size(vecReps));
				for intOri=find(vecReps(:)'>1)
					vecDur = vecHz_Dur(cellIDs{intOri});
					vecPre = vecHz_Pre(cellIDs{intOri});
					
					[h,vecP(intOri)] = ttest(vecDur,vecPre);
				end
				
				dblCorrP = min(vecP*intUseComps);
				
				% assign data
				matUseP(intNeuron,intSubsampleTrialIdx,:) = vecReps(:)'>1;
				matAllP(intNeuron,intSubsampleTrialIdx,:) = vecP;
				matCorrP(intNeuron,intSubsampleTrialIdx) = dblCorrP;
			end
			
		end
		if boolSave
			save([strDataTargetPath 'ZetaDataTrialNum2' strRunType strRunStim '.mat' ],...
				'vecSubsampleTrials','matUseP','matAllP','matCorrP');
		end
	end
end
