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

strAnalysisType = 'A'; %A=500ms, B=300ms; 1=pooled,2=split,3=split,uncorr
strDisk = 'F:';
strDataSourcePath = [strDisk '\Data\Processed\Neuropixels\'];
strDataTargetPath = [strDisk '\Data\Processed\ZETA\TrialNum\'];
strFigPath = [strDisk '\Data\Results\ZETA\Examples\'];
intMakePlots =0; %0=none, 1=normal plot, 2=including raster
intLatencyPeaks = 0;
vecRandTypes = [1 2];%1=normal,2=rand
vecRestrictRange = [0 inf];
boolSave = true;
vecSubsampleTrials = 100:100:4000;
intSubsampleNum = numel(vecSubsampleTrials);
intResampleNum = 100;%10:10:90;%[10:10:100];
vecRunAreas = 7:16%[8];%[1 8];%[7:24];%[1:4];%1:6;%1:5;
if contains(strAnalysisType,'A')
    dblPreDur = 0.5;
elseif contains(strAnalysisType,'B')
    dblPreDur = 0.3;
end
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
		vecUseRunStim = 1;
	else
		vecUseRunStim = vecRunStim;
	end
	for intRunStim=vecUseRunStim
		for intRandType=vecRandTypes
			%reset vars
			clearvars -except strAnalysisType dblPreDur intLatencyPeaks intSubsampleNum vecSubsampleTrials vecTrialNum intResampleNum vecRestrictRange cellRepStr intRandType vecRandTypes intRunStim vecRunStim cellRunStim intArea vecRunAreas cellUniqueAreas boolSave vecResamples strDataSourcePath strDataTargetPath strFigPath intMakePlots vecRunTypes
			strArea = cellUniqueAreas{intArea};
			strRunStim = cellRunStim{intRunStim};
			
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
			matZeta = nan(intNeurons,intSubsampleNum);
			matZetaP = nan(intNeurons,intSubsampleNum);
			matHzD = nan(intNeurons,intSubsampleNum);
			matHzP = nan(intNeurons,intSubsampleNum);
			
			%% analyze
			for intNeuron=[1:intNeurons]%31
				for intSubsampleTrialIdx=1:intSubsampleNum
					
					%% message
					if toc(hTic) > 5
						fprintf('Processing neuron %d/%d, %d trials [%s]\n',intNeuron,intNeurons,vecSubsampleTrials(intSubsampleTrialIdx),getTime);
						hTic=tic;
					end
					clear vecTrialStarts;
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
					for intRec=1:numel(sThisRec.cellStim)
						vecStimOnTime = cat(2,vecStimOnTime,sThisRec.cellStim{intRec}.structEP.vecStimOnTime);
						vecStimOffTime = cat(2,vecStimOffTime,sThisRec.cellStim{intRec}.structEP.vecStimOffTime);
					end
					
					vecTrialStarts = [];
					vecTrialStarts(:,1) = vecStimOnTime;
					vecTrialStarts(:,2) = vecStimOffTime;
					
					%% get visual responsiveness
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
					
					vecTrialNum = unique([vecTrialNum size(matEventTimes,1)]);
					% subsample trials
					intTrialNum = vecSubsampleTrials(intSubsampleTrialIdx);
					if intTrialNum > size(matEventTimes,1),break;end
					matUseEventTimes = matEventTimes(sort(randperm(size(matEventTimes,1),intTrialNum)),:);
					
					%if size(matEventTimes,1) > 0,continue;end
					close;close;
					[dblZetaP,vecLatencies,sZETA,sRate] = getZeta(vecSpikeTimes,matUseEventTimes,dblUseMaxDur,intResampleNum,intMakePlots,intLatencyPeaks,vecRestrictRange);
					
					intSpikeNum = numel(vecSpikeTimes);
					
					% assign data
					matZeta(intNeuron,intSubsampleTrialIdx) = sZETA.dblZETA;
					matZetaP(intNeuron,intSubsampleTrialIdx) = dblZetaP;
					matHzD(intNeuron,intSubsampleTrialIdx) = sZETA.dblMeanD;
					matHzP(intNeuron,intSubsampleTrialIdx) = sZETA.dblMeanP;
					
					%continue;
					%% build vector for cells to plot
					%% save plot
					if intMakePlots
						%plot
						vecHandles = get(gcf,'children');
						ptrFirstSubplot = vecHandles(find(contains(get(vecHandles,'type'),'axes'),1,'last'));
						axes(ptrFirstSubplot);
						title(sprintf('%s-N%d, %s-%s,U%d/C%d',strArea,intNeuron,strDate,strBlock,intSU,intClust));
						drawnow;
						
						strFileName = sprintf('%s%s-%sSU%dC%d',strArea,strDate,strBlock,intSU,intClust);
						%export_fig([strFigPath strFileName 'Zeta' strArea 'Resamp' num2str(intResampleNum) '.tif']);
						%print(gcf, '-dpdf', [strFigPath strFileName 'Zeta' strRunType 'Resamp' num2str(intResampleNum) '.pdf']);
						%export_fig([strFigPath strFileName 'Zeta_EF_' strArea 'Resamp' num2str(intResampleNum) '.pdf']);
						pause
					end
				end
				%save
				%vecNumSpikes = nan(1,intNeurons);
				%vecZeta = nan(1,intNeurons);
				%vecP = nan(1,intNeurons);
				%vecHzD = nan(1,intNeurons);
				%vecHzP = nan(1,intNeurons);
				%cellZeta = cell(1,intNeurons);
				%cellArea = cell(1,intNeurons);
				
			end
			if boolSave
				save([strDataTargetPath 'ZetaDataTrialNum' strRunType strRunStim '.mat' ],...
					'vecSubsampleTrials','matZetaP','matZeta','matHzD','matHzP');
			end
		end
	end
end
