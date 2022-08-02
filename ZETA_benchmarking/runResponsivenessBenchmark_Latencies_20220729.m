%% load spiking data & plot tuning curves
%{
Columns 1 through 30

     1     4     7     8    10    13    14    17    18    19    20    21    22    24    27    38    41    46    49    50    65    75    79    81    82    85    87   106   111   112
%}

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
strDataTargetPath = [strDisk '\Data\Processed\ZETA\Latencies\'];
strFigPath = [strDisk '\Data\Results\ZETA\Latencies\'];
intMakePlots =0; %0=none, 1=normal plot, 2=including raster
vecRandTypes = [1 2];%1=normal,2=rand
vecRestrictRange = [0 inf];
boolSave = true;
vecRunAreas = 8;%[10:16];%[7:24];%[1:4];%1:6;%1:5;
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
			clearvars -except strDataSourcePath vecBinDurs vecRestrictRange cellRepStr intRandType vecRandTypes intRunStim vecRunStim cellRunStim intArea vecRunAreas cellUniqueAreas boolSave vecResamples strDataMasterPath strDataTargetPath strFigPath intMakePlots vecRunTypes
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
				cellFiles(contains(cellFiles,'MA'))=[];
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
							indUseStims = ismember(cellfun(@(x) x.strExpType,sAP.cellBlock,'uniformoutput',false),strRunStim);
							if isempty(indUseStims) || ~any(indUseStims)
								continue;
							end
							%add data
							if intNeurons == 0
								intNewFile = 0;
								sAggNeuron(1) = sAP.sCluster(intClust);
								sAggStim(1).cellBlock = sAP.cellBlock(indUseStims);
								sAggStim(1).Rec = sAggNeuron(end).Rec;
							elseif ~isempty(indUseStims) && any(indUseStims)
								sAggNeuron(end+1) = sAP.sCluster(intClust);
							end
							if intNewFile
								sAggStim(end+1).cellBlock = sAP.cellBlock(indUseStims);
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
			
			
			%% pre-allocate output variables
			vecNumSpikes = nan(1,intNeurons);
			vecRateDiffAtPeak = nan(1,intNeurons);
			vecSU = nan(1,intNeurons);
			vecOnsetLatencies = nan(1,intNeurons);
			vecPeakLatencies = nan(1,intNeurons);
			vecZetaP = nan(1,intNeurons);
				
			%% message
			fprintf('Processing %s [%s]\n',strRunType,getTime);
			hTic=tic;
			
			%% analyze
			for intNeuron=1:intNeurons%[33:intNeurons]%31 [33 53][6 9] 103
				
				%% message
				if toc(hTic) > 5
					fprintf('Processing neuron %d/%d [%s]\n',intNeuron,intNeurons,getTime);
					hTic=tic;
				end
				clear vecTrialStarts;
				%% load or generate data
				if contains(strRunType,cellUniqueAreas(7:end),'IgnoreCase',true)
					%% get neuronal data
					sThisNeuron = sAggNeuron(intNeuron);
					vecSpikeTimes = sThisNeuron.SpikeTimes;
					strRecIdx = sThisNeuron.Rec;
					strMouse = sThisNeuron.Subject;
					strBlock = '';
					strArea = strName;
					strDate = sThisNeuron.Exp;
					intSU = sThisNeuron.Cluster;
					intClust = sThisNeuron.IdxClust;
					
					%% get matching recording data
					sThisRec = sAggStim(strcmpi(strRecIdx,cellRecIdx));
					vecStimOnTime = [];
					vecStimOffTime = [];
					for intRec=1:numel(sThisRec.cellBlock)
						vecStimOnTime = cat(2,vecStimOnTime,sThisRec.cellBlock{intRec}.vecStimOnTime);
						vecStimOffTime = cat(2,vecStimOffTime,sThisRec.cellBlock{intRec}.vecStimOffTime);
					end
					
					vecTrialStarts = [];
					vecTrialStarts(:,1) = vecStimOnTime;
					vecTrialStarts(:,2) = vecStimOffTime;
				end
				
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
				close;close;
				
				%zeta
				%{
				sOpt = struct;
		sOpt.handleFig =-1;
		[vecMean,vecSEM,vecWindowBinCenters] = doPEP(vecSpikeTimes,0:0.025:dblUseMaxDur,matEventTimes(:,1),sOpt);
		errorbar(vecWindowBinCenters,vecMean,vecSEM);
		ylim([0 max(get(gca,'ylim'))]);
		title(sprintf('Mean spiking over trials'));
		xlabel('Time from event (s)');
		ylabel('Mean spiking rate (Hz)');
		fixfig
		pause
		continue
				%}
				%[dblZetaP,sZETA,sRate,vecLatencies] = zetatest(vecSpikeTimes,matEventTimes,dblUseMaxDur,100,intMakePlots,4,vecRestrictRange);
				intMakePlots = 0;
				[dblZetaP,sZETA,sRate,vecLatencies] = zetatest(vecSpikeTimes,matEventTimes,dblUseMaxDur,100,intMakePlots,4,vecRestrictRange);
				vecOnsetLatencies(intNeuron) = vecLatencies(4);
				vecPeakLatencies(intNeuron) = vecLatencies(3);
				vecZetaP(intNeuron) = dblZetaP;
				vecSU(intNeuron) = intSU;
				if intMakePlots > 0
					strTit = sprintf('%s-N%dSU%d',strRunType,intNeuron,intSU);
					title(subplot(2,3,2),strTit);
					drawnow;
					pause
					export_fig([strFigPath strTit '.tif']);
					export_fig([strFigPath strTit '.pdf']);
					boolSave = false;
					continue;
				end
				
				%% save output
				% assign data
				vecNumSpikes(intNeuron) = numel(vecSpikeTimes);
				if isempty(sRate) || isempty(sRate.vecT) || all(isnan(vecLatencies))
					dblRateDiff = nan;
				else
					dblRateDiff = sRate.vecRate(find(vecLatencies(4) < sRate.vecT,1)-1) - mean(sRate.vecRate);
				end
				if ~isempty(dblRateDiff)
					vecRateDiffAtPeak(intNeuron) = dblRateDiff;
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
			if boolSave
				save([strDataTargetPath 'NpxLatencies' strRunType strRunStim '.mat' ],...
					'vecSU','vecOnsetLatencies','vecPeakLatencies','vecZetaP','vecRateDiffAtPeak','vecNumSpikes');
			end
		end
	end
end