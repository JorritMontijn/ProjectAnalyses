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
strDataTargetPath = [strDisk '\Data\Processed\ZETA\Latencies\'];
strFigPath = [strDisk '\Data\Results\ZETA\Latencies\'];
intMakePlots =4; %0=none, 1=normal plot, 2=including raster
vecRandTypes = [1 2];%1=normal,2=rand
vecRestrictRange = [0 inf];
boolSave = true;
vecBinDurs = [(2.^(-10:0))*0.1];
vecRunAreas = 8;%[7:16];%[7:24];%[1:4];%1:6;%1:5;
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
			intBinNum = numel(vecBinDurs);
			vecNumSpikes = nan(1,intNeurons);
			matBinLatencies = nan(intBinNum,intNeurons);
			vecZetaLatencies = nan(1,intNeurons);
			vecSU = nan(1,intNeurons);
			
			%% message
			fprintf('Processing %s, # of bins = %d [%s]\n',strRunType,intBinNum,getTime);
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
				
				[dblZetaP,sZETA,sRate,vecLatencies] = zetatest(vecSpikeTimes,matEventTimes,dblUseMaxDur,100,intMakePlots,4,vecRestrictRange);
				vecZetaLatencies(intNeuron) = vecLatencies(4);
				vecSU(intNeuron) = intSU;
				if 0%intMakePlots > 0
					strTit = sprintf('%s-N%dSU%d',strRunType,intNeuron,intSU);
					title(subplot(2,3,2),strTit);
					drawnow;
					pause
					export_fig([strFigPath strTit '.tif']);
					export_fig([strFigPath strTit '.pdf']);
					boolSave = false;
					continue;
				end
				%%{
				%% calculate mean-rate difference
				%pre-allocate
				intMaxRep = size(matEventTimes,1);
				vecEventStops = matEventTimes(:,2);
				vecStimHz = zeros(intMaxRep,1);
				vecBaseHz = zeros(intMaxRep,1);
				dblMedianBaseDur = median(matEventTimes(2:end,1) - matEventTimes(1:(end-1),2));
				
				%go through trials to build spike time vector
				for intEvent=1:intMaxRep
					%get times
					dblStartT = matEventTimes(intEvent,1);
					dblStopT = dblStartT + dblUseMaxDur;
					dblPreT = dblStartT - dblMedianBaseDur;
					
					% build trial assignment
					vecStimHz(intEvent) = sum(vecSpikeTimes < dblStopT & vecSpikeTimes > dblStartT)/(dblStopT - dblStartT);
					vecBaseHz(intEvent) = sum(vecSpikeTimes < dblStartT & vecSpikeTimes > dblPreT)/dblMedianBaseDur;
				end
				
				%get metrics
				dblMeanD = mean(vecStimHz - vecBaseHz) / ( (std(vecStimHz) + std(vecBaseHz))/2);
				[h,dblMeanP]=ttest(vecStimHz,vecBaseHz);
				
				figure
				bplot([vecBaseHz vecStimHz])
				
				xlabel('Base/stim')
				ylabel('Firing rate (Hz)')
				xlim([0 3]);
				fixfig;
				title(sprintf('p=%.3f',dblMeanP))
				%}
				%% get bin-wise approach
				%get data
				for intBinIdx=1:intBinNum
					
					dblFrameDur = vecBinDurs(intBinIdx);
					dblStimDur = median(diff(vecStimOnTime));
					vecBinEdges = 0:dblFrameDur:1.5;
					vecBinCenters = vecBinEdges(2:end)-dblFrameDur/2;
					intBins = numel(vecBinCenters);
					intTrials = numel(vecStimOnTime);
					matResp = nan(intTrials,intBins);
					for intTrial=1:intTrials
						matResp(intTrial,:) = histcounts(vecSpikeTimes,vecBinEdges+matEventTimes(intTrial,1));
					end
					
					%test
					vecR = mean(matResp,1);
					%get MSD peak
					[dblPeakRate,dblPeakTime,dblPeakWidth,vecPeakStartStop,intPeakLoc,vecPeakStartStopIdx] = getPeak(vecR,vecBinCenters,vecRestrictRange);
					
					
					if ~isnan(dblPeakTime)
						%get onset
						[dblOnset,dblOnsetVal] = getOnset(vecR,vecBinCenters,dblPeakTime,vecRestrictRange);
					end
					
					matBinLatencies(intBinIdx,intNeuron) = dblOnset;
				end
				
				%% MIMI
				%vecCoeffs0 = sMIMI.FitCoeffs;
				vecCoeffs0 = [];
				[dblMIMI_P,vecLatencies,sMIMI,sRate] = getMIMI(vecSpikeTimes,matEventTimes,dblUseMaxDur,intMakePlots,2,vecRestrictRange,[],[],[],vecCoeffs0);
				vecMimiLatencies(intNeuron) = vecLatencies(2);
				
				if intMakePlots > 0
					strTit = sprintf('MIMI_%s-N%dSU%d',strRunType,intNeuron,intSU);
					title(subplot(2,3,2),[strTit sprintf('; onset=%.2f',sRate.dblOnset*1000)]);
					hold(subplot(2,3,2),'on');
					scatter(subplot(2,3,2),sRate.dblOnset,sRate.dblPeakRate,'x')
					drawnow;
					
					export_fig([strFigPath strTit '.tif']);
					export_fig([strFigPath strTit '.pdf']);
					boolSave = false;
					%continue;
				end
				
				%% save output
				% assign data
				vecNumSpikes(intNeuron) = numel(vecSpikeTimes);
				
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
				save([strDataTargetPath 'ZetaDataBinsLatencies2' strRunType strRunStim '.mat' ],...
					'vecSU','vecBinDurs','matBinLatencies','vecZetaLatencies','vecNumSpikes');
			end
		end
	end
end