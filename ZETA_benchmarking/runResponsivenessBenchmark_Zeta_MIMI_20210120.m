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

strDrive = 'D:';
strDataMasterPath = [strDrive '\Data\Processed\ePhys\'];
strDataTargetPath = [strDrive '\Data\Processed\ZETA\Inclusion\'];
strFigPath = [strDrive '\Data\Results\ZETA\Examples\'];
intMakePlots =0; %0=none, 1=normal plot, 2=including raster
vecRandTypes = [1 2];%1=normal,2=rand
vecRestrictRange = [0 inf];
boolSave = true;
vecResamples = 100;%10:10:90;%[10:10:100];
vecRunAreas = [16:24];%[8];%[1 8];%[7:24];%[1:4];%1:6;%1:5;
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
			clearvars -except vecTrialNum vecRestrictRange cellRepStr intRandType vecRandTypes intRunStim vecRunStim cellRunStim intArea vecRunAreas cellUniqueAreas boolSave vecResamples strDataMasterPath strDataTargetPath strFigPath intMakePlots vecRunTypes
			strArea = cellUniqueAreas{intArea};
			strRunStim = cellRunStim{intRunStim};
			
			if intRandType == 1
				strRunType = strArea;
				fprintf('Prepping normal... [%s]\n',getTime);
			elseif intRandType ==2
				strRunType = [strArea '-Rand'];
				fprintf('Prepping random... [%s]\n',getTime);
			end
			
			%% load data
			if contains(strRunType,cellUniqueAreas(7:end),'IgnoreCase',true)
				[sAggStim,sAggNeuron]=loadDataNpx(strArea,strRunStim);
				if isempty(sAggStim),continue;end
				cellRecIdx = {sAggStim.Rec};
				intNeurons = numel(sAggNeuron);
				strArea = replace([lower(strArea) strRunStim],lower(cellRepStr(:,1)),cellRepStr(:,2));
				if intNeurons < 10,continue;end
			end
			
			for intResampleIdx = 1:numel(vecResamples)
				intResampleNum = vecResamples(intResampleIdx);
				%% message
				fprintf('Processing %s, resampling %d (%d/%d) [%s]\n',strRunType,intResampleNum,intResampleIdx,numel(vecResamples),getTime);
				hTicMessage=tic;
				
				%% pre-allocate output variables
				vecComputTimeZETA = nan(1,intNeurons);
				vecComputTimeMIMI = nan(1,intNeurons);
				vecNumSpikes = nan(1,intNeurons);
				vecMIMIP = nan(1,intNeurons);
				vecZetaP = nan(1,intNeurons);
				vecHzD = nan(1,intNeurons);
				vecHzP = nan(1,intNeurons);
				cellArea = cell(1,intNeurons);
					
				%% analyze
				parfor intNeuron=[1:intNeurons]%31
					%% load or generate data
					if contains(strRunType,cellUniqueAreas(7:end),'IgnoreCase',true)
						%% get neuronal data
						sThisNeuron = sAggNeuron(intNeuron);
						vecSpikeTimes = sThisNeuron.SpikeTimes;
						strRecIdx = sThisNeuron.Rec;
						strMouse = sThisNeuron.Mouse;
						strBlock = '';
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
					
					end
					
					%% get visual responsiveness
					%get trial dur
					dblUseMaxDur = round(median(diff(vecTrialStarts(:,1)))*2)/2;
					%set derivative params
					if contains(strRunType,'Rand')
						dblDur = dblUseMaxDur;
						vecJitter = 4*dblDur*rand([numel(vecStimOnTime) 1])-dblDur;
						matEventTimes = bsxfun(@plus,vecTrialStarts,vecJitter);
					else
						matEventTimes = vecTrialStarts;
					end
					
					%vecTrialNum = unique([vecTrialNum size(matEventTimes,1)]);
					
					%ZETA
					hTicZ=tic;
					intPlot = 0;
					[dblZetaP,vecLatencies,sZETA] = getZeta(vecSpikeTimes,matEventTimes,dblUseMaxDur,intResampleNum,intPlot,0);
					dblComputTimeZETA = toc(hTicZ);
					%MIMI
					hTicM = tic;
					[dblMIMI_P,vecLatenciesMIMI,sMIMI] = getMIMI(vecSpikeTimes,matEventTimes(:,1),dblUseMaxDur,intPlot,0);
					dblComputTimeMIMI = toc(hTicM);
					
					intSpikeNum = numel(vecSpikeTimes);
					
					% assign data
					vecNumSpikes(intNeuron) = intSpikeNum;
					vecMIMIP(intNeuron) = dblMIMI_P;
					vecZetaP(intNeuron) = dblZetaP;
					vecHzD(intNeuron) = sZETA.dblMeanD;
					vecHzP(intNeuron) = sZETA.dblMeanP;
					cellArea{intNeuron} = strArea;
					vecComputTimeZETA(intNeuron) = dblComputTimeZETA;
					vecComputTimeMIMI(intNeuron) = dblComputTimeMIMI;
					
					%% message
					%if toc(hTicMessage) > 5 && intNeuron > 1
						fprintf('Processed neuron %d/%d, ZETA took %.1fs, MIMI took %.1fs [%s]\n',intNeuron,intNeurons,...
							vecComputTimeZETA(intNeuron),vecComputTimeMIMI(intNeuron),getTime);
						%hTicMessage=tic;
					%end
					%clear vecTrialStarts;
					
				end
				%save
				if boolSave
					save([strDataTargetPath 'ZetaMIMI' strRunType strRunStim 'Resamp' num2str(intResampleNum) '.mat' ],...
						'vecNumSpikes','vecMIMIP','vecZetaP','vecHzD','vecHzP','cellArea','vecComputTimeZETA','vecComputTimeMIMI');
				end
			end
		end
	end
end