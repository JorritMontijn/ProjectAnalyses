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


strDataMasterPath = 'F:\Data\Processed\ePhys\';
strDataTargetPath = 'F:\Data\Processed\ZETA\Inclusion\';
strFigPath = 'F:\Data\Results\ZETA\Examples\';
intMakePlots =0; %0=none, 1=normal plot, 2=including raster
vecRandTypes = 2;%[1 2];%1=normal,2=rand
vecRestrictRange = [0 inf];
boolSave = true;
vecResamples = 100;%100;%10:10:90;%[10:10:100];
vecRunAreas = 8;%10:16;%4;%6;%14:16;%7:16%[8];%[1 8];%[7:24];%[1:4];%1:6;%1:5;
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
				cellRecIdx = {sAggStim.Rec};
				intNeurons = numel(sAggNeuron);
				strArea = replace([lower(strArea) strRunStim],lower(cellRepStr(:,1)),cellRepStr(:,2));
				
			elseif contains(strRunType,'V1') || contains(strRunType,'SC')
				%% find data
				strDataSourcePath = [strDataMasterPath 'DriftingGratings\'];
				sFiles = dir([strDataSourcePath '*' strArea '_*.mat']);
				cellFiles = {sFiles(:).name}';
				
				%% go through files
				sAggRec = struct;
				sAggNeuron = struct;
				for intFile=1:numel(cellFiles)
					%% load
					sLoad = load([strDataSourcePath cellFiles{intFile}]);
					sNeuron = sLoad.sNeuron;
					sRecording = sLoad.sRecording;
					
					%% aggregate data
					if intFile == 1
						sAggRec= sRecording;
						sAggNeuron = sNeuron;
					else
						sAggRec(end+1) = sRecording;
						sAggNeuron((end+1):(end+numel(sNeuron))) = sNeuron;
					end
				end
				intNeurons = numel(sAggNeuron);
				cellRecIdx = {sAggRec.RecIdx};
			elseif contains(strRunType,'Retina')
				%% find data
				strDataSourcePath = strDataMasterPath;
				sFiles = dir([strDataSourcePath '*' strArea '_*.mat']);
				cellFiles = {sFiles(:).name}';
				sLoad = load([strDataSourcePath cellFiles{1}]);
				cellRetData = sLoad.LightFlash_For_J;
				intNeurons = size(cellRetData,1)-2;
			elseif contains(strRunType,'Poisson')
				intNeurons = 10000;
			elseif contains(strRunType,'CaNM')
				%% find data
				strDataSourcePath = 'F:\Data\Processed\imagingGCaMP\';
				sFiles = dir([strDataSourcePath '20150511xyt02_ses.mat']);
				cellFiles = {sFiles(:).name}';
				sLoad = load([strDataSourcePath cellFiles{1}]);
				ses = sLoad.ses;
				intNeurons = numel(ses.neuron);
				vecOn = ses.structStim.FrameOn;
				cellData = {ses.neuron(:).dFoF};
			elseif contains(strRunType,'CaDG')
				strDataSourcePath = 'F:\Data\Processed\imagingGCaMP\';
				sFiles = dir([strDataSourcePath '20150511xyt01_ses.mat']);
				cellFiles = {sFiles(:).name}';
				sLoad = load([strDataSourcePath cellFiles{1}]);
				ses = sLoad.ses;
				intNeurons = numel(ses.neuron);
				vecOn = ses.structStim.FrameOn;
				cellData = {ses.neuron(:).dFoF};
			end
			
			for intResampleIdx = 1:numel(vecResamples)
				intResampleNum = vecResamples(intResampleIdx);
				%% message
				fprintf('Processing %s, resampling %d (%d/%d) [%s]\n',strRunType,intResampleNum,intResampleIdx,numel(vecResamples),getTime);
				hTicMessage=tic;
				
				%% pre-allocate output variables
				vecNumSpikes = nan(1,intNeurons);
				vecZetaP = nan(1,intNeurons);
				vecHzP = nan(1,intNeurons);
				
				vecPoissP = nan(1,intNeurons);
				
				vecISIP_ks = nan(1,intNeurons);
				vecISIP_IntG = nan(1,intNeurons);
				vecISIP_g = nan(1,intNeurons);
				
				vecBISIP_g = nan(1,intNeurons);
				vecBISIP_ks = nan(1,intNeurons);
				
				vecComputTimeZETA = nan(1,intNeurons);
				vecComputTimePoiss = nan(1,intNeurons);
				vecComputTimeISI = nan(1,intNeurons);
				vecComputTimeBISI = nan(1,intNeurons);
				
				
				cellArea = cell(1,intNeurons);
				
				%% analyze
				for intNeuron=[1:intNeurons]%31
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
					
					elseif contains(strRunType,'Retina')
						%% generate
						strCellID = cellRetData{intNeuron+2,1};
						cellStr = strsplit(strCellID,'_');
						cellStr2 = strsplit(cellStr{2},'-');
						strRecIdx = cellStr2{1};
						strMouse = 'Kamermans';
						strBlock = cellStr2{1};
						strDate = cellStr{1};
						intSU = intNeuron;
						intClust = str2double(cellStr2{2}(2));
						
						%get data
						cellData = cellRetData(3:end,2);
						matN = cellData{intNeuron};
						
						%get trial timing
						intOn = find(cellRetData{3,3}(:,2)==1,1);
						dblStimOn = cellRetData{3,3}(intOn,1);
						intOff = find(cellRetData{3,3}(intOn:end,2)==0,1) + intOn - 1;
						dblStimOff = cellRetData{3,3}(intOff,1);
						dblTrialDur=cellRetData{3,3}(end-1,1);
						%build vectors
						vecAllTrialStarts = (matN(:,2)-1)*dblTrialDur;
						vecTrialStartTime = ((1:max(matN(:,2)))-1)'*dblTrialDur;
						vecStimOnTime = vecTrialStartTime + dblStimOn;
						vecStimOffTime = vecTrialStartTime + dblStimOff;
						
						%put in output
						vecSpikeTimes = matN(:,1) + vecAllTrialStarts;
						vecTrialStarts(:,1) = vecStimOnTime;
						vecTrialStarts(:,2) = vecStimOffTime;
					elseif contains(strRunType,'Poisson')
						%% generate
						strRecIdx = 'x';
						strMouse = 'Artificial';
						strBlock = '1';
						strDate = getDate();
						intSU = intNeuron;
						intClust = intNeuron;
						
						%set parameters
						dblBaseRate = exprnd(5);
						dblPrefRate = dblBaseRate+exprnd(20);
						dblKappa = rand(1)*5+5;
						vecTrialAngles=repmat([0:45:359],[1 20]);
						dblTrialDur=2;
						vecStimOnTime = dblTrialDur*(1:numel(vecTrialAngles))';
						vecStimOffTime = vecStimOnTime + 1;
						
						vecTrialStarts(:,1) = vecStimOnTime;
						vecTrialStarts(:,2) = vecStimOffTime;
						[vecSpikeTimes,dblPrefOri] = getGeneratedSpikingData(vecTrialAngles,vecTrialStarts,dblBaseRate,dblPrefRate,dblKappa,true);
					end
					
					%% get visual responsiveness
					%get trial dur
					dblUseMaxDur = round(median(diff(vecTrialStarts(:,1)))*2)/2;
					%set derivative params
					if contains(strRunType,'Rand')
						vecISI = diff(vecSpikeTimes);
						dblT0 = vecSpikeTimes(1);
						vecRandISI = vecISI(randperm(numel(vecISI)));
						vecSpikeTimes = dblT0 + cumsum(vecSpikeTimes) - vecSpikeTimes(1);
					end
					matEventTimes = vecTrialStarts;
					%vecTrialNum = unique([vecTrialNum size(matEventTimes,1)]);
					
					%vecTrialNum = unique([vecTrialNum size(matEventTimes,1)]);
					%pre-alloc
					sZETA = [];
					sZETA.dblMeanP = 1;
					dblZetaP = 1;
					dblPoissP = 1;
					dblP_KS = 1;
					dblP_G = 1;
					dblIntP_G = 1;
					dblP_BISI = 1;
					dblP_BISI_KS = 1;
					%ZETA
					hTic=tic;
					%[dblZetaP,vecLatencies,sZETA] = getZeta(vecSpikeTimes,matEventTimes,dblUseMaxDur,intResampleNum,0,0);
					dblComputTimeZETA = toc(hTic);
					
					%Poisson
					hTic=tic;
					%dblPoissP = getPoissonTest(vecSpikeTimes,matEventTimes(:,1),dblUseMaxDur,intResampleNum);
					dblComputTimePoiss= toc(hTic);
					
					%ISI varieties
					hTic = tic;
					[dblP_KS,dblP_G,dblIntP_G] = getISItest(vecSpikeTimes,matEventTimes(:,1),dblUseMaxDur,intResampleNum);
					dblComputTimeISI = toc(hTic);
					
					%BISI
					hTic = tic;
					%[dblP_BISI,dblP_BISI_KS] = getBISI(vecSpikeTimes,matEventTimes(:,1),dblUseMaxDur,intResampleNum);
					dblComputTimeBISI = toc(hTic);
					
					intSpikeNum = numel(vecSpikeTimes);
					
					% assign data
					vecNumSpikes(intNeuron) = intSpikeNum;
					vecZetaP(intNeuron) = dblZetaP;
					vecHzP(intNeuron) = sZETA.dblMeanP;
					
					vecPoissP(intNeuron) = dblPoissP;
					
					vecISIP_ks(intNeuron) = dblP_KS;
					vecISIP_IntG(intNeuron) = dblIntP_G;
					vecISIP_g(intNeuron) = dblP_G;
					
					vecBISIP_g(intNeuron) = dblP_BISI;
					vecBISIP_ks(intNeuron) = dblP_BISI_KS;
					
					cellArea{intNeuron} = strArea;
					vecComputTimeZETA(intNeuron) = dblComputTimeZETA;
					vecComputTimePoiss(intNeuron) = dblComputTimePoiss;
					vecComputTimeISI(intNeuron) = dblComputTimeISI;
					vecComputTimeBISI(intNeuron) = dblComputTimeBISI;
					
					%% message
					%if toc(hTicMessage) > 5 && intNeuron > 1
						fprintf('Processed neuron %d/%d, Compute time: ZETA=%.1fs, Poiss=%.1fs, ISI=%.1fs, BISI=%.1fs [%s]\n',intNeuron,intNeurons,...
							vecComputTimeZETA(intNeuron),...
							vecComputTimePoiss(intNeuron),...
							vecComputTimeISI(intNeuron),...
							vecComputTimeBISI(intNeuron),getTime);
						%hTicMessage=tic;
					%end
					%clear vecTrialStarts;
					
				end
				%save
				if boolSave
					save([strDataTargetPath 'ZetaVars2B' strRunType strRunStim 'Resamp' num2str(intResampleNum) '.mat' ],...
						'vecNumSpikes','cellArea',...
						'vecZetaP',...
						'vecHzP',...
						'vecPoissP',...
						'vecISIP_ks',...
						'vecISIP_IntG',...
						'vecISIP_g',...
						'vecBISIP_g',...
						'vecBISIP_ks',...
						'vecComputTimeZETA',...
						'vecComputTimePoiss',...
						'vecComputTimeISI',...
						'vecComputTimeBISI');
				end
			end
		end
	end
end