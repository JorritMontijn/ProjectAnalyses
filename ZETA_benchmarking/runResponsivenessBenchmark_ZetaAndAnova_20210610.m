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

strDisk = 'F:\';
strDataMasterPath = fullfile(strDisk,'\Data\Processed\ePhys\');
strDataTargetPath = fullfile(strDisk,'\Data\Processed\ZETA\Inclusion\');
strFigPath = fullfile(strDisk,'\Data\Results\ZETA\Examples\');
intMakePlots =0; %0=none, 1=normal plot, 2=including raster
vecRandTypes = [1 2];%[1 2];%1=normal,2=rand
vecRestrictRange = [0 inf];
boolSave = true;
vecResamples = 100;%10:10:90;%[10:10:100];
vecRunAreas = 8;%[1 2 4 6 7:14];%[1 2 4 6]% 7:14];%7:16;%4;%6;%14:16;%7:16%[8];%[1 8];%[7:24];%[1:4];%1:6;%1:5;
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
				strName = replace([lower(strArea) strRunStim],lower(cellRepStr(:,1)),cellRepStr(:,2));
				[sAggStim,sAggNeuron]=loadDataNpx(strArea,strRunStim);
				if isempty(sAggStim),continue;end
				cellRecIdx = {sAggStim.Rec};
				intNeurons = numel(sAggNeuron);
				
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
				cellData = {ses.neuron(:).apFrames};
			elseif contains(strRunType,'CaDG')
				strDataSourcePath = 'F:\Data\Processed\imagingGCaMP\';
				sFiles = dir([strDataSourcePath '20150511xyt01_ses.mat']);
				cellFiles = {sFiles(:).name}';
				sLoad = load([strDataSourcePath cellFiles{1}]);
				ses = sLoad.ses;
				intNeurons = numel(ses.neuron);
				vecOn = ses.structStim.FrameOn;
				cellData = {ses.neuron(:).apFrames};
			end
			
			for intResampleIdx = 1:numel(vecResamples)
				intResampleNum = vecResamples(intResampleIdx);
				%% message
				fprintf('Processing %s, resampling %d (%d/%d) [%s]\n',strRunType,intResampleNum,intResampleIdx,numel(vecResamples),getTime);
				hTic=tic;
				
				%% pre-allocate output variables
				cellNeuron = cell(1,intNeurons);
				vecNumSpikes = ones(1,intNeurons);
				vecZetaP = ones(1,intNeurons);
				vecAnovaP = ones(1,intNeurons);
				vecZetaTime = ones(1,intNeurons);
				vecAnovaTime = ones(1,intNeurons);
				vecTtestP = ones(1,intNeurons);
				vecTtestTime = ones(1,intNeurons);
				load([strDataTargetPath 'ZetaDataAnova' strRunType strRunStim '.mat']);
				
				%% analyze
				for intNeuron=1:intNeurons%31
					%% message
					if toc(hTic) > 5
						fprintf('Processing neuron %d/%d [%s]\n',intNeuron,intNeurons,getTime);
						hTic=tic;
					end
					clear vecTrialStarts;
					
					%% prep data
					if contains(strRunType,cellUniqueAreas(7:end),'IgnoreCase',true)
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
						for intRec=1%:numel(sThisRec.cellStim)
							vecStimOnTime = cat(2,vecStimOnTime,sThisRec.cellStim{intRec}.structEP.vecStimOnTime);
							vecStimOffTime = cat(2,vecStimOffTime,sThisRec.cellStim{intRec}.structEP.vecStimOffTime);
						end
						
						vecTrialStarts = [];
						vecTrialStarts(:,1) = vecStimOnTime;
						vecTrialStarts(:,2) = vecStimOffTime;
					elseif contains(strRunType,'V1') || contains(strRunType,'SC')
						%% get neuronal data
						sThisNeuron = sAggNeuron(intNeuron);
						vecSpikeTimes = sThisNeuron.SpikeTimes;
						strRecIdx = sThisNeuron.RecIdx;
						strMouse = sThisNeuron.Mouse;
						strBlock = sThisNeuron.Block;
						strArea = sThisNeuron.Area;
						strDate = sThisNeuron.Date;
						intSU = sThisNeuron.IdxSU;
						intClust = sThisNeuron.IdxClust;
						
						%% get matching recording data
						sThisRec = sAggRec(strcmpi(strRecIdx,cellRecIdx));
						vecStimOnTime = sThisRec.vecStimOnTime;
						vecStimOffTime = sThisRec.vecStimOffTime;
						vecStimOriDegrees = sThisRec.vecStimOriDegrees;
						vecStimOriRads = deg2rad(vecStimOriDegrees);
						vecEyeTimestamps = sThisRec.vecEyeTimestamps;
						matEyeData = sThisRec.matEyeData;
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
					elseif contains(strRunType,'CaNM') || contains(strRunType,'CaDG')
						%% generate
						strRecIdx = ses.session;
						strMouse = getFlankedBy(ses.experiment,'_','.mat','last');
						strBlock = sprintf('xyt%02d',ses.recording);
						strDate = ses.session;
						intSU = intNeuron;
						intClust = intNeuron;
						
						%np sub
						vecSoma_dFoF = calcdFoF(ses.neuron(intNeuron).F, ses.samplingFreq);
						vecNp_dFoF = calcdFoF(ses.neuron(intNeuron).npF, ses.samplingFreq);
						
						% calculate dFoF
						dblNpFactor = 0.7;
						vec_dFoF = vecSoma_dFoF - vecNp_dFoF*dblNpFactor;
						
						%params
						dblSpikeTau = 0.7; %0.7
						dblThresholdFactor = 0.25; %0.25
						intBlockSize = 997;
						
						%spike detection
						[vecFramesAP, vecNumberAP, vecSpikes, vecExpFit, vecSpikeTimes] = doDetectSpikes(vec_dFoF,ses.samplingFreq,dblSpikeTau,intBlockSize,dblThresholdFactor);
						
						%get data
						vecOn = ses.structStim.FrameOn;
						vecOff = ses.structStim.FrameOff;
						%{
						vecAP_frames = cellData{intNeuron};
						vecAP_T = vecAP_frames ./ ses.samplingFreq;
						%}
						vecStimOnTime = vecOn(:) ./ ses.samplingFreq;
						if vecOff(1) == vecOn(2)
							vecStimOffTime = vecStimOnTime + median(diff(vecStimOnTime))/2;
						else
							vecStimOffTime = vecOff(:) ./ ses.samplingFreq;
						end
						
						%put in output
						%vecSpikeTimes = vecAP_T + (rand(size(vecAP_T))-0.5)/ses.samplingFreq;
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
						dblDur = dblUseMaxDur;
						vecJitter = (2*dblDur*rand([numel(vecStimOnTime) 1])-dblDur);
						matEventTimes = bsxfun(@plus,vecTrialStarts,vecJitter);
					else
						matEventTimes = vecTrialStarts;
					end
					vecTrialNum = unique([vecTrialNum size(matEventTimes,1)]);
					intTrials = size(matEventTimes,1);
					intSpikeNum = numel(vecSpikeTimes);
					if intSpikeNum>50000,continue;end
					
					%if size(matEventTimes,1) > 0,continue;end
					%%{
					intGetLatencies = 0;
					hTic1 = tic;
					dblZetaP = getZeta(vecSpikeTimes,matEventTimes(:,1),dblUseMaxDur,intResampleNum,0,intGetLatencies);
					dblZetaDur = toc(hTic1);
					
					%% ANOVA
					hTic2 = tic;
					[vecTrialPerSpike,vecTimePerSpike] = getSpikesInTrial(vecSpikeTimes,matEventTimes(:,1),dblUseMaxDur);
					if numel(vecTimePerSpike) < 3,continue;end
					N0=100;
					%[optN, dblC] = opthist(vecTimePerSpike,round(N0));
					[optN, dblC, allN, allC] = opthist(vecTimePerSpike);
					if optN==1,optN=2;end %at least 2 bins
					dblBinWidth = dblUseMaxDur/optN;
					vecBins = 0:dblBinWidth:dblUseMaxDur;
					matPSTH = nan(intTrials,numel(vecBins)-1);
					for intTrial=1:intTrials
						matPSTH(intTrial,:) = histcounts(vecTimePerSpike(vecTrialPerSpike==intTrial),vecBins);
					end
					dblAnovaP=anova1(matPSTH,[],'off');
					dblAnovaDur = toc(hTic2);
					
					%%
					[vecX,vecReorder]=sort(dblUseMaxDur./allN);
					vecY = allC(vecReorder)-min(allC)+2;
					[dblMinY,intMinIdx] = min(vecY);
					figure
					hold on
					plot(vecX,vecY,'x-')
					scatter(vecX(intMinIdx),dblMinY,'xr');
					set(gca,'xscale','log')
					set(gca,'yscale','log')
					ylabel('Shimazaki-Shinomoto loss (a.u.)')
					xlabel('Bin width (s)')
					title(sprintf('Neuron %d',intNeuron))
					fixfig;grid off
					%%
					maxfig;
					subplot(2,3,1)
					plotRaster(vecSpikeTimes,matEventTimes(:,1),dblUseMaxDur,50000)
					title(sprintf('Neuron %d',intNeuron))
					grid off
					
					subplot(2,3,2)
					errorbar(vecBins(2:end)-dblBinWidth/2,mean(matPSTH./dblBinWidth,1),std(matPSTH./dblBinWidth,[],1)/sqrt(intTrials))
					ylabel('Binned spiking rate (Hz)')
					xlabel('Time after trial start (s)')
					title(sprintf('Bin width %.3f s',dblBinWidth));
					fixfig;
					
					subplot(2,3,3)
					[vecTime,vecRate,sIFR] = getIFR(vecSpikeTimes,matEventTimes(:,1),dblUseMaxDur,[],[],[],1);
					xlabel('Time after trial start (s)')
					ylabel('Instantaneous firing rate (Hz)')
					ylim([0 max(get(gca,'ylim'))]);
					title('No binning');
					
					pause;
					continue;
					%%}
					%% t-test
					%'vecTtestP','vecTtestTime'
					hTic3 = tic;
					vecRespBinsDur = sort(flat([matEventTimes(:,1) matEventTimes(:,2)]));
					vecR = histcounts(vecSpikeTimes,vecRespBinsDur);
					vecD = diff(vecRespBinsDur)';
					vecMu_Dur = vecR(1:2:end)./vecD(1:2:end);
					dblStart1 = min(vecRespBinsDur);
					dblFirstPreDur = dblStart1 - max([0 dblStart1 - median(vecD(2:2:end))]);
					dblR1 = sum(vecSpikeTimes > (dblStart1 - dblFirstPreDur) & vecSpikeTimes < dblStart1);
					vecMu_Pre = [dblR1 vecR(2:2:end)]./[dblFirstPreDur vecD(2:2:end)];
					
					%get metrics
					dblMeanD = mean(vecMu_Dur - vecMu_Pre) / ( (std(vecMu_Dur) + std(vecMu_Pre))/2);
					[h,dblTtestP]=ttest(vecMu_Dur,vecMu_Pre);
					dblTtestDur = toc(hTic3);
					
					%%
					% assign data
					cellNeuron{intNeuron} = [strArea strDate 'N' num2str(intSU)];
					vecNumSpikes(intNeuron) = intSpikeNum;
					vecZetaP(intNeuron) = dblZetaP;
					vecAnovaP(intNeuron) = dblAnovaP;
					vecZetaTime(intNeuron) = dblZetaDur;
					vecAnovaTime(intNeuron) = dblAnovaDur;
					vecTtestP(intNeuron) = dblTtestP;
					vecTtestTime(intNeuron) = dblTtestDur;
					
					%continue;
					%% build vector for cells to plot
					%% save plot
					if intMakePlots && (intResampleNum == vecResamples(end))% && ~(exist('vecTraceAct','var') && ~isempty(vecTraceAct)))
						%plot
						vecHandles = get(gcf,'children');
						ptrFirstSubplot = vecHandles(find(contains(get(vecHandles,'type'),'axes'),1,'last'));
						axes(ptrFirstSubplot);
						title(sprintf('%s-N%d, %s-%s,U%d/C%d',strName,intNeuron,strDate,strBlock,intSU,intClust));
						drawnow;
						
						strFileName = sprintf('%s%s-%sSU%dC%d',strName,strDate,strBlock,intSU,intClust);
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
				if boolSave
					save([strDataTargetPath 'ZetaDataAnova' strRunType strRunStim '.mat' ],...
						'cellNeuron','vecNumSpikes','vecZetaP','vecAnovaP','vecZetaTime','vecAnovaTime','vecTtestP','vecTtestTime');
				end
			end
		end
	end
end