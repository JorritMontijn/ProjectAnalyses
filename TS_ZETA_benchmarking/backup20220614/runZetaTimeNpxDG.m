%% load spiking data & plot tuning curves

%% set recording
close all;
clear all;
strDataSourcePath = 'F:\Data\Processed\Neuropixels\';
strDataTargetPath = 'D:\Data\Processed\TraceZeta\';
strFigPath = 'D:\Data\Results\TraceZeta\';
vecRandTypes = [1 2];%1=normal,2=rand
vecTestResamps = 100;
intResampNum = 100;
boolSave = true;

vecRunAreas = 8;%[7:24];%[8];%[1 8];%[7:24];%[1:4];%1:6;%1:5;
cellRunStim = {'','RunDriftingGratings','RunNaturalMovie'};
vecRunStim = 2;%2:3;

%% set variables
strRec = sprintf('NpxResamp%d',intResampNum);

%% set recording
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
				[sAggStim,sAggNeuron]=loadDataNpx(strArea,strRunStim,strDataSourcePath);
				if isempty(sAggStim),continue;end
				%take only BL6
				indUseNeurons = contains({sAggNeuron.SubjectType},'BL6','IgnoreCase',true);
				sAggNeuron = sAggNeuron(indUseNeurons);
				cellRecIdx = {sAggStim.Rec};
				intNeurons = numel(sAggNeuron);
				strArea = replace([lower(strArea) strRunStim],lower(cellRepStr(:,1)),cellRepStr(:,2));
				if intNeurons < 10,continue;end
			end
			
			if intRandType == min(vecRandTypes)
				%% pre-allocate output variables
				matTtest = nan(intNeurons,2);
				matZeta = nan(intNeurons,2);
				matAnova = nan(intNeurons,2);
				matZetaTime = nan(intNeurons,2,numel(vecTestResamps));
				matCompTimeZeta = nan(intNeurons,2);
				matCompTimeZetaTime = nan(intNeurons,2);
				matCompTimeAnova = nan(intNeurons,2);
			end
			
			%% analyze
			hTicN = tic;
			for intNeuron=[1:intNeurons]%31 %20    24    62    79    97   117
				%% message
				if toc(hTicN) > 5
					fprintf('Processing %s neuron %d/%d [%s]\n',strArea,intNeuron,intNeurons,getTime);
					hTicN=tic;
				end
				
				%% load or generate data
				if contains(strRunType,cellUniqueAreas(7:end),'IgnoreCase',true)
					%% get neuronal data
					sThisNeuron = sAggNeuron(intNeuron);
					vecSpikeTimes = sThisNeuron.SpikeTimes;
					strRecIdx = sThisNeuron.Rec;
					strBlock = '';
					strDate = sThisNeuron.Date;
					intSU = sThisNeuron.Cluster;
					intClust = sThisNeuron.IdxClust;
					
					%% get matching recording data
					sThisRec = sAggStim(strcmpi(strRecIdx,cellRecIdx));
					vecStimOnTime = [];
					vecStimOffTime = [];
					for intRec=1:numel(sThisRec.cellBlock)
						if ~isempty(vecStimOnTime) || numel(sThisRec.cellBlock{intRec}.vecStimOnTime) ~= 480,continue;end
						vecStimOnTime = cat(2,vecStimOnTime,sThisRec.cellBlock{intRec}.vecStimOnTime);
						vecStimOffTime = cat(2,vecStimOffTime,sThisRec.cellBlock{intRec}.vecStimOffTime);
					end
					
					matTrialStart = [];
					matTrialStart(:,1) = vecStimOnTime;
					matTrialStart(:,2) = vecStimOffTime;
					
				end
				%subsample trials
				matTrialStart = matTrialStart(1:120,:);
				
				%get trial dur
				dblUseTrialDur = round(median(diff(matTrialStart(:,1)))*2)/2;
				%set derivative params
				if intRandType ==2
					vecJitter = 2*dblUseTrialDur*((rand(size(matTrialStart(:,1)))-1/2)*2);
					matUseTrialT = bsxfun(@plus,matTrialStart,vecJitter);
					strRand = 'Rand';
				else
					matUseTrialT = matTrialStart;
					strRand = 'Real';
				end
				
				intPlot = 0;
				hTicZ = tic;
				[dblZetaP,sZETA] = zetatest(vecSpikeTimes,matUseTrialT,dblUseTrialDur,intResampNum,0);
				dblTimeZ = toc(hTicZ);
				dblMeanP = sZETA.dblMeanP;
				
				%% anova
				hTicA=tic;
				[vecTrialPerSpike,vecTimePerSpike] = getSpikesInTrial(vecSpikeTimes,matUseTrialT(:,1),dblUseTrialDur);
				intTrials = size(matUseTrialT,1);
				if numel(vecTimePerSpike) < 3,continue;end
				[optN, dblC, allN, allC] = opthist(vecTimePerSpike);
				if optN==1,optN=2;end %at least 2 bins
				dblBinWidth = dblUseTrialDur/optN;
				vecBins = 0:dblBinWidth:dblUseTrialDur;
				matPSTH = nan(intTrials,numel(vecBins)-1);
				for intTrial=1:intTrials
					matPSTH(intTrial,:) = histcounts(vecTimePerSpike(vecTrialPerSpike==intTrial),vecBins);
				end
				dblAnovaP=anova1(matPSTH,[],'off');
				dblAnovaDur = toc(hTicA);
				
				%% calculate IFR of random jitters to set threshold for significant epochs
				%get average of multi-scale derivatives, and rescaled to instantaneous spiking rate
				hTicTZ = tic;
				intRandIters = numel(sZETA.cellRandT);
				matRandRate = nan(numel(sZETA.vecSpikeT),intRandIters);
				vecT = sZETA.vecSpikeT;
				%figure
				%hold on
				for intRandIter=1:intRandIters
					vecRandT = sZETA.cellRandT{intRandIter};
					vecRandD = sZETA.cellRandDiff{intRandIter};
					
					
					dblMeanRandRate = (numel(vecRandT)/(dblUseTrialDur*size(matUseTrialT,1)));
					dblBase = 3;
					[vecRandRate,sRate] = getMultiScaleDeriv(vecRandT,vecRandD,[],[],dblBase,0,dblMeanRandRate,dblUseTrialDur);
					vecInterpR = interp1(vecRandT,vecRandRate,vecT);
					matRandRate(:,intRandIter) = vecInterpR;
					%plot(vecT,vecInterpR,'color',[0.5 0.5 0.5]);
				end
				%get real rate
				dblMeanRealRate = (numel(vecT)/(dblUseTrialDur*size(matUseTrialT,1)));
				[vecRealRate,sRate] = getMultiScaleDeriv(vecT,sZETA.vecD,[],[],dblBase,0,dblMeanRealRate,dblUseTrialDur);
				dblTimeTZ = toc(hTicTZ);
				
				%% test multiple resamp nums
				for intResampIdx=1:numel(vecTestResamps)
					intUseResamps = vecTestResamps(intResampIdx);
					vecRandRate = flat(matRandRate(:,randperm(intResampNum,intUseResamps)));
					
					[vecQ,dblZETA] = getZetaP(vecRealRate,vecRandRate,false);
					vecP = 1-abs(1-(vecQ)*2);
					vecP_corr = (1 - (1 - vecP).^(numel(vecQ))); %sidak
					
					%new
					vecP_corr3 = nan(size(vecQ));
					indP_down = vecQ<0.5;
					vecP_corr3(indP_down) = (1 - (1 - vecQ(indP_down)).^(numel(vecQ)));
					%vecP_corr3(indP_down) = (1 - (1 - 2*vecQ(indP_down)).^(numel(vecRandRate)));
					indP_up = vecQ>=0.5;
					vecP_corr3(indP_up) = 1 - ((1 - (1 - vecQ(indP_up))).^numel(vecQ));
					%vecP_corr3(indP_up) = 1 - ((1 - (1 - vecQ(indP_up))*2).^numel(vecRandRate));
					
					%extra
					intType = 2;
					if intType==2
						vecZ = -norminv(vecP_corr); %double the statistic for two-sided test
						%vecZ(vecZ<0)=0;
						vecP_corr2 = (1-normcdf(vecZ))*2; %transform back to p-value
					end
					vecZ = -norminv(vecP_corr/2); %double the statistic for two-sided test
					vecZ2 = -norminv(vecP_corr2/2); %double the statistic for two-sided test
					vecZ3 = -norminv(vecP_corr3/2); %double the statistic for two-sided test
					
					%save
					dblZetaTimeP = min(vecP_corr);
					dblZetaTimeP2 = min(vecP_corr2);
					dblZetaTimeP3 = min(vecP_corr3);
					
					%save
					matZetaTime(intNeuron,intRandType,intResampIdx) = dblZetaTimeP3;
				end
				vecZ_corr = -norminv(vecP_corr/2);
				%real spiking
				if 0
					%% plot
					figure
					hold on
					for intRandIter=1:2;%intRandIters
						plot(vecT,matRandRate,'color',[0.5 0.5 0.5]);
					end
					hold off
					colormap(gca,parula);
					set(gca,'clim',[min(vecZ_corr),max(vecZ_corr)]);
					hold on
					vecH = cline(gca,sZETA.vecSpikeT,vecRealRate,[],vecZ_corr,true);
					indSign = vecP_corr<0.05;
					scatter(sZETA.vecSpikeT(indSign),vecRealRate(indSign),'k.');
					hold off
					h=colorbar;
					clabel(h,'Z-score');
					xlabel('Time (s)')
					ylabel('Spiking rate (Hz)')
					title(sprintf('ZETA+, p=%.3f',dblZetaTimeP3))
					fixfig;
					return
				end
				
				%% save data
				matTtest(intNeuron,intRandType) = dblMeanP;
				matZeta(intNeuron,intRandType) = dblZetaP;
				matAnova(intNeuron,intRandType) = dblAnovaP;
				matCompTimeZeta(intNeuron,intRandType) = dblTimeZ;
				matCompTimeZetaTime(intNeuron,intRandType) = dblTimeTZ;
				matCompTimeAnova(intNeuron,intRandType) = dblAnovaDur;
				
			end
			
			
			
		end
		%% save
		if boolSave
			save([strDataTargetPath 'ZetaTime' strArea strRec '.mat' ],...
				'matAnova','matCompTimeAnova','vecTestResamps','matTtest','matZeta','matZetaTime','matCompTimeZeta','matCompTimeZetaTime','intResampNum','strRec','strArea');
		end
		
	end
end
