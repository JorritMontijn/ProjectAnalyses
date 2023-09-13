%% load spiking data & plot tuning curves

%% set recording
%close all;
clear all;
cellUniqueAreas = {...
	'HeteroPoissonPeak',...Area 1
	'TriPhasic',...Area 2
	'',...Area 3
	'',...Area 4
	'',...Area 5
	'',...Area 6
	'lateral geniculate',...Area 7
	'Primary visual',...Area 8
	};

strDisk = 'F:\';
strDataMasterPath = fullfile(strDisk,'\Data\Processed\ePhys\');
strDataTargetPath = fullfile(strDisk,'\Data\Processed\ZETAv2\Inclusion\');
strFigPath = fullfile(strDisk,'\Data\Results\ZETAv2\Examples\');
intMakePlots =0; %0=none, 1=normal plot, 2=including raster
vecRandTypes = [1 2];%[1 2];%1=normal,2=rand
vecRestrictRange = [0 inf];
boolSave = true;
vecResamples = 250;%250;%10:10:90;%[10:10:100];
intUseGenN = 100;
vecRunAreas = 2;%[1 8]
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
			clearvars -except intUseGenN vecTrialNum vecRestrictRange cellRepStr intRandType vecRandTypes intRunStim vecRunStim cellRunStim intArea vecRunAreas cellUniqueAreas boolSave vecResamples strDataMasterPath strDataTargetPath strFigPath intMakePlots vecRunTypes
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
				cellRecIdx = {sAggStim.Exp};
				intNeurons = numel(sAggNeuron);
				
			elseif contains(strRunType,'Poisson')
				intNeurons = intUseGenN;
			elseif contains(strRunType,'TriPhasic')
			    intNeurons = intUseGenN;
            end
			
			for intResampleIdx = 1:numel(vecResamples)
				intResampleNum = vecResamples(intResampleIdx);
				%% message
				fprintf('Processing %s, resampling %d (%d/%d) [%s]\n',strRunType,intResampleNum,intResampleIdx,numel(vecResamples),getTime);
				hTic=tic;
				
				%% pre-allocate output variables
				cellNeuron = cell(1,intNeurons);
				vecNumSpikes = zeros(1,intNeurons);
				vecZetaP_old = ones(1,intNeurons);
				vecZetaP_new = ones(1,intNeurons);
				vecAnovaP = ones(1,intNeurons);
				vecZetaTime = ones(1,intNeurons);
				vecAnovaTime = ones(1,intNeurons);
				vecTtestP = ones(1,intNeurons);
				vecTtestTime = ones(1,intNeurons);
				%load([strDataTargetPath 'ZetaDataAnova' strRunType strRunStim '.mat']);
				
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
						strRecIdx = sThisNeuron.Exp;
						strMouse = sThisNeuron.Subject;
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
							vecStimOnTime = cat(2,vecStimOnTime,sThisRec.cellBlock{intRec}.vecStimOnTime);
							vecStimOffTime = cat(2,vecStimOffTime,sThisRec.cellBlock{intRec}.vecStimOffTime);
						end
						
						vecTrialStarts = [];
						vecTrialStarts(:,1) = vecStimOnTime;
						vecTrialStarts(:,2) = vecStimOffTime;
					    dblUseMaxDur = round(median(diff(vecTrialStarts(:,1)))*2)/2;
					elseif contains(strRunType,'HeteroPoissonPeak')
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
						vecTrialDur=linspace(0.2,2,numel(vecTrialAngles));
						vecStimOnTime = cumsum(vecTrialDur)+1;
						vecStimOffTime = vecStimOnTime + 1;
						dblJitter = 0.05;
                        boolDoublePeaked = true;
                        intAddSpikes = 0;%numel(vecTrialDur)/10;

						vecTrialStarts(:,1) = vecStimOnTime;
						vecTrialStarts(:,2) = vecStimOffTime;
						dblUseMaxDur = round(median(diff(vecTrialStarts(:,1)))*2)/2;
					    [vecSpikeTimes,dblPrefOri] = getGeneratedSpikingDataWithPeak(vecTrialAngles,vecTrialStarts,dblBaseRate,dblPrefRate,dblJitter,dblKappa,boolDoublePeaked,[],intAddSpikes);
					elseif contains(strRunType,'TriPhasic')
                        %% generate
						strRecIdx = 'x';
						strMouse = 'Artificial';
						strBlock = '1';
						strDate = getDate();
						intSU = intNeuron;
						intClust = intNeuron;
						
						%set parameters
						%m: intTrialNum
                        %Tr: response dur
                        %T: trial dur
                        %tau: total dur
                        %L_b: base rate
                    	%L_s: stim rate
                        %rate is 0 between Tr and T

                        intNumT = 160;
                        vecTrialDur=linspace(1.1,3,intNumT);
						
                        dblBaseRate = exprnd(1);
						dblFirstRate = exprnd(1);
						dblThirdRate = exprnd(10);
						dblUseMaxDur = 1;
                        dblFirstDur = 0.1;
                        m=1;
                        Tr = 0.5;
                        T = 1;
                        
                        L_b=dblThirdRate;
                        L_s=dblFirstRate;
                        vecTrialStarts = nan(size(vecTrialDur))';
                        %generate preceding period
                        dblPreT=5;
                        cellSpikeT = cell(1,intNumT+3);
                        dblNextT = dblPreT;
                        for intTrial=1:intNumT
                            tau = vecTrialDur(intTrial);
                            cellSpikeT{intTrial} = dblNextT+getGeneratedTriPhasicR(m,dblFirstDur,dblUseMaxDur,tau,L_b,L_s*dblFirstDur);
                            vecTrialStarts(intTrial) = dblNextT;
                            dblNextT = dblNextT + tau;
                        end
                        vecTrialStarts(:,2) = vecTrialStarts(:,1) + dblUseMaxDur;
                        %generate end period
                        dblPostT=5;
                        %add baseline spikes
                        dblTotT=dblNextT+dblPostT;
                        vecSpikeTimesBase = getGeneratedSpikingData(0,[0;dblTotT],dblBaseRate,dblBaseRate,10);
                        cellSpikeT{intNumT+1} = vecSpikeTimesBase;
                        vecSpikeTimes = sort(cell2vec(cellSpikeT));
                        pNewStitch=zetatest(vecSpikeTimes,vecTrialStarts,dblUseMaxDur,[],2,[],[],[],true)
                        pNewNoStitch=zetatest(vecSpikeTimes,vecTrialStarts,dblUseMaxDur,[],2,[],[],[],false)
                        pOld=getZeta(vecSpikeTimes,vecTrialStarts,dblUseMaxDur,[],2)
                    end
					
					%% get visual responsiveness
					%get trial dur
					%set derivative params
					if contains(strRunType,'Rand')
						dblDur = dblUseMaxDur;
						vecJitter = (2*dblDur*rand([size(vecTrialStarts,1) 1])-dblDur);
						matEventTimes = bsxfun(@plus,vecTrialStarts,vecJitter);
					else
						matEventTimes = vecTrialStarts;
					end
					vecTrialNum = unique([vecTrialNum size(matEventTimes,1)]);
					intTrials = size(matEventTimes,1);
					intSpikeNum = numel(vecSpikeTimes);
					if intSpikeNum>50000 || intSpikeNum<3,continue;end
					
					%if size(matEventTimes,1) > 0,continue;end
					%%{
					intGetLatencies = 0;
					hTic1 = tic;
					dblZetaP_old = getZeta(vecSpikeTimes,matEventTimes(:,1),dblUseMaxDur,intResampleNum,0,intGetLatencies);
					dblZetaDur_old = toc(hTic1);
					hTic1b = tic;
					dblZetaP_new = zetatest(vecSpikeTimes,matEventTimes(:,1),dblUseMaxDur,intResampleNum,0,intGetLatencies);
					dblZetaDur_new = toc(hTic1b);
					
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
					vecZetaP_old(intNeuron) = dblZetaP_old;
					vecZetaP_new(intNeuron) = dblZetaP_new;
					vecAnovaP(intNeuron) = dblAnovaP;
					vecZetaTime(intNeuron) = dblZetaDur_new;
					vecAnovaTime(intNeuron) = dblAnovaDur;
					vecTtestP(intNeuron) = dblTtestP;
					vecTtestTime(intNeuron) = dblTtestDur;
                end
                
				if boolSave
					save([strDataTargetPath 'ZetaDataAnova' strRunType strRunStim 'Resamp' num2str(intResampleNum) '.mat' ],...
						'cellNeuron','vecNumSpikes','vecZetaP_old','vecZetaP_new','vecAnovaP','vecZetaTime','vecAnovaTime','vecTtestP','vecTtestTime');
				end
			end
		end
	end
end