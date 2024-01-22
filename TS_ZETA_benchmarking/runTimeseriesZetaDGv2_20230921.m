%% load spiking data & plot tuning curves

%% set recording
%close all;
clear all;
if isfolder('F:\Drive\MontijnHeimel_TimeseriesZeta')
	strPath = 'F:\Drive\MontijnHeimel_TimeseriesZeta';
else
	strPath = 'C:\Drive\MontijnHeimel_TimeseriesZeta';
end
strDataPath = fullfile(strPath,'\Data\');
strFigPath = fullfile(strPath,'\Figs\');

vecRunTypes = [1 2];
intResampNum = 500;
boolSave = true;%true;
dblUseDur = 8;
boolDirectQuantile = false;
intUseTrials = inf; %limit number of used trials to reduce performance saturation
dblSuperResFactor = 100; %1 or 100
warning('off','zetatstest:InsufficientDataLength');

%% load data
for boolDoOGB = false%[false true]
	%% load data
	if boolDoOGB
		strIndicator = 'OGB';
		strDataSourcePath = 'D:\Data\Processed\PlaidsAndGratings\Gratings\';
		% load database
		sDir=dir(fullpath(strDataSourcePath,'*.mat'));
		cellRunRecs = {sDir.name};
		
	else
		strIndicator = 'GCaMP';
		strDataSourcePath = 'D:\Data\Processed\imagingGCaMP\';
		
		% load database
		sLoad=load([strDataSourcePath 'SessionDatabase.mat']);
		vecCheck = find(~cellfun(@isempty,sLoad.cellRecGratings));
		cellRunRecs = cell(1,numel(vecCheck));
		for i=1:numel(vecCheck)
			intEntry = vecCheck(i);
			cellRecs = sLoad.cellRecGratings{intEntry};
			cellRunRecs{i} = cellRecs{floor(numel(cellRecs)/2)};
		end
	end
	
	for intRunType=vecRunTypes
		%% load data
		for intFile=1:numel(cellRunRecs)
			strFile = cellRunRecs{intFile};
			[dummy,strRec,strExt]=fileparts(strFile);
			sLoad = load([strDataSourcePath strFile]);
			ses = sLoad.ses;
			intNeurons = numel(ses.neuron);
			cellData = {ses.neuron(:).dFoF};
			
			if intRunType == 1
				strRunType = strRec;
			elseif intRunType ==2
				strRunType = [strRec '-Rand'];
			end
			
			%% pre-allocate output variables
			vecZetaP = nan(1,intNeurons);
			vecMeanP = nan(1,intNeurons);
			vecAnovaP = nan(1,intNeurons);
			vecAnovaDur = nan(1,intNeurons);
			vecZetaDur = nan(1,intNeurons);
			vecKsP = nan(1,intNeurons);
			vecWilcoxP = nan(1,intNeurons);
			vecKsDur = nan(1,intNeurons);
			vecWilcoxDur = nan(1,intNeurons);
		
			%% analyze
			hTic = tic;
			for intNeuron=1:intNeurons%[1:intNeurons]%43%1:27, 2:69
				%% message
				if toc(hTic) > 5
					fprintf('Processing neuron %d/%d [%s]\n',intNeuron,intNeurons,getTime);
					hTic=tic;
				end
				clear vecTrialStarts;
				%% load or generate data
				strRecIdx = ses.session;
				strMouse = getFlankedBy(ses.experiment,'_','.mat','last');
				strBlock = sprintf('xyt%02d',ses.recording);
				strDate = ses.session;
				intSU = intNeuron;
				intClust = intNeuron;
				
				%get data
				intTotTrials = numel(ses.structStim.FrameOn);
				if isnan(intUseTrials) || intUseTrials > intTotTrials
					intUseTrials = intTotTrials;
				end
				vecOn = ses.structStim.FrameOn(1:intUseTrials);
				vecOff = ses.structStim.FrameOff(1:intUseTrials);
				vecdFoF = cellData{intNeuron};
				vecTraceT = (1:numel(vecdFoF)) ./ ses.samplingFreq;
				vecStimOnTime = vecOn(:) ./ ses.samplingFreq;
				if vecOff(1) == vecOn(2)
					vecStimOffTime = vecStimOnTime + median(diff(vecStimOnTime))/2;
				else
					vecStimOffTime = vecOff(:) ./ ses.samplingFreq;
				end
				%put in output
				vecTraceAct = vecdFoF;
				matTrialStarts = [];
				matTrialStarts(:,1) = vecStimOnTime-3;
				matTrialStarts(:,2) = vecStimOnTime;
				dblBaseDur = roundi(mean(vecStimOnTime(2:end) - vecStimOffTime(1:(end-1))),0);
				
				%% remove superfluous data
				%sampling interval
				dblSamplingInterval = median(diff(vecTraceT));
				dblUseMaxDur = dblUseDur;
				vecUseDur = [dblUseMaxDur 6];
				strUseDur = ['Dur' num2str(dblUseDur)];
				dblJitterSize = 1;
				vecEventStarts = matTrialStarts(:,1);
				%discard leading/lagging data
				dblPreUse = -dblUseMaxDur*((dblJitterSize-1)/2);
				dblPostUse = dblUseMaxDur*((dblJitterSize+1)/2);
				
				dblStartT = min(vecEventStarts) + dblPreUse*2 - dblUseMaxDur;
				dblStopT = max(vecEventStarts) + dblPostUse*2 + dblUseMaxDur;
				indRemoveEntries = (vecTraceT < dblStartT) | (vecTraceT > dblStopT);
				vecTraceT(indRemoveEntries) = [];
				vecTraceAct(indRemoveEntries) = [];
				intTrials = numel(vecEventStarts);
				return
				%% get visual responsiveness
				%set derivative params
				if contains(strRunType,'Rand')
					vecJitterPerTrial = dblJitterSize*linspace(dblUseMaxDur/intTrials,dblUseMaxDur,intTrials)';
					vecStartT =  vecEventStarts + vecJitterPerTrial(randperm(numel(vecJitterPerTrial)));
					matEventTimes = cat(2,vecStartT,vecStartT+matTrialStarts(:,2)-matTrialStarts(:,1));
				else
					matEventTimes = cat(2,vecEventStarts,vecEventStarts+matTrialStarts(:,2)-matTrialStarts(:,1));
				end
				
				%check length
				dblMaxEndDur = max(vecTraceT) - matEventTimes(end,1);
				if dblUseMaxDur > dblMaxEndDur
					matEventTimes(end,:) = [];
				end
				
				%ANOVA
				hTicA = tic;
				[vecRefT2,matTracePerTrial] = getTraceInTrial(vecTraceT,vecTraceAct,matEventTimes(:,1),dblSamplingInterval,dblUseMaxDur);
				dblBinWidth = median(diff(vecRefT2));
				dblAnovaP=anova1(matTracePerTrial,[],'off');
				dblAnovaDur = toc(hTicA);
				vecAnovaP(intNeuron) = dblAnovaP;
				vecAnovaDur(intNeuron) = dblAnovaDur;
				
				%TS-ZETA
				hTicZ = tic;
				intPlot = 0;
				%continue;
				[dblZetaP,sZETA] = zetatstest(vecTraceT,vecTraceAct,matEventTimes,vecUseDur,intResampNum,intPlot,boolDirectQuantile,dblJitterSize,dblSuperResFactor);
				%pause
				% assign data
				dblMeanP = sZETA.dblMeanP;
				dblMeanZ = -norminv(dblMeanP/2);
				dblZetaZ = sZETA.dblZETA;
				vecZetaP(intNeuron) = dblZetaP;
				vecMeanP(intNeuron) = dblMeanP;
				dblZetaDur = toc(hTicZ);
				vecZetaDur(intNeuron) = dblZetaDur;
				
				%kruskal-wallis
				hTicK = tic;
				dblKsP = kruskalwallis(matTracePerTrial,[],'off');
				vecKsP(intNeuron) = dblKsP;
				vecKsDur(intNeuron) = toc(hTicK);
				
				%wilcoxon ranksum
				if ~isempty(sZETA.vecMu_Dur) && ~isempty(sZETA.vecMu_Pre)
					hTicW = tic;
					vecWilcoxP(intNeuron) = ranksum(sZETA.vecMu_Dur,sZETA.vecMu_Pre);
					vecWilcoxDur(intNeuron) = toc(hTicW);
				end
				
				if 0%dblZetaZ > 3 && dblMeanZ < 1.5
					[dblZetaP,sZETA] = zetatstest(vecTraceT,vecTraceAct,matEventTimes,vecUseDur,intResampNum,2); %16
					
					subplot(2,3,5)
					hold on
					h1=bplot(sZETA.vecMeanBase,1);
					h2=bplot(sZETA.vecMeanStim,2);
					ylabel('Mean dF/F0');
					set(gca,'xtick',[1 2],'xticklabel',{'Base','Stim'});
					title(sprintf('%s - neuron %d',strRec,intNeuron),'interpreter','none');fixfig;
					pause
					
					if 0
						%%
						export_fig(fullpath(strFigPath,sprintf('Example_%s_Cell%d.tif',strRec,intNeuron)));
						export_fig(fullpath(strFigPath,sprintf('Example_%s_Cell%d.pdf',strRec,intNeuron)));
						
					end
				end
			end
			fprintf('%s; Mean comp time per neuron was %.3fs\n',strRec,mean(vecZetaDur));
			
			if boolSave
				if dblSuperResFactor == 1
					strSR = 'SR1';
				else
					strSR = '';
				end
				save([strDataPath 'TsZeta' strIndicator '_Q' num2str(boolDirectQuantile) '_' strRunType strUseDur 'T' num2str(intUseTrials) 'Resamp' num2str(intResampNum) strSR '.mat' ],...
					'intUseTrials','vecWilcoxP','vecWilcoxDur','vecKsP','vecKsDur',...
					'vecAnovaP','vecZetaP','vecMeanP','vecAnovaDur','vecZetaDur','strRunType','strRecIdx');
			end
		end
	end
end