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
boolDirectQuantile = false;
intUseTrials = inf; %limit number of used trials to reduce performance saturation
dblSuperResFactor = 100; %1 or 100
warning('off','zetatstest:InsufficientDataLength');

%% load data
strIndicator = 'GCaMP';
strDataSourcePath = 'D:\Data\Processed\imagingGCaMP\';

% load database
sLoad=load([strDataSourcePath 'SessionDatabase.mat']);
vecCheck = find(~cellfun(@isempty,sLoad.cellRecNatMovs));
cellRunRecs = cell(1,numel(vecCheck));
for i=1:numel(vecCheck)
	intEntry = vecCheck(i);
	cellRecs = sLoad.cellRecNatMovs{intEntry};
	cellRunRecs{i} = cellRecs{floor(numel(cellRecs)/2)};
end

for intRandType=vecRunTypes
	%% load data
	for intFile=1:numel(cellRunRecs)
		strFile = cellRunRecs{intFile};
		[dummy,strRec,strExt]=fileparts(strFile);
		sLoad = load([strDataSourcePath strFile]);
		ses = sLoad.ses;
		intNeurons = numel(ses.neuron);
		cellData = {ses.neuron(:).dFoF};
		
		if intRandType == 1
			strRunType = ['NM-' strRec];
		elseif intRandType ==2
			strRunType = ['NM-' strRec '-Rand'];
		end
		
		%% find transitions
		boolFindTransitionsInMovie = false;
		if boolFindTransitionsInMovie
			vecDiff = nan(1,500);
			for intFrame1=1:500
				intFrame2 = modx(intFrame1+1,500);
				vecDiff(intFrame2) = sum(flat(Earthflight_WingedPlanet__CondorFlightSchool_NarratedByDavidTen(intFrame1).cdata...
					- Earthflight_WingedPlanet__CondorFlightSchool_NarratedByDavidTen(intFrame2).cdata));
			end
			vecFrameTransitions = find(vecDiff>6e7); %1    91   251   340
		else
			vecFrameTransitions = [1 91 251 340];
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
			if all(isnan(vecdFoF)),continue;end
			vecTraceT = (1:numel(vecdFoF)) ./ ses.samplingFreq;
			vecTraceAct = vecdFoF;
			vecStimOnTime = vecOn(:) ./ ses.samplingFreq;
			vecStimOffTime = vecOff(:) ./ ses.samplingFreq;
			
			%% get transitions
			vecStimDur = vecStimOffTime-vecStimOnTime;
			vecTransition1 = vecStimOnTime+((vecFrameTransitions(1)-1)/500)*vecStimDur;
			vecTransition2 = vecStimOnTime+((vecFrameTransitions(2)-1)/500)*vecStimDur;
			vecTransition3 = vecStimOnTime+((vecFrameTransitions(3)-1)/500)*vecStimDur;
			vecTransition4 = vecStimOnTime+((vecFrameTransitions(4)-1)/500)*vecStimDur;
			dblTotDur = 500/25;
			dblFrameDur = dblTotDur/500;
			dblDur1 = median(vecTransition2-vecTransition1);
			dblDur2 = median(vecTransition3-vecTransition2);
			dblDur3 = median(vecTransition4-vecTransition3);
			dblDur4 = median(vecStimOffTime-vecTransition4);
			%dblShiftBy = -1;
			%dblUseMaxDur = 6;
			dblShiftBy = 0;
			dblUseMaxDur = dblTotDur;
			vecUseDur = [dblUseMaxDur 10];
			dblJitterSize = 1;
			
			matEventTimes = cat(2,vecTransition1(:),vecTransition2(:));
			dblMaxJitterFirstEvent = matEventTimes(1)-vecTraceT(1);
			%matEventTimes = cat(2,vecTransition2(:),vecTransition3(:));
			if intRandType == 2
				dblDur = dblUseMaxDur;
				vecJitter = (2*dblDur*rand([numel(vecTransition2) 1])-dblDur);
				while vecJitter(1) < -dblMaxJitterFirstEvent
					vecJitter(1) = 2*dblDur*rand(1)-dblDur;
				end
				matEventTimes = bsxfun(@plus,matEventTimes,vecJitter);
			else
				matEventTimes = matEventTimes;
			end
			
			%sampling interval
			dblSamplingInterval = median(diff(vecTraceT));
			
			%% get visual responsiveness
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
				[dblZetaP,sZETA] = zetatstest(vecTraceT,vecTraceAct,matEventTimes,vecUseDur,intResampNum,2,boolDirectQuantile,dblJitterSize,dblSuperResFactor); %16
				
				subplot(2,3,5)
				hold on
				h1=bplot(sZETA.vecMu_Pre,1);
				h2=bplot(sZETA.vecMu_Dur,2);
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
			save([strDataPath 'TsZeta' strRunType 'Resamp' num2str(intResampNum) strSR '.mat' ],...
				'intUseTrials','vecWilcoxP','vecWilcoxDur','vecKsP','vecKsDur',...
				'vecAnovaP','vecZetaP','vecMeanP','vecAnovaDur','vecZetaDur','strRunType','strRecIdx');
		end
	end
end