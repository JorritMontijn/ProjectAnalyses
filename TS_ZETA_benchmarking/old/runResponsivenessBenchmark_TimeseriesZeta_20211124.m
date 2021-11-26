%% load spiking data & plot tuning curves

%% set recording
close all;
clear all;
strDataSourcePath = 'D:\Data\Processed\imagingGCaMP\';
strDataTargetPath = 'D:\Data\Processed\TraceZeta\';
vecRunTypes = [1 2];
intResampNum = 100;
boolSave = true;%true;
strFigPath = 'D:\Data\Results\TraceZeta\';

%% load database
sLoad=load([strDataSourcePath 'SessionDatabase.mat']);
vecCheck = find(~cellfun(@isempty,sLoad.cellRecGratings));
cellRunRecs = cell(1,numel(vecCheck));
for i=1:numel(vecCheck)
	intEntry = vecCheck(i);
	cellRecs = sLoad.cellRecGratings{intEntry};
	cellRunRecs{i} = cellRecs{floor(numel(cellRecs)/2)};
end


%% load data
for intRunType=vecRunTypes
	%% load data
	for intFile=1:numel(cellRunRecs)
		strFile = cellRunRecs{intFile};
		[dummy,strRec,strExt]=fileparts(strFile);
		sLoad = load([strDataSourcePath strFile]);
		ses = sLoad.ses;
		intNeurons = numel(ses.neuron);
		vecOn = ses.structStim.FrameOn;
		cellData = {ses.neuron(:).dFoF};
		
		if intRunType == 1
			strRunType = strRec;
		elseif intRunType ==2
			strRunType = [strRec '-Rand'];
		end
		
		%% pre-allocate output variables
		vecZetaP = nan(1,intNeurons);
		vecMeanP = nan(1,intNeurons);
		
		%% analyze
		hTic = tic;
		for intNeuron=1:intNeurons%[1:intNeurons]%1:27, 2:69
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
			vecOn = ses.structStim.FrameOn;
			vecOff = ses.structStim.FrameOff;
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
			
			%% stitch
			%sampling interval
			dblSamplingInterval = median(diff(vecTraceT));
			dblDur = 6;
			dblUseMaxDur = dblDur;
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
			
			%build pseudo vector
			intTrials = numel(vecEventStarts);
			vecTimestamps = vecTraceT;
			vecData = vecTraceAct;
			vecEventT = vecEventStarts+dblPreUse;
			dblWindowDur = dblPostUse-dblPreUse;
			[vecPseudoTime,vecPseudoData,vecPseudoStartT] = getPseudoTimeSeries(vecTimestamps,vecData,vecEventT,dblWindowDur);
			vecPseudoStartT = vecPseudoStartT(:);
			
			%% get visual responsiveness
			%set derivative params
			if contains(strRunType,'Rand')
				vecJitterPerTrial = dblJitterSize*linspace(dblUseMaxDur/intTrials,dblUseMaxDur,intTrials)';
				vecStartT =  vecPseudoStartT + vecJitterPerTrial(randperm(numel(vecJitterPerTrial)));
				matEventTimes = cat(2,vecStartT,vecStartT+matTrialStarts(:,2)-matTrialStarts(:,1));
			else
				matEventTimes = cat(2,vecPseudoStartT,vecPseudoStartT+matTrialStarts(:,2)-matTrialStarts(:,1));
			end
			intPlot = 0;
			%continue;
			[dblZetaP,sZETA] = getTraceZeta(vecPseudoTime,vecPseudoData,matEventTimes,dblUseMaxDur,intResampNum,intPlot); %16
			%pause
			% assign data
			dblMeanP = sZETA.dblMeanP;
			dblMeanZ = -norminv(dblMeanP/2);
			dblZetaZ = sZETA.dblZETA;
			vecZetaP(intNeuron) = dblZetaP;
			vecMeanP(intNeuron) = dblMeanP;
			
			if 0%dblZetaZ > 3 && dblMeanZ < 1.5
				[dblZetaP,sZETA] = getTraceZeta(vecPseudoTime,vecPseudoData,matEventTimes,dblUseMaxDur,intResampNum,2); %16
				
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
		
		if boolSave
			save([strDataTargetPath 'TraceZeta3GCaMP_' strRunType 'Resamp' num2str(intResampNum) '.mat' ],...
				'vecZetaP','vecMeanP','strRecIdx');
		end
	end
end
