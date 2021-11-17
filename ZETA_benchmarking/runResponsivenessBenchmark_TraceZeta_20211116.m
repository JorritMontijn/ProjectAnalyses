%% load spiking data & plot tuning curves

%% set recording
close all;
clear all;
strDataSourcePath = 'D:\Data\Processed\imagingGCaMP\';
strDataTargetPath = 'D:\Data\Processed\TraceZeta\';
vecRunTypes = [1 2];
intResampNum = 250;
boolSave = true;%true;

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
		for intNeuron=1:intNeurons%[1:intNeurons]%31
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
			vecTrialStarts(:,1) = vecStimOnTime-3;
			vecTrialStarts(:,2) = vecStimOnTime;
			dblBaseDur = roundi(mean(vecStimOnTime(2:end) - vecStimOffTime(1:(end-1))),0);
			
			%% get visual responsiveness
			%set derivative params
			dblDur = 6;
			if contains(strRunType,'Rand')
				vecJitter = 2*dblDur*rand(size(vecStimOnTime))-dblDur;
				matEventTimes = bsxfun(@plus,vecTrialStarts,vecJitter);
			else
				matEventTimes = vecTrialStarts;
			end
			dblUseMaxDur = dblDur;
			intPlot = 0;
			%continue;
			[dblZetaP,sZETA] = getTraceZeta(vecTraceT,vecTraceAct,matEventTimes,dblUseMaxDur,intResampNum,intPlot); %16
			%pause
			% assign data
			dblMeanP = sZETA.dblMeanP;
			vecZetaP(intNeuron) = dblZetaP;
			vecMeanP(intNeuron) = dblMeanP;
		end
		if boolSave
			save([strDataTargetPath 'TraceZeta' strRunType 'Resamp' num2str(intResampNum) '.mat' ],...
				'vecZetaP','vecMeanP','strRecIdx');
		end
	end
end
