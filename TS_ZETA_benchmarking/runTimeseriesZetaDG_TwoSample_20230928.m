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
intResampNum = 250;
boolSave = true;%true;
dblUseDur = 8;
boolDirectQuantile = false;
intUseTrials = 80; %limit number of used trials to reduce performance saturation
intSuperResFactor = 1; %1 or 100
warning('off','zetatstest:InsufficientDataLength');

%% load data
for boolDoOGB = [false true]
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
			matResp = getRespMatSes(ses);
			indPreRem = any(isnan(matResp),2);
			ses.neuron(indPreRem) = [];
			
			% get tuning & discard non-tuned
			vecOrientation = ses.structStim.Orientation;
			matResp = getRespMatSes(ses);
			sOut = getTuningCurves(matResp,vecOrientation);
			dblCutOff = 0.05;
			indUseNeurons = sOut.vecFitP < dblCutOff;
			ses.neuron(~indUseNeurons) = [];
			intNeurons = numel(ses.neuron);
			cellData = {ses.neuron(:).dFoF};
			dblOriStep = median(diff(sort(sOut.vecUniqueDegs)));
			vecPrefOri = mod(round(rad2deg(sOut.matFittedParams(indUseNeurons,1))/dblOriStep)*dblOriStep,180);
			[N,edges] = histcounts(vecPrefOri,[(-dblOriStep/2):dblOriStep:180]);
			
			if intRunType == 1
				strRunType = strRec;
			elseif intRunType ==2
				strRunType = [strRec '-Rand'];
			end
			
			%% pre-allocate output variables
			vecTsZetaP = nan(1,intNeurons);
			vecTtestP = nan(1,intNeurons);
			vecAnovaP = nan(1,intNeurons);
			vecAnovaDur = nan(1,intNeurons);
			vecZetaDur = nan(1,intNeurons);
			
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
				
				%get pref ori trials
				vecOrientation = ses.structStim.Orientation;
				vecPrefOriTrials = find(vecPrefOri(intNeuron)==vecOrientation);
				vecOrthOriTrials = find(mod(vecPrefOri(intNeuron)+90,180)==vecOrientation);
				intMaxReps = numel(vecPrefOriTrials);
				intUseReps = intMaxReps/2;
				
				indSelectPref = false(1,intMaxReps);
				if intRunType == 1
					%select trials from pref and orth
					vecTrials1 = vecPrefOriTrials(randperm(intMaxReps,intUseReps));
					vecTrials2 = vecOrthOriTrials(randperm(intMaxReps,intUseReps));
				elseif intRunType ==2
					%select trials from pref and pref
					indSelectPref(randperm(intMaxReps,intUseReps)) = true;
					vecSubFirstSet = find(indSelectPref==true);
					vecSubSecondSet = find(indSelectPref==false);
					vecSubSecondSet((intUseReps+1):end)=[];
					vecTrials1 = vecPrefOriTrials(vecSubFirstSet);
					vecTrials2 = vecPrefOriTrials(vecSubSecondSet);
				end
				
				%get data
				dblSamplingFreq = ses.samplingFreq;
				vecdFoF = cellData{intNeuron};
				vecTraceT = (1:numel(vecdFoF)) ./ ses.samplingFreq;
				vecTraceAct = vecdFoF;
				dblUseMaxDur = dblUseDur;
				vecUseDur = [dblUseMaxDur 6];
				
				%set 1
				vecOn1 = ses.structStim.FrameOn(vecTrials1);
				vecOff1 = ses.structStim.FrameOff(vecTrials1);
				vecStimOnTime1 = vecOn1(:) ./ ses.samplingFreq;
				if vecOff1(1) == vecOn1(2)
					vecStimOffTime1 = vecStimOnTime1 + median(diff(vecStimOnTime1))/2;
				else
					vecStimOffTime1 = vecOff1(:) ./ ses.samplingFreq;
				end
				matTrialT1 = [];
				matTrialT1(:,1) = vecStimOnTime1;
				matTrialT1(:,2) = vecStimOffTime1;
				
				%set 2
				vecOn2 = ses.structStim.FrameOn(vecTrials2);
				vecOff2 = ses.structStim.FrameOff(vecTrials2);
				vecStimOnTime2 = vecOn2(:) ./ ses.samplingFreq;
				if vecOff2(1) == vecOn2(2)
					vecStimOffTime2 = vecStimOnTime2 + median(diff(vecStimOnTime2))/2;
				else
					vecStimOffTime2 = vecOff2(:) ./ ses.samplingFreq;
				end
				matTrialT2 = [];
				matTrialT2(:,1) = vecStimOnTime2;
				matTrialT2(:,2) = vecStimOffTime2;
				
				%% remove superfluous data
				intPlot = 0;
				boolDirectQuantile=0;
				[dblZeta2P,sZETA] = zetatstest2b(vecTraceT,vecTraceAct,matTrialT1,vecTraceT,vecTraceAct,matTrialT2,dblUseMaxDur,intResampNum,intPlot,boolDirectQuantile,intSuperResFactor);
				
				vecTtestP(intNeuron) = sZETA.dblMeanP;
				vecTsZetaP(intNeuron) = dblZeta2P;
				
				%ANOVA
				hTicA = tic;
				[vecRefT1,matTracePerTrial1] = getTraceInTrial(vecTraceT,vecTraceAct,matTrialT1(:,1),1/dblSamplingFreq,dblUseMaxDur);
				[vecRefT2,matTracePerTrial2] = getTraceInTrial(vecTraceT,vecTraceAct,matTrialT2(:,1),1/dblSamplingFreq,dblUseMaxDur);
				
				vecBin1 = flat(repmat(1:size(matTracePerTrial1,2),[size(matTracePerTrial1,1) 1]));
				vecBin2 = flat(repmat(1:size(matTracePerTrial2,2),[size(matTracePerTrial2,1) 1]));
				vecR1 = matTracePerTrial1(:);
				vecR2 = matTracePerTrial2(:);
				
				%two-sample
				g1 = cat(1,ones(size(vecR1)),2*ones(size(vecR2)));
				g2 = cat(1,vecBin1,vecBin2);
				vecP=anovan(cat(1,vecR1,vecR2),{g1,g2},'model','interaction','display','off');%,'varnames',{'g1','g2'})%
				[h crit_p adj_p]=fdr_bh(vecP([1 3]));
				dblAnova2P = min(adj_p);
		
				%one-sample diff
				vecAnovaP(intNeuron) = dblAnova2P;

			end
			if boolSave
				save([strDataPath 'TsZeta2' strIndicator '_Q' num2str(boolDirectQuantile) '_' strRunType 'Resamp' num2str(intResampNum) '.mat' ],...
					'vecAnovaP','vecTsZetaP','vecTtestP','strRunType','strRecIdx');
			end
		end
	end
end