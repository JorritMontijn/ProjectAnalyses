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

vecRandTypes = [1 2];
intResampNum = 500;
boolSave = true;%true;
dblUseDur = 8;
boolDirectQuantile = false;
intSuperResFactor = 1; %1 or 100
warning('off','zetatstest:InsufficientDataLength');
vecCutOffs = 0.05./(2.^(6:-2:-4));%0.0008    0.0031    0.0125    0.0500    0.2000    0.8000
intReps = 1500;

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

for intRandType=vecRandTypes
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
		vecTsZetaP = nan(1,intNeurons);
		vecTtestP = nan(1,intNeurons);
		vecAnovaP = nan(1,intNeurons);
		matClustP = nan(1,intNeurons);
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
			
			%get data
			intTotTrials = numel(ses.structStim.FrameOn);
			vecOn = ses.structStim.FrameOn;
			vecOff = ses.structStim.FrameOff;
			vecdFoF = cellData{intNeuron};
			vecTraceT = (1:numel(vecdFoF)) ./ ses.samplingFreq;
			vecTraceAct = vecdFoF;
			vecStimOnTime = vecOn(:) ./ ses.samplingFreq;
			vecStimOffTime = vecOff(:) ./ ses.samplingFreq;
			
			%sampling interval
			dblSamplingInterval = median(diff(vecTraceT));
			
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
			dblShiftBy = -1;
			dblUseMaxDur = 6;
			
			%stim 1
			matTrialTS1 = [];
			matTrialTS1(:,1) = vecTransition2+dblShiftBy;
			matTrialTS1(:,2) = vecTransition2+dblUseMaxDur+dblShiftBy;
			intTrialsS1 = numel(vecTransition2);

			%stim 2
			matTrialTS2 = [];
			matTrialTS2(:,1) = vecTransition4+dblShiftBy;
			matTrialTS2(:,2) = vecTransition4+dblUseMaxDur+dblShiftBy;
			intTrialsS2 = numel(vecTransition4);
			
			%% randomize
			if intRandType == 2
				%if random, use 50% of trials from s1 and 50 from s2 for both sets
				vecUseS1for1 = randperm(intTrialsS1,round(intTrialsS1/2));
				vecUseS2for1 = randperm(intTrialsS2,round(intTrialsS2/2));
				vecUseS1for2 = find(~ismember(1:intTrialsS1,vecUseS1for1));
				vecUseS2for2 = find(~ismember(1:intTrialsS2,vecUseS2for1));
				
				matTrialT1 = cat(1,matTrialTS1(vecUseS1for1,:),matTrialTS2(vecUseS2for1,:));
				matTrialT2 = cat(1,matTrialTS1(vecUseS1for2,:),matTrialTS2(vecUseS2for2,:));
			else
				%if not random, use normal data
				matTrialT1 = matTrialTS1;
				matTrialT2 = matTrialTS2;
			end
			
			%% run tests
			boolDirectQuantile=0;
			intPlot = 0;
			[dblZeta2P,sZETA] =  zetatstest2(vecTraceT,vecTraceAct,matTrialT1,vecTraceT,vecTraceAct,matTrialT2,dblUseMaxDur,intResampNum,intPlot,boolDirectQuantile,intSuperResFactor);
			if 0%sZETA.dblMeanP > 0.1 && dblZeta2P < 0.01
				intPlot = 4;
				[dblZeta2P,sZETA] = zetatstest2(vecTraceT,vecTraceAct,matTrialT1,vecTraceT,vecTraceAct,matTrialT2,dblUseMaxDur,intResampNum,intPlot,boolDirectQuantile,intSuperResFactor);
				subplot(2,3,5)
				bplot([sZETA.vecMu1' sZETA.vecMu2'],'outliers');
				ylabel('dF/F0');
				set(gca,'xtick',[1 2],'xticklabel',{'Pref stim','Orth stim'});
				title(sprintf('%s N%d',strRecIdx,intNeuron1));
				fixfig;
				pause
			end
			
			vecTsZetaP(intNeuron) = dblZeta2P;
			
			%ANOVA
			hTicA = tic;
			[vecRefT1,matTracePerTrial1] = getTraceInTrial(vecTraceT,vecTraceAct,matTrialT1(:,1),dblSamplingInterval,dblUseMaxDur);
			[vecRefT2,matTracePerTrial2] = getTraceInTrial(vecTraceT,vecTraceAct,matTrialT2(:,1),dblSamplingInterval,dblUseMaxDur);
			
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
			
			%t-test
			vecMu1 = mean(matTracePerTrial1,2);
			vecMu2 = mean(matTracePerTrial2,2);
			[h,dblMeanP]=ttest2(vecMu1,vecMu2);
			vecTtestP(intNeuron) = dblMeanP;
			
			%cluster analysis
			for intCutOffIdx=1:numel(vecCutOffs)
				dblCutOff = vecCutOffs(intCutOffIdx);
				[dblClustP,sClustPos,sClustNeg] = clustertest(matTracePerTrial1,matTracePerTrial2,intReps,[],dblCutOff);
				matClustP(intNeuron,intCutOffIdx) = dblClustP;
			end
			
		end
		if boolSave
			save([strDataPath 'TsZeta2NM_' strRunType 'Resamp' num2str(intResampNum) '.mat' ],...
				'vecCutOffs','matClustP','vecAnovaP','vecTsZetaP','vecTtestP','strRunType','strRecIdx');
		end
	end
end