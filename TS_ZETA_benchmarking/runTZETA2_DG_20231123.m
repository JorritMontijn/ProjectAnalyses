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
intSuperResFactor = 1; %1 or 100
warning('off','zetatstest:InsufficientDataLength');
vecCutOffs = 0.05./(2.^(6:-2:-4));%0.0008    0.0031    0.0125    0.0500    0.2000    0.8000
intReps = 1500;

%% load data
for intCompType=2%:2
	if intCompType == 1
		strCompType = 'DiffNeurons';
	else
		strCompType = 'DiffStims';
	end
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
	
	
	for intRunType=vecRunTypes
		%% load data
		for intFile=1:numel(cellRunRecs)
			strFile = cellRunRecs{intFile};
			fprintf('Running file %d/%d: %s [%s]\n',intFile,numel(cellRunRecs),strFile,getTime);
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
			dblCutOff = inf;%0.05
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
			matClustP = nan(1,intNeurons);
			vecAnovaDur = nan(1,intNeurons);
			vecZetaDur = nan(1,intNeurons);
			
			%% analyze
			hTic = tic;
			for intNeuron1=1:intNeurons%[1:intNeurons]%43%1:27, 2:69
				%% message
				if toc(hTic) > 5
					fprintf('Processing neuron %d/%d [%s]\n',intNeuron1,intNeurons,getTime);
					hTic=tic;
				end
				clear vecTrialStarts;
				%% load or generate data
				strRecIdx = ses.session;
				strMouse = getFlankedBy(ses.experiment,'_','.mat','last');
				strBlock = sprintf('xyt%02d',ses.recording);
				strDate = ses.session;
				intSU = intNeuron1;
				intClust = intNeuron1;
				
				%get pref ori trials
				vecOrientation = ses.structStim.Orientation;
				vecPrefOriTrials = find(vecPrefOri(intNeuron1)==vecOrientation);
				vecOrthOriTrials = find(mod(vecPrefOri(intNeuron1)+90,180)==vecOrientation);
				intMaxReps = numel(vecPrefOriTrials);
				intUseReps = intMaxReps/2;
				intTrialNum = numel(ses.structStim.Orientation);
				
				indSelectPref = false(1,intMaxReps);
				if intCompType == 1
					%diff neurons
					
					%select random 50% of trials
					vecTrials1 = sort(randperm(intTrialNum,round(intTrialNum/2)));
					if intRunType == 1
						vecTrials2=vecTrials1;
						intNeuron2 = randi(intNeurons);
						while intNeuron1==intNeuron2
							intNeuron2 = randi(intNeurons);
						end
					elseif intRunType ==2
						vecTrials2 = find(~ismember(1:intTrialNum,vecTrials1));
						intNeuron2 = intNeuron1; %FPR: same neuron
					end
					
				else
					%diff stims
					intNeuron2 = intNeuron1;
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
				end
				
				%get data
				dblSamplingFreq = ses.samplingFreq;
				vecTraceAct1 = cellData{intNeuron1};
				vecTraceAct2 = cellData{intNeuron2};
				vecTraceT = (1:numel(vecTraceAct1)) ./ ses.samplingFreq;
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
				
				%% run tests
				boolDirectQuantile=0;
				intPlot = 0;
				[dblZeta2P,sZETA] =  zetatstest2(vecTraceT,vecTraceAct1,matTrialT1,vecTraceT,vecTraceAct2,matTrialT2,dblUseMaxDur,intResampNum,intPlot,boolDirectQuantile,intSuperResFactor);
				if 0%sZETA.dblMeanP > 0.1 && dblZeta2P < 0.01
					intPlot = 4;
					[dblZeta2P,sZETA] = zetatstest2(vecTraceT,vecTraceAct1,matTrialT1,vecTraceT,vecTraceAct2,matTrialT2,dblUseMaxDur,intResampNum,intPlot,boolDirectQuantile,intSuperResFactor);
					subplot(2,3,5)
					bplot([sZETA.vecMu1' sZETA.vecMu2'],'outliers');
					ylabel('dF/F0');
					set(gca,'xtick',[1 2],'xticklabel',{'Pref stim','Orth stim'});
					title(sprintf('%s N%d',strRecIdx,intNeuron1));
					fixfig;
					pause
				end
				
				vecTsZetaP(intNeuron1) = dblZeta2P;
				
				%ANOVA
				hTicA = tic;
				[vecRefT1,matTracePerTrial1] = getTraceInTrial(vecTraceT,vecTraceAct1,matTrialT1(:,1),1/dblSamplingFreq,dblUseMaxDur);
				[vecRefT2,matTracePerTrial2] = getTraceInTrial(vecTraceT,vecTraceAct2,matTrialT2(:,1),1/dblSamplingFreq,dblUseMaxDur);
				
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
				vecAnovaP(intNeuron1) = dblAnova2P;
				
				%t-test
				vecMu1 = mean(matTracePerTrial1,2);
				vecMu2 = mean(matTracePerTrial2,2);
				[h,dblMeanP]=ttest2(vecMu1,vecMu2);
				vecTtestP(intNeuron1) = dblMeanP;
				
				%cluster analysis
				for intCutOffIdx=1:numel(vecCutOffs)
					dblCutOff = vecCutOffs(intCutOffIdx);
					[dblClustP,sClustPos,sClustNeg] = clustertest(matTracePerTrial1,matTracePerTrial2,intReps,[],dblCutOff);
					matClustP(intNeuron1,intCutOffIdx) = dblClustP;
				end
				
			end
			if boolSave
				save([strDataPath 'TsZeta2_' strCompType '_' strRunType 'Resamp' num2str(intResampNum) '.mat' ],...
					'vecCutOffs','matClustP','vecAnovaP','vecTsZetaP','vecTtestP','strRunType','strRecIdx');
			end
		end
	end
end