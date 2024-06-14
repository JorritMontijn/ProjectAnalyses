%% aim
%{
https://www.nature.com/articles/s41598-021-95037-z
https://www.nature.com/articles/s42003-021-02437-y

q1: Progression of orientation information over time: how does trial activity evolve? what is the function of the onset peak?
> is decoding better when matched for stimulus phase? => no.

q2: "Spike time and rate coding can be represented within a single model of spiking probability as a
function of time: rate codes are uniform over a certain period tau, while spike time codes are
temporally localized peaks"
> is this true?

q3: Rate codes do not exist; a rate code is simply a subset of spike time codes where the temporal
integration window is very large. But what about multi dim codes? Those are all rate based. Can we
formulate a multidimensional spike-time code? I.e., can we make a taxonomy of neural codes?

q4: How does information evolve over time, is initial peak indeed less tuned? Is pop activity
rhythmic? Are stimuli encoded invariant to brain state? Eg, high arousal, low arousal. Or is
stimulus manifold dynamic over time? Does manifold scale with arousal? => How does manifold depend
on binning size? What is the optimal time window?

%}

%% set parameters
clear all;%close all;
cellTypes = {'Real','ShuffTid','Uniform'};%'Real','ShuffTid'};
cellRunTypes = {'RecTopo','SimDG18'};
intRunType = 2; %topo or sim

boolFixedSpikeGroupSize = false;
intOnsetType = 0; %normal (0) or rem onset (1)
cellConnType = {'Least','Most'};

%% define qualifying areas
boolSaveFig = false;
cellUseAreas = {...
	'Primary visual area',...
	...'posteromedial visual area',...
	};
if isfolder('F:\Drive\PopTimeCoding') && isfolder('F:\Data\Processed\Neuropixels\')
	strDataPath = 'F:\Data\Processed\Neuropixels\';
	strFigurePathSR = 'F:\Drive\PopTimeCoding\single_recs';
	strFigurePath = 'F:\Drive\PopTimeCoding\figures\';
	strTargetDataPath = 'F:\Drive\PopTimeCoding\data\';
else
	strDataPath = 'E:\DataPreProcessed\';
	strFigurePathSR = 'C:\Drive\PopTimeCoding\single_recs';
	strFigurePath = 'C:\Drive\PopTimeCoding\figures\';
	strTargetDataPath = 'C:\Drive\PopTimeCoding\data\';
end

%% load data
strArea = 'V1';
vecMaxT = [1 0.5];
vecStimNum = [12 9];
vecRemOnsets = [0.25 0.125];
strRunType = cellRunTypes{intRunType}; %Sim or ABI or Npx?
dblMaxT = vecMaxT(intRunType);
intStimNum = vecStimNum(intRunType);
intQuantiles = 10;

vecBinsDur = 0:0.01:1;
vecBinsDurC = vecBinsDur(2:end) - diff(vecBinsDur(1:2))/2;
vecColGrey = [0.7 0.7 0.7];
cellSummaryMatrices = {};

%for boolRemOnset=true%[false true]
if intOnsetType == 1
	dblRemOnset = vecRemOnsets(intRunType);
else
	dblRemOnset = 0;
end
%define variables
matR = nan(8,4,2); %[r rci_lo r_ci_hi r_se]
cellSuffices = {'Dur','Conf','IFR'};
cellSufficesLong = {'Duration','Confidence','AvgRate'};
cellVars = {'Num','Rate','Tune','IFR','Conf','FrIn','PrDi','Bndw'};
cellVarsShort = {'#','R','T','V1','C','%I','PrD','B','%S','%G','%I'};
cellVarsLong = {'NumOfCells','AvgRateOfCells','AvgTuningOfCells','AvgRate','Confidence','FractionInterneurons',...
	'AvgPrefDistToStim','AvgBandwidthOfCells'};
cellVarsLegend = {'# of cells','Cell Rate','Cell Tune','V1 Rate','Conf','% InterN',...
	'Pref-d','BndW'};
intVarNum = numel(cellVars);
for intConnType=1:numel(cellConnType)
	strConnType = cellConnType{intConnType};
	
	for intType=1:numel(cellTypes)
		strType = cellTypes{intType};
		if boolFixedSpikeGroupSize
			sFiles = dir ([strTargetDataPath 'M1Data_' strRunType '*' strType '*Fixed*' strConnType '*.mat']);
			strSGS = 'FixedSGS';
		else
			sFiles = dir ([strTargetDataPath 'M1Data_' strRunType '*' strType '*Var*' strConnType '*.mat']);
			strSGS = 'VarSGS';
		end
		strSelectOnset = sprintf('%.2f',vecRemOnsets(intRunType));
		indWithOnset = cellfun(@(x) strcmp(x((end-7):end),[strSelectOnset '.mat']), {sFiles.name}');
		if dblRemOnset == 0
			indWithOnset = cellfun(@(x) strcmp(x((end-7):end),[strSelectOnset '.mat']), {sFiles.name}');
			sFiles(indWithOnset) = [];
		else
			sFiles(~indWithOnset) = [];
		end
		strOnset = sprintf('%.2f',dblRemOnset);
		intRecNum = numel(sFiles);
		cellQuantileDur = cell(intRecNum,intQuantiles);
		cellQuantileConf = cell(intRecNum,intQuantiles);
		matLatConf = nan(intRecNum,numel(vecBinsDurC));
		vecSpikeGroupSize = nan(intRecNum,1);
		
		%% plot
		figure;maxfig;
		h1=subplot(2,3,1);
		h2=subplot(2,3,2);
		h3=subplot(2,3,3);
		vecColChance = [0.3 0.3 0.3];
		plot(h1,[0 dblMaxT],([1 1]./intStimNum),'--','color',vecColChance);
		hold(h1,'on')
		title(h1,sprintf('%s, %s',strType,strSGS));
		
		plot(h2,[0 100],([1 1]./intStimNum),'--','color',vecColChance);
		hold(h2,'on')
		
		%plot(h3,[0 10],([1 1]./intStimNum),'--','color',vecColChance);
		hold(h3,'on')
		
		for intFile=1:intRecNum
			%% load
			sLoad=load(fullpath(sFiles(intFile).folder,sFiles(intFile).name));
			sSpikeGroup = sLoad.sSpikeGroup;
			strRec = sLoad.strRec;
			fprintf('Loading %s (%d/%d) [%s]\n',strRec,intFile,intRecNum,getTime);
			strType = sLoad.strType;
			vecOri180 = sLoad.vecOri180;
			vecOrientation = sLoad.vecOrientation;
			vecStimOffTime = sLoad.vecStimOffTime;
			vecStimOnTime = sLoad.vecStimOnTime;
			intSpikeGroupSize = sLoad.intSpikeGroupSize;
			vecTuningPerCell = real(sLoad.sTuningV2.vecFitT);
			vecRatePerCell = mean(sLoad.sTuningV2.matMeanResp,2);
			
			%% create derived variables
			%put confidence in deciles per recording
			vecSpikeGroupDuration = [sSpikeGroup.Duration]';
			vecSpikeGroupCorrect = [sSpikeGroup.Correct]';
			vecSpikeGroupConfidence = [sSpikeGroup.Confidence]';
			vecSpikeGroupLatency = dblRemOnset+[sSpikeGroup.Latency]';
			vecSpikeGroupNumOfCells = [sSpikeGroup.NumOfCells]';
			vecSpikeGroupAvgRateOfCells = [sSpikeGroup.AvgRateOfCells]';
			vecSpikeGroupAvgTuningOfCells = [sSpikeGroup.AvgTuningOfCells]';
			vecSpikeGroupAvgRate = [sSpikeGroup.AvgRate]';
			
			vecSpikeGroupFractionInterneurons = [sSpikeGroup.FractionInterneurons]';
			vecSpikeGroupAvgPrefDistToStim = [sSpikeGroup.AvgPrefDistToStim]';
			vecSpikeGroupAvgBandwidthOfCells = [sSpikeGroup.AvgBandwidthOfCells]';
			
			%sort by duration
			[vecSortedDur,vecSort]=sort(vecSpikeGroupDuration);
			vecSortedCorr = vecSpikeGroupCorrect(vecSort);
			vecSortedConf = vecSpikeGroupConfidence(vecSort);
			vecSortedLat = vecSpikeGroupLatency(vecSort);
			vecSortedIFR = vecSpikeGroupAvgRate(vecSort);
			
			%calculate fraction correct and confidence per bin of equal size
			intQuantileNum = 10;
			[vecMeanDur,vecSemDur,vecMeanConf,vecSemConf,vecQuantile,cellValsDur,cellValsConf]=getQuantiles(vecSortedDur,vecSortedConf,intQuantileNum);
			
			%conf with dur
			errorbar(h2,vecMeanDur*1000,vecMeanConf,vecSemConf,vecSemConf,vecSemDur,vecSemDur,'color',lines(1));
			
			%conf with time
			[vecCounts,vecMeans,vecSDs]=makeBins(vecSpikeGroupLatency,vecSpikeGroupConfidence,vecBinsDur);
			plot(h1,vecBinsDurC,vecMeans,'color',vecColGrey);%,vecSDs./sqrt(vecCounts),vecSDs./sqrt(vecCounts),'color',lines(1));
			
			%% single-group based
			%predict dur
			for intVar=1:intVarNum
				strVarShort = cellVars{intVar};
				strVarLong = cellVarsLong{intVar};
				strVarSource = ['vecSpikeGroup' strVarLong];
				vecSource = eval(strVarSource);
				
				%corrs
				[r_Var_Dur,p_Var_Dur]=corr(vecSource,vecSpikeGroupDuration);
				eval(['r_' strVarShort '_Dur=r_Var_Dur;']);
				eval(['p_' strVarShort '_Dur=p_Var_Dur;']);
				eval(['vecR_' strVarShort '_Dur(intFile) = r_Var_Dur;']);
			end
			
			tbl = table(vecSpikeGroupNumOfCells,vecSpikeGroupAvgRateOfCells,vecSpikeGroupAvgTuningOfCells,vecSpikeGroupAvgRate,vecSpikeGroupConfidence,vecSpikeGroupDuration,...
				'VariableNames',{'CellNumber','AvgCellRate','AvgCellTuning','AvgPopIFR','Conf','Dur'});
			mdl_predDur = fitlm(tbl,'linear');
			
			%predict Conf
			for intVar=1:intVarNum
				strVarShort = cellVars{intVar};
				strVarLong = cellVarsLong{intVar};
				strVarSource = ['vecSpikeGroup' strVarLong];
				vecSource = eval(strVarSource);
				
				%corrs
				[r_Var_Conf,p_Var_Conf]=corr(vecSource,vecSpikeGroupConfidence);
				eval(['r_' strVarShort '_Conf=r_Var_Conf;']);
				eval(['p_' strVarShort '_Conf=p_Var_Conf;']);
				eval(['vecR_' strVarShort '_Conf(intFile) = r_Var_Conf;']);
			end
			
			tbl = table(vecSpikeGroupNumOfCells,vecSpikeGroupAvgRateOfCells,vecSpikeGroupAvgTuningOfCells,vecSpikeGroupAvgRate,vecSpikeGroupConfidence,...
				'VariableNames',{'CellNumber','AvgCellRate','AvgCellTuning','AvgPopIFR','Conf'});
			mdl_predConf = fitlm(tbl,'linear');
			
			%predict IFR
			for intVar=1:intVarNum
				strVarShort = cellVars{intVar};
				strVarLong = cellVarsLong{intVar};
				strVarSource = ['vecSpikeGroup' strVarLong];
				vecSource = eval(strVarSource);
				
				%corrs
				[r_Var_IFR,p_Var_IFR]=corr(vecSource,vecSpikeGroupAvgRate);
				eval(['r_' strVarShort '_IFR=r_Var_IFR;']);
				eval(['p_' strVarShort '_IFR=p_Var_IFR;']);
				eval(['vecR_' strVarShort '_IFR(intFile) = r_Var_IFR;']);
			end
			
			tbl = table(vecSpikeGroupAvgRateOfCells,vecSpikeGroupAvgTuningOfCells,vecSpikeGroupAvgRate,...
				'VariableNames',{'AvgCellRate','AvgCellTuning','AvgPopIFR'});
			mdl_predIFR = fitlm(tbl,'linear');
			
			%rate vs tuning
			[rGroup_Rate_Tune,pGroup_Rate_Tune]=corr(vecTuningPerCell,vecRatePerCell);
			%}
			
			%% save
			%result: during high ifr epochs, cells with low firing rate and low tuning are also active,
			%while during low ifr epochs, only cells with high firing rate or high tuning are active
			
			%vecR_CombIFR_Rate(intFile) = table2array(mdl_predifr.Coefficients(2,1));
			%vecR_CombIFR_Tune(intFile) = table2array(mdl_predifr.Coefficients(3,1));
			
			%vecR_Group_Rate_Tune(intFile) = rGroup_Rate_Tune;
			
			cellQuantileDur(intFile,:) = cellValsDur;
			cellQuantileConf(intFile,:) = cellValsConf;
			matLatConf(intFile,:) = vecMeans;
			vecSpikeGroupSize(intFile) = intSpikeGroupSize;
			
			%% add quantile data
			intQuantileNum = 10;
			for intSuffixIdx=1:numel(cellSuffices)
				strSuffix = cellSuffices{intSuffixIdx}; %Dur/Conf/IFR
				strSuffixLong = cellSufficesLong{intSuffixIdx};
				vecSuffixSource = [sSpikeGroup.(strSuffixLong)]';
				for intVar=1:intVarNum
					strVarShort = cellVars{intVar};
					strVarLong = cellVarsLong{intVar};
					vecSource = [sSpikeGroup.(strVarLong)]';
					%strVarSource = ['vecSpikeGroup' strVarLong];
					%vecSource = eval(strVarSource);
					
					[~,~,~,~,~,cellValsSuff_Suff,cellValsVar_Suff]=getQuantiles(vecSuffixSource,vecSource,intQuantileNum);
					
					%assign output
					eval(['cellVals' strSuffix '_' strSuffix '_' strType ' = cellValsSuff_Suff;']);
					eval(['cellVals' strVarShort '_' strSuffix '_' strType ' = cellValsVar_Suff;']);
					eval(['cellQuantile' strSuffix '_' strSuffix '_' strType '(intFile,:) = cellValsSuff_Suff;']);
					eval(['cellQuantile' strVarShort '_' strSuffix '_' strType '(intFile,:) = cellValsVar_Suff;']);
				end
			end
		end
		hold(h2,'off');
		xlabel(h2,'Duration of n-spike block (ms)');
		ylabel(h2,'Decoder confidence');
		title(h2,sprintf('Deciles'));
		ylim(h2,[0 max(get(h2,'ylim'))]);
		
		plot(h1,vecBinsDurC,mean(matLatConf,1),'color',lines(1));%,vecSDs./sqrt(vecCounts),vecSDs./sqrt(vecCounts),'color',lines(1));
		hold(h1,'off')
		xlabel(h1,'Latency of n-spike block after stim onset (s)');
		ylabel(h1,'Decoder confidence');
		ylim(h1,[0 max(get(h1,'ylim'))]);
		
		
		%%
		matQuantDur = cellfun(@mean,cellQuantileDur);
		matQuantConf = cellfun(@mean,cellQuantileConf);
		matQuantConf = zscore(matQuantConf,[],2);
		matX = repmat(1:10,[intRecNum 1]);
		mdl = fitlm(matX(:),matQuantConf(:));
		r=mdl.Coefficients.Estimate(2);
		p=mdl.Coefficients.pValue(2);
		
		hold(h3,'on')
		plot(h3,matX',matQuantConf','color',vecColGrey);
		plot(h3,mean(matX,1),mean(matQuantConf,1),'color',lines(1));
		hold(h3,'off')
		xlabel(h3,'Duration decile of n-spike block');
		ylabel(h3,sprintf('Decoder confidence, z-scored per rec (%s)',getGreek('sigma')));
		%ylim(h3,[0 max(get(h3,'ylim'))]);
		title(h3,sprintf('Lin reg, r=%.3f %s/decile, p=%.3e',r,getGreek('sigma'),p));
		fixfig;
		
		%% plot corrs
		%test vs 0
		matP = nan(3,intVarNum);
		for i=1:numel(cellSuffices)
			%transform variable names and calculate p
			strSuffix = cellSuffices{i};
			for intVar=1:intVarNum
				strVar = cellVars{intVar};
				if strcmp(strVar,strSuffix)
					eval(['vecR_' strVar ' = nan(size(vecR_' cellVars{1} '_' strSuffix '));']);
				else
					eval(['[h,matP(i,intVar)] = ttest(vecR_' strVar '_' strSuffix ');']);
					eval(['vecR_' strVar ' = vecR_' strVar '_' strSuffix ';']);
				end
			end
			vecP=matP(i,:);
			[h, crit_p, adj_p]=fdr_bh(vecP,0.05);
			
			%% plot
			subplot(2,3,3+i);cla
			hold on;
			plot([0 intVarNum*2-1],[0 0],'--','color',[0.5 0.5 0.5]);
			matCol = lines(intVarNum);
			matCol = matCol([2 1 3:size(matCol,1)],:);
			matCol(3,:) = 0;
			matColP = 1-((1-matCol)./2);
			
			intRecNum = numel(vecR_Num);
			vecX = ones(size(vecR_Num));
			strTitle = 'p: ';
			for intVar=1:intVarNum
				strVar = cellVars{intVar};
				intPlotLoc = intVar*2-1;
				vecSource = eval(['vecR_' strVar]);
				swarmchart(intPlotLoc*vecX,vecSource,[],matColP(intVar,:),'filled');
				errorbar(intPlotLoc,mean(vecSource),std(vecSource)./sqrt(intRecNum),'x','color',matCol(intVar,:));
				strTitle = [strTitle sprintf('%s=%.2e,',cellVarsShort{intVar},adj_p(intVar))];
				if intVar==round(intVarNum/2)
					strTitle = [strTitle sprintf('\n')];
				end
			end
			hold off
			set(gca,'xtick',[1:2:(intVarNum*2)],'xticklabel',cellVarsLegend);
			xtickangle(gca,45);
			ylabel(['Correlation with group ' strSuffix ' (r)']);
			title(strTitle)
		end
		fixfig;
		
		%%
		if boolSaveFig
			drawnow;
			export_fig(fullpath(strFigurePath,sprintf('M1_SpikeBlockDecoding%s%s%s%s%s.tif',strRunType,strType,strSGS,strConnType,strOnset)));
			export_fig(fullpath(strFigurePath,sprintf('M1_SpikeBlockDecoding%s%s%s%s%s.pdf',strRunType,strType,strSGS,strConnType,strOnset)));
		end
		
		%% plot deciles per var
		%test vs 0
		matP = nan(3,numel(cellVars));
		for intSuffixIdx=1:numel(cellSuffices)
			%prep plot
			figure;maxfig;
			hSummary=subplot(2,6,12);
			cla(hSummary);
			hold(hSummary,'on');
			plot(hSummary,[0 (2*numel(cellVars))-1],[0 0],'--','color',[0.5 0.5 0.5]);
			title(hSummary,sprintf('%s-%s, Bonferroni-corrected p & CI',strRunType,strType));
			
			%transform variable names and calculate p
			strSuffix = cellSuffices{intSuffixIdx};
			for intVar=1:numel(cellVars)
				strVar = cellVars{intVar};
				eval(['cellY = cellQuantile' strVar '_' strSuffix '_' strType ';']);
				
				
				subplot(2,6,intVar)
				matQuantY = cellfun(@mean,cellY);
				cellSummaryMatrices{intType,intVar,intSuffixIdx,intConnType} = matQuantY;
				
				matQuantY_Z = zscore(matQuantY,[],2);
				matX = repmat(1:10,[intRecNum 1]);
				mdl = fitlm(matX(:),matQuantY_Z(:));
				r=mdl.Coefficients.Estimate(2);
				r_SE=mdl.Coefficients.SE(2);
				dblAlpha = 0.05/numel(cellVars);
				r_CI = coefCI(mdl,dblAlpha);
				r_CI = r_CI(2,:);
				p=mdl.Coefficients.pValue(2)*5;
				
				hold('on')
				plot(matX',matQuantY_Z','color',vecColGrey);
				plot(mean(matX,1),mean(matQuantY_Z,1),'color',matCol(intVar,:));
				hold('off')
				xlabel([strSuffix ' decile of n-spike block']);
				ylabel(sprintf('%s, z-scored per rec (%s)',strVar,getGreek('sigma')));
				title(sprintf('OLS, y=%.2f %s/x, p=%.2e',r,getGreek('sigma'),p));
				
				%plot in summary
				errorbar(hSummary,(intVar*2-1),r,r-r_CI(1),r-r_CI(2),'x','color',matCol(intVar,:));
			end
			hold(hSummary,'off');
			set(hSummary,'xtick',[1:2:(2*numel(cellVars))],'xticklabel',cellVars);
			xtickangle(hSummary,45);
			ylabel(hSummary,sprintf('Lin reg slope, mean +/- 95 CI (%s/decile)',getGreek('sigma')));
			ylim(hSummary,[-0.4 0.4]);
			fixfig;
			
			if boolSaveFig
				%%
				export_fig(fullpath(strFigurePath,sprintf('M1_SpikeBlockCorr%s%s%s%s%s%s.tif',strType,strRunType,strSuffix,strSGS,strConnType,strOnset)));
				export_fig(fullpath(strFigurePath,sprintf('M1_SpikeBlockCorr%s%s%s%s%s%s.pdf',strType,strRunType,strSuffix,strSGS,strConnType,strOnset)));
			end
		end
	end
end

%% plot real and shufftid normalized to uniform
intPlotType = 3; %dur/conf/ifr => ifr is v1, others are only within v2
for intConnType=1:numel(cellConnType)
	strConnType = cellConnType{intConnType};
	
	%test vs 0
	intTestVarNum = 6;
	matP = nan(2,numel(cellVars));
	figure;maxfig;
	for i=1:numel(cellSuffices) %dur, conf, ifr
		%prep plot
		strSuffix = cellSuffices{i};
		hSummary=[1 2];
		cellStrP = {'',''};
		for intRealShuff=1:2
			if intRealShuff == 1
				hSummary(intRealShuff)=subplot(2,3,i);
			else
				hSummary(intRealShuff)=subplot(2,3,i+3);
			end
			cla(hSummary(intRealShuff));
			hold(hSummary(intRealShuff),'on');
			plot(hSummary(intRealShuff),[0 (2*numel(cellVars))-1],[0 0],'--','color',[0.5 0.5 0.5]);
		end
		
		%transform variable names and calculate p
		for intVar=1:numel(cellVars)
			strVar = cellVars{intVar};
			matQuantY_Real = cellSummaryMatrices{1,intVar,i,intConnType};
			matQuantY_ShuffTid = cellSummaryMatrices{2,intVar,i,intConnType};
			matQuantY_Uniform = cellSummaryMatrices{3,intVar,i,intConnType};
			
			for intRealShuff=1:2
				if intRealShuff == 1
					matRelY = matQuantY_Real ./ matQuantY_ShuffTid;
					strPlot = 'x';
				else
					matRelY = matQuantY_Uniform ./ matQuantY_ShuffTid;
					strPlot = 'o--';
				end
				matQuantY_Z = zscore(matRelY,[],2);
				matX = repmat(1:10,[intRecNum 1]);
				mdl = fitlm(matX(:),matQuantY_Z(:));
				r=mdl.Coefficients.Estimate(2);
				r_SE=mdl.Coefficients.SE(2);
				dblAlpha = 0.05/intTestVarNum;
				r_CI = coefCI(mdl,dblAlpha);
				r_CI = r_CI(2,:);
				p=mdl.Coefficients.pValue(2)*numel(cellVars);
				matP(intRealShuff,intVar) = p;
				cellStrP{intRealShuff} = [cellStrP{intRealShuff} sprintf('%s=%.2f',cellVars{intVar}(1),p)];
				%plot in summary
				errorbar(hSummary(intRealShuff),(intVar*2-1),r,r-r_CI(1),r-r_CI(2),strPlot,'color',matCol(intVar,:));
				
				%save
				if i == intPlotType && intRealShuff==1
					matR(intVar,:,intConnType) = [r r_CI r_SE];
				end
			end
		end
		%finish
		cellRealShuff = {'Real','Shuff'};
		for intRealShuff=1:2
			hold(hSummary(intRealShuff),'off');
			set(hSummary(intRealShuff),'xtick',[1:2:(2*numel(cellVars))],'xticklabel',cellVars);
			xtickangle(hSummary(intRealShuff),45);
			ylabel(hSummary(intRealShuff),sprintf('ShuffTid-norm. r, mean +/- 95 CI (%s/decile)',getGreek('sigma')));
			ylim(hSummary(intRealShuff),[-0.4 0.4]);
			strP = cellStrP{intRealShuff};
			strTitleType = [strSuffix ', ' cellRealShuff{intRealShuff} ', ' strRunType ', ' strSGS ', O' strOnset];
			title(hSummary(intRealShuff),sprintf('%s, Bonf-corr p & CI\n%s',strTitleType,strP),'interpreter','none');
		end
		fixfig;
		
		if boolSaveFig
			%%
			export_fig(fullpath(strFigurePath,sprintf('M1_SummaryShuffTidNormalized%s%s%s%s.tif',strRunType,strSGS,strConnType,strOnset)));
			export_fig(fullpath(strFigurePath,sprintf('M1_SummaryShuffTidNormalized%s%s%s%s.pdf',strRunType,strSGS,strConnType,strOnset)));
		end
	end
	
end

% plot most/least comparison
matLeastR = matR(:,:,1); %[r rci_lo r_ci_hi r_se]
matMostR = matR(:,:,2); %[r rci_lo r_ci_hi r_se]

%p-value
vecP = nan(1,intVarNum);
figure;hold on
plot([1 intVar*2+1],[0 0],'--','color',[0.5 0.5 0.5]);
for intVar=1:intVarNum
	vecLeast = matLeastR(intVar,:);
	vecMost = matMostR(intVar,:);
	
	mDiff = vecLeast(1)-vecMost(1);
	sDiff = sqrt(0.5*((range(vecLeast(2:3))/2).^2+(range(vecMost(2:3))/2).^2));
	zDiff = abs(mDiff/sDiff);
	pDiff = normcdf(zDiff,'upper')*2;
	
	vecP(intVar) = pDiff;
	
	vecM = [vecMost(1) vecLeast(1) ];
	errorbar(intVar*2 + [-0.5 0.5],vecM,vecM-[vecMost(2) vecLeast(2) ],vecM-[vecMost(3) vecLeast(3) ],'x-','color',matCol(intVar,:));
end
%vecP_bh = bonf_holm(vecP);%CIs are already corrected
strTitle = '';
for intVar=1:intVarNum
	strTitle = [strTitle sprintf(';%s=%.2e',cellVars{intVar}(1),vecP(intVar))];
	if intVar == 4
		strTitle = [strTitle newline];
	end
end

vecX = 2:2:(intVar*2);
set(gca,'xtick',vecX,'xticklabel',cellVars);
title(strTitle);
ylim([ -0.4 0.4]);
ylabel(sprintf('ShuffTid-norm. r, mean +/- 95 CI (%s/decile)',getGreek('sigma')));
fixfig

if boolSaveFig
	%%
	export_fig(fullpath(strFigurePath,sprintf('M1_SummaryMostVsLeast%s%s%s.tif',strRunType,strSGS,strOnset)));
	export_fig(fullpath(strFigurePath,sprintf('M1_SummaryMostVsLeast%s%s%s.pdf',strRunType,strSGS,strOnset)));
end