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
intGroupOrDecile = 1;
intOnsetType = 0; %normal (0) or rem onset (1)

%% define qualifying areas
boolSaveFig = true;
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
cellSuffices = {'Dur','Conf','IFR'};
cellSufficesLong = {'Duration','Confidence','AvgRate'};
cellVars = {'Num','Rate','Tune','IFR','Conf','FrIn','PrDi','Bndw','FrSupra','FrGranu','FrInfra'};
cellVarsShort = {'#','R','T','Ifr','C','%I','PrD','B','%S','%G','%I'};
cellVarsLong = {'NumOfCells','AvgRateOfCells','AvgTuningOfCells','AvgRate','Confidence','FractionInterneurons',...
	'AvgPrefDistToStim','AvgBandwidthOfCells','FractionSupra','FractionGranu','FractionInfra'};
cellVarsLegend = {'# of cells','Cell Rate','Cell Tune','Pop Rate','Conf','% InterN',...
	'Pref-d','BndW','% Supra','% Granu','% Infra'};

%remove layers if model
if intRunType == 2
	cellVars = cellVars(1:8);
	cellVarsShort = cellVarsShort(1:numel(cellVars));
	cellVarsLong = cellVarsLong(1:numel(cellVars));
	cellVarsLegend = cellVarsLegend(1:numel(cellVars));
end

intVarNum = numel(cellVars);
for intType=1:numel(cellTypes)
	%%
	strType = cellTypes{intType};
	if boolFixedSpikeGroupSize
		sFiles = dir ([strTargetDataPath 'Q1cData_' strRunType '*' strType '*Fixed*.mat']);
		strSGS = 'FixedSGS';
	else
		sFiles = dir ([strTargetDataPath 'Q1cData_' strRunType '*' strType '*Var*.mat']);
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
	
	%plot(h2,[0 100],([1 1]./intStimNum),'--','color',vecColChance);
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
		%rSpike_IFR_Rate = sLoad.rSpike_IFR_Rate;
		%pSpike_IFR_Rate = sLoad.pSpike_IFR_Rate;
		%rSpike_IFR_Tune = sLoad.rSpike_IFR_Tune;
		%pSpike_IFR_Tune = sLoad.pSpike_IFR_Tune;
		vecTuningPerCell = sLoad.vecTuningPerCell;
		vecRatePerCell = sLoad.vecRatePerCell;
		
		%delete outliers
		[vecSortedRate,vecSort]=sort(intSpikeGroupSize./[sSpikeGroup.Duration]');
		vecRem = vecSort([1:round(numel(vecSortedRate)*0.025) round(numel(vecSortedRate)*0.975):numel(vecSortedRate)]);
		sSpikeGroup(vecRem) = [];
		
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
		if isfield(sSpikeGroup,'FractionSupra')
			vecSpikeGroupFractionSupra = [sSpikeGroup.FractionSupra]';
			vecSpikeGroupFractionGranu = [sSpikeGroup.FractionGranu]';
			vecSpikeGroupFractionInfra = [sSpikeGroup.FractionInfra]';
		else
			vecSpikeGroupFractionSupra = nan*[sSpikeGroup.AvgBandwidthOfCells]';
			vecSpikeGroupFractionGranu = nan*[sSpikeGroup.AvgBandwidthOfCells]';
			vecSpikeGroupFractionInfra = nan*[sSpikeGroup.AvgBandwidthOfCells]';
		end
		%sort by duration
		[vecSortedDur,vecSort]=sort(vecSpikeGroupDuration);
		vecSortedCorr = vecSpikeGroupCorrect(vecSort);
		vecSortedConf = vecSpikeGroupConfidence(vecSort);
		vecSortedLat = vecSpikeGroupLatency(vecSort);
		vecSortedNum = vecSpikeGroupNumOfCells(vecSort);
		vecSortedRate = vecSpikeGroupAvgRateOfCells(vecSort);
		vecSortedTune = vecSpikeGroupAvgTuningOfCells(vecSort);
		vecSortedIFR = vecSpikeGroupAvgRate(vecSort);
		
		%calculate fraction correct and confidence per bin of equal size
		intQuantileNum = 10;
		%[vecMeanDur,vecSemDur,vecMeanConf,vecSemConf,vecQuantile,cellValsDur,cellValsConf]=getQuantiles(vecSortedDur,vecSortedConf,intQuantileNum);
		[vecMeanIfr,vecSemIfr,vecMeanTune,vecSemCTune,vecQuantile,cellValsIfr,cellValsTune]=getQuantiles(vecSortedIFR,vecSortedTune-mean(vecSortedTune),intQuantileNum);
		
		%conf with dur
		%errorbar(h2,vecMeanDur*1000,vecMeanConf,vecSemConf,vecSemConf,vecSemDur,vecSemDur,'color',lines(1));
		errorbar(h2,1:intQuantileNum,vecMeanTune,vecSemCTune,'color',lines(1));
		
		%conf with time
		[vecCounts,vecMeans,vecSDs]=makeBins(vecSpikeGroupLatency,vecSpikeGroupConfidence,vecBinsDur);
		plot(h1,vecBinsDurC,vecMeans,'color',vecColGrey);%,vecSDs./sqrt(vecCounts),vecSDs./sqrt(vecCounts),'color',lines(1));
		pout=binotest(sum(vecSpikeGroupCorrect),numel(vecSpikeGroupCorrect),1/numel(unique(vecOri180)));
		%if mean(vecSpikeGroupCorrect)<0.3,continue;end
		
		%% group or decile
		if intGroupOrDecile == 2
			%% deciles of dur
			for intVar=1:intVarNum
				strVarShort = cellVars{intVar};
				strVarLong = cellVarsLong{intVar};
				strVarMean = ['vecMean' strVarShort];
				strVarSem = ['vecSem' strVarShort];
				
				strVarSource = ['vecSpikeGroup' strVarLong];
				vecSource = eval(strVarSource);
				
				[vecMeanDurDeciles,vecSemDurDeciles,vecMeanOut,vecSemOut]=getQuantiles(vecSpikeGroupDuration,vecSource,intQuantileNum);
				vecMeanDurDeciles = vecMeanDurDeciles';
				
				%assign output
				eval([strVarMean '=vecMeanOut'';']);
				eval([strVarSem '=vecSemOut'';']);
				
				%corrs
				[r_Var_Dur,p_Var_Dur]=corr(vecMeanOut,vecMeanDurDeciles);
				eval(['r_' strVarShort '_Dur=r_Var_Dur;']);
				eval(['p_' strVarShort '_Dur=p_Var_Dur;']);
				
			end
			
			%prediction with dur
			tbl = table(vecMeanNum,vecMeanRate,vecMeanTune,vecMeanIFR,vecMeanTune,vecMeanDurDeciles,...
				'VariableNames',{'CellNumber','AvgCellRate','AvgCellTuning','AvgPopIFR','Conf','Dur'});
			mdl_predDur = fitlm(tbl,'linear');
			
			%% deciles of conf
			for intVar=1:intVarNum
				strVarShort = cellVars{intVar};
				strVarLong = cellVarsLong{intVar};
				strVarMean = ['vecMean' strVarShort];
				strVarSem = ['vecSem' strVarShort];
				
				strVarSource = ['vecSpikeGroup' strVarLong];
				vecSource = eval(strVarSource);
				
				[vecMeanConfDeciles,vecSemConfDeciles,vecMeanOut,vecSemOut]=getQuantiles(vecSpikeGroupConfidence,vecSource,intQuantileNum);
				vecMeanConfDeciles = vecMeanConfDeciles';
				
				%assign output
				eval([strVarMean '=vecMeanOut'';']);
				eval([strVarSem '=vecSemOut'';']);
				
				%corrs
				[r_Var_Conf,p_Var_Conf]=corr(vecMeanOut,vecMeanConfDeciles);
				eval(['r_' strVarShort '_Conf=r_Var_Conf;']);
				eval(['p_' strVarShort '_Conf=p_Var_Conf;']);
			end
			
			tbl = table(vecMeanNum,vecMeanRate,vecMeanTune,vecMeanIFR,vecMeanConfDeciles,...
				'VariableNames',{'CellNumber','AvgCellRate','AvgCellTuning','AvgPopIFR','Conf'});
			mdl_predConf = fitlm(tbl,'linear');
			
			%% deciles of ifr
			for intVar=1:intVarNum
				strVarShort = cellVars{intVar};
				strVarLong = cellVarsLong{intVar};
				strVarMean = ['vecMean' strVarShort];
				strVarSem = ['vecSem' strVarShort];
				
				strVarSource = ['vecSpikeGroup' strVarLong];
				vecSource = eval(strVarSource);
				
				[vecMeanIFRDeciles,vecSemIFRDeciles,vecMeanOut,vecSemOut]=getQuantiles(vecSpikeGroupAvgRate,vecSource,intQuantileNum);
				vecMeanIFRDeciles = vecMeanIFRDeciles';
				
				%assign output
				eval([strVarMean '=vecMeanOut'';']);
				eval([strVarSem '=vecSemOut'';']);
				
				%corrs
				[r_Var_IFR,p_Var_IFR]=corr(vecMeanOut,vecMeanIFRDeciles);
				eval(['r_' strVarShort '_IFR=r_Var_IFR;']);
				eval(['p_' strVarShort '_IFR=p_Var_IFR;']);
			end
			
			%predict ifr
			[r_Num_IFR,p_Num_IFR]=corr(vecMeanNum,vecMeanIFRDeciles);
			[r_Rate_IFR,p_Rate_IFR]=corr(vecMeanRate,vecMeanIFRDeciles);
			[r_Tune_IFR,p_Tune_IFR]=corr(vecMeanTune,vecMeanIFRDeciles);
			tbl = table(vecMeanRate,vecMeanTune,vecMeanIFRDeciles,...
				'VariableNames',{'AvgCellRate','AvgCellTuning','AvgPopIFR'});
			mdl_predIFR = fitlm(tbl,'linear');
			
			%rate vs tuning
			[rGroup_Rate_Tune,pGroup_Rate_Tune]=corr(vecTuningPerCell,vecRatePerCell);
			
		else
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
		end
		%% save
		%result: during high ifr epochs, cells with low firing rate and low tuning are also active,
		%while during low ifr epochs, only cells with high firing rate or high tuning are active
		
		%vecR_CombIFR_Rate(intFile) = table2array(mdl_predifr.Coefficients(2,1));
		%vecR_CombIFR_Tune(intFile) = table2array(mdl_predifr.Coefficients(3,1));
		
		%vecR_Group_Rate_Tune(intFile) = rGroup_Rate_Tune;
		
		cellQuantileDur(intFile,:) = cellValsIfr;
		cellQuantileConf(intFile,:) = cellValsTune;
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
				eval(['cellVals' strSuffix '_' strSuffix ' = cellValsSuff_Suff;']);
				eval(['cellVals' strVarShort '_' strSuffix ' = cellValsVar_Suff;']);
				eval(['cellQuantile' strSuffix '_' strSuffix '(intFile,:) = cellValsSuff_Suff;']);
				eval(['cellQuantile' strVarShort '_' strSuffix '(intFile,:) = cellValsVar_Suff;']);
			end
		end
	end
	hold(h2,'off');
	%xlabel(h2,'Duration of n-spike block (ms)');
	%ylabel(h2,'Decoder confidence');
	xlabel(h2,'Avg pop rate of n-spike block (decile)');
	ylabel(h2,sprintf('Avg cell tuning (%st-stat)',getGreek('Delta')));
	title(h2,sprintf('Deciles'));
	%ylim(h2,[0 max(get(h2,'ylim'))]);
	
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
			matR(intType,i,intVar,:) = vecSource; %[Real/Shuff/Unif x Dur/Conf/IFR x Var x Rec]
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
		export_fig(fullpath(strFigurePath,sprintf('Q1c_SpikeBlockDecoding%s%s%s%s.tif',strRunType,strType,strSGS,strOnset)));
		export_fig(fullpath(strFigurePath,sprintf('Q1c_SpikeBlockDecoding%s%s%s%s.pdf',strRunType,strType,strSGS,strOnset)));
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
			eval(['cellY = cellQuantile' strVar '_' strSuffix ';']);
			
			
			subplot(2,6,intVar)
			matQuantY = cellfun(@mean,cellY);
			cellSummaryMatrices{intType,intVar,intSuffixIdx} = matQuantY;
			
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
			export_fig(fullpath(strFigurePath,sprintf('Q1c_SpikeBlockCorr%s%s%s%s%s.tif',strType,strRunType,strSuffix,strSGS,strOnset)));
			export_fig(fullpath(strFigurePath,sprintf('Q1c_SpikeBlockCorr%s%s%s%s%s.pdf',strType,strRunType,strSuffix,strSGS,strOnset)));
		end
	end
end

%% plot real and shufftid normalized to uniform
%test vs 0
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
		eval(['cellY = cellQuantile' strVar '_' strSuffix ';']);
		
		
		matQuantY_Real = cellSummaryMatrices{1,intVar,i};
		matQuantY_ShuffTid = cellSummaryMatrices{2,intVar,i};
		matQuantY_Uniform = cellSummaryMatrices{3,intVar,i};
		
		for intRealShuff=1:2
			if intRealShuff == 1
				matRelY = matQuantY_Real ./ matQuantY_Uniform;
				strPlot = 'x';
			else
				matRelY = matQuantY_ShuffTid ./ matQuantY_Uniform;
				strPlot = 'o--';
			end
			matQuantY_Z = zscore(matRelY,[],2);
			matX = repmat(1:10,[intRecNum 1]);
			mdl = fitlm(matX(:),matQuantY_Z(:));
			r=mdl.Coefficients.Estimate(2);
			r_SE=mdl.Coefficients.SE(2);
			dblAlpha = 0.05/numel(cellVars);
			r_CI = coefCI(mdl,dblAlpha);
			r_CI = r_CI(2,:);
			p=mdl.Coefficients.pValue(2)*numel(cellVars);
			matP(intRealShuff,intVar) = p;
			cellStrP{intRealShuff} = [cellStrP{intRealShuff} sprintf('%s=%.2f',cellVars{intVar}(1),p)];
			%plot in summary
			errorbar(hSummary(intRealShuff),(intVar*2-1),r,r-r_CI(1),r-r_CI(2),strPlot,'color',matCol(intVar,:));
		end
	end
	%finish
	cellRealShuff = {'Real','Shuff'};
	for intRealShuff=1:2
		hold(hSummary(intRealShuff),'off');
		set(hSummary(intRealShuff),'xtick',[1:2:(2*numel(cellVars))],'xticklabel',cellVars);
		xtickangle(hSummary(intRealShuff),45);
		ylabel(hSummary(intRealShuff),sprintf('Uniform-norm. r, mean +/- 95 CI (%s/decile)',getGreek('sigma')));
		ylim(hSummary(intRealShuff),[-0.4 0.4]);
		strP = cellStrP{intRealShuff};
		strTitleType = [strSuffix ', ' cellRealShuff{intRealShuff} ', ' strRunType ', ' strSGS ', O' strOnset];
		title(hSummary(intRealShuff),sprintf('%s, Bonf-corr p & CI\n%s',strTitleType,strP),'interpreter','none');
	end
	fixfig;
	
	if boolSaveFig
		%%
		export_fig(fullpath(strFigurePath,sprintf('Q1c_SummaryUniformNormalized%s%s%s.tif',strRunType,strSGS,strOnset)));
		export_fig(fullpath(strFigurePath,sprintf('Q1c_SummaryUniformNormalized%s%s%s.pdf',strRunType,strSGS,strOnset)));
	end
end

%% plot real and shufftid normalized to uniform
%test vs 0
intVarNum = 11;
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
		eval(['cellY = cellQuantile' strVar '_' strSuffix ';']);
		
		
		matQuantY_Real = cellSummaryMatrices{1,intVar,i};
		matQuantY_ShuffTid = cellSummaryMatrices{2,intVar,i};
		matQuantY_Uniform = cellSummaryMatrices{3,intVar,i};
		
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
			dblAlpha = 0.05/intVarNum;
			r_CI = coefCI(mdl,dblAlpha);
			r_CI = r_CI(2,:);
			p=mdl.Coefficients.pValue(2)*intVarNum;
			matP(intRealShuff,intVar) = p;
			cellStrP{intRealShuff} = [cellStrP{intRealShuff} sprintf('%s=%.2e',cellVars{intVar}(1),p)];
			%plot in summary
			errorbar(hSummary(intRealShuff),(intVar*2-1),r,r-r_CI(1),r-r_CI(2),strPlot,'color',matCol(intVar,:));
		end
	end
	%finish
	cellRealShuff = {'Real','Uniform'};
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
		export_fig(fullpath(strFigurePath,sprintf('Q1c_SummaryShuffTidNormalized%s%s%s.tif',strRunType,strSGS,strOnset)));
		export_fig(fullpath(strFigurePath,sprintf('Q1c_SummaryShuffTidNormalized%s%s%s.pdf',strRunType,strSGS,strOnset)));
	end
end

%% plot real and shufftid pairwise
%test vs 0
vecP = nan(1,numel(cellVars));
figure;maxfig;
for i=1:numel(cellSuffices) %dur, conf, ifr
	%prep plot
	h=subplot(2,3,i);
	strSuffix = cellSuffices{i};
	hold('on');
	plot([0 (numel(cellVars))-1],[0 0],'--','color',[0.5 0.5 0.5]);
	
	%transform variable names and calculate p
	for intVar=1:numel(cellVars)
		strVar = cellVars{intVar};
		vecR_Real = squeeze(matR(1,i,intVar,:));
		vecR_Shuff = squeeze(matR(2,i,intVar,:));
		vecR_Unif = squeeze(matR(3,i,intVar,:));
		indRem = vecR_Real==0;
		vecR_Real(indRem) = [];
		vecR_Shuff(indRem) = [];
		vecR_Unif(indRem) = [];
		intRecNum = numel(vecR_Real);
		matX = repmat(intVar+[-0.2 0.2],[numel(vecR_Unif) 1]);
		plot(matX', [vecR_Real vecR_Shuff]','color', [0.5 0.5 0.5]);
		
		errorbar(intVar+[-0.2 0.2],mean([vecR_Real vecR_Shuff]),std([vecR_Real vecR_Shuff])./sqrt(intRecNum),'color',matCol(intVar,:));
		
	end
	
	%plot diff
	subplot(2,3,i+3);hold on
	plot([0 (2*numel(cellVars))-1],[0 0],'--','color',[0.5 0.5 0.5]);
	
	
	%transform variable names and calculate p
	vecP = nan(1,numel(cellVars));
	for intVar=1:numel(cellVars)
		strVar = cellVars{intVar};
		vecR_Real = squeeze(matR(1,i,intVar,:));
		vecR_Shuff = squeeze(matR(2,i,intVar,:));
		vecR_Unif = squeeze(matR(3,i,intVar,:));
		indRem = vecR_Real==0;
		vecR_Real(indRem) = [];
		vecR_Shuff(indRem) = [];
		vecR_Unif(indRem) = [];
		intRecNum = numel(vecR_Real);
		
		vecDiff = vecR_Real - vecR_Shuff;
		swarmchart(repmat(intVar*2,[intRecNum 1]),vecDiff,[],matColP(intVar,:),'filled');
		errorbar(intVar*2,mean(vecDiff),std(vecDiff)./sqrt(intRecNum),'x','color',matCol(intVar,:));
		[h,p]=ttest(vecDiff);
		vecP(intVar) = p;
	end
	
	hold('off');
	set(gca,'xtick',[2:2:(2*numel(cellVars))],'xticklabel',cellVars);
	xtickangle(gca,45);
	ylabel(gca,sprintf('ShuffTid-norm. %sr, mean +/- sem',getGreek('delta')));
	ylim(gca,[-0.4 0.4]);
	strTitleType = [strSuffix ', ' strRunType ', ' strSGS ', O' strOnset];
	
	indUsePs = ~isnan(vecP);
	vecP_corr2 = vecP;
	vecP_corr = vecP;
	[d,d2,vecPs_corr] = fdr_bh(vecP(indUsePs));
	[vecPs_corr2] = bonf_holm(vecP(indUsePs));
	vecP_corr(indUsePs) = vecPs_corr;
	vecP_corr2(indUsePs) = vecPs_corr2;
	
	strP= '';
	for intVar=1:numel(cellVars)
		strP = [strP sprintf('%s=%.2e',cellVars{intVar}(1),vecP_corr2(intVar))];
	end
	title(gca,sprintf('%s, Bonf-holm corr p\n%s',strTitleType,strP),'interpreter','none');
	ylim([ -0.5 0.5]);
end
fixfig;
%
if boolSaveFig
	%%
	export_fig(fullpath(strFigurePath,sprintf('Q1c_SummaryPairwise%s%s%s.tif',strRunType,strSGS,strOnset)));
	export_fig(fullpath(strFigurePath,sprintf('Q1c_SummaryPairwise%s%s%s.pdf',strRunType,strSGS,strOnset)));
end