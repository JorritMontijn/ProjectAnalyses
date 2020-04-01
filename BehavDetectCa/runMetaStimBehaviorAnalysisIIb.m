%{
figures to save and take mean over animals
- cellSavePreStim: dF/F detect/no-detect preceding stimulus frames [5 (c) x 2 (d/nd) cell array] x n (animals) with every cell = vector of values
- cellSaveBehavDetect: behavioral detection [5 (c) x n (animals)] cell array with every cell = vector of detect/no detect values (1/0)
- cellSaveBehavRT: behavioral RT [5 (c) x n (animals)] cell array with every cell = vector of RTs
- cellSaveBaseSubtrStimAct: dF/F detect/no-detect baseline-subtracted activity [5 (c) x 2 (d/nd) cell array] x n (animals) with every cell = vector of values
- cellSaveDuringStim: dF/F detect/no-detect during stimulus frames [5 (c) x 2 (d/nd) cell array] x n (animals) with every cell = vector of values
- cellSaveTraceAct: dF/F trace detect/no-detect pref/non-pref [2 (pref/no-pref) x 2 (detect/no-detect) x 6 (c)] x n (animals) with every cell = [Nreps x Nframes] matrix of trace dF/F values
-  AD: cellSaveMatrices: blob matrices in [1 x 4 (detect/no-detect/diff/significance)] cell array with [5 (c) x 31 (bins)] matrix
-  AD: cellSaveHCARNeuronConsistency: [1 x i (neuronal population)] cell array with vector of shuffled correlation values, where first value is actual data
-  AD: cellSaveSignalCorrs: signal correlations [3 (within-HCAR, between, within-non-HCAR) x n (animals/blocks)] with every cell = vector of pairwise correlation values
-  AD: cellSaveNoiseCorrs: noise correlations [3 (within-HCAR, between, within-non-HCAR) x n (animals/blocks)] with every cell = vector of pairwise correlation values
-  AD: cellSaveHCAR: normalized increase in dF/F for detection trials [5 (c) x n (animals/blocks)] cell array with each cell = vector of neurons with normalized dF/F increase for detection trials
-  AD: cellSaveSignalCorrsBiDir: signal correlations [3 (within-HCAR, between, within-non-HCAR) x n (animals/blocks)] with every cell = vector of pairwise correlation values
-  AD: cellSaveNoiseCorrsBiDir: noise correlations [3 (within-HCAR, between, within-non-HCAR) x n (animals/blocks)] with every cell = vector of pairwise correlation values
-  AD: cellSaveNormActDissim: within-group z-scored activation dissimilarity [2 (high/low HCAR) x 2 (detect/no-detect)] x n (animals/blocks) with every cell = vector of dissimilarity values
-  AD: cellSaveDissimCorrITD: inter-trial-distance dependence of assembly consistency [2 (high/low HCAR) x 2 (hits/misses) x n (animals/blocks)] with every cell = [p x 2] matrix with intertrial distance and correlation values for each trial pair
	- extra graph of above with means comparison of hits/misses high/low HCAR

%}
close all
clear all
%% parameters
strFigDir = 'D:\Data\ResultsStimDetectionCa\_meta';
strDataDir = 'D:\Data\ResultsStimDetectionCa';
cellInclude = {'20140207','20140314','20140425','20140507','20140530','20140604','20140711','20140715'};%exclude:20140430
%cellInclude = {'20140207','20140507','20140711','20140715'};%exclude:20140430
boolSavePlots = true;
boolNPS = true; 
intFigCounter = 0;
if boolNPS
	cellInclude = cellfun(@strcat,cellInclude,cellfill('NPS',size(cellInclude)),'UniformOutput',false);
	strFigDir = [strFigDir 'NPS'];
end

%% pre-allocate
intAnimal = 0;
vecAnimalSource = [];
vecSubPopulationSource = [];
cellPreStim = {};%: dF/F detect/no-detect preceding stimulus frames [5 (c) x 2 (d/nd) cell array] x n (animals) with every cell = vector of values
cellBehavDetect = {};%: behavioral detection [5 (c) x n (animals)] cell array with every cell = vector of detect/no detect values (1/0)
cellBehavRT = {};%: behavioral RT [5 (c) x n (animals)] cell array with every cell = vector of RTs
cellBaseSubtrStimAct = {};%: dF/F detect/no-detect baseline-subtracted activity [5 (c) x 2 (d/nd) cell array] x n (animals) with every cell = vector of values
cellDuringStim = {};%: dF/F detect/no-detect during stimulus frames [5 (c) x 2 (d/nd) cell array] x n (animals) with every cell = vector of values
cellTraceAct = {};%: dF/F trace detect/no-detect pref/non-pref [2 (pref/no-pref) x 2 (detect/no-detect) x 6 (c)] x n (animals) with every cell = [Nreps x Nframes] matrix of trace dF/F values

matActCohensD = [];
matActPPCohensD = [];
matActNPPCohensD = [];
matHetCohensD = [];
matHetPPCohensD = [];
matHetNPPCohensD = [];

matRelActValsR2 = [];
matRelActSDR2 = [];
matCorrs = [];
cellNormActDissim = {};%: within-group z-scored activation dissimilarity [2 (high/low HCAR) x 2 (detect/no-detect)] x n (animals/blocks) with every cell = vector of dissimilarity values
matITC = [];%: inter-trial-distance dependence of assembly consistency [2 (hit/miss) x 3 (ITD/HCAR/non-HCAR) x n (animals/blocks)] with every cell = vector with correlation values
cellRTDependency = {};%: RT dependence of activity dissimilarity [3 (RT/HCAR/non-HCAR) x n (animals/blocks)] with every cell = vector of trials
matMixBootstrappedDecodingOutput = [];
matContHetero = [];
matContAct = [];
cellHetTime = {};
cellActTime = {};
matRespDecode = [];
matRespDecodeDist = [];
cellHetPoprespCorr = {};
cellAcrossRepCorr = {};

%% create aggregate data structures
sDir = dir(strDataDir);
strOldDir = cd(strDataDir);
intCounter = 0;

for intFile=1:numel(sDir)
	%check if data part 2
	strFile = sDir(intFile).name;
	strRec = getFlankedBy(strFile,'data_aggregate','_');
	intRec = find(strcmp(strRec,cellInclude),1);
	if ~isempty(intRec)
		%load file
		load(strFile);
		fprintf('Loaded %s [%s]\n',strFile,getTime);
		
		%update population source
		intNrPops = length(cellSaveRelActValsR2);
		intAnimal = intAnimal + 1;
		vecAnimalSource(end+1) = intAnimal;
		vecSubPopulationSource(end+1) = intAnimal;
		if intNrPops > 1 %group populations
			%take mean over cohen's d's
			cellSaveActCohensD = {nanmean(cat(1,cellSaveActCohensD{1},cellSaveActCohensD{2}),1)};
			cellSaveActPPCohensD = {nanmean(cat(1,cellSaveActPPCohensD{1},cellSaveActPPCohensD{2}),1)};
			cellSaveActNPPCohensD = {nanmean(cat(1,cellSaveActNPPCohensD{1},cellSaveActNPPCohensD{2}),1)};
			cellSaveHetCohensD = {nanmean(cat(1,cellSaveHetCohensD{1},cellSaveHetCohensD{2}),1)};
			cellSaveHetPPCohensD = {nanmean(cat(1,cellSaveHetPPCohensD{1},cellSaveHetPPCohensD{2}),1)};
			cellSaveHetNPPCohensD = {nanmean(cat(1,cellSaveHetNPPCohensD{1},cellSaveHetNPPCohensD{2}),1)};
		
			%take mean over matrices
			cellSaveRelActValsR2 = {mean(cat(1,cellSaveRelActValsR2{1},cellSaveRelActValsR2{2}),1)};
			cellSaveRelActSDR2 = {mean(cat(1,cellSaveRelActSDR2{1},cellSaveRelActSDR2{2}),1)};
			
			%mean of correlations
			matSelect1 = tril(true(size(cellSaveStimCorrs{1,1})),-1);
			matSelect2 = tril(true(size(cellSaveStimCorrs{1,2})),-1);
			
			vecCorrs = nan(size(cellSaveStimCorrs,1),1);
			for intType=1:size(cellSaveStimCorrs,1)
				vecCorrs(intType) = mean([mean(cellSaveStimCorrs{intType,1}(matSelect1)) mean(cellSaveStimCorrs{intType,2}(matSelect2))]);
			end
			
			%concatenate heterogeneity measures
			for intDim1=1:size(cellSaveNormActDissim,1)
				for intDim2=1:size(cellSaveNormActDissim,2)
					for intDim4=1:size(cellSaveNormActDissim,4)
						cellSaveNormActDissim{intDim1,intDim2,1,intDim4} = [cellSaveNormActDissim{intDim1,intDim2,1,intDim4};cellSaveNormActDissim{intDim1,intDim2,2,intDim4}];
					end
				end
			end
			cellSaveNormActDissim(:,:,2,:) = [];
			
			%concatenate IT correlations
			cellSaveITC = {mean(cat(3,cellSaveITC{1},cellSaveITC{2}),3)};
		
			%concatenate RT dependency
			cellSaveRTDependency{1,1} = [cellSaveRTDependency{1,1} cellSaveRTDependency{1,2}];
			for intDim1=2:size(cellSaveRTDependency,1)
				cellSaveRTDependency{intDim1,1} = [cellSaveRTDependency{intDim1,1}; cellSaveRTDependency{intDim1,2}];
			end
			cellSaveRTDependency(:,2) = [];
			
			%take mean over heterogen
			matTempHet = cellSaveMatContHetero{1};
			matTempHet(:,:,2) = cellSaveMatContHetero{2};
			cellSaveMatContHetero = {mean(matTempHet,3)};
			
			%take mean over dF/F
			cellSaveActSigmoid{1,1} = nanmean(cat(3,cellSaveActSigmoid{1,1},cellSaveActSigmoid{2,1}),3);
			cellSaveActSigmoid{1,2} = nanmean(cat(3,cellSaveActSigmoid{1,2},cellSaveActSigmoid{2,2}),3);
			cellSaveActSigmoid{1,3} = nanmean(cat(3,cellSaveActSigmoid{1,3},cellSaveActSigmoid{2,3}),3);
			cellSaveActSigmoid(2,:) = [];
			
			%concatenate heterogen-time
			for intC=1:size(cellSaveHetTime,1)
				cellSaveHetTime{intC,1} = cat(1,cellSaveHetTime{intC,1},cellSaveHetTime{intC,2});
				cellSaveActTime{intC,1} = cat(1,cellSaveActTime{intC,1},cellSaveActTime{intC,2});
			end
			cellSaveHetTime = cellSaveHetTime(:,1);
			cellSaveActTime = cellSaveActTime(:,1);
			
			%concatenate resp type decoding
			cellSaveRespDecode{1,1} = mean(cat(4,cellSaveRespDecode{1,1},cellSaveRespDecode{2,1}),4);
			cellSaveRespDecode{1,2} = mean(cat(4,cellSaveRespDecode{1,2},cellSaveRespDecode{2,2}),4);
			cellSaveRespDecode(2,:) = [];
			
			%concatenate cellHetPoprespCorr
			cellSaveHetPoprespCorr{1,1} = [cellSaveHetPoprespCorr{1,1};cellSaveHetPoprespCorr{1,2}];
			cellSaveHetPoprespCorr{2,1} = [cellSaveHetPoprespCorr{2,1};cellSaveHetPoprespCorr{2,2}];
			cellSaveHetPoprespCorr{3,1} = mean([cellSaveHetPoprespCorr{3,1} cellSaveHetPoprespCorr{3,2}],2);
			
			cellSaveHetPoprespCorr(:,2) = [];
			
		else
			%mean of correlations
			matSelect = tril(true(size(cellSaveStimCorrs{1,1})),-1);
			
			vecCorrs = nan(size(cellSaveStimCorrs,1),1);
			for intType=1:size(cellSaveStimCorrs,1)
				vecCorrs(intType) = mean(cellSaveStimCorrs{intType,1}(matSelect));
			end
		end
		%if sum(cellBehavDetect{6,intAnimal}==0)<2
		%	cellSaveMatrices{2,1}(:,6)=nan;
		%end
		matActSigmoid = cat(3,cellSaveActSigmoid{1,1},cellSaveActSigmoid{1,2},cellSaveActSigmoid{1,3});
			
		%assign data
		matActCohensD = cat(1,matActCohensD,cellSaveActCohensD{1});
		matActPPCohensD = cat(1,matActPPCohensD,cellSaveActPPCohensD{1});
		matActNPPCohensD = cat(1,matActNPPCohensD,cellSaveActNPPCohensD{1});
		matHetCohensD = cat(1,matHetCohensD,cellSaveHetCohensD{1});
		matHetPPCohensD = cat(1,matHetPPCohensD,cellSaveHetPPCohensD{1});
		matHetNPPCohensD = cat(1,matHetNPPCohensD,cellSaveHetNPPCohensD{1});
		matRelActValsR2 = cat(1,matRelActValsR2,cellSaveRelActValsR2{1});
		matRelActSDR2 = cat(1,matRelActSDR2,cellSaveRelActSDR2{1});
		matCorrs = cat(2,matCorrs,vecCorrs);%: signal correlations [3 (within-HCAR, between, within-non-HCAR) x n (animals/blocks)] with every cell = vector of pairwise correlation values
		cellNormActDissim = cat(3,cellNormActDissim,cellSaveNormActDissim);%: within-group z-scored activation dissimilarity [2 (high/low HCAR) x 2 (detect/no-detect)] x n (animals/blocks) with every cell = vector of dissimilarity values
		matITC = cat(3,matITC,cellSaveITC{1});%: inter-trial-distance dependence of assembly consistency [2 (high/low HCAR) x 2 (hits/misses) x n (animals/blocks)] with every cell = [p x 2] matrix with intertrial distance and correlation values for each trial pair
		cellRTDependency = cat(2,cellRTDependency,cellSaveRTDependency);%: RT dependence of activity dissimilarity [3 (RT/HCAR/non-HCAR) x n (animals/blocks)] with every cell = vector of trials
		matContHetero = cat(3,matContHetero,cellSaveMatContHetero{1});
		matContAct = cat(4,matContAct,matActSigmoid);
		cellHetTime = cat(2,cellHetTime,cellSaveHetTime);
		cellActTime = cat(2,cellActTime,cellSaveActTime);
		matRespDecode = cat(4,matRespDecode,cellSaveRespDecode{1,1});
		matRespDecodeDist = cat(3,matRespDecodeDist,cellSaveRespDecode{1,2});
		cellHetPoprespCorr = cat(2,cellHetPoprespCorr,cellSaveHetPoprespCorr);
		intCounter = intCounter + 1;
	end
end

%% meta analyses part 1
cd(strFigDir);
%close all;
intAnimals = intCounter;
vecContrasts = [0 0.5 2 8 32 100];

%% meta analyses part 2


%% cohen's d
vecActD = nanmean(matActCohensD(:,2:5),2);
vecActPPD = nanmean(matActPPCohensD(:,2:5),2);
vecActNPPD = nanmean(matActNPPCohensD(:,2:5),2);
vecHetD = nanmean(matHetCohensD(:,2:5),2);
vecHetPPD = nanmean(matHetPPCohensD(:,2:5),2);
vecHetNPPD = nanmean(matHetNPPCohensD(:,2:5),2);
[h,p1]=ttest(vecActD,vecHetD);
[h,p2]=ttest(vecActPPD,vecHetPPD);
[h,p3]=ttest(vecActNPPD,vecHetNPPD);
[h crit_p vecCohenP_adj]=fdr_bh([p1 p2 p3],0.05,'pdep');

figure
hold on
errorbar(1/6,mean(vecActD),std(vecActD)/sqrt(length(vecActD)),'bx')
errorbar(2/6,mean(vecHetD),std(vecHetD)/sqrt(length(vecHetD)),'kx')
errorbar(3/6,mean(vecActPPD),std(vecActPPD)/sqrt(length(vecActPPD)),'bx')
errorbar(4/6,mean(vecHetPPD),std(vecHetPPD)/sqrt(length(vecHetPPD)),'kx')
errorbar(5/6,mean(vecActNPPD),std(vecActNPPD)/sqrt(length(vecActNPPD)),'bx')
errorbar(6/6,mean(vecHetNPPD),std(vecHetNPPD)/sqrt(length(vecHetNPPD)),'kx')
hold off
xlim([-0.1 1.1])
ylabel('Cohen''s D hit-miss difference')
set(gca,'xtick',(1/6):(1/6):1,'xticklabel',{'dF/F0 All','Het All','Act Pref-pop','Het Pref-pop','Act NPP','Het NPP'})
title(sprintf('Behavioral correlates test contrasts, paired t-tests het vs dF/F0, all=%.3f, pref=%.3f, np=%.3f',p1,p2,p3));

%{
vecHetD-vecActD
vecHetPPD-vecActPPD
vecHetNPPD-vecActNPPD

vecActD
vecHetD
vecHetPPD
vecHetNPPD
%}

if boolSavePlots
	intFigCounter = intFigCounter + 1;
	figure(gcf);
	drawnow;
	strFig = sprintf('Meta%d_cohensD_raw',intFigCounter);
	export_fig([strFig '.tif']);
	export_fig([strFig '.pdf']);
end


%% subgroup analysis
%cellSaveRelActValsR2 [dblPercExplainedNeurons dblPercExplainedNeuronsAdj dblPercExplainedTrials dblPercExplainedTrialsAdj dblPercExplainedAll intNrSig];
%cellSaveRelActSDR2 [dblNrSDNeurons dblNrSDTrials dblNrSDAll dblNrSDSig];
		
vecR2Neurons = matRelActValsR2(:,1);
vecR2NeuronsSD = matRelActSDR2(:,1);
vecR2Trials = matRelActValsR2(:,3);
vecR2TrialsSD = matRelActSDR2(:,2);
vecR2All = matRelActValsR2(:,5);
vecR2AllSD = matRelActSDR2(:,3);
vecSigNeurons = matRelActValsR2(:,6);
vecSigNeuronsSD = matRelActSDR2(:,4);

%CI nr of significant neurons
figure
subplot(1,4,1)
%errorbar(0.6,mean(vecR2Neurons),std(vecR2Neurons)*2,'bx');
%hold on
scatter(0.5*ones(size(vecR2Neurons)),vecR2Neurons,'rx')
%hold off
set(gca,'xtick',0.5,'xticklabel','R^2 neurons')
ylabel('Explained variance by neuron ID')
title(sprintf('Significant; %d/%d',sum(~(vecR2NeuronsSD>-2 & vecR2NeuronsSD<2)),length(vecR2NeuronsSD)))
set(gca,'ylim',[0 0.4]);

subplot(1,4,2)
%errorbar(0.6,mean(vecR2Trials),std(vecR2Trials)*2,'bx');
%hold on
scatter(0.5*ones(size(vecR2Trials)),vecR2Trials,'rx')
%hold off
set(gca,'xtick',0.5,'xticklabel','Nr sig. neurons')
ylabel('Explained variance by trial ID')
title(sprintf('Significant; %d/%d',sum(~(vecR2TrialsSD>-2 & vecR2TrialsSD<2)),length(vecR2TrialsSD)))
set(gca,'ylim',[0 0.4]);

subplot(1,4,3)
%errorbar(0.6,mean(vecR2All),std(vecR2All)*2,'bx');
%hold on
scatter(0.5*ones(size(vecR2All)),vecR2All,'rx')
%hold off
set(gca,'xtick',0.5,'xticklabel','R^2 trials+neuron')
ylabel('Explained variance by trial+neuron ID')
title(sprintf('Significant; %d/%d',sum(~(vecR2AllSD>-2 & vecR2AllSD<2)),length(vecR2AllSD)))
set(gca,'ylim',[0 0.4]);

subplot(1,4,4)
%errorbar(0.6,mean(vecR2All),std(vecR2All)*2,'bx');
%hold on
scatter(0.5*ones(size(vecSigNeurons)),vecSigNeurons,'rx')
%hold off
set(gca,'xtick',0.5,'xticklabel','Sign. Neurons')
ylabel('Fraction of sign. hit-modulated neurons')
title(sprintf('Significant; %d/%d',sum(~(vecSigNeuronsSD>-2 & vecSigNeuronsSD<2)),length(vecSigNeuronsSD)))
set(gca,'ylim',[0 max(get(gca,'ylim'))]);

if boolSavePlots
	intFigCounter = intFigCounter + 1;
	figure(gcf);
	drawnow;
	strFig = sprintf('Meta%d_detectcorrelated_predictability_raw',intFigCounter);
	export_fig([strFig '.tif']);
	export_fig([strFig '.pdf']);
end
%return
%% signal correlations
hSignalCorrs = figure;
intPopulations = size(matCorrs,2);

vecSC_Miss = matCorrs(1,:);
vecSC_Hit = matCorrs(3,:);
vecSC_Fast = matCorrs(5,:);
vecSC_Slow = matCorrs(7,:);

errorbar(1:4,[mean(vecSC_Miss) mean(vecSC_Hit) mean(vecSC_Fast) mean(vecSC_Slow)],[std(vecSC_Miss) std(vecSC_Hit) std(vecSC_Fast) std(vecSC_Slow)]/sqrt(intPopulations),'Linestyle','none','Marker','x','Color','k');
set(gca,'XTick',1:4,'XTickLabel',{'Miss','Hit-All','Hit-Fast','Hit-Slow'})


[h,dblP_MH] = ttest2(vecSC_Miss,vecSC_Hit); %assembly/non-assembly
[h,dblP_MF] = ttest2(vecSC_Miss,vecSC_Fast);%assembly/between
[h,dblP_MS] = ttest2(vecSC_Miss,vecSC_Slow);%between/non-assembly
[h,dblP_FS] = ttest2(vecSC_Fast,vecSC_Slow);%between/non-assembly

title(sprintf('Signal correlations;p-diff; MH=%.3f; MF=%.3f; MS=%.3f; FS=%.3f',dblP_MH,dblP_MF,dblP_MS,dblP_FS))


if boolSavePlots
	intFigCounter = intFigCounter + 1;
	drawnow;
	strFig = sprintf('Meta%d_signalcorrelations_raw',intFigCounter);
	export_fig([strFig '.tif']);
	export_fig([strFig '.pdf']);
end

%% noise correlations
hNoiseCorrs = figure;
intPopulations = size(matCorrs,2);

vecNC_Miss = matCorrs(2,:);
vecNC_Hit = matCorrs(4,:);
vecNC_Fast = matCorrs(6,:);
vecNC_Slow = matCorrs(8,:);

errorbar(1:4,[mean(vecNC_Miss) mean(vecNC_Hit) mean(vecNC_Fast) mean(vecNC_Slow)],[std(vecNC_Miss) std(vecNC_Hit) std(vecNC_Fast) std(vecNC_Slow)]/sqrt(intPopulations),'Linestyle','none','Marker','x','Color','k');
set(gca,'XTick',1:4,'XTickLabel',{'Miss','Hit-All','Hit-Fast','Hit-Slow'})
ylim([0.2 0.4])

[h,dblP_MH] = ttest2(vecNC_Miss,vecNC_Hit); %assembly/non-assembly
[h,dblP_MF] = ttest2(vecNC_Miss,vecNC_Fast);%assembly/between
[h,dblP_MS] = ttest2(vecNC_Miss,vecNC_Slow);%between/non-assembly
[h,dblP_FS] = ttest2(vecNC_Fast,vecNC_Slow);%between/non-assembly

title(sprintf('Noise correlations;p-diff; MH=%.3f; MF=%.3f; MS=%.3f; FS=%.3f',dblP_MH,dblP_MF,dblP_MS,dblP_FS))


if boolSavePlots
	intFigCounter = intFigCounter + 1;
	drawnow;
	strFig = sprintf('Meta%d_noisecorrelations_raw',intFigCounter);
	export_fig([strFig '.tif']);
	export_fig([strFig '.pdf']);
end

%% relative change in correlations
%diff fast/slow OR relative fast and slow vs miss

hDiffCorrs = figure;
scatter(vecNC_Fast-vecNC_Slow,vecSC_Fast-vecSC_Slow,'bx')
xlim([-0.6 0.6])
ylim([-0.6 0.6])
hold on
plot(get(gca,'XLim'),[0 0],'k--')
plot([0 0],get(gca,'YLim'),'k--')

[hSC,pSC,ciSC] = ttest(vecSC_Fast-vecSC_Slow);
[hNC,pNC,ciNC] = ttest(vecNC_Fast-vecNC_Slow);

plot([mean(vecNC_Fast-vecNC_Slow) mean(vecNC_Fast-vecNC_Slow)],ciSC,'r')
plot(ciNC,[mean(vecSC_Fast-vecSC_Slow) mean(vecSC_Fast-vecSC_Slow)],'r')
hold off

xlabel('d(Noise correlation) Fast-Slow')
ylabel('d(Signal correlation) Fast-Slow')

title(sprintf('Difference in corrs; SC, p=%.3f; NC, p=%.3f; cross is 95%% CI',pSC,pNC))

if boolSavePlots
	intFigCounter = intFigCounter + 1;
	drawnow;
	strFig = sprintf('Meta%d_diffcorrelations_raw',intFigCounter);
	export_fig([strFig '.tif']);
	export_fig([strFig '.pdf']);
end


%% within-group z-scored activation dissimilarity
%% reaction-time dependency
vecX = [0 3];
hRTDependenceFig = figure;

vecPs = nan(1,8);
for intActType=1:4
	if intActType == 1
		%dissimilarity
		intHCAR = 2;
		intNonHCAR = 3;
		strLabelY = 'Mean population activity heterogeneity';
	elseif intActType == 2
		%z-scored act
		strLabelY = 'Mean population z-scored activity';
	elseif intActType == 3
		%dF/F act
		strLabelY = 'Mean population dF/F activity';
	elseif intActType == 4
		%dF/F act
		strLabelY = 'Variance of population activity';
	end
	
	matReg = zeros(intPopulations,length(vecX));
	vecSlopes = zeros(intPopulations,1);
	
	subplot(2,4,intActType);
	hold on;
	vecAggRTs = [];
	vecAggAct = [];
	for intPopulation=1:intPopulations
		vecAggRTs = [vecAggRTs cellRTDependency{1,intPopulation}];
		vecAggAct = [vecAggAct cellRTDependency{intActType+1,intPopulation}'];
		
		%do regressions per population
		sStatsC=regstats(cellRTDependency{intActType+1,intPopulation}',cellRTDependency{1,intPopulation},'linear');
		vecY = polyval(sStatsC.beta([2 1]),vecX);
		matReg(intPopulation,:) = vecY;
		vecSlopes(intPopulation) = sStatsC.beta(2);
		plot(vecX,vecY,'Color',[0.5 0.5 1.0]);
	end
	
	%plot means
	vecMeanY = mean(matReg,1);
	plot(vecX,vecMeanY,'Color',[1 0 0],'LineWidth',2)
	
	%ttest
	[h,dblPSlope] = ttest(vecSlopes);
	[corrected_p, h]=bonf_holm([ones([1 7]) dblPSlope],0.05);
	dblPSlopeCorr = corrected_p(end);
	vecPs((intActType-1)*2+1) = dblPSlope;
	title(sprintf('Mean of linear regressions over animals; slope-p=%.3f',dblPSlopeCorr))
	xlabel('Reaction Time (s)')
	ylabel(strLabelY)
	
	%plot aggregate
	subplot(2,4,intActType+4);
	scatter(vecAggRTs,vecAggAct,'bx')
	hold on
	drawnow
	
	%perform regressions
	sStatsC=regstats(vecAggAct,vecAggRTs,'linear');
	vecX = [0 3];
	vecY = polyval(sStatsC.beta([2 1]),vecX);
	plot(vecX,vecY,'r')
	
	hold off
	
	title(sprintf('Aggregate data set; slope=%.3f; p=%.3f; R^2=%.3f',...
		sStatsC.beta(2),sStatsC.tstat.pval(2),sStatsC.rsquare))
	vecPs((intActType-1)*2+2) = sStatsC.tstat.pval(2);
	xlabel('Reaction Time (s)')
	ylabel(strLabelY)
end

[corrected_p, h]=bonf_holm(vecPs,0.05);
fprintf('\nCorrected p-values RT dependence: %.3f\n',corrected_p)

if boolSavePlots
	intFigCounter = intFigCounter + 1;
	drawnow;
	jFig = get(handle(gcf), 'JavaFrame');
	jFig.setMaximized(true);
	figure(gcf);
	drawnow;
	strFig = sprintf('Meta%d_activationdissimilarity_RTDependency_raw',intFigCounter);
	export_fig([strFig '.tif']);
	export_fig([strFig '.pdf']);
end

%% inter-trial-distance dependence of assembly consistency
%for whole pop or only preferred pop
matMean = median(matITC,3);
matErr = std(matITC,[],3)/sqrt(intAnimals);
hPopCorrPlot = figure;

matErrLow = (matMean-quantile(matITC,normcdf(-1,0,1),3))/sqrt(intAnimals);
matErrHigh = (quantile(matITC,normcdf(1,0,1),3)-matMean)/sqrt(intAnimals);

subplot(2,2,1)
cellColor={'r','m','g'};
hold on;
for intResp=1:3
	errorbar(intResp,matMean(1,intResp),matErr(1,intResp),['x' cellColor{intResp}]);
end
hold off
title('Non-preferred population');
cellLabels = {'Miss','Slow','Fast'};
set(gca,'xtick',1:3,'xticklabel',cellLabels)
ylabel('Inter-trial correlation');
ylim([0 0.15]);

subplot(2,2,2)
hold on;
for intResp=1:3
	errorbar(intResp,matMean(2,intResp),matErr(2,intResp),['x' cellColor{intResp}]);
end
hold off
title('Preferred population');
set(gca,'xtick',1:3,'xticklabel',cellLabels)
ylabel('Inter-trial correlation');
ylim([0 0.15]);

[h,dblP_NP_MS] = ttest(matITC(1,1,:),matITC(1,2,:));
[h,dblP_NP_MF] =  ttest(matITC(1,1,:),matITC(1,3,:));
[h,dblP_NP_SF] =  ttest(matITC(1,2,:),matITC(1,3,:));

[h,dblP_P_MS] = ttest(matITC(2,1,:),matITC(2,2,:));
[h,dblP_P_MF] =  ttest(matITC(2,1,:),matITC(2,3,:));
[h,dblP_P_SF] =  ttest(matITC(2,2,:),matITC(2,3,:));

[h,pNPM]=ttest(matITC(1,1,:));
[h,pNPS]=ttest(matITC(1,2,:));
[h,pNPF]=ttest(matITC(1,3,:));
[h,pPM]=ttest(matITC(2,1,:));
[h,pPS]=ttest(matITC(2,2,:));
[h,pPF]=ttest(matITC(2,3,:));

vecP_adjP = bonf_holm([pPM pPS pPF]);
vecP_adjNP = bonf_holm([pNPM pNPS pNPF]);

if boolSavePlots
	intFigCounter = intFigCounter + 1;
	figure(hPopCorrPlot);
	drawnow;
	jFig = get(handle(hPopCorrPlot), 'JavaFrame');
	jFig.setMaximized(true);
	figure(hPopCorrPlot);
	drawnow;
	strFig = sprintf('Meta%d_intertrialcorr_%s_pop%d_raw',intFigCounter);
	export_fig([strFig '.tif']);
	export_fig([strFig '.pdf']);
end

%% heterogeneity over contrasts
%normalize heterogeneity per animal
matContHeteroNorm = matContHetero;%zeros(size(matContHetero)); %1=all-hit,2=all-miss,3=pref-hit,4=pref-miss,5=np-hit,6=np-miss
for intAnimal=1:size(matContHetero,3)
	matTempHet = matContHetero(:,:,intAnimal);
	%matContHeteroNorm(:,:,intAnimal) = (matTempHet ./ mean(matTempHet(:)))-0.1;
end

%pre-compute variables
vecContrasts = [0.2 0.5 2 8 32 100];
vecWindow = [1 length(vecContrasts)];
vecWindowSelect = vecWindow(1):vecWindow(end);
intWL = length(vecWindowSelect);
vecWindowInv = intWL:-1:1;
vecWindowPlotInv = vecWindow(end):-1:vecWindow(1);
vecX = [vecContrasts(vecWindowSelect) vecContrasts(vecWindowPlotInv)];
vecLineX = vecContrasts(vecWindowSelect);

%get data
hHetCon=figure;
cellStrPopType = {'Pop-all','Pref-pop','Non-pref-pop'};
for intPopType=1:3
	subplot(2,2,intPopType)
	
	%get data
	matHitHet = matContHeteroNorm(:,1+(intPopType-1)*2,:);
	matMissHet = matContHeteroNorm(:,2+(intPopType-1)*2,:);
	vecHitY = nanmean(matHitHet,3)';
	vecHitE = nanstd(matHitHet,[],3)'/sqrt(size(matContHeteroNorm,3));
	vecMissY = nanmean(matMissHet,3)';
	vecMissE = nanstd(matMissHet,[],3)'/sqrt(size(matContHeteroNorm,3));
	
	%get sig
	vecP = zeros(1,size(matContHeteroNorm,1));
	for intC=1:size(matContHeteroNorm,1)
		[h,p,ci] = ttest(matHitHet(intC,1,:),matMissHet(intC,1,:));
		
		%put in vector
		vecP(intC) = p;
	end
	
	%overall t-test
	matHetHits = matHitHet(2:5,1,:);
	matHetMisses = matMissHet(2:5,1,:);
	[h,dblP_All,ci] = ttest(matHetHits(:),matHetMisses(:));
	
	for intResp=[0 1]
		if intResp == 1
			vecMeanTrace = vecHitY(vecWindowSelect);
			vecSE = vecHitE(vecWindowSelect);
			vecColorFill = [0.7 1 0.7];
			vecColorLine = [0 1 0];
		else
			vecColorLine = [1 0 0];
			vecColorFill = [1 0.7 0.7];
			vecMeanTrace = vecMissY(vecWindowSelect);
			vecSE = vecMissE(vecWindowSelect);
		end
		
		vecMinTrace = vecMeanTrace-vecSE;
		vecMaxTrace = vecMeanTrace+vecSE;
		vecY = [vecMinTrace vecMaxTrace(vecWindowInv)];
		
		%plot
		hold on
		fill(vecX,vecY,vecColorFill,'EdgeColor','none');
		plot(vecLineX,vecMeanTrace,'-','LineWidth',2,'Color',vecColorLine);
		hold off
	end
	set(gca,'XScale','log','YScale','linear')
	set(gca,'XTick',vecLineX,'XTickLabel',[0 vecLineX(2:end)])
	
	title(sprintf('%s het during stim; p-vals: 0=%.3f; 0.5=%.3f; 2=%.3f; 8=%.3f; 32=%.3f; 100=%.3f; Overall,p=%.3f',cellStrPopType{intPopType},[vecP dblP_All]))
	grid on
	xlabel('Stimulus contrast')
	ylabel(['Heterogeneity ' cellStrPopType{intPopType}])
	xlim(vecContrasts(vecWindow))
	if intPopType==2
		ylim([0.6-eps 1.6+eps])
	else
		ylim([0.6-eps 1.1+eps])
	end
	legend({'SEM','Miss','SEM','Hit'},'Location','Best')
	drawnow;
end
strFigTitle = 'hetero_over_contrasts';

if boolSavePlots
	intFigCounter = intFigCounter + 1;
	drawnow;
	jFig = get(handle(gcf), 'JavaFrame');
	jFig.setMaximized(true);
	figure(gcf);
	drawnow;
	strFig = sprintf('Meta%d_%s_raw',intFigCounter,strFigTitle);
	export_fig([strFig '.tif']);
	export_fig([strFig '.pdf']);
end

%% dF/F0 over contrasts
%normalize act per animal
matContActNorm = matContAct;
for intAnimal=1:size(matContActNorm,4)
	matTempAct = matContActNorm(:,:,:,intAnimal);
	matContActNorm(:,:,:,intAnimal) = matTempAct - min(matTempAct(:));
end

%pre-compute variables
vecContrasts = [0.2 0.5 2 8 32 100];
vecWindow = [1 length(vecContrasts)];
vecWindowSelect = vecWindow(1):vecWindow(end);
intWL = length(vecWindowSelect);
vecWindowInv = intWL:-1:1;
vecWindowPlotInv = vecWindow(end):-1:vecWindow(1);
vecX = [vecContrasts(vecWindowSelect) vecContrasts(vecWindowPlotInv)];
vecLineX = vecContrasts(vecWindowSelect);

%get data
hActCon=figure;
cellStrPopType = {'Pref-pop','Pop-all','Non-pref-pop'};
for intPopType=1:3 %pref;all;non-pref
	subplot(2,2,intPopType)
	
	%get data
	matHitAct = matContActNorm(2,:,intPopType,:);
	matMissAct = matContActNorm(1,:,intPopType,:);
	vecHitY = nanmean(matHitAct,4);
	vecHitE = nanstd(matHitAct,[],4)/sqrt(size(matContActNorm,4));
	vecMissY = nanmean(matMissAct,4);
	vecMissE = nanstd(matMissAct,[],4)/sqrt(size(matContActNorm,4));
	
	%get sig
	vecP = zeros(1,size(matContActNorm,2));
	for intC=1:size(matContActNorm,2)
		[h,p,ci] = ttest(matHitAct(1,intC,1,:),matMissAct(1,intC,1,:));
		
		%put in vector
		vecP(intC) = p;
	end
	
	%overall t-test
	matActHits = matHitAct(1,2:5,1,:);
	matActMisses = matMissAct(1,2:5,1,:);
	[h,dblP_All,ci] = ttest(matActHits(:),matActMisses(:));
	
	for intResp=[0 1]
		if intResp == 1
			vecMeanTrace = vecHitY(vecWindowSelect);
			vecSE = vecHitE(vecWindowSelect);
			vecColorFill = [0.7 1 0.7];
			vecColorLine = [0 1 0];
		else
			vecColorLine = [1 0 0];
			vecColorFill = [1 0.7 0.7];
			vecMeanTrace = vecMissY(vecWindowSelect);
			vecSE = vecMissE(vecWindowSelect);
		end
		
		vecMinTrace = vecMeanTrace-vecSE;
		vecMaxTrace = vecMeanTrace+vecSE;
		vecY = [vecMinTrace vecMaxTrace(vecWindowInv)];
		
		%plot
		hold on
		fill(vecX,vecY,vecColorFill,'EdgeColor','none');
		plot(vecLineX,vecMeanTrace,'-','LineWidth',2,'Color',vecColorLine);
		hold off
	end
	set(gca,'XScale','log','YScale','linear')
	set(gca,'XTick',vecLineX,'XTickLabel',[0 vecLineX(2:end)])
	
	title(sprintf('%s Act during stim; p-vals: 0=%.3f; 0.5=%.3f; 2=%.3f; 8=%.3f; 32=%.3f; 100=%.3f; Overall,p=%.3f',cellStrPopType{intPopType},[vecP dblP_All]))
	grid on
	xlabel('Stimulus contrast')
	ylabel(['Act (dfof) ' cellStrPopType{intPopType}])
	xlim(vecContrasts(vecWindow))
	ylim([0-eps 0.06+eps])
	legend({'SEM','Miss','SEM','Hit'},'Location','Best')
	drawnow;
end
strFigTitle = 'act_over_contrasts';

if boolSavePlots
	intFigCounter = intFigCounter + 1;
	drawnow;
	jFig = get(handle(gcf), 'JavaFrame');
	jFig.setMaximized(true);
	figure(gcf);
	drawnow;
	strFig = sprintf('Meta%d_%s_raw',intFigCounter,strFigTitle);
	export_fig([strFig '.tif']);
	export_fig([strFig '.pdf']);
end

%% heterogeneity+dF/F over time
%cellHetTime{intC,intAnimal} = {trial x time}[1=miss] {trial x time}[2=slow] {trial x time}[3=fast]
%normalize heterogeneity per animal

%pre-allocate for later usage
matDiffFS = nan(2,size(cellHetTime,1),size(cellHetTime,2),size(cellHetTime{1,1}{1},2));
for intAct = [0 1];
	if intAct == 0
		cellDataTime = cellHetTime;
	else
		cellDataTime = cellActTime;
	end
	
	cellRawDataMiss = {};
	cellRawDataSlow = {};
	cellRawDataFast = {};
	
	
	for intAnimal=1:size(cellDataTime,2)
		if intAct == 0
			matTempHet = [];
			for intC=1:size(cellDataTime,1)
				matTempHet = [matTempHet;cellDataTime{intC,intAnimal}{1};cellDataTime{intC,intAnimal}{2};cellDataTime{intC,intAnimal}{3}];
			end
			dblMean = mean(matTempHet(:)); %mean over all values of animal
		end
		%apply normalization
		for intC=1:size(cellDataTime,1)
			if intAct == 0
				cellDataTime{intC,intAnimal}{1} = cellDataTime{intC,intAnimal}{1} ./ dblMean;
				cellDataTime{intC,intAnimal}{2} = cellDataTime{intC,intAnimal}{2} ./ dblMean;
				cellDataTime{intC,intAnimal}{3} = cellDataTime{intC,intAnimal}{3} ./ dblMean;
			end
			cellRawDataMiss{intC}(intAnimal,:) = mean(cellDataTime{intC,intAnimal}{1},1);
			cellRawDataSlow{intC}(intAnimal,:) = mean(cellDataTime{intC,intAnimal}{2},1);
			cellRawDataFast{intC}(intAnimal,:) = mean(cellDataTime{intC,intAnimal}{3},1);
			
			%put in output
			matDiffFS(intAct+1,intC,intAnimal,:) =  cellRawDataFast{intC}(intAnimal,:) - cellRawDataSlow{intC}(intAnimal,:);
		end
	end
	
	%general
	vecC = [0 0.5 2 8 32 100];
	vecWindowSecs = [-3 5];
	vecWindow = round(vecWindowSecs*25.4);
	vecLineX = (vecWindow(1):vecWindow(end))/25.4;
	vecWindowInv = length(vecLineX):-1:1;
	vecX = [vecLineX vecLineX(vecWindowInv)];
	
	%plot
	hHetTime = figure;
	for intC=1:length(vecC)
		subplot(2,3,intC);
		for intType=1:3
			if intType == 1 %miss
				matRawData = cellRawDataMiss{intC};
				strType = 'Miss';
				vecColorLine = [1 0 0];
				vecColorFill = [1 0.7 0.7];
			elseif intType == 2 %slow
				matRawData = cellRawDataSlow{intC};
				strType = 'Slow';
				vecColorLine = [1 1 0];
				vecColorFill = [1 1 0.7];
			else %fast
				matRawData = cellRawDataFast{intC};
				strType = 'Fast';
				vecColorLine = [0 1 0];
				vecColorFill = [0.7 1 0.7];
			end
			
			
			%get data
			vecMeanTrace = mean(matRawData,1);
			vecSE = std(matRawData,[],1)./sqrt(size(cellDataTime,2));
			vecMinTrace = vecMeanTrace-vecSE;
			vecMaxTrace = vecMeanTrace+vecSE;
			vecY = [vecMinTrace vecMaxTrace(vecWindowInv)];
			
			%plot
			hold on
			fill(vecX,vecY,vecColorFill,'EdgeColor','none');
			plot(vecLineX,vecMeanTrace,'-','LineWidth',2,'Color',vecColorLine);
			hold off
		end
		
		%labels
		if intAct == 0
			title(sprintf('Pop heterogeneity over time; Contrast %.1f',vecC(intC)))
			ylabel('Normalized population response heterogeneity')
			ylim([0.8 1.6])
			strFigTitle = 'hetero_over_time';
		else
			title(sprintf('Pop dF/F0 over time; Contrast %.1f',vecC(intC)))
			ylabel('Normalized population dF/F0')
			ylim([-0.02 0.06])
			strFigTitle = 'activity_over_time';
		end
		grid on
		xlabel('Time after stim onset (s)')
		xlim(vecWindowSecs)
		drawnow;
	end
	
	if boolSavePlots
		intFigCounter = intFigCounter + 1;
		drawnow;
		jFig = get(handle(gcf), 'JavaFrame');
		jFig.setMaximized(true);
		figure(gcf);
		drawnow;
		strFig = sprintf('Meta%d_%s_raw',intFigCounter,strFigTitle);
		export_fig([strFig '.tif']);
		export_fig([strFig '.pdf']);
	end
end

%% plot resp type predictability
vecHetMiss = [];
vecHetSlow = [];
vecHetFast = [];
vecActMiss = [];
vecActSlow = [];
vecActFast = [];

%which data points?
vecWindowSecs = [-3 5];
vecWindow = round(vecWindowSecs*25.4);
vecLineX = (vecWindow(1):vecWindow(end))/25.4;
vecSelectSecs = [-3 0];
vecSelect = round(vecSelectSecs*25.4);
intStart = find(vecLineX>=vecSelectSecs(1),1);
intStop = find(vecLineX>=vecSelectSecs(2),1);

%get data
for intAnimal=1:size(cellHetTime,2)
	matTemp = [];
	for intC=2:size(cellHetTime,1)
		matTemp = [matTemp;cellHetTime{intC,intAnimal}{1};cellHetTime{intC,intAnimal}{2};cellHetTime{intC,intAnimal}{3}];
	end
	dblMean = mean(matTemp(:)); %mean over all values of animal
	
	%apply normalization
	for intC=2:size(cellHetTime,1)
		cellHetTime{intC,intAnimal}{1} = cellHetTime{intC,intAnimal}{1} ./ dblMean;
		cellHetTime{intC,intAnimal}{2} = cellHetTime{intC,intAnimal}{2} ./ dblMean;
		cellHetTime{intC,intAnimal}{3} = cellHetTime{intC,intAnimal}{3} ./ dblMean;
	end
	
	%prep data
	matTempHetM = [];
	matTempHetS = [];
	matTempHetF = [];
	matTempActM = [];
	matTempActS = [];
	matTempActF = [];
	for intC=2:size(cellHetTime,1)
		matTempHetM = [matTempHetM;cellHetTime{intC,intAnimal}{1}];
		matTempHetS = [matTempHetS;cellHetTime{intC,intAnimal}{2}];
		matTempHetF = [matTempHetF;cellHetTime{intC,intAnimal}{3}];
		matTempActM = [matTempActM;cellActTime{intC,intAnimal}{1}];
		matTempActS = [matTempActS;cellActTime{intC,intAnimal}{2}];
		matTempActF = [matTempActF;cellActTime{intC,intAnimal}{3}];
	end
	vecHetMiss(intAnimal) = mean(mean(matTempHetM(:,intStart:intStop),1));
	vecHetSlow(intAnimal) = mean(mean(matTempHetS(:,intStart:intStop),1));
	vecHetFast(intAnimal) = mean(mean(matTempHetF(:,intStart:intStop),1));
	vecActMiss(intAnimal) = mean(mean(matTempActM(:,intStart:intStop),1));
	vecActSlow(intAnimal) = mean(mean(matTempActS(:,intStart:intStop),1));
	vecActFast(intAnimal) = mean(mean(matTempActF(:,intStart:intStop),1));
end
% check separability
%{
calculate per point intra-cluster vs inter-cluster difference (ratio) for
dF/F and for heterogeneity dimension; will give distribution of
separability for each cluster; then perform ANOVA/ttests for
between-cluster distribution comparison
%}

%get cluster means
dblMeanActMiss = mean(vecActMiss);
dblMeanActSlow = mean(vecActSlow);
dblMeanActFast = mean(vecActFast);
dblMeanHetMiss = mean(vecHetMiss);
dblMeanHetSlow = mean(vecHetSlow);
dblMeanHetFast = mean(vecHetFast);
intAnimals = length(vecActMiss);

%separability; inter vs intra cluster distance
vecSepRatioActMiss_w_Slow = abs(vecActMiss - dblMeanActSlow) ./ (abs(vecActMiss - dblMeanActMiss) + abs(vecActMiss - dblMeanActSlow));%inter / (inter+intra) distance
vecSepRatioActMiss_w_Fast = abs(vecActMiss - dblMeanActFast) ./ (abs(vecActMiss - dblMeanActMiss) + abs(vecActMiss - dblMeanActFast));
vecSepRatioActSlow_w_Miss = abs(vecActSlow - dblMeanActMiss) ./ (abs(vecActSlow - dblMeanActSlow) + abs(vecActSlow - dblMeanActMiss));
vecSepRatioActSlow_w_Fast = abs(vecActSlow - dblMeanActFast) ./ (abs(vecActSlow - dblMeanActSlow) + abs(vecActSlow - dblMeanActFast));
vecSepRatioActFast_w_Slow = abs(vecActFast - dblMeanActSlow) ./ (abs(vecActFast - dblMeanActFast) + abs(vecActFast - dblMeanActSlow));
vecSepRatioActFast_w_Miss = abs(vecActFast - dblMeanActMiss) ./ (abs(vecActFast - dblMeanActFast) + abs(vecActFast - dblMeanActMiss));
vecSepRatioHetMiss_w_Slow = abs(vecHetMiss - dblMeanHetSlow) ./ (abs(vecHetMiss - dblMeanHetMiss) + abs(vecHetMiss - dblMeanHetSlow));
vecSepRatioHetMiss_w_Fast = abs(vecHetMiss - dblMeanHetFast) ./ (abs(vecHetMiss - dblMeanHetMiss) + abs(vecHetMiss - dblMeanHetFast));
vecSepRatioHetSlow_w_Miss = abs(vecHetSlow - dblMeanHetMiss) ./ (abs(vecHetSlow - dblMeanHetSlow) + abs(vecHetSlow - dblMeanHetMiss));
vecSepRatioHetSlow_w_Fast = abs(vecHetSlow - dblMeanHetFast) ./ (abs(vecHetSlow - dblMeanHetSlow) + abs(vecHetSlow - dblMeanHetFast));
vecSepRatioHetFast_w_Slow = abs(vecHetFast - dblMeanHetSlow) ./ (abs(vecHetFast - dblMeanHetFast) + abs(vecHetFast - dblMeanHetSlow));
vecSepRatioHetFast_w_Miss = abs(vecHetFast - dblMeanHetMiss) ./ (abs(vecHetFast - dblMeanHetFast) + abs(vecHetFast - dblMeanHetMiss));

%act agg
vecActSlowMiss = [vecSepRatioActMiss_w_Slow vecSepRatioActSlow_w_Miss];
vecActSlowFast = [vecSepRatioActFast_w_Slow vecSepRatioActSlow_w_Fast];
vecActFastMiss = [vecSepRatioActMiss_w_Fast vecSepRatioActFast_w_Miss];

%het agg
vecHetSlowMiss = [vecSepRatioHetMiss_w_Slow vecSepRatioHetSlow_w_Miss];
vecHetSlowFast = [vecSepRatioHetFast_w_Slow vecSepRatioHetSlow_w_Fast];
vecHetFastMiss = [vecSepRatioHetMiss_w_Fast vecSepRatioHetFast_w_Miss];

%fig
hFigRTP = figure;

%plot
subplot(2,2,1)
scatter(vecHetMiss,vecActMiss,'rx')
hold on
scatter(vecHetSlow,vecActSlow,'yx')
scatter(vecHetFast,vecActFast,'gx')
hold off
title('Stim; Response type predictability')
xlabel('Heterogeneity during 1s preceding stimulus')
ylabel('dF/F during 1s preceding stimulus')

subplot(2,2,2)

dblOffset = 0.1;
plot([0.5 3.5],[0.5 0.5],'k--')
%plot act
hold on
errorbar(1-dblOffset,mean(vecActSlowMiss),std(vecActSlowMiss)/sqrt(length(vecActSlowMiss)),'bx')
errorbar(2-dblOffset,mean(vecActSlowFast),std(vecActSlowFast)/sqrt(length(vecActSlowFast)),'bx')
errorbar(3-dblOffset,mean(vecActFastMiss),std(vecActFastMiss)/sqrt(length(vecActFastMiss)),'bx')
%plot het
errorbar(1+dblOffset,mean(vecHetSlowMiss),std(vecHetSlowMiss)/sqrt(length(vecHetSlowMiss)),'kx')
errorbar(2+dblOffset,mean(vecHetSlowFast),std(vecHetSlowFast)/sqrt(length(vecHetSlowFast)),'kx')
errorbar(3+dblOffset,mean(vecHetFastMiss),std(vecHetFastMiss)/sqrt(length(vecHetFastMiss)),'kx')
hold off
set(gca,'XTick',1:3,'XTickLabel',{'S-M','S-F','F-M'})
ylabel(sprintf('Separability \n(inter/(intra+inter) cluster distance)'))
[h,pActSM] = ttest(vecActSlowMiss,0.5);
[h,pActSF] = ttest(vecActSlowFast,0.5);
[h,pActFM] = ttest(vecActFastMiss,0.5);
[h,pHetSM] = ttest(vecHetSlowMiss,0.5);
[h,pHetSF] = ttest(vecHetSlowFast,0.5);
[h,pHetFM] = ttest(vecHetFastMiss,0.5);
[h,pSM] = ttest(vecActSlowMiss,vecHetSlowMiss);
[h,pSF] = ttest(vecActSlowFast,vecHetSlowFast);
[h,pFM] = ttest(vecActFastMiss,vecHetFastMiss);
title(sprintf(['Stim; Sep act; SM,p=%.3f; SF,p=%.3f; FM,p=%.3f\n'...
	'Sep het; SM,p=%.3f; SF,p=%.3f; FM,p=%.3f\n',...
	'Diff act-het; SM,p=%.3f; SF,p=%.3f; FM,p=%.3f'],...
	pActSM,pActSF,pActFM,...
	pHetSM,pHetSF,pHetFM,...
	pSM,pSF,pFM));
ylim([0 1])
xlim([0.5 3.5])
legend('dF/F','Heterogeneity')




% NOW PERFORM FOR PROBE TRIALS
vecHetMiss = [];
vecHetSlow = [];
vecHetFast = [];
vecActMiss = [];
vecActSlow = [];
vecActFast = [];

%get data
for intAnimal=1:size(cellHetTime,2)
	matTemp = [];
	for intC=1
		matTemp = [matTemp;cellHetTime{intC,intAnimal}{1};cellHetTime{intC,intAnimal}{2};cellHetTime{intC,intAnimal}{3}];
	end
	dblMean = mean(matTemp(:)); %mean over all values of animal
	
	%apply normalization
	for intC=1
		cellHetTime{intC,intAnimal}{1} = cellHetTime{intC,intAnimal}{1} ./ dblMean;
		cellHetTime{intC,intAnimal}{2} = cellHetTime{intC,intAnimal}{2} ./ dblMean;
		cellHetTime{intC,intAnimal}{3} = cellHetTime{intC,intAnimal}{3} ./ dblMean;
	end
	
	%prep data
	matTempHetM = [];
	matTempHetS = [];
	matTempHetF = [];
	matTempActM = [];
	matTempActS = [];
	matTempActF = [];
	for intC=1
		matTempHetM = [matTempHetM;cellHetTime{intC,intAnimal}{1}];
		matTempHetS = [matTempHetS;cellHetTime{intC,intAnimal}{2}];
		matTempHetF = [matTempHetF;cellHetTime{intC,intAnimal}{3}];
		matTempActM = [matTempActM;cellActTime{intC,intAnimal}{1}];
		matTempActS = [matTempActS;cellActTime{intC,intAnimal}{2}];
		matTempActF = [matTempActF;cellActTime{intC,intAnimal}{3}];
	end
	vecHetMiss(intAnimal) = mean(mean(matTempHetM(:,intStart:intStop),1));
	vecHetSlow(intAnimal) = mean(mean(matTempHetS(:,intStart:intStop),1));
	vecHetFast(intAnimal) = mean(mean(matTempHetF(:,intStart:intStop),1));
	vecActMiss(intAnimal) = mean(mean(matTempActM(:,intStart:intStop),1));
	vecActSlow(intAnimal) = mean(mean(matTempActS(:,intStart:intStop),1));
	vecActFast(intAnimal) = mean(mean(matTempActF(:,intStart:intStop),1));
end
% check separability
%{
calculate per point intra-cluster vs inter-cluster difference (ratio) for
dF/F and for heterogeneity dimension; will give distribution of
separability for each cluster; then perform ANOVA/ttests for
between-cluster distribution comparison
%}

%get cluster means
dblMeanActMiss = mean(vecActMiss);
dblMeanActSlow = mean(vecActSlow);
dblMeanActFast = mean(vecActFast);
dblMeanHetMiss = mean(vecHetMiss);
dblMeanHetSlow = mean(vecHetSlow);
dblMeanHetFast = mean(vecHetFast);
intAnimals = length(vecActMiss);

%separability; inter vs intra cluster distance
vecSepRatioActMiss_w_Slow = abs(vecActMiss - dblMeanActSlow) ./ (abs(vecActMiss - dblMeanActMiss) + abs(vecActMiss - dblMeanActSlow));%inter / (inter+intra) distance
vecSepRatioActMiss_w_Fast = abs(vecActMiss - dblMeanActFast) ./ (abs(vecActMiss - dblMeanActMiss) + abs(vecActMiss - dblMeanActFast));
vecSepRatioActSlow_w_Miss = abs(vecActSlow - dblMeanActMiss) ./ (abs(vecActSlow - dblMeanActSlow) + abs(vecActSlow - dblMeanActMiss));
vecSepRatioActSlow_w_Fast = abs(vecActSlow - dblMeanActFast) ./ (abs(vecActSlow - dblMeanActSlow) + abs(vecActSlow - dblMeanActFast));
vecSepRatioActFast_w_Slow = abs(vecActFast - dblMeanActSlow) ./ (abs(vecActFast - dblMeanActFast) + abs(vecActFast - dblMeanActSlow));
vecSepRatioActFast_w_Miss = abs(vecActFast - dblMeanActMiss) ./ (abs(vecActFast - dblMeanActFast) + abs(vecActFast - dblMeanActMiss));
vecSepRatioHetMiss_w_Slow = abs(vecHetMiss - dblMeanHetSlow) ./ (abs(vecHetMiss - dblMeanHetMiss) + abs(vecHetMiss - dblMeanHetSlow));
vecSepRatioHetMiss_w_Fast = abs(vecHetMiss - dblMeanHetFast) ./ (abs(vecHetMiss - dblMeanHetMiss) + abs(vecHetMiss - dblMeanHetFast));
vecSepRatioHetSlow_w_Miss = abs(vecHetSlow - dblMeanHetMiss) ./ (abs(vecHetSlow - dblMeanHetSlow) + abs(vecHetSlow - dblMeanHetMiss));
vecSepRatioHetSlow_w_Fast = abs(vecHetSlow - dblMeanHetFast) ./ (abs(vecHetSlow - dblMeanHetSlow) + abs(vecHetSlow - dblMeanHetFast));
vecSepRatioHetFast_w_Slow = abs(vecHetFast - dblMeanHetSlow) ./ (abs(vecHetFast - dblMeanHetFast) + abs(vecHetFast - dblMeanHetSlow));
vecSepRatioHetFast_w_Miss = abs(vecHetFast - dblMeanHetMiss) ./ (abs(vecHetFast - dblMeanHetFast) + abs(vecHetFast - dblMeanHetMiss));

%act agg
vecActSlowMiss = [vecSepRatioActMiss_w_Slow vecSepRatioActSlow_w_Miss];
vecActSlowFast = [vecSepRatioActFast_w_Slow vecSepRatioActSlow_w_Fast];
vecActFastMiss = [vecSepRatioActMiss_w_Fast vecSepRatioActFast_w_Miss];

%het agg
vecHetSlowMiss = [vecSepRatioHetMiss_w_Slow vecSepRatioHetSlow_w_Miss];
vecHetSlowFast = [vecSepRatioHetFast_w_Slow vecSepRatioHetSlow_w_Fast];
vecHetFastMiss = [vecSepRatioHetMiss_w_Fast vecSepRatioHetFast_w_Miss];

%plot
subplot(2,2,3)
scatter(vecHetMiss,vecActMiss,'rx')
hold on
scatter(vecHetSlow,vecActSlow,'yx')
scatter(vecHetFast,vecActFast,'gx')
hold off
title('Probe; Response type predictability')
xlabel('Heterogeneity during 1s preceding stimulus')
ylabel('dF/F during 1s preceding stimulus')

subplot(2,2,4)

dblOffset = 0.1;
plot([0.5 3.5],[0.5 0.5],'k--')
%plot act
hold on
errorbar(1-dblOffset,mean(vecActSlowMiss),std(vecActSlowMiss)/sqrt(length(vecActSlowMiss)),'bx')
errorbar(2-dblOffset,mean(vecActSlowFast),std(vecActSlowFast)/sqrt(length(vecActSlowFast)),'bx')
errorbar(3-dblOffset,mean(vecActFastMiss),std(vecActFastMiss)/sqrt(length(vecActFastMiss)),'bx')
%plot het
errorbar(1+dblOffset,mean(vecHetSlowMiss),std(vecHetSlowMiss)/sqrt(length(vecHetSlowMiss)),'kx')
errorbar(2+dblOffset,mean(vecHetSlowFast),std(vecHetSlowFast)/sqrt(length(vecHetSlowFast)),'kx')
errorbar(3+dblOffset,mean(vecHetFastMiss),std(vecHetFastMiss)/sqrt(length(vecHetFastMiss)),'kx')
hold off
set(gca,'XTick',1:3,'XTickLabel',{'S-M','S-F','F-M'})
ylabel(sprintf('Separability \n(inter/(intra+inter) cluster distance)'))
[h,pActSM] = ttest(vecActSlowMiss,0.5);
[h,pActSF] = ttest(vecActSlowFast,0.5);
[h,pActFM] = ttest(vecActFastMiss,0.5);
[h,pHetSM] = ttest(vecHetSlowMiss,0.5);
[h,pHetSF] = ttest(vecHetSlowFast,0.5);
[h,pHetFM] = ttest(vecHetFastMiss,0.5);
[h,pSM] = ttest(vecActSlowMiss,vecHetSlowMiss);
[h,pSF] = ttest(vecActSlowFast,vecHetSlowFast);
[h,pFM] = ttest(vecActFastMiss,vecHetFastMiss);
title(sprintf(['Probe; Sep act; SM,p=%.3f; SF,p=%.3f; FM,p=%.3f\n'...
	'Sep het; SM,p=%.3f; SF,p=%.3f; FM,p=%.3f\n',...
	'Diff act-het; SM,p=%.3f; SF,p=%.3f; FM,p=%.3f'],...
	pActSM,pActSF,pActFM,...
	pHetSM,pHetSF,pHetFM,...
	pSM,pSF,pFM));
ylim([0 1])
xlim([0.5 3.5])
legend('dF/F','Heterogeneity')

if boolSavePlots
	strFigTitle = 'RespTypePredictability';
	intFigCounter = intFigCounter + 1;
	drawnow;
	jFig = get(handle(gcf), 'JavaFrame');
	jFig.setMaximized(true);
	figure(gcf);
	drawnow;
	strFig = sprintf('Meta%d_%s_raw',intFigCounter,strFigTitle);
	export_fig([strFig '.tif']);
	export_fig([strFig '.pdf']);
end


%% resp type decoding
hFigRTD = figure;
for intRespType=1:3
	if intRespType == 1
		strTitle = 'Miss';
		cellColor = {'r','k','k'};
	elseif intRespType == 2
		strTitle = 'Slow';
		cellColor = {'k','y','k'};
	elseif intRespType == 3
		strTitle = 'Fast';
		cellColor = {'k','k','g'};
	end
	subplot(2,2,intRespType)
	
	%make polar background
	[x,y] = pol2cart((4/3)*pi,0.9);%miss
	text(x,y,'Miss','Color',cellColor{1});
	hold on
	[x,y] = pol2cart((0/3)*pi,0.7);%slow
	text(x,y,'Slow','Color',cellColor{2});
	[x,y] = pol2cart((2/3)*pi,0.9);%fast
	text(x,y,'Fast','Color',cellColor{3});
	[x,y] = pol2cart((1/3)*pi,1);%slow/fast
	plot([0 x],[0 y],'k');
	[x,y] = pol2cart((3/3)*pi,1);%miss/fast
	plot([0 x],[0 y],'k');
	[x,y] = pol2cart((5/3)*pi,1);%slow/miss
	plot([0 x],[0 y],'k');
	ang=0:0.001:2*pi;
	x=cos(ang);
	y=sin(ang);
	plot(x,y, 'k');
	
	%plot data
	for intAnimal=1:size(matRespDecode,4)
		thetaAct = matRespDecode(1,1,intRespType,intAnimal);
		rhoAct = matRespDecode(1,2,intRespType,intAnimal);
		thetaHet = matRespDecode(2,1,intRespType,intAnimal);
		rhoHet = matRespDecode(2,2,intRespType,intAnimal);
		
		[x,y]=pol2cart(thetaAct,rhoAct);
		scatter(x,y,'bx');
		[x,y]=pol2cart(thetaHet,rhoHet);
		scatter(x,y,'kx');
	end
	hold off;
	xlim([-1 1]);
	ylim([-1 1]);
	title(sprintf('%s trials',strTitle))
end

%plot summary
subplot(2,2,4);
matPerformance = 1-abs(matRespDecodeDist);
vecActPerfomance = matPerformance(1,:);
vecHetPerfomance = matPerformance(2,:);

%statistics
[hH,pH]=ttest(vecHetPerfomance);
[hA,pA]=ttest(vecActPerfomance);
[hHA,pHA]=ttest(vecActPerfomance,vecHetPerfomance);
vecSubHetP = zeros(1,3);
vecSubActP = zeros(1,3);
vecSubDiffP = zeros(1,3);

%plot
dblOffset=0.1;
dblMeanAct=mean(vecActPerfomance);
dblSDAct = std(vecActPerfomance);
errorbar((1+dblOffset),dblMeanAct,dblSDAct/sqrt(length(vecActPerfomance)),'xb')
hold on
dblMeanHet=mean(vecHetPerfomance);
dblSDHet = std(vecHetPerfomance);
errorbar((2+dblOffset),dblMeanHet,dblSDHet/sqrt(length(vecHetPerfomance)),'xk')
plot([0.5 2.5],[0 0],'k--')
for intRespType=1:3
	if intRespType == 1
		strColor = 'r';
	elseif intRespType == 2
		strColor = 'y';
	elseif intRespType == 3
		strColor = 'g';
	end
	vecThisAct = vecActPerfomance(intRespType:3:end);
	vecThisHet = vecHetPerfomance(intRespType:3:end);
	scatter((1-dblOffset)*ones(size(vecThisAct)),vecThisAct,[strColor 'x']);
	scatter((2-dblOffset)*ones(size(vecThisHet)),vecThisHet,[strColor 'x']);
	
	[dummy,vecSubHetP(intRespType)] = ttest(vecThisHet);
	[dummy,vecSubActP(intRespType)] = ttest(vecThisAct);
	[dummy,vecSubDiffP(intRespType)] = ttest(vecThisAct,vecThisHet);
end
hold off
title(sprintf('Normalized resp type decoding; pA=%.3f;pH=%.3f;pHA=%.3f;\nM:pA=%.3f;pH=%.3f;pHA=%.3f S:pA=%.3f;pH=%.3f;pHA=%.3f F:pA=%.3f;pH=%.3f;pHA=%.3f',...
	pA,pH,pHA,vecSubActP(1),vecSubHetP(1),vecSubDiffP(1),vecSubActP(2),vecSubHetP(2),vecSubDiffP(2),vecSubActP(3),vecSubHetP(3),vecSubDiffP(3)))
ylabel('Norm decod perf')
set(gca,'XTick',[1 2],'XTickLabel',{'dF/F0', 'Heterogen'})
ylim([-1 1])
xlim([0.5 2.5])
if boolSavePlots
	strFigTitle = 'RespTypeDecoding';
	intFigCounter = intFigCounter + 1;
	drawnow;
	strFig = sprintf('Meta%d_%s_raw',intFigCounter,strFigTitle);
	export_fig([strFig '.tif']);
	export_fig([strFig '.pdf']);
end


%% difference fast/slow response trials for dF/F + heterogeneity over time
%intAct=1 = het; intAct=2 = dF/F
%matDiffFS(intAct,intC,intAnimal,intT)


%vars
vecC = [0 0.5 2 8 32 100];
vecWindowSecs = [-3 5];
vecWindow = round(vecWindowSecs*25.4);
vecLineX = (vecWindow(1):vecWindow(end))/25.4;
vecWindowInv = length(vecLineX):-1:1;
vecX = [vecLineX vecLineX(vecWindowInv)];

%loop
for intAct=1:2
	if intAct == 1 %miss
		strType = 'Heterogeneity';
		vecColorLine = [0 0 0];
		vecColorFill = [0.7 0.7 0.7];
	else %fast
		strType = 'dF/F0';
		vecColorLine = [0 0 1];
		vecColorFill = [0.7 0.7 1];
	end
	
	%plot
	hFigDFST = figure;
	for intC=1:length(vecC)
		subplot(2,3,intC);
		
		%pre-allocate
		intMaxT = size(matDiffFS,4);
		vecSignificant = zeros(1,intMaxT);
		vecP = zeros(1,intMaxT);
		matCI = zeros(2,intMaxT);
		for intT=1:intMaxT
			%get data
			vecData = matDiffFS(intAct,intC,:,intT);
			
			%test
			[boolH,dblP,vecCI] = ttest(vecData);
			matCI(:,intT) = vecCI;
			vecP(intT) = dblP;
			vecSignificant(intT) = boolH;
		end
		vecMeanTrace = squeeze(mean(matDiffFS(intAct,intC,:,:),3));
		
		%plot curves
		vecMinTrace = vecMeanTrace - squeeze(std(matDiffFS(intAct,intC,:,:),[],3))/sqrt(size(matDiffFS,3));
		vecMaxTrace = vecMeanTrace + squeeze(std(matDiffFS(intAct,intC,:,:),[],3))/sqrt(size(matDiffFS,3));
		vecY = [vecMinTrace' vecMaxTrace(vecWindowInv)'];
		
		%plot
		hold on
		fill(vecX,vecY,vecColorFill,'EdgeColor','none');
		plot(vecLineX,vecMeanTrace,'-','LineWidth',2,'Color',vecColorLine);
		
		%plot significance
		vecLimY = get(gca,'YLim');
		dblPlotY = vecLimY(2) - (vecLimY(2) - vecLimY(1))*(1/20);
		boolLastPlotted = false;
		intCounterT = 1;
		while ~boolLastPlotted
			if vecSignificant(1)
				intNextT = find(vecSignificant==0,1,'first');
				if isempty(intNextT),boolLastPlotted=true;break;end
				plot([vecLineX(intCounterT) vecLineX(intNextT+intCounterT-1)],[dblPlotY dblPlotY]);
			else
				intNextT = find(vecSignificant==1,1,'first');
				if isempty(intNextT),boolLastPlotted=true;break;end
			end
			
			%update counters
			intCounterT = intCounterT + intNextT -1;
			vecSignificant = vecSignificant(intNextT:end);
		end
		hold off
		
		%labels
		if intAct == 1
			title(sprintf('Diff Fast-slow heterogeneity over time; Contrast %.1f',vecC(intC)))
			ylabel('Difference in heterogeneity')
			ylim([-0.1 0.4])
			strFigTitle = 'diffFS_hetero_over_time';
		else
			title(sprintf('Diff Fast-slow dF/F0 over time; Contrast %.1f',vecC(intC)))
			ylabel('Difference in dF/F0')
			ylim([-0.01 0.04])
			strFigTitle = 'diffFS_activity_over_time';
		end
		grid on
		xlabel('Time after stim onset (s)')
		xlim(vecWindowSecs)
		drawnow;
	end
	
	if boolSavePlots
		intFigCounter = intFigCounter + 1;
		drawnow;
		jFig = get(handle(gcf), 'JavaFrame');
		jFig.setMaximized(true);
		figure(gcf);
		drawnow;
		strFig = sprintf('Meta%d_%s_raw',intFigCounter,strFigTitle);
		export_fig([strFig '.tif']);
		export_fig([strFig '.pdf']);
	end
end

%% plot heterogeneity vs correlation ideal/actual pop response (pref stim)
hHetCorrIdActPopRe = figure;
vecX = [0 3];
vecSlopes = nan(1,intAnimals);
vecIntercepts = nan(1,intAnimals);
vecHetPoints = [];
vecPrefRespCorrPoints = [];
hold on
cMap=colormap(lines(intAnimals));
for intAnimal=1:intAnimals
	vecSlopes(intAnimal) = cellHetPoprespCorr{3,intAnimal}(2);
	vecIntercepts(intAnimal) = cellHetPoprespCorr{3,intAnimal}(1);
	
	vecY = polyval([vecSlopes(intAnimal) vecIntercepts(intAnimal)],vecX);
	plot(vecX,vecY,'Color',[0.5 0.5 0.5]);
	
	vecHetPoints = [vecHetPoints; cellHetPoprespCorr{1,intAnimal}];
	vecPrefRespCorrPoints = [vecPrefRespCorrPoints; cellHetPoprespCorr{2,intAnimal}];
	
	scatter(cellHetPoprespCorr{1,intAnimal},cellHetPoprespCorr{2,intAnimal},[],cMap(intAnimal,:),'x');
end
dblMeanSlope = mean(vecSlopes);
dblMeanIntercept = mean(vecIntercepts);
vecY = polyval([dblMeanSlope dblMeanIntercept],vecX);
plot(vecX,vecY,'b');
hold off
[h,p] = ttest(vecSlopes);
title(sprintf('T-test slopes; mean=%.3f; p=%.3f;',dblMeanSlope,p));

xlabel('Heterogeneity')
ylabel('Correlation ideal/actual population response')

ylim([-0.8 0.8])
		

strFigTitle = 'hetcorridealactualpopresp';
if boolSavePlots
	intFigCounter = intFigCounter + 1;
	figure(hHetCorrIdActPopRe);
	drawnow;
	strFig = sprintf('Meta%d_%s_raw',intFigCounter,strFigTitle);
	export_fig([strFig '.tif']);
	export_fig([strFig '.pdf']);
end

%% clean up
cd(strOldDir);