%close all
clear all
%% parameters
strFigDir = 'E:\UvA_Backup\Data\Results\spikeAnalysis\meta2';
strDataDir = 'E:\UvA_Backup\Data\Results\spikeAnalysis';
cellInclude = {'20140207','20140314','20140425','20140507','20140530','20140604','20140711','20140715'};%exclude:20140430
boolSavePlots = false;
boolUseNPS = false;
boolOnlyTuned = false;
intFigCounter = 0;

%% pre-allocate
intAnimal = 0;
vecAnimalSource = [];
vecSubPopulationSource = [];

matTempSeqStab = [];
matDecodingAnimals = [];
matAccOff = [];
matPowerSpectraSelf = [];
matPowerSpectraOthers = [];
matOccurrenceT = [];
matOccurrenceSelf = [];
matOccurrenceOther = [];
matOccurrenceStim = [];
vecOccurrenceBase = [];
matOccurrenceOffset = [];
vecOccurrenceOffsetBase = [];
matHMMDecoding = [];
matOriSimOrder = [];
matOriSimOverlap = [];
matRespSimOrder = [];
matRespSimOverlap = [];
matOriSimOrderS = [];
matOriSimOverlapS = [];
matRespSimOrderS = [];
matRespSimOverlapS = [];
vecHitPresRatioZ = [];
vecMissPresRatioZ = [];
		
%2
cellClusterZ = [];
		
if boolUseNPS,strNPS='NPS';else strNPS='';end
if boolOnlyTuned,strOT='OT';else strOT='';end

%% create aggregate data structures
sDir = dir(strDataDir);
strOldDir = cd(strDataDir);
intCounterF1 = 0;
intCounterF2 = 0;
intCounterF3 = 0;
intCounterF4 = 0;
vecNumClusts = [];
for intFile=1:numel(sDir)
	%check if data part 1
	strFile = sDir(intFile).name;
	strRec = getFlankedBy(strFile,['dataAssAnalFiveAA' strNPS],[strOT '_']);
	intRec = find(strcmp(strRec,cellInclude),1);
	if ~isempty(intRec)
		%load file
		load(strFile);
		
		%update population source
		intNrPops = size(cellSaveEnsembleSequences,1);
		if intNrPops > 1 %group populations
			if exist('cellSaveTemporalSequenceStability','var')
				cellSaveTemporalSequenceStability{1} = mean(cat(3,cellSaveTemporalSequenceStability{1},cellSaveTemporalSequenceStability{2}),3);
			end
			if exist('cellSaveOverallOriRespDecoding','var')
				cellSaveOverallOriRespDecoding{1} = mean(cat(3,cellSaveOverallOriRespDecoding{1},cellSaveOverallOriRespDecoding{2}),3);
				cellSaveOverallOriRespDecoding(2) = [];
			end
			
			if exist('cellSaveDecodingAT','var')
				cellSaveDecodingAT{1,4} = mean(cat(3,cellSaveDecodingAT{1,4},cellSaveDecodingAT{2,4}),3);
			end
			if exist('cellSaveConsistencyPerTrial','var')
				for intProp=1:size(cellSaveConsistencyPerTrial,2)
					cellSaveConsistencyPerTrial{1,intProp} = [cellSaveConsistencyPerTrial{1,intProp} cellSaveConsistencyPerTrial{2,intProp}];
				end
			end
			for intProp=1:size(cellSaveEnsembleSequences,2)
				cellSaveEnsembleSequences{1,intProp} = nanmean(cat(3,cellSaveEnsembleSequences{1,intProp},cellSaveEnsembleSequences{2,intProp}),3);
			end
			cellSaveEnsembleSequences(2,:) = [];
			if exist('cellSavePresRatioZ','var')
				cellSavePresRatioZ{1,1} = [cellSavePresRatioZ{1,1}; cellSavePresRatioZ{2,1}];
				cellSavePresRatioZ{1,2} = [cellSavePresRatioZ{1,2}; cellSavePresRatioZ{2,2}];
				cellSavePresRatioZ(2,:) = [];
			end
		end
		if exist('cellSavePresRatioZ','var')
			vecHitPresRatioZ = [vecHitPresRatioZ; cellSavePresRatioZ{1,1}];
			vecMissPresRatioZ = [vecMissPresRatioZ; cellSavePresRatioZ{1,2}];
		end
		if exist('cellSaveConsistencyPerTrial','var')
			matConsistencyTemp = cellfun(@mean,cellSaveConsistencyPerTrial);
			matDecodingAnimals(end+1,:) = matConsistencyTemp(1,:);
		end
		
		%HMM
		if exist('cellSaveOverallOriRespDecoding','var')
			matHMMDecoding(:,:,intCounterF3+1) = cellSaveOverallOriRespDecoding{1};
		end
		
		%cluster sequences
		matOriSimOrder(:,:,intCounterF3+1) = cellSaveEnsembleSequences{1,1};
		matOriSimOverlap(:,:,intCounterF3+1) = cellSaveEnsembleSequences{1,2};
		matRespSimOrder(:,:,intCounterF3+1) = cellSaveEnsembleSequences{1,3};
		matRespSimOverlap(:,:,intCounterF3+1) = cellSaveEnsembleSequences{1,4};
		matOriSimOrderS(:,:,intCounterF3+1) = cellSaveEnsembleSequences{1,5};
		matOriSimOverlapS(:,:,intCounterF3+1) = cellSaveEnsembleSequences{1,6};
		matRespSimOrderS(:,:,intCounterF3+1) = cellSaveEnsembleSequences{1,7};
		matRespSimOverlapS(:,:,intCounterF3+1) = cellSaveEnsembleSequences{1,8};
		
		
		%occurrences in time
		if exist('cellSaveAssemblyRecurrenceMinDur','var')
			cellCat = cellSaveAssemblyRecurrenceMinDur(:,:,4);
			matOccurrenceT = cat(1,matOccurrenceT,cell2mat(cellCat(:)));
			cellCat = cellSaveAssemblyRecurrenceMinDur(:,:,5);
			matOccurrenceSelf = cat(1,matOccurrenceSelf,cell2mat(cellCat(:)));
			cellCat = cellSaveAssemblyRecurrenceMinDur(:,:,6);
			matOccurrenceOther = cat(1,matOccurrenceOther,cell2mat(cellCat(:)));
			cellCat = cellSaveAssemblyRecurrenceMinDur(:,:,7);
			matOccurrenceStim = cat(1,matOccurrenceStim,cell2mat(cellCat(:)));
			cellCat = cellSaveAssemblyRecurrenceMinDur(:,:,8);
			vecOccurrenceBase = cat(1,vecOccurrenceBase,cell2mat(cellCat(:)));
			cellCat = cellSaveAssemblyRecurrenceMinDur(:,:,9);
			matOccurrenceOffset = cat(1,matOccurrenceOffset,cell2mat(cellCat(:)));
			cellCat = cellSaveAssemblyRecurrenceMinDur(:,:,10);
			vecOccurrenceOffsetBase = cat(1,vecOccurrenceOffsetBase,cell2mat(cellCat(:)));
		end
		
		%save data
		if exist('cellSaveDecodingAT','var')
			vecTimeOff = cellSaveDecodingAT{1,3};
			matAccOff(end+1,:) = cellSaveDecodingAT{1,4};
		end
		
		%assign data
		if exist('cellSaveTemporalSequenceStability','var')
			matTempSeqStab = [matTempSeqStab; cellSaveTemporalSequenceStability{1}];
		end
		
		intCounterF3 = intCounterF3 + 1;
	end
	
	%check if data part 2
	strFile = sDir(intFile).name;
	strRec = getFlankedBy(strFile,['dataAssAnalSixAA' strNPS],[strOT '_']);
	intRec = find(strcmp(strRec,cellInclude),1);
	if ~isempty(intRec)
		%load file
		load(strFile);
		
		%update population source
		intNrPops = size(cellSaveClusterZ,1);
		if intNrPops > 1 %group populations
			for intMinSize=1:size(cellSaveClusterZ,2)
				for intType=1:size(cellSaveClusterZ,3)
					cellSaveClusterZ{1,intMinSize,intType} = [cellSaveClusterZ{1,intMinSize,intType} cellSaveClusterZ{2,intMinSize,intType}];
				end
			end
			cellSaveClusterZ(2,:,:) = [];
		end
		cellSaveClusterZ = squeeze(cellSaveClusterZ);
		
		%concatenate
		if isempty(cellClusterZ)
			cellClusterZ=cellSaveClusterZ;
		else
			cellClusterZ = cellfun(@cat,cellfill(2,size(cellSaveClusterZ)),cellClusterZ,cellSaveClusterZ,'UniformOutput',false);
		end
		
		intCounterF2 = intCounterF2 + 1;
	end

end

%check for inconsistencies
%if ~((intCounterF1 == intCounterF2)% && (intCounterF2 == intCounterF3) && (intCounterF3 == intCounterF4) && (intCounterF4 == intCounterF5))
%	warning([mfilename ':InconsistentFileNumbers'],'File counters inconsistent; F1=%d; F2=%d',...; F3=%d; F4=%d; F5=%d',
%		intCounterF1,intCounterF2);%,intCounterF3,intCounterF4,intCounterF5)
%end

%% meta analyses
cd(strFigDir);
close all;
intAnimals = size(matTempSeqStab,1);
vecContrasts = [0 0.5 2 8 32 100];

%% cluster presence
dblStep = 1;
vecHistBins = -(4-dblStep/2):dblStep:(4-dblStep/2);

vecDistHit = hist(vecHitPresRatioZ,vecHistBins);
vecDistMiss = hist(vecMissPresRatioZ,vecHistBins);
stairs([vecHistBins-dblStep/2 4],[vecDistMiss vecDistMiss(end)],'color','r');
hold on
stairs([vecHistBins-dblStep/2 4],[vecDistHit vecDistHit(end)],'color','g');
hold off
xlim([-4 4])
[h,pT]=ttest2(vecHitPresRatioZ,vecMissPresRatioZ);
[h,pF]=vartest2(vecHitPresRatioZ,vecMissPresRatioZ);
title(sprintf('Hit (mu=%.3f,sd=%.3f) vs Miss (mu=%.3f,sd=%.3f), t-test, p=%.3f; F-test, p=%.3f',mean(vecHitPresRatioZ),std(vecHitPresRatioZ),mean(vecMissPresRatioZ),std(vecMissPresRatioZ),pT,pF));
xlabel('Mean trial presence ratio normalized to shuffled distro (sd)');
ylabel('Nr of clusters (count)');

if boolSavePlots
	drawnow;
	intFigCounter = intFigCounter + 1;
	strFig = sprintf('Meta%s%d_trial_pres_ratio_raw',strNPS,intFigCounter);
	export_fig([strFig '.tif']);
	export_fig([strFig '.pdf']);
end

%% recurring cluster properties
vecHistBins = -3.5:1:3.5;

%{
cellDistro = cellfun(@hist,cellClusterZ,cellfill(vecHistBins,size(cellClusterZ)),'UniformOutput',false);
matDistroDistZ = bsxfun(@rdivide,cell2mat(cellDistro(:,1)),sum(cell2mat(cellDistro(:,1)),2));
[h,matP] = cellfun(@ttest,cellClusterZ);
matMean = cellfun(@nanmean,cellClusterZ);
matClusterN = cellfun(@numel,cellClusterZ)-5;

matDistrodPAZ = bsxfun(@rdivide,cell2mat(cellDistro(:,2)),sum(cell2mat(cellDistro(:,2)),2));
matDistroOSIZ = bsxfun(@rdivide,cell2mat(cellDistro(:,3)),sum(cell2mat(cellDistro(:,3)),2));
%}

vecDistZ = cellClusterZ{2,1};
vecOSIZ = cellClusterZ{2,3};
vecdPAZ = cellClusterZ{2,2};

figure
subplot(2,2,2);
hist(vecDistZ,vecHistBins)
xlim([-4 4])
[h,p]=ttest(vecDistZ);
title(sprintf('Inter-member distance, t-test vs 0, p=%.3f; mean=%.3f sd; frac >2sd: %.3f',p,mean(vecDistZ),sum(abs(vecDistZ)>2)/numel(vecDistZ)));
xlabel('Mean distance normalized to shuffled distro (sd)');
ylabel('Nr of recurring clusters (count)');

subplot(2,2,3);
hist(vecOSIZ,vecHistBins)
xlim([-4 4])
[h,p]=ttest(vecOSIZ);
title(sprintf('OSI clusters, t-test vs 0, p=%.3f; mean=%.3f sd; frac >2sd: %.3f',p,mean(vecOSIZ),sum(abs(vecOSIZ)>2)/numel(vecOSIZ)));
xlabel('Mean OSI normalized to shuffled distro (sd)');
ylabel('Nr of recurring clusters (count)');

subplot(2,2,4);
hist(vecdPAZ,vecHistBins)
xlim([-4 4])
[h,p]=ttest(vecdPAZ);
title(sprintf('Pref ori diff assembly members, t-test vs 0, p=%.3f; mean=%.3f sd; frac >2sd: %.3f',p,mean(vecdPAZ),sum(abs(vecdPAZ)>2)/numel(vecdPAZ)));
xlabel('Mean dPO normalized to shuffled distro (sd)');
ylabel('Nr of recurring clusters (count)');



%% ensemble sequences
figure
%orientation-based
subplot(2,2,1)
imagesc(nanmean(matOriSimOrder,3),[0.1 0.3]);
colormap(hot(256));
subplot(2,2,2)
imagesc(nanmean(matOriSimOrderS,3),[0.1 0.3]);
colormap(hot(256));


vecSameOri = [];
vecOtherOri = [];

vecSameOriS = [];
vecOtherOriS = [];
for intAnimal=1:8
	matT = matOriSimOrder(:,:,intAnimal);
	vecSameOri(intAnimal) = mean(matT(tril(triu(true(4)))));
	vecOtherOri(intAnimal) = mean(matT(tril(true(4),-1)));
	
	matT = matOriSimOrderS(:,:,intAnimal);
	vecSameOriS(intAnimal) = mean(matT(tril(triu(true(4)))));
	vecOtherOriS(intAnimal) = mean(matT(tril(true(4),-1)));
	
end

%t-tests
[h,pSameShuffled]=ttest2(vecSameOri,vecSameOriS);
[h,pOtherShuffled]=ttest2(vecOtherOri,vecOtherOriS);
[h,pSameOther]=ttest2(vecSameOri,vecOtherOri);

vecMeansO = [mean(vecOtherOriS) mean(vecOtherOri) mean(vecSameOriS) mean(vecSameOri)];
vecSEMsO = [std(vecOtherOriS)/sqrt(numel(vecOtherOriS)) std(vecOtherOri)/sqrt(numel(vecOtherOri)) std(vecSameOriS)/sqrt(numel(vecSameOriS)) std(vecSameOri)/sqrt(numel(vecSameOri))];

subplot(2,2,3)
errorbar(1:4,vecMeansO,vecSEMsO,'xb');
set(gca,'xtick',1:4,'xticklabel',{'Other shuf','Other','Same Shuf','Same'});
title(sprintf('T-tests, other-shuf,p=%.3f; same-shuf,p=%.3f;same-other,p=%.3f',pOtherShuffled,pSameShuffled,pSameOther))
ylim([0.1 0.3+eps])
ylabel('Sequence similarity')

%stim detection
vecSeqMissesShuffled = matRespSimOrderS(1,1,:);
vecSeqMisses = matRespSimOrder(1,1,:);
vecSeqHitsShuffled = matRespSimOrderS(2,2,:);
vecSeqHits = matRespSimOrder(2,2,:);

vecMeans = [mean(vecSeqMissesShuffled) mean(vecSeqMisses) mean(vecSeqHitsShuffled) mean(vecSeqHits)];
vecSEMs = [std(vecSeqMissesShuffled)/sqrt(numel(vecSeqMissesShuffled)) std(vecSeqMisses)/sqrt(numel(vecSeqMisses)) std(vecSeqHitsShuffled)/sqrt(numel(vecSeqHitsShuffled)) std(vecSeqHits)/sqrt(numel(vecSeqHits))];

%t-tests
[h,pMissShuffled]=ttest(vecSeqMissesShuffled,vecSeqMisses);
[h,pHitShuffled]=ttest(vecSeqHitsShuffled,vecSeqHits);
[h,pMissHit]=ttest(vecSeqHits,vecSeqMisses);

subplot(2,2,4)
errorbar(1:4,vecMeans,vecSEMs,'xb');
set(gca,'xtick',1:4,'xticklabel',{'Miss shuf','Miss','Hit Shuf','Hit'});
title(sprintf('T-tests, miss-shuf,p=%.3f; hit-shuf,p=%.3f;hit-miss,p=%.3f',pMissShuffled,pHitShuffled,pMissHit))
ylim([0.1 0.3+eps])
ylabel('Sequence similarity')

%%
if boolSavePlots
	drawnow;
	jFig = get(handle(gcf), 'JavaFrame');
	jFig.setMaximized(true);
	figure(gcf);
	drawnow;
	intFigCounter = intFigCounter + 1;
	strFig = sprintf('Meta%s%d_ensemble_seq_hitmiss_raw',strNPS,intFigCounter);
	export_fig([strFig '.tif']);
	export_fig([strFig '.pdf']);
end

%% HMM decoding
matMean = mean(matHMMDecoding,3);
matSEM = std(matHMMDecoding,[],3)/sqrt(8);
[h,p]=ttest(matHMMDecoding,0.25,[],[],3);
vecPlotC = vecContrasts;
vecPlotC(1) = 0.2;
vecHits = matHMMDecoding(2,2:6,:);
vecMiss = matHMMDecoding(1,2:6,:);
[h,pHC]=ttest(vecHits(:),0.25);
[h,pMC]=ttest(vecMiss(:),0.25);
[h,pHM]=ttest(vecMiss(:),vecHits(:));
errorfill(vecPlotC,matMean(1,:),matSEM(1,:),[0.7 0 0],[1 0.5 0.5])
hold on
errorfill(vecPlotC,matMean(2,:),matSEM(2,:),[0 0.7 0],[0.5 1 0.5])
hold off
set(gca,'XScale','log')
set(gca,'Xtick',vecPlotC,'xticklabel',vecContrasts)
ylabel('Orientation decoding accuracy')
xlabel('Stimulus contrast (%)')
title(['t-test vs 0.25, MH, p=' sprintf('%.3f,',p) sprintf('\n0.5-100, H-M: p=%.3f, H-chance, p=%.3f, M-chance, p=%.3f',pHM,pHC,pMC)])

if boolSavePlots
	drawnow;
	jFig = get(handle(gcf), 'JavaFrame');
	jFig.setMaximized(true);
	figure(gcf);
	drawnow;
	intFigCounter = intFigCounter + 1;
	strFig = sprintf('Meta%s%d_HMM_ori_decoding_raw',strNPS,intFigCounter);
	export_fig([strFig '.tif']);
	export_fig([strFig '.pdf']);
end

%% occurrence prob in time
vecBlur = normpdf(-10:10,0,5);
vecBlur = vecBlur./sum(vecBlur);
intClusters = size(matOccurrenceSelf,1);
vecOccurrenceT = mean(matOccurrenceT,1);
vecMeanOccSelf = mean(matOccurrenceSelf,1);
vecSEMOccSelf = std(matOccurrenceSelf,[],1)/sqrt(intClusters);
[h,vecSelfP]=ttest(matOccurrenceSelf,1);
[a,b,vecCorrSelfP]=fdr_bh(vecSelfP);
vecFiltCorrSelfP = imfiltreflectpad(vecCorrSelfP,vecBlur);
indSignificantSelf = vecFiltCorrSelfP < 0.05;

vecMeanOccOther = mean(matOccurrenceOther,1);
vecMeanOccOther(vecMeanOccOther<0.5) = 0.5;
vecSEMOccOther = std(matOccurrenceOther,[],1)/sqrt(intClusters);
[h,vecOtherP]=ttest(matOccurrenceOther,1);
[a,b,vecCorrOtherP]=fdr_bh(vecOtherP);
vecFiltCorrOtherP = imfiltreflectpad(vecCorrOtherP,vecBlur);
indSignificantOther = vecFiltCorrOtherP < 0.05;

vecBlur2 = normpdf(-2:2,0,1);
vecBlur2 = vecBlur2./sum(vecBlur2);
matOccurrenceStimNorm = bsxfun(@rdivide,matOccurrenceStim,vecOccurrenceBase);
vecMeanOccStim = mean(matOccurrenceStimNorm,1);
vecSEMOccStim = std(matOccurrenceStimNorm,[],1)/sqrt(intClusters);
[h,vecStimP]=ttest(matOccurrenceStimNorm,1);
[a,b,vecCorrStimP]=fdr_bh(vecStimP);
vecFiltCorrStimP = imfiltreflectpad(vecCorrStimP,vecBlur2);
indSignificantStim = vecFiltCorrStimP < 0.05;

matOccurrenceOffsetNorm = bsxfun(@rdivide,matOccurrenceOffset,vecOccurrenceOffsetBase);
vecMeanOccOffset = mean(matOccurrenceOffsetNorm,1);
vecSEMOccOffset = std(matOccurrenceOffsetNorm,[],1)/sqrt(intClusters);
[h,vecStimP]=ttest(matOccurrenceOffsetNorm,1);
[a,b,vecCorrStimP]=fdr_bh(vecStimP);
vecFiltCorrStimP = imfiltreflectpad(vecCorrStimP,vecBlur2);
indSignificantOffset = vecFiltCorrStimP < 0.05;

%plot
figure
subplot(2,2,1)
stairs(vecOccurrenceT,9+0.5*indSignificantSelf,'r');
hold on;
stairs(vecOccurrenceT,8.9+0.5*indSignificantOther,'b');
plot([min(vecOccurrenceT) max(vecOccurrenceT)],[1 1],'k--');
errorfill(vecOccurrenceT,vecMeanOccSelf,vecSEMOccSelf,[1 0 0],[1 0.7 0.7]);
errorfill(vecOccurrenceT,vecMeanOccOther,vecSEMOccOther,[0 0 1],[0.7 0.7 1]);
hold off;
title(sprintf('Blue=others, red=self'));
set(gca,'yscale','log');
xlim([min(vecOccurrenceT) max(vecOccurrenceT)]);
ylim([0.2 10]);
xlabel('Time (s)');
ylabel('Occ. prob. relative to random');

subplot(2,2,2)
dblThresh = mean(vecMeanOccStim)+2*std(vecMeanOccStim);
stairs(vecOccurrenceT,2.7+0.1*(vecMeanOccStim>dblThresh),'r');
hold on
errorfill([-5 10],[1 1]*mean(vecMeanOccStim),[2 2]*std(vecMeanOccStim),[0 0 0],[0.7 0.7 0.7]);
errorfill(vecOccurrenceT,vecMeanOccStim,vecSEMOccStim,[1 0 0],[1 0.7 0.7]);
hold off;
title(sprintf('centered on stim onset'));
%set(gca,'yscale','log');
xlim([-5 10])
ylim([0 3])
xlabel('Time (s)')
ylabel('Occ. prob. relative to random');


subplot(2,2,3)
vecLimX = [-8 8];
dblThresh = mean(vecMeanOccOffset)+2*std(vecMeanOccOffset);
stairs(vecOccurrenceT,4.5+0.2*(vecMeanOccOffset>dblThresh),'r');
hold on
errorfill(vecLimX,[1 1]*mean(vecMeanOccOffset),[2 2]*std(vecMeanOccOffset),[0 0 0],[0.7 0.7 0.7]);
errorfill(vecOccurrenceT,vecMeanOccOffset,vecSEMOccOffset,[1 0 0],[1 0.7 0.7]);
hold off;
title(sprintf('centered on licking response'));
%set(gca,'yscale','log');
xlim(vecLimX)
ylim([0 5])
xlabel('Time (s)')
ylabel('Occ. prob. relative to random');



	
if boolSavePlots
	drawnow;
	jFig = get(handle(gcf), 'JavaFrame');
	jFig.setMaximized(true);
	figure(gcf);
	drawnow;
	intFigCounter = intFigCounter + 1;
	strFig = sprintf('Meta%s%d_temporal_dependence_occurrences_MinDur_raw',strNPS,intFigCounter);
	export_fig([strFig '.tif']);
	export_fig([strFig '.pdf']);
end

return

%% consistency
subplot(2,2,1)
errorfill(vecTimeOff,mean(matAccOff,1),std(matAccOff,[],1)/sqrt(8))
ylim([0.2 0.5])
ylabel('Mean orientation decoding accuracy');
xlabel('Time from stimulus offset (s)');
title('One-second sliding window')
%}

figure
subplot(2,2,2);
vecMissTSS = matTempSeqStab(:,1);
vecHitsTSS = matTempSeqStab(:,2);
vecBaseTSS = matTempSeqStab(:,3);
[h,dblP_HM]=ttest2(vecHitsTSS,vecMissTSS);
[h,dblP_HB]=ttest2(vecHitsTSS,vecBaseTSS);
[h,dblP_MB]=ttest2(vecMissTSS,vecBaseTSS);
errorbar(1,nanmean(vecMissTSS),nanstd(vecMissTSS)/sqrt(length(vecMissTSS)),'rx');
hold on
errorbar(2,nanmean(vecHitsTSS),nanstd(vecHitsTSS)/sqrt(length(vecHitsTSS)),'gx');
errorbar(3,nanmean(vecBaseTSS),nanstd(vecBaseTSS)/sqrt(length(vecBaseTSS)),'bx');
hold off
ylim([0.91 0.96])
ylabel('Mean variability of temporal sequence (s)');
xlabel('Behavioral response speed');
set(gca,'xtick',1:3,'xticklabel',{'Miss','Hit','Baseline'});
title(sprintf('Paired t-test h-m=%.3f, h-b=%.3f, m-b, p=%.3f',dblP_HM,dblP_HB,dblP_MB));

%decode acc; [low-consist high-consist miss hit]
%matDecodingAnimals 

subplot(2,2,3);
vecMissDec = matDecodingAnimals(:,3);
vecHitsDec = matDecodingAnimals(:,4);
[h,dblP_HM]=ttest(vecHitsDec,vecMissDec);
errorbar(1,nanmean(vecMissDec),nanstd(vecMissDec)/sqrt(length(vecMissDec)),'rx');
hold on
errorbar(2,nanmean(vecHitsDec),nanstd(vecHitsDec)/sqrt(length(vecHitsDec)),'gx');
hold off
ylim([0 1+eps]);
ylabel('Mean decoding accuracy');
xlabel('Behavioral response speed');
set(gca,'xtick',1:3,'xticklabel',{'Miss','Hit'});
title(sprintf('Paired t-test h-m=%.3f',dblP_HM));

subplot(2,2,4);
vecLoCoDec = matDecodingAnimals(:,1);
vecHiCoDec = matDecodingAnimals(:,2);
[h,dblP_HM]=ttest(vecLoCoDec,vecHiCoDec);
errorbar(1,nanmean(vecLoCoDec),nanstd(vecLoCoDec)/sqrt(length(vecLoCoDec)),'rx');
hold on
errorbar(2,nanmean(vecHiCoDec),nanstd(vecHiCoDec)/sqrt(length(vecHiCoDec)),'gx');
plot([1.2 1.8],[vecLoCoDec vecHiCoDec],'Color',[0.5 0.5 0.5])
hold off
ylim([0 1+eps]);
ylabel('Mean decoding accuracy');
xlabel('Sequence consistency');
set(gca,'xtick',1:3,'xticklabel',{'Lowest 50%','Highest 50%'});
title(sprintf('Paired t-test h-m=%.3f',dblP_HM));


if boolSavePlots
	drawnow;
	jFig = get(handle(gcf), 'JavaFrame');
	jFig.setMaximized(true);
	figure(gcf);
	drawnow;
	intFigCounter = intFigCounter + 1;
	strFig = sprintf('Meta%s%d_temporal_sequence stability_raw',strNPS,intFigCounter);
	export_fig([strFig '.tif']);
	export_fig([strFig '.pdf']);
end
%}

%% clean up
cd(strOldDir);