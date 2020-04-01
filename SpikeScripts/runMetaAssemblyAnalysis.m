%close all
clear all
%% parameters
strFigDir = 'D:\Data\Results\spikeAnalysis\meta';
strDataDir = 'D:\Data\Results\spikeAnalysis';
cellInclude = {'20140207','20140314','20140425','20140507','20140530','20140604','20140711','20140715'};%exclude:20140430
boolSavePlots = true;
boolUseNPS = false;
boolOnlyTuned = false;
intFigCounter = 0;

%% pre-allocate
intAnimal = 0;
vecAnimalSource = [];
vecSubPopulationSource = [];
matDecAcc = [];
matCont_dFoF = [];
matCont_AE = [];
matCont_PE = [];

cellAssSigmoids = [];

cellRTDependency = [];

matMeanConsistencyOverall = [];
vecSlowAssemblyConsistencies = [];
vecFastAssemblyConsistencies = [];
matSlowFastOverall = [];
matPowerSpectraSelf = [];
matPowerSpectraOthers = [];

matTempSeqStab = [];

matAssProps = [];
vecAssPropAnimal = [];
matAssPropShuff = [];
matOSIShuf = [];

matAssDist = [];
vecAssDistAnimal = [];
matAssDistShuff = [];

vecAssemblyOSIs = [];
vecNeuronOSIs = [];

vecDecAccLowConsist = [];
vecDecAccHighConsist = [];
vecDecAccMiss = [];
vecDecAccHit = [];
		
		
matPE_AccMean = [];
matBehavConsistency = [];

matOccurrenceT = [];
matOccurrenceSelf = [];
matOccurrenceOther = [];
matOccurrenceStim = [];
vecOccurrenceBase = [];
matOccurrenceOffset = [];
vecOccurrenceOffsetBase = [];
		
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
	strRec = getFlankedBy(strFile,['dataAssAnalAA' strNPS],[strOT '_']);
	intRec = find(strcmp(strRec,cellInclude),1);
	if ~isempty(intRec)
		%load file
		load(strFile);
		
		%update population source
		intPopDim = 5;
		intNrPops = size(cellSaveAssAnal,intPopDim);
		intAnimal = intAnimal + 1;
		vecAnimalSource(end+1) = intAnimal;
		vecSubPopulationSource(end+1) = intAnimal;
		if intNrPops > 1 %group populations
			%take mean over matrices
			intEndPos = intPopDim-1;
			vecOrigSize = size(cellSaveAssAnal);
			vecOrigSize(end) = [];
			vecMatDims = ones(size(vecOrigSize));
			boolDontStop = true;
			while boolDontStop
				%concatenate
				strExec = ['cellSaveAssAnal{' sprintf('%d,',vecMatDims) '1} = nanmean(cat(3,cellSaveAssAnal{' sprintf('%d,',vecMatDims) '1},cellSaveAssAnal{' sprintf('%d,',vecMatDims) '2}),3);'];
				eval(strExec);
				
				%loop to go through all possible locations
				intPos = intEndPos;
				vecMatDims(intPos) = vecMatDims(intPos) + 1;
				while any(vecMatDims>vecOrigSize)
					intPos = intPos - 1;
					if intPos == 0
						boolDontStop = false;
						break
					end
					vecMatDims(intPos) = vecMatDims(intPos) + 1;
					for intTempPos=intPos:intEndPos-1
						vecMatDims(intTempPos+1) = 1;
					end
				end
			end
			%remove last population dimension
			strExec = ['cellSaveAssAnal(' repmat(':,',[1 intEndPos]) '2) = [];'];
			eval(strExec);
		end
		
		%assign data ass sigmoids
		cellAssSigmoids = cat(intPopDim,cellAssSigmoids,cellSaveAssAnal);
		
		%assign data cont anal
		matDecAcc = cat(1,matDecAcc,cellSaveContAnal{1});
		matCont_dFoF = cat(3,matCont_dFoF,cellSaveContAnal{2});
		matCont_AE = cat(3,matCont_AE,cellSaveContAnal{3});
		matCont_PE = cat(3,matCont_PE,cellSaveContAnal{4});
		
		cellRTDependency = cat(2,cellRTDependency,cellSaveRTDependency);
		intCounterF1 = intCounterF1 + 1;
	end
	
	%check if data part 2
	strFile = sDir(intFile).name;
	strRec = getFlankedBy(strFile,['dataAssAnalTOAA' strNPS],[strOT '_']);
	intRec = find(strcmp(strRec,cellInclude),1);
	if ~isempty(intRec)
		%load file
		load(strFile);
		
		%update population source
		intNrPops = size(cellSaveTemporalOrder,1);
		if intNrPops > 1 %group populations
			%take mean over matrices
			cellSaveTemporalOrder{1,1} = nanmean(cat(3,cellSaveTemporalOrder{1,1},cellSaveTemporalOrder{2,1}),3);
			for intType=2:3
				cellSaveTemporalOrder{1,intType} = [cellSaveTemporalOrder{1,intType} cellSaveTemporalOrder{2,intType}];
			end
			cellSaveTemporalOrder{1,4} = nanmean(cat(3,cellSaveTemporalOrder{1,4},cellSaveTemporalOrder{2,4}),3);
			cellSaveTemporalOrder(2,:) = [];
			
			%cellSaveTemporalSequenceStability{1} = mean(cat(3,cellSaveTemporalSequenceStability{1},cellSaveTemporalSequenceStability{2}),3);
			
			for intType=1:size(cellSaveAssemblyProperties,2)
				cellSaveAssemblyProperties{1,intType} = [cellSaveAssemblyProperties{1,intType} cellSaveAssemblyProperties{2,intType}];
			end
			cellSaveAssemblyProperties(2,:) = [];
			
		end
		matAssPropShufTemp = cellSaveAssemblyProperties{2};
		matOSIShufTemp = cellSaveAssemblyProperties{10};
		vecProps = true(1,length(cellSaveAssemblyProperties));
		vecProps([2 10]) = false;
		matAssPropsTemp = cell2mat(cellSaveAssemblyProperties(vecProps)');
		
		
		%power spectra
		intTotNum = 0;
		for intPopulation = 1:size(cellSaveAssemblyRecurrence,1)
			for intAssembly = 1:size(cellSaveAssemblyRecurrence,2)
				vecSpectrumSelf = cellSaveAssemblyRecurrence{intPopulation,intAssembly,2}';
				if isempty(vecSpectrumSelf),continue;end
				if length(vecSpectrumSelf) < 255,vecSpectrumSelf = [vecSpectrumSelf nan(1,255-length(vecSpectrumSelf))];end
				matPowerSpectraSelf = [matPowerSpectraSelf; vecSpectrumSelf];
				intTotNum = intTotNum + 1;
				
				vecSpectrumOthers = cellSaveAssemblyRecurrence{intPopulation,intAssembly,3}';
				if isempty(vecSpectrumOthers),continue;end
				if length(vecSpectrumOthers) < 255,vecSpectrumOthers = [vecSpectrumOthers nan(1,255-length(vecSpectrumOthers))];end
				matPowerSpectraOthers = [matPowerSpectraOthers; vecSpectrumOthers];
			end
		end
		vecFrequencies = cellSaveAssemblyRecurrence{1,1,1}';
		%fprintf('Session %s; properties %d assemblies, recurrence %d assemblies\n',strRec,size(matAssPropsTemp,2),intTotNum);
		
		%assign data
		matMeanConsistencyOverall = [matMeanConsistencyOverall; cellSaveTemporalOrder{1}]; %correlation temporal order within assemblies with temporal order overall
		vecSlowAssemblyConsistencies = [vecSlowAssemblyConsistencies nanmean(cellSaveTemporalOrder{2})]; %correlation temporal order of assembly sequence during slow trials with overall assembly sequence
		vecFastAssemblyConsistencies = [vecFastAssemblyConsistencies nanmean(cellSaveTemporalOrder{3})]; %correlation temporal order of assembly sequence during fast trials with overall assembly sequence
		matSlowFastOverall = [matSlowFastOverall; cellSaveTemporalOrder{4}]; %correlation temporal order of assembly sequence during fast trials with overall assembly sequence
		%matTempSeqStab = [matTempSeqStab; cellSaveTemporalSequenceStability{1}];
		matAssProps = [matAssProps matAssPropsTemp];
		vecAssPropAnimal = [vecAssPropAnimal ones(1,size(matAssPropsTemp,2))*(intCounterF2+1)];
		matAssPropShuff = [matAssPropShuff matAssPropShufTemp];
		matOSIShuf = [matOSIShuf matOSIShufTemp];
		
		intCounterF2 = intCounterF2 + 1;
		
		
		vecNumClusts(end+1) = sum(logical(matAssPropsTemp(7,:)) & ~any(isnan(matAssPropsTemp) | isinf(matAssPropsTemp),1));
	end
	
	%check if data part 3
	strFile = sDir(intFile).name;
	strRec = getFlankedBy(strFile,['dataAssAnalThreeAA' strNPS],[strOT '_']);
	intRec = find(strcmp(strRec,cellInclude),1);
	if ~isempty(intRec)
		%load file
		load(strFile);
		
		%update population source
		intNrPops = size(cellSaveAssemblyDistances,1);
		if intNrPops > 1 %group populations
			for intType=1:size(cellSaveAssemblyDistances,2)
				cellSaveAssemblyDistances{1,intType} = [cellSaveAssemblyDistances{1,intType} cellSaveAssemblyDistances{2,intType}];
			end
			cellSaveAssemblyDistances(2,:) = [];
			
			cellSaveTemporalSequenceStability{1} = mean(cat(3,cellSaveTemporalSequenceStability{1},cellSaveTemporalSequenceStability{2}),3);
			
			for intType=1:size(cellSaveOSIs,2)
				cellSaveOSIs{1,intType} = [cellSaveOSIs{1,intType} cellSaveOSIs{2,intType}];
			end
			cellSaveOSIs(2,:) = [];
		end
		matAssDistShufTemp = cellSaveAssemblyDistances{2};
		vecProps = true(1,length(cellSaveAssemblyDistances));
		vecProps(2) = false;
		matAssDistTemp = cell2mat(cellSaveAssemblyDistances(vecProps)');
		
		%assign data
		matAssDist = [matAssDist matAssDistTemp];
		vecAssDistAnimal = [vecAssDistAnimal ones(1,size(matAssDistTemp,2))*(intCounterF3+1)];
		matAssDistShuff = [matAssDistShuff matAssDistShufTemp];
		%matTempSeqStab = [matTempSeqStab; cellSaveTemporalSequenceStability{1}];
		
		% save data
		vecAssemblyOSIs = [vecAssemblyOSIs cellSaveOSIs{1}];
		vecNeuronOSIs = [vecNeuronOSIs cellSaveOSIs{2}];
		
		intCounterF3 = intCounterF3 + 1;
	end
	
	%check if data part 4
	strFile = sDir(intFile).name;
	strRec = getFlankedBy(strFile,['dataAssAnalFourAA' strNPS],[strOT '_']);
	intRec = find(strcmp(strRec,cellInclude),1);
	if ~isempty(intRec)
		%load file
		load(strFile);
		intCounterF4 = intCounterF4 + 1;
		
		%{
		%assign data
		if numel(cellSaveTemporalSequenceStability) > 1,cellSaveTemporalSequenceStability{1} = mean(cat(3,cellSaveTemporalSequenceStability{1},cellSaveTemporalSequenceStability{2}),3);end
		matTempSeqStab = [matTempSeqStab; cellSaveTemporalSequenceStability{1}];
		
		%get consistency for hit/miss/slow/fast
		matBehavConsistency(intCounterF4,:) = nanmean(cat(3,cellSaveConsistencyPerRespType{1},cellSaveConsistencyPerRespType{end}),3);
		
		%get decoding accuracy for 50% high/low temporal consistency
		vecDecAccLowHigh = nanmean(cellfun(@mean,cellSaveConsistencyPerTrial),1);
		vecDecAccLowConsist(intCounterF4) = vecDecAccLowHigh(1);
		vecDecAccHighConsist(intCounterF4) = vecDecAccLowHigh(2);
		vecDecAccMiss(intCounterF4) = vecDecAccLowHigh(3);
		vecDecAccHit(intCounterF4) = vecDecAccLowHigh(4);
		
		%get decoding accuracy per trial as function of number of PEs
		matPE_AccMean(:,:,intCounterF4) = nanmean(cat(3,cellSaveDASPE{1},cellSaveDASPE{end}),3);
		%}
		
		%occurrences in time
		cellCat = cellSaveAssemblyRecurrence(:,:,4);
		matOccurrenceT = cat(1,matOccurrenceT,cell2mat(cellCat(:)));
		cellCat = cellSaveAssemblyRecurrence(:,:,5);
		matOccurrenceSelf = cat(1,matOccurrenceSelf,cell2mat(cellCat(:)));
		cellCat = cellSaveAssemblyRecurrence(:,:,6);
		matOccurrenceOther = cat(1,matOccurrenceOther,cell2mat(cellCat(:)));
		cellCat = cellSaveAssemblyRecurrence(:,:,7);
		matOccurrenceStim = cat(1,matOccurrenceStim,cell2mat(cellCat(:)));
		cellCat = cellSaveAssemblyRecurrence(:,:,8);
		vecOccurrenceBase = cat(1,vecOccurrenceBase,cell2mat(cellCat(:)));
		cellCat = cellSaveAssemblyRecurrence(:,:,9);
		matOccurrenceOffset = cat(1,matOccurrenceOffset,cell2mat(cellCat(:)));
		cellCat = cellSaveAssemblyRecurrence(:,:,10);
		vecOccurrenceOffsetBase = cat(1,vecOccurrenceOffsetBase,cell2mat(cellCat(:)));
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
intAnimals = size(matDecAcc,1);
vecContrasts = [0 0.5 2 8 32 100];

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
indSignificantSelf = vecFiltCorrSelfP < 0.01;

vecMeanOccOther = mean(matOccurrenceOther,1);
vecMeanOccOther(vecMeanOccOther<0.5) = 0.5;
vecSEMOccOther = std(matOccurrenceOther,[],1)/sqrt(intClusters);
[h,vecOtherP]=ttest(matOccurrenceOther,1);
[a,b,vecCorrOtherP]=fdr_bh(vecOtherP);
vecFiltCorrOtherP = imfiltreflectpad(vecCorrOtherP,vecBlur);
indSignificantOther = vecFiltCorrOtherP < 0.01;

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
ylim([0.5 10]);
xlabel('Time (s)');
ylabel('Occ. prob. relative to random');

subplot(2,2,2)
dblThresh = mean(vecMeanOccStim)+2*std(vecMeanOccStim);
stairs(vecOccurrenceT,2.2+0.1*(vecMeanOccStim>dblThresh),'r');
hold on
errorfill([-5 10],[1 1]*mean(vecMeanOccStim),[2 2]*std(vecMeanOccStim),[0 0 0],[0.7 0.7 0.7]);
errorfill(vecOccurrenceT,vecMeanOccStim,vecSEMOccStim,[1 0 0],[1 0.7 0.7]);
hold off;
title(sprintf('centered on stim onset'));
%set(gca,'yscale','log');
xlim([-5 10])
ylim([0.5 2.5])
xlabel('Time (s)')
ylabel('Occ. prob. relative to random');


subplot(2,2,3)
vecLimX = [-8 8];
dblThresh = mean(vecMeanOccOffset)+2*std(vecMeanOccOffset);
stairs(vecOccurrenceT,3.5+0.2*(vecMeanOccOffset>dblThresh),'r');
hold on
errorfill(vecLimX,[1 1]*mean(vecMeanOccOffset),[2 2]*std(vecMeanOccOffset),[0 0 0],[0.7 0.7 0.7]);
errorfill(vecOccurrenceT,vecMeanOccOffset,vecSEMOccOffset,[1 0 0],[1 0.7 0.7]);
hold off;
title(sprintf('centered on licking response'));
%set(gca,'yscale','log');
xlim(vecLimX)
ylim([0 4])
xlabel('Time (s)')
ylabel('Occ. prob. relative to random');



	
if boolSavePlots
	drawnow;
	jFig = get(handle(gcf), 'JavaFrame');
	jFig.setMaximized(true);
	figure(gcf);
	drawnow;
	intFigCounter = intFigCounter + 1;
	strFig = sprintf('Meta%s%d_temporal_dependence_occurrences_raw',strNPS,intFigCounter);
	export_fig([strFig '.tif']);
	export_fig([strFig '.pdf']);
end

%% temporal order
figure
subplot(2,2,1)
boxplot(matMeanConsistencyOverall)
ylabel('Pearson correlation of temporal order of neuronal activation (muCC)');
%ylim([-1-eps 1+eps]*max(get(gca,'ylim')));
ylim([-1-eps 1+eps]);
[h,pAvsO]=ttest2(matMeanConsistencyOverall(:,1),matMeanConsistencyOverall(:,2));
set(gca,'xtick',[1 2],'xticklabel',{'Real assemblies','Other pop. events'});
[h,pAvs0]=ttest(matMeanConsistencyOverall(:,1));
[h,pOvs0]=ttest(matMeanConsistencyOverall(:,2));
title(sprintf('t-tests: assemblies vs 0;p=%.3f, other vs. 0;p=%.3f, A vs. O; p=%.3f',pAvs0,pOvs0,pAvsO));

subplot(2,2,2);
vecMissTSS = matTempSeqStab(:,1);
vecSlowTSS = matTempSeqStab(:,3);
vecFastTSS = matTempSeqStab(:,4);
vecBaseTSS = matTempSeqStab(:,5);
vecHitTSS = mean(cat(2,vecSlowTSS,vecFastTSS),2);
[h,dblP_HM]=ttest2(vecHitTSS,vecMissTSS);
[h,dblP_HB]=ttest2(vecHitTSS,vecBaseTSS);
[h,dblP_MB]=ttest2(vecMissTSS,vecBaseTSS);
errorbar(1,nanmean(vecMissTSS),nanstd(vecMissTSS)/sqrt(length(vecMissTSS)),'rx');
hold on
errorbar(2,nanmean(vecHitTSS),nanstd(vecHitTSS)/sqrt(length(vecHitTSS)),'gx');
errorbar(3,nanmean(vecBaseTSS),nanstd(vecBaseTSS)/sqrt(length(vecBaseTSS)),'bx');
hold off
%ylim([0 1+eps]*max(get(gca,'ylim')));
ylabel('Mean variability of temporal sequence (s)');
xlabel('Behavioral response speed');
set(gca,'xtick',1:3,'xticklabel',{'Miss','Hit','Baseline'});
title(sprintf('Paired t-test h-m=%.3f, h-b=%.3f, m-b, p=%.3f',dblP_HM,dblP_HB,dblP_MB));

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

%% decoding accuracy analyses
%get decoding accuracy for 50% high/low temporal consistency
%vecDecAccLowHigh = nanmean(cellfun(@mean,cellSaveConsistencyPerTrial),1);
%vecDecAccLowConsist(intCounterF4) = vecDecAccLowHigh(1);
%vecDecAccHighConsist(intCounterF4) = vecDecAccLowHigh(2);

%get decoding accuracy per trial as function of number of PEs
%matPE_AccMean(:,:,intCounterF4) = nanmean(cat(3,cellSaveDASPE{1},cellSaveDASPE{end}),3);
%[vecRecurAccMean;vecNonRecurAccMean;vecAllAccMean]

%prep
matMeanDecAccPE = nanmean(matPE_AccMean,3);
matSEMDecAccPE = nanstd(matPE_AccMean,[],3)/sqrt(intAnimals);

intMaxPEs = size(matMeanDecAccPE,2)-2;
vecX = 0:(intMaxPEs+1);
vecLimX = [-0.5 intMaxPEs+1.5];
cellTicks = num2cell(num2str(vecX'));
cellTicks{end} = ['>' cellTicks{end-1}];
cellTicks = cellTicks';


figure
subplot(2,2,1)
plot(vecLimX,[1 1]/4,'k--');
hold on
errorbar(vecX,matMeanDecAccPE(1,:),matSEMDecAccPE(1,:));
hold off
ylim([0 1]);
xlim(vecLimX);
set(gca,'xtick',0:7,'xticklabel',cellTicks);
xlabel('Number of PEs during stim. pres.');
ylabel('Decoding accuracy');
title('Recurring PEs');

subplot(2,2,2)
plot(vecLimX,[1 1]/4,'k--');
hold on
errorbar(vecX,matMeanDecAccPE(2,:),matSEMDecAccPE(2,:));
hold off
ylim([0 1]);
xlim(vecLimX);
set(gca,'xtick',vecX,'xticklabel',cellTicks);
xlabel('Number of PEs during stim. pres.');
ylabel('Decoding accuracy');
title('Non-recurring PEs');

subplot(2,2,3)
plot(vecLimX,[1 1]/4,'k--');
hold on
errorbar(vecX,matMeanDecAccPE(3,:),matSEMDecAccPE(3,:));
hold off
ylim([0 1]);
xlim(vecLimX);
set(gca,'xtick',0:7,'xticklabel',cellTicks);
xlabel('Number of PEs during stim. pres.');
ylabel('Decoding accuracy');
title('All PEs');

if boolSavePlots
	drawnow;
	jFig = get(handle(gcf), 'JavaFrame');
	jFig.setMaximized(true);
	figure(gcf);
	drawnow;
	intFigCounter = intFigCounter + 1;
	strFig = sprintf('Meta%s%d_Decoding_PEs',strNPS,intFigCounter);
	export_fig([strFig '.tif']);
	export_fig([strFig '.pdf']);
end



figure
subplot(2,2,1)
plot([0 1],[1 1]/4,'k--');
hold on
plot([zeros(1,intAnimals)+0.3; ones(1,intAnimals)-0.3],[vecDecAccMiss; vecDecAccHit],'Color',[0.5 0.5 0.5]);
errorbar([0.2 0.8],[mean(vecDecAccMiss) mean(vecDecAccHit)],[std(vecDecAccMiss) std(vecDecAccHit)]/sqrt(intAnimals),'bx');
hold off
ylim([0 1]);
ylabel('Orientation decoding accuracy');
set(gca,'xtick',[0.2 0.8],'xticklabel',{'Misses','Hits'});
xlabel('Behavioral response');
[h,p]=ttest(vecDecAccMiss,vecDecAccHit);
title(sprintf('Decoding during hits vs misses, p=%.3f',p));
	
		
subplot(2,2,2)
vecPlotTypes = [1 2 5];
plot([zeros(1,intAnimals)+0.3; ones(1,intAnimals)-0.3],[matBehavConsistency(:,1) matBehavConsistency(:,2)]','Color',[0.5 0.5 0.5]);
hold on
errorbar([0.2 0.8 0.5],mean(matBehavConsistency(:,vecPlotTypes),1),std(matBehavConsistency(:,vecPlotTypes),[],1)/sqrt(intAnimals),'xb');
hold off
set(gca,'xtick',[0.2 0.5 0.8],'xticklabel',{'Hit','Pre','Miss'});
xlabel('Behavioral response')
ylabel('Mean temporal consistency (r)');
[h,pHM]=ttest2(matBehavConsistency(:,1),matBehavConsistency(:,2)); %hit/miss
[h,pHP]=ttest2(matBehavConsistency(:,1),matBehavConsistency(:,5)); %hit/pre
[h,pMP]=ttest2(matBehavConsistency(:,2),matBehavConsistency(:,5)); %miss/pre
title(sprintf('t-tests,hit-miss,p=%.3f;hit-pre,p=%.3f;miss-pre,p=%.3f',[pHM pHP pMP]));

subplot(2,2,3)
plot([0 1],[1 1]/4,'k--');
hold on
plot([zeros(1,intAnimals)+0.3; ones(1,intAnimals)-0.3],[vecDecAccLowConsist; vecDecAccHighConsist],'Color',[0.5 0.5 0.5]);
errorbar([0.2 0.8],[mean(vecDecAccLowConsist) mean(vecDecAccHighConsist)],[std(vecDecAccLowConsist) std(vecDecAccHighConsist)]/sqrt(intAnimals),'bx');
hold off
ylim([0 1]);
ylabel('Orientation decoding accuracy');
set(gca,'xtick',[0.2 0.8],'xticklabel',{'Lowest 50%','Highest 50%'});
xlabel('Mean temporal sequence consistency during stim. pres.');
[h,p]=ttest(vecDecAccLowConsist,vecDecAccHighConsist);
title(sprintf('Decoding during high vs low temporal consistency, p=%.3f',p));


if boolSavePlots
	drawnow;
	jFig = get(handle(gcf), 'JavaFrame');
	jFig.setMaximized(true);
	figure(gcf);
	drawnow;
	intFigCounter = intFigCounter + 1;
	strFig = sprintf('Meta%s%d_Decoding_behav_and_consistency',strNPS,intFigCounter);
	export_fig([strFig '.tif']);
	export_fig([strFig '.pdf']);
end

%% get assembly properties
%dist
vecAssemblyMeanDist = matAssDist(1,:);

%props
vecAssemblyMean_dPA = matAssProps(1,:);
vecLaggingRatio = matAssProps(2,:);
vecHitMissRatioPresence = matAssProps(3,:);
return
vecSlowFastRatioPresence = matAssProps(4,:);
vecAssOSI = matAssProps(5,:);
vecStimDrivenRatio = matAssProps(6,:);
vecAssPropShuffSD = std(matAssPropShuff,[],1);
vecAssPropShuffMean = mean(matAssPropShuff,1);
vecdPA_zScored = (vecAssemblyMean_dPA-vecAssPropShuffMean)./vecAssPropShuffSD;

vecAssemblyMean_OSI = matAssProps(8,:);
vecOSIShuffSD = std(matOSIShuf,[],1);
vecOSIShuffMean = mean(matOSIShuf,1);
vecdOSI_zScored = (vecAssemblyMean_OSI-vecOSIShuffMean)./vecOSIShuffSD;

indRealAssemblies = logical(matAssProps(7,:)) & ~any(isnan(matAssProps) | isinf(matAssProps),1);


%% hit/miss occurrence ratio permutation test
vecCorrHM = log(vecHitMissRatioPresence(indRealAssemblies));
histx(vecCorrHM)

dblZ = mean(vecCorrHM)/std(vecCorrHM);
dblP = 1 - (normcdf(dblZ,0,1) - normcdf(-dblZ,0,1));

%% OSIs & recurrence
%OSIs
vecNeuronOSIs = vecNeuronOSIs(~isnan(vecNeuronOSIs));
[h,p]=ttest2(vecNeuronOSIs,vecAssemblyOSIs);
vecBinsOSI = 0.05:0.1:0.95;
vecN_N=hist(vecNeuronOSIs,vecBinsOSI);
vecN_A=hist(vecAssemblyOSIs,vecBinsOSI);

figure
subplot(2,2,1)
plot(vecBinsOSI,vecN_N/sum(vecN_N),'gx-');
hold on;
plot(vecBinsOSI,vecN_A/sum(vecN_A),'rx-');
plot([1 1]*mean(vecNeuronOSIs),get(gca,'ylim'),'g--')
plot([1 1]*mean(vecAssemblyOSIs),get(gca,'ylim'),'r--')
hold off;
title(sprintf('OSIs neurons (green; %.3f) and assemblies (red; %.3f), p=%.3f',mean(vecNeuronOSIs),mean(vecAssemblyOSIs),p));
ylim([0 0.25]);
xlabel('OSI');
ylabel('Fraction of total')

%recurrence
vecMeanSelf = nanmean(matPowerSpectraSelf,1);
vecMeanOthers = nanmean(matPowerSpectraOthers,1);
vecNormalize = vecFrequencies > 8 & vecFrequencies < 12;
dblNormSelf = nanmean(vecMeanSelf(vecNormalize));
dblNormOthers = nanmean(vecMeanOthers(vecNormalize));
matPowerSpectraSelfNorm = matPowerSpectraSelf/dblNormSelf;
matPowerSpectraOthersNorm = matPowerSpectraOthers/dblNormOthers;
[h,p]=ttest(matPowerSpectraSelfNorm,matPowerSpectraOthersNorm);
[h,x,pc]=fdr_bh(p);
dblMax = 10^6;
dblStep = dblMax/2;

subplot(2,2,2);
stairs(vecFrequencies,(dblMax-dblStep)+(pc<0.05)*(dblStep/2));
hold on
errorfill(vecFrequencies,nanmean(matPowerSpectraSelfNorm,1),nanstd(matPowerSpectraSelfNorm,[],1)/sqrt(size(matPowerSpectraSelfNorm,1)),[1 0 0],[1 0.5 0.5]);
errorfill(vecFrequencies,nanmean(matPowerSpectraOthersNorm,1),nanstd(matPowerSpectraOthersNorm,[],1)/sqrt(size(matPowerSpectraOthersNorm,1)),[0 0 1],[0.5 0.5 1]);
hold off
set(gca,'yscale','log')
ylabel('Power (dB/Hz)');
xlabel('Frequency (Hz)');
title('Assembly recurrence Welch power spectrum; red=self, blue=others');

if boolSavePlots
	drawnow;
	jFig = get(handle(gcf), 'JavaFrame');
	jFig.setMaximized(true);
	figure(gcf);
	drawnow;
	intFigCounter = intFigCounter + 1;
	strFig = sprintf('Meta%s%d_OSI_and_recurrence',strNPS,intFigCounter);
	export_fig([strFig '.tif']);
	export_fig([strFig '.pdf']);
end

%% distances
%indRealAssemblies = logical(matAssProps(7,:)) & ~(any(isnan(matAssDist),1) | any(isinf(matAssProps),1));

vecAssDistShuffSD = std(matAssDistShuff,[],1);
vecAssDistShuffMean = mean(matAssDistShuff,1);

vecDist_zScored = (vecAssemblyMeanDist-vecAssDistShuffMean)./vecAssDistShuffSD;
vecDistReal = vecAssemblyMeanDist(indRealAssemblies) - vecAssDistShuffMean(indRealAssemblies);
figure
subplot(2,2,1);
vecDist_zScoredReal = vecDist_zScored(indRealAssemblies);
vecDistPlot = vecDist_zScoredReal;
vecHistBins = -3.5:1:3.5;
%vecHistBins = -3.75:0.5:3.75;
hist(vecDistPlot,vecHistBins)
xlim([-4 4])
[h,p]=ttest(vecDistPlot);
title(sprintf('Spatial distance assembly members, t-test vs 0, p=%.3f; mean=%.3f sd; frac >2sd: %.3f',p,nanmean(vecDistPlot),sum(abs(vecDistPlot)>2)/numel(vecDistPlot)));
xlabel('Mean distance normalized to shuffled distro (sd)');
ylabel('Nr of assemblies (count)');

%{
subplot(2,2,2);
vecDist_zScoredNR = vecDist_zScored(~indRealAssemblies);
vecDistPlot = vecDist_zScoredNR(~isnan(vecDist_zScoredNR));
vecHistBins = -3.5:1:3.5;
%vecHistBins = -3.75:0.5:3.75;
hist(vecDistPlot,vecHistBins)
xlim([-4 4])
[h,p]=ttest(vecDistPlot);
title(sprintf('Spatial distance assembly members, t-test vs 0, p=%.3f; mean=%.3f sd; frac >2sd: %.3f',p,mean(vecDistPlot),sum(abs(vecDistPlot)>2)/numel(vecDistPlot)));
xlabel('Mean distance normalized to shuffled distro (sd)');
ylabel('Nr of assemblies (count)');
%}

subplot(2,2,2);
scatter(vecDistPlot,vecSlowFastRatioPresence(indRealAssemblies));
[dblCorr,dblP] = corr(vecDistPlot',vecSlowFastRatioPresence(indRealAssemblies)','rows','pairwise');
xlabel('Mean distance normalized to shuffled distro (sd)');
ylabel('Fast/Slow presence ratio');
title(sprintf('Correlation; r=%.3f; p=%.3f',dblCorr,dblP));
xlim([-4 4])

subplot(2,2,3);
scatter(vecDistPlot,vecAssOSI(indRealAssemblies));
[dblCorr,dblP] = corr(vecDistPlot',vecAssOSI(indRealAssemblies)','rows','pairwise');
xlabel('Mean distance normalized to shuffled distro (sd)');
ylabel('Assembly OSI');
title(sprintf('Correlation; r=%.3f; p=%.3f',dblCorr,dblP));
xlim([-4 4])

subplot(2,2,4);
scatter(vecDistPlot,vecLaggingRatio(indRealAssemblies));
[dblCorr,dblP] = corr(vecDistPlot',vecLaggingRatio(indRealAssemblies)','rows','pairwise');
xlabel('Mean distance normalized to shuffled distro (sd)');
ylabel('Lagging ratio');
title(sprintf('Correlation; r=%.3f; p=%.3f',dblCorr,dblP));
xlim([-4 4])

if boolSavePlots
	drawnow;
	jFig = get(handle(gcf), 'JavaFrame');
	jFig.setMaximized(true);
	figure(gcf);
	drawnow;
	intFigCounter = intFigCounter + 1;
	strFig = sprintf('Meta%s%d_distance_dependence_assemblies',strNPS,intFigCounter);
	export_fig([strFig '.tif']);
	export_fig([strFig '.pdf']);
end

%%  assembly properties
%{
matAssProps(1,intAssembly) = vecAssemblyMean_dPA;
matAssProps(2,intAssembly) =  vecLaggingRatio;
matAssProps(3,intAssembly) = vecHitMissRatioPresence;
matAssProps(4,intAssembly) = vecSlowFastRatioPresence;
matAssProps(5,intAssembly) =  vecAssOSI;
matAssProps(6,intAssembly) =  vecStimDrivenRatio;
matAssProps(7,intAssembly) =  double(indRealAssemblies);
%}

figure
subplot(2,2,1);
vecdPA_zScoredReal = vecdPA_zScored(indRealAssemblies);
hist(vecdPA_zScoredReal,vecHistBins)
xlim([-4 4])
[h,p]=ttest(vecdPA_zScoredReal);
title(sprintf('Pref ori diff assembly members, t-test vs 0, p=%.3f; mean=%.3f sd; frac >2sd: %.3f',p,mean(vecdPA_zScoredReal),sum(abs(vecdPA_zScoredReal)>2)/numel(vecdPA_zScoredReal)));
xlabel('Mean dPO normalized to shuffled distro (sd)');
ylabel('Nr of assemblies (count)');


subplot(2,2,2);
vec_dPA = rad2ang(vecAssemblyMean_dPA(indRealAssemblies));
scatter(vec_dPA,vecSlowFastRatioPresence(indRealAssemblies));
[dblCorr,dblP] = corr(vecAssemblyMean_dPA(indRealAssemblies)',vecSlowFastRatioPresence(indRealAssemblies)');
xlabel('Mean dPO (deg)');
ylabel('Fast/Slow presence ratio');
title(sprintf('Correlation; r=%.3f; p=%.3f',dblCorr,dblP));

subplot(2,2,3);
vecOSI_zScoredReal = vecdOSI_zScored(indRealAssemblies);
hist(vecOSI_zScoredReal,vecHistBins)
xlim([-4 4])
[h,p]=ttest(vecOSI_zScoredReal);
title(sprintf('Pref ori diff assembly members, t-test vs 0, p=%.3f; mean=%.3f sd; frac >2sd: %.3f',p,mean(vecOSI_zScoredReal),sum(abs(vecOSI_zScoredReal)>2)/numel(vecOSI_zScoredReal)));
xlabel('Mean OSI normalized to shuffled distro (sd)');
ylabel('Nr of assemblies (count)');



%{
subplot(2,2,3);
scatter(vec_dPA,vecStimDrivenRatio(indRealAssemblies));
[dblCorr,dblP] = corr(vecAssemblyMean_dPA(indRealAssemblies)',vecStimDrivenRatio(indRealAssemblies)');
xlabel('Mean dPO (deg)');
ylabel('Stim-driven ratio');
title(sprintf('Correlation; r=%.3f; p=%.3f',dblCorr,dblP));
%}
subplot(2,2,4);
scatter(vec_dPA,vecLaggingRatio(indRealAssemblies));
[dblCorr,dblP] = corr(vecAssemblyMean_dPA(indRealAssemblies)',vecLaggingRatio(indRealAssemblies)');
xlabel('Mean dPO (deg)');
ylabel('Lagging ratio');
title(sprintf('Correlation; r=%.3f; p=%.3f',dblCorr,dblP));

if boolSavePlots
	drawnow;
	jFig = get(handle(gcf), 'JavaFrame');
	jFig.setMaximized(true);
	figure(gcf);
	drawnow;
	intFigCounter = intFigCounter + 1;
	strFig = sprintf('Meta%s%d_orientation_tuning_assemblies',strNPS,intFigCounter);
	export_fig([strFig '.tif']);
	export_fig([strFig '.pdf']);
end

% other properties
figure
subplot(2,2,1);
scatter(vecStimDrivenRatio(indRealAssemblies),vecAssOSI(indRealAssemblies));
[dblCorr,dblP] = corr(vecStimDrivenRatio(indRealAssemblies)',vecAssOSI(indRealAssemblies)');
xlabel('Stim-driven ratio');
ylabel('OSI');
title(sprintf('Correlation; r=%.3f; p=%.3f',dblCorr,dblP));

subplot(2,2,2);
scatter(vecLaggingRatio(indRealAssemblies),vecSlowFastRatioPresence(indRealAssemblies));
[dblCorr,dblP] = corr(vecLaggingRatio(indRealAssemblies)',vecSlowFastRatioPresence(indRealAssemblies)');
xlabel('Lagging ratio');
ylabel('Fast/Slow presence ratio');
title(sprintf('Correlation; r=%.3f; p=%.3f',dblCorr,dblP));

subplot(2,2,3);
scatter(vecLaggingRatio(indRealAssemblies),vecStimDrivenRatio(indRealAssemblies));
[dblCorr,dblP] = corr(vecLaggingRatio(indRealAssemblies)',vecStimDrivenRatio(indRealAssemblies)');
xlabel('Lagging ratio');
ylabel('Stim-driven ratio');
title(sprintf('Correlation; r=%.3f; p=%.3f',dblCorr,dblP));

subplot(2,2,4);
scatter(vecSlowFastRatioPresence(indRealAssemblies),vecStimDrivenRatio(indRealAssemblies));
[dblCorr,dblP] = corr(vecStimDrivenRatio(indRealAssemblies)',vecSlowFastRatioPresence(indRealAssemblies)');
xlabel('Fast/Slow presence ratio');
ylabel('Stim-driven ratio');
title(sprintf('Correlation; r=%.3f; p=%.3f',dblCorr,dblP));

if boolSavePlots
	drawnow;
	jFig = get(handle(gcf), 'JavaFrame');
	jFig.setMaximized(true);
	figure(gcf);
	drawnow;
	intFigCounter = intFigCounter + 1;
	strFig = sprintf('Meta%s%d_assembly_property_correlations_raw',strNPS,intFigCounter);
	export_fig([strFig '.tif']);
	export_fig([strFig '.pdf']);
end

%% decoding with neurons & contrast response functions

%save data
%cellSaveAssAnal{1,intPopulation} = vecDecodingPerformance; %dFoF AE Ass
%cellSaveAssAnal{2,intPopulation} = [vecMeanHit;vecMeanMiss]; %hit-miss assembly presence (time) during trials
%cellSaveAssAnal{3,intPopulation} = [vecTrialOccurrences_Hit;vecTrialOccurrences_Miss]/intTrialsPerContrast; %hit-miss mean assemblies per trial
%cellSaveAssAnal{4,intPopulation} = [vecTrialMeanSize_Hit;vecTrialMeanSize_Miss]; %hit-miss mean neurons per assembly


vecY = mean(matDecAcc(:,1:2),1);
vecE = std(matDecAcc(:,1:2),[],1)/sqrt(intAnimals);
[h,dblP_dec_Act_AE]=ttest(matDecAcc(:,1),matDecAcc(:,2));

figure
%decoding
subplot(2,2,1)
plot([0 1],[0.25 0.25],'k--')
hold on
errorbar([0.3 0.7],vecY,vecE,'linestyle','none','marker','x')
hold off
xlim([0 1])
ylim([0 0.5])
ylabel('CV Decoding accuracy')
set(gca,'xtick',[0.3 0.7],'xticklabel',{'dF/F0','Activation Events'})
title(sprintf('Mean decoding accuracy over animals (n=%d) +/- st. err., p=%.3f',intAnimals,dblP_dec_Act_AE))

%contrast response curves
%dF/F
matMean_dFoF = nanmean(matCont_dFoF,3);
matSD_dFoF = nanstd(matCont_dFoF,[],3);

subplot(2,2,2)
vecPlotC = vecContrasts;
vecPlotC(1) = 0.2;
hold on
errorfill(vecPlotC,matMean_dFoF(2,:),matSD_dFoF(2,:)/sqrt(intAnimals),[1 0 0],[1 0.5 0.5]);
errorfill(vecPlotC,matMean_dFoF(1,:),matSD_dFoF(1,:)/sqrt(intAnimals),[0 1 0],[0.5 1 0.5]);
hold off
set(gca,'XScale','log')
set(gca,'Xtick',vecPlotC,'xticklabel',vecContrasts)
ylabel('Mean dF/F0')
xlabel('Stimulus contrast (%)')
[h,pC0]=ttest(matCont_dFoF(1,1,:),matCont_dFoF(2,1,:));
[h,pC05]=ttest(matCont_dFoF(1,2,:),matCont_dFoF(2,2,:));
[h,pC2]=ttest(matCont_dFoF(1,3,:),matCont_dFoF(2,3,:));
[h,pC8]=ttest(matCont_dFoF(1,4,:),matCont_dFoF(2,4,:));
[h,pC32]=ttest(matCont_dFoF(1,5,:),matCont_dFoF(2,5,:));
[h,pC100]=ttest(matCont_dFoF(1,6,:),matCont_dFoF(2,6,:));


%AE
matMean_AE = nanmean(matCont_AE,3);
matSD_AE = nanstd(matCont_AE,[],3);

subplot(2,2,3)
vecPlotC = vecContrasts;
vecPlotC(1) = 0.2;
hold on
errorfill(vecPlotC,matMean_AE(2,:),matSD_AE(2,:)/sqrt(intAnimals),[1 0 0],[1 0.5 0.5]);
errorfill(vecPlotC,matMean_AE(1,:),matSD_AE(1,:)/sqrt(intAnimals),[0 1 0],[0.5 1 0.5]);
hold off
set(gca,'XScale','log')
set(gca,'Xtick',vecPlotC,'xticklabel',vecContrasts)
ylabel('Mean AEs')
xlabel('Stimulus contrast (%)')
[h,pC0]=ttest(matCont_AE(1,1,:),matCont_AE(2,1,:));
[h,pC05]=ttest(matCont_AE(1,2,:),matCont_AE(2,2,:));
[h,pC2]=ttest(matCont_AE(1,3,:),matCont_AE(2,3,:));
[h,pC8]=ttest(matCont_AE(1,4,:),matCont_AE(2,4,:));
[h,pC32]=ttest(matCont_AE(1,5,:),matCont_AE(2,5,:));
[h,pC100]=ttest(matCont_AE(1,6,:),matCont_AE(2,6,:));

%PE
matMean_PE = nanmean(matCont_PE,3);
matSD_PE = nanstd(matCont_PE,[],3);

subplot(2,2,4)
vecPlotC = vecContrasts;
vecPlotC(1) = 0.2;
hold on
errorfill(vecPlotC,matMean_PE(2,:),matSD_PE(2,:)/sqrt(intAnimals),[1 0 0],[1 0.5 0.5]);
errorfill(vecPlotC,matMean_PE(1,:),matSD_PE(1,:)/sqrt(intAnimals),[0 1 0],[0.5 1 0.5]);
hold off
set(gca,'XScale','log')
set(gca,'Xtick',vecPlotC,'xticklabel',vecContrasts)
ylabel('Mean PEs')
xlabel('Stimulus contrast (%)')
[h,pC0]=ttest(matCont_PE(1,1,:),matCont_PE(2,1,:));
[h,pC05]=ttest(matCont_PE(1,2,:),matCont_PE(2,2,:));
[h,pC2]=ttest(matCont_PE(1,3,:),matCont_PE(2,3,:));
[h,pC8]=ttest(matCont_PE(1,4,:),matCont_PE(2,4,:));
[h,pC32]=ttest(matCont_PE(1,5,:),matCont_PE(2,5,:));
[h,pC100]=ttest(matCont_PE(1,6,:),matCont_PE(2,6,:));


if boolSavePlots
	drawnow;
	jFig = get(handle(gcf), 'JavaFrame');
	jFig.setMaximized(true);
	figure(gcf);
	drawnow;
	intFigCounter = intFigCounter + 1;
	strFig = sprintf('Meta%s%d_contrast_responses_raw',strNPS,intFigCounter);
	export_fig([strFig '.tif']);
	export_fig([strFig '.pdf']);
end

%% assembly stuff
%cellSaveAssAnal{1,1,1,1,intAnimal} = [vecMeanHit;vecMeanMiss]; %hit-miss assembly presence (time) during trials
%cellSaveAssAnal{1,2,1,1,intAnimal} = [vecMeanHitNR;vecMeanMissNR]; %hit-miss assembly presence (time) during trials
%cellSaveAssAnal{1,1,2,1,intAnimal} = [vecMeanFast;vecMeanSlow]; %hit-miss assembly presence (time) during trials
%cellSaveAssAnal{1,2,2,1,intAnimal} = [vecMeanFastNR;vecMeanSlowNR]; %hit-miss assembly presence (time) during trials
%cellSaveAssAnal{1,1,1,2,intAnimal} = [mean(vecMeanPreHitRealAssAct);mean(vecMeanPreMissRealAssAct)]; %hit-miss assembly presence (time) preceding trials
%cellSaveAssAnal{1,2,1,2,intAnimal} = [mean(vecMeanPreHitNRAssAct);mean(vecMeanPreMissNRAssAct)]; %hit-miss assembly presence (time) during trials
%cellSaveAssAnal{1,1,2,2,intAnimal} = [mean(vecMeanPreFastRealAssAct);mean(vecMeanPreSlowRealAssAct)]; %hit-miss assembly presence (time) during trials
%cellSaveAssAnal{1,2,2,2,intAnimal} = [mean(vecMeanPreFastNRAssAct);mean(vecMeanPreSlowNRAssAct)]; %hit-miss assembly presence (time) during trials

%% create plotting variables and perform plotting
hNonRecAssemblyCorrelates = figure;
hRealAssemblyCorrelates = figure;
hNonRecAssemblyCorrelatesSlowFast = figure;
hRealAssemblyCorrelatesSlowFast = figure;

%pre-allocate selection index vectors
cellRespTypeNames = {'Hit','Miss','Fast','Slow'};
cellRecurTypeNames = {'Real','NonRec'};
cellPrePostTypeNames = {'','Pre'};

%plot
for intRespPlotType=1%:2
	for intRecurType=2
		%shorten name
		strRecur = cellRecurTypeNames{intRecurType};
		
		%get target figure
		if intRespPlotType == 1 %hit/miss
			strGood = cellRespTypeNames{1}; %hit
			strBad = cellRespTypeNames{2}; %miss
			if intRecurType == 1 %recurring
				figure(hRealAssemblyCorrelates);
			elseif intRecurType == 2 %non-recurring
				figure(hNonRecAssemblyCorrelates);
			end
		elseif intRespPlotType == 2 %fast/slow
			strGood = cellRespTypeNames{3}; %fast
			strBad = cellRespTypeNames{4}; %slow
			if intRecurType == 1 %recurring
				figure(hRealAssemblyCorrelatesSlowFast);
			elseif intRecurType == 2 %non-recurring
				figure(hNonRecAssemblyCorrelatesSlowFast);
			end
		end
		matEffectSize = nan(2,6,intAnimals);
		%{
		%% plot presence
		matPresence = cell2mat(cellAssSigmoids(1,intRecurType,intRespPlotType,1,:));
		matMean = nanmean(matPresence,5);
		matSE = nanstd(matPresence,[],5)/sqrt(intAnimals);
		
		matPresPre = cell2mat(cellAssSigmoids(1,intRecurType,intRespPlotType,2,:));
		vecMeanPre = nanmean(matPresPre,5);
		vecSEPre = nanstd(matPresPre,[],5)/sqrt(intAnimals);
		
		subplot(2,2,1)
		errorbar([vecPlotC(1) vecPlotC(end)],vecMeanPre(1)*[1 1],vecSEPre(1)*[1 1],'color',[0 1 0],'linestyle','--');
		hold on
		errorbar([vecPlotC(1) vecPlotC(end)],vecMeanPre(2)*[1 1],vecSEPre(2)*[1 1],'color',[1 0 0],'linestyle','--');
		errorfill(vecPlotC,matMean(2,:),matSE(2,:),[1 0 0],[1 0.5 0.5]);
		errorfill(vecPlotC,matMean(1,:),matSE(1,:),[0 1 0],[0.5 1 0.5]);
		hold off
		set(gca,'XScale','log')
		set(gca,'Xtick',vecPlotC,'xticklabel',vecContrasts)
		ylabel('Mean assembly presence per trial (s)')
		xlabel('Stimulus contrast (%)')
		ylim([0 max(get(gca,'ylim'))]);
		title([strRecur ' ' strGood '/' strBad]);
		%}
		%%{
		%plot occurrence
		matOccurrence = cell2mat(cellAssSigmoids(2,intRecurType,intRespPlotType,1,:));
		matMean = nanmean(matOccurrence,5);
		matSE = nanstd(matOccurrence,[],5)/sqrt(intAnimals);
		
		matOccPre = cell2mat(cellAssSigmoids(2,intRecurType,intRespPlotType,2,:));
		vecMeanPre = nanmean(matOccPre,5);
		vecSEPre = nanstd(matOccPre,[],5)/sqrt(intAnimals);
		
		matOcc = squeeze(matOccurrence);
		matOPS = squeeze(matOccPre);
		matOccP = nan(2,6);
		vecRegP = nan(1,2);
		vecRegB = nan(1,2);
		for intR=1:2
			matThisOcc = bsxfun(@minus,squeeze(matOcc(intR,:,:)),matOPS(intR,:));
			%matThisOcc = squeeze(matOcc(intR,:,:));
			matC = repmat((1:6)',[1 intAnimals]);
			sStats=regstats(matThisOcc(:),matC(:),'linear',{'beta','rsquare','tstat'});
			vecRegB(intR) = sStats.beta(2);
			vecRegP(intR) = sStats.tstat.pval(2);
			
			for intC=1:6
				[h,matOccP(intR,intC)] = ttest(squeeze(matOPS(intR,:)),squeeze(matOcc(intR,intC,:))');
			end
			vec0 = squeeze(matOcc(intR,1,:))' - matOPS(intR,:);
			vec100 = squeeze(matOcc(intR,6,:))' - matOPS(intR,:);
			matEffectSize(intR,2,:) = vec100./(vec100+vec0);
		end
		[h,x,matOccP_corr]=fdr_bh(matOccP);
		%matOccP_corr
		
		subplot(2,2,1)
		errorbar([vecPlotC(1) vecPlotC(end)],vecMeanPre(1)*[1 1],vecSEPre(1)*[1 1],'color',[0 1 0],'linestyle','--');
		hold on
		errorbar([vecPlotC(1) vecPlotC(end)],vecMeanPre(2)*[1 1],vecSEPre(2)*[1 1],'color',[1 0 0],'linestyle','--');
		errorfill(vecPlotC,matMean(2,:),matSE(2,:),[1 0 0],[1 0.5 0.5]);
		errorfill(vecPlotC,matMean(1,:),matSE(1,:),[0 1 0],[0.5 1 0.5]);
		hold off
		set(gca,'XScale','log')
		set(gca,'Xtick',vecPlotC,'xticklabel',vecContrasts)
		ylim([0 max(get(gca,'ylim'))]);
		ylabel('Mean assembly occurrences per trial')
		xlabel('Stimulus contrast (%)')
		%title([strRecur ' ' strGood '/' strBad sprintf('; regP=%.3f,%.3f',vecRegP)]);
		title([strRecur ' ' strGood '/' strBad]);
		%%}
		
		
		% plot duration
		matDur = cell2mat(cellAssSigmoids(4,intRecurType,intRespPlotType,1,:));
		matMean = nanmean(matDur,5);
		matSE = nanstd(matDur,[],5)/sqrt(intAnimals);
		
		matDurPre = cell2mat(cellAssSigmoids(4,intRecurType,intRespPlotType,2,:));
		vecMeanPre = nanmean(matDurPre,5);
		vecSEPre = nanstd(matDurPre,[],5)/sqrt(intAnimals);
		
		matDDS = squeeze(matDur);
		matPDS = squeeze(matDurPre);
		matDP = nan(2,6);
		for intR=1:2
			matThisDur = bsxfun(@minus,squeeze(matDur(intR,:,:)),matPDS(intR,:));
			%matThisDur = squeeze(matDur(intR,:,:));
			matC = repmat((1:6)',[1 intAnimals]);
			sStats=regstats(matThisDur(:),matC(:),'linear',{'beta','rsquare','tstat'});
			vecRegB(intR) = sStats.beta(2);
			vecRegP(intR) = sStats.tstat.pval(2);
			
			for intC=1:6
				[h,matDP(intR,intC)] = ttest(squeeze(matPDS(intR,:)),squeeze(matDDS(intR,intC,:))');
			end
			vec0 = squeeze(matDDS(intR,1,:))' - matPDS(intR,:);
			vec100 = squeeze(matDDS(intR,6,:))' - matPDS(intR,:);
			matEffectSize(intR,4,:) = vec100./(vec100+vec0);
		end
		[h,x,matDP_corr]=fdr_bh(matDP);
		%matDP_corr
		
		subplot(2,2,2)
		errorbar([vecPlotC(1) vecPlotC(end)],vecMeanPre(1)*[1 1],vecSEPre(1)*[1 1],'color',[0 1 0],'linestyle','--');
		hold on
		errorbar([vecPlotC(1) vecPlotC(end)],vecMeanPre(2)*[1 1],vecSEPre(2)*[1 1],'color',[1 0 0],'linestyle','--');
		errorfill(vecPlotC,matMean(2,:),matSE(2,:),[1 0 0],[1 0.5 0.5]);
		errorfill(vecPlotC,matMean(1,:),matSE(1,:),[0 1 0],[0.5 1 0.5]);
		hold off
		set(gca,'XScale','log')
		set(gca,'Xtick',vecPlotC,'xticklabel',vecContrasts)
		ylabel('Mean duration per assembly occurrence (s)')
		xlabel('Stimulus contrast (%)')
		ylim([0 max(get(gca,'ylim'))]);
		%title([strRecur ' ' strGood '/' strBad sprintf('; regP=%.3f,%.3f',vecRegP)]);
		title([strRecur ' ' strGood '/' strBad]);
		
		
		%plot size
		matSize = cell2mat(cellAssSigmoids(3,intRecurType,intRespPlotType,1,:))./matDur;
		matMean = nanmean(matSize,5);
		matSE = nanstd(matSize,[],5)/sqrt(intAnimals);
		
		matSizePre = cell2mat(cellAssSigmoids(3,intRecurType,intRespPlotType,2,:))./matDurPre;
		vecMeanPre = nanmean(matSizePre,5);
		vecSEPre = nanstd(matSizePre,[],5)/sqrt(intAnimals);
		
		matDSS = squeeze(matSize);
		matPSS = squeeze(matSizePre);
		matSSP = nan(2,6);
		for intR=1:2
			matThisSize = bsxfun(@minus,squeeze(matSize(intR,:,:)),matPSS(intR,:));
			%matThisSize = squeeze(matSize(intR,:,:));
			matC = repmat((1:6)',[1 intAnimals]);
			sStats=regstats(matThisSize(:),matC(:),'linear',{'beta','rsquare','tstat'});
			vecRegB(intR) = sStats.beta(2);
			vecRegP(intR) = sStats.tstat.pval(2);
			
			for intC=1:6
				[h,matSSP(intR,intC)] = ttest(squeeze(matPSS(intR,:)),squeeze(matDSS(intR,intC,:))');
			end
			vec0 = squeeze(matDSS(intR,1,:))' - matPSS(intR,:);
			vec100 = squeeze(matDSS(intR,6,:))' - matPSS(intR,:);
			matEffectSize(intR,3,:) = vec100./(vec100+vec0);
		end
		[h,x,matSSP_corr]=fdr_bh(matSSP);
		%matSSP_corr
		
		subplot(2,2,3)
		errorbar([vecPlotC(1) vecPlotC(end)],vecMeanPre(1)*[1 1],vecSEPre(1)*[1 1],'color',[0 1 0],'linestyle','--');
		hold on
		errorbar([vecPlotC(1) vecPlotC(end)],vecMeanPre(2)*[1 1],vecSEPre(2)*[1 1],'color',[1 0 0],'linestyle','--');
		errorfill(vecPlotC,matMean(2,:),matSE(2,:),[1 0 0],[1 0.5 0.5]);
		errorfill(vecPlotC,matMean(1,:),matSE(1,:),[0 1 0],[0.5 1 0.5]);
		hold off
		set(gca,'XScale','log')
		set(gca,'Xtick',vecPlotC,'xticklabel',vecContrasts)
		ylabel('Mean # of neurons per sec assembly occurrence')
		xlabel('Stimulus contrast (%)')
		ylim([0 max(get(gca,'ylim'))]);
		%title([strRecur ' ' strGood '/' strBad sprintf('; regP=%.3f,%.3f',vecRegP)]);
		title([strRecur ' ' strGood '/' strBad]);
		
		
		
		%{
		%plot AEs
		matAEs = cell2mat(cellAssSigmoids(5,intRecurType,intRespPlotType,1,:));
		matMean = nanmean(matAEs,5);
		matSE = nanstd(matAEs,[],5)/sqrt(intAnimals);
		
		matAEsPre = cell2mat(cellAssSigmoids(5,intRecurType,intRespPlotType,2,:));
		vecMeanPre = nanmean(matAEsPre,5);
		vecSEPre = nanstd(matAEsPre,[],5)/sqrt(intAnimals);
		
		subplot(2,2,5)
		errorbar([vecPlotC(1) vecPlotC(end)],vecMeanPre(1)*[1 1],vecSEPre(1)*[1 1],'color',[0 1 0],'linestyle','--');
		hold on
		errorbar([vecPlotC(1) vecPlotC(end)],vecMeanPre(2)*[1 1],vecSEPre(2)*[1 1],'color',[1 0 0],'linestyle','--');
		errorfill(vecPlotC,matMean(2,:),matSE(2,:),[1 0 0],[1 0.5 0.5]);
		errorfill(vecPlotC,matMean(1,:),matSE(1,:),[0 1 0],[0.5 1 0.5]);
		hold off
		set(gca,'XScale','log')
		set(gca,'Xtick',vecPlotC,'xticklabel',vecContrasts)
		ylabel('Mean AEs per assembly occurrence')
		xlabel('Stimulus contrast (%)')
		ylim([0 max(get(gca,'ylim'))]);
		title([strRecur ' ' strGood '/' strBad]);
		%}
		%plot AEs / sec
		matAEsPerSec = cell2mat(cellAssSigmoids(5,intRecurType,intRespPlotType,1,:)) ./ matDur;
		matMean = nanmean(matAEsPerSec,5);
		matSE = nanstd(matAEsPerSec,[],5)/sqrt(intAnimals);
		
		matAEsPerSecPre = cell2mat(cellAssSigmoids(5,intRecurType,intRespPlotType,2,:)) ./ matDurPre;
		vecMeanPre = nanmean(matAEsPerSecPre,5);
		vecSEPre = nanstd(matAEsPerSecPre,[],5)/sqrt(intAnimals);
		
		matADS = squeeze(matAEsPerSec);
		matAPS = squeeze(matAEsPerSecPre);
		matAP = nan(2,6);
		for intR=1:2
			matThisA = bsxfun(@minus,squeeze(matAEsPerSec(intR,:,:)),matAPS(intR,:));
			%matThisA = squeeze(matAEsPerSec(intR,:,:));
			matC = repmat((1:6)',[1 intAnimals]);
			sStats=regstats(matThisA(:),matC(:),'linear',{'beta','rsquare','tstat'});
			vecRegB(intR) = sStats.beta(2);
			vecRegP(intR) = sStats.tstat.pval(2);
			
			for intC=1:6
				[h,matAP(intR,intC)] = ttest(squeeze(matAPS(intR,:)),squeeze(matADS(intR,intC,:))');
			end
			vec0 = squeeze(matADS(intR,1,:))' - matAPS(intR,:);
			vec100 = squeeze(matADS(intR,6,:))' - matAPS(intR,:);
			matEffectSize(intR,6,:) = vec100./(vec100+vec0);
		end
		[h,x,matAP_corr]=fdr_bh(matAP);
		%matAP_corr
		
		subplot(2,2,4)
		errorbar([vecPlotC(1) vecPlotC(end)],vecMeanPre(1)*[1 1],vecSEPre(1)*[1 1],'color',[0 1 0],'linestyle','--');
		hold on
		errorbar([vecPlotC(1) vecPlotC(end)],vecMeanPre(2)*[1 1],vecSEPre(2)*[1 1],'color',[1 0 0],'linestyle','--');
		errorfill(vecPlotC,matMean(2,:),matSE(2,:),[1 0 0],[1 0.5 0.5]);
		errorfill(vecPlotC,matMean(1,:),matSE(1,:),[0 1 0],[0.5 1 0.5]);
		hold off
		set(gca,'XScale','log')
		set(gca,'Xtick',vecPlotC,'xticklabel',vecContrasts)
		ylabel('Mean AEs per second assembly occurrence')
		xlabel('Stimulus contrast (%)')
		ylim([0 max(get(gca,'ylim'))]);
		%title([strRecur ' ' strGood '/' strBad sprintf('; regP=%.3f,%.3f',vecRegP)]);
		title([strRecur ' ' strGood '/' strBad]);
		
		
	end
end

%% save figs
if boolSavePlots
	figure(hRealAssemblyCorrelates);
	drawnow;
	jFig = get(handle(gcf), 'JavaFrame');
	jFig.setMaximized(true);
	figure(gcf);
	drawnow;
	intFigCounter = intFigCounter + 1;
	strFig = sprintf('Meta%s%d_assembly_sigmoids_RecurHitMiss_raw',strNPS,intFigCounter);
	export_fig([strFig '.tif']);
	export_fig([strFig '.pdf']);
	
	figure(hNonRecAssemblyCorrelates);
	drawnow;
	jFig = get(handle(gcf), 'JavaFrame');
	jFig.setMaximized(true);
	figure(gcf);
	drawnow;
	intFigCounter = intFigCounter + 1;
	strFig = sprintf('Meta%s%d_assembly_sigmoids_NonRecurHitMiss_raw',strNPS,intFigCounter);
	export_fig([strFig '.tif']);
	export_fig([strFig '.pdf']);
	
	figure(hRealAssemblyCorrelatesSlowFast);
	drawnow;
	jFig = get(handle(gcf), 'JavaFrame');
	jFig.setMaximized(true);
	figure(gcf);
	drawnow;
	intFigCounter = intFigCounter + 1;
	strFig = sprintf('Meta%s%d_assembly_sigmoids_RecurFastSlow_raw',strNPS,intFigCounter);
	export_fig([strFig '.tif']);
	export_fig([strFig '.pdf']);
	
	figure(hNonRecAssemblyCorrelatesSlowFast);
	drawnow;
	jFig = get(handle(gcf), 'JavaFrame');
	jFig.setMaximized(true);
	figure(gcf);
	drawnow;
	intFigCounter = intFigCounter + 1;
	strFig = sprintf('Meta%s%d_assembly_sigmoids_NonRecurFastSlow_raw',strNPS,intFigCounter);
	export_fig([strFig '.tif']);
	export_fig([strFig '.pdf']);
end

%% reaction-time dependency
vecX = [0 3];
hRTDependenceFig = figure;
intActTypes = size(cellRTDependency,1)-1;
matR2 = nan(intActTypes,intAnimals);
vecPs = nan(1,size(cellRTDependency,1)-1);
for intActType=1:intActTypes
	figure(hRTDependenceFig);
	if intActType == 1
		%heterogeneity
		strLabelY = 'dF/F0';
		strC = 'k';
	elseif intActType == 2
		%z-scored act
		strLabelY = 'Activation events';
		strC = 'b';
	elseif intActType == 3
		%dF/F act
		strLabelY = 'Recurring assembly presence';
		strC = 'r';
	elseif intActType == 4
		%dF/F act
		strLabelY = 'Non-recurring assembly presence';
		strC = 'r';
	end
	
	matReg = zeros(intAnimals,length(vecX));
	vecSlopes = zeros(intAnimals,1);
	vecR2 = zeros(intAnimals,1);
	subplot(2,intActTypes,intActType);
	hold on;
	vecAggRTs = [];
	vecAggAct = [];
	for intPopulation=1:intAnimals
		vecAggRTs = [vecAggRTs cellRTDependency{1,intPopulation}];
		vecAggAct = [vecAggAct cellRTDependency{intActType+1,intPopulation}];
		
		%do regressions per population
		sStatsC=regstats(cellRTDependency{intActType+1,intPopulation}',cellRTDependency{1,intPopulation},'linear');
		vecY = polyval(sStatsC.beta([2 1]),vecX);
		matReg(intPopulation,:) = vecY;
		vecSlopes(intPopulation) = sStatsC.beta(2);
		vecR2(intPopulation) = sStatsC.rsquare;
		plot(vecX,vecY,'Color',[0.5 0.5 1.0]);
	end
	matR2(intActType,:) = vecR2;
	
	%plot means
	vecMeanY = mean(matReg,1);
	plot(vecX,vecMeanY,'Color',strC,'LineWidth',2)
	
	%ttest
	[h,dblPSlope] = ttest(vecSlopes);
	[corrected_p, h]=bonf_holm([ones([1 intActTypes*2-1]) dblPSlope],0.05);
	dblPSlopeCorr = corrected_p(end);
	vecPs(intActType) = dblPSlope;
	title(sprintf('Mean of linear regressions over animals; slope-p=%.3f',dblPSlope))
	xlabel('Reaction Time (s)')
	ylabel(strLabelY)
	xlim([0 3]);
end

matP = nan(intActTypes,intActTypes);
for intActType1=1:(intActTypes-1)
	for intActType2=(intActType1+1):intActTypes
		[h,matP(intActType2,intActType1)] = ttest(matR2(intActType1,:),matR2(intActType2,:));
	end
end
vecP = matP(~isnan(matP));
[h,crit,vecFDR]=fdr_bh(matP(~isnan(matP)),0.05);

matV = ~isnan(matP);
matP_FDR = nan(size(matP));
intCounter=0;
for intActType1=1:(intActTypes-1)
	for intActType2=(intActType1+1):intActTypes
		if matV(intActType2,intActType1) == 1
			intCounter = intCounter + 1;
			matP_FDR(intActType2,intActType1) = vecFDR(intCounter);
		end
	end
end
matP_FDR=matP_FDR';
matP = matP';
%return

if boolSavePlots
	figure(hRTDependenceFig);
	intFigCounter = intFigCounter + 1;
	drawnow;
	jFig = get(handle(gcf), 'JavaFrame');
	jFig.setMaximized(true);
	figure(gcf);
	drawnow;
	strFig = sprintf('Meta%d_assembly_activity_RTDependency_raw',intFigCounter);
	export_fig([strFig '.tif']);
	export_fig([strFig '.pdf']);
end

%% clean up
cd(strOldDir);