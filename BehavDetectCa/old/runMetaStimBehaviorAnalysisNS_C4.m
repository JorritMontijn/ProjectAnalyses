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
strFigDir = 'D:\Data\Results\stimdetection\metaNS_Supp2';
strDataDir = 'D:\Data\Results\stimdetection';
cellInclude = {'20140207','20140314','20140425','20140507','20140530','20140604','20140711','20140715'};%exclude:20140430
boolSavePlots = true;
intFigCounter = 0;

%% pre-allocate
intAnimal = 0;
vecAnimalSource = [];
vecSubPopulationSource = [];

matContAct = [];
cellHetTime = {};
cellActTime = {};
matQuintileRemoval = []; %contrast x quintile x hit/miss x animal
cellHetAll = {};
cellPupAll = {};
cellRTAll = {};
cellContAll = {};
		
%% create aggregate data structures
sDir = dir(strDataDir);
strOldDir = cd(strDataDir);
intCounterF1 = 0;
intCounterF2 = 0;

for intFile=1:numel(sDir)
	%check if data part 1
	strFile = sDir(intFile).name;
	strRec = getFlankedBy(strFile,'data_aggregateNS_Supp2_','_');
	intRec = find(strcmp(strRec,cellInclude),1);
	if ~isempty(intRec)
		%load file
		load(strFile);
		
		%update population source
		intNrPops = size(cellSaveMatContAct,1);
		intAnimal = intAnimal + 1;
		vecAnimalSource(end+1) = intAnimal;
		vecSubPopulationSource(end+1) = intAnimal;
		if intNrPops > 1 %group populations
			%take mean over matrices
			cellSaveMatContAct = {nanmean(cat(3,cellSaveMatContAct{1},cellSaveMatContAct{2}),3)};
			cellSaveMatDetect = {nanmean(cat(3,cellSaveMatDetect{1},cellSaveMatDetect{2}),3)};
			
			cellSaveHetQs = {cat(2,cellSaveHetQs{1},cellSaveHetQs{2})}; 
			cellSaveHetAll = {cat(2,cellSaveHetAll{1},cellSaveHetAll{2})};
			cellSavePupAll = {cat(2,cellSavePupAll{1},cellSavePupAll{2})};
			cellSaveRespTimeTrials = {cat(2,cellSaveRespTimeTrials{1},cellSaveRespTimeTrials{2})};
			cellSaveContrastTrials = {cat(2,cellSaveContrastTrials{1},cellSaveContrastTrials{2})};
		end
		
		%quintile pre-pro
		matHet = [cellSaveHetQs{1};cellSaveHetAll{1}];
		matQ = nan(6,6,2);%contrast x quintile x hit/miss
		vecC = unique(cellSaveContrastTrials{1});
		for intContrast=1:6
			indC = cellSaveContrastTrials{1} == vecC(intContrast);
			for intQuintile=1:6
				for intHitMiss=1:2
					if intHitMiss == 1,indR = ~isnan(cellSaveRespTimeTrials{1});
					else indR = isnan(cellSaveRespTimeTrials{1});end
					if sum(indC&indR) == 0
						vecR=find(indC);
						[dummy,vecI]=findmin(cellSaveRespTimeTrials{1}(vecR),2);
						cellSaveRespTimeTrials{1}(vecR(vecI)) = nan;
						if intHitMiss == 1,indR = ~isnan(cellSaveRespTimeTrials{1});
						else indR = isnan(cellSaveRespTimeTrials{1});end
					end
					matQ(intContrast,intQuintile,intHitMiss) = mean(matHet(intQuintile,indC&indR));
				end
			end
		end
		
		%assign data
		intCounterF1 = intCounterF1 + 1;
		cellHetAll{intCounterF1} = cellSaveHetAll{1};
		cellPupAll{intCounterF1} = cellSavePupAll{1};
		cellRespTimeTrials{intCounterF1} = cellSaveRespTimeTrials{1};
		cellContrastTrials{intCounterF1} = cellSaveContrastTrials{1};
		
		matContAct = cat(3,matContAct,cellSaveMatContAct{1});
		matQuintileRemoval = cat(4,matQuintileRemoval,matQ);
	end
	
	%check if data part 2
	strFile = sDir(intFile).name;
	strRec = getFlankedBy(strFile,'data_aggregateNS','_');
	intRec = find(strcmp(strRec,cellInclude),1);
	if ~isempty(intRec)
		%load file
		load(strFile);
		
		%update population source
		intNrPops = size(cellSaveHetTime,2);
		if intNrPops > 1 %group populations
			%concatenate heterogen-time
			for intC=1:size(cellSaveHetTime,1)
				cellSaveHetTime{intC,1} = cat(1,cellSaveHetTime{intC,1},cellSaveHetTime{intC,2});
				cellSaveActTime{intC,1} = cat(1,cellSaveActTime{intC,1},cellSaveActTime{intC,2});
			end
			cellSaveHetTime = cellSaveHetTime(:,1);
			cellSaveActTime = cellSaveActTime(:,1);
		end
		
		%assign data
		cellHetTime = cat(2,cellHetTime,cellSaveHetTime);
		cellActTime = cat(2,cellActTime,cellSaveActTime);
		intCounterF2 = intCounterF2 + 1;
	end
end

%check for inconsistencies
%if ~((intCounterF1 == intCounterF2) && (intCounterF2 == intCounterF3) && (intCounterF3 == intCounterF4) && (intCounterF4 == intCounterF5))
%	warning([mfilename ':InconsistentFileNumbers'],'File counters inconsistent; F1=%d; F2=%d; F3=%d; F4=%d; F5=%d',intCounterF1,intCounterF2,intCounterF3,intCounterF4,intCounterF5)
%end

%% meta analyses part 1
cd(strFigDir);
%close all;
intAnimals = intCounterF1;

%% dF/F over contrast for hit/miss and still/move
vecAnimals = [1:8];
matContActNorm = nan(size(matContAct));
for intAnimal=vecAnimals
	matTempAct = matContAct(:,:,intAnimal);
	matContActNorm(:,:,intAnimal) = matTempAct;% ./ mean(matTempAct(:));
end
hActCon=figure;

%pre-compute variables
vecContrasts = [0.2 0.5 2 8 32 100];
vecWindow = [1 length(vecContrasts)];
vecWindowSelect = vecWindow(1):vecWindow(end);
intWL = length(vecWindowSelect);
vecWindowInv = intWL:-1:1;
vecWindowPlotInv = vecWindow(end):-1:vecWindow(1);
vecX = [vecContrasts(vecWindowSelect) vecContrasts(vecWindowPlotInv)];
vecLineX = vecContrasts(vecWindowSelect);

%plot dF/F over contrasts
matMES = nan(3,intAnimals);
for intActType=1:3
	subplot(2,2,intActType)
	if intActType==1
		strActType = 'dF/F';
		matHitAct = matContActNorm(:,1,:);
		matMissAct = matContActNorm(:,2,:);
		vecLimY = [-0.01 0.06];
	elseif intActType==2
		strActType = 'Heterogeneity';
		matHitAct = matContActNorm(:,3,:);
		matMissAct = matContActNorm(:,4,:);
		%for intAnimal=1:size(matHitAct,3)
		%	matTempActH = matHitAct(:,:,intAnimal);
		%	matTempActM = matMissAct(:,:,intAnimal);
		%	matHitAct(:,:,intAnimal) = matTempActH ./ mean([matTempActH(:);matTempActM(:)]);
		%	matMissAct(:,:,intAnimal) = matTempActM ./ mean([matTempActH(:);matTempActM(:)]);
		%end
		vecLimY = [0.7 1.1];
	elseif intActType==3
		strActType = 'PupilSize';
		matHitAct = matContActNorm(:,5,vecAnimals);
		matMissAct = matContActNorm(:,6,vecAnimals);
		%for intAnimal=1:size(matHitAct,3)
		%	matTempActH = matHitAct(:,:,intAnimal);
		%	matTempActM = matMissAct(:,:,intAnimal);
		%	matHitAct(:,:,intAnimal) = matTempActH ./ mean([matTempActH(:);matTempActM(:)]);
		%	matMissAct(:,:,intAnimal) = matTempActM ./ mean([matTempActH(:);matTempActM(:)]);
		%end
		vecLimY = 'auto';
	end
	
	vecHitY = nanmean(matHitAct,3)';
	vecHitE = nanstd(matHitAct,[],3)'/sqrt(size(matContActNorm,3));
	vecMissY = nanmean(matMissAct,3)';
	vecMissE = nanstd(matMissAct,[],3)'/sqrt(size(matContActNorm,3));
	
	%get sig
	vecP = zeros(1,size(matContActNorm,1));
	for intC=1:size(matContActNorm,1)
		[h,p,ci] = ttest(matHitAct(intC,1,:),matMissAct(intC,1,:));
		
		%put in vector
		vecP(intC) = p;
	end
	
	%overall t-test
	matActHits = matHitAct(2:5,1,:);
	matActMisses = matMissAct(2:5,1,:);
	[h,dblP_All,ci] = ttest(matActHits(:),matActMisses(:));
	
	%get data
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
	dblFracMoved = mean(mean(matContAct(:,6,:)));
	title(sprintf('Pop act during stim; p-vals: 0: %.3f; 0.5: %.3f; 2: %.3f; 8: %.3f; 32: %.3f; 100: %.3f; Overall,p=%.3f',[vecP dblP_All]))
	grid on
	xlabel('Stimulus contrast (%)');
	ylabel(sprintf('%s during stimulus',strActType));
	xlim(vecContrasts(vecWindow))
	ylim(vecLimY);
	legend({'SEM','Miss','SEM','Hit'},'Location','Best')
end

vecMES_dFoF = mean(squeeze(matContActNorm(2:5,7,:)),1); %dFoF
vecMES_Het = mean(squeeze(matContActNorm(2:5,8,:)),1); %het
[h,p] = ttest(vecMES_dFoF,vecMES_Het);
fprintf('Measure of effect size dFoF (mean=%.3f) vs. heterogeneity (mean=%.3f); p=%.3f\n',mean(vecMES_dFoF),mean(vecMES_Het),p);

%plot heterogeneity/pupil size dependency
subplot(2,2,4);
vecIntercepts = nan(1,intAnimals);
vecSlopes = nan(1,intAnimals);
for intAnimal=vecAnimals
	indResp = ~isnan(cellRespTimeTrials{intAnimal});
	
	vecRespPupSize = cellPupAll{intAnimal}(indResp);
	vecNoRespPupSize = cellPupAll{intAnimal}(~indResp);
	
	vecRespHetero = cellHetAll{intAnimal}(indResp);
	vecNoRespHetero = cellHetAll{intAnimal}(~indResp);
	
	hold on
	vecDataY = [vecRespPupSize; vecNoRespPupSize];
	vecDataX = [vecRespHetero vecNoRespHetero];
	sStats=regstats(vecDataY,vecDataX,'linear');
	vecX = [min(vecDataX(:)) max(vecDataX(:))];
	vecY = polyval(sStats.beta([2 1]),vecX);
	plot(vecX,vecY,'Color',[0.5 0.5 0.5])
	
	hold off
	vecIntercepts(intAnimal) = sStats.beta(1);
	vecSlopes(intAnimal) = sStats.beta(2);
end
hold on
vecX = [0.5 2];
vecY = polyval([nanmean(vecSlopes) nanmean(vecIntercepts)],vecX);
plot(vecX,vecY,'Color','b');
hold off
ylim([-0.6 0.6])
[h,p]=ttest(vecSlopes);
title(sprintf('Pupil size vs heterogeneity; p-val slopes=%.3f',p))
xlabel('Heterogeneity  during stimulus')
ylabel('Normalized pupil size')

drawnow;
strFigTitle = 'eyetracking';
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

%% clean up
cd(strOldDir);