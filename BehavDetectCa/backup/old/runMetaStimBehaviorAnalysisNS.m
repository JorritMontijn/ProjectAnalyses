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
strFigDir = 'D:\Data\Results\stimdetection\metaNS';
strDataDir = 'D:\Data\Results\stimdetection';
cellInclude = {'20140207','20140314','20140425','20140507','20140530','20140604','20140711','20140715'};%exclude:20140430
boolSavePlots = true;
intFigCounter = 0;

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

cellMatrices = {};%: blob matrices in [1 x 4 (detect/no-detect/diff/significance)] cell array with [5 (c) x 31 (bins)] matrix
matCorrs = [];
cellNormActDissim = {};%: within-group z-scored activation dissimilarity [2 (high/low HCAR) x 2 (detect/no-detect)] x n (animals/blocks) with every cell = vector of dissimilarity values
cellCorrITD = {};%: inter-trial-distance dependence of assembly consistency [2 (hit/miss) x 3 (ITD/HCAR/non-HCAR) x n (animals/blocks)] with every cell = vector with correlation values
cellRTDependency = {};%: RT dependence of activity dissimilarity [3 (RT/HCAR/non-HCAR) x n (animals/blocks)] with every cell = vector of trials
matMixBootstrappedDecodingOutput = [];
matContHetero = [];
cellHetTime = {};
cellActTime = {};
matRespDecode = [];
matRespDecodeDist = [];
cellHetPoprespCorr = {};
cellAcrossRepCorr = {};

%% create aggregate data structures
sDir = dir(strDataDir);
strOldDir = cd(strDataDir);
intCounterF1 = 0;
intCounterF2 = 0;
intCounterF3 = 0;
intCounterF4 = 0;
intCounterF5 = 0;

for intFile=1:numel(sDir)
	%check if data part 1
	strFile = sDir(intFile).name;
	strRec = getFlankedBy(strFile,'data_aggregate','_');
	intRec = find(strcmp(strRec,cellInclude),1);
	if ~isempty(intRec)
		%load file
		load(strFile);
		
		%assign data
		cellBehavDetect = cat(2,cellBehavDetect,cellSaveBehavDetect);%: behavioral detection [6 (c) x n (animals)] cell array with every cell = vector of detect/no detect values (1/0)
		cellBehavRT = cat(2,cellBehavRT,cellSaveBehavRT);%: behavioral RT [5 (c) x n (animals)] cell array with every cell = vector of RTs
		cellPreStim = cat(3,cellPreStim,cellSavePreStim);%: dF/F detect/no-detect preceding stimulus frames [5 (c) x 2 (d/nd) cell array] x n (animals) with every cell = vector of values
		cellBaseSubtrStimAct = cat(3,cellBaseSubtrStimAct,cellSaveBaseSubtrStimAct);%: dF/F detect/no-detect baseline-subtracted activity [5 (c) x 2 (d/nd) cell array] x n (animals) with every cell = vector of values
		cellDuringStim = cat(3,cellDuringStim,cellSaveDuringStim);%: dF/F detect/no-detect during stimulus frames [5 (c) x 2 (d/nd) cell array] x n (animals) with every cell = vector of values
		cellTraceAct = cat(4,cellTraceAct,cellSaveTraceAct);%: dF/F trace detect/no-detect pref/non-pref [2 (pref/no-pref) x 2 (detect/no-detect) x 6 (c)] x n (animals) with every cell = [Nreps x Nframes] matrix of trace dF/F values
		intCounterF1 = intCounterF1 + 1;
	end
	
	%check if data part 2
	strFile = sDir(intFile).name;
	strRec = getFlankedBy(strFile,'data_aggregateNS','_');
	intRec = find(strcmp(strRec,cellInclude),1);
	if ~isempty(intRec)
		%load file
		load(strFile);
		
		%update population source
		intNrPops = size(cellSaveMatrices,2);
		intAnimal = intAnimal + 1;
		vecAnimalSource(end+1) = intAnimal;
		vecSubPopulationSource(end+1) = intAnimal;
		if intNrPops > 1 %group populations
			%take mean over matrices
			cellSaveMatrices{1,1} = nanmean(cat(3,cellSaveMatrices{1,1},cellSaveMatrices{1,2}),3);
			cellSaveMatrices{2,1} = nanmean(cat(3,cellSaveMatrices{2,1},cellSaveMatrices{2,2}),3);
			cellSaveMatrices{3,1} = nanmean(cat(3,cellSaveMatrices{3,1},cellSaveMatrices{3,2}),3);
			cellSaveMatrices{4,1} = cellSaveMatrices{4,1} | cellSaveMatrices{4,2};
			cellSaveMatrices(:,2) = [];
			
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
			
			%concatenate ITD correlations
			for intDim1=1:size(cellSaveCorrITD,1)
				for intDim2=1:size(cellSaveCorrITD,2)
					for intDim3=1:size(cellSaveCorrITD,3)
						cellSaveCorrITD{intDim1,intDim2,intDim3,1} = [cellSaveCorrITD{intDim1,intDim2,intDim3,1} cellSaveCorrITD{intDim1,intDim2,intDim3,2}];
					end
				end
			end
			cellSaveCorrITD(:,:,:,2) = [];
			
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
			
			%concatenate across rep corr
			cellSaveAcrossRepCorr{3,1} = [cellSaveAcrossRepCorr{3,1} cellSaveAcrossRepCorr{3,2}];
			cellSaveAcrossRepCorr{4,1} = [cellSaveAcrossRepCorr{4,1} cellSaveAcrossRepCorr{4,2}];
			cellSaveAcrossRepCorr{5,1} = [cellSaveAcrossRepCorr{5,1} cellSaveAcrossRepCorr{5,2}];
			cellSaveAcrossRepCorr{6,1}.beta = mean([cellSaveAcrossRepCorr{6,1}.beta cellSaveAcrossRepCorr{6,2}.beta],2);

		else
			%mean of correlations
			matSelect = tril(true(size(cellSaveStimCorrs{1,1})),-1);
			
			vecCorrs = nan(size(cellSaveStimCorrs,1),1);
			for intType=1:size(cellSaveStimCorrs,1)
				vecCorrs(intType) = mean(cellSaveStimCorrs{intType,1}(matSelect));
			end
		end
		%retrieve rep corr data
		cellSaveAcrossRepCorr2 = {};
		cellSaveAcrossRepCorr2{1,1} = cellSaveAcrossRepCorr{3,1};
		cellSaveAcrossRepCorr2{2,1} = cellSaveAcrossRepCorr{4,1};
		cellSaveAcrossRepCorr2{3,1} = cellSaveAcrossRepCorr{5,1};
		cellSaveAcrossRepCorr2{4,1} = cellSaveAcrossRepCorr{6,1}.beta;

		
		%assign data
		cellMatrices = cat(2,cellMatrices,cellSaveMatrices);%: blob matrices in [1 x 4 (detect/no-detect/diff/significance)] cell array with [5 (c) x 31 (bins)] matrix
		matCorrs = cat(2,matCorrs,vecCorrs);%: signal correlations [3 (within-HCAR, between, within-non-HCAR) x n (animals/blocks)] with every cell = vector of pairwise correlation values
		cellNormActDissim = cat(3,cellNormActDissim,cellSaveNormActDissim);%: within-group z-scored activation dissimilarity [2 (high/low HCAR) x 2 (detect/no-detect)] x n (animals/blocks) with every cell = vector of dissimilarity values
		cellCorrITD = cat(4,cellCorrITD,cellSaveCorrITD);%: inter-trial-distance dependence of assembly consistency [2 (high/low HCAR) x 2 (hits/misses) x n (animals/blocks)] with every cell = [p x 2] matrix with intertrial distance and correlation values for each trial pair
		cellRTDependency = cat(2,cellRTDependency,cellSaveRTDependency);%: RT dependence of activity dissimilarity [3 (RT/HCAR/non-HCAR) x n (animals/blocks)] with every cell = vector of trials
		matContHetero = cat(3,matContHetero,cellSaveMatContHetero{1});
		cellHetTime = cat(2,cellHetTime,cellSaveHetTime);
		cellActTime = cat(2,cellActTime,cellSaveActTime);
		matRespDecode = cat(4,matRespDecode,cellSaveRespDecode{1,1});
		matRespDecodeDist = cat(3,matRespDecodeDist,cellSaveRespDecode{1,2});
		cellHetPoprespCorr = cat(2,cellHetPoprespCorr,cellSaveHetPoprespCorr);
		cellAcrossRepCorr = cat(2,cellAcrossRepCorr,cellSaveAcrossRepCorr2);
		intCounterF2 = intCounterF2 + 1;
	end
	
	%check for decoding
	strFile = sDir(intFile).name;
	strRec = getFlankedBy(strFile,'data_aggregateNSSD','_');
	intRec = find(strcmp(strRec,cellInclude),1);
	if ~isempty(intRec) && length(strRec) == 8
		%load file
		load(strFile);
		
		%update population source
		intNrPops = numel(cellSaveDecoding)/2;
		if intNrPops > 1 %group populations
			%take mean over matrices
			cellSaveDecoding(1).matMixBootstrappedDecodingOutput = nanmean(cat(4,cellSaveDecoding(1).matMixBootstrappedDecodingOutput,cellSaveDecoding(2).matMixBootstrappedDecodingOutput),4);
			cellSaveDecoding(2) = [];
		end
		
		%assign data
		matMixBootstrappedDecodingOutput = cat(1,matMixBootstrappedDecodingOutput,nanmean(cellSaveDecoding(1).matMixBootstrappedDecodingOutput,1));
		intCounterF3 = intCounterF3 + 1;
	end
	
	%check for decoding 2
	strFile = sDir(intFile).name;
	strRec = getFlankedBy(strFile,'data_aggregateNSSD2_','_');
	intRec = find(strcmp(strRec,cellInclude),1);
	if ~isempty(intRec)
		%load file
		load(strFile);
		
		%update population source
		intNrPops = numel(cellSaveDecoding2);
		if intNrPops > 1 %group populations
			%take mean over matrices
			cellSaveDecoding2(1).matHD_Orientation = cat(1,cellSaveDecoding2(1).matHD_Orientation,cellSaveDecoding2(2).matHD_Orientation);
			cellSaveDecoding2(1).matHD_Contrast = cat(1,cellSaveDecoding2(1).matHD_Contrast,cellSaveDecoding2(2).matHD_Contrast);
			cellSaveDecoding2(1).matHD_ReactionTime = cat(1,cellSaveDecoding2(1).matHD_ReactionTime,cellSaveDecoding2(2).matHD_ReactionTime);
			cellSaveDecoding2(1).matHD_DecodingAccuracy = cat(1,cellSaveDecoding2(1).matHD_DecodingAccuracy,cellSaveDecoding2(2).matHD_DecodingAccuracy);
			cellSaveDecoding2(1).matHD_HeterogeneityPre = cat(1,cellSaveDecoding2(1).matHD_HeterogeneityPre,cellSaveDecoding2(2).matHD_HeterogeneityPre);
			cellSaveDecoding2(1).matHD_HeterogeneityDuring = cat(1,cellSaveDecoding2(1).matHD_HeterogeneityDuring,cellSaveDecoding2(2).matHD_HeterogeneityDuring);
			cellSaveDecoding2(1).matHD_ActivityPre = cat(1,cellSaveDecoding2(1).matHD_ActivityPre,cellSaveDecoding2(2).matHD_ActivityPre);
			cellSaveDecoding2(1).matHD_ActivityDuring = cat(1,cellSaveDecoding2(1).matHD_ActivityDuring,cellSaveDecoding2(2).matHD_ActivityDuring);
			cellSaveDecoding2(2) = [];
		end
		
		%assign data
		cellHD_Orientation{intCounterF4+1} = cellSaveDecoding2(1).matHD_Orientation;
		cellHD_Contrast{intCounterF4+1} = cellSaveDecoding2(1).matHD_Contrast;
		cellHD_ReactionTime{intCounterF4+1} = cellSaveDecoding2(1).matHD_ReactionTime;
		cellHD_DecodingAccuracy{intCounterF4+1} = cellSaveDecoding2(1).matHD_DecodingAccuracy;
		cellHD_HeterogeneityPre{intCounterF4+1} = cellSaveDecoding2(1).matHD_HeterogeneityPre;
		cellHD_HeterogeneityDuring{intCounterF4+1} = cellSaveDecoding2(1).matHD_HeterogeneityDuring;
		cellHD_ActivityPre{intCounterF4+1} = cellSaveDecoding2(1).matHD_ActivityPre;
		cellHD_ActivityDuring{intCounterF4+1} = cellSaveDecoding2(1).matHD_ActivityDuring;
		
		intCounterF4 = intCounterF4 + 1;
	end
	
	%check for decoding 3
	strFile = sDir(intFile).name;
	strRec = getFlankedBy(strFile,'data_aggregateNSSD3_','_');
	intRec = find(strcmp(strRec,cellInclude),1);
	if ~isempty(intRec)
		%load file
		%{
		load(strFile);
		
		%update population source
		intNrPops = size(cellSaveDecoding3,2);
		if intNrPops > 1 %group populations
			%take mean over matrices
			cellSaveDecoding3(1,1).matHD_Orientation = cat(1,cellSaveDecoding3(1,1).matHD_Orientation,cellSaveDecoding3(1,2).matHD_Orientation);
			cellSaveDecoding3(1,1).matHD_Contrast = cat(1,cellSaveDecoding3(1,1).matHD_Contrast,cellSaveDecoding3(1,2).matHD_Contrast);
			cellSaveDecoding3(1,1).matHD_ReactionTime = cat(1,cellSaveDecoding3(1,1).matHD_ReactionTime,cellSaveDecoding3(1,2).matHD_ReactionTime);
			cellSaveDecoding3(1,1).matHD_DecodingAccuracy = cat(1,cellSaveDecoding3(1,1).matHD_DecodingAccuracy,cellSaveDecoding3(1,2).matHD_DecodingAccuracy);
			cellSaveDecoding3(1,1).matHD_HeterogeneityPre = cat(1,cellSaveDecoding3(1,1).matHD_HeterogeneityPre,cellSaveDecoding3(1,2).matHD_HeterogeneityPre);
			cellSaveDecoding3(1,1).matHD_HeterogeneityDuring = cat(1,cellSaveDecoding3(1,1).matHD_HeterogeneityDuring,cellSaveDecoding3(1,2).matHD_HeterogeneityDuring);
			cellSaveDecoding3(1,1).matHD_ActivityPre = cat(1,cellSaveDecoding3(1,1).matHD_ActivityPre,cellSaveDecoding3(1,2).matHD_ActivityPre);
			cellSaveDecoding3(1,1).matHD_ActivityDuring = cat(1,cellSaveDecoding3(1,1).matHD_ActivityDuring,cellSaveDecoding3(1,2).matHD_ActivityDuring);
			
			cellSaveDecoding3(2,1).matHD_Orientation = cat(1,cellSaveDecoding3(2,1).matHD_Orientation,cellSaveDecoding3(2,2).matHD_Orientation);
			cellSaveDecoding3(2,1).matHD_Contrast = cat(1,cellSaveDecoding3(2,1).matHD_Contrast,cellSaveDecoding3(2,2).matHD_Contrast);
			cellSaveDecoding3(2,1).matHD_ReactionTime = cat(1,cellSaveDecoding3(2,1).matHD_ReactionTime,cellSaveDecoding3(2,2).matHD_ReactionTime);
			cellSaveDecoding3(2,1).matHD_DecodingAccuracy = cat(1,cellSaveDecoding3(2,1).matHD_DecodingAccuracy,cellSaveDecoding3(2,2).matHD_DecodingAccuracy);
			cellSaveDecoding3(2,1).matHD_HeterogeneityPre = cat(1,cellSaveDecoding3(2,1).matHD_HeterogeneityPre,cellSaveDecoding3(2,2).matHD_HeterogeneityPre);
			cellSaveDecoding3(2,1).matHD_HeterogeneityDuring = cat(1,cellSaveDecoding3(2,1).matHD_HeterogeneityDuring,cellSaveDecoding3(2,2).matHD_HeterogeneityDuring);
			cellSaveDecoding3(2,1).matHD_ActivityPre = cat(1,cellSaveDecoding3(2,1).matHD_ActivityPre,cellSaveDecoding3(2,2).matHD_ActivityPre);
			cellSaveDecoding3(2,1).matHD_ActivityDuring = cat(1,cellSaveDecoding3(2,1).matHD_ActivityDuring,cellSaveDecoding3(2,2).matHD_ActivityDuring);
			cellSaveDecoding3(:,2) = [];
		end
		
		%assign data
		cellHD_Orientation2{1,intCounterF5+1} = cellSaveDecoding3(1).matHD_Orientation;
		cellHD_Contrast2{1,intCounterF5+1} = cellSaveDecoding3(1).matHD_Contrast;
		cellHD_ReactionTime2{1,intCounterF5+1} = cellSaveDecoding3(1).matHD_ReactionTime;
		cellHD_DecodingAccuracy2{1,intCounterF5+1} = cellSaveDecoding3(1).matHD_DecodingAccuracy;
		cellHD_HeterogeneityPre2{1,intCounterF5+1} = cellSaveDecoding3(1).matHD_HeterogeneityPre;
		cellHD_HeterogeneityDuring2{1,intCounterF5+1} = cellSaveDecoding3(1).matHD_HeterogeneityDuring;
		cellHD_ActivityPre2{1,intCounterF5+1} = cellSaveDecoding3(1).matHD_ActivityPre;
		cellHD_ActivityDuring2{1,intCounterF5+1} = cellSaveDecoding3(1).matHD_ActivityDuring;
		
		cellHD_Orientation2{2,intCounterF5+1} = cellSaveDecoding3(2).matHD_Orientation;
		cellHD_Contrast2{2,intCounterF5+1} = cellSaveDecoding3(2).matHD_Contrast;
		cellHD_ReactionTime2{2,intCounterF5+1} = cellSaveDecoding3(2).matHD_ReactionTime;
		cellHD_DecodingAccuracy2{2,intCounterF5+1} = cellSaveDecoding3(2).matHD_DecodingAccuracy;
		cellHD_HeterogeneityPre2{2,intCounterF5+1} = cellSaveDecoding3(2).matHD_HeterogeneityPre;
		cellHD_HeterogeneityDuring2{2,intCounterF5+1} = cellSaveDecoding3(2).matHD_HeterogeneityDuring;
		cellHD_ActivityPre2{2,intCounterF5+1} = cellSaveDecoding3(2).matHD_ActivityPre;
		cellHD_ActivityDuring2{2,intCounterF5+1} = cellSaveDecoding3(2).matHD_ActivityDuring;
		%}
		intCounterF5 = intCounterF5 + 1;
	end
end

%check for inconsistencies
if ~((intCounterF1 == intCounterF2) && (intCounterF2 == intCounterF3) && (intCounterF3 == intCounterF4) && (intCounterF4 == intCounterF5))
	warning([mfilename ':InconsistentFileNumbers'],'File counters inconsistent; F1=%d; F2=%d; F3=%d; F4=%d; F5=%d',intCounterF1,intCounterF2,intCounterF3,intCounterF4,intCounterF5)
end

%% meta analyses part 1
cd(strFigDir);
close all;
intAnimals = size(cellBehavDetect,2);
vecContrasts = [0 0.5 2 8 32 100];
%{
%% behavioral analyses
%% overview graph 0%/100% detection rates with significance calculation
matRTAgg = ones(5,intAnimals);
matDetectAgg = ones(6,intAnimals);
matDetect_C0 = ones(3,intAnimals); %upper, mean, lower
matDetect_C100 = ones(3,intAnimals);
for intAnimal=1:intAnimals
	%get data
	vecC0 = cellBehavDetect{1,intAnimal};
	vecC100 = cellBehavDetect{6,intAnimal};
	for intC=1:6
		matDetectAgg(intC,intAnimal) = mean(cellBehavDetect{intC,intAnimal});
		if intC>1,matRTAgg(intC-1,intAnimal) = mean(cellBehavRT{intC-1,intAnimal});end
	end
	
	%calculate confidence intervals using Clopper-Pearson method
	[dblP_C0,dblCI_C0] = binofit(sum(vecC0),length(vecC0));
	matDetect_C0(1,intAnimal) = dblCI_C0(2); %upper
	matDetect_C0(2,intAnimal) = dblP_C0; %mean
	matDetect_C0(3,intAnimal) = dblCI_C0(1); %lower
	[dblP_C100,dblCI_C100] = binofit(sum(vecC100),length(vecC100));
	matDetect_C100(1,intAnimal) = dblCI_C100(2); %upper
	matDetect_C100(2,intAnimal) = dblP_C100; %mean
	matDetect_C100(3,intAnimal) = dblCI_C100(1); %lower
end
hDetect_C100_C0 = figure;
errorbar(0.4:1:intAnimals,matDetect_C0(2,:),matDetect_C0(3,:)-matDetect_C0(2,:),matDetect_C0(1,:)-matDetect_C0(2,:),'LineStyle','none','Marker','x','Color','b')
hold on
errorbar(0.6:1:intAnimals,matDetect_C100(2,:),matDetect_C100(3,:)-matDetect_C100(2,:),matDetect_C100(1,:)-matDetect_C100(2,:),'LineStyle','none','Marker','x','Color','r')
hold off
legend('0% stimulus contrast','100% stimulus contrast')
ylim([0 1])
ylabel('Behavioral response probability')
xlabel('Animal')
set(gca,'XTick',0.5:1:intAnimals,'XTickLabel',1:intAnimals)
title('Response proportion with 95% confidence intervals using Clopper-Pearson method')
if boolSavePlots
	intFigCounter = intFigCounter + 1;
	drawnow;
	strFig = sprintf('Meta%d_behavPerformance_raw',intFigCounter);
	export_fig([strFig '.tif']);
	export_fig([strFig '.pdf']);
end

%% mean detection rate over contrasts (mean = blue, single animals = grey lines)
hDetect_AllContrasts = figure;
set(gca,'XScale','log','YScale','linear')
vecX = vecContrasts;
vecLabelsX = vecX;
vecX(1) = 0.25;
vecLabelsX(1) = 0;
semilogx(vecX(2:end),matDetectAgg(2:end,:),'-x','Color',[0.5 0.5 0.5])
hold on
semilogx(vecX(1),matDetectAgg(1,:),'-x','Color',[0.5 0.5 0.5])
semilogx(vecX(1),mean(matDetectAgg(1,:),2),'-x','Color',[0 0 1])
semilogx(vecX(2:end),mean(matDetectAgg(2:end,:),2),'-x','Color',[0 0 1])
hold off
set(gca,'XScale','log','YScale','linear')
set(gca,'XTick',vecX,'XTickLabel',vecLabelsX)
ylabel('Fraction correct')
ylim([0 1])
xlabel('Stimulus contrast (%)')
grid on


%regression analysis
%regression
matFitY = matDetectAgg(2:end,:)';
vecFitX = log(vecX(2:end));
matFitX = repmat(vecFitX,[size(matFitY,1) 1]);
stats = regstats(matFitY(:),matFitX(:),'linear',{'tstat'});

%nlinfit
[beta,resid,J,COVB,mse] = nlinfit(matFitX(:),matFitY(:),@linfit,[stats.tstat.beta(2) stats.tstat.beta(1)]);
[ypred,delta] = nlpredci(@linfit,vecFitX,beta,resid,'covar',COVB);

hold on
vecFillX = [vecX(2:end) vecX(end:-1:2)];
vecFillY = [ypred+delta' ypred(end:-1:1)-delta(end:-1:1)'];
fill(vecFillX,vecFillY,[0.5 0.5 0.5]);
semilogx(vecX(2:end),vecFitX * stats.tstat.beta(2) + stats.tstat.beta(1),'k');
hold off
title(sprintf('Behavioral stimulus detection per animal (grey) and mean (blue);regress p=%.9f',stats.tstat.pval(2)))
if boolSavePlots
	intFigCounter = intFigCounter + 1;
	drawnow;
	strFig = sprintf('Meta%d_behavDetect_raw',intFigCounter);
	export_fig([strFig '.tif']);
	export_fig([strFig '.pdf']);
end

%%  mean RT over contrasts (mean = blue, single animals = grey lines)
matRTs = matRTAgg*1000;
hRT_AllContrasts = figure;
set(gca,'XScale','log','YScale','linear')
vecX = vecContrasts(2:end);
vecLabelsX = vecX;
semilogx(vecX,matRTs,'-x','Color',[0.5 0.5 0.5])
hold on
semilogx(vecX,mean(matRTs,2),'-x','Color',[0 0 1])
hold off
set(gca,'XScale','log','YScale','linear')
set(gca,'XTick',vecX,'XTickLabel',vecLabelsX)
ylabel('Mean Reaction Time (ms)')
xlabel('Stimulus contrast (%)')
grid on

%regression analysis
%regression
matFitY = matRTs';
vecFitX = log(vecX);
matFitX = repmat(vecFitX,[size(matFitY,1) 1]);
stats = regstats(matFitY(:),matFitX(:),'linear',{'tstat'});

%nlinfit
[beta,resid,J,COVB,mse] = nlinfit(matFitX(:),matFitY(:),@linfit,[stats.tstat.beta(2) stats.tstat.beta(1)]);
[ypred,delta] = nlpredci(@linfit,vecFitX,beta,resid,'covar',COVB);

hold on
vecFillX = [vecX vecX(end:-1:1)];
vecFillY = [ypred+delta' ypred(end:-1:1)-delta(end:-1:1)'];
fill(vecFillX,vecFillY,[0.5 0.5 0.5]);
semilogx(vecX,vecFitX * stats.tstat.beta(2) + stats.tstat.beta(1),'k');
hold off
title(sprintf('Behavioral RT for hit trials per animal (grey) and mean (blue);regress p=%.9f',stats.tstat.pval(2)))

if boolSavePlots
	intFigCounter = intFigCounter + 1;
	drawnow;
	strFig = sprintf('Meta%d_behavRT_raw',intFigCounter);
	export_fig([strFig '.tif']);
	export_fig([strFig '.pdf']);
end

%% dF/F over contrasts for detect/no detect
%for pre-stim, during-stim, pre-subtracted-stim
%mean dF/F +/- st err over animals for detect/nodetect trials
for intType =1:3
	%get type
	if intType==1
		hActPreStim = figure;
		strFigTitle = 'ActPreStim';
		strLabelY = 'dF/F preceding stimulus';
		cellData = cellPreStim;
	elseif intType==2
		hActDurStim = figure;
		strFigTitle = 'ActDurStim';
		strLabelY = 'dF/F during stimulus';
		cellData = cellDuringStim;
	elseif intType==3
		hActPBSStim = figure;
		strFigTitle = 'ActBaSuStim';
		strLabelY = 'Baseline-subtracted dF/F during stimulus';
		cellData = cellBaseSubtrStimAct;
	end
	
	%get data
	matActHit = ones(5,intAnimals);
	matActMiss = ones(5,intAnimals);
	vecP = ones(5,1);
	
	for intC=1:size(cellPreStim,1)
		for intAnimal=1:intAnimals
			matActHit(intC,intAnimal) = mean(cellData{intC,1,intAnimal});
			matActMiss(intC,intAnimal) = mean(cellData{intC,2,intAnimal});
		end
		
		%ttest
		[h,vecP(intC),ci] = ttest(matActHit(intC,:),matActMiss(intC,:));
	end
	
	%overall t-test
	matHits = matActHit(1:4,:);
	matMisses = matActMiss(1:4,:);
	[h,dblP_All,ci] = ttest(matHits(:),matMisses(:));
	
	%pre-compute variables
	vecWindow = [2 6];
	vecWindowSelect = vecWindow(1):vecWindow(end);
	intWL = length(vecWindowSelect);
	vecWindowInv = intWL:-1:1;
	vecWindowPlotInv = vecWindow(end):-1:vecWindow(1);
	vecX = [vecContrasts(vecWindowSelect) vecContrasts(vecWindowPlotInv)];
	vecLineX = vecContrasts(vecWindowSelect);
	
	%get data
	for intDetect=[0 1]
		if intDetect == 1
			vecMeanTrace = mean(matActHit,2)';
			vecSE = std(matActHit,[],2)'/sqrt(size(matActHit,2));
			vecColorFill = [0.7 1.0 0.7];
			vecColorLine = [0 1 0];
		else
			vecColorLine = [1 0 0];
			vecColorFill = [1 0.7 0.7];
			vecMeanTrace = mean(matActMiss,2)';
			vecSE = std(matActMiss,[],2)'/sqrt(size(matActMiss,2));
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
	
	[h crit_p vecP_Adj]=fdr_bh(vecP,0.05,'pdep');
	
	set(gca,'XScale','log','YScale','linear')
	set(gca,'XTick',vecLineX,'XTickLabel',vecLineX)
	strTitle = 'T-test p: ';
	for intC=1:length(vecP_Adj)
		strTitle = [strTitle sprintf('C%.1f: ',vecLineX(intC)) sprintf('%.3f; ',vecP(intC))];
	end
	strTitle = [strTitle '; overall, p=' sprintf('%.3f',dblP_All)];
	
	title(strTitle)
	grid on
	xlabel('Contrast')
	ylabel(strLabelY)
	xlim(vecContrasts(vecWindow))
	ylim([-0.02 0.08])
	legend({'StErr','Miss','StErr','Hit'},'Location','Best')
	
	if boolSavePlots
		intFigCounter = intFigCounter + 1;
		drawnow;
		strFig = sprintf('Meta%d_%s_raw',intFigCounter,strFigTitle);
		export_fig([strFig '.tif']);
		export_fig([strFig '.pdf']);
	end
end

%% trace activity plots
%2 (detect/nodetect) x 2 (pref neurons/non-pref neurons) plot for each contrast
%mean (blue line) and single animals (grey lines)
hTraceAct = zeros(1,6);
cellStrDetect = {'Hit','Miss'};
cellStrPopType = {'Pref','Non-pref'};
for intC=1:6
	hTraceAct(intC) = figure;
	intPlotCounter = 0;
	for intDetect=1:2
		strDetect = cellStrDetect{intDetect};
		for intPrefPop=1:2
			strPopType = cellStrPopType{intPrefPop};
			%get data
			intFrames = size(cellTraceAct{intPrefPop,intDetect,intC,1},2);
			matAct = zeros(intAnimals,intFrames);
			for intAnimal = 1:intAnimals
				matAct(intAnimal,:) = nanmean(cellTraceAct{intPrefPop,intDetect,intC,intAnimal},1);
			end
			
			%plot
			intPlotCounter = intPlotCounter + 1;
			subplot(2,2,intPlotCounter);
			
			%pre-compute variables
			dblSamplingFreq = 25.36;
			vecWindow = round([-3*dblSamplingFreq 5*dblSamplingFreq]);
			vecWindowSelect = vecWindow(1):vecWindow(end);
			intWL = length(vecWindowSelect);
			vecWindowInv = intWL:-1:1;
			vecWindowPlotInv = length(vecWindowSelect):-1:1;
			vecX = [vecWindowSelect vecWindowSelect(vecWindowPlotInv)]/dblSamplingFreq;
			vecLineX = vecWindowSelect/dblSamplingFreq;
			
			
			vecColorLine = [0 0 1];
			vecColorFill = [0.7 0.7 1];
			vecMeanTrace = mean(matAct,1);
			vecSE = std(matAct,[],1)/sqrt(size(matAct,1));
			
			
			vecMinTrace = vecMeanTrace-vecSE;
			vecMaxTrace = vecMeanTrace+vecSE;
			vecY = [vecMinTrace vecMaxTrace(vecWindowInv)];
			
			%plot
			hold on
			
			fill(vecX,vecY,vecColorFill,'EdgeColor','none');
			plot(vecLineX,vecMeanTrace,'-','LineWidth',2,'Color',vecColorLine);
			hold off
			xlim([-3 5])
			ylim([-0.02 0.08])
			title(sprintf('Contrast %.1f%%; %s %s',vecContrasts(intC),strDetect,strPopType))
			grid on
			xlabel('Time after stimulus onset (s)')
			ylabel('Mean pop dF/F over animals')
		end
	end
	if boolSavePlots
		intFigCounter = intFigCounter + 1;
		drawnow;
		strFig = sprintf('Meta%d_traceAct_C%.1f_raw',intFigCounter,vecContrasts(intC));
		export_fig([strFig '.tif']);
		export_fig([strFig '.pdf']);
	end
end
%}

%% meta analyses part 2
%% blob matrices
%2 x 2 heat map plot; detect/no-detect/diff/significance
%mean over animals
%cellMatrices = cat(2,cellMatrices,cellSaveMatrices');
%: blob matrices in [1 x 4 (detect/no-detect/diff/significance)] cell array with [5 (c) x 31 (bins)] matrix
intSwitchZ = 1;
if intSwitchZ == 1
	vecBins = -1:0.1:2;
	vecContrasts = [0 0.5 2 8 32 100];
	vecTickY = [1 11 21 31];
	vecTickLabelY = vecBins(vecTickY);
	vecTickLabelX = vecContrasts;
	strLabelY = sprintf('Z-scored dF/F0 distribution');
else
	vecBins = -0.05:0.01:0.25;
	vecContrasts = [0 0.5 2 8 32 100];
	vecTickY = [6 16 26];%[6 16 26 36];%[1 6 11 16 20];
	vecTickLabelY = vecBins(vecTickY);
	vecTickLabelX = vecContrasts;
	strLabelY = sprintf('Normalized dF/F0 distribution');
end

hFigHeatMaps = figure;
for intType = 1:4
	if intType == 1
		strTitle = sprintf('Hit; z=normalized count; mean over animals (n=%d)',size(cellMatrices,2));
	elseif intType == 2
		strTitle = sprintf('Miss; z=normalized count; mean over animals (n=%d)',size(cellMatrices,2));
	elseif intType == 3
		strTitle = sprintf('Difference; z=normalized count; mean over animals (n=%d)',size(cellMatrices,2));
	elseif intType == 4
		strTitle = sprintf('Significance; z=normalized count; mean over animals (n=%d)',size(cellMatrices,2));
	end
	
	subplot(2,2,intType)
	matAgg = zeros(size(cellMatrices{intType,1},1),size(cellMatrices{intType,1},2),size(cellMatrices,2));
	for intAnimal=1:size(cellMatrices,2)
		matAgg(:,:,intAnimal) = double(cellMatrices{intType,intAnimal});
	end
	cellMatMean{intType} = nanmean(matAgg,3);
	if intType == 4
		intTestType = 3;
		
		if intTestType == 2
			matAgg = zeros(size(cellMatrices{3,1},1),size(cellMatrices{3,1},2),size(cellMatrices,2));
			for intAnimal=1:size(cellMatrices,2)
				matAgg(:,:,intAnimal) = double(cellMatrices{3,intAnimal});
			end
			[matH,matP]=ttest(matAgg,0.5,[],[],3);
			
			[matH_adj crit_p matP_Adj]=fdr_bh(matP,0.05,'pdep');
			
			cellMatMean{intType} = matH;
			
		elseif intTestType == 1
			cellMatMean{intType} = cellMatMean{intType}*size(cellMatrices,2);
			matDiffT = matDiff;
			matDiffT(isnan(matDiffT)) = nanmean(matDiff(:));
			matP = reshape(zscore(matDiffT(:)),size(matDiff));
			cellMatMean{intType} = matP>2;
		elseif intTestType == 3
			matAgg = zeros(size(cellMatrices{3,1},1),size(cellMatrices{3,1},2),size(cellMatrices,2));
			for intAnimal=1:size(cellMatrices,2)
				matAgg(:,:,intAnimal) = double(cellMatrices{3,intAnimal});
			end
			matRealMean = matDiff+0.5;
			
			%shuffle values per animal
			intIters = 1000;
			matShuffledMean = zeros(size(matAgg,1),size(matAgg,2),intIters);
			for intIter=1:1000
				matShuffled = zeros(size(matAgg));
				for intAnimal=1:size(matAgg,3)
					matTemp = matAgg(:,:,intAnimal);
					matShuffled(:,:,intAnimal) = reshape(matTemp(randperm(numel(matTemp))),size(matTemp));
				end
				matShuffledMean(:,:,intIter) = mean(matShuffled,3);
			end
			
			%calc 95% CI
			matCutOff=mean(matShuffledMean,3)+2*std(matShuffledMean,[],3);
			cellMatMean{intType} = matRealMean>matCutOff;
		end
	end
	if intType < 3,imagesc(cellMatMean{intType},[0 0.18]);colormap(hot);colorbar;drawnow;freezeColors;cbfreeze;
	elseif intType == 3
		matA = cellMatMean{1} ./ (cellMatMean{1} + cellMatMean{2});
		%matDiff = buildConvolutedMatrix(matA,1,(1/3)*[1;1;1]);
		matDiff = conv2(matA-0.5,(1/3)*[1;1;1],'same');
		imagesc(matDiff+0.5,[0.35-eps 0.65+eps]);
		colormap(redblue);colorbar;drawnow;freezeColors;cbfreeze;
	else imagesc(cellMatMean{intType});colormap(hot);colorbar;drawnow;freezeColors;cbfreeze;
	end
	
	%plot
	axis xy
	title(strTitle)
	set(gca,'YTick',vecTickY,'YTickLabel',vecTickLabelY)
	set(gca,'XTick',1:length(vecTickLabelX),'XTickLabel',vecTickLabelX)
	ylabel(strLabelY)
	xlabel('Contrast (%)')
end
if boolSavePlots
	intFigCounter = intFigCounter + 1;
	drawnow;
	jFig = get(handle(gcf), 'JavaFrame');
	jFig.setMaximized(true);
	figure(gcf);
	drawnow;
	strFig = sprintf('Meta%d_heatmaps_raw',intFigCounter);
	export_fig([strFig '.tif']);
	export_fig([strFig '.pdf']);
end

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

%% population activity heterogeneity (similarity distance) with removal most/least active active
hHeterogeneity = figure;
intPopulations = size(cellNormActDissim,3);
for intRemType=1:5
	if intRemType == 1
		%no rem
		strRemType = 'No';
		intPlotNr = 1;
	elseif intRemType == 2
		%rem hi
		strRemType = 'High Z';
		intPlotNr = 2;
	elseif intRemType == 3
		%rem lo
		strRemType = 'Low Z';
		intPlotNr = 3;
	elseif intRemType == 4
		%rem hi
		strRemType = 'High A';
		intPlotNr = 5;
	elseif intRemType == 5
		%rem lo
		strRemType = 'Low A';
		intPlotNr = 6;
	end
	matGroup = nan(2,intPopulations);
	matData = nan(2,intPopulations);
	for intPopulation = 1:intPopulations
		matGroup(:,intPopulation) = 1:2;
		matData(2,intPopulation) = mean(cellNormActDissim{1,1,intPopulation,intRemType}); %hit
		matData(1,intPopulation) = mean(cellNormActDissim{1,2,intPopulation,intRemType});%miss
	end
	
	% perform ttests
	[h,dblP] = ttest(matData(1,:),matData(2,:)); %assembly/non-assembly
	
	%normalize
	matData = matData./repmat(mean(matData,1),[2 1]);
	vecMean = nanmean(matData,2);
	vecErr = nanstd(matData,[],2)/sqrt(intPopulations);
	
	%plot
	subplot(2,3,intPlotNr)
	errorbar((1:2)+0.1,vecMean,vecErr,'ob','LineStyle','none');
	hold on
	scatter(matGroup(:)-0.1,matData(:),'rx')
	hold off
	set(gca,'XTick',1:2,'XTickLabel',{'Miss','Hit'})
	ylabel('Normalized heterogeneity')
	title(sprintf('%s rem; T-test p-value=%.3f',strRemType,dblP))
	ylim([0.94 1.06]);
end

if boolSavePlots
	intFigCounter = intFigCounter + 1;
	drawnow;
	jFig = get(handle(gcf), 'JavaFrame');
	jFig.setMaximized(true);
	figure(gcf);
	drawnow;
	strFig = sprintf('Meta%d_activationdissimilarity_raw',intFigCounter);
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
%{
cellSaveCorrITD{1,1,intPopType,intPopulation} = vecDistMiss;
cellSaveCorrITD{1,2,intPopType,intPopulation} = vecCorrMiss;
cellSaveCorrITD{1,3,intPopType,intPopulation} = vecHetMiss;
cellSaveCorrITD{1,4,intPopType,intPopulation} = vecActMiss;

cellSaveCorrITD{2,1,intPopType,intPopulation} = vecDistHit;
cellSaveCorrITD{2,2,intPopType,intPopulation} = vecCorrHit;
cellSaveCorrITD{2,3,intPopType,intPopulation} = vecHetHit;
cellSaveCorrITD{2,4,intPopType,intPopulation} = vecActHit;

cellSaveCorrITD{3,1,intPopType,intPopulation} = vecDistSlow;
cellSaveCorrITD{3,2,intPopType,intPopulation} = vecCorrSlow;
cellSaveCorrITD{3,3,intPopType,intPopulation} = vecHetSlow;
cellSaveCorrITD{3,4,intPopType,intPopulation} = vecActSlow;

cellSaveCorrITD{4,1,intPopType,intPopulation} = vecDistFast;
cellSaveCorrITD{4,2,intPopType,intPopulation} = vecCorrFast;
cellSaveCorrITD{4,3,intPopType,intPopulation} = vecHetFast;
cellSaveCorrITD{4,4,intPopType,intPopulation} = vecActFast;
%}
%create figs
hPopAct = figure;
hPopCorr = figure;

%for whole pop or only preferred pop
matHetMean = nan(4,3,intAnimals);
matActMean = nan(4,3,intAnimals);
matCorrMean = nan(4,3,intAnimals);
matFitsMean = nan(4,3,intAnimals,2);

for intPopType=1:3
	if intPopType == 1
		strPopType = 'AllNeurons';
	elseif intPopType == 2
		strPopType = 'PrefNeurons';
	elseif intPopType == 3
		strPopType = 'NonPrefNeurons';
	end
	
	%get data
	for intAnimal=1:intAnimals
		vecDistMiss = cellCorrITD{1,1,intPopType,intAnimal};
		vecCorrMiss = cellCorrITD{1,2,intPopType,intAnimal};
		vecHetMiss = cellCorrITD{1,3,intPopType,intAnimal};
		vecActMiss = cellCorrITD{1,4,intPopType,intAnimal};
		
		vecDistHit = cellCorrITD{2,1,intPopType,intAnimal};
		vecCorrHit = cellCorrITD{2,2,intPopType,intAnimal};
		vecHetHit = cellCorrITD{2,3,intPopType,intAnimal};
		vecActHit = cellCorrITD{2,4,intPopType,intAnimal};
		
		vecDistSlow = cellCorrITD{3,1,intPopType,intAnimal};
		vecCorrSlow = cellCorrITD{3,2,intPopType,intAnimal};
		vecHetSlow = cellCorrITD{3,3,intPopType,intAnimal};
		vecActSlow = cellCorrITD{3,4,intPopType,intAnimal};
		
		vecDistFast = cellCorrITD{4,1,intPopType,intAnimal};
		vecCorrFast = cellCorrITD{4,2,intPopType,intAnimal};
		vecHetFast = cellCorrITD{4,3,intPopType,intAnimal};
		vecActFast = cellCorrITD{4,4,intPopType,intAnimal};
		
		%inter-trial correlation
		dblMeanCorrHit = mean(vecCorrHit);
		dblMeanCorrMiss = mean(vecCorrMiss);
		dblMeanCorrFast = mean(vecCorrFast);
		dblMeanCorrSlow = mean(vecCorrSlow);
		
		matCorrMean(1,intPopType,intAnimal) = dblMeanCorrHit;
		matCorrMean(2,intPopType,intAnimal) = dblMeanCorrMiss;
		matCorrMean(3,intPopType,intAnimal) = dblMeanCorrSlow;
		matCorrMean(4,intPopType,intAnimal) = dblMeanCorrFast;
		
		%heterogeneity
		matHetMean(1,intPopType,intAnimal) = mean(vecHetHit);
		matHetMean(2,intPopType,intAnimal) = mean(vecHetMiss);
		matHetMean(3,intPopType,intAnimal) = mean(vecHetSlow);
		matHetMean(4,intPopType,intAnimal) = mean(vecHetFast);
		
		%activity
		matActMean(1,intPopType,intAnimal) = mean(vecActHit);
		matActMean(2,intPopType,intAnimal) = mean(vecActMiss);
		matActMean(3,intPopType,intAnimal) = mean(vecActSlow);
		matActMean(4,intPopType,intAnimal) = mean(vecActFast);
		
		%fits
		sStats=regstats(vecCorrHit,vecDistHit,'linear',{'beta'});
		matFitsMean(1,intPopType,intAnimal,1) = sStats.beta(1);
		matFitsMean(1,intPopType,intAnimal,2) = sStats.beta(2);
		
		sStats=regstats(vecCorrMiss,vecDistMiss,'linear',{'beta'});
		matFitsMean(2,intPopType,intAnimal,1) = sStats.beta(1);
		matFitsMean(2,intPopType,intAnimal,2) = sStats.beta(2);
		
		sStats=regstats(vecCorrSlow,vecDistSlow,'linear',{'beta'});
		matFitsMean(3,intPopType,intAnimal,1) = sStats.beta(1);
		matFitsMean(3,intPopType,intAnimal,2) = sStats.beta(2);
		
		sStats=regstats(vecCorrFast,vecDistFast,'linear',{'beta'});
		matFitsMean(4,intPopType,intAnimal,1) = sStats.beta(1);
		matFitsMean(4,intPopType,intAnimal,2) = sStats.beta(2);
	end
	
	%get data
	vecHetHit = matHetMean(1,intPopType,:);
	vecActHit = matActMean(1,intPopType,:);
	
	vecHetMiss = matHetMean(2,intPopType,:);
	vecActMiss = matActMean(2,intPopType,:);
	
	vecHetSlow = matHetMean(3,intPopType,:);
	vecActSlow = matActMean(3,intPopType,:);
	
	vecHetFast = matHetMean(4,intPopType,:);
	vecActFast = matActMean(4,intPopType,:);
	
	
	%plot heterogeneity + activity
	figure(hPopAct);
	
	subplot(3,2,1+(intPopType-1)*2)
	
	matHetMean(1,intPopType,intAnimal)
	
	[h,pHM] = ttest2(vecHetHit,vecHetMiss);
	[h,pMF] = ttest2(vecHetMiss,vecHetFast);
	[h,pMS] = ttest2(vecHetMiss,vecHetSlow);
	[h,pFS] = ttest2(vecHetFast,vecHetSlow);
	
	errorbar(1:4,[mean(vecHetMiss) mean(vecHetHit) mean(vecHetFast) mean(vecHetSlow)],...
		[std(vecHetMiss)/sqrt(length(vecHetMiss)) std(vecHetHit)/sqrt(length(vecHetHit)) std(vecHetFast)/sqrt(length(vecHetFast)) std(vecHetSlow)/sqrt(length(vecHetSlow))],...
		'Linestyle','none','Marker','x','Color','k');
	set(gca,'XTick',1:4,'XTickLabel',{'Miss','Hit-All','Hit-Fast','Hit-Slow'})
	title(sprintf('Mean heterogeneity of %s population\nT-tests; Hit-miss,p=%.3f; Miss-Fast,p=%.3f; Miss-Slow,p=%.3f; Fast-Slow,p=%.3f',strPopType,pHM,pMF,pMS,pFS))
	ylabel('Heterogeneity')
	
	subplot(3,2,2+(intPopType-1)*2)
	
	[h,pHM] = ttest2(vecActHit,vecActMiss);
	[h,pMF] = ttest2(vecActMiss,vecActFast);
	[h,pMS] = ttest2(vecActMiss,vecActSlow);
	[h,pFS] = ttest2(vecActFast,vecActSlow);
	
	errorbar(1:4,[mean(vecActMiss) mean(vecActHit) mean(vecActFast) mean(vecActSlow)],...
		[std(vecActMiss)/sqrt(length(vecActMiss)) std(vecActHit)/sqrt(length(vecActHit)) std(vecActFast)/sqrt(length(vecActFast)) std(vecActSlow)/sqrt(length(vecActSlow))],...
		'Linestyle','none','Marker','x','Color','k');
	set(gca,'XTick',1:4,'XTickLabel',{'Miss','Hit-All','Hit-Fast','Hit-Slow'})
	title(sprintf('Mean dF/F0 activity of %s population\nT-tests; Hit-miss,p=%.3f; Miss-Fast,p=%.3f; Miss-Slow,p=%.3f; Fast-Slow,p=%.3f',strPopType,pHM,pMF,pMS,pFS))
	ylabel('Mean dF/F0 activity')
	
	
	%% plot inter-trial correlation vs distance; + relation dF/F0 vs heterogeneity
	figure(hPopCorr);
	
	vecCorrHit = matCorrMean(1,intPopType,:);
	vecCorrMiss = matCorrMean(2,intPopType,:);
	vecCorrFast = matCorrMean(3,intPopType,:);
	vecCorrSlow = matCorrMean(4,intPopType,:);
	
	subplot(3,2,2+(intPopType-1)*2)
	dblMeanCorrHit = mean(vecCorrHit);
	dblMeanCorrMiss = mean(vecCorrMiss);
	dblMeanCorrFast = mean(vecCorrFast);
	dblMeanCorrSlow = mean(vecCorrSlow);
	dblErrCorrHit = std(vecCorrHit)/sqrt(length(vecCorrHit));
	dblErrCorrMiss = mean(vecCorrMiss)/sqrt(length(vecCorrMiss));
	dblErrCorrFast = mean(vecCorrFast)/sqrt(length(vecCorrFast));
	dblErrCorrSlow = mean(vecCorrSlow)/sqrt(length(vecCorrSlow));
	
	
	[h,pM] = ttest(vecCorrMiss);
	[h,pH] = ttest(vecCorrHit);
	[h,pF] = ttest(vecCorrFast);
	[h,pS] = ttest(vecCorrSlow);
	
	errorbar(1:4,[dblMeanCorrMiss dblMeanCorrHit dblMeanCorrFast dblMeanCorrSlow],[dblErrCorrMiss dblErrCorrHit dblErrCorrFast dblErrCorrSlow],'Linestyle','none','Marker','x','Color','k');
	hold on;
	plot(get(gca,'XLim'),[0 0],'k--');
	hold off;
	set(gca,'XTick',1:4,'XTickLabel',{'Miss','Hit-All','Hit-Fast','Hit-Slow'})
	title(sprintf('Mean inter-trial correlation of %s population activity\nT-tests; Miss,p=%.3f; Hit,p=%.3f; Fast,p=%.3f; Slow,p=%.3f',strPopType,pM,pH,pF,pS))
	ylabel('Pearson correlation of population activity')
	ylim([-0.01 0.05])
	
	%{
	%plot ITD
	subplot(3,2,1+(intPopType-1)*2)
	
	vecHitB2 = matFitsMean(1,intPopType,:,2);
	vecMissB2 = matFitsMean(2,intPopType,:,2);
	vecFastB2 = matFitsMean(3,intPopType,:,2);
	vecSlowB2 = matFitsMean(4,intPopType,:,2);
	
	[h,pH] = ttest(vecHitB2);
	[h,pM] = ttest(vecMissB2);
	[h,pF] = ttest(vecFastB2);
	[h,pS] = ttest(vecSlowB2);
	
	errorbar(1:4,[mean(vecMissB2) mean(vecHitB2) mean(vecFastB2) mean(vecSlowB2)],[std(vecMissB2) std(vecHitB2) std(vecFastB2) std(vecSlowB2)]/sqrt(intAnimals),'Linestyle','none','Marker','x','Color','k');
	set(gca,'XTick',1:4,'XTickLabel',{'Miss','Hit-All','Hit-Fast','Hit-Slow'})
	title(sprintf('Mean slopes of inter-trial correlation drop-off of %s population activity\nT-tests; Miss,p=%.3f; Hit,p=%.3f; Fast,p=%.3f; Slow,p=%.3f',strPopType,pH,pM,pF,pS))
	ylabel('Slopes of inter-trial correlation')
	%}
end


if boolSavePlots
	intFigCounter = intFigCounter + 1;
	figure(hPopCorr);
	drawnow;
	jFig = get(handle(hPopCorr), 'JavaFrame');
	jFig.setMaximized(true);
	figure(hPopCorr);
	drawnow;
	strFig = sprintf('Meta%d_intertrialcorr_%s_pop%d_raw',intFigCounter);
	export_fig([strFig '.tif']);
	export_fig([strFig '.pdf']);
	
	intFigCounter = intFigCounter + 1;
	figure(hPopAct);
	drawnow;
	jFig = get(handle(hPopAct), 'JavaFrame');
	jFig.setMaximized(true);
	figure(hPopAct);
	drawnow;
	strFig = sprintf('Meta%d_HetAct__raw',intFigCounter);
	export_fig([strFig '.tif']);
	export_fig([strFig '.pdf']);
end

%% decoding plot (hit/miss HCAR/random)
matMeanCorrectDecoded = nan(2,6); %rand/miss rand/hit HCAR/miss HCAR/hit
matSdCorrectDecoded = nan(2,6);

%calc mean
matMeanCorrectDecoded(2,:) = nanmean(matMixBootstrappedDecodingOutput(:,:,1),1); %HCAR hit
matMeanCorrectDecoded(1,:) = nanmean(matMixBootstrappedDecodingOutput(:,:,2),1); %HCAR miss

%calc sd
matSdCorrectDecoded(2,:) = nanstd(matMixBootstrappedDecodingOutput(:,:,1),[],1); %HCAR hit
matSdCorrectDecoded(1,:) = nanstd(matMixBootstrappedDecodingOutput(:,:,2),[],1); %HCAR miss

%calc over contrasts
matMeanOverContrasts(2,:) = nanmean(matMixBootstrappedDecodingOutput(:,:,1),2);
matMeanOverContrasts(1,:) = nanmean(matMixBootstrappedDecodingOutput(:,:,2),2);

%set vars
vecLineX = 1:6;
vecWindowInv = 6:-1:1;
vecX = [vecLineX vecWindowInv];

%plot
intN = size(matMixBootstrappedDecodingOutput,1);
h=figure;
for intPlotType=1:2
	
	vecMeanTrace = matMeanCorrectDecoded(intPlotType,:);
	vecSE = matSdCorrectDecoded(intPlotType,:)./sqrt(intN);
	
	if intPlotType == 1
		vecColorFill = [1.0 0.7 0.7];
		vecColorLine = [1 0 0];
	elseif intPlotType == 2
		vecColorFill = [0.7 1.0 0.7];
		vecColorLine = [0 1 0];
	end
	
	vecMinTrace = vecMeanTrace-vecSE;
	vecMaxTrace = vecMeanTrace+vecSE;
	vecY = [vecMinTrace vecMaxTrace(vecWindowInv)];
	
	%plot
	hold on
	fill(vecX,vecY,vecColorFill,'EdgeColor','none');
	plot(vecLineX,vecMeanTrace,'-','LineWidth',2,'Color',vecColorLine);
	
	%scatter mean
	dblX = max(vecX)+0.5+0.2*intPlotType;
	dblY = mean(matMeanOverContrasts(intPlotType,:),2);
	scatter(dblX,dblY,[],vecColorLine);
	dblErr = std(matMeanOverContrasts(intPlotType,:),[],2)./sqrt(intN);
	errorbar(dblX,dblY,dblErr,'Color',vecColorLine);
	
	hold off
	%end
end
legend(gca,'Miss','Hit','Location', 'Best')
set(gca,'XTick',[1 2 3 4 5 6],'XTickLabel',[0 0.5 2 8 32 100]);
xlabel('Stimulus Contrast (%)')
ylabel('Decoding performance (proportion correct)')
ylim([0 0.3])
xlim([0.5 7.5])

hold on
plot(get(gca,'XLim'),[1/21 1/21],'k--')
hold off

[h,dblHitMiss]=ttest(matMeanOverContrasts(1,:),matMeanOverContrasts(2,:));

%ttests
[h,p0]=ttest(matMixBootstrappedDecodingOutput(:,1,1),matMixBootstrappedDecodingOutput(:,1,2));
[h,p05]=ttest(matMixBootstrappedDecodingOutput(:,2,1),matMixBootstrappedDecodingOutput(:,2,2));
[h,p2]=ttest(matMixBootstrappedDecodingOutput(:,3,1),matMixBootstrappedDecodingOutput(:,3,2));
[h,p8]=ttest(matMixBootstrappedDecodingOutput(:,4,1),matMixBootstrappedDecodingOutput(:,4,2));
[h,p32]=ttest(matMixBootstrappedDecodingOutput(:,5,1),matMixBootstrappedDecodingOutput(:,5,2));
[h,p100]=ttest(matMixBootstrappedDecodingOutput(:,6,1),matMixBootstrappedDecodingOutput(:,6,2));

title(sprintf('(1)Red=Miss;(2)Green=Hit;ttest p-vals;hit-miss: %.3f;Hit/miss p-vals per contrast; 0=%.3f; 0.5=%.3f; 2=%.3f; 8=%.3f; 32=%.3f; 100=%.3f'...
	,dblHitMiss,p0,p05,p2,p8,p32,p100))
strFigTitle = 'decoding_performance';
if boolSavePlots
	intFigCounter = intFigCounter + 1;
	drawnow;
	strFig = sprintf('Meta%d_%s_raw',intFigCounter,strFigTitle);
	export_fig([strFig '.tif']);
	export_fig([strFig '.pdf']);
end

%% heterogeneity over contrasts
%normalize heterogeneity per animal
matContHeteroNorm = zeros(size(matContHetero));
for intAnimal=1:size(matContHetero,3)
	matTempHet = matContHetero(:,:,intAnimal);
	matContHeteroNorm(:,:,intAnimal) = matTempHet ./ mean(matTempHet(:));
end
matHitHet = matContHeteroNorm(:,1,:);
matMissHet = matContHeteroNorm(:,2,:);
vecHitY = mean(matHitHet,3)';
vecHitE = std(matHitHet,[],3)'/sqrt(size(matContHeteroNorm,3));
vecMissY = mean(matMissHet,3)';
vecMissE = std(matMissHet,[],3)'/sqrt(size(matContHeteroNorm,3));

%get sig
vecP = zeros(1,size(matContHeteroNorm,1));
for intC=1:size(matContHeteroNorm,1)
	[h,p,ci] = ttest(matHitHet(intC,1,:),matMissHet(intC,1,:));
	
	%put in vector
	vecP(intC) = p;
end

%overall t-test
matHetHits = matHitHet(1:4,1,:);
matHetMisses = matMissHet(1:4,1,:);
[h,dblP_All,ci] = ttest(matHetHits(:),matHetMisses(:));
	
%pre-compute variables
vecContrasts = [0.5 2 8 32 100];
vecWindow = [1 5];
vecWindowSelect = vecWindow(1):vecWindow(end);
intWL = length(vecWindowSelect);
vecWindowInv = intWL:-1:1;
vecWindowPlotInv = vecWindow(end):-1:vecWindow(1);
vecX = [vecContrasts(vecWindowSelect) vecContrasts(vecWindowPlotInv)];
vecLineX = vecContrasts(vecWindowSelect);

%get data
hHetCon=figure;
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
set(gca,'XTick',vecLineX,'XTickLabel',vecLineX)
title(sprintf('Pop het during stim; p-vals: 0.5: %.3f; 2: %.3f; 8: %.3f; 32: %.3f; 100: %.3f; Overall,p=%.3f',[vecP dblP_All]))
grid on
xlabel('Contrast')
ylabel('Normalized population response heterogeneity')
xlim(vecContrasts(vecWindow))
%ylim([-0.01 0.06])
legend({'SEM','Miss','SEM','Hit'},'Location','Best')
drawnow;
strFigTitle = 'hetero_over_contrasts';
if boolSavePlots
	intFigCounter = intFigCounter + 1;
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

%% decoding vs heterogeneity
%decoding vs contrast; decoding vs RT; decoding vs hetero-pre; decoding vs
%hetero-dur; decoding vs act-pre; decoding vs act-dur
hDecode2 = figure;

%contrast
subplot(2,2,1);
vecContrasts = unique(cellHD_Contrast{1});
vecContrastsPlot = [0.0025; vecContrasts(2:end)];
intContrasts = length(vecContrasts);
matDecCon = nan(intAnimals,intContrasts);
vecDecConP = nan(intContrasts,1);
strTitle = sprintf('\n');
for intContrast=1:intContrasts
	for intAnimal=1:intAnimals
		indContrast = cellHD_Contrast{intAnimal}(:) == vecContrasts(intContrast);
		matDecCon(intAnimal,intContrast) = mean(cellHD_DecodingAccuracy{intAnimal}(indContrast));
	end
	[h,vecDecConP(intContrast)] = ttest(matDecCon(:,intContrast),0.25);
end
vecDecConMean = mean(matDecCon(:,:,1),1);
vecDecConSD = std(matDecCon(:,:,1),[],1);

[h crit_p vecDecConP_Adj]=fdr_bh(vecDecConP,0.05,'pdep');


for intContrast=1:intContrasts
	strTitle = [strTitle sprintf(';%.1f:p=%.3f',vecContrasts(intContrast)*100,vecDecConP_Adj(intContrast))];
end

errorbar(vecContrastsPlot,vecDecConMean,vecDecConSD/sqrt(intAnimals),'r');
hold on
plot([vecContrastsPlot(1) vecContrastsPlot(end)],[0.25 0.25],'--k')
hold off
title(['Cross-validated orientation decoding; adjusted for MCs using FDR correction' strTitle])
xlabel('Stimulus Contrast')
ylabel('Orientation decoding performance')
set(gca,'XScale','log','YScale','linear')

%{
%split for hit/miss
subplot(2,3,4);

vecContrasts = unique(cellHD_Contrast{1});
vecContrastsPlot = [0.0025; vecContrasts(2:end)];
intContrasts = length(vecContrasts);
matDecCon_M = nan(intAnimals,intContrasts);
vecDecConP_M = nan(intContrasts,1);
matDecCon_H = nan(intAnimals,intContrasts);
vecDecConP_H = nan(intContrasts,1);
strTitle = sprintf('\n');
for intContrast=1:intContrasts
	for intAnimal=1:intAnimals
		indHit = ~isnan(cellHD_ReactionTime{intAnimal}(:));
		indMiss = isnan(cellHD_ReactionTime{intAnimal}(:));
		indContrast = cellHD_Contrast{intAnimal}(:) == vecContrasts(intContrast);
		matDecCon_M(intAnimal,intContrast) = mean(cellHD_DecodingAccuracy{intAnimal}(indContrast&indMiss));
		
		matDecCon_H(intAnimal,intContrast) = mean(cellHD_DecodingAccuracy{intAnimal}(indContrast&indHit));
	end
	[h,vecDecConP_M(intContrast)] = ttest(matDecCon_M(:,intContrast),0.25);
	[h,vecDecConP_H(intContrast)] = ttest(matDecCon_H(:,intContrast),0.25);
end
vecDecConMean_M = nanmean(matDecCon_M(:,:,1),1);
vecDecConSD_M = nanstd(matDecCon_M(:,:,1),[],1);

vecDecConMean_H = nanmean(matDecCon_H(:,:,1),1);
vecDecConSD_H = nanstd(matDecCon_H(:,:,1),[],1);

[h crit_p vecDecConP_M_Adj]=fdr_bh(vecDecConP_M,0.05,'pdep');
[h crit_p vecDecConP_H_Adj]=fdr_bh(vecDecConP_H,0.05,'pdep');

for intContrast=1:intContrasts
	strTitle = [strTitle sprintf(';%.1f:M-p=%.3f;H-p=%.3f',vecContrasts(intContrast)*100,vecDecConP_M_Adj(intContrast),vecDecConP_H_Adj(intContrast))];
end

errorbar(vecContrastsPlot,vecDecConMean_M,vecDecConSD_M/sqrt(intAnimals),'r');
hold on
errorbar(vecContrastsPlot,vecDecConMean_H,vecDecConSD_H/sqrt(intAnimals),'g');
plot([vecContrastsPlot(1) vecContrastsPlot(end)],[0.25 0.25],'--k')
hold off
title(['Cross-validated orientation decoding; adjusted for MCs using FDR correction' strTitle])
xlabel('Stimulus Contrast')
ylabel('Orientation decoding performance')
set(gca,'XScale','log','YScale','linear')
%}
%{
%split for hit/miss
subplot(2,3,4);

vecContrasts = unique(cellHD_Contrast2{1,1});
vecContrastsPlot = [0.0025; vecContrasts(2:end)];
intContrasts = length(vecContrasts);
matDecCon_M = nan(intAnimals,intContrasts);
vecDecConP_M = nan(intContrasts,1);
matDecCon_H = nan(intAnimals,intContrasts);
vecDecConP_H = nan(intContrasts,1);
strTitle = sprintf('\n');
for intContrast=1:intContrasts
	for intAnimal=1:intAnimals
		indContrast = cellHD_Contrast2{1,intAnimal}(:) == vecContrasts(intContrast);
		matDecCon_M(intAnimal,intContrast) = mean(cellHD_DecodingAccuracy2{1,intAnimal}(indContrast));
		
		indContrast = cellHD_Contrast2{2,intAnimal}(:) == vecContrasts(intContrast);
		matDecCon_H(intAnimal,intContrast) = mean(cellHD_DecodingAccuracy2{2,intAnimal}(indContrast));
	end
	[h,vecDecConP_M(intContrast)] = ttest(matDecCon_M(:,intContrast),0.25);
	[h,vecDecConP_H(intContrast)] = ttest(matDecCon_H(:,intContrast),0.25);
end
vecDecConMean_M = mean(matDecCon_M(:,:,1),1);
vecDecConSD_M = std(matDecCon_M(:,:,1),[],1);

vecDecConMean_H = mean(matDecCon_H(:,:,1),1);
vecDecConSD_H = std(matDecCon_H(:,:,1),[],1);

[h crit_p vecDecConP_M_Adj]=fdr_bh(vecDecConP_M,0.05,'pdep');
[h crit_p vecDecConP_H_Adj]=fdr_bh(vecDecConP_H,0.05,'pdep');

for intContrast=1:intContrasts
	strTitle = [strTitle sprintf(';%.1f:M-p=%.3f;H-p=%.3f',vecContrasts(intContrast)*100,vecDecConP_M_Adj(intContrast),vecDecConP_H_Adj(intContrast))];
end

errorbar(vecContrastsPlot,vecDecConMean_M,vecDecConSD_M/sqrt(intAnimals),'r');
hold on
errorbar(vecContrastsPlot,vecDecConMean_H,vecDecConSD_H/sqrt(intAnimals),'g');
plot([vecContrastsPlot(1) vecContrastsPlot(end)],[0.25 0.25],'--k')
hold off
title(['Cross-validated orientation decoding; adjusted for MCs using FDR correction' strTitle])
xlabel('Stimulus Contrast')
ylabel('Orientation decoding performance')
set(gca,'XScale','log','YScale','linear')
%}

%take highest/lowest 50% act/het per contrast
matDecAcc_ActHigh = zeros(intContrasts,intAnimals);
matDecAcc_ActLow = zeros(intContrasts,intAnimals);
matDecAcc_HetHigh = zeros(intContrasts,intAnimals);
matDecAcc_HetLow = zeros(intContrasts,intAnimals);

%remove catch trials & normalize activity
dblFrac = 1/2;
for intAnimal=1:intAnimals
	for intC=1:intContrasts
		indContrast = cellHD_Contrast{intAnimal} == vecContrasts(intC);
		
		vecHet = cellHD_HeterogeneityDuring{intAnimal}(indContrast);
		vecAct = cellHD_ActivityDuring{intAnimal}(indContrast);
		vecDec = cellHD_DecodingAccuracy{intAnimal}(indContrast);
		
		%sort het
		intVals = floor(length(vecHet)*dblFrac);
		[vecSort,vecIndex] = sort(vecHet,'ascend');
		dblDecLowHet = mean(vecDec(vecIndex(1:intVals)));
		dblDecHighHet = mean(vecDec(vecIndex((end-intVals+1):end)));
		
		%sort act
		intVals = floor(length(vecAct)*dblFrac);
		[vecSort,vecIndex] = sort(vecAct,'ascend');
		dblDecLowAct = mean(vecDec(vecIndex(1:intVals)));
		dblDecHighAct = mean(vecDec(vecIndex((end-intVals+1):end)));
		
		%output
		matDecAcc_ActHigh(intC,intAnimal) = dblDecHighAct;
		matDecAcc_ActLow(intC,intAnimal) = dblDecLowAct;
		matDecAcc_HetHigh(intC,intAnimal) = dblDecHighHet;
		matDecAcc_HetLow(intC,intAnimal) = dblDecLowHet;
	end
end

%diff
matDecAccInc_Act = (matDecAcc_ActHigh - matDecAcc_ActLow) ./ (matDecAcc_ActHigh + matDecAcc_ActLow);
matDecAccInc_Het = (matDecAcc_HetHigh - matDecAcc_HetLow) ./ (matDecAcc_HetHigh + matDecAcc_HetLow);

%remove probe trials
matDecAccInc_Het = matDecAccInc_Het(2:end,:);
matDecAccInc_Act = matDecAccInc_Act(2:end,:);

%t-test over animals (first mean over contrasts per animal and convert to %)
vecDecHet = mean(matDecAccInc_Het,1)*100;
vecDecAct = mean(matDecAccInc_Act,1)*100;
[h,pH]=ttest(vecDecHet);
[h,pA]=ttest(vecDecAct);
[h,pHA]=ttest(vecDecHet,vecDecAct);

%dur-hetero
subplot(2,2,3);
plot([0 3],[0 0],'k--')
hold on
scatter(ones(size(vecDecHet)),vecDecHet,'x','MarkerEdgeColor',[0.5 0.5 0.5])
errorbar(2,mean(vecDecHet),std(vecDecHet)/sqrt(intAnimals),'Linestyle','none','Marker','x','Color','k');
hold off
title(sprintf('Heterogeneity during stim;T-test,p=%.3f',pH));
ylabel('Decoding accuracy increase for high heterogeneity trials (%)')
xlim([0 3])
ylim([-20 20])

%dur act
subplot(2,2,4);
plot([0 3],[0 0],'k--')
hold on
scatter(ones(size(vecDecAct)),vecDecAct,'x','MarkerEdgeColor',[0.5 0.5 0.5])
errorbar(2,mean(vecDecAct),std(vecDecAct)/sqrt(intAnimals),'Linestyle','none','Marker','x','Color','b');
hold off
title(sprintf('Activity during stim;T-test,p=%.3f',pA));
ylabel('Accuracy accuracy increase for high dF/F0 trials (%)')
xlim([0 3])
ylim([-20 20])

strFigTitle = 'decoding_dependence';
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

%% across-repetition correlation
hRepCorr = figure;
vecSlopes = nan(1,intAnimals);
vecIntercepts = nan(1,intAnimals);
dblMin = inf;
dblMax = -inf;
hold on
cMap=colormap(lines(intAnimals));
for intAnimal=1:intAnimals
	vecCorr = cellAcrossRepCorr{1,intAnimal};
	vecTime = cellAcrossRepCorr{3,intAnimal};
	vecBeta = cellAcrossRepCorr{4,intAnimal};
	
	%scatter raw points
	scatter(vecTime/60,vecCorr,[],cMap(intAnimal,:),'x');
	
	vecSlopes(intAnimal) = vecBeta(2);
	vecIntercepts(intAnimal) = vecBeta(1);
	vecX = [min(vecTime/60) max(vecTime/60)];
	vecY = polyval([vecBeta(2) vecBeta(1)],vecX);
	plot(vecX,vecY,'Color',cMap(intAnimal,:));
	dblMin = min([dblMin vecTime/60]);
	dblMax = max([dblMax vecTime/60]);
end

vecX = [dblMin dblMax];
vecBeta = [mean(vecIntercepts) mean(vecSlopes)];
vecY = polyval([vecBeta(2) vecBeta(1)],vecX);
plot(vecX,vecY,'k-','LineWidth',2);
vecLimX = get(gca,'XLim');
xlim([0 vecLimX(2)]);
plot([0 vecLimX(2)],[0 0],'k--');
hold off

[h,p] = ttest(vecSlopes);

ylabel('Across-repetition correlation of correlation structure')
xlabel('Time between stimulus repetitions (minutes)')

title(sprintf('Across-repetition correlation vs distance in time; mean slope=%.3f; p=%.3f; ',mean(vecSlopes),p));

strFigTitle = 'AcrossRepetitionCorrelation';
if boolSavePlots
	intFigCounter = intFigCounter + 1;
	figure(hRepCorr);
	drawnow;
	strFig = sprintf('Meta%d_%s_raw',intFigCounter,strFigTitle);
	export_fig([strFig '.tif']);
	export_fig([strFig '.pdf']);
end
%% clean up
cd(strOldDir);