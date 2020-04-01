%% analyze input strength dependency for dimensionality of pop responses

%% initialize
clearvars;
boolLoad = false;
if ~exist('matModelResp','var')
	clearvars -except bool*;%close all;
	intLoadSim = 13; %% SET WHICH SIM TO LOAD HERE; 41 vs 42; 2 vs. 5
	boolLoad = true;
end
if intLoadSim == 11 && boolLoad
	strSimulation = 'xAreaDistributed_OriFull_2017-06-22';
	%strSimulation = 'xAreaDistributed_ContFull_2017-06-13';
	%strSimulation = 'xAreaDistributed_494423_2017-06-13'; %contrast
elseif intLoadSim == 12 && boolLoad
	strSimulation = 'xAreaDistributed_OriInhOnly_2017-06-23'; %contrast
elseif intLoadSim == 13 && boolLoad
	strSimulation = 'xAreaDistributed_OriExcOnly_2017-06-23'; %contrast
elseif intLoadSim == 14 && boolLoad
	strSimulation = 'xAreaDistributed_OriNone_2017-06-22'; %contrast
	
	
	
elseif intLoadSim == 21 && boolLoad
	strSimulation = 'xAreaDistributed_ContFull_2017-05-04'; %contrast
	
elseif intLoadSim == 31 && boolLoad
	strSimulation = 'xAreaDistributed_LumFull_2017-05-03'; %luminance

end

%% RUN: #header
strFigDir = 'D:\Data\Results\InfoDecoding\';
if boolLoad
	runModelHeader;
end
strTag = [strStimType strSimulation(end-9:end)];

%% set parameters
boolSaveFigs = false;
boolSaveData = false;
dblCutOff = 0.90;
dblLambdaInfo = 1;
dblLambdaPred = 0;
dblNoiseLevel = 0;
boolAnalI = true;
boolAnalPred = true;
boolAnalPopM = true;
intMaxDimAnal = 30;
%intWithinArea = 1; %1, within V1; 2, within V2; 3, across V1-V2;
intIters = 1;

%% calculate which parameter to use
[d,intUseParam]=max(range(matStimTypeCombos,2));
if intUseParam == 1,vecParamStimType = vecStimTypeOris;strParam='Ori';
elseif intUseParam == 2,vecParamStimType = vecStimTypeSFs;strParam='SF';
elseif intUseParam == 3,vecParamStimType = vecStimTypeTFs;strParam='TF';
elseif intUseParam == 4,vecParamStimType = vecStimTypeContrasts;strParam='Contrast';
elseif intUseParam == 5,vecParamStimType = vecStimTypeLuminance;strParam='Luminance';
end
vecParamIdx = matStimTypeCombos(intUseParam,:);
vecParamVals = unique(vecParamStimType);
intParamNum=numel(vecParamVals);
mapC = redbluepurple(intParamNum);

%% add noise to improve numerical stability
matModelResp = double(matModelResp);
intNeurons = size(matModelResp,1);
intTrials = size(matModelResp,2);
vecNeuronSD = xstd(matModelResp,2);
matNoise = normrnd(zeros(intNeurons,intTrials),repmat(vecNeuronSD*dblNoiseLevel,[1 intTrials]));
matModelRespP = matModelResp + matNoise;
%matModelRespP = zscore(matModelRespP,[],1);
vecGroupSizes = 2.^(1:8);
intGroupSizes = numel(vecGroupSizes);
%dblStep = intCellsV1/intPopSize;
%vecNeurons = round((1:dblStep:intCellsV1) + randi(dblStep) - 1);

%% build pairwise stim type combinations & pre-allocate output matrices
matCombos = nchoosek(vecParamIdx,2);
matCombosReduced = nchoosek(vecParamIdx(1:floor(numel(vecParamIdx)/2)),2);
matI_Dir = nan(intParamNum,intParamNum,intIters,intGroupSizes);
matI_Dir_shuff = nan(intParamNum,intParamNum,intIters,intGroupSizes);
matA_ML = nan(intParamNum,intParamNum,intIters,intGroupSizes);
matA_ML_shuff = nan(intParamNum,intParamNum,intIters,intGroupSizes);

cellArea = {'V1','V2'};
for intArea=1:2
	if intArea == 1
		matData = matModelRespP(1:intCellsV1,:);
	else
		matData = matModelRespP((intCellsV1+1):end,:);
	end
	% prep figure
	h=figure;
	
	% fisher anal
	for intCombo=1:size(matCombosReduced,1)
		vecUseStimTypes = matCombosReduced(intCombo,:);
		vecUseStimTypes(1) = vecUseStimTypes(1)*2;
		vecUseStimTypes(2) = vecUseStimTypes(2)*2-1;
		
		%build parameter structure
		sParamsAnal = struct;
		sParamsAnal.vecUseStimTypes = vecUseStimTypes;
		sParamsAnal.dblLambda = dblLambdaInfo;
		sParamsAnal.intIters = intIters;
		sParamsAnal.vecGroupSizes = vecGroupSizes;
		sParamsAnal.boolVerbose = true;
		sParamsAnal.dblDiffTheta = abs(rad2ang(circ_dist(ang2rad(vecTrialOris(vecUseStimTypes)*2)))/2);
		sParamsAnal.boolDirectI = true;
		sParamsAnal.boolDecodeLR = false;
		sParamsAnal.boolDecodeML = true;
		sParamsAnal.boolDecodeMD = false;
		
		%do analysis
		sOut = doFisherAnal(matData,vecTrialStimType,sParamsAnal);
		
		%% assign to output
		matI_DirTemp = nanmean(sOut.matI_Direct_bc,3); %[iter x group_size x fold]
		matI_Dir(vecUseStimTypes(1),vecUseStimTypes(2),:,:) = matI_DirTemp;
		matI_Dir(vecUseStimTypes(2),vecUseStimTypes(1),:,:) = matI_DirTemp;
		matI_DirTemp_shuff= nanmean(sOut.matI_Direct_bc_shuff,3);
		matI_Dir_shuff(vecUseStimTypes(1),vecUseStimTypes(2),:,:) = matI_DirTemp_shuff;
		matI_Dir_shuff(vecUseStimTypes(2),vecUseStimTypes(1),:,:) = matI_DirTemp_shuff;
		
		matATemp_ML = nanmean(sOut.matA_ML,3);
		matA_ML(vecUseStimTypes(1),vecUseStimTypes(2),:,:) = matATemp_ML;
		matA_ML(vecUseStimTypes(2),vecUseStimTypes(1),:,:) = matATemp_ML;
		matATemp_ML_shuff= nanmean(sOut.matA_ML_shuff,3);
		matA_ML_shuff(vecUseStimTypes(1),vecUseStimTypes(2),:,:) = matATemp_ML_shuff;
		matA_ML_shuff(vecUseStimTypes(2),vecUseStimTypes(1),:,:) = matATemp_ML_shuff;
		
		%plot
		doPlotInfo(vecGroupSizes,matI_Dir,matI_Dir_shuff,matA_ML,matA_ML_shuff);
		
	end
	
	%% save figure
	cellSplit = strsplit(strSimulation,'_');
	strFile = sprintf('%s%s_InfoDecoding_%s',cellSplit{2},cellArea{intArea},cellSplit{3});
	drawnow;
	jFig = get(handle(gcf), 'JavaFrame');
	jFig.setMaximized(true);
	figure(gcf);
	drawnow;
	export_fig([strFigDir strFile '.tif']);
	export_fig([strFigDir strFile '.pdf']);
	fprintf('Saved figure to %s [%s]\n',strFile,getTime);
end
return

%% plot covar
for intStimType=vecParamIdx
	indTrials = vecTrialStimType==intStimType;
	matCorr = corr(matData(:,indTrials)');
end

%% rotate vectors in matrix to align 0
matMeanLRICent = nan(size(matMeanLRI));
for intPos=1:intParamNum
	matMeanLRICent(intPos,:) = matMeanLRI(intPos,mod(((1:intParamNum)-1)-(intParamNum/2)+intPos,intParamNum)+1);
end
