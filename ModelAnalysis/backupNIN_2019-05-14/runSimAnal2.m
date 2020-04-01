%% analyze input strength dependency for dimensionality of pop responses

%% initialize
clearvars;
boolLoad = false;
if ~exist('matModelResp','var')
	clearvars -except bool*;%close all;
	intLoadSim = 12; %% SET WHICH SIM TO LOAD HERE; 41 vs 42; 2 vs. 5
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
intIters = 3;
intPopSize = 3;

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

%dblStep = intCellsV1/intPopSize;
%vecNeurons = round((1:dblStep:intCellsV1) + randi(dblStep) - 1);

%% build pairwise stim type combinations & pre-allocate output matrices
matCombos = nchoosek(vecParamIdx,2);
matI_LR = nan(intParamNum,intParamNum,intIters);
matI_LR_shuff = nan(intParamNum,intParamNum,intIters);
matI_Dir = nan(intParamNum,intParamNum,intIters);
matI_Dir_shuff = nan(intParamNum,intParamNum,intIters);
matA_LR = nan(intParamNum,intParamNum,intIters);
matA_LR_shuff = nan(intParamNum,intParamNum,intIters);
matA_ML = nan(intParamNum,intParamNum,intIters);
matA_ML_shuff = nan(intParamNum,intParamNum,intIters);
matA_MD = nan(intParamNum,intParamNum,intIters);
matA_MD_shuff = nan(intParamNum,intParamNum,intIters);
matData = matModelRespP(1:intCellsV1,:);

%% prep figure
h=figure;

%% fisher anal
for intCombo=1:size(matCombos,1)
	vecUseStimTypes = matCombos(intCombo,:);
	
	%build parameter structure
	sParamsAnal = struct;
	sParamsAnal.vecUseStimTypes = vecUseStimTypes;
	sParamsAnal.dblLambda = dblLambdaInfo;
	sParamsAnal.boolDirectI = true;
	sParamsAnal.intIters = intIters;
	sParamsAnal.vecGroupSizes = intPopSize;
	sParamsAnal.boolVerbose = true;
	sParamsAnal.dblDiffTheta = abs(rad2ang(circ_dist(ang2rad(vecTrialOris(vecUseStimTypes)*2)))/2);
	sParamsAnal.boolDecodeLR = true;
	sParamsAnal.boolDecodeML = true;
	sParamsAnal.boolDecodeMD = true;
	
	%do analysis
	sOut = doFisherAnal(matData,vecTrialStimType,sParamsAnal);
	
	%% assign to output
	vecI_LR = nanmean(sOut.matI_LogReg_bc_CV,3);
	matI_LR(vecUseStimTypes(1),vecUseStimTypes(2),:) = vecI_LR;
	matI_LR(vecUseStimTypes(2),vecUseStimTypes(1),:) = vecI_LR;
	
	vecI_LR_shuff = nanmean(sOut.matI_LogReg_bc_CV_shuff,3);
	matI_LR_shuff(vecUseStimTypes(1),vecUseStimTypes(2),:) = vecI_LR_shuff;
	matI_LR_shuff(vecUseStimTypes(2),vecUseStimTypes(1),:) = vecI_LR_shuff;
	
	vecI_Dir = nanmean(sOut.matI_Direct_bc,3);
	matI_Dir(vecUseStimTypes(1),vecUseStimTypes(2),:) = vecI_Dir;
	matI_Dir(vecUseStimTypes(2),vecUseStimTypes(1),:) = vecI_Dir;
	vecI_Dir_shuff= nanmean(sOut.matI_Direct_bc_shuff,3);
	matI_Dir_shuff(vecUseStimTypes(1),vecUseStimTypes(2),:) = vecI_Dir_shuff;
	matI_Dir_shuff(vecUseStimTypes(2),vecUseStimTypes(1),:) = vecI_Dir_shuff;
	
	
	vecA_LR = nanmean(sOut.matA_LR,3);
	matA_LR(vecUseStimTypes(1),vecUseStimTypes(2),:) = vecA_LR;
	matA_LR(vecUseStimTypes(2),vecUseStimTypes(1),:) = vecA_LR;
	vecA_LR_shuff= nanmean(sOut.matA_LR_shuff,3);
	matA_LR_shuff(vecUseStimTypes(1),vecUseStimTypes(2),:) = vecA_LR_shuff;
	matA_LR_shuff(vecUseStimTypes(2),vecUseStimTypes(1),:) = vecA_LR_shuff;
	
	vecA_ML = nanmean(sOut.matA_ML,3);
	matA_ML(vecUseStimTypes(1),vecUseStimTypes(2),:) = vecA_ML;
	matA_ML(vecUseStimTypes(2),vecUseStimTypes(1),:) = vecA_ML;
	vecA_ML_shuff= nanmean(sOut.matA_ML_shuff,3);
	matA_ML_shuff(vecUseStimTypes(1),vecUseStimTypes(2),:) = vecA_ML_shuff;
	matA_ML_shuff(vecUseStimTypes(2),vecUseStimTypes(1),:) = vecA_ML_shuff;
	
	vecA_MD = nanmean(sOut.matA_MD,3);
	matA_MD(vecUseStimTypes(1),vecUseStimTypes(2),:) = vecA_MD;
	matA_MD(vecUseStimTypes(2),vecUseStimTypes(1),:) = vecA_MD;
	vecA_MD_shuff= nanmean(sOut.matA_MD_shuff,3);
	matA_MD_shuff(vecUseStimTypes(1),vecUseStimTypes(2),:) = vecA_MD_shuff;
	matA_MD_shuff(vecUseStimTypes(2),vecUseStimTypes(1),:) = vecA_MD_shuff;
	
	%% plot
	matMeanLRI = nanmean(matI_LR,3);

	
	subplot(2,5,1);
	imagesc(matMeanLRI);
	title(sprintf('Fisher LR; range: %.2f %.2f',get(gca,'clim')))
	subplot(2,5,6);
	imagesc(nanmean(matI_LR_shuff,3));
	title(sprintf('Fisher LR shuffled; range: %.2f %.2f',get(gca,'clim')))
	
	subplot(2,5,2);
	imagesc(nanmean(matI_Dir,3));
	title(sprintf('Fisher Direct; range: %.2f %.2f',get(gca,'clim')))
	subplot(2,5,7);
	imagesc(nanmean(matI_Dir_shuff,3));
	title(sprintf('Fisher Direct shuffled; range: %.2f %.2f',get(gca,'clim')))
	
	subplot(2,5,3);
	imagesc(nanmean(matA_LR,3));
	title(sprintf('Accuracy LR; range: %.2f %.2f',get(gca,'clim')))
	subplot(2,5,8);
	imagesc(nanmean(matA_LR_shuff,3));
	title(sprintf('Accuracy LR shuffled; range: %.2f %.2f',get(gca,'clim')))
	
	subplot(2,5,4);
	imagesc(nanmean(matA_ML,3));
	title(sprintf('Accuracy ML; range: %.2f %.2f',get(gca,'clim')))
	subplot(2,5,9);
	imagesc(nanmean(matA_ML_shuff,3));
	title(sprintf('Accuracy ML shuffled; range: %.2f %.2f',get(gca,'clim')))
	
	subplot(2,5,5);
	imagesc(nanmean(matA_MD,3));
	title(sprintf('Accuracy MD; range: %.2f %.2f',get(gca,'clim')))
	subplot(2,5,10);
	imagesc(nanmean(matA_MD_shuff,3));
	title(sprintf('Accuracy MD shuffled; range: %.2f %.2f',get(gca,'clim')))
	drawnow;
end
%% save figure
cellSplit = strsplit(strSimulation,'_');
strFile = sprintf('%s_InfoDecoding_PopSize%d_%s.tif',cellSplit{2},intPopSize,cellSplit{3});
drawnow;
jFig = get(handle(gcf), 'JavaFrame');
jFig.setMaximized(true);
figure(gcf);
drawnow;
export_fig([strFigDir strFile '.tif']);
export_fig([strFigDir strFile '.pdf']);
fprintf('Saved figure to %s [%s]\n',strFile,getTime);
return
%% rotate vectors in matrix to align 0
	matMeanLRICent = nan(size(matMeanLRI));
	for intPos=1:intParamNum
		matMeanLRICent(intPos,:) = matMeanLRI(intPos,mod(((1:intParamNum)-1)-(intParamNum/2)+intPos,intParamNum)+1);
	end
