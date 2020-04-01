
%% initialize
clearvars;
boolLoad = false;
if ~exist('matModelResp','var')
	clearvars -except bool*;%close all;
	intLoadSim = 24; %% SET WHICH SIM TO LOAD HERE; 41 vs 42; 2 vs. 5
	boolLoad = true;
end


vecRunSims = [41 42];
for intLoadSim=vecRunSims
	clearvars -except vecRunSims intLoadSim
	boolLoad = true;
	if intLoadSim == 11 && boolLoad
		strSimulation = 'xAreaDistributed_OriFull_2017-06-15';
	elseif intLoadSim == 12 && boolLoad
		strSimulation = 'xAreaDistributed_ContExcOnly_2017-06-08'; %contrast
	elseif intLoadSim == 13 && boolLoad
		strSimulation = 'xAreaDistributed_ContInhOnly_2017-06-12'; %contrast
	elseif intLoadSim == 14 && boolLoad
		strSimulation = 'xAreaDistributed_ContNone_2017-06-12'; %contrast
		
	elseif intLoadSim == 21 && boolLoad
		strSimulation = 'xAreaDistributed_LumFull_2017-06-26'; %luminance full
	elseif intLoadSim == 22 && boolLoad
		strSimulation = ''; %lum exc
	elseif intLoadSim == 23 && boolLoad
		strSimulation = ''; %lum inh
	elseif intLoadSim == 24 && boolLoad
		strSimulation = 'xAreaDistributed_LumNone_2017-06-26'; %lum none
		
		
	elseif intLoadSim == 31 && boolLoad
		strSimulation = 'xAreaDistributed_ContFull_2017-06-26'; %contrast full
	elseif intLoadSim == 32 && boolLoad
		strSimulation = ''; %contrast exc
	elseif intLoadSim == 33 && boolLoad
		strSimulation = ''; %contrast inh
	elseif intLoadSim == 34 && boolLoad
		strSimulation = 'xAreaDistributed_ContNone_2017-06-26'; %contrast none
		
	elseif intLoadSim == 41 && boolLoad
		strSimulation = 'xAreaDistributed_Ori2OldConn_2017-09-18'; %old connectivity
	elseif intLoadSim == 42 && boolLoad
		strSimulation = 'xAreaDistributed_Ori2NewConn_2017-09-18'; %new connectivity
	
	end
	
	%% RUN: #header
	strFigDir = 'D:\Data\Results\Block3\';
	if boolLoad
		runModelHeader;
	end
	cellIn = strsplit(strSimulation,'_');
	strFiller = cellIn{1};
	strType = cellIn{2};
	strDate = cellIn{3};
	strTag = [strType strDate];
	
	%% set parameters
	boolSaveFigs = true;
	boolSaveData = true;
	dblCutOff = 0.90;
	dblLambdaInfo = 1;
	vecGroupSizes = 2.^(2:9);
	intResamplings = 10;
	intIters = 5;
	dblNoiseLevel = 0;
	boolAnalFprimeVar = true;
	%boolAnalPred = true;
	%boolAnalPopM = true;
	
	%% calculate which parameter to use
	vecRanges = range(matStimTypeCombos(2:end,:),2);
	[d,intUseParam]=max(vecRanges); %exclude orientation
	intUseParam = intUseParam + 1;
	if all(vecRanges==0), intUseParam = 0;vecParamStimType = [0];strParam='None';
		vecParamIdx = ones(1,size(matStimTypeCombos,2));
		vecParamVals = 1;
		intParamNum=1;
		mapC = redbluepurple(intParamNum);
	else
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
	end
	
	%% add noise to improve numerical stability
	matModelResp = double(matModelResp);
	intNeurons = size(matModelResp,1);
	intTrials = size(matModelResp,2);
	vecNeuronSD = xstd(matModelResp,2);
	matNoise = normrnd(zeros(intNeurons,intTrials),repmat(vecNeuronSD*dblNoiseLevel,[1 intTrials]));
	matModelRespP = matModelResp + matNoise;
	%matModelRespP = zscore(matModelRespP,[],1);
	cellArea = {'V1','V2','V1-V2'};
	
	
	intWithinArea=1;
	
	matMeanCorr = nan(numel(vecGroupSizes),intIters,intParamNum);
	matSDCorr = nan(numel(vecGroupSizes),intIters,intParamNum);
	matCovCorr = nan(numel(vecGroupSizes),intIters,intParamNum,3);
	for intSizeIdx=1:numel(vecGroupSizes)
		intSize = vecGroupSizes(intSizeIdx);
		fprintf('Starting pop size %d (%d/%d) [%s]\n',intSize,intSizeIdx,numel(vecGroupSizes),getTime);
		for intIter=1:intIters
			%% compare different stimulus intensity levels; from blank white noise to structured stimulus
			for intStimTypeIdx=1:intParamNum
				%% start
				dblParamVal = vecParamVals(intStimTypeIdx);
				cellC{intStimTypeIdx} = sprintf('%03d',dblParamVal);
				strC = cellC{intStimTypeIdx};
				
				vecUseStimTypes = find(sort(vecParamIdx,'ascend')==intStimTypeIdx);
				
				%check if no spikes
				if (sum(sum(matModelRespP(1:intCellsV1,vecTrialStimType==vecUseStimTypes(1)))) + sum(sum(matModelRespP(1:intCellsV1,vecTrialStimType==vecUseStimTypes(end))))) == 0
					boolAllZero = true;else boolAllZero = false;end
				
				%set parameters
				sParams = struct;
				sParams.vecUseStimTypes = vecUseStimTypes;
				sParams.dblLambda = dblLambdaInfo;
				sParams.dblDiffTheta = range(vecTrialOris);
				
				if intWithinArea == 1
					vecCells = 1:intCellsV1;
				elseif intWithinArea == 2
					vecCells = (intCellsV1+1):size(matModelRespP,1);
				end
				vecCells = vecCells(randperm(numel(vecCells),intSize));
				
				if boolAllZero
					%get output
				else
					%do analysis
					%% calculate mahalanobis distances (inherently heteroscadastic, so no need for stimulus switching)
					%[(real class) x (target class) x (trial) x (iteration) x (group size)]
					
					vecUseStimTypes=vecUseStimTypes([2 1]);
					matData = matModelRespP(vecCells,:);
					matR = matData(:,vecTrialStimType==vecUseStimTypes(1));
					
					intReps = size(matR,2);
					intUseReps = floor(intReps/2);
					matR1 = matR(:,1:intUseReps);
					matR2 = matR(:,(intUseReps+1):(intUseReps*2));
					%%
					matCovar1 = cov(matR1');
					matCovar2 = cov(matR2');
					
					[dblR,d,dblUp,dblDown]=corrcoef(matCovar1(:),matCovar2(:));
					dblR = dblR(1,end);
					dblUp = dblUp(1,end);
					dblDown = dblDown(1,end);
					
					matCovCorr(intSizeIdx,intIter,intStimTypeIdx,1) = dblR;
					matCovCorr(intSizeIdx,intIter,intStimTypeIdx,2) = dblUp;
					matCovCorr(intSizeIdx,intIter,intStimTypeIdx,3) = dblDown;
					
					
					vecMean1 = xmean(matR1,2);
					vecMean2 = xmean(matR2,2);
						
					vecSD1 = xstd(matR1,2);
					vecSD2 = xstd(matR2,2);
					
					matMeanCorr(intSizeIdx,intIter,intStimTypeIdx) = corr(vecMean1,vecMean2);
					matSDCorr(intSizeIdx,intIter,intStimTypeIdx) =  corr(vecSD1,vecSD2);
					%}
				end
			end
		end
	end
	%% plot
	figure
	matMeanCV = squeeze(mean(matCovCorr(:,:,:,1),2));
	cellGroups = cellfun(@num2str,num2cell(vecGroupSizes),'UniformOutput',false);
	intGroups = numel(vecGroupSizes);
	mapC = winter(intGroups);
	subplot(2,2,1)
	colormap(mapC)
	hold on
	for intSizeIdx=1:intGroups
		plot(vecParamVals,matMeanCV(intSizeIdx,:),'x','Color',mapC(intSizeIdx,:))
	end
	legend(cellGroups,'Location','Best')
	xlabel(strParam)
	ylabel('Correlation Covar of 1st vs 2nd half')
	title('Colors: different group sizes of random neurons')
	fixfig
	
	matMeanMu = squeeze(mean(matMeanCorr(:,:,:),2));
	subplot(2,2,2)
	colormap(mapC)
	hold on
	for intSizeIdx=1:intGroups
		plot(vecParamVals,matMeanMu(intSizeIdx,:),'x','Color',mapC(intSizeIdx,:))
	end
	legend(cellGroups,'Location','Best')
	xlabel(strParam)
	ylabel('Correlation Mean of 1st vs 2nd half')
	title('Colors: different group sizes of random neurons')
	fixfig
	
	matMeanSD = squeeze(mean(matSDCorr(:,:,:),2));
	subplot(2,2,3)
	colormap(mapC)
	hold on
	for intSizeIdx=1:intGroups
		plot(vecParamVals,matMeanSD(intSizeIdx,:),'x','Color',mapC(intSizeIdx,:))
	end
	legend(cellGroups,'Location','Best')
	xlabel(strParam)
	ylabel('Correlation SD of 1st vs 2nd half')
	title('Colors: different group sizes of random neurons')
	fixfig
	
	%legend(vecGroupSizes);
	
	return
	%% fuse stimuli
	matPercCosSimMean = cat(3,matPercCosSimMean(1,:,:),matPercCosSimMean(2,:,:));
	matPercCosSimOrthMean = cat(3,matPercCosSimOrthMean(1,:,:),matPercCosSimOrthMean(2,:,:));
	matVarianceOnFprime = cat(3,matVarianceOnFprime(1,:,:),matVarianceOnFprime(2,:,:));
	matDistanceBetweenClassMeans = cat(3,matDistanceBetweenClassMeans(1,:,:),matDistanceBetweenClassMeans(2,:,:));
	matIncreaseInDiffCorrsMean = cat(3,matIncreaseInDiffCorrsMean(1,:,:),matIncreaseInDiffCorrsMean(2,:,:));
	cellProjectedPointsFprime = cat(3,cellProjectedPointsFprime(1,:,:),cellProjectedPointsFprime(2,:,:));
	cellProjectedPointsFprimeShuffled = cat(3,cellProjectedPointsFprimeShuffled(1,:,:),cellProjectedPointsFprimeShuffled(2,:,:));
	
	%% plot
	%extract
	vecPercCosSimMean = xmean(matPercCosSimMean,3);
	vecPercCosSimSE = xstd(matPercCosSimMean,3)/sqrt(intIters);
	vecPercCosSimOrthMean = xmean(matPercCosSimOrthMean,3);
	vecPercCosSimOrthSE = xstd(matPercCosSimOrthMean,3)/sqrt(intIters);
	
	vecVarianceOnFprimeMean = xmean(matVarianceOnFprime.^2,3);
	vecVarianceOnFprimeSE = xstd(matVarianceOnFprime.^2,3)/sqrt(intIters);
	vecDistanceBetweenClassMeansMean = xmean(matDistanceBetweenClassMeans.^2,3);
	vecDistanceBetweenClassMeansSE = xstd(matDistanceBetweenClassMeans.^2,3)/sqrt(intIters);
	
	% plot
	strArea = cellArea{intWithinArea};
	subplot(2,2,(intWithinArea-1)*2+1);
	plot([0 100],[100 100],'k--')
	hold on;
	errorbar(vecParamVals,vecPercCosSimMean,vecPercCosSimSE,'r-')
	errorbar(vecParamVals,vecPercCosSimOrthMean,vecPercCosSimOrthSE,'g-')
	hold off;
	xlim([-1 101])
	ylim([0 max(get(gca,'ylim'))+1]);
	xlabel(strParam);
	ylabel('Cos sim relative to shuffled');
	legend({'','F''','Orthogonal to F'''},'Location','Best')
	title([strType ' ' strArea '; ' num2str(vecGroupSizes) ' cells'])
	fixfig
	
	subplot(2,4,(intWithinArea-1)*4+3);
	mapC = redbluepurple(intParamNum);
	
	errorbar(vecVarianceOnFprimeMean,vecDistanceBetweenClassMeansMean,vecDistanceBetweenClassMeansSE,'k.');
	hold on
	herrorbar(vecVarianceOnFprimeMean,vecDistanceBetweenClassMeansMean,vecVarianceOnFprimeSE,'k.');
	cline(vecVarianceOnFprimeMean,vecDistanceBetweenClassMeansMean,vecParamVals);colormap(redbluepurple)
	hold off
	ylim([0 max(get(gca,'ylim'))]);
	xlim([0 max(get(gca,'xlim'))]);
	ylabel('Distance between class means (spikes)');
	xlabel('Variance along F'' (spikes)');
	title([strType ' ' strArea])
	fixfig
	
	vecI = (vecDistanceBetweenClassMeansMean ./vecVarianceOnFprimeMean);
	subplot(2,4,(intWithinArea-1)*4+4);
	plot(vecParamVals,vecI)
	xlabel(strParam);
	ylabel('Predicted Fisher information (d''^2)');
	ylim([0 max(get(gca,'ylim'))]);
	xlim([-1 101])
	title([strType ' ' strArea])
	fixfig
	drawnow;
	return
	%% save data
	if boolSaveData
		strDataFile = ['D:\Data\Processed\V1_LIFmodel\Simulation_' strTag '_AB3_Area' num2str(intWithinArea) '.mat'];
		if exist(strDataFile,'file')
			sSave = load(strDataFile);
		else
			sSave = struct;
		end
		if boolAnalFprimeVar
			sSave = struct;
			sSave.strParam = strParam;
			sSave.strArea = cellArea{intWithinArea};
			sSave.strType = strType;
			sSave.vecParamVals = vecParamVals;
			sSave.sParams=sParams;
			
			sSave.vecPercCosSimMean = vecPercCosSimMean;
			sSave.vecPercCosSimSE = vecPercCosSimSE;
			sSave.vecPercCosSimOrthMean = vecPercCosSimOrthMean;
			sSave.vecPercCosSimOrthSE = vecPercCosSimOrthSE;
			
			sSave.cellProjectedPointsFprime = cellProjectedPointsFprime;
			sSave.cellProjectedPointsFprimeShuffled = cellProjectedPointsFprimeShuffled;
			
			sSave.vecVarianceOnFprimeMean = vecVarianceOnFprimeMean;
			sSave.vecVarianceOnFprimeSE = vecVarianceOnFprimeSE;
			sSave.vecDistanceBetweenClassMeansMean = vecDistanceBetweenClassMeansMean;
			sSave.vecDistanceBetweenClassMeansSE = vecDistanceBetweenClassMeansSE;
			
		end
		
		save(strDataFile,'-struct','sSave');
	end
end
