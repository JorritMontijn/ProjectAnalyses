
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
	intUseCells = 256;
	intResamplings = 10;
	intIters = 10;
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
	if intUseCells > intCellsV1,intIters = 1;end
	
	%% add noise to improve numerical stability
	matModelResp = double(matModelResp);
	intNeurons = size(matModelResp,1);
	intTrials = size(matModelResp,2);
	vecNeuronSD = xstd(matModelResp,2);
	matNoise = normrnd(zeros(intNeurons,intTrials),repmat(vecNeuronSD*dblNoiseLevel,[1 intTrials]));
	matModelRespP = matModelResp + matNoise;
	%matModelRespP = zscore(matModelRespP,[],1);
	cellArea = {'V1','V2','V1-V2'};
	figure
	
	for intWithinArea=[1 2]
		%% pre-allocate output matrices
		matPercCosSimMean = nan(2,intParamNum,intIters);
		%matPercCosSimCI = nan(2,intParamNum,intIters);
		matPercCosSimOrthMean = nan(2,intParamNum,intIters);
		%matPercCosSimOrthCI = nan(2,intParamNum,intIters);
		
		%get projection onto F'; what is distribution compared to cluster separation?
		cellProjectedPointsFprime = cell(2,intParamNum,intIters);
		cellProjectedPointsFprimeShuffled = cell(2,intParamNum,intIters);
		
		matIncreaseInDiffCorrsMean = nan(2,intParamNum,intIters);
		%matIncreaseInDiffCorrsCI = nan(2,intParamNum,intIters);
		
		matDistanceBetweenClassMeans = nan(2,intParamNum,intIters);
		matVarianceOnFprime = nan(2,intParamNum,intIters);
		matFisherI = nan(2,intParamNum,intIters);
		for intIter=1:intIters
			%% compare different stimulus intensity levels; from blank white noise to structured stimulus
			for intStimTypeIdx=1:intParamNum
				%% start
				dblParamVal = vecParamVals(intStimTypeIdx);
				cellC{intStimTypeIdx} = sprintf('%03d',dblParamVal);
				strC = cellC{intStimTypeIdx};
				fprintf('Starting %s %s (%d/%d) [%s]\n',strParam,strC,intStimTypeIdx,intParamNum,getTime);
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
				if intUseCells < numel(vecCells),vecCells = vecCells(randperm(numel(vecCells),intUseCells));end
				
				if boolAllZero
					%get output
				else
					%% switch R1/R2
					for intStimSwitch=1:2
						if intStimSwitch==2,vecUseStimTypes=vecUseStimTypes([2 1]);end
						matData = matModelRespP(vecCells,:);
						matR1 = matData(:,vecTrialStimType==vecUseStimTypes(1));
						matR2 = matData(:,vecTrialStimType==vecUseStimTypes(2));
						vecMu1=mean(matR1,2);
						matR1 = bsxfun(@minus,matR1,vecMu1);
						matR2 = bsxfun(@minus,matR2,vecMu1);
						vecFprime = mean(matR2,2);
						
						%build orthogonal vector (or as close as we can get)
						vecF2 = vecFprime;
						intN = numel(vecCells);
						for intEl = 1:intN
							vecF2(intEl) = -vecFprime(intEl)*sum(vecFprime(~ismember(1:intN,1)).^2);
							if vecF2'*vecFprime < 0,break;end
						end
						[dblPercVarInDirOverShuffled,vecCIPCIDOS,dblCohensD] = getVarInDir(matR1',vecFprime',intResamplings);
						[dblPercVarInDirOverShuffledOrthogonal,vecCIPCIDOSO,dblCohensDOrth] = getVarInDir(matR1',vecF2',intResamplings);
						
						%get projection onto F'; what is distribution compared to
						%cluster separation?
						[vecProjectedPoints1D,dblRefNorm,matProjectedPointsShuffled,vecProjectedPointsClass2] = getProjOnLine(matR1',matR2',intResamplings);
						vecIncreaseInDiffCorrs = (xstd(vecProjectedPoints1D,1) ./ xstd(matProjectedPointsShuffled,1)) * 100;
						dblIncreaseInDiffCorrs = mean(vecIncreaseInDiffCorrs);
						vecIncreaseInDiffCorrCI = getCI(vecIncreaseInDiffCorrs);
						
						%% get actual Fisher information
						%get corr facs
						dblLambda=sParams.dblLambda;
						dblDiffTheta=sParams.dblDiffTheta;
						intGroupSize = size(matR1,1);
						intTrials12 = size(matR1,2); %check if this one
						dblSubFac =(2*intGroupSize)/(intTrials12*(dblDiffTheta.^2));
						dblProdFacRaw = ((2*intTrials12-intGroupSize-3)/(2*intTrials12-2));
						
						%get logistic regression output
						[vecWeightsLogReg, dblLLH] = doBinLogReg([matR1 matR2], [zeros(1,size(matR1,2)) ones(1,size(matR2,2))], dblLambda);
						
						%non-CV
						vecClass1NonCV = vecWeightsLogReg'*[matR1;ones(1,size(matR1,2))];
						vecClass2NonCV = vecWeightsLogReg'*[matR2;ones(1,size(matR2,2))];
						dblDprimeLogReg = getdprime2(vecClass1NonCV,vecClass2NonCV);
						%CV
						vecClass1 = vecWeightsLogReg'*[matR1;ones(1,size(matR1,2))];
						vecClass2 = vecWeightsLogReg'*[matR2;ones(1,size(matR2,2))];
						dblDprimeLogRegCV = getdprime2(vecClass1,vecClass2);
						dblI_bc = (dblDprimeLogRegCV.^2)*dblProdFacRaw-dblSubFac;
						
						%% assign output
						matPercCosSimMean(intStimSwitch,intStimTypeIdx,intIter) = dblPercVarInDirOverShuffled;
						%matPercCosSimCI(:,intStimTypeIdx,intIter) = vecCIPCIDOS;
						matPercCosSimOrthMean(intStimSwitch,intStimTypeIdx,intIter) = dblPercVarInDirOverShuffledOrthogonal;
						%matPercCosSimOrthCI(:,intStimTypeIdx,intIter) = vecCIPCIDOSO;
						
						%get projection onto F'; what is distribution compared to cluster separation?
						cellProjectedPointsFprime{intStimSwitch,intStimTypeIdx,intIter} = vecProjectedPoints1D;
						cellProjectedPointsFprimeShuffled{intStimSwitch,intStimTypeIdx,intIter} = matProjectedPointsShuffled;
						
						matIncreaseInDiffCorrsMean(intStimSwitch,intStimTypeIdx,intIter) = dblIncreaseInDiffCorrs;
						%matIncreaseInDiffCorrsCI(:,intStimTypeIdx,intIter) = vecIncreaseInDiffCorrCI;
						
						matDistanceBetweenClassMeans(intStimSwitch,intStimTypeIdx,intIter) = (norm(vecFprime) - mean(vecProjectedPoints1D));
						matVarianceOnFprime(intStimSwitch,intStimTypeIdx,intIter) = var(vecProjectedPoints1D);
						matFisherI(intStimSwitch,intStimTypeIdx,intIter) = dblI_bc;
					end
					%{
					%%
					figure
					subplot(2,2,1)
					[vecN,vecBins] = histx(vecProjectedPoints1D);
					plot(vecBins,vecN,'b');
					dblM = mean(vecProjectedPoints1D);
					dblS = std(vecProjectedPoints1D);
					vecG = normpdf(vecBins,dblM,dblS);
					vecNormG = (vecG/max(vecG))*max(vecN);
					hold on
					plot(vecBins,vecNormG,'r--')
					hold off
					title(sprintf('Class 2 mean at %.3f',dblRefNorm));
					
					subplot(2,2,2)
					[vecN2,vecBins2] = histx(vecProjectedPointsClass2);
					plot(vecBins2,vecN2,'b');
					dblM = mean(vecProjectedPointsClass2);
					dblS = std(vecProjectedPointsClass2);
					vecG = normpdf(vecBins2,dblM,dblS);
					vecNormG = (vecG/max(vecG))*max(vecN);
					hold on
					plot(vecBins2,vecNormG,'r--')
					hold off
					title(sprintf('Class 1 mean at %.3f',mean(vecProjectedPoints1D)));
					
					subplot(2,2,3)
					dblStep = mean([mean(diff(vecBins)) mean(diff(vecBins2))]);
					dblMin = min([min(vecBins) min(vecBins2)])-dblStep;
					dblMax = max([max(vecBins) max(vecBins2)])+dblStep;
					vecBinsBoth = dblMin:dblStep:dblMax;
					vecBinsBothPlot = vecBinsBoth(2:end)-(dblStep/2);
					vecNs = histcounts([vecProjectedPoints1D vecProjectedPointsClass2],vecBinsBoth);
					plot(vecBinsBothPlot,vecNs,'b');
					title(sprintf('Class 1 mean at %.3f',mean(vecProjectedPoints1D)));
					
					dblDprime = getdprime2(vecProjectedPoints1D,vecProjectedPointsClass2);
					
					%%
					matMahal11 = squeeze(matMahalDists(1,1,:,:,end));
					matMahal12 = squeeze(matMahalDists(1,2,:,:,end));
					
					vecBins=0:5:200;
					vecEdges = -2.5:5:202.5;
					mapC=redbluepurple(10);
					figure
					intI = 1;
					%for intI=1:9
					%subplot(3,3,intI)
					dblMax = max([max(matMahal11(:)) max(matMahal12(:))]);
					plot([0 dblMax],[0 dblMax],'k--')
					hold on
					scatter(matMahal11(:,intI),matMahal12(:,intI))
					hold off
					xlabel('Mahalanobis distance to correct class')
					ylabel('Distance to incorrect class')
					%end
					
					
					
					return
					%}
				end
			end
		end
		%% fuse stimuli
		matPercCosSimMean = cat(3,matPercCosSimMean(1,:,:),matPercCosSimMean(2,:,:));
		matPercCosSimOrthMean = cat(3,matPercCosSimOrthMean(1,:,:),matPercCosSimOrthMean(2,:,:));
		matVarianceOnFprime = cat(3,matVarianceOnFprime(1,:,:),matVarianceOnFprime(2,:,:));
		matDistanceBetweenClassMeans = cat(3,matDistanceBetweenClassMeans(1,:,:),matDistanceBetweenClassMeans(2,:,:));	
		matIncreaseInDiffCorrsMean = cat(3,matIncreaseInDiffCorrsMean(1,:,:),matIncreaseInDiffCorrsMean(2,:,:));	
		cellProjectedPointsFprime = cat(3,cellProjectedPointsFprime(1,:,:),cellProjectedPointsFprime(2,:,:));	
		cellProjectedPointsFprimeShuffled = cat(3,cellProjectedPointsFprimeShuffled(1,:,:),cellProjectedPointsFprimeShuffled(2,:,:));	
		
		matFisherI = cat(3,matFisherI(1,:,:),matFisherI(2,:,:));
		
						
		%% plot
		%extract
		vecPercCosSimMean = xmean(matPercCosSimMean,3);
		vecPercCosSimSE = xstd(matPercCosSimMean,3)/sqrt(intIters);
		vecPercCosSimOrthMean = xmean(matPercCosSimOrthMean,3);
		vecPercCosSimOrthSE = xstd(matPercCosSimOrthMean,3)/sqrt(intIters);
		
		vecVarianceOnFprimeMean = xmean(matVarianceOnFprime,3);
		vecVarianceOnFprimeSE = xstd(matVarianceOnFprime,3)/sqrt(intIters);
		vecDistanceBetweenClassMeansMean = xmean(matDistanceBetweenClassMeans,3);
		vecDistanceBetweenClassMeansSE = xstd(matDistanceBetweenClassMeans,3)/sqrt(intIters);
		
		vecFisherMean = xmean(matFisherI,3);
		vecFisherSE = xstd(matFisherI,3)/sqrt(intIters);
		
		
		
		% plot
		intNrCols = 4;
		strArea = cellArea{intWithinArea};
		subplot(2,intNrCols,(intWithinArea-1)*intNrCols+1);
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
		title([strType ' ' strArea '; ' num2str(intUseCells) ' cells'])
		fixfig
		
		subplot(2,intNrCols,(intWithinArea-1)*intNrCols+2);
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
		title(sprintf('Corr (d(mu) <-> Fisher I); %.3f)',nancorr(vecDistanceBetweenClassMeansMean',vecFisherMean')))
		fixfig
		
		vecPredI = (vecDistanceBetweenClassMeansMean ./sqrt(vecVarianceOnFprimeMean)).^2;
		subplot(2,intNrCols,(intWithinArea-1)*intNrCols+3);
		plot(vecParamVals,vecPredI)
		xlabel(strParam);
		ylabel('Predicted Fisher information (mu/sd)^2');
		ylim([0 max(get(gca,'ylim'))]);
		xlim([-1 101])
		title(sprintf('Corr (pred I <-> Fisher I); %.3f)',nancorr(vecPredI',vecFisherMean')))
		fixfig
		drawnow;
		
		subplot(2,intNrCols,(intWithinArea-1)*intNrCols+4);
		errorbar(vecParamVals,vecFisherMean,vecFisherSE);
		xlabel(strParam);
		ylabel('Real Fisher information (d''^2) (LogReg)');
		ylim([0 max(get(gca,'ylim'))]);
		xlim([-1 101])
		title('Note: not CV')
		fixfig
		drawnow;
		
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
	
	%full screen
	drawnow;
	jFig = get(handle(gcf), 'JavaFrame');
	jFig.setMaximized(true);
	figure(gcf);
	drawnow;
	if boolSaveFigs
		%save figure
		export_fig([strFigDir  'AnalBlock3FprimeVar_' strTag '.tif']);
		export_fig([strFigDir  'AnalBlock3FprimeVar_' strTag '.pdf']);
	end
end
