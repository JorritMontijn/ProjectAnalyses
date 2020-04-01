%{
Date		Type	BiDir			Plaid			Plane	Notes
20130612	SR101	xyt02			xyt01			1
20130625	PV		xyt02			xyt01			1
20131016	PV		xyt01/03		xyt02[/04:LQ]	1		xyt04 low qual
20131016	PV		xyt05							2
20131016	PV		xyt06							3
20131022	PV		xyt01			xyt02			1
20131022	PV		xyt03							2
20140129	None	xyt01/03/04		xyt02			1		No SR101
20140129	None	xyt05							2		No SR101
20140314	SR101	xyt08			xyt09			1
20140423	Both	xyt01/03		xyt02			1		PV and SR101!
20140425	SR101	xyt09			xyt10			1
%}
%stimuli:
%{
		Dur1 (3s)	Dur2 (3s) 	Dur3 (3s)					[+3 ITI]
	1)	Low on		High on		high off	low off			12s
	2)	High on		high off								6s
	3)	Low on								low off			12s
	4)	High on		low on		low off		high off		12s
	5)	low on		low off									6s
	6)	High on								high off		12s
	7)	H+L on		high off	low off						9s
	8)	H+L on		low off		high off					9s
	9)	Low on		low on		low off		low off			12s
	10)	High on		High on		high off	high off		12s
	11)	L+L on		L+L off									6s
	12)	H+H on		H+H off									6s
%}

%{
%step 1
load data

%step 2

%reformat data; output:
sTrialType(intTrialType).sPrimOri(intOriType).sEpoch(intEpoch).matFrames(intRepetition,intFrame)

intTrialType = [1:12]; trial types as above
intOriType = [1 2]; 1=0; 2=90
intEpoch = [1:3]; 1-3: short, mid, long
intFrame = frame number (1 - max in epoch)
intRepetition = repetition number (1-4)

matFrames: matrix containing corresponding original frame numbers
%}
%% set data environment
clear all;
close all;
strMasterPath = 'D:\Data\Processed\imagingdata\';
cellRec =		{'20130625','20131016','20131022','20140423'};
cellGrating =	{'xyt02',	'xyt01',	'xyt01',	'xyt02'};
cellPlaid =		{'xyt01',	'xyt02',	'xyt02',	'xyt01'};
vecUse = [1 2 3 4];

%cellRec = {'20130625','20140423'};
%indGratingFirst = [0 0];
matAggContRespNormPyr = [];
matAggContRespNormPV = [];

%% run
for intRec=vecUse
	% get data
	strDate = cellRec{intRec};
	figure;
	strRecGratings = cellGrating{intRec};
	strRecPlaids = cellPlaid{intRec};
	
	%set inclusion fraction
	dblIncludeFraction = 0.5;
	
	%load gratings
	strTargetDirGratings = [strMasterPath strDate filesep strRecGratings filesep];
	strTargetFileGratings = [strDate strRecGratings '_ses.mat'];
	sLoad = load([strTargetDirGratings strTargetFileGratings]);
	sesGratings = sLoad.ses;
	
	%load plaids
	strTargetDirPlaids = [strMasterPath strDate filesep strRecPlaids filesep];
	strTargetFilePlaids = [strDate strRecPlaids '_ses.mat'];
	sLoad = load([strTargetDirPlaids strTargetFilePlaids]);
	sesPlaids = sLoad.ses;
	
	for intSwitchPV = [0 1]
		if intSwitchPV == 1
			% switch to PV
			sesGratings.neuron = sesGratings.PV;
			sesPlaids.neuron = sesPlaids.PV;
			vecColor = [1 0.1 0.1];
			strTitle = 'PV';
			intPlot = 2;
		else
			vecColor = [0.1 0.1 1];
			strTitle = 'Non-PV';
			intPlot = 1;
		end
		
		%% remove absent cells
		intNeurons = numel(sesGratings.neuron);
		vecPresence = nan(1,intNeurons);
		for intNeuron=1:intNeurons
			vecPresence(intNeuron) = sesGratings.neuron(intNeuron).intPresence & sesPlaids.neuron(intNeuron).intPresence;
		end
		sesGratings.neuron(~vecPresence) = [];
		sesPlaids.neuron(~vecPresence) = [];
		
		%% transform gratings to resp mat
		intPBR = 75;
		sP.intPreBaselineRemoval = intPBR;
		matRespG = getRespMat(sesGratings,[],[],sP); %G=Gratings
		%matRespG = getRespMat(sesGratings); %G=Gratings
		intNeurons = size(matRespG,1);
		vecOriG = getOriListFromTrials(sesGratings.structStim.Orientation);
		sTypesG = getStimulusTypes(sesGratings.structStim);
		cellSelectG = getSelectionVectors(sesGratings.structStim,sTypesG);
		sTuningG = calcTuningRespMat(matRespG,cellSelectG,vecOriG);
		vecSort = sort(sTuningG.vecOSI,'descend');
		if intSwitchPV
			dblThresholdOSI = 0;
		else
			dblThresholdOSI = vecSort(round(intNeurons*dblIncludeFraction));
		end
		indUseNeurons = sTuningG.vecOSI >= dblThresholdOSI;
		indPref0 = sTuningG.vecPrefAngle < abs(sTuningG.vecPrefAngle-90);
		vecBinOriPrefIdx = (indPref0==0)*find(vecOriG==0) + (indPref0==1)*find(vecOriG==90);
		sP2.intPreBaselineRemoval = intPBR;
		matRespC0 = getRespMat(sesGratings,[],-2,sP2);
		
		%% transform resp mat to stimresp
		matStimResponse = doMatRespTransform(matRespG,cellSelectG);
		matMeanResp = squeeze(xmean(matStimResponse,2));
		intOriReps = size(matStimResponse,2);
		matRespCont100 = nan(intOriReps,intNeurons);
		matRespCont0 = nan(intOriReps,intNeurons);
		for intNeuron=1:intNeurons
			matRespCont100(:,intNeuron) = matStimResponse(vecBinOriPrefIdx(intNeuron),:,intNeuron);
			matRespCont0(:,intNeuron) = matRespC0(intNeuron,cellSelectG{vecBinOriPrefIdx(intNeuron)});
		end
		matRespCont100(:,~indUseNeurons) = [];
		matRespCont0(:,~indUseNeurons) = [];
		
		%% transform plaids to resp mat
		sTrialType = doReformatPlaidData(sesPlaids);
		intReps = size(sTrialType(1).sPrimOri(1).sEpoch(1).matFrames,1);
		
		vecOriTypes = [1 2];
		vecLowTrials = [1, 3, 5, 9];
		vecHighTrials = [2, 4, 6, 10];
		%sTrialType([1:6,9,10]).sPrimOri([1,2]).sEpoch(1).matFrames([1:6],[1:77])
		
		matMeanActNoneLow{1} = nan(intNeurons,numel(vecLowTrials)*intReps);
		matMeanActNoneLow{2} = nan(intNeurons,numel(vecLowTrials)*intReps);
		matMeanActNoneHigh{1} = nan(intNeurons,numel(vecHighTrials)*intReps);
		matMeanActNoneHigh{2} = nan(intNeurons,numel(vecHighTrials)*intReps);
		matMeanActLow{1} = nan(intNeurons,numel(vecLowTrials)*intReps);
		matMeanActLow{2} = nan(intNeurons,numel(vecLowTrials)*intReps);
		matMeanActHigh{1} = nan(intNeurons,numel(vecHighTrials)*intReps);
		matMeanActHigh{2} = nan(intNeurons,numel(vecHighTrials)*intReps);
		matActContrastRepNeuron = nan(2,numel(vecLowTrials)*intReps,sum(indUseNeurons));
		matActNoneRepNeuron = nan(2,numel(vecLowTrials)*intReps,sum(indUseNeurons));
		for	intNeuron=1:intNeurons
			for intHorVertIdx = 1:numel(vecOriTypes)
				intRepIdxLow = 0;
				intRepIdxHigh = 0;
				for intLowTrialIdx=1:numel(vecLowTrials)
					matRepByFrame = sTrialType(vecLowTrials(intLowTrialIdx)).sPrimOri(intHorVertIdx).sEpoch(1).matFrames;
					for intRep=1:intReps
						intRepIdxLow = intRepIdxLow + 1;
						vecPreFrames = round(matRepByFrame(intRep,1)-sesPlaids.samplingFreq):(matRepByFrame(intRep,1)-1);
						matMeanActNoneLow{intHorVertIdx}(intNeuron,intRepIdxLow) = mean(sesPlaids.neuron(intNeuron).dFoF(vecPreFrames));
						matMeanActLow{intHorVertIdx}(intNeuron,intRepIdxLow) = mean(sesPlaids.neuron(intNeuron).dFoF(matRepByFrame(intRep,:)));
					end
				end
				for intHighTrialIdx=1:numel(vecHighTrials)
					matRepByFrame = sTrialType(vecHighTrials(intHighTrialIdx)).sPrimOri(intHorVertIdx).sEpoch(1).matFrames;
					for intRep=1:intReps
						intRepIdxHigh = intRepIdxHigh + 1;
						vecPreFrames = round(matRepByFrame(intRep,1)-sesPlaids.samplingFreq):(matRepByFrame(intRep,1)-1);
						matMeanActNoneHigh{intHorVertIdx}(intNeuron,intRepIdxHigh) = mean(sesPlaids.neuron(intNeuron).dFoF(vecPreFrames));
						matMeanActHigh{intHorVertIdx}(intNeuron,intRepIdxHigh) = mean(sesPlaids.neuron(intNeuron).dFoF(matRepByFrame(intRep,:)));
					end
				end
			end
			matActContrastRepNeuron(1,:,intNeuron) = matMeanActLow{indPref0(intNeuron)+1}(intNeuron,:);
			matActContrastRepNeuron(2,:,intNeuron) = matMeanActHigh{indPref0(intNeuron)+1}(intNeuron,:);
			
			matActNoneRepNeuron(1,:,intNeuron) = matMeanActNoneLow{indPref0(intNeuron)+1}(intNeuron,:);
			matActNoneRepNeuron(2,:,intNeuron) = matMeanActNoneHigh{indPref0(intNeuron)+1}(intNeuron,:);
		end
		matRespContLowHigh = matActContrastRepNeuron(:,:,indUseNeurons);%-matActNoneRepNeuron(:,:,indUseNeurons);
		matRespContNone = matActNoneRepNeuron(:,:,indUseNeurons);
		
		%% aggregate
		vecMeanC0 = zeros(1,sum(indUseNeurons));
		vecSDC0 = zeros(1,sum(indUseNeurons));
		vecMeanC25 = squeeze(xmean(matRespContLowHigh(1,:,:),2))';
		vecSDC25 = squeeze(xstd(matRespContLowHigh(1,:,:),2))';
		vecMeanC50 = squeeze(xmean(matRespContLowHigh(2,:,:),2))';
		vecSDC50 = squeeze(xstd(matRespContLowHigh(2,:,:),2))';
		vecMeanC100 = xmean(matRespCont100-matRespCont0,1);
		vecSDC100 = xstd(matRespCont100-matRespCont0,1);
		
		matMeanContResp = cat(1,vecMeanC0,vecMeanC25,vecMeanC50,vecMeanC100);
		matSDContResp = cat(1,vecSDC0,vecSDC25,vecSDC50,vecSDC100);
		
		%normalize
		matMeanContRespNorm = matMeanContResp - matMeanContResp(1,:);
		matMeanContRespNorm = matMeanContRespNorm./max(matMeanContRespNorm);
		indPostRem = matMeanContRespNorm(4,:) < matMeanContRespNorm(1,:) | any(matMeanContRespNorm < -0.4,1);
		matMeanContRespNorm(:,indPostRem) = [];
		
		%% plot
		vecContrasts = [0 25 50 100];
		subplot(2,2,intPlot);
		cla;
		plot(vecContrasts,matMeanContRespNorm,'Color',[0.5 0.5 0.5])
		hold on
		errorbar(vecContrasts,xmean(matMeanContRespNorm,2),xstd(matMeanContRespNorm,2)/sqrt(sum(indUseNeurons)),'Color',vecColor)
		ylabel('Norm. Resp.');
		xlabel('Stim contrast (%)');
		set(gca,'xtick',vecContrasts);
		title([strTitle '; ' strDate]);
		fixfig
		drawnow;
		
		%% save data
		if intSwitchPV == 0
			matAggContRespNormPyr = cat(2,matAggContRespNormPyr,matMeanContRespNorm);
		elseif intSwitchPV == 1
			matAggContRespNormPV = cat(2,matAggContRespNormPV,matMeanContRespNorm);
		end
	end
end

%% plot aggregate
figure
for intSwitchAggPV = [0 1]
	if intSwitchAggPV == 1
		% switch to PV
		matPlotData = matAggContRespNormPV;
		vecColor = [1 0.1 0.1];
		strTitle = 'PV';
		intPlot = 2;
		
	elseif intSwitchAggPV == 0
		matPlotData = matAggContRespNormPyr;
		vecColor = [0.1 0.1 1];
		strTitle = 'Non-PV';
		intPlot = 1;
	end
	subplot(2,2,intPlot);
	
	
	cla;
	plot(vecContrasts,matPlotData,'Color',[0.5 0.5 0.5])
	hold on
	errorbar(vecContrasts,xmean(matPlotData,2),xstd(matPlotData,2)/sqrt(size(matPlotData,2)),'Color',vecColor)
	ylabel('Norm. Resp.');
	xlabel('Stim contrast (%)');
	set(gca,'xtick',vecContrasts);
	title(['Aggregate; ' strTitle]);
	fixfig
end
[h,vecP]=ttest2(matAggContRespNormPyr',matAggContRespNormPV');

subplot(2,2,3);
cla;
hold on
errorbar(vecContrasts,xmean(matAggContRespNormPV,2),xstd(matAggContRespNormPV,2)/sqrt(size(matPlotData,2)),'Color',[1 0.1 0.1])
errorbar(vecContrasts,xmean(matAggContRespNormPyr,2),xstd(matAggContRespNormPyr,2)/sqrt(size(matPlotData,2)),'Color',[0.1 0.1 1])
ylabel('Norm. Resp.');
xlabel('Stim contrast (%)');
set(gca,'xtick',vecContrasts);
title([sprintf('p=%.3f;',vecP(2:end))]);
fixfig

%% plot 2
matAggContRespNormNormPyr = bsxfun(@rdivide,matAggContRespNormPyr,mean(matAggContRespNormPyr(end,:)));
matAggContRespNormNormPV = bsxfun(@rdivide,matAggContRespNormPV,mean(matAggContRespNormPV(end,:)));
figure
for intSwitchAggPV = [0 1]
	if intSwitchAggPV == 1
		% switch to PV
		matPlotData = matAggContRespNormNormPV;
		vecColor = [1 0.1 0.1];
		strTitle = 'PV';
		intPlot = 2;
		
	elseif intSwitchAggPV == 0
		matPlotData = matAggContRespNormNormPyr;
		vecColor = [0.1 0.1 1];
		strTitle = 'Non-PV';
		intPlot = 1;
	end
	subplot(2,2,intPlot);
	
	
	cla;
	plot(vecContrasts,matPlotData,'Color',[0.5 0.5 0.5])
	hold on
	errorbar(vecContrasts,xmean(matPlotData,2),xstd(matPlotData,2)/sqrt(size(matPlotData,2)),'Color',vecColor)
	ylabel('Norm. Resp.');
	xlabel('Stim contrast (%)');
	set(gca,'xtick',vecContrasts);
	title(['Aggregate; ' strTitle]);
	fixfig
end
[h,vecP]=ttest2(matAggContRespNormNormPyr',matAggContRespNormNormPV');

subplot(2,2,3);
cla;
hold on
errorbar(vecContrasts,xmean(matAggContRespNormNormPV,2),xstd(matAggContRespNormNormPV,2)/sqrt(size(matPlotData,2)),'Color',[1 0.1 0.1])
errorbar(vecContrasts,xmean(matAggContRespNormNormPyr,2),xstd(matAggContRespNormNormPyr,2)/sqrt(size(matPlotData,2)),'Color',[0.1 0.1 1])
ylabel('Norm. Resp.');
xlabel('Stim contrast (%)');
set(gca,'xtick',vecContrasts);
title([sprintf('p=%.3f;',vecP(2:end))]);
fixfig

vecAUC_PV = sum(matAggContRespNormNormPV,1);
vecAUC_PV = vecAUC_PV;
vecAUC_Pyr = sum(matAggContRespNormNormPyr,1);
vecAUC_Pyr = vecAUC_Pyr;
[h,pAUC]=ttest2(vecAUC_PV,vecAUC_Pyr);
subplot(2,2,4);
vecX = [1/3 2/3];
errorbar(vecX,[mean(vecAUC_Pyr) mean(vecAUC_PV)],[std(vecAUC_Pyr)/numel(vecAUC_Pyr) std(vecAUC_PV)/numel(vecAUC_PV)],'Linestyle','none');
set(gca,'xtick',vecX,'xticklabel',{'Non-PV','PV'});
ylim([0 3])
xlim([0 1])
title(sprintf('Mean AUC Non-PV=%.3f,PV=%.3f,p=%.3f',mean(vecAUC_Pyr),mean(vecAUC_PV),pAUC));
fixfig
