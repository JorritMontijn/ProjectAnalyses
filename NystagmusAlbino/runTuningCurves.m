%% load spiking data & plot tuning curves

%% set recording
%close all;
cellUniqueAreas = {'V1','SC','NOT'};
strArea = ''; %V1, SC, NOT
strDataSourcePath = 'D:\Data\Processed\ePhys\DriftingGratings\';
strDataTargetPath = 'D:\Data\ResultsNystagmus\TuningMetrics\';
strFigPath = 'D:\Data\ResultsNystagmus\TuningCurves\';
boolMakePlots = true;

%% find data
sFiles = dir([strDataSourcePath '*' strArea '_*.mat']);
cellFiles = {sFiles(:).name}';

%% go through files
sAggRec = struct;
sAggNeuron = struct;
for intFile=1:numel(cellFiles)
	%% load
	sLoad = load([strDataSourcePath cellFiles{intFile}]);
	sNeuron = sLoad.sNeuron;
	sRecording = sLoad.sRecording;
	
	%% aggregate data
	if intFile == 1
		sAggRec= sRecording;
		sAggNeuron = sNeuron;
	else
		sAggRec(end+1) = sRecording;
		sAggNeuron((end+1):(end+numel(sNeuron))) = sNeuron;
	end
end

%% pre-allocate output variables
intNeurons = numel(sAggNeuron);
vecPrefDir = nan(1,intNeurons);
cellArea = cell(1,intNeurons);
matFitParams = nan(5,intNeurons);
vecFitOriDegs = 1:360;
matFitResp = nan(numel(vecFitOriDegs),intNeurons);
vecRho = nan(1,intNeurons);
vecDeltaPrime = nan(1,intNeurons);
cellRawOri = cell(1,intNeurons);
cellRawMean = cell(1,intNeurons);
cellRawSD = cell(1,intNeurons);
vecPeaksCV = nan(2,intNeurons);
vecPeaksBW = nan(2,intNeurons);
vecVisZ = nan(1,intNeurons);
vecVisHzD = nan(1,intNeurons);
cellVisZ = cell(1,intNeurons);
vecNonStatIdx = nan(1,intNeurons);
vecViolIdx = nan(1,intNeurons);
	
%% analyze
cellRecIdx = {sAggRec.RecIdx};
for intNeuron=[1:intNeurons]%31
	%% message
	fprintf('Processing neuron %d/%d [%s]\n',intNeuron,intNeurons,getTime);
	
	%% get neuronal data
	sThisNeuron = sAggNeuron(intNeuron);
	vecSpikeTimes = sThisNeuron.SpikeTimes;
	strRecIdx = sThisNeuron.RecIdx;
	strMouse = sThisNeuron.Mouse;
	strBlock = sThisNeuron.Block;
	strArea = sThisNeuron.Area;
	strDate = sThisNeuron.Date;
	intSU = sThisNeuron.IdxSU;
	intClust = sThisNeuron.IdxClust;
	
	%% get matching recording data
	sThisRec = sAggRec(strcmpi(strRecIdx,cellRecIdx));
	vecStimOnTime = sThisRec.vecStimOnTime;
	vecStimOffTime = sThisRec.vecStimOffTime;
	vecStimOriDegrees = sThisRec.vecStimOriDegrees;
	vecStimOriRads = deg2rad(vecStimOriDegrees);
	vecEyeTimestamps = sThisRec.vecEyeTimestamps;
	matEyeData = sThisRec.matEyeData;
	
	%% get trial data
	vecStimCounts = getSpikeCounts(vecSpikeTimes,vecStimOnTime,vecStimOffTime);
	vecStimResp = bsxfun(@rdivide,vecStimCounts,(vecStimOffTime-vecStimOnTime)'); %transform to Hz
	dblBaseDur = 0.5;
	vecBaseCounts = getSpikeCounts(vecSpikeTimes,vecStimOnTime-dblBaseDur,vecStimOnTime);
	vecBaseResp = bsxfun(@rdivide,vecBaseCounts,dblBaseDur); %transform to Hz
	if 0
	%% get tuning curves & parameters
	dblRho = getTuningRho(vecStimResp,vecStimOriRads);
	dblDeltaPrime = getDeltaPrime(vecStimResp,vecStimOriRads);
	sTuning = getTuningCurves(vecStimResp,vecStimOriDegrees);
	[matRespNSR,vecStimTypes,vecUniqueDegs] = getStimulusResponses(vecStimResp,vecStimOriDegrees);
		
	%% save data
	cellArea{intNeuron} = strArea;
	vecParams = sTuning.matFittedParams;
	matFitParams(:,intNeuron) = vecParams;
	vecPrefDir(intNeuron) = vecParams(1);
	matFitResp(:,intNeuron) = feval(sTuning.funcFit,vecParams,deg2rad(vecFitOriDegs));
	%vecParams(1) = pi/2;
	%rho
	vecRho(intNeuron) = dblRho;
	%d'
	vecDeltaPrime(intNeuron) = dblDeltaPrime;
	%raw ori, mean per ori, and sd per ori
	cellRawOri{intNeuron} = sTuning.vecUniqueRads;
	cellRawMean{intNeuron} = sTuning.matMeanResp;
	cellRawSD{intNeuron} = sTuning.matSDResp;
	%peak CV + BW
	vecPeaksCV(:,intNeuron) = sTuning.matVariance;
	vecPeaksBW(:,intNeuron) = sTuning.matBandwidth;
	
	%%% plot raw tuning curves
	%vecPlotTheta = [cellRawOri{intNeuron} cellRawOri{intNeuron}(1)];
	%vecPlotRho = cellRawMean{intNeuron}/max(cellRawMean{intNeuron});
	%vecPlotRho(end+1) = vecPlotRho(1);
	%intPlot = find(ismember(cellUniqueAreas,cellArea{intNeuron}));
	%axes(vecH(intPlot));
	%polarplot(vecPlotTheta,vecPlotRho);
	
	
	%% get cluster quality
	%build figure name
	strFileName = sprintf('%s%sB%sSU%dC%d',strArea,strDate,strBlock,intSU,intClust);
	boolMakePlotsCQ = true;
	sClustQual = getClusterQuality(vecSpikeTimes,vecStimOnTime,boolMakePlotsCQ);
	vecNonStatIdx(intNeuron) = sClustQual.dblNonstationarityIndex;
	vecViolIdx(intNeuron) = sClustQual.dblViolIdx1ms;
	
	%% plot neuron
	if boolMakePlotsCQ
	
	%add orientation tuning curve
	subplot(2,2,4)
	vecMean = nanmean(matRespNSR,3);
	vecSD = nanstd(matRespNSR,[],3);
	[dblMaxResp,intMaxIdx] = max(vecMean);
	vecShiftDegs = circshift(1:numel(vecUniqueDegs),[0 (round(numel(vecUniqueDegs)/2))-intMaxIdx]);
	errorbar(vecUniqueDegs-180,vecMean(vecShiftDegs),vecSD(vecShiftDegs))
	xlim([-180 180]);
	set(gca,'xtick',-180:90:180);
	ylim([0 max(get(gca,'ylim'))]);
	xlabel('Angle from strongest response (deg)');
	ylabel('Firing rate (Hz)');
	title(sprintf('%s %sB%s, SU %d (C%d), %s''=%.3f, %s=%.3f',strArea,strDate,strBlock,intSU,intClust,getGreek('delta','lower'),dblDeltaPrime,getGreek('rho','lower'),dblRho));
	fixfig
	
	%save plot
	drawnow;
	export_fig([strFigPath strFileName 'ClustQual.tif']);
	export_fig([strFigPath strFileName 'ClustQual.pdf']);
	end
	end
	%% get visual responsiveness
	matStimOnOff = cat(2,vecStimOnTime,vecStimOffTime);
	[dblZ,sOptOut] = getZeta(vecSpikeTimes,matStimOnOff,1);%boolMakePlots);
	dblHzD = sOptOut.dblHzD;
	%[dblZ,vecInterpT,vecZ,matDiffTest,dblHzD,dblP] = getVisualResponsivenessCurv(vecSpikeTimes,matStimOnOff,1);%boolMakePlots);
	%[dblZ,vecInterpT,vecZ,matDiffTest,dblHzD,dblP] = getVisualResponsiveness(vecSpikeTimes,matStimOnOff,1);%boolMakePlots);
	hAx=get(gcf,'Children');
	axes(hAx(end));
	dblD = 0;
	title(sprintf('%s %sB%s, SU %d (C%d), d=%.3f, d(Hz)=%.3f',strArea,strDate,strBlock,intSU,intClust,dblZ,dblHzD));
	% assign data
	vecVisZ(intNeuron) = dblZ;
	vecVisHzD(intNeuron) = dblHzD;
	cellVisZ{intNeuron} = sOptOut.vecZ;
	
	%% save plot
	if boolMakePlots
	drawnow;
	
	strRecIdx = sThisNeuron.RecIdx;
	strMouse = sThisNeuron.Mouse;
	strBlock = sThisNeuron.Block;
	strArea = sThisNeuron.Area;
	strDate = sThisNeuron.Date;
	intSU = sThisNeuron.IdxSU;
	intClust = sThisNeuron.IdxClust;
	
	strFileName = sprintf('%s%sB%sSU%dC%d',strArea,strDate,strBlock,intSU,intClust);
	export_fig([strFigPath strFileName 'VisResp3.tif']);
	export_fig([strFigPath strFileName 'VisResp3.pdf']);
	
	pause(0.1);
	close
	close
	end
	%return
end
return
%% save data
sNeuron = rmfield(sAggNeuron,'SpikeTimes');
strTargetFile = ['AllData' getDate '.mat'];
save([strDataTargetPath strTargetFile],'sNeuron','vecPrefDir','cellArea',...
	'matFitParams','vecRho','vecDeltaPrime','cellRawOri','cellRawMean',...
	'cellRawSD','vecPeaksCV','vecPeaksBW','vecVisD','vecVisHzD',...
	'vecNonStatIdx','vecViolIdx');

%% plot metrics
vecSpikeTot = cellfun(@numel,{sAggNeuron(:).SpikeTimes});
vecVisD_Plot = abs(vecVisZ);
%vecVisD_Plot = real(sqrt(vecVisD_Plot)) - imag(sqrt(vecVisD_Plot));
figure
subplot(2,2,1)
scatter(vecVisHzD,vecVisD_Plot)
%set(gca,'xscale','log','yscale','log')
xlabel('Vis Hz D')
ylabel('Vis D')

indInclude = vecDeltaPrime > 0.5 | vecVisD_Plot > 0.5;
subplot(2,2,2)
scatter(vecDeltaPrime,vecVisD_Plot,[],indInclude)
%set(gca,'xscale','log','yscale','log')
xlabel('delta-prime')
ylabel('Vis D')

subplot(2,2,3)
scatter(vecVisHzD,vecVisD_Plot)
%set(gca,'xscale','log','yscale','log')
xlabel('Vis Hz D')
ylabel('Vis D')
vecLim = [-2 2];
xlim(vecLim);
ylim(vecLim);

indInclude = vecDeltaPrime > 0.5 | vecVisD_Plot > 0.5;
subplot(2,2,4)
scatter(vecDeltaPrime,vecVisD_Plot,[],indInclude)
%set(gca,'xscale','log','yscale','log')
xlabel('delta-prime')
ylabel('Vis D')
xlim(vecLim);
ylim(vecLim);

return
indReview = vecSpikeTot < 1500 & ~indInclude;
for intNeuron=find(indReview)
	%% get neuronal data
	sThisNeuron = sAggNeuron(intNeuron);
	vecSpikeTimes = sThisNeuron.SpikeTimes;
	strRecIdx = sThisNeuron.RecIdx;
	strMouse = sThisNeuron.Mouse;
	strBlock = sThisNeuron.Block;
	strArea = sThisNeuron.Area;
	strDate = sThisNeuron.Date;
	intSU = sThisNeuron.IdxSU;
	intClust = sThisNeuron.IdxClust;
	
	%% get matching recording data
	sThisRec = sAggRec(strcmpi(strRecIdx,cellRecIdx));
	vecStimOnTime = sThisRec.vecStimOnTime;
	sClustQual = getClusterQuality(vecSpikeTimes,vecStimOnTime,true);
	title(sprintf('%d; %s %sB%s, SU %d (C%d)',intNeuron,strArea,strDate,strBlock,intSU,intClust));
	pause
end

%% manually add edge cases
indInclude([56 67 114 115 121]) = true;
return
%% plot d' and pref angle distros
%num
intN_V1 = sum(vecDprimeV1);
intN_SC = sum(vecDprimeSC);
intN_NOT = sum(vecDprimeNOT);

%d'
dblStep = 0.5;
vecBinEdgesDprime = -1:dblStep:3;
vecBinCentersDprime = vecBinEdgesDprime(2:end) - dblStep/2;
vecDprimeNOT = histcounts(vecDeltaPrime(ismember(cellArea,'NOT')),vecBinEdgesDprime);
vecDprimeSC = histcounts(vecDeltaPrime(ismember(cellArea,'SC')),vecBinEdgesDprime);
vecDprimeV1 = histcounts(vecDeltaPrime(ismember(cellArea,'V1')),vecBinEdgesDprime);

%pref ori
dblStep = pi/6;
vecBinEdgesPrefOri = 0:dblStep:(2*pi);
vecBinCentersPrefOri = vecBinEdgesPrefOri(2:end) - dblStep/2;
vecPrefAnglesNOT = vecPrefDir(ismember(cellArea,'NOT'));
vecPrefAnglesSC = vecPrefDir(ismember(cellArea,'SC'));
vecPrefAnglesV1 = vecPrefDir(ismember(cellArea,'V1'));

%plot
figure
%v1
vecLimY = [0 0.6];
subplot(2,3,1)
bar(vecBinCentersDprime,vecDprimeV1/intN_V1)
title(sprintf('V1; %d cells',intN_V1));
xlabel(sprintf('Orientation tuning (%s'')',getGreek(4,'lower')));
ylabel('Fraction of V1 cells');
ylim(vecLimY);
fixfig

subplot(2,3,4)
polarhistogram(vecPrefAnglesV1,vecBinEdgesPrefOri);
title('Preferred angle, V1')
fixfig

%SC
subplot(2,3,2)
bar(vecBinCentersDprime,vecDprimeSC/sum(vecDprimeSC))
title(sprintf('SC; %d cells',intN_SC));
xlabel(sprintf('Orientation tuning (%s'')',getGreek(4,'lower')));
ylabel('Fraction of SC cells');
ylim(vecLimY);
fixfig

subplot(2,3,5)
polarhistogram(vecPrefAnglesSC,vecBinEdgesPrefOri);
title('Preferred angle, SC')
fixfig

%NOT
subplot(2,3,3)
bar(vecBinCentersDprime,vecDprimeNOT/sum(vecDprimeNOT))
title(sprintf('NOT; %d cells',intN_NOT));
xlabel(sprintf('Orientation tuning (%s'')',getGreek(4,'lower')));
ylabel('Fraction of NOT cells');
ylim(vecLimY);
fixfig

subplot(2,3,6)
polarhistogram(vecPrefAnglesNOT,vecBinEdgesPrefOri);
title('Preferred angle, NOT')
fixfig

%% plot bandwidth and mean tuning curves
%FWHM=2*arccos(1- [(1/kappa) * ln(2)] )
%mean tuning curves
figure;
vecOri = cellRawOri{1};
matRawResp = cell2mat(cellRawMean');
dblRadCenterTo = pi;
intCenterTo = find(vecPlotTheta>=dblRadCenterTo,1);
for intArea=1:numel(cellUniqueAreas)
	strArea = cellUniqueAreas{intArea};
matRaw = matRawResp(ismember(cellArea,strArea),:);
[vecMaxVal,vecMaxIdx] = max(matRaw,[],2);
matNorm = bsxfun(@rdivide,matRaw,vecMaxVal);
for intN=1:size(matNorm,1)
	matNorm(intN,:) = circshift(matNorm(intN,:),[0 intCenterTo-vecMaxIdx(intN)]);
end

%plot
subplot(2,3,intArea+3)
hold on
plot(vecOri-dblRadCenterTo,matNorm,'color',[0.7 0.7 0.7]);
errorbar(vecOri-dblRadCenterTo,mean(matNorm,1),std(matNorm,[],1)/sqrt(size(matNorm,1)),'b');
hold off
ylabel('Normalized response (1/max)');
xlabel('Radians');
fixfig
end

%BW
vecPeaksBW= real(vecPeaksBW);
dblStep = pi/6;
vecBinEdgesBW = 0:dblStep:(2*pi);
vecBinCentersBW = vecBinEdgesBW(2:end) - dblStep/2;
vecBW_NOT = histcounts(vecPeaksBW(1,ismember(cellArea,'NOT')),vecBinEdgesBW);
vecBW_SC = histcounts(vecPeaksBW(1,ismember(cellArea,'SC')),vecBinEdgesBW);
vecBW_V1 = histcounts(vecPeaksBW(1,ismember(cellArea,'V1')),vecBinEdgesBW);


%plot
subplot(2,3,1)
bar(vecBinCentersBW,vecBW_V1/intN_V1)
title(sprintf('V1; %d cells',intN_V1));
xlabel(sprintf('Bandwidth (rad)'));
ylabel('Fraction of V1 cells');
ylim(vecLimY);
fixfig

%SC
subplot(2,3,2)
bar(vecBinCentersBW,vecBW_SC/intN_SC)
title(sprintf('SC; %d cells',intN_SC));
xlabel(sprintf('Bandwidth (rad)'));
ylabel('Fraction of SC cells');
ylim(vecLimY);
fixfig

%NOT
subplot(2,3,3)
bar(vecBinCentersBW,vecBW_NOT/intN_NOT)
title(sprintf('NOT; %d cells',intN_NOT));
xlabel(sprintf('Bandwidth (rad)'));
ylabel('Fraction of NOT cells');
ylim(vecLimY);
fixfig


%% save