%% set recording
clear all;%close all;
strSelectMouseType = 'WT';
strSelectArea = 'V1';
strSelectStim = 'DG'; %DG=drifting grating, RF=receptive field mapping, RF-DG=merged block
intSelectPopulation	= 3;
strFigPath = 'D:\Data\ResultsNystagmus\CSD\';
boolSavePlots = true;

%% load data
loadNystagmusData;

%% run header
runNystagmusHeader;

%% remove sample outliers
dblSampleOutlierSD = 6;
matTempLFP = matLFPData;
vecSD = std(matTempLFP,[],2);
vecMean = mean(matTempLFP,2);
matLFPZ = bsxfun(@rdivide,bsxfun(@minus,matTempLFP,vecMean),vecSD);
matLFPZ(:,any(abs(matLFPZ)>dblSampleOutlierSD,1)) = nan;

%% build trial LFP & MUA
sOpt = struct;
sOpt.handleFig = -1;
sOpt.vecWindow = [0 0.5];
intChannels= 32;
intEvents = numel(vecStimOnTime);
intTimepointsLFP = dblResampFreq*range(sOpt.vecWindow)+1;
matTrialLFP = nan(intChannels,intTimepointsLFP,intEvents);
intTimepointsEnv = ceil((1/mean(diff(vecEnvTimestamps)))*range(sOpt.vecWindow)+1);
matTrialEnv = nan(intChannels,intTimepointsEnv,intEvents);
for intCh=1:32
	%LFP
	[vecTrialLFP,vecSEM,vecWindowBinCenters,matPET_LFP] = doPEP(vecLFPTimestamps,matLFPZ(intCh,:),vecStimOnTime,sOpt);
	matTrialLFP(intCh,:,:) = matPET_LFP';
	%MUA
	[vecTrialEnv,vecSEM,vecWindowBinCenters,matPET_Env] = doPEP(vecEnvTimestamps,matEnvData(intCh,:),vecStimOnTime,sOpt);
	matTrialEnv(intCh,:,:) = matPET_Env';
end

%% remove trial outliers
dblTrialOutlierSD = 3;
matAvTrialLFP = squeeze(nanmean(matTrialLFP,2));
vecSD = std(matAvTrialLFP,[],2);
vecMean = mean(matAvTrialLFP,2);
matTrialLFPZ = bsxfun(@rdivide,bsxfun(@minus,matAvTrialLFP,vecMean),vecSD);
indRemTrials = (sum(abs(matTrialLFPZ)>dblTrialOutlierSD,1)>1);
matTrialLFP = matTrialLFP(:,:,~indRemTrials);
matMeanTrialLFP = nanmean(matTrialLFP,3);
matTrialEnv = matTrialEnv(:,:,~indRemTrials);
matMeanTrialEnv = nanmean(matTrialEnv,3);

%% build CSD
close all;
cellStr = {'Odd','Even'};
for intOddEven = 1:2
	figure
	%get data
	matUseLFP = nanmean(matMeanTrialLFP(intOddEven:2:end,:),3);
	matUseMUA = nanmean(matMeanTrialEnv(intOddEven:2:end,:),3);

% remove channels
vecChRem = [11 14];
for intChRem=vecChRem
	matUseLFP(intChRem,:) = (matUseLFP(intChRem-1,:) + matUseLFP(intChRem+1,:))/2;
end
%plot MUA
subplot(2,2,1)
imagesc(vecWindowBinCenters,vecChannelDepth((intOddEven):2:(end)),matUseMUA);colormap(hot);freezeColors;
xlabel('Time after stim onset (s)');
ylabel('Approx. depth (micron)');
fixfig;grid off;
title(sprintf('MUA %s, %sB%s %s',strSelectArea,strDate,strBlock,cellStr{intOddEven}));


%plot LFP
subplot(2,2,3)
imagesc(vecWindowBinCenters,vecChannelDepth((intOddEven):2:(end)),matUseLFP);colormap(parula);freezeColors;
xlabel('Time after stim onset (s)');
ylabel('Approx. depth (micron)');fixfig;grid off;title('LFP');

%loop
c = 0;
sep = 2;
matTrialCSD = nan(intChannels/2-sep*2,intTimepointsLFP);
intStart = intOddEven+sep;
intEnd = intChannels-sep;
vecUseChannels = (1+sep):(size(matUseLFP,1)-sep);
for ch = vecUseChannels
	c = c+1;
	matTempCSD = -0.4.*(matUseLFP(ch-sep,:)-(2.*matUseLFP(ch,:))+matUseLFP(ch+sep,:))./(((25.*10.^-6).*sep).^2);
	matTrialCSD(c,:) = nanmean(matTempCSD,3);
end
%smooth
vecFiltX = normpdf(-2:2,0,1);
vecFiltY = normpdf(-2:2,0,1);
matFilt = vecFiltY' * vecFiltX;
matFilt = matFilt ./ sum(matFilt(:));
dblMean = nanmean(matTrialCSD(:));
matTrialCSD(isnan(matTrialCSD))=dblMean;
matTrialCSD = matTrialCSD - dblMean;
matTrialCSD = conv2(matTrialCSD,matFilt,'same');
matTrialCSD = matTrialCSD + dblMean;
subplot(2,2,4)
imagesc(vecWindowBinCenters,vecChannelDepth(vecUseChannels*2),matTrialCSD);colormap(parula);freezeColors;
xlabel('Time after stim onset (s)');
ylabel('Approx. depth (micron)');fixfig;grid off;title('CSD');

%plot xcorr
subplot(2,2,2)
matCorr=corr(matLFPData');
vecUseCh = intOddEven:2:intChannels;
imagesc(vecUseCh,vecUseCh,matCorr(vecUseCh,vecUseCh),[-1 1]);colormap(redblue);freezeColors;colorbar
xlabel('Channel #');
ylabel('Channel #');fixfig;grid off;title('x-corr');

% save plot
if boolSavePlots
	drawnow;
	strFileName = sprintf('%s%sB%s_%s',strSelectArea,strDate,strBlock,cellStr{intOddEven});
	%new maximized figure
	
	drawnow;
	jFig = get(handle(gcf), 'JavaFrame');
	jFig.setMaximized(true);
	figure(gcf);
	drawnow;
	
	export_fig([strFigPath strFileName 'CSD.tif']);
	export_fig([strFigPath strFileName 'CSD.pdf']);
	
	pause(0.1);
end
end