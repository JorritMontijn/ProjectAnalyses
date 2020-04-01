%get block data
close all
strOutputDir = 'D:\Data\Results\stimdetection\raw_movie\';
strMasterPath = 'D:\Data\Processed\imagingdata\';
strSession = '20140530';
intRecording = 1;
strRecording = sprintf('xyt%02d',intRecording);
strRecPath = [strMasterPath strSession filesep strRecording filesep];
strFilename = sprintf('%s%s_ses.mat',strSession,strRecording);
load(['D:\Data\Results\stimdetection\dataPreProAggregate' strSession '.mat']);

vecTrials = 1:48;
%vecFrames = 6701:10500;
vecFrames = 1:10000;
dblSamplingFreq = 25.4;
vecTimeSecs = vecFrames/dblSamplingFreq;
intFrames = length(vecFrames);

%get reaction times
vecRTsPerC = nan(1,6);
vecC = [0 0.5 2 8 32 100]/100;
vecAllRTs = cellMultiSes{1}.structStim.vecTrialRespSecs - cellMultiSes{1}.structStim.SecsOn; %get RTs for hit trials and selected contrasts
for intC=1:6
	vecTheseTrials = cellMultiSes{1}.structStim.Contrast == vecC(intC);
	vecRTsPerC(intC) = nanmedian(vecAllRTs(vecTheseTrials));
end

%%
%pre-allocate
vecTrialContrastsHit = [];
vecTrialContrastsMiss = [];
vecTrialActHit = [];
vecTrialActMiss = [];
vecTrialHetHit = [];
vecTrialHetMiss = [];

%perform analysis
h=figure;
set(gcf,'Color',[0 0 0])
drawnow;
jFig = get(handle(gcf), 'JavaFrame');
jFig.setMaximized(true);
figure(gcf);

dblMarkerSize = 30;
dblRangeT = 8; %secs
intFrameRateStimulus = 60; %hz

%prepare neural data
[indTuned,vecNeuronPrefStim,vecOSI,cellAggregate] = getTunedStimDetectionNeurons(cellMultiSes{1});
%remove non-tuned neurons
cellMultiSes{1}.neuron(~indTuned) = [];
vecNeuronPrefStim(~indTuned) = [];
intNeurons = numel(cellMultiSes{1}.neuron);
matAct = zeros(intNeurons,12000);
for intNeuron=1:intNeurons
	matAct(intNeuron,:) = cellMultiSes{1}.neuron(intNeuron).dFoF(1:12000);
end
[vecHeterogeneity,vecActivity] = calcMatRespHeteroGen(matAct);
vecHeterogeneity(1:5) = nan;
vecActivity(1:5) = nan;
vecHetLimY = [0.6 1.6+eps];
dblHetLowStimBound = 1.4;
vecActLimY = [-0.06 0.1+eps];
dblActLowStimBound = 0.07;

%create vector stimulus features
vecOri = nan(1,intFrames);
vecContrast = ones(1,intFrames)*-1;
vecLicking = cell2mat(cellMultiSes{1}.structStim.cellRespPulses);
vecLicking(isnan(vecLicking)) = [];
indLicking = false(1,intFrames);
indLicking(vecLicking-1) = true;
for intTrial=vecTrials
	vecStimFrames = (cellMultiSes{1}.structStim.FrameOn(intTrial):cellMultiSes{1}.structStim.FrameOff(intTrial))-vecFrames(1);
	intOri = cellMultiSes{1}.structStim.Orientation(intTrial);
	intC = cellMultiSes{1}.structStim.Contrast(intTrial);
	vecContrast(vecStimFrames) = intC;
	vecOri(vecStimFrames) = intOri;
end
indStimulus = [false(1,vecFrames(1)-1) vecContrast>-1 false(1,500)];
boolStimulus = false;

drawnow;
intStimFrame = 1;
matGrating = zeros(1024,1024,1,'uint8')+128;
intSwitchTrialActive = 0;


%% run loop
for intT=1:intFrames
	%
	%get raw data
	intFrame = vecFrames(intT);
	if mod(intT,10) == 0,fprintf('Creating frame %d [%d/%d] [%s]\n',intFrame,intT,intFrames,getTime);end
	matRGB = zeros(512,512,3,'uint8');
	matRGB(:,:,1) = imread([strRecPath 'images' filesep sprintf('t%05d',intFrame-1) '_ch00.tif']);
	matRGB(:,:,2) = imread([strRecPath 'images' filesep sprintf('t%05d',intFrame-1) '_ch01.tif']);
	
	%get stimulus+response data
	if vecContrast(intT)> -1 && ~boolStimulus %start stim
		boolStimulus = true;
		intOri = vecOri(intT);
		dblContrast = vecContrast(intT);
		dblPhaseRand = rand(1);
		dblStartT = vecTimeSecs(intT);
		intStimFrame = ceil(dblPhaseRand * intFrameRateStimulus);
		
		%get grating file
		strGratingFile = sprintf('%sgratingmovie_SQUARE_%dfps_Ori%03d_Speed1.mat','D:\Acquisition\LaserComputer\experiment_scripts\GratingMovies\',60,intOri);
		sLoad = load(strGratingFile);
		matGrating = imcontrast(sLoad.matGrating,dblContrast);
	elseif vecContrast(intT)> -1 && boolStimulus %continue stim
		dblSecsT = vecTimeSecs(intT);
		tStamp = mod(dblSecsT-dblStartT+dblPhaseRand,1);
		intStimFrame = ceil(tStamp * intFrameRateStimulus);
	elseif vecContrast(intT) == -1 && boolStimulus %stop stim
		intStimFrame = 1;
		matGrating = zeros(1024,1024,1,'uint8')+128;
		boolStimulus = false;
	end
	
	%get neural data
	intPreFrame = intFrame-round(dblRangeT*dblSamplingFreq);
	vecUseFrames = intPreFrame:intFrame;
	vecUseFrames(vecUseFrames<1) = 1;
	intPlotFrames = length(vecUseFrames);
	vecThisAct = vecActivity(vecUseFrames);
	vecThisHet = vecHeterogeneity(vecUseFrames);
	vecPreLicking = indLicking(vecUseFrames);
	vecPlotT = linspace(0,dblRangeT,length(vecThisAct));
	
	%get stimulus stretches
	indPlotStim = indStimulus(vecUseFrames);
	intStart = find(indPlotStim,1,'first')-1;
	intStop = find(indPlotStim,1,'last')-1;
	vecPlotStim = intStart:intStop;
	
	%display raw movie frame
	clf
	subplot(1,2,1)
	imshow(matRGB);
	title(sprintf('Mouse 5; t=%03.3fs [frame %d]',intFrame/dblSamplingFreq,intFrame),'FontSize',16,'Color','w');
	
	%display stimulus on screen
	subplot(7,4,3)
	imshow(matGrating(:,:,intStimFrame));
	title('Screen','FontSize',12,'Color','w');
	
	%display behavior
	subplot(7,4,4)
	set(gca,'color','none')
	if any(indLicking(max([1 (intFrame-10)]):intFrame))
		text(0,0.5,'Response: Licking','FontSize',20,'Color','r');
	else
		text(0,0.5,'Response: No licking','FontSize',20,'Color','w');
	end
	
	%display mean dF/F0 so far
	subplot(7,4,7)
	hold on
	set(gca,'color','none')
	if ~isempty(vecTrialContrastsHit)
		scatter(vecTrialContrastsHit,vecTrialActHit,dblMarkerSize,'go','filled');
	end
	if ~isempty(vecTrialContrastsMiss)
		scatter(vecTrialContrastsMiss,vecTrialActMiss,dblMarkerSize,'ro','filled');
	end
	set(gca,'xscale','log');
	grid on
	xlim([0.15 110]);
	ylim([-0.02 0.06+eps]);
	set(gca,'xtick',[0.2 0.5 2 8 32 100],'xticklabel',[0 0.5 2 8 32 100]);
	set(gca,'color','none','xcolor','w','ycolor','w ');
	xlabel('Stimulus Contrast (%)','FontSize',12,'Color','w');
	ylabel('dF/F0','FontSize',12,'Color','w');
	set(gca,'FontSize',12,'Linewidth',1); %set grid line width and change font size of x/y ticks
	
	%display heterogeneity so far
	subplot(7,4,8)
	hold on
	set(gca,'color','none')
	if ~isempty(vecTrialContrastsHit)
		scatter(vecTrialContrastsHit,vecTrialHetHit,dblMarkerSize,'go','filled');
	end
	if ~isempty(vecTrialContrastsMiss)
		scatter(vecTrialContrastsMiss,vecTrialHetMiss,dblMarkerSize,'ro','filled');
	end
	set(gca,'xscale','log');
	grid on
	xlim([0.15 110]);
	ylim([0.8 1.401]);
	set(gca,'xtick',[0.2 0.5 2 8 32 100],'xticklabel',[0 0.5 2 8 32 100]);
	set(gca,'color','none','xcolor','w','ycolor','w ');
	xlabel('Stimulus Contrast (%)','FontSize',12,'Color','w');
	ylabel('Heterogeneity','FontSize',12,'Color','w');
	set(gca,'FontSize',12,'Linewidth',1); %set grid line width and change font size of x/y ticks
	
	
	%create tick variables
	vecTickX = 0:1:dblRangeT;
	cellTickLabelX = mat2cell(num2str((vecTickX-dblRangeT)'),ones(1,length(vecTickX)),2)';
	cellTickLabelX{end} = [cellTickLabelX{end}(2) ' (now)'];
	intStimStartFrame = intT-(intPlotFrames-intStart)+1;
	
	%plot dF/F0
	subplot(7,2,8)
	plot(vecPlotT,vecThisAct,'Color',[1 1 1]*0.6,'linewidth',2);
	hold on;
	if ~isempty(intStart)
		patch([intStart intStop intStop intStart]/dblSamplingFreq,[dblActLowStimBound dblActLowStimBound vecActLimY(2) vecActLimY(2)],[0.5 0.5 0.5]);
		plot(vecPlotStim/dblSamplingFreq,vecThisAct(vecPlotStim+1),'Color',[0.2 0.2 1],'linewidth',3);
		if intStart < (intPlotFrames-20) && intStop > 20
			text((intStart+2)/dblSamplingFreq,mean([dblActLowStimBound vecActLimY(2)]),sprintf('%.1f%%',vecContrast(intStimStartFrame)*100),'FontSize',12,'Color','w');
		end
	end
	for intLick=find(vecPreLicking)
		plot([1 1]*(intLick/dblSamplingFreq),vecActLimY,'--','Color',[1 0 0],'linewidth',2);
	end
	hold off;
	title('Mean population dF/F0','FontSize',12,'Color','w');
	xlabel('Time (s)','FontSize',12,'Color','w');
	ylabel('dF/F0','FontSize',12,'Color','w');
	set(gca,'color','none','xcolor','w','ycolor','w ');
	set(gca,'xtick',vecTickX,'xticklabel',cellTickLabelX)
	set(gca,'FontSize',12,'Linewidth',2); %set grid line width and change font size of x/y ticks
	grid on
	xlim([0 dblRangeT+eps])
	ylim(vecActLimY);
	
	%plot heterogeneity
	subplot(7,2,12)
	plot(vecPlotT,vecThisHet,'Color',[1 1 1]*0.6,'linewidth',2);
	hold on;
	if ~isempty(intStart)
		patch([intStart intStop intStop intStart]/dblSamplingFreq,[dblHetLowStimBound dblHetLowStimBound vecHetLimY(2) vecHetLimY(2)],[0.5 0.5 0.5]);
		plot(vecPlotStim/dblSamplingFreq,vecThisHet(vecPlotStim+1),'Color',[1 0.2 0.2],'linewidth',3);
		if intStart < (intPlotFrames-20) && intStop > 20
			text((intStart+2)/dblSamplingFreq,mean([dblHetLowStimBound vecHetLimY(2)]),sprintf('%.1f%%',vecContrast(intStimStartFrame)*100),'FontSize',12,'Color','w');
		end
	end
	for intLick=find(vecPreLicking)
		plot([1 1]*(intLick/dblSamplingFreq),vecHetLimY,'--','Color',[1 0 0],'linewidth',2);
	end
	hold off;
	title('Mean population heterogeneity','FontSize',12,'Color','w');
	xlabel('Time (s)','FontSize',12,'Color','w');
	ylabel('Heterogeneity','FontSize',12,'Color','w');
	set(gca,'color','none','xcolor','w','ycolor','w ');
	set(gca,'xtick',vecTickX,'xticklabel',cellTickLabelX)
	set(gca,'FontSize',12,'Linewidth',2); %set grid line width and change font size of x/y ticks
	grid on
	xlim([0 dblRangeT+eps])
	ylim(vecHetLimY);
	
	%get trial data
	if intSwitchTrialActive == 1 && ~boolStimulus
		intSwitchTrialActive = 0;
		
		dblStimContrast = vecContrast(intStimStartFrame)*100;
		
		if vecPreLicking(intStop+1) == 1 %hit
			dblHet = mean(vecThisHet(vecPlotStim+1));
			dblAct = mean(vecThisAct(vecPlotStim+1));
			
			if dblStimContrast == 0,dblStimContrast=0.2;end
			vecTrialContrastsHit(end+1) = dblStimContrast;
			vecTrialActHit(end+1) = dblAct;
			vecTrialHetHit(end+1) = dblHet;
		else %miss
			[dummy,intC]=find(vecC*100 == dblStimContrast,1,'first');
			dblMedianRT = vecRTsPerC(intC);
			
			intSelectNrFrames = round(dblMedianRT*dblSamplingFreq);
			
			dblHet = mean(vecThisHet(vecPlotStim(1:intSelectNrFrames)+1));
			dblAct = mean(vecThisAct(vecPlotStim(1:intSelectNrFrames)+1));
			
			if dblStimContrast == 0,dblStimContrast=0.2;end
			vecTrialContrastsMiss(end+1) = dblStimContrast;
			vecTrialActMiss(end+1) = dblAct;
			vecTrialHetMiss(end+1) = dblHet;
		end
	elseif intSwitchTrialActive == 0 && boolStimulus
		intSwitchTrialActive = 1;
	end
	
	strFilenameMovie = sprintf('frame%04d.tif',intT);
	drawnow;
	export_fig([strOutputDir strFilenameMovie],'-nocrop');
end