%make pop coding fig

%{
load('D:\Data\Processed\imagingdata\20140530\xyt01\20140530xyt01_ses.mat')
ses = doRecalcdFoF(ses,3,[],'astrocyte');
ses = doRecalcdFoF(ses,3,[],'neuron');
%}

%get stim info
%ses.structStim.Orientation = mod(ses.structStim.Orientation,180);
vecStimTypes = getUniqueVals(ses.structStim.Orientation);
intStimTypes = length(vecStimTypes);
cmap = colormap(jet(intStimTypes));

%set plotting vars
intStart = 12750;
intStop = 16000;
dblMinY = -0.05;
dblMaxY = 0.25;
vecKernel = normpdf(-2:0.3:2,0,1);
vecKernel = vecKernel/sum(vecKernel);

%{
%plot heterogeneity
[vecHeterogeneity,vecActivity] = calcSlidingHeteroGen(ses);

figure
hold on
for intStim=1:numel(ses.structStim.FrameOn)
	intStimStart = ses.structStim.FrameOn(intStim);
	intStimStop = ses.structStim.FrameOff(intStim);
	dblOri = ses.structStim.Orientation(intStim);
	intStimType = find(vecStimTypes==dblOri,1);
	fill([intStimStart intStimStart intStimStop intStimStop],[0.8 1.6 1.6 0.8],cmap(intStimType,:),'EdgeColor',cmap(intStimType,:));
	
end
vecHeterogeneity = conv(vecHeterogeneity,vecKernel,'same');
plot(vecHeterogeneity,'k-','LineWidth',2)
hold off

set(gca,'XTick',ses.structStim.FrameOn)
set(gca,'XTickLabel',roundi(ses.structStim.FrameOn/ses.samplingFreq,1))
grid on
xlim([intStart intStop])
title(sprintf('heterogeneity'))
drawnow
export_fig(sprintf('heterogeneity.tif'))
export_fig(sprintf('heterogeneity.pdf'))
pause
close all
%}
%{
%plot mean act
figure
hold on
for intStim=1:numel(ses.structStim.FrameOn)
	intStimStart = ses.structStim.FrameOn(intStim);
	intStimStop = ses.structStim.FrameOff(intStim);
	dblOri = ses.structStim.Orientation(intStim);
	intStimType = find(vecStimTypes==dblOri,1);
	fill([intStimStart intStimStart intStimStop intStimStop],[dblMinY dblMaxY dblMaxY dblMinY],cmap(intStimType,:),'EdgeColor',cmap(intStimType,:));
	
end
vecActivity = conv(vecActivity,vecKernel,'same');
plot(vecActivity,'b-','LineWidth',2)
hold off
ylim([dblMinY dblMaxY])
set(gca,'XTick',ses.structStim.FrameOn)
set(gca,'XTickLabel',roundi(ses.structStim.FrameOn/ses.samplingFreq,1))
grid on
xlim([intStart intStop])
title(sprintf('mean df/f'))
drawnow
export_fig(sprintf('meanPopAct.tif'))
export_fig(sprintf('meanPopAct.pdf'))
pause
close all
%}




%plot licks
vecLickFramesOn = [];
for intStim=1:length(ses.structStim.cellRespPulses)
	vecLickFramesOn = [vecLickFramesOn ses.structStim.cellRespPulses{intStim}];
end
vecDetections = ses.structStim.vecTrialRespPulses(~isnan(ses.structStim.vecTrialRespPulses));
vecLickFramesOn(isnan(vecLickFramesOn)) = [];
vecDetectTrace = zeros(size(ses.neuron(1).dFoF));
vecDetectTrace(vecDetections) = 1;
vecLickTrace = zeros(size(ses.neuron(1).dFoF));
vecLickTrace(vecLickFramesOn) = 1;
figure
hold on
for intStim=1:numel(ses.structStim.FrameOn)
	intStimStart = ses.structStim.FrameOn(intStim);
	intStimStop = ses.structStim.FrameOff(intStim);
	dblOri = ses.structStim.Orientation(intStim);
	intStimType = find(vecStimTypes==dblOri,1);
	fill([intStimStart intStimStart intStimStop intStimStop],[dblMinY dblMaxY dblMaxY dblMinY],cmap(intStimType,:),'EdgeColor',cmap(intStimType,:));
	
end

stairs(vecLickTrace,'k-','LineWidth',2)
stairs(vecDetectTrace,'g-','LineWidth',2)
hold off

set(gca,'XTick',ses.structStim.FrameOn)
set(gca,'XTickLabel',roundi(ses.structStim.FrameOn/ses.samplingFreq,1))
grid on
xlim([intStart intStop])
title(sprintf('licks'))
drawnow
export_fig(sprintf('licks.tif'))
export_fig(sprintf('licks.pdf'))
pause
close all


37-46
%return
%{
%plot moves
vecMoveFramesOn = [];
for intStim=1:length(ses.structStim.cellMovePulses)
	vecMoveFramesOn = [vecMoveFramesOn ses.structStim.cellMovePulses{intStim}];
end
vecMoveTrace = zeros(size(ses.neuron(1).dFoF));
vecMoveTrace(vecMoveFramesOn) = 1;
figure
hold on
for intStim=1:numel(ses.structStim.FrameOn)
	intStimStart = ses.structStim.FrameOn(intStim);
	intStimStop = ses.structStim.FrameOff(intStim);
	dblOri = ses.structStim.Orientation(intStim);
	intStimType = find(vecStimTypes==dblOri,1);
	fill([intStimStart intStimStart intStimStop intStimStop],[dblMinY dblMaxY dblMaxY dblMinY],cmap(intStimType,:),'EdgeColor',cmap(intStimType,:));
	
end

stairs(vecMoveTrace,'k-','LineWidth',2)
hold off

set(gca,'XTick',ses.structStim.FrameOn)
set(gca,'XTickLabel',roundi(ses.structStim.FrameOn/ses.samplingFreq,1))
grid on
xlim([intStart intStop])
title(sprintf('movements'))
drawnow
export_fig(sprintf('moves.tif'))
export_fig(sprintf('moves.pdf'))
pause
close all
%}
%{
%% plot astrocytes
for intAstro=1:numel(ses.astrocyte)%
	vecTrace = ses.astrocyte(intAstro).dFoF;
	vecTrace = conv(vecTrace,vecKernel,'same');
	
	figure
	hold on
	for intStim=1:numel(ses.structStim.FrameOn)
		intStimStart = ses.structStim.FrameOn(intStim);
		intStimStop = ses.structStim.FrameOff(intStim);
		dblOri = ses.structStim.Orientation(intStim);
		intStimType = find(vecStimTypes==dblOri,1);
		fill([intStimStart intStimStart intStimStop intStimStop],[dblMinY dblMaxY dblMaxY dblMinY],cmap(intStimType,:),'EdgeColor',cmap(intStimType,:));
		
	end
	
	plot(vecTrace,'k-','LineWidth',2)
	hold off
	set(gca,'XTick',ses.structStim.FrameOn)
	set(gca,'XTickLabel',roundi(ses.structStim.FrameOn/ses.samplingFreq,1))
	grid on
	xlim([intStart intStop])
	ylim([dblMinY dblMaxY])
	title(sprintf('astrocyte  %d',intAstro))
	drawnow
	export_fig(sprintf('astro%d.tif',intAstro))
	export_fig(sprintf('astro%d.pdf',intAstro))
	pause
	close
end

%8 32 33 37 72 75 96? 158 183
%[8 32 33 37 72 75 158]
%21 45 81
%}

%%
intNeurons = numel(ses.neuron);
for intNeuron=[2 3 4 6 29 40 50 79 83 112 117 123]
	vecTrace = ses.neuron(intNeuron).dFoF;
	vecTrace = conv(vecTrace,vecKernel,'same');

	figure
	hold on
	for intStim=1:numel(ses.structStim.FrameOn)
		intStimStart = ses.structStim.FrameOn(intStim);
		intStimStop = ses.structStim.FrameOff(intStim);
		dblOri = ses.structStim.Orientation(intStim);
		intStimType = find(vecStimTypes==dblOri,1);
		fill([intStimStart intStimStart intStimStop intStimStop],[dblMinY dblMaxY dblMaxY dblMinY],cmap(intStimType,:),'EdgeColor',cmap(intStimType,:));
		
	end
	
	plot(vecTrace,'k-','LineWidth',2)
	hold off
	set(gca,'XTick',ses.structStim.FrameOn)
	set(gca,'XTickLabel',roundi(ses.structStim.FrameOn/ses.samplingFreq,1))
	grid on
	xlim([intStart intStop])
	ylim([dblMinY*2 dblMaxY*2])
	title(sprintf('neuron %d',intNeuron))
	drawnow
	export_fig(sprintf('neuron%d.tif',intNeuron))
	export_fig(sprintf('neuron%d.pdf',intNeuron))
	pause
	close
end
%{
PV:
1, 3

neuron:
2[/], 4?, 6, 8, 12, 15!, 17, 19
%}