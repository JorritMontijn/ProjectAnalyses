% Script to plot nice neurons with opto:

% ChannelSel = 7; %PV ChR2 neuron
% ChannelSel = 20; %Ìnhibited multi-unit activity

%% Get Data
% [Data] = MOL_SelectData();

%% Get arguments:
sessionData     = Data.sessionData;
trialData       = Data.trialData;
spikeData       = Data.spikeData;

RawDataDir      = 'F:\Data\PMAL\RawData\P3.16\2017-12-13_11-54-15\RawData';

%Save figure as given by user input
OutputDir       = 'E:\Documents\PhD\Figures\Project PMAL\Results\';

% Alpha = 1;
Alpha = 0.3;

TimeWindow = [sessionData.t_start sessionData.t_stop];

%% parameters for plotting example piece:
exampletrials           = [1 6];
example_tstart          = trialData.photostimStart(exampletrials(1))-1e6;
example_tstop           = trialData.photostimStart(exampletrials(2))+2e6;

%High pass filter:
hp_butter           = 600;  %Low pass filter (Hz) (Butter)
ord_butter          = 4;    %Butterworth filter order
[B_butt,A_butt]     = butter(ord_butter,hp_butter/(32000/2),'high');

%% First for the channel that is excited:

%Get the data:
lfpData             = MOL_extractLFP(RawDataDir,7,TimeWindow);
lfpData.hpsignal    = filtfilt(B_butt,A_butt,lfpData.signal);

%% Plot for a few trials:
figure; 
set(gcf,'units','normalized','Position',[0.2 0.1 0.5 0.7],'color','w')
selectedpiece           = lfpData.ts>example_tstart & lfpData.ts<example_tstop;
plot(lfpData.ts(selectedpiece),lfpData.hpsignal(selectedpiece),'k'); hold all;
ymax    = max(lfpData.hpsignal(selectedpiece))*1.1;
ymin    = min(lfpData.hpsignal(selectedpiece))*1.1;

for iT = exampletrials(1):exampletrials(2)
    optostart   = trialData.photostimStart(iT);
    optoend     = trialData.photostimEnd(iT);
    X = [optostart optoend optoend optostart];
    Y = [ymin ymin ymax ymax];
    patch(X, Y,'blue','EdgeColor','blue','FaceAlpha',Alpha,'EdgeAlpha',0); hold on;
end

ylim([ymin*1.2 ymax*1.2]);
xlim([example_tstart example_tstop]);
print(gcf,fullfile(OutputDir,strcat('Ch7_Raw_Trace','.pdf')),'-dpdf','-bestfit');
% export_fig(fullfile(OutputDir,strcat('Ch7_Raw_Trace','.pdf')),'-eps');

%%
figure;
set(gcf,'units','normalized','Position',[0.2 0.1 0.5 0.7],'color','w')

preposttime         = 0.5e6;
exampletrial        = 2;
selectedpiece       = lfpData.ts>trialData.photostimStart(exampletrial)-preposttime & lfpData.ts<trialData.photostimEnd(exampletrial)+preposttime;

plot(lfpData.ts(selectedpiece),lfpData.hpsignal(selectedpiece),'k'); hold all;
ymax = max(lfpData.hpsignal(selectedpiece))*1.1;
ymin = min(lfpData.hpsignal(selectedpiece))*1.1;

for i = 1:1000/50
    pulsestart = trialData.photostimStart(exampletrial)+50e3*(i-1)-8e3;
    X = [pulsestart pulsestart+10e3 pulsestart+10e3 pulsestart];
    Y = [ymin ymin ymax ymax];
    patch(X, Y,'blue','EdgeColor','blue','FaceAlpha',Alpha,'EdgeAlpha',0); hold on;
end

ylim([ymin*1.1 ymax*1.1]);
set(gca,'linewidth',3)
print(gcf,fullfile(OutputDir,strcat('Ch7_Pulse_Closeup','.pdf')),'-dpdf','-bestfit');

%% Plot the average waveform of the neuron:

for iSpike=1:size(spikeData.wf{7},2)
    spikeData.wf{7}(:,iSpike) = spikeData.wf{7}(:,iSpike)-repmat(mean(spikeData.wf{7}(1:10,iSpike)),64,1);
end

figure;
set(gcf,'units','normalized','Position',[0.2 0.1 0.5 0.7],'color','w')
shadedErrorBar(1:64,mean(spikeData.wf{7},2),std(spikeData.wf{7},[],2))
set(gca,'linewidth',3)
print(gcf,fullfile(OutputDir,strcat('Ch7_Spike_Waveform','.pdf')),'-dpdf','-bestfit');

%% Get the data for the multiUnit figures:
lfpData             = MOL_extractLFP(RawDataDir,20,TimeWindow);
lfpData.hpsignal    = filtfilt(B_butt,A_butt,lfpData.signal);

%% Plot for a few trials:
figure;
set(gcf,'units','normalized','Position',[0.2 0.1 0.5 0.7],'color','w')
selectedpiece           = lfpData.ts>example_tstart & lfpData.ts<example_tstop;

plot(lfpData.ts(selectedpiece),lfpData.hpsignal(selectedpiece),'k'); hold all;
ymax    = max(lfpData.hpsignal(selectedpiece))*1.1;
ymin    = min(lfpData.hpsignal(selectedpiece))*1.1;

for iT = exampletrials(1):exampletrials(2)
    optostart   = trialData.photostimStart(iT);
    optoend     = trialData.photostimEnd(iT);
    X = [optostart optoend optoend optostart];
    Y = [ymin ymin ymax ymax];
    patch(X, Y,'blue','EdgeColor','blue','FaceAlpha',Alpha,'EdgeAlpha',0); hold on;
end

ylim([ymin*1.2 ymax*1.2]);
xlim([example_tstart example_tstop]);
set(gca,'linewidth',3)

print(gcf,fullfile(OutputDir,strcat('Ch20_Raw_Trace','.pdf')),'-dpdf','-bestfit');

%% plot MUA over time:
nthresholds = 4;

idx                 = find(abs(lfpData.hpsignal)>nthresholds*std(lfpData.hpsignal));
idx                 = idx(diff(idx)~=1);

params              = params_histresponse;
params.conv_sigma           = 0.1e6;        %sd of gaussian window for smoothing

edges               = lfpData.t_start:params.binsize:lfpData.t_end;
hist_mat            = histc(lfpData.ts(idx),edges) * 1e6/params.binsize;

N                   = params.conv_twin/params.binsize;
alpha               = ((N-1)/(params.conv_sigma/params.binsize))/2; %? = (N – 1)/(2?)
win                 = gausswin(N,alpha); %convolution with gaussian
win                 = win/sum(win); %normalized

%Smooth either the total or the individual trials:
hist_mat            = padarray(hist_mat,[0 round(length(win)/2)],'symmetric','both'); %pad the array on both sides for convolution
hist_mat            = conv(hist_mat,win,'valid'); %Take only the valid overlapping center of convolution
hist_mat            = hist_mat(1:length(edges)); %slight correction to get same size (edge vs convolution)

figure;
set(gcf,'units','normalized','Position',[0.2 0.1 0.5 0.7],'color','w')
selectedpiece           = edges>example_tstart & edges<example_tstop;

plot(edges(selectedpiece),hist_mat(selectedpiece),'k'); hold all;
ymax    = max(hist_mat(selectedpiece))*1.1;
ymin    = min(hist_mat(selectedpiece))*1.1;

for iT = exampletrials(1):exampletrials(2)
    optostart   = trialData.photostimStart(iT);
    optoend     = trialData.photostimEnd(iT);
    X = [optostart optoend optoend optostart];
    Y = [ymin ymin ymax ymax];
    patch(X, Y,'blue','EdgeColor','blue','FaceAlpha',Alpha,'EdgeAlpha',0); hold on;
end

ylim([ymin*1.2 ymax*1.2]);
xlim([example_tstart example_tstop]);
set(gca,'linewidth',3)
print(gcf,fullfile(OutputDir,strcat('Ch20_MUA_Trace','.pdf')),'-dpdf','-bestfit');

%% 
% 
% figure;
% plot(randn(1,500));
% 
% 
% print(gcf,fullfile(OutputDir,strcat('Ch20_MUA_Trace','.pdf')),'-dpdf','-bestfit');
% 
% 
% export_fig(fullfile(OutputDir,strcat('Ch20_MUA_Trace','.eps')),'-eps');
% 
% 
% f = figure;
% plot(rand(4))
% % Specify that the image in the PDF should be 2x2 
% % My PaperUnits is inches
% f.PaperPosition(3:4) = [2 2];
% print('-dpdf',fullfile(OutputDir,strcat('Ch20_MUA_Trace','.pdf')))

