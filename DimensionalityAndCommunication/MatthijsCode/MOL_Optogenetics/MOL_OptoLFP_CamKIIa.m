% Script to plot modulation by CamKIIa neurons with opto:

% ChannelSel = ; % neuron
% ChannelSel = ; %Ìnhibited multi-unit activity

%% Get Data
[Data] = MOL_SelectData();

%% Get arguments:
sessionData     = Data.sessionData;
trialData       = Data.trialData;
% spikeData       = Data.spikeData;

RootDir     = 'F:\';
RawDataDir  = fullfile(RootDir,'Data','CHDET','RawData',sessionData.mousename{1},sessionData.Rec_datetime{1},'RawData');
% RawDataDir  = 'H:\Data\CHDET\RawData\2009\2018-08-24_11-09-42\RawData';

%Save figure as given by user input
% OutputDir = 'E:\Documents\PhD\Figures\Project CHDET\Optogenetics\CaMKIIa_excitation\';
OutputDir = 'E:\Documents\PhD\Figures\Project CHDET\AuNoiseReponses\';

Alpha = 0.3;

TimeWindow = [sessionData.t_start sessionData.t_stop];

%% parameters for plotting example piece:
exampletrials           = [1 4];
example_tstart          = trialData.photostimStart(exampletrials(1))-1e6;
example_tstop           = trialData.photostimStart(exampletrials(2))+2e6;

pulselength             = 50e3;

selectedchannel         = 60;

%High pass filter:
hp_butter           = 600;  %Low pass filter (Hz) (Butter)
ord_butter          = 4;    %Butterworth filter order
[B_butt,A_butt]     = butter(ord_butter,hp_butter/(32000/2),'high');

%% First for the channel that is excited:

%Get the data:
lfpData             = MOL_extractLFP(RawDataDir,selectedchannel,TimeWindow);
lfpData.hpsignal    = filtfilt(B_butt,A_butt,lfpData.signal);
lfpData.ts          = lfpData.t_start:(1/32000*1e6):lfpData.t_end;

% lfpData.hpsignal    = lfpData.signal;

%% Plot for a few trials:
figure; 
set(gcf,'units','normalized','Position',[0.2 0.1 0.5 0.7],'color','w')
selectedpiece           = lfpData.ts>example_tstart & lfpData.ts<example_tstop;
plot(lfpData.ts(selectedpiece),lfpData.hpsignal(selectedpiece),'k'); hold all;
ymax    = max(lfpData.hpsignal(selectedpiece))*1.1;
ymin    = min(lfpData.hpsignal(selectedpiece))*1.1;

% ymax = 3e-4;
% ymin = -3e-4;

for iT = exampletrials(1):exampletrials(2)
    optostart   = trialData.photostimStart(iT);
    optoend     = trialData.photostimEnd(iT);
    X = [optostart optoend optoend optostart];
    Y = [ymin ymin ymax ymax];
    patch(X, Y,'blue','EdgeColor','blue','FaceAlpha',Alpha,'EdgeAlpha',0); hold on;
end

ylim([ymin*1.2 ymax*1.2]);
xlim([example_tstart example_tstop]);
print(gcf,fullfile(OutputDir,strcat('Ch',num2str(selectedchannel),'_Raw_Trace_','.pdf')),'-dpdf','-bestfit');
export_fig(fullfile(OutputDir,strcat('Ch',num2str(selectedchannel),'_Raw_Trace_','.png')),'-png');
% print(gcf,fullfile(OutputDir,strcat('Ch',num2str(selectedchannel),'_Raw_Trace_',num2str(sessionData.Photostimpower),'mW','.pdf')),'-dpdf','-bestfit');
% export_fig(fullfile(OutputDir,strcat('Ch',num2str(selectedchannel),'_Raw_Trace_',num2str(sessionData.Photostimpower),'mW','.png')),'-png');

%%
figure;
set(gcf,'units','normalized','Position',[0.2 0.1 0.5 0.7],'color','w')

preposttime         = 0.5e6;
exampletrial        = 2;
selectedpiece       = lfpData.ts>trialData.photostimStart(exampletrial)-preposttime & lfpData.ts<trialData.photostimEnd(exampletrial)+preposttime;

plot(lfpData.ts(selectedpiece),lfpData.hpsignal(selectedpiece),'k'); hold all;
ymax = max(lfpData.hpsignal(selectedpiece))*1.1;
ymin = min(lfpData.hpsignal(selectedpiece))*1.1;

% ymax = 3e-4;
% ymin = -3e-4;

for i = 1:sessionData.PhotostimFreq
    pulsestart = trialData.photostimStart(exampletrial)+1e6/sessionData.PhotostimFreq*(i-1);
    X = [pulsestart pulsestart+pulselength pulsestart+pulselength pulsestart];
    Y = [ymin ymin ymax ymax];
    patch(X, Y,'blue','EdgeColor','blue','FaceAlpha',Alpha,'EdgeAlpha',0); hold on;
end

% xlim([example_tstart example_tstop]);
ylim([ymin*1.1 ymax*1.1]);
set(gca,'linewidth',3)
print(gcf,fullfile(OutputDir,strcat('Ch',num2str(selectedchannel),'_Pulse_Closeup_',num2str(sessionData.Photostimpower),'mW','.pdf')),'-dpdf','-bestfit');
% print(gcf,fullfile(OutputDir,strcat('Ch',num2str(selectedchannel),'_Pulse_Closeup','.bmp')),'-bestfit');
export_fig(fullfile(OutputDir,strcat('Ch',num2str(selectedchannel),'_Pulse_Closeup_',num2str(sessionData.Photostimpower),'mW','.png')),'-png');

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
% area(edges(selectedpiece),hist_mat(selectedpiece)); hold all;

area(edges(selectedpiece),hist_mat(selectedpiece),'FaceColor','k'); hold all;
% plot(edges(selectedpiece),hist_mat(selectedpiece),'k'); hold all;
ymax    = max(hist_mat(selectedpiece))*1.1;
ymin    = min(hist_mat(selectedpiece))*1.1;

% ymax    = 120;

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
print(gcf,fullfile(OutputDir,strcat('Ch',num2str(selectedchannel),'_MUA_Trace','.pdf')),'-dpdf','-bestfit');
export_fig(fullfile(OutputDir,strcat('Ch',num2str(selectedchannel),'_MUA_Trace','.png')),'-png');
% print(gcf,fullfile(OutputDir,strcat('Ch',num2str(selectedchannel),'_MUA_Trace_',num2str(sessionData.Photostimpower),'mW','.pdf')),'-dpdf','-bestfit');
% export_fig(fullfile(OutputDir,strcat('Ch',num2str(selectedchannel),'_MUA_Trace_',num2str(sessionData.Photostimpower),'mW','.png')),'-png');

%% 
preposttime         = 0.025e6;
counter             = 1;
nSamples            = 2*preposttime*lfpData.fs*1e-6;
edges               = -preposttime*1e-6:1/32000:preposttime*1e-6-1/32000;

pulseresponse = NaN(sessionData.PhotostimFreq*numel(trialData.trialNum),nSamples);
for iT = 1:length(trialData.photostimStart)
    for iPu = 1:sessionData.PhotostimFreq
        pulsestart = trialData.photostimStart(iT)+1e6/sessionData.PhotostimFreq*(iPu-1);
        pulseresponse(counter,:) = lfpData.signal(lfpData.ts>(pulsestart-preposttime) & lfpData.ts<(pulsestart+preposttime));
        counter = counter+1;
    end
end

figure;
set(gcf,'units','normalized','Position',[0.2 0.1 0.5 0.7],'color','w')
pulseresponse = pulseresponse-repmat(nanmean(pulseresponse(:,edges<0),2),1,nSamples);
plot(edges,nanmean(pulseresponse,1),'k','LineWidth',3); hold all;
X = [0 0+pulselength*1e-6 0+pulselength*1e-6 0];
ymax = max(nanmean(pulseresponse,1))*1.1;
ymin = min(nanmean(pulseresponse,1))*1.1;
Y = [ymin ymin ymax ymax];
patch(X, Y,'blue','EdgeColor','blue','FaceAlpha',Alpha,'EdgeAlpha',0); hold on;

print(gcf,fullfile(OutputDir,strcat('Ch',num2str(selectedchannel),'_AveragePulse_',num2str(sessionData.Photostimpower),'mW','.pdf')),'-dpdf','-bestfit');
export_fig(fullfile(OutputDir,strcat('Ch',num2str(selectedchannel),'_AveragePulse_',num2str(sessionData.Photostimpower),'mW','.png')),'-png');

