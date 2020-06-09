function MOL_UpdateQualOverview()

global mainfig params lfpData trialData

figure(mainfig)

%% Panel 1: raw traces:
subplot(1,3,1);
cla;

% If there is a task, then take a random lick, extract lfp surrounding it
% to investigate movement artefacts:
dirfiles = dir(fullfile(params.sessionRootDir,'ChangeDetection*'));
if ~isempty(dirfiles)
    eventfilename = fullfile(params.sessionRootDir,dirfiles.name,'Events.nev');
    %Settings for getting events
    FieldSelection      = [1 0 1 0 1]; %Fields to be read: TimeStamps, TTLs, EventStrings
    ExtractHeader       = 0;
    ExtractMode         = 1; %Extract all data, from beginning to end
    [TimeStamps, TTLs, EventStrings] =   Nlx2MatEV(eventfilename,FieldSelection,ExtractHeader,ExtractMode);
    TTL_leftLick        = 1; %B0001xxxx
    TTL_rightLick       = 2; %B0010xxxx
    Port_leftLick      	= 2;
    Port_rightLick      = 2;
    EventStrings = EventStrings';
    EventPort = NaN(1,length(EventStrings));
    for i = 1:length(EventStrings)
        EventPort(i) = str2double(EventStrings{i}(strfind(EventStrings{i},'port ')+5));
    end
    TS_Licks            = TimeStamps((EventPort == Port_leftLick & TTLs == TTL_leftLick) | (EventPort == Port_rightLick & TTLs == TTL_rightLick));
    ChosenLick          = TS_Licks(round(length(TS_Licks)/randi(9)));
    TimeWindow          = [ChosenLick + params.t_pre ChosenLick + params.t_post];
    
    %Parameters for filters and resampling:
    prm.UseButter               = 0;
    prm.UseKaiser               = 0;
    prm.UseResampling           = 0;    %Resample data to lower frequency
    %Extract the lfp in a subfunction:
    lfpData_lick                = MOL_extractLFP(fullfile(params.sessionRootDir,'RawData'),params.OrigChannelSel,TimeWindow,prm);
    
    if params.flipsignLFP
        lfpData_lick.signal = cellfun(@(x) -x, lfpData_lick.signal, 'UniformOutput',false);
    end
    
    lfpData_lick.ts             = repmat(lfpData_lick.t_start(1):1/lfpData_lick.fs(1)*1e6:lfpData_lick.t_end(1),1,1);
    
    timeselec                   = lfpData_lick.ts> TimeWindow(1) & lfpData_lick.ts<TimeWindow(2);
    offsetmat                   = repmat(lfpData.sortedChannelDepth(lfpData.sortedisgood),1,length(params.xtime));

    lfpData_lick.sortedSignal                = cell2mat(lfpData_lick.signal);
    lfpData_lick.sortedSignal                = lfpData_lick.sortedSignal(lfpData.sortedChannelidx,:);
    
    %Plot:
    baseline = repmat(mean(lfpData_lick.sortedSignal(lfpData.sortedisgood,timeselec),2),1,sum(timeselec));
    handles = plot(params.xtime,(lfpData_lick.sortedSignal(lfpData.sortedisgood,timeselec)-baseline)*20e4 + offsetmat,'LineWidth',0.2); hold on;
    %add legend to identify channels by their colors:
    legend(flipud(handles),num2cell(num2str(flipud(lfpData.sortedChannelNum(lfpData.sortedisgood))),2),'Location','NorthEast','FontSize',350/params.nChannels);
else
    iTrial                      = randi(length(trialData.(params.AlignOnEvent)));
    timeselec                   = lfpData.ts>trialData.(params.AlignOnEvent)(iTrial)+params.t_pre & lfpData.ts<trialData.(params.AlignOnEvent)(iTrial)+params.t_post;
    offsetmat                   = repmat(lfpData.sortedChannelDepth(lfpData.sortedisgood),1,length(params.xtime));
    
    %Plot:
    baseline = repmat(mean(lfpData.sortedSignal(lfpData.sortedisgood,timeselec),2),1,sum(timeselec));
    handles = plot(params.xtime,(lfpData.sortedSignal(lfpData.sortedisgood,timeselec)-baseline)*20e4 + offsetmat,'LineWidth',0.2); hold on;
    %add legend to identify channels by their colors:
    legend(flipud(handles),num2cell(num2str(flipud(lfpData.sortedChannelNum(lfpData.sortedisgood))),2),'Location','NorthEast','FontSize',350/params.nChannels);
end

%Make up:
set(gca,'XTick',params.xticks,'XTickLabels',params.xticklabels)
xlabel('Time from stimulus (ms)','FontSize', 15);
ylabel('Channel','FontSize', 15)
plot([0 0],params.Ylim_all,'k','LineWidth',2);
xlim([params.xtime(1) params.xtime(end)]);
title('Raw Traces','FontSize',15)
ylabel('Channel Depth (in um from dura)','FontSize', 15)

%%  MUA power:
subplot(1,3,2);
cla; hold all;

%Get power estimate for different frequency bands:
%                         [pwr_out,pwr_f] = pwelch(sortedlfpmat,[],[],[],lfpData.fs(1,1),'power');
[pwr_out,pwr_f]         = pwelch(lfpData.sortedSignal(:,1:20000)',[],[],[],lfpData.fs(1,1),'power');

hf_pwr                  = sum(pwr_out(pwr_f>params.hf_lower & pwr_f<params.hf_upper,:),1);
hf_pwr                  = hf_pwr/max(hf_pwr(lfpData.sortedisgood)); %normalize to maximum
plot(hf_pwr(lfpData.sortedisgood),lfpData.sortedChannelDepth(lfpData.sortedisgood),'-.k*','LineWidth',3)

linenoise_pwr                  = sum(pwr_out(pwr_f>49 & pwr_f<51,:),1);
linenoise_pwr                  = linenoise_pwr/max(linenoise_pwr(lfpData.sortedisgood)); %normalize to maximum
plot(linenoise_pwr(lfpData.sortedisgood),lfpData.sortedChannelDepth(lfpData.sortedisgood),'-.r','LineWidth',1)

theta_pwr                  = sum(pwr_out(pwr_f>8 & pwr_f<12,:),1);
theta_pwr                  = theta_pwr/max(theta_pwr(lfpData.sortedisgood)); %normalize to maximum
plot(theta_pwr(lfpData.sortedisgood),lfpData.sortedChannelDepth(lfpData.sortedisgood),'-.g','LineWidth',1)

xlabel('Normalized MUA power')
legend({'MUA' '50 Hz' 'Theta'})

%% Panel 3: CSD power
subplot(1,3,3)
cla;
[meancsd,meanresp]                  = calcCSD(params,lfpData,trialData);
figure(mainfig);
params.cscale                       = [-max(max(meancsd))*0.9 max(max(meancsd))*0.9];

imagesc(params.xtime,lfpData.sortedChannelDepth,meancsd,params.cscale); hold on;
xlabel('Time from stimulus (s)','FontSize', 15)
set(gca,'YDir','normal')

switch params.colormap
    case 'redblue'
        h       = coolwarm + repmat(0.1-abs(linspace(-0.1,.1,64))',1,3);
        h(h>1)  = 1;
        if params.sinkhot
            h = flipud(h);
        end
        colormap(h);
    case 'parula'
        if params.sinkhot
            colormap(flipud(parula));
        else                     colormap(parula);
        end
end

offsetmat = repmat(lfpData.sortedChannelDepth,1,length(params.xtime));
plot(params.xtime,meanresp*20e4 + offsetmat,'k','LineWidth',0.5); hold on;

ylabel('Channel','FontSize', 15)
plot([0 0],ylim,'k','LineWidth',1);
xlim([params.xtime(1) params.xtime(end)]);
set(gca,'XTick',params.xticks,'XTickLabels',params.xticklabels)
xlabel('Time from stimulus (ms)','FontSize', 15)

%% General make up of figures
sepfraction = 0.08;
for i = 1:3
    subplot(1,3,i);
    ylim(params.Ylim_all);
    x_width = (params.xmax-params.xmin)/3.1;
    set(gca,'Position',[params.xmin+(i-1)*x_width params.ymin x_width-sepfraction params.ymax-params.ymin])
    if i==3
        set(gca, 'YTick', lfpData.sortedChannelDepth, 'YTickLabels', lfpData.sortedChannelNum,'FontSize', 15)
    else
        Vert_Ticks = -2000:100:400;
        Vert_Ticks = Vert_Ticks(Vert_Ticks>min(lfpData.sortedChannelDepth)-25 & Vert_Ticks<max(lfpData.sortedChannelDepth)+25);
        set(gca, 'YTick', Vert_Ticks, 'YTickLabels', Vert_Ticks,'FontSize', 15)
    end
end

lfpData.CheckerCSD                  = meancsd;
lfpData.CheckerERP                  = meanresp;
lfpData.HF_PWR                      = hf_pwr';

end

function [meancsd,meanresp] = calcCSD(params,lfpData,trialData)
%make butter filter:
[params.B_butt,params.A_butt]   = butter(params.ord_butter,params.lp_butter/(lfpData.fs(1)/2));

if params.UseKaiser
    [N,Wn,beta,ftype]   = kaiserord([params.hp_kaiser-params.dp_kaiser params.hp_kaiser params.lp_kaiser params.lp_kaiser+params.dp_kaiser],[0 1 0],[1e-6 1 1e-6],lfpData.fs(1));
    B_kai               = fir1(N,Wn,ftype,kaiser(N+1,beta),'noscale');
    A_kai               = 1;
end
    
nEvents                         = length(trialData.(params.AlignOnEvent));
lfpmat                          = zeros(params.nChannels,length(params.xtime),nEvents); %init lfp matrix

%Get lfp per trial and apply low-pass filter:
for ev = 1:nEvents
    
    tempsignal      = lfpData.sortedSignal(:,lfpData.ts>trialData.(params.AlignOnEvent)(ev)+params.t_pre & lfpData.ts<trialData.(params.AlignOnEvent)(ev)+params.t_post);
    lfpmat(:,:,ev)  = tempsignal(:,1:length(params.xtime));
    
    if params.UseButter
        for iCh = 1:params.nChannels
            lfpmat(iCh,:,ev)            = filtfilt(params.B_butt,params.A_butt,lfpmat(iCh,:,ev));
        end
    end
    
    if params.UseKaiser
        for iCh = 1:params.nChannels
            lfpmat(iCh,:,ev)            = filtfilt(B_kai,A_kai,lfpmat(iCh,:,ev));
        end
    end
       
end

% Replace bad channels with mean of adjacent channels
for iCh = find(~lfpData.sortedisgood)'
    chbelow          = find(lfpData.sortedChannelDepth > lfpData.sortedChannelDepth(iCh) & lfpData.sortedisgood,1,'first');
    chabove          = find(lfpData.sortedChannelDepth < lfpData.sortedChannelDepth(iCh) & lfpData.sortedisgood,1,'last');
    lfpmat(iCh,:,:)  = mean(lfpmat([chbelow chabove],:,:),1);
end

%Loop over the different shanks (CSD does not work well with
%considering no horizontal spacing)
for iShank = 1:params.nShanks
    ShankChannelSel                     = params.ChannelSel([1:params.nChannelsShank] + params.nChannelsShank*(iShank-1));
    
    idx                                 = ismember(lfpData.sortedChannelNum,ShankChannelSel);
    sortedShankChannelNum               = lfpData.sortedChannelNum(idx);
    sortedShankChannelDepth             = lfpData.sortedChannelDepth(idx);
    %     sortedShankChannelIdx       = lfpData.sortedChannelNum(ismember(lfpData.sortedChannelNum,ShankChannelSel));
    
    shanklfpmat                         = lfpmat(idx,:,:);

    switch params.csdmethod
        case 'DoubleDiff'
            meancsd                             = csdfromlfp_doublediff(mean(shanklfpmat,3));
        case 'NicholsonFreeman'
            meancsd                             = csdfromlfp_NicholsonFreeman(mean(shanklfpmat,3),lfpData);
    end
   
    %Baseline subtraction of csd:
    meancsd_shank{iShank}               = meancsd - repmat(mean(meancsd(:,params.xtime<-0.025),2),1,length(params.xtime));
    ChannelY_shank{iShank}              = sortedShankChannelDepth;
end

meanresp                = mean(lfpmat(:,:,:),3);
baseline                = repmat(mean(meanresp(:,params.xtime<0),2),1,length(params.xtime));
meanresp                = meanresp - baseline;

meancsd_shank{1}        = tointerpol2(meancsd_shank{1},ChannelY_shank{1}',lfpData.sortedChannelDepth);
meancsd_shank{2}        = tointerpol2(meancsd_shank{2},ChannelY_shank{2}',lfpData.sortedChannelDepth);

meancsd_shank{1}        = meancsd_shank{1}/max(max(meancsd_shank{1}));
meancsd_shank{2}        = meancsd_shank{2}/max(max(meancsd_shank{2}));

meancsd                 = cat(3,meancsd_shank{1},meancsd_shank{2});
meancsd                 = nanmean(meancsd,3);

end

function csd = csdfromlfp_doublediff(lfp)

% Vaknin transform: (Vaknin 1985: add channels above and below assuming isoconductivity):
lfp_vaknin      = [lfp(1,:); lfp(1,:); lfp; lfp(end,:); lfp(end,:)];

% Filter temporally and spatially: is key to this
spat_filter     = fspecial('gaussian',[3 5],1.1);
lfp_vaknin      = conv2(lfp_vaknin,spat_filter,'same');

% compute the CSD map
csd             = diff(lfp_vaknin,2,1); %Taking the double derivative from up to down

% Remove first and last channels of csd map to ensure same size as original
% lfp matrix and set outer channels to zero (nonsensical csd values)
csd             = [zeros(1,size(csd,2));csd(3:end-2,:);zeros(1,size(csd,2))];

end


function csd = csdfromlfp_NicholsonFreeman(lfp,lfpData)

% Vaknin transform: (Vaknin 1985: add channels above and below assuming isoconductivity):
lfp_vaknin      = [lfp(1,:); lfp; lfp(end,:);];

% Filter temporally and spatially: is key to obtaining a smooth csd profile (filter over 3 channels and 5 time points)
spat_filter     = fspecial('gaussian',[3 5],1.1);
lfp_vaknin      = conv2(lfp_vaknin,spat_filter,'same');

% compute the CSD map:
csd             = CSD(lfp_vaknin',lfpData.fs(1),50*1e-6,'conductivity',0.4);
csd             = csd'; %Transpose

% Remove first and last 2 channels of csd map to ensure same size as original
% lfp matrix and set outer channels to zero (nonsensical csd values)
csd             = [zeros(1,size(csd,2));csd(3:end-2,:);zeros(1,size(csd,2))];

end