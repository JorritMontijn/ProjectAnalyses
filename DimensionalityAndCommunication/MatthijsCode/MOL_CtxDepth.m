function MOL_CtxDepth()

%% TODO:
%- check sink source naming
%- remove noise channels (how to assess which noisy?

RootDir                 = 'E:\';
ExternalDir             = 'J:\';

PreprocessProject       = 'PMAL';
PreprocessAllAnimals    = 0;        %Whether to preprocess all mice, instead of selection below

switch PreprocessProject
    case 'CHDET'
        PreprocessAnimals       = {'2011'}; %Which mouse to preprocess
    case 'MODDISCR'
        PreprocessAnimals       = {'P3.18'}; %Which mouse to preprocess
        %         PreprocessAnimals       = {'2.14' 'P3.3'}; %Which mouse to preprocess
    case 'PMAL'
        PreprocessAnimals       = {'2005' '2006' '2008' 'P3.13' 'P3.15' 'P3.16'}; %Which mouse to preprocess
        %         PreprocessAnimals       = {'2.7' '2.8' '1.1' '1.2' '1.3' '1.4' '1.37' '1.38' '1.40'}; %Which mouse to preprocess
end

%% Parameters
params.t_pre                = -0.1e6; %All time in microseconds
params.t_post               = .25e6;  %All time in microseconds
% params.AlignOnEvent         = 'stimChange';
params.AlignOnEvent         = 'stimStart';
params.interpolate          = 10; %times the spatial resolution increases
% params.cscale               = [-16 16]*1e-6; %Colorbar scale
params.colormap             = 'parula'; %redblue or parula
params.showFig              = 1;

params.sinkhot              = 1;
% There are two common conventions in CSD world. Here the convention is held
% that when net positive current enters the cell, the CSD will be negative and is termed a sink.
% As this is excitatory input, this will be represent by 'hot' colors;

%% Parameters for Butterworth filter
params.UseButter           = 1;
params.lp_butter           = 80;  %Low pass filter (Hz) (Butter)
params.ord_butter          = 4;    %Butterworth filter order

%Parameters for Kaiser filter
params.UseKaiser           = 0;
params.lp_kaiser           = 80;  %Low pass filter (Hz) (Kaiser)
params.hp_kaiser           = 0.5;    %High pass filter (Hz) (Kaiser)
params.dp_kaiser           = 1;    %Transition bandwidth (Kaiser only)

%% Configurations from excel overview of recordings per project:
MainDataDir             = fullfile(RootDir,'Data',PreprocessProject);
[~,~,RawLibrary]        = xlsread(fullfile(MainDataDir,sprintf('%s_Recordings_Overview.xlsx',PreprocessProject)));
RawLibrary              = RawLibrary(~cellfun(@isnumeric,RawLibrary(:,1)),:); %Trim the RawLibrary to entered values
AllMice                 = unique(RawLibrary(2:end,1));
[XlsRows,XlsColumns]    = size(RawLibrary);
if PreprocessAllAnimals
    PreprocessAnimals = AllMice;
end

&& any(strcmp({'V1','PM','AL'},sessionData.(sprintf('Probe%d_Area',iprobe))))

%% Loop over experiments %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for iAnimal = 1:length(PreprocessAnimals)
    % Find all experiments for this animal:
    ExpIdx               = find(strcmp(RawLibrary(:,strcmp(RawLibrary(1,:),'mousename')),PreprocessAnimals{iAnimal}));
    
    for iExp = ExpIdx' %Loop over experiments
        clear sessionData trialData spikeData lfpData pupilData
        
        %% sessionData (taking info from the excel recordings overview):
        for nCol = 1:XlsColumns
            sessionData.(RawLibrary{1,nCol}) = RawLibrary{iExp,nCol};
        end
        
        if strcmp(sessionData.Experiment,'CheckerboardReversals')
            sessionData.session_ID          = MOL_giveSessionID(sessionData.mousename,sessionData.Rec_datetime,iExp);
            
            %Set the directories for raw and eventually preprocessed data:
            sessionRootDir              = fullfile(MainDataDir,'RawData',sessionData.mousename,sessionData.Rec_datetime);
            RawDataDir                  = fullfile(sessionRootDir,'RawData');
            PreprocessDataDir           = fullfile(MainDataDir,sessionData.Experiment,sessionData.mousename,sessionData.Rec_datetime);
            
            [sessionData, trialData]    = MOL_Preprocessing_Checkerboard(sessionData,sessionRootDir);
            
            %The time window that the LFP is extracted (from start to end of this protocol)
            TimeWindow = [sessionData.t_start sessionData.t_stop];
            if sessionData.AdapterFlipped
                [ChannelSel] = MOL_AdapterChFlip(ChannelSel);
                error('script not designed to work with flipped adapter')
            end
            
            nPreviousChannels = 0;
            for iprobe = 1:3
                ChannelSel = [];             %Make a list of the channels to be extracted to LFP (often list of all channels):
                AREAS = {'V1' 'PM' 'AL' 'PPC' 'CG1'};
                if sessionData.(sprintf('Probe%d',iprobe)) && ~strcmp(sessionData.(sprintf('Probe%d_Config',iprobe)),'A4x8-5mm-100-200-177-CM32')...
                    switch sessionData.(sprintf('Probe%d_Config',iprobe))
                        case 'A4x8-5mm-100-200-177-CM32'
                            nChannels           = 32;
                            nShanks             = 4;
                            nChannelsShank      = 8;
                            intersitedistance   = 100;
                        case 'A1x32-Poly2-5mm-50s-177-CM32'
                            nChannels           = 32;
                            nShanks             = 2;
                            nChannelsShank      = 16;
                            intersitedistance   = 50;
                        case 'A1x64-Poly2-6mm-23s-160'
                            nChannels           = 64;
                            nShanks             = 2;
                            nChannelsShank      = 32;
                            intersitedistance   = 46;
                    end
                    ChannelSel                  = 1:nChannels;
                    
                    
                        CSDfile         = fullfile(sessionRootDir,'CSD','csdData.mat');
                        loadstruct      = load(CSDfile,'csdData');
                        csdData         = loadstruct.csdData;
                        
                        save(fullfile(sessionRootDir,'CSD',sprintf('CSD_%d_%s.mat',lfpData.session_ID(1),lfpData.Area{1})),'csdData','-v7.3');
                        
                        
                        xtime                           = (params.t_pre:1e6/lfpData.fs:params.t_post-0.01e6) * 1e-6;
                        
                        
                        %                     %Interpolate values for visualization purposes:
                        %                     newChannelY = linspace(min(min([ChannelY_shank{:}])),max(max([ChannelY_shank{:}])),100);
                        %
                        %                     if params.interpolate
                        %                         meancsd_shank{1}                            = tointerpol2(meancsd_shank{1},ChannelY_shank{1}',newChannelY);
                        %                         meancsd_shank{2}                            = tointerpol2(meancsd_shank{2},ChannelY_shank{2}',newChannelY);
                        %                     end
                        %                     meancsd_shank{1}        = meancsd_shank{1}/max(max(meancsd_shank{1}));
                        %                     meancsd_shank{2}        = meancsd_shank{2}/max(max(meancsd_shank{2}));
                        %
                        %                     meancsd                 = cat(3,meancsd_shank{1},meancsd_shank{2});
                        %                     meancsd                 = nanmean(meancsd,3);
                        %
                        %                     params.cscale           = [-max(max(meancsd))*0.9 max(max(meancsd))*0.9];
                        
                        csdfig = figure; set(gcf,'units','normalized','Position',[0.05 0.4 0.9 0.4],'color','w');
                        set(gcf,'defaultAxesFontSize',15)
                        suptitle(sprintf('CSD %s - %s',sessionData.mousename,sessionData.(sprintf('Probe%d_Area',iprobe))));
                        
                        subplot(1,3,1)
                        offsetmat = repmat(ChannelY_meanLFP,1,size(lfpmat(:,:,1),2));
                        plot(xtime,mean(lfpmat(:,:,:),3)*20e4 + offsetmat); hold on;
                        xlabel('Time from stimulus (s)','FontSize', 15)
                        ylabel('Channel','FontSize', 15)
                        plot([0 0],ylim,'k','LineWidth',2);
                        xlim([xtime(1) xtime(end)]);
                        %         ylim([min(ChannelY_meanLFP)-50 max(ChannelY_meanLFP)+50]);
                        title('Mean LFP (<50Hz) to Checker reversal','FontSize',15)
                        ylabel('Channel Depth (in um from dura)','FontSize', 15)
                        
                        subplot(1,3,2)
                        set(gcf,'defaultAxesFontSize',15)
                        %         imagesc(xtime,flipud(newChannelY),flipud(meancsd),params.cscale)
                        imagesc(xtime,newChannelY,meancsd,params.cscale)
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
                        ticks = linspace(params.cscale(1),params.cscale(2),9);
                        c = colorbar('Ticks',ticks,'TickLabels',num2cell(ticks*1e6),'FontSize',10);
                        c.Label.String = 'Sink                    (uA/mm3)                 Source';
                        
                        
                        for i = 1:3
                            subplot(1,3,i)
                            ylim([min(ChannelY_meanLFP)-50 max(ChannelY_meanLFP)+50]);
                            
                            Vert_Ticks = -2000:100:400;
                            Vert_Ticks = Vert_Ticks(Vert_Ticks>min(ChannelY_meanLFP)-50 & Vert_Ticks<max(ChannelY_meanLFP)+50);
                            
                            set(gca, 'YTick', Vert_Ticks, 'YTickLabels', Vert_Ticks,'FontSize', 15)
                        end
                        
                        csdData                 = struct();
                        csdData.meanCSD         = meancsd;
                        csdData.meanERP         = mean(lfpmat(:,:,:),3);
                        csdData.time            = xtime;
                        csdData.hfPower         = hf_pwr';
                        csdData.ChannelY_sorted = ChannelY_meanLFP;
                        csdData.channelnum      = sortidx+nPreviousChannels; %#ok<STRNU>
                        
                        fprintf('Saving CSD data...\n')
                        mkdir(fullfile(sessionRootDir,'CSD'))
                        save(fullfile(sessionRootDir,'CSD','csdData.mat'),'csdData','-v7.3');
                        %                         figtitle = fullfile(sessionRootDir,'CSD',sprintf('%d_%s_CSD.pdf',lfpData.session_ID(1),lfpData.Area{1}));
                        figtitle = fullfile(sessionRootDir,'CSD',sprintf('%d_%s_CSD.bmp',lfpData.session_ID(1),lfpData.Area{1}));
                        %                         print(csdfig,figtitle,'-dpdf','-bestfit');
                        saveas(csdfig,figtitle, 'bmp')
                        figtitle = fullfile('E:\Documents\PhD\Figures\CSD',sprintf('%d_%s_CSD.bmp',lfpData.session_ID(1),lfpData.Area{1}));
                        saveas(csdfig,figtitle, 'bmp')
                        
                    end
                    
                end
                nPreviousChannels           = nPreviousChannels + nChannels;
                %             end
            end
        end
    end
    
end


    function calcCSD(lfpData)
        
    end

    function out = calcMUA(lfpData)
        
        %Get power estimate for different frequency bands:
        [pwr_out,pwr_f] = pwelch(lfpData.signal(ChannelSel,:)',[],[],[],lfpData.fs(1,1),'power');
        
        hf_pwr = sum(pwr_out(pwr_f>500 & pwr_f<5000,:),1);
        hf_pwr = hf_pwr/max(hf_pwr); %normalize to maximum
        
        subplot(1,3,3)
        [ChannelY_meanLFP,sortidx]  = sort(-lfpData.ChannelY(ChannelSel));
        hf_pwr = hf_pwr(sortidx);
        plot(hf_pwr,ChannelY_meanLFP)
        
    end

