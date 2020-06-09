function MOL_CSD(varargin)

%% Get input arguments:
sessionData     = varargin{1};
trialData       = varargin{2};
% spikeData       = varargin{3};
lfpData         = varargin{4};


%% DFT filter on 50 Hz:
% L = size(lfpData.signal,2);
% f = lfpData.fs*(0:(L/2))/L;
% for iCh = 1:length(lfpData.ts)
%     fprintf('Filtering 50Hz noise channel %d\n',iCh)
%     h = fft(lfpData.signal(iCh,:));
%     h(find(f>49 & f<51)) = complex(zeros(numel(find(f>49 & f<51)),1));
%     lfpData.signal(iCh,:) = ifft(h);
% end

%%
nPreviousChannels = 0; %counter that keeps track of cumulative channel number
for iProbe = 1:3
    if sessionData.(sprintf('Probe%d',iProbe))
        switch sessionData.(sprintf('Probe%d_Config',iProbe))
            case 'A1x32-Poly2-5mm-50s-177-CM32'
                nChannels           = 32;
                nShanks             = 2;
                nChannelsShank      = 16;
                intersitedistance   = 50;
            case 'A4x8-5mm-100-200-177-CM32'
                nChannels           = 32;
                nShanks             = 4;
                nChannelsShank      = 8;
                intersitedistance   = 100;
            case 'A1x64-Poly2-6mm-23s-160'
                nChannels           = 64;
                nShanks             = 2;
                nChannelsShank      = 32;
                intersitedistance   = 46;
        end
        
        ChannelSel                      = [1:nChannels] + nPreviousChannels;
        nPreviousChannels               = nPreviousChannels + nChannels;
        xtime                           = (params.t_pre:1e6/lfpData.fs:params.t_post-0.01e6) * 1e-6;
        nEvents                         = length(trialData.(params.AlignOnEvent));
        lfpmat                          = zeros(nChannels,length(xtime),nEvents); %init lfp matrix
        
        %%Get lfp per trial
        %make butter filter:
        [params.B_butt,params.A_butt]     = butter(params.ord_butter,params.lp_butter/(lfpData.fs(1)/2));

        for ev = 1:nEvents
            tempsignal      = lfpData.signal(ChannelSel,lfpData.ts{1,:}>trialData.(params.AlignOnEvent)(ev)+params.t_pre & lfpData.ts{1,:}<trialData.(params.AlignOnEvent)(ev)+params.t_post);
            lfpmat(:,:,ev)  = tempsignal(:,1:length(xtime));
            
            if params.UseButter
                for iCh = 1:nChannels
                    lfpmat(iCh,:,ev)            = filtfilt(params.B_butt,params.A_butt,lfpmat(iCh,:,ev));
                end
            end
        end
        
        %Loop over the different shanks (CSD does not work well with
        %considering no horizontal spacing)
        for iShank = 1:nShanks
            ShankChannelSel                     = ChannelSel([1:nChannelsShank] + nChannelsShank*(iShank-1));
            
            [ShankChannelY,sortidx]             = sort(-lfpData.ChannelY(ShankChannelSel),'ascend');
            sortidx                             = sortidx + nChannelsShank*(iShank-1);
            
            shanklfpmat                         = NaN(nChannelsShank,size(lfpmat,2),size(lfpmat,3));
            shanklfpmat(:,:,:)                  = lfpmat(sortidx,:,:);
            
            %Init csd matrix:
            csd                                 = zeros(size(shanklfpmat));
            for ev = 1:size(shanklfpmat,3)
                csd(:,:,ev)                     = csdfromlfp(shanklfpmat(:,:,ev));
            end
            meancsd                             = mean(csd,3);
            meancsd                             = meancsd / (intersitedistance/100);

            %Baseline subtraction of csd:
            meancsd_shank{iShank}               = meancsd - repmat(mean(meancsd(:,xtime<-0.025),2),1,length(xtime));
            ChannelY_shank{iShank}              = ShankChannelY;
        end
        
        %Interpolate values for visualization purposes:
        newChannelY = linspace(min(min([ChannelY_shank{:}])),max(max([ChannelY_shank{:}])),100);

        if params.interpolate
            meancsd_shank{1}                            = tointerpol2(meancsd_shank{1},ChannelY_shank{1}',newChannelY);
            meancsd_shank{2}                            = tointerpol2(meancsd_shank{2},ChannelY_shank{2}',newChannelY);
        end
        
        meancsd = cat(3,meancsd_shank{1},meancsd_shank{2});
        meancsd = nanmean(meancsd,3);
        
        params.cscale = [-max(max(meancsd))*0.9 max(max(meancsd))*0.9];
        
        [ChannelY_meanLFP,sortidx]  = sort(-lfpData.ChannelY(ChannelSel));
        lfpmat                      = lfpmat(sortidx,:,:);
        
        csdfig = figure; set(gcf,'units','normalized','Position',[0.05 0.4 0.9 0.4],'color','w');
        set(gcf,'defaultAxesFontSize',15)
        suptitle(sprintf('CSD %s - %s',sessionData.mousename,sessionData.(sprintf('Probe%d_Area',iProbe))));
        
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
        
        %Get power estimate for different frequency bands:
        [pwr_out,pwr_f] = pwelch(lfpData.signal(ChannelSel,:)',[],[],[],lfpData.fs(1,1),'power');
        
        hf_pwr = sum(pwr_out(pwr_f>500 & pwr_f<5000,:),1);
        hf_pwr = hf_pwr/max(hf_pwr); %normalize to maximum
        
        subplot(1,3,3)
        [ChannelY_meanLFP,sortidx]  = sort(-lfpData.ChannelY(ChannelSel));
        hf_pwr = hf_pwr(sortidx);
        plot(hf_pwr,ChannelY_meanLFP)

        for i = 1:3
            subplot(1,3,i)
            ylim([min(ChannelY_meanLFP)-50 max(ChannelY_meanLFP)+50]);
            
            Vert_Ticks = -2000:100:400;
            Vert_Ticks = Vert_Ticks(Vert_Ticks>min(ChannelY_meanLFP)-50 & Vert_Ticks<max(ChannelY_meanLFP)+50);
            
            set(gca, 'YTick', Vert_Ticks, 'YTickLabels', Vert_Ticks,'FontSize', 15)
        end
                
    end
end

end


function csd = csdfromlfp(lfp)

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