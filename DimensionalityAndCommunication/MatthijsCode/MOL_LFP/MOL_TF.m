function MOL_TF(varargin)

params                  = struct();
params.RootDataDir      = 'E:\Data\';
params.Project          = 'CHDET';
params.Protocol         = 'ChangeDetectionConflict';
params.Animals          = {'2009'}; %Which mouse to preprocess

params.loadsessionData  = 1; 
params.loadtrialData    = 1; 
params.loadspikeData    = 0; 
params.loadlfpData      = 1; 
params.loadpupilData    = 0; 

[sessionData,trialData,~,lfpData,~] = MOL_getDataInputs(params);

%% Parameters for time-frequency analysis
% params = params_timefreq();
%% parameters for window of interest:
params.t_pre                = -0.5e6;
params.t_post               = .5e6;

%% Parameters for analysis:
params.eventofinterest      = 'stimStart';

%% 
for sesid = unique(sessionData.session_ID)'
    %% Get the relevant data for each session individually:
    [tempsessionData,temptrialData,templfpData] = MOL_getTempPerSes(sesid,sessionData,trialData,lfpData);
    
    figure; hold all;
    
    for iCh = 50%2%:32; %length(templfpData.session_ID)
        if ischar(templfpData.signal{iCh})
            datadir = fullfile(params.RootDataDir,params.Project,params.Protocol,tempsessionData.mousename{1},tempsessionData.Rec_datetime{1});
            filename = fullfile(datadir,sprintf('LFP_%d_CSC%d.dat',tempsessionData.session_ID,iCh));

            m       = memmapfile(filename,'Format','double');
            
%             trialselec = find(strcmp(temptrialData.trialType,'X') & abs(temptrialData.visualOriChange)==90);
            trialselec = find(strcmp(temptrialData.trialType,'X')); %& temptrialData.noResponse==1);
            Sxx_alltrials = NaN(2,length(trialselec),params.t_post/1e3/2+1);
            for iTr = 1:length(trialselec)
                %                 timeselec = lfpData.ts{iCh}>temptrialData.stimChange(i)+params.t_pre & lfpData.ts{iCh}<temptrialData.stimChange(i)+params.t_post;
                trial = trialselec(iTr);
                timeselec(1,:) = lfpData.ts{iCh}>temptrialData.stimChange(trial)+params.t_pre & lfpData.ts{iCh}<temptrialData.stimChange(trial);
                timeselec(2,:) = lfpData.ts{iCh}>temptrialData.stimChange(trial) & lfpData.ts{iCh}<temptrialData.stimChange(trial)+params.t_post;
                colorset = {[0 0 0],[0.5 0 .5]};
                for iB = 1:2
                    
                    data    = m.Data(timeselec(iB,:));
                    
                    dt = 1/templfpData.fs(1);
                    T = (templfpData.ts{iCh}(timeselec(iB,:))-templfpData.ts{iCh}(find(timeselec(iB,:),1,'first')))*1e-6+1e-3;
                    x = data';
                    xf = fft(x); %1. Compute the Fourier transform of x.
                    Sxx = (2*dt^2./T) .* (xf.*conj(xf)); %2. Compute the power spectrum.
                    Sxx = Sxx(1:length(x)/2+1); %3. Ignore negative frequencies.
                    Sxx = conv(Sxx,ones(5,1),'same');
                    df = 1/max(T); %4. Determine the frequency resolution.
                    fNQ=1/dt/2; %5. Determine the Nyquist     frequency.
                    faxis = (0:df:fNQ); %6. Construct the frequency axis.
%                     faxis = log10(faxis);  Sxx = log10(Sxx);
%                     plot(faxis, Sxx,'Color',colorset{iB},'LineWidth',0.3) %7. Plot power versus frequency.
                    
                    Sxx_alltrials(iB,iTr,:) = Sxx;
                    plot(data,'Color',colorset{iB});
                end
            end
            
        end
        
%         plot(faxis, squeeze(mean(Sxx_alltrials(1,:,:),2)),'Color',colorset{1},'LineWidth',3) %7. Plot power versus frequency.
%         plot(faxis, squeeze(mean(Sxx_alltrials(2,:,:),2)),'Color',colorset{2},'LineWidth',3) %7. Plot power versus frequency.
    end
    
    
    xlim([3 100]) %8. Select frequency range.
    xlabel('Frequency [Hz]');
    ylabel('Power')
    set(gca,'xscale','log','yscale','log')
    
    resol = 10000;
    chunks = 1:resol:length(data);
    for iChunk = 1:length(chunks)-1
        [pxx(iChunk,:),f] = pmtm(data(chunks(iChunk):chunks(iChunk+1)),3,resol,templfpData.fs(iCh)); %returns a PSD computed as a function of physical frequency.
    end
    
    figure; plot(f(f<100),mean(pxx(:,f<100),1))
    
    end
    data        = data(1:10000);
    pxx         = pmtm(data);
    
    

    [pxx,f] = pmtm(data,3,length(data),templfpData.fs(iCh)); %returns a PSD computed as a function of physical frequency.

    figure; plot(f(f<50),pxx(f<50));
    fft
    pspectrum
    coloraxis = caxis;
    caxis(coloraxis);
    
    %% 
%     figurefilename = char(strcat(cfg.Outputroot, strtok(cfg.Animalroot,'\'), '.', strtok(cfg.Dateroot,'\'), '.', strtok(cfg.Bundle), '.TFR_Stimulus','.jpg'));
%     saveas(gcf, figurefilename);
%     fprintf('Figure is saved as %s\n', figurefilename);
end
end
