
%MOL comments:

% - perhaps you can make the savitzky golay filter (and script in general)
% flexible to multiple input sampling frequencies? Although we will not
% have different values in practice, it might help with defining
% and understanding what the value of parameters mean.
% windowL = binL*2000; %why 2000? add comments... doesn't the window become
% twice the length? 1.2 seconds?

% what happened to SelectData?
% - try to make script work with lfp signal in cell array
% LF coherence? 


% - Many more additional comments are needed
% delate duplications
% the length of pupilbinned and spectrum did not match, error
%make it function: function[allChanData] = MOL_PupilLFP_coherence(Data)
%appending struct: field=fieldnames(allChanData), vertcat

%%% - sParams.deafultLagchange name to default lag
%%% - sParams.coh_bound change name to hf bound
%%% singleChannel change to iCh (shorter and more standardized?)
%%% - these three are not parameters: sParams.pupilFs,sParams.lfpFs,sParams.binnedFs   
%%% You apply a low-pass filter at the end of pupil preprocessing, but then
% again here, the savitzky golay filter. Is it possible to do this in one
% go, and at which step then?
%%% you have to round xL in correlation lag, script error



%This script is superclose to what we want: to make this script a function, to
%give it one session and it returns for every channel the optimal coherence
%and lag for low and high frequencies. Then if you run this for multiple
%session you can start to make averages for different areas :)
% - Overall very nice and concise script: the loops are appropriate and
% script is fast.
%%

%OUTPUT
%struct 'allChanData' with fields:
% corr       - correlation of LFP band and pupil data (for particular lag)
% coh_mean   - mean coherence between signals 
% coh_max    - max coherence between signals
% coh_max_Fr - what frequency gives the max coherence
% coh_meanLF - mean coherence for low frequencies between signals
% coh_maxLF  - max coherence for low frequencies between signals 
% coh_maxLF_Fr - what frequency gives the max coherence for low frequencies 

%if the option for finding the lag is chosen it will give 2 vectors
% lagmax - misalignment of LFP and pupil that gives maximum correlation (in sec)
% lagmin - misalignment of LFP and pupil that gives minimum correlation (in sec)

%% GET DATA (pupil, session, lfp)
[Data] = MOL_SelectData(); 
% frequency of lfp data MUST be 1000Hz, frequency of pupil data preferably 25Hz

%% PARAMETERS
%choose LFP frequency band for which you want to find correlations/coherence
%with pupil data
sParams.freq            = [3,10; 30,100];       %boundaries for frequency band (in Hz)
sParams.binL            = [0.6; 0.2];           %a bin size for each frequency band(in sec)
bandN                   = length(sParams.freq); %how many frequencies bands do you want

%Do you want to find lags for highest and lowest correlation?
sParams.checkLag        = false;         %if this is true script will only get vectors with lags    
sParams.checkLag_sec    = 4;             %up to how many seconds of lag between pupil and LFP will you check (in sec)
%If not, what is your lag?
sParams.deafultLag      = 0.8;           %(in sec)

sParams.coh_boundHF     = [0.1,0.4];     % which frequencies should you consider for coherence (in Hz)
sParams.coh_boundLF     = [0.01,0.1];    % which low frequencies should you consider for coherence (in Hz)

%choose method of computing length of window for coherence
% 'seconds' - give the exact length of window in seconds
% 'bins' - give number of bins you want to divide your data
sParams.cohWind         = 'bins'; %'seconds', 'bins'
sParams.cohWind_sec     = 50;    
sParams.cohWind_bin     = 100;  

sParams.sgFilter        = true;    %Whether to apply a Savitzky–Golay filter
sParams.sgFilterOrder   = 4;       %Polynomial fit (Savitzky–Golay filter)
sParams.sgFilterWindow  = 15;      %Size of the window length (Savitzky–Golay filter)

%sParams.lfpTs         = Data.lfpData.ts(1,:);

%% prepare pupil, lfp and final struct

numChannels     = length(Data.lfpData.signal);
allChanData     = struct();
for band = 1:bandN
    allChanData(band).area        = Data.lfpData.Area;
    allChanData(band).ChannelY    = Data.lfpData.ChannelY';
    allChanData(band).freq        = sParams.freq(band,:);
    allChanData(band).binL        = sParams.binL(band);
end

lagmax          = zeros(bandN,numChannels);
lagmin          = zeros(bandN,numChannels);

pupil           = cell2mat(Data.pupilData.area(:));
lfp             = cell2mat(Data.lfpData.signal(1));

% apply filter
if sParams.sgFilter == true
     pupil = sgolayfilt(pupil,sParams.sgFilterOrder,sParams.sgFilterWindow);
end

%cope information to the other struct that will change its properties with
%each freqency band
fStruct = sParams;

fStruct.pupilFs        = Data.pupilData.fs(1);
fStruct.lfpFs          = Data.lfpData.fs(1); % must be 1000Hz
fStruct.binnedFs       = 1./sParams.binL;

%% LOOP for each freqency band and each channel to get lag or correlation and coherence
for band = 1:bandN 
    
    fStruct.freq        = sParams.freq(band,:);
    fStruct.binL        = sParams.binL(band);
    fStruct.binnedFs    = fStruct.binnedFs(band);
    
    %% BIN DATA (do this for each band if bin lengths are different)
    if ~all(sParams.binL==sParams.binL(1)) ||...
          (all(sParams.binL==sParams.binL(1)) && band == 1)
      fStruct= MOL_PupilLFP_binning(pupil,lfp, fStruct);
      fStruct.binnedPupil = cell2mat(fStruct.binnedPupil);
    end
%     fStruct.binnedPupil = fStruct.binnedPupil(1:end-1);
    
    for iCh = 1:numChannels %for each channel
        %% TIME-FREQUENCY ANALYSIS
        lfp         = Data.lfpData.signal{iCh}(1:fStruct.lfpL);
        spect       = MOL_PupilLFP_spectrum(lfp, fStruct); 
        
        %% FIND LAG (if sParams.checkLag=true script won't proceed further)
        if fStruct.checkLag
            [lagmax(band ,iCh), lagmin(band ,iCh)] = MOL_PupilLFP_correlationLAG(spect, fStruct);
            if iCh ==1 || rem(iCh, 5)==0
                fprintf('Channel %d was processed [%s] \n', iCh,getTime)
            end
            continue
        end
        
        %% CORRELATIONS
        corrValue   = MOL_PupilLFP_correlation(spect, fStruct);
        
        %% COHERENCE
        [cxy, fcxy] = MOL_PupilLFP_coherence(spect, fStruct);
        
        %% FINAL STRUCT
        allChanData(band).coh_meanHF(iCh) = mean(cxy(fcxy>sParams.coh_boundHF(1) & fcxy<sParams.coh_boundHF(2)));
        allChanData(band).coh_maxHF(iCh) = max(cxy(fcxy>sParams.coh_boundHF(1) & fcxy<sParams.coh_boundHF(2)));
        allChanData(band).coh_maxHF_Fr(iCh) = fcxy(cxy == allChanData(band).coh_max(iCh));
        allChanData(band).coh_meanLF(iCh)= mean(cxy(fcxy>=sParams.coh_boundLF(1) & fcxy<sParams.coh_boundLF(2)));
        allChanData(band).coh_maxLF(iCh) = max(cxy(fcxy>=sParams.coh_boundLF(1)& fcxy<sParams.coh_boundLF(2)));
        allChanData(band).coh_maxLF_Fr(iCh)= fcxy(cxy == allChanData(band).coh_maxLF(iCh));
        
        allChanData(band).corr(iCh) = corrValue;
        
        if iCh ==1 || rem(iCh, 5)==0
            fprintf('Channel %d was processed [%s] \n', iCh,getTime)
        end
    end
end

%% VISUALISATION
band = 2; % which band do you want to visualize
recarea = {'V1', 'PPC', 'CG1'}; %which area? you can choose multiple
relation1 = allChanData(band).corr;
relation2 = allChanData(band).coh_max;

% not finished - no labels etc.
for areaL = 1: length(recarea)
    area = strcmp(allChanData(band).area, recarea(areaL));
    area_depth = allChanData(band).ChannelY(area);
    figure; scatter(relation1(area), relation2(area));
    figure; scatter(relation1(area), area_depth)
end

% not finished; plot pupil and lfp
channel = 10;
secN = 60;

lfp = cell2mat(Data.lfpData.signal(channel));
sessN = find(Data.pupilData.session_ID == Data.lfpData.session_ID(channel));
pupil = cell2mat(Data.pupilData.area(sessN));
zPupil = (pupil - nanmean(pupil))/nanstd(pupil);
plot(lfp(1:secN*Data.lfpData.fs(1)))
hold on; plot(zPupil(1:secN*Data.pupilData.fs(1)))