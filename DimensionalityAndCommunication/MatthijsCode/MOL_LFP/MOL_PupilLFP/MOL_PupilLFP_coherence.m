function [cxy, fcxy] = MOL_PupilLFP_coherence(spect, fParams)
%% MOL_PupilLFP_coherence finds the magnitude-squared coherence estimate
%  for pupil and sum of frequency magnitudes (lfp)

%You need in fParams struct: binnedPupil, binnedFs, freq, cohWind, cohWind_bin

    % add magnitudes of all frequencies
    sumFreq = sum(abs(spect(fParams.freq(1):fParams.freq(2),1:end)));
    
    % take z-score from pupil and time-sParams.frequency data
    zPupil = (fParams.binnedPupil - nanmean(fParams.binnedPupil))/nanstd(fParams.binnedPupil);
    zSumFreq = (sumFreq - (nanmean(sumFreq)))/nanstd(sumFreq);
     
    zPupil(isnan(zPupil)) = 0;
    
    %choose method for calculating window for coherence
    if strcmp(fParams.cohWind, 'seconds')  
        coh_window = sParams.cohWind_sec*fParams.binnedFs;
        
    elseif strcmp(fParams.cohWind, 'bins') 
        coh_window = (round(length(fParams.binnedPupil)/fParams.cohWind_bin));   
    end
    
    noverlap = round(coh_window/fParams.binnedFs);
    f = [];
    
    [cxy,fcxy] = mscohere(zPupil,zSumFreq,coh_window,noverlap,f,fParams.binnedFs);
   
%     figure; plot(fcxy, cxy) 

end
