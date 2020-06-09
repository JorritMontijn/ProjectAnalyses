function [spect] = MOL_PupilLFP_spectrum(lfp, fStruct)
% You need in fStruct struct: lfpFs, binL

%% FREQUENCY ANALYSIS
    binL  = fStruct.binL;
    
    windowL = binL*2000;
    window  = hamming(windowL, 'periodic');
%   window = kaiser(windowL,5 );
    noverlap = 0.5*windowL;
    nfft = fStruct.lfpFs;
    Fs = fStruct.lfpFs;
    
    s = spectrogram(lfp,window,noverlap,nfft,Fs);
%   figure; spectrogram(lfp(1,1:10000),window,noverlap,nfft,Fs);
    
    spect = s(2:101, :);
    
end