function [corrData] = MOL_PupilLFP_correlation(spect,fStruct)
%% MOL_PupilLFP_correlation finds correlation between pupil and sum of frequency magnitudes (lfp) for given time lag
% You need in fStruct struct: binnedPupil, binnedFs, freq, lag

pupil = fStruct.binnedPupil;
jj = round(fStruct.lag*fStruct.binnedFs);
kk = jj-1;
pupilShifted = fStruct.binnedPupil(jj:end);
sumFreq0 = sum(abs(spect(fStruct.freq(1):fStruct.freq(2),1:end-kk)));
R = corrcoef(pupilShifted,sumFreq0,'rows','complete');

corrData  = max(R(1,2));
end