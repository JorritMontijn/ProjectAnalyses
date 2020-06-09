function [lagmax, lagmin] = MOL_PupilLFP_correlationLAG(spectrum,fStruct)

% You need in fStruct struct: binnedPupil, binnedFs, check_lag_sec, freq

pupil       = fStruct.binnedPupil;

xL          = round(fStruct.binnedFs*fStruct.checkLag_sec);
Rmat        = zeros(1,xL);
Pmat        = zeros(1,xL);

for jj = 1:xL
    kk = jj-1;
    pupilShifted = pupil(jj:end);
    sumFreq0 = sum(abs(spectrum(fStruct.freq(1):fStruct.freq(2),1:end-kk)));
    
    [R,P]=corrcoef(pupilShifted,sumFreq0,'rows','complete');
    Rmat(jj) = R(2);
    Pmat(jj) = P(2);
end

%find lowest and highest correlation and corresponding time lag

lagmax = find(Rmat == max(Rmat))*fStruct.binL;
lagmin = find(Rmat == min(Rmat))*fStruct.binL;

end