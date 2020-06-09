function [fStruct, transLFP, transTs] = MOL_PupilLFP_binning(pupil, lfp, fStruct)
%% MOL_PupilLFP_binning bin pupil data so it will match amount of time-frequency bins
% You need in fStruct struct: pupilFs, lfpFs, binL

pupilFs     = [fStruct.pupilFs];    
lfpFs       = [fStruct.lfpFs];

%take exactly the same amount of sec from pupil and lfp data
lfp0         = lfp;
lfpL         = floor(((length(lfp0(1,:))/lfpFs)))*lfpFs;

if (lfpL/lfpFs)*pupilFs > length(pupil)
    lfpL = lfpL-lfpFs;
end

%lfp = lfp0(:,1:lfpL);

fStruct.lfpL = lfpL;

pupil = pupil(1, 1: (lfpL/lfpFs)*pupilFs);
%pupil_ts_nlx = pupil_ts_nlx(1, 1: (lfpL/lfpFs)*pupilFs);
pupilL = length(pupil);

%% loop for low frequency (1) and high frequency (2) magnitude
    
    binL  = fStruct.binL;
    
    %create vectors that divide the data into chunks of the lenght definded in
    %Params
    
    vecPupil      = 1: binL*pupilFs: pupilL;
    vecSignal     = 1: binL*lfpFs: lfpL;
    
    chunksAmt   = length(vecPupil)-1; %how many chunks
    
% generate vector with time stamps for new sampling frequency
% to activate it add vector with lfp times tamps as an input (lfpTs)
transP = zeros(1,chunksAmt);

% transTs = zeros(1,chunksAmt);
% transLFP = zeros(1,chunksAmt);

% tsTemp = cell2mat(fStruct.lfpTs);
% tsTemp = tsTemp(1,:);

warning('off','all')

for jj = 1:chunksAmt
    binnedP(jj) = nanmean(pupil(vecPupil(jj):vecPupil(jj+1)));
         binnedLFP(jj) = nanmean(lfp(1, vecSignal(jj):vecSignal(jj+1)));
%         transTs(jj) = tsTemp(1, vecSignal(jj+1)-((vecSignal(jj+1)-vecSignal(jj))/2));
end

fStruct.binnedPupil = {binnedP};
end

