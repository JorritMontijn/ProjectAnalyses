%% description
%{
    - The lack of a 2nd bump in the visual miss trials was expected. I find
    it very interesting that there may also be a difference between hits
    and misses for what pertains the first bump. Do you think it would be
    possible to look at B1_hit / B1_miss alignment? Of course we could also
    try a simpler approach (just look at whether the trajectories in a
    latent space differ between hits and misses).

When you look at the absolute numbers (in the non-normalized figures) in
figures 3B/4B, there actually doesn't seem to be that much difference in B1
dynamics between hits and misses: it's mostly the presence/absence of B2
that changes the relative size of the B1 trajectories. But of course the
B1-B2 are different logistic regressors, so we can't say for certain based
on these graphs.     
I think the right analysis to answer your question here would be to do a
reduced rank regression like before to compare the performance of
B1_hit(B1_hit)->B1_hit with B1_hit(B1_miss)->B1_hit (and the other way
around). This would allow us to quantify the difference in neural subspaces
between B1 hits and misses.
%}

%% clear
clear all;

%% set variables
%On which timestamp to align as t=0
strAlignOn = 'Change';
intFullIters = 10;
intResamplings = 10;
boolShuffle = false;
dblLambda = 1/10;
boolSavePlots = true;
strFigDir = 'D:\Data\Results\BumpsMatthijs\';

%get neuronal activity from t=0 - t=0.4
dblBinSizeSecs = 100/1000;%25/1000;
vecBinEdges = [0 dblBinSizeSecs];
dblBinOffsetStep = (dblBinSizeSecs/10);
vecBinOffsets = -2:(dblBinSizeSecs/10):2;
vecBinStartT = vecBinOffsets;
intBins = numel(vecBinStartT);
vecRemoveSessions = [3 6 9]; %due to bad performance

%% HEADER, load data
MOL_Header;

%% analysis
%is B1 activity hit-specific?
%predict B1_hit(B1_hit)->B1_hit and compare with B1_hit(B1_miss)->B1_hit => x% generalizable
%is B1 hit-specificity larger than B2 hit-specificity?
%predict B2_hit(B2_hit)->B2_hit and compare with B2_hit(B2_miss)->B2_hit  => y% generalizable
%x < y? p=?

