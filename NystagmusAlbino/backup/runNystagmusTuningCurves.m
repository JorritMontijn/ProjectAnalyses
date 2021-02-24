%% load spiking data & plot tuning curves
%{
%20190315:
2 (loc1, cort): some possibly visual responses, not very convincing
4 (loc1, subcort), awesome cells, very strong visual responses
7 (loc1, subcort), several strongly visual cells
10 (loc2, cort), multiple strongly visual cells
12 (loc2, subcort, end error), possibly some visual cells, not very strong
%}

%% set recording
clear all;close all;
strSelectMouseType = 'WT';
strSelectArea = 'V1';
strSelectStim = 'DG'; %DG=drifting grating, RF=receptive field mapping, RF-DG=merged block
intSelectPopulation	= 1;

%% load data
strPath = 'D:\Data\Processed\ePhys\DriftingGratings\';
sFiles = dir([strPath '*.mat']);
for intFile=1:numel(sFiles)
sLoad = load(sFiles(intFile).name);
end
sLoad.sNeuron
sLoad.sRecording