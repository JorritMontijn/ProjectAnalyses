%% set variables
% parameters for window of interest:
dblStartBaseT = -0.4;
dblStopBaseT = -0.2;
dblStartEp1T = 0;
dblStopEp1T = 0.2;
dblStartEp2T = 0.2;
dblStopEp2T = 0.4;

%% load data
strFigDir = 'F:\Data\Results\BumpsMatthijs\';
fprintf('Loading data... [%s]\n',getTime);
[Data] = MOL_GetData('F:\Data\Processed\BumpsMatthijs','CHDET',{'ChangeDetectionConflict'},{'2003' '2009' '2010' '2011' '2012' '2013' '2019' '2020' '2021'},[],{'sessionData' 'trialData' 'spikeData'});
sessionData     = Data.sessionData;
trialData       = Data.trialData;
spikeData       = Data.spikeData;
