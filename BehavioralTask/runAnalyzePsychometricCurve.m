%runAnalyzePsychometricCurve

%define variables
strSes = '20140207';
strMasterPath = 'D:\Data\Processed\imagingdata';
vecRecordings = 1:8;

%put in structure
sIn.strSes = strSes;
sIn.strMasterPath = strMasterPath;
sIn.vecRecordings = vecRecordings;

%get stim aggregate
sStimAggregate = buildStimAggregate(sIn);

%calculate psychometric curve on % responded and RT
sOut = calcPsychometricCurve(sStimAggregate);

%plot
doPlotPsychoMetricCurve(sStimAggregate,sOut);