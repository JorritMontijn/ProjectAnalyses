%% set paths
strPath2 = 'C:\Code\Toolboxes\zetatest\dependencies\';
strPath1 = 'C:\Code\Acquisition\UniversalProbeFinder\zetatest\dependencies\';

%% gen data
dblTau = 2;
dblNoise = 0.025;
dblSamplingFreq = 25;
boolQuick = false;

%set indicator properties
sIndicatorProps = struct;
sIndicatorProps.dblTimescale = dblTau;
sIndicatorProps.dblNoise = dblNoise;

% stimulus data
dblFactor1=50;
dblFactor2=50;
dblBaseRate = exprnd(0.1)+0.1;
vecDurs = [0.1 0.9 0.1];
vecRates = [dblFactor1*exprnd(1)+dblFactor1/10 exprnd(0.2)+0.1 dblFactor2*exprnd(1)+dblFactor2/10];
intNumT = 100;
vecTrialDur=linspace(0.5,10,intNumT);
vecRepStarts = 5+cumsum(vecTrialDur);
dblEndT = vecRepStarts(end)+5;
vecSpikeTimes = getGeneratedMultiPhasicR(dblBaseRate,vecRates,vecDurs,vecRepStarts,dblEndT);

dblUseMaxDur = 1;
vecTrialStarts = vecRepStarts(:);
matTrialT1 = cat(2,vecTrialStarts,vecTrialStarts+dblUseMaxDur);
matTrialT1 = matTrialT1 + 0.05*rand(size(matTrialT1));

% generate dfof
[vecTimestamps,vecdFoF] = getGeneratedFluorescence(vecSpikeTimes,dblSamplingFreq,sIndicatorProps,boolQuick);

%% test
vecEventStartT = matTrialT1(:,1);
%test
%vecTimestamps = [0.1:0.011:7];
%vecdFoF = ones(size(vecTimestamps));
%vecEventStartT = [1 3 5];
%dblUseMaxDur = 1;

% getTsRefT.m
cd(strPath1)
rng(1,'mt19937ar')
vecRefT = getTsRefT(vecTimestamps,vecEventStartT,dblUseMaxDur);
cd(strPath2)
rng(1,'mt19937ar')
vecRefT2 = getTsRefT(vecTimestamps,vecEventStartT,dblUseMaxDur);
sum(abs(vecRefT-vecRefT2))


% getInterpolatedTimeSeries.m
cd(strPath1)
rng(1,'mt19937ar')
[vecRefT2_a,matTracePerTrial] = getInterpolatedTimeSeries(vecTimestamps,vecdFoF,vecEventStartT,dblUseMaxDur,vecRefT);
cd(strPath2)
rng(1,'mt19937ar')
[vecRefT2_b,matTracePerTrial2] = getInterpolatedTimeSeries(vecTimestamps,vecdFoF,vecEventStartT,dblUseMaxDur,vecRefT2);
sum(abs(matTracePerTrial(:)-matTracePerTrial2(:)))

%calcTsZetaOne
intResampNum = 100;boolDirectQuantile=false;dblJitterSize=2;boolUseParallel=false;intJitterDistro=2;
cd(strPath1)
rng(1,'mt19937ar')
tic
[vecRefT,vecRealDiff,vecRealFrac,vecRealFracLinear,cellRandT,cellRandDiff,dblZetaP,dblZETA,intZETALoc] = ...
 	calcTsZetaOne(vecTimestamps,vecdFoF,vecEventStartT,dblUseMaxDur,intResampNum,boolDirectQuantile,dblJitterSize,intJitterDistro);
toc
cd(strPath2)
rng(1,'mt19937ar')
tic
[vecRefT2,vecRealDiff2,vecRealFrac2,vecRealFracLinear2,cellRandT2,cellRandDiff2,dblZetaP2,dblZETA2,intZETALoc2] = ...
 	calcTsZetaOne(vecTimestamps,vecdFoF,vecEventStartT,dblUseMaxDur,intResampNum,boolDirectQuantile,dblJitterSize,boolUseParallel);
toc
any(cellfun(@numel,cellRandT)~=cellfun(@numel,cellRandT2))

% getTraceOffsetOne.m
cd(strPath1)
rng(1,'mt19937ar')
[vecRandDiff,vecThisFrac,vecThisFracLinear,vecRandT] = getTraceOffsetOne(vecTimestamps,vecdFoF,vecEventStartT,dblUseMaxDur);
cd(strPath2)
rng(1,'mt19937ar')
[vecRandDiff2,vecThisFrac2,vecThisFracLinear2,vecRandT2] = getTraceOffsetOne(vecTimestamps,vecdFoF,vecEventStartT,dblUseMaxDur);
sum(abs(vecRandDiff-vecRandDiff2))
sum(abs(vecThisFrac-vecThisFrac2))
sum(abs(vecThisFracLinear-vecThisFracLinear2))
sum(abs(vecRandT-vecRandT2))

%% test interp
vecRefT = getTsRefT(vecTimestamps,vecEventStartT,dblUseMaxDur);
% get original times
for intTrial=1:numel(vecEventStartT)
dblStartT = vecEventStartT(intTrial);
intStartT = max([1 find(vecTimestamps > (dblStartT + vecRefT(1)),1) - 1]);
intStopT = min([numel(vecTimestamps) find(vecTimestamps > (dblStartT + vecRefT(end)),1) + 1]);
vecSelectSamples = intStartT:intStopT;

%% get data
vecUseTimes = vecTimestamps(vecSelectSamples);
vecUseTrace = vecdFoF(vecSelectSamples);
vecUseInterpT = vecRefT+dblStartT;
		
vecInterpTrace1 = interp1(vecUseTimes,vecUseTrace,vecUseInterpT);
F = griddedInterpolant(vecUseTimes,vecUseTrace,'linear');
vecInterpTrace2 = F(vecUseInterpT);

all(vecInterpTrace1==vecInterpTrace2)
end