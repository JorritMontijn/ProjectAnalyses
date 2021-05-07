%% params
m = 100;
intExtraEdge = 2;
Tr=0.8;
T=1;
tau=2;
vecHz = 1:10;
intHzNum = numel(vecHz);
intReps = 100;

%% run
matTtestCompTime = nan(intHzNum,intReps);
matZetaCompTime = nan(intHzNum,intReps);
matTtestP = nan(intHzNum,intReps);
matZetaP = nan(intHzNum,intReps);
hTic = tic;
for intHzIdx=1:intHzNum
	L_b=vecHz(intHzIdx);
	L_s =L_b;
	m=100;
	if toc(hTic) > 5
		fprintf('Processing %d/%d [%s]\n',intHzIdx,intHzNum,getTime);
		hTic = tic;
	end
	
	for intRep=1:intReps
		% gen data
		[vecSpikeTimes,vecEventStarts] = getGeneratedTriPhasicR(m+intExtraEdge,Tr,T,tau,L_b,L_s);
		vecEventStarts = vecEventStarts((1+intExtraEdge/2):(end-(intExtraEdge/2)));
		
		% zeta
		intResampNum = 100;
		hTicZ=tic;
		[dblZetaP,dummy,sZETA] = getZeta(vecSpikeTimes,vecEventStarts,tau,intResampNum);
		dblCTZ = toc(hTicZ);
		matZetaP(intHzIdx,intRep) = dblZetaP;
		matZetaCompTime(intHzIdx,intRep) = dblCTZ;
		
		% t-test
		hTicT=tic;
		vecStimHz = getSpikeCounts(vecSpikeTimes,vecEventStarts,vecEventStarts+tau/2);
		vecBaseHz = getSpikeCounts(vecSpikeTimes,vecEventStarts+tau/2,vecEventStarts+tau);
		[h,dblTtestP] = ttest(vecStimHz,vecBaseHz);
		dblCTT = toc(hTicT);
		matTtestP(intHzIdx,intRep) = dblTtestP;
		matTtestCompTime(intHzIdx,intRep) = dblCTT;
	end
end

%% plot
vecBins = []

figure
matLowHighT = getCI(matTtestCompTime,2,0.01,true);
vecMeanCTT = median(matTtestCompTime,2);
matLowHighZ = getCI(matZetaCompTime,2,0.01,true);
vecMeanCTZ = median(matZetaCompTime,2);

subplot(2,3,1)

hold on
plot(vecHz,vecMeanCTZ./vecMeanCTT,'b');
plot(vecHz,matLowHighZ(:,1)./matLowHighT(:,1),'b--');
plot(vecHz,matLowHighZ(:,2)./matLowHighT(:,2),'b--');
hold off
%set(gca,'yscale','log')

subplot(2,3,2)

hold on
plot(vecHz,vecMeanCTZ,'b');
plot(vecHz,vecMeanCTT*10,'k');
hold off