%{
Ik heb iig nu een scriptje om met simulaties de AUC van de KS-test, t-test
en zeta te vergelijken, zoals jij ook al had laten zien. Ik vind daar ook
dat de t-test beter is dan zeta wanneer je perfecte overlap hebt van T met
de stimulus tijd, maar dat zeta bij <75% overlap het beter doet. De KS-test
doet het slecht, maar verrassend genoeg opeens beter dan beiden wanneer de
overlap 10% is. Ik ga ook nog een plotje maken wat er gebeurt als je een
mean-rate verschil + piek hebt, dan kan ik een 2D plot maken van %
improvement over t-test met duration-matching en peak size. Ik kan dan ook
nog een versie maken voor de theoretical predictions; hopelijk laat die
iets vergelijkbaars zien met de simulaties. Alles samen geeft dan volgens
mij wel een goed beeld van waar de t-test en zeta tot vergelijkbare
performance leiden; en dat geeft dan ook meteen aan dat het gebied waar de
t-test beter werkt dan zeta eigenlijk alleen bij strict poisson-spiking
neuronen is.

%}

%% shrink stim time (while keeping same number of spikes?)
T = 1; %stim duration
tau = 2; %total trial duration (s+b)
m = 100;%number of trials
intExtraEdge = 8;
intReps = 1000;

dblBaseHz = 1;
vecDiffHz = 0:0.1:1;
vecTr = 1:-0.025:0.5;%1:-0.1:0.1;
intNumHz = numel(vecDiffHz);
intNumT = numel(vecTr);
matTtestP = nan(intNumHz,intNumT,intReps);
matZetaP = nan(intNumHz,intNumT,intReps);
matKsTestP = nan(intNumHz,intNumT,intReps);

matTtestP_FA = nan(intNumHz,intNumT,intReps);
matZetaP_FA = nan(intNumHz,intNumT,intReps);
matKsTestP_FA = nan(intNumHz,intNumT,intReps);


hTic=tic;
tic
for intT=1:intNumT
	%% prep
	Tr = vecTr(intT);
	for intDiffHz=1:intNumHz
		dblDiffHz = vecDiffHz(intDiffHz);
		L_b = dblBaseHz;
		L_s = dblBaseHz + dblDiffHz;
		if toc(hTic) > 5
			fprintf('Processing T=%.3f,Hz=%.3f [%s]\n',Tr,dblDiffHz,getTime);
			hTic = tic;
		end
		parfor intRep=1:intReps
			%% real
			% gen data
			[vecSpikeTimes,vecEventStarts] = getGeneratedTriPhasicR(m+intExtraEdge,Tr,T,tau,L_b,L_s);
			vecEventStarts = vecEventStarts((1+intExtraEdge/2):(end-(intExtraEdge/2)));
			
			% zeta
			intResampNum = 100;
			hTicZ=tic;
			[dblZetaP,dummy,sZETA] = getZeta(vecSpikeTimes,vecEventStarts,tau,intResampNum);
			dblCTZ = toc(hTicZ);
			matZetaP(intDiffHz,intT,intRep) = dblZetaP;
			matZetaCompTime(intDiffHz,intT,intRep) = dblCTZ;
			
			% t-test
			hTicT=tic;
			vecStimHz = getSpikeCounts(vecSpikeTimes,vecEventStarts,vecEventStarts+tau/2);
			vecBaseHz = getSpikeCounts(vecSpikeTimes,vecEventStarts+tau/2,vecEventStarts+tau);
			[h,dblTtestP] = ttest(vecStimHz,vecBaseHz);
			dblCTT = toc(hTicT);
			matTtestP(intDiffHz,intT,intRep) = dblTtestP;
			matTtestCompTime(intDiffHz,intT,intRep) = dblCTT;
			
			% k-s
			hTicKs=tic;
			[vecDelta] = getDelta(vecSpikeTimes,vecEventStarts,tau);
			vecISIs = diff(vecSpikeTimes);
			vecRandSpikeTimes = cat(1,vecSpikeTimes(1),cumsum(vecISIs(randperm(numel(vecISIs)))));
			[vecDeltaRand] = getDelta(vecRandSpikeTimes,vecEventStarts,tau);
			[h,pKs] = kstest2(vecDelta-mean(vecDelta),vecDeltaRand-mean(vecDeltaRand));
			dblCTK = toc(hTicKs);
			matKsTestP(intDiffHz,intT,intRep) = pKs;
			matKsTestCompTime(intDiffHz,intT,intRep) = pKs;
			
			%% shuffled
			% gen data
			L_mean = numel(vecSpikeTimes)/((m+intExtraEdge)*tau);
			%L_mean = (L_s + L_b*(tau-T))/tau;
			%L_mean = (L_s + L_b)/tau;
			[vecSpikeTimes_FA,vecEventStarts_FA]= getGeneratedSpikingPoisson(m+intExtraEdge,T,tau,L_mean,L_mean);
			vecISIs_FA = diff(vecSpikeTimes_FA);
			vecSpikeTimes_FA = cat(1,vecISIs_FA(randi(numel(vecISIs_FA))),cumsum(vecISIs_FA(randperm(numel(vecISIs_FA)))));
			vecEventStarts_FA = vecEventStarts_FA((1+intExtraEdge/2):(end-(intExtraEdge/2)));
			
			% zeta
			intResampNum = 100;
			[dblZetaP_FA,vecLatencies,sZETA_FA] = getZeta(vecSpikeTimes_FA,vecEventStarts_FA,tau,intResampNum);
			matZetaP_FA(intDiffHz,intT,intRep) = dblZetaP_FA;
			% t-test
			vecStimHz_FA = getSpikeCounts(vecSpikeTimes_FA,vecEventStarts_FA,vecEventStarts_FA+tau/2);
			vecBaseHz_FA = getSpikeCounts(vecSpikeTimes_FA,vecEventStarts_FA+tau/2,vecEventStarts_FA+tau);
			[h,dblTtestP_FA] = ttest(vecStimHz_FA,vecBaseHz_FA);
			matTtestP_FA(intDiffHz,intT,intRep) = dblTtestP_FA;
			% k-s
			[vecDelta_FA] = getDelta(vecSpikeTimes_FA,vecEventStarts_FA,tau);
			vecRandSpikeTimes_FA = cat(1,vecISIs_FA(randi(numel(vecISIs_FA))),cumsum(vecISIs_FA(randperm(numel(vecISIs_FA)))));
			[vecDeltaRand_FA] = getDelta(vecRandSpikeTimes_FA,vecEventStarts_FA,tau);
			[h,pKs_FA] = kstest2(vecDelta_FA-mean(vecDelta_FA),vecDeltaRand_FA-mean(vecDeltaRand_FA));
			matKsTestP_FA(intDiffHz,intT,intRep) = pKs_FA;
			
		end
	end
end
toc

%% statistical test AUC
matMeanAUC_Z = nan(intNumHz,intNumT);
matMeanAUC_T = nan(intNumHz,intNumT);
matAUC_TZ_P = nan(intNumHz,intNumT);
matSeAUC_Z = nan(intNumHz,intNumT,2);
matSeAUC_T = nan(intNumHz,intNumT,2);
matMeanAUC_K = nan(intNumHz,intNumT);
matSeAUC_K = nan(intNumHz,intNumT,2);
for intT=1:intNumT
	for intHz=1:intNumHz
		vecP_TP_Z=matZetaP(intHz,intT,:);
		vecP_FP_Z=matZetaP_FA(intHz,intT,:);
		
		vecP_TP_T=matTtestP(intHz,intT,:);
		vecP_FP_T=matTtestP_FA(intHz,intT,:);
		
		vecP_TP_K=matKsTestP(intHz,intT,:);
		vecP_FP_K=matKsTestP_FA(intHz,intT,:);
		
		vecAllT = cat(1,vecP_TP_T(:),vecP_FP_T(:));
		vecAllZ = cat(1,vecP_TP_Z(:),vecP_FP_Z(:));
		vecAllK = cat(1,vecP_TP_K(:),vecP_FP_K(:));
		
		vecClass = cat(1,0*vecP_TP_T(:),0*vecP_FP_T(:)+1);
		[AUC_T,AUC_T_ci,AUC_T_se] = auc([vecClass vecAllT],0.05,'mann-whitney');
		[AUC_Z,AUC_Z_ci,AUC_Z_se] = auc([vecClass vecAllZ],0.05,'mann-whitney');
		[AUC_K,AUC_K_ci,AUC_K_se] = auc([vecClass vecAllK],0.05,'mann-whitney');
		
		% Observed data
		m0 = AUC_T - AUC_Z;
		s0 = (AUC_T_se + AUC_Z_se)/2;
		z = m0/s0;
		AUC_p = 1 - abs(normcdf(z)-normcdf(-z));
		
		%output
		matMeanAUC_Z(intHz,intT) = AUC_Z;
		matMeanAUC_T(intHz,intT) = AUC_T;
		matSeAUC_Z(intHz,intT,:) = AUC_Z_ci;
		matSeAUC_T(intHz,intT,:) = AUC_T_ci;
		matAUC_TZ_P(intHz,intT) = AUC_p;
		
		matMeanAUC_K(intHz,intT) = AUC_K;
		matSeAUC_K(intHz,intT,:) = AUC_K_ci;
	end
end
%% plot
figure
maxfig;
subplot(3,5,6)
vecPlotTr = vecTr(end:-1:1);
imagesc(vecDiffHz,vecTr,matMeanAUC_Z',[0.5 1]);
title(sprintf('ZETA'));
h=colorbar;
xlabel('d(Hz), stim vs ITI');
ylabel('Response duration (s)');
clabel(h,'AUC');
fixfig;
grid off;

subplot(3,5,1)
imagesc(vecDiffHz,vecTr,matMeanAUC_T',[0.5 1]);
title(sprintf('T-test'));
h=colorbar;
clabel(h,'AUC');
xlabel('d(Hz), stim vs ITI');
ylabel('Response duration (s)');
fixfig;
grid off;
%{
subplot(3,5,6)
imagesc(vecDiffHz,vecTr,matMeanAUC_K',[0.5 1]);
title(sprintf('KS-test'));
h=colorbar;
clabel(h,'AUC');
xlabel('d(Hz), stim vs ITI');
ylabel('Response duration (s)');
fixfig;
grid off;
%}
%plot examples
vecExampleTr = [1 1 0.9 0.8];
vecExampleDHz = [0 0.5 0 0];
for intExample=1:4
	Tr = vecExampleTr(intExample);
	dblDiffHz = vecExampleDHz(intExample);
	
	L_b = dblBaseHz;
	L_s = dblBaseHz + dblDiffHz;
		
	% gen data
	[vecSpikeTimes,vecEventStarts] = getGeneratedTriPhasicR(m+intExtraEdge,Tr,T,tau,L_b,L_s);
	vecEventStarts = vecEventStarts((1+intExtraEdge/2):(end-(intExtraEdge/2)));
	dblStep = 0.1;
	vecBins = 0:dblStep:tau;
	h=subplot(3,5,1+intExample);
	[vecMean,vecSEM,vecWindowBinCenters] = doPEP(vecSpikeTimes,vecBins,vecEventStarts,-1);
	hold on
	plot(vecWindowBinCenters,vecMean,'b');
	plot(vecWindowBinCenters,vecMean-vecSEM,'b--');
	plot(vecWindowBinCenters,vecMean+vecSEM,'b--');
	hold off
	ylabel('Firing rate (Hz)');
	xlabel('Time after onset (s)');
	title(sprintf('Tr=%.3f s, r_s=%.1f Hz, r_b=%.1f Hz',Tr,L_s,L_b));
	ylim([0 3]);
	%ylim([0 max(get(gca,'ylim'))]);
	fixfig;
	
	%AUCs
	h=subplot(3,5,6+intExample);
	vecP_TP_Z=sort(matZetaP(vecDiffHz==dblDiffHz,vecTr==Tr,:));
	vecP_FP_Z=sort(matZetaP_FA(vecDiffHz==dblDiffHz,vecTr==Tr,:));
	vecAllPZ = cat(1,vecP_TP_Z(:),vecP_FP_Z(:));
	vecTPRZ = sum(vecAllPZ>=vecP_TP_Z(:)',2)/numel(vecP_TP_Z);
	vecFPRZ = sum(vecAllPZ>=vecP_FP_Z(:)',2)/numel(vecP_FP_Z);
	
	vecP_TP_T=matTtestP(vecDiffHz==dblDiffHz,vecTr==Tr,:);
	vecP_FP_T=matTtestP_FA(vecDiffHz==dblDiffHz,vecTr==Tr,:);
	vecAllPZ = cat(1,vecP_TP_T(:),vecP_FP_T(:));
	vecTPRT = sum(vecAllPZ>=vecP_TP_T(:)',2)/numel(vecP_TP_T);
	vecFPRT = sum(vecAllPZ>=vecP_FP_T(:)',2)/numel(vecP_FP_T);
	
	hold on
	plot(sort(vecFPRZ),sort(vecTPRZ),'b');
	plot(sort(vecFPRT),sort(vecTPRT),'k');
	hold off
	xlabel('FP rate');
	ylabel('TP rate');
	fixfig;
end
return
%% save
export_fig(sprintf('ExpFiringComparisonZetaTtestKs%s.tif',getDate));
export_fig(sprintf('ExpFiringComparisonZetaTtestKs%s.pdf',getDate));


save(sprintf('ExpFiringComparisonZetaTtestKs%s',getDate),'T','tau','m','intExtraEdge','dblBaseHz','vecDiffHz','vecTr',...
	'matTtestP','matZetaP','matKsTestP','matTtestP_FA','matZetaP_FA','matKsTestP_FA');


%% rest
return
%subplot(2,3,5)
[M,c] = contour(vecDiffHz,vecTr,matMeanAUC_T',0.6:0.1:0.9);
clabel(M,c,'Color','black');
axis ij

%subplot(2,3,6)
[M,c] = contour(vecDiffHz,vecTr,matMeanAUC_K',0.6:0.1:0.9);
clabel(M,c,'Color','red');
axis ij
return
figure
hold on
errorbar(vecTr,matMeanAUC_Z,matSeAUC_Z(1,:)-matMeanAUC_Z,matSeAUC_Z(2,:)-matMeanAUC_Z);
errorbar(vecTr,matMeanAUC_T,matSeAUC_T(1,:)-matMeanAUC_T,matSeAUC_T(2,:)-matMeanAUC_T,'k');
errorbar(vecTr,matMeanAUC_K,matSeAUC_K(1,:)-matMeanAUC_K,matSeAUC_K(2,:)-matMeanAUC_K,'r');
hold off
legend({'Zeta','T-test','KS-test'},'Location','Best');
strTitP = sprintf('%.1f=%.3f;',flat(cat(1,vecTr,matAUC_TZ_P)));
title(strTitP)
xlabel('Response duration (s)');
ylabel('AUC');
ylim([0.4 1]);
fixfig
return
%% %%
%% sliding window
%% %%

%% shrink stim time (while keeping same number of spikes?)
vecShiftT = -0.5:0.1:0.5;
intNumShiftT = numel(vecShiftT);
vecRev=1;
intRevN = numel(vecRev);
intReps=100;
matTtestP = nan(intReps,intNumShiftT,intRevN);
matZetaP = nan(intReps,intNumShiftT,intRevN);
matKsTestP = nan(intReps,intNumShiftT,intRevN);
matTtestP_FA = nan(intReps,intNumShiftT);
matZetaP_FA = nan(intReps,intNumShiftT);
matKsTestP_FA = nan(intReps,intNumShiftT);
hTic=tic;
tic

T = tau/2;
for intT=1:intNumShiftT
	%% prep
	dblShiftT = vecShiftT(intT);
	if toc(hTic) > 5
		fprintf('Processing shiftT=%.1f [%s]\n',dblShiftT,getTime);
		hTic = tic;
	end
	
	parfor intRep=1:intReps
		%% real
		% gen data
		[vecSpikeTimes,vecEventStarts]= getGeneratedSpikingPoisson(m+intExtraEdge,T,tau,L_b,L_s);
		vecEventStarts = vecEventStarts((1+intExtraEdge/2):(end-(intExtraEdge/2)));
		dblEndT = (m+intExtraEdge)*tau;
		
		for intReverse=vecRev
			if intReverse == 2
				vecEventStarts = sort(dblEndT - vecEventStarts);
				vecSpikeTimes = sort(dblEndT - vecSpikeTimes);
				vecEventStartsShifted = vecEventStarts + dblShiftT;
			else
				vecEventStartsShifted = vecEventStarts + dblShiftT;
			end
			% zeta
			intResampNum = 100;
			[dblZetaP,dummy,sZETA] = getZeta(vecSpikeTimes,vecEventStartsShifted,tau,intResampNum);
			matZetaP(intRep,intT,intReverse) = dblZetaP;
			% t-test
			vecStimHz = getSpikeCounts(vecSpikeTimes,vecEventStartsShifted,vecEventStartsShifted+T);
			vecBaseHz = getSpikeCounts(vecSpikeTimes,vecEventStartsShifted+T,vecEventStartsShifted+tau);
			[h,dblTtestP] = ttest(vecStimHz,vecBaseHz);
			matTtestP(intRep,intT,intReverse) = dblTtestP;
			% k-s
			%[h,pKs] = kstest2(sZETA.vecD,sZETA.matRandD(:));
			%matKsTestP(intRep,intT) = pKs;
			% k-s
			[vecDelta] = getDelta(vecSpikeTimes,vecEventStartsShifted,tau);
			vecISIs = diff(vecSpikeTimes);
			vecRandSpikeTimes = cat(1,vecSpikeTimes(1),cumsum(vecISIs(randperm(numel(vecISIs)))));
			[vecDeltaRand] = getDelta(vecRandSpikeTimes,vecEventStartsShifted,tau);
			%calc difference
			[h,pKs] = kstest2(vecDelta-mean(vecDelta),vecDeltaRand-mean(vecDeltaRand));
			matKsTestP(intRep,intT,intReverse) = pKs;
			%[h,pKs2] = kstest2(sZETA.vecD,sZETA.matRandD(:));
		end
		%error check what happens when we reverse spike times
		
		%% shuffled
		% gen data
		L_mean = (L_s*T + L_b*(tau-T))/tau;
		[vecSpikeTimes_FA,vecEventStarts_FA]= getGeneratedSpikingPoisson(m+intExtraEdge,T,tau,L_mean,L_mean);
		vecISIs_FA = diff(vecSpikeTimes_FA);
		vecSpikeTimes_FA = cat(1,vecISIs_FA(randi(numel(vecISIs_FA))),cumsum(vecISIs_FA(randperm(numel(vecISIs_FA)))));
		vecEventStarts_FA = vecEventStarts_FA((1+intExtraEdge/2):(end-(intExtraEdge/2)));
		vecEventStarts_FA = vecEventStarts_FA + dblShiftT;
		% zeta
		intResampNum = 100;
		[dblZetaP_FA,vecLatencies,sZETA_FA] = getZeta(vecSpikeTimes_FA,vecEventStarts_FA,tau,intResampNum);
		matZetaP_FA(intRep,intT) = dblZetaP_FA;
		% t-test
		vecStimHz_FA = getSpikeCounts(vecSpikeTimes_FA,vecEventStarts_FA,vecEventStarts_FA+1);
		vecBaseHz_FA = getSpikeCounts(vecSpikeTimes_FA,vecEventStarts_FA+1,vecEventStarts_FA+2);
		[h,dblTtestP_FA] = ttest(vecStimHz_FA,vecBaseHz_FA);
		matTtestP_FA(intRep,intT) = dblTtestP_FA;
		
		% k-s
		[vecDelta_FA] = getDelta(vecSpikeTimes_FA,vecEventStarts_FA,tau);
		vecISIs_FA = diff(vecSpikeTimes_FA);
		vecRandSpikeTimes_FA = cat(1,vecSpikeTimes_FA(1),cumsum(vecISIs_FA(randperm(numel(vecISIs_FA)))));
		[vecDeltaRand_FA] = getDelta(vecRandSpikeTimes_FA,vecEventStarts_FA,tau);
		
		%calc difference
		[h,pKs_FA] = kstest2(vecDelta_FA-mean(vecDelta_FA),vecDeltaRand_FA-mean(vecDeltaRand_FA));
		matKsTestP_FA(intRep,intT) = pKs_FA;
		
		%[h,pKs2_FA] = kstest2(sZETA_FA.vecD,sZETA_FA.matRandD(:));
		
	end
end
toc

%% statistical test AUC
matMeanAUC_Z = nan(intNumShiftT,intRevN);
matMeanAUC_T = nan(intNumShiftT,intRevN);
matAUC_TZ_P = nan(intNumShiftT,intRevN);
matSeAUC_Z = nan(2,intNumShiftT,intRevN);
matSeAUC_T = nan(2,intNumShiftT,intRevN);
matMeanAUC_K = nan(intNumShiftT,intRevN);
matSeAUC_K = nan(2,intNumShiftT,intRevN);
for intT=1:intNumShiftT
	for intRev=vecRev
		vecP_TP_Z=matZetaP(:,intT,intRev);
		vecP_TP_T=matTtestP(:,intT,intRev);
		vecP_TP_K=matKsTestP(:,intT,intRev);
		
		vecP_FP_Z=matZetaP_FA(:,intT);
		vecP_FP_T=matTtestP_FA(:,intT);
		vecP_FP_K=matKsTestP_FA(:,intT);
		
		vecAllT = cat(1,vecP_TP_T,vecP_FP_T);
		vecAllZ = cat(1,vecP_TP_Z,vecP_FP_Z);
		vecAllK = cat(1,vecP_TP_K,vecP_FP_K);
		
		vecClass = cat(1,0*vecP_TP_T,0*vecP_FP_T+1);
		[AUC_T,AUC_T_ci,AUC_T_se] = auc([vecClass vecAllT],0.05,'mann-whitney');
		[AUC_Z,AUC_Z_ci,AUC_Z_se] = auc([vecClass vecAllZ],0.05,'mann-whitney');
		[AUC_K,AUC_K_ci,AUC_K_se] = auc([vecClass vecAllK],0.05,'mann-whitney');
		
		% Observed data
		m0 = AUC_T - AUC_Z;
		s0 = (AUC_T_se + AUC_Z_se)/2;
		z = m0/s0;
		AUC_p = 1 - abs(normcdf(z)-normcdf(-z));
		
		%output
		matMeanAUC_Z(intT,intRev) = AUC_Z;
		matMeanAUC_T(intT,intRev) = AUC_T;
		matSeAUC_Z(:,intT,intRev) = AUC_Z_ci;
		matSeAUC_T(:,intT,intRev) = AUC_T_ci;
		matAUC_TZ_P(intT,intRev) = AUC_p;
		
		matMeanAUC_K(intT,intRev) = AUC_K;
		matSeAUC_K(:,intT,intRev) = AUC_K_ci;
	end
end

%% move stim time
figure
hold on
for intRev=vecRev
	if intRev==2
		strEx = 'x--';
	else
		strEx = '';
	end
	errorbar(vecShiftT,matMeanAUC_Z(:,intRev),flat(matSeAUC_Z(1,:,intRev))-matMeanAUC_Z(:,intRev),flat(matSeAUC_Z(2,:,intRev))-matMeanAUC_Z(:,intRev),['b' strEx]);
	errorbar(vecShiftT,matMeanAUC_T(:,intRev),flat(matSeAUC_T(1,:,intRev))-matMeanAUC_T(:,intRev),flat(matSeAUC_T(2,:,intRev))-matMeanAUC_T(:,intRev),['k' strEx]);
	errorbar(vecShiftT,matMeanAUC_K(:,intRev),flat(matSeAUC_K(1,:,intRev))-matMeanAUC_K(:,intRev),flat(matSeAUC_K(2,:,intRev))-matMeanAUC_K(:,intRev),['r' strEx]);
end
hold off
legend({'Zeta','T-test','KS-test'},'Location','Best');
%strTitP = sprintf('%.1f=%.3f;',flat(cat(1,vecShiftT,vecAUC_TZ_P)));
%title(strTitP)
xlabel('Response shift (s)');
ylabel('AUC');
ylim([0.4 1]);
fixfig