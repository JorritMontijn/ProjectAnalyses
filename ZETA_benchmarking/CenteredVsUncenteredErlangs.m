%% build erlang from mean-centered delta-dots and compare with non-mean-centered delta-dots erlang

L_s=10; %rate during stimulus
L_b=10; %rate during baseline
T = 1; %stim duration
tau = 2; %total trial duration (s+b)
m = 100;%number of trials
intReps=100;

n_s = m*T*L_s;
n_b = m*(tau-T)*L_b;
n=n_s + n_b;

intSavePoints = round(n);
matDeltadot_i = nan(intReps,intSavePoints);
matDelta_i = nan(intReps,intSavePoints);
matDeltaU_i = nan(intReps,intSavePoints);

for intRep=1:intReps
%d
[vecSpikeTimes,vecEventStarts]= getGeneratedSpikingPoisson(m,T,tau,L_b,L_s);


% prepare interpolation points
dblUseMaxDur = tau;
vecSpikeT = getSpikeT(vecSpikeTimes,vecEventStarts,dblUseMaxDur);

% run normal
%pre-allocate
vecStimUseOnTime = vecEventStarts(:,1);
vecThisSpikeTimes = unique(getSpikeT(vecSpikeTimes,vecStimUseOnTime,dblUseMaxDur));
vecThisSpikeFracs = linspace(0,1,numel(vecThisSpikeTimes))';
vecThisFrac = interp1(vecThisSpikeTimes,vecThisSpikeFracs,vecSpikeT);

%get linear fractions
vecThisFracLinear = (vecSpikeT./dblUseMaxDur);

%calc difference
vecDelta = vecThisFrac - vecThisFracLinear;
vecDeltaDot = diff(vecDelta);
vecDelta = [0; cumsum(vecDeltaDot)];
varDeltaDotEmp = var(vecDeltaDot);
varDeltaDotThe = 1/n^2;

%% delta_i is not erlang-distributed over repetitions, but S_i is (or should be).
n_e = numel(vecThisSpikeTimes);
S = ((1/(n_e-1)) - vecDeltaDot); %v_i normalized to [0 1]

varDeltaDotEmp2 = nanvar(S); %variance is the same, but mean of S is 1/n

vecDelta_uncentered = [0; cumsum(S)]; %v_i normalized to [0 1] => vecThisFracLinear

%matDeltadot_i(intRep,1:intSavePoints) = S(1:intSavePoints); %all iid exponentials

vecDelta = interp1(linspace(0,1,n_e),vecDelta,linspace(0,1,n));
vecDelta_uncentered = interp1(linspace(0,1,n_e),vecDelta_uncentered,linspace(0,1,n));
matDelta_i(intRep,1:intSavePoints) = vecDelta(1:intSavePoints);
matDeltaU_i(intRep,1:intSavePoints) = vecDelta_uncentered(1:intSavePoints);
end

% predict brownian bridge
t=1:n;
W = (1:n)./(n^2);
W_sub=(t./(n))*(1/n);
%x=((n-t)./sqrt(n))*W_alt

%B Klar (2014/2015), gamma difference distributions:
%?=?1/?1 - ?2/?2;
%?^2=?1/?1^2 + ?2/?2^2;

%moreover,
%if X~Gamma(a,b)
%then cX~Gamma(a,b/c)
%
vecE = nan(1,n);
vecS = nan(1,n);

vecS_diff = nan(1,n);
vecE_diff = nan(1,n);
vecM_T = nan(1,n);
vecV_T = nan(1,n);
for intI=1:n
i = intI;
a1 = i;%shape (i)
b1 = 1/n;%rate (lambda)
deltaU_i_reps = gamrnd(a1,b1,[1 intReps]);

vecE(intI) = mean(deltaU_i_reps);
vecS(intI) = var(deltaU_i_reps);

B1 = (n-i)/(n);
B2 = 1-B1;%4*(0.25-((n/2-i+1)/n)^2);

a2 = n-i;%shape (i)
b2 = 1/n;%1/i;rate (lambda)
brown_rem = gamrnd(a2,b2,[1 intReps]);

delta_i_reps = B1*deltaU_i_reps-B2*brown_rem;
%delta_i_reps = 2.5*(1-dblBrownFactor)*brown_rem;
vecE_diff(intI) = mean(delta_i_reps);
vecS_diff(intI) = var(delta_i_reps);

vecM_T(intI) = B1*a1*b1 - B2*a2*b2;
vecV_T(intI) = (B1*(a1/((1/b1)^2)) + B2*(a2/((1/b2)^2)))/2;
end

%
figure
subplot(2,3,1)
%imagesc(matDelta_i-nanmean(matDelta_i))
plot(matDelta_i(1:10,:)')
ylabel('delta_i');
xlabel('i');
fixfig;

subplot(2,3,2)
%imagesc(matDelta_i-nanmean(matDelta_i))
hold on
plot(vecE_diff)
plot(nanmean(matDelta_i))
plot(vecM_T)
hold off
ylabel('E(delta_i)');
xlabel('i');
legend({'Sampled','Emp., 1 run','Theory'});
fixfig;

subplot(2,3,3)
hold on
plot(vecS_diff)
plot(nanvar(matDelta_i))
plot(vecV_T)
hold off
ylabel('Var(delta_i)');
xlabel('i');
legend({'Sampled','Emp., 1 run','Theory'});
fixfig;

subplot(2,3,4)
%imagesc(matDelta_i-nanmean(matDelta_i))
plot(matDeltaU_i(1:10,:)')
ylabel('deltaU_i');
xlabel('i');
fixfig;

subplot(2,3,5)
%imagesc(matDeltaU_i-nanmean(matDeltaU_i))
hold on
plot(vecE)
plot(nanmean(matDeltaU_i))
hold off
ylabel('E(deltaU_i)');
xlabel('i');
legend({'Sampled','Emp., 1 run'});
fixfig;

subplot(2,3,6)
hold on
plot(vecS)
plot(nanvar(matDeltaU_i))
hold off
ylabel('Var(deltaU_i)');
xlabel('i');
legend({'Sampled','Emp., 1 run'});
fixfig;
