%% normal
mu = 10;
sd = 3;
n = 10;
figure;
hold on
dblStep=10;
vecN = dblStep:dblStep:(dblStep*10);
intNumN = numel(vecN);
meanEmp = nan(1,intNumN);
meanThe = nan(1,intNumN);
meanElf = nan(1,intNumN);
for intN=1:intNumN
	n = vecN(intN);
reps = 1000;
vecMaxX = nan(1,reps);
for i=1:reps
	x=normrnd(mu,sd,[1 n]);
	vecMaxX(i) = max(x);
end
fQuantile = @(mu,sd,p) mu + sd*sqrt(2)*erfinv(2*p - 1);
%p=1 - (1/(sqrt(2)*n+1)); %why sqrt(2)?
p=(n-alpha)./(n-2*alpha+1);
dblPredMaxX = fQuantile(mu,sd,p);
bplot(vecMaxX,n);
plot([-0.1 0.1]*dblStep+n,dblPredMaxX*[1 1],'b--')

%elfving
alpha = pi/8;
E_maxX = norminv((n-alpha)./(n-2*alpha+1)); %expected maximum of a single d' sample with X~N(0,1)
E_maxDprime = mu + sd.*E_maxX;

plot([-0.1 0.1]*dblStep+n,E_maxDprime*[1 1],'r--')

meanEmp(intN) = mean(vecMaxX);
meanThe(intN) = mean(dblPredMaxX);
meanElf(intN) = mean(E_maxDprime);
end
mean(abs(meanEmp-meanElf))
mean(abs(meanEmp- meanThe))
%% exponential
lambda=13;
n = 10;
figure;
hold on
dblStep=10;
for n=dblStep:dblStep:100
reps = 1000;
vecMaxX = nan(1,reps);
for i=1:reps
	x=exprnd(lambda,[1 n]);
	vecMaxX(i) = max(x);
end
fQuantile = @(lambda,p) -log(1-p)*lambda;
%p=1 - (1/(sqrt(2)*n+1)); %why sqrt(2)?
p=(n-alpha)./(n-2*alpha+1);
dblPredMaxX = fQuantile(lambda,p);
bplot(vecMaxX,n);
plot([-0.1 0.1]*dblStep+n,dblPredMaxX*[1 1],'b--')
plot([-0.1 0.1]*dblStep+n,mean(vecMaxX)*[1 1],'k--')
end

%mean(abs(meanEmp- meanThe))


%% beta
n = 1000;
k = 0:n; %# of successes
p = 0.5;

%dblA = k + 1;
%dblB = n - k + 1;

dblA = 2;
dblB = 2;
X = 0:0.01:1;
%pdf
y1=betapdf(X,dblA,dblB);
y2=((X.^(dblA-1)).*((1-X).^(dblB-1)))/beta(dblA,dblB);

%cdf
I = betainc(X,dblA,dblB); %already regularized
%B = beta(dblA,dblB);
Y1 = I;% ./ B;
Y2 = betacdf(X,dblA,dblB);

%inv cdf (q?)
vecP=0:0.01:1;
Yprime1 = betainv(vecP,dblA,dblB);


figure;
hold on
dblStep=100;
for n=dblStep:dblStep:(dblStep*10)
reps = 1000;
vecMaxX = nan(1,reps);
for i=1:reps
	x=betarnd(dblA,dblB,[1 n]);
	vecMaxX(i) = max(x);
end
%p=1 - (1/(n+1)); %why sqrt(2)?
alpha = pi/8;
%p=(n -(n/4) -alpha)./(n -(n/4) -2*alpha+1); %why n - n/4?
p=(n -alpha)./(n -2*alpha+1); %why n - n/4?
dblPredMaxX =  betainv(p,dblA,dblB);
bplot(vecMaxX,n);
plot([-0.1 0.1]*dblStep+n,dblPredMaxX*[1 1],'b--')
plot([-0.1 0.1]*dblStep+n,mean(vecMaxX)*[1 1],'k--')

%binom
k = dblA - 1;
n = dblA + dblB;

X2 = betapdf(p,dblA,dblB)
X3 = binopdf(k,n,p)/(n+1)
end

%error try product of quantile functions to approximate E[max(d')]

%%
n = 1000; %number of spikes
vecK=0:n;
vecP = vecK(2:(end-1))./n;
vecI = 1:(n-1);

vecX = vecP;
subplot(2,3,1);cla;
hold on
vecProbK = log(ones(size(vecK))/(n+1));
matProbK = nan(n-1,n+1);
matProbDelta = nan(n-1,2*n+2);
%vecProbK = ones(size(vecK))/(n+1);
for intI=vecI
dblX = vecX(intI);
vecY = binopdf(vecK,n,dblX);
matProbK(intI,:) = vecY;
vecDelta = zeros(1,2*n+2);
vecDelta(round(n-intI):round(2*n-intI)) = vecY;
matProbDelta(intI,:) = vecDelta;
%plot(vecY)
vecLogY = log(vecY);
vecLogY(isinf(vecLogY)) = -flintmax;
vecProbK = vecProbK + vecLogY;
vecProbK = vecProbK - log(sum(exp(vecProbK)));
%vecProbK = vecProbK + vecY;
end

subplot(2,3,2)
imagesc(vecI,vecK,matProbK')
axis xy

subplot(2,3,3)
vecDi = (-(n+1):(n+1))/(n+1);
imagesc(vecI,vecDi,matProbDelta')
axis xy