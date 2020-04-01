clear all;
%% generate data
vecS = [10:10:100 200:100:1000 2000:1000:10000];
intNumS = numel(vecS);
intMaxIter = 10000;
matZ = nan(intNumS,intMaxIter);
for intIdxS= 1:intNumS
	intS = vecS(intIdxS);
	for intIter=1:intMaxIter
		vecRand = normrnd(0,1,[1 intS]);
		matZ(intIdxS,intIter) = max(abs(vecRand));
	end
end

%% test distros
%Elfving (1947), https://doi.org/10.1093/biomet/34.1-2.111
%Royston (1982), DOI: 10.2307/2347982 
%https://stats.stackexchange.com/questions/394960/variance-of-normal-order-statistics
%https://stats.stackexchange.com/questions/9001/approximate-order-statistics-for-normal-random-variables
%https://en.wikipedia.org/wiki/Extreme_value_theory
%https://en.wikipedia.org/wiki/Gumbel_distribution
figure
vecUseS = vecS;
dblAlpha=pi/8;
fMaxE = @(N) norminv((1-dblAlpha)./((N*2)-2*dblAlpha+1));
dblE = -fMaxE(vecUseS);
dblRealMean = mean(matZ,2)';

subplot(2,2,1)
plot([0.5 4.5],[0.5 4.5],'k--');
hold on
scatter(dblE,dblRealMean)
hold off;
xlabel('Real mean')
ylabel('Approx mean')
title('Approximation of mean');

fMaxV = @(N) ((N) ./ (((N+1).^2).*(N+2))).*(1./((normpdf(norminv(N./(N+1)))).^2));
dblV = fMaxV(vecUseS);

dblRealVar = var(matZ,[],2)';

subplot(2,2,2)
plot([0 0.35],[0 0.35],'k--');
hold on
scatter(dblRealVar,dblV)
hold off;
xlabel('Real var')
ylabel('Approx var')
title('Approximation of var');
%
%define beta
intNofS = 10;
intS = vecUseS(intNofS);
dblVar = dblV(intNofS);
dblBeta = (sqrt(6)*sqrt(dblVar))./(pi);

dblEulerMascheroni = vpa(eulergamma);
dblEulerMascheroni = 0.5772156649015328606065120900824;

dblMode = dblE(intNofS) - dblBeta.*dblEulerMascheroni;

dblApproxMean = dblMode + dblBeta*dblEulerMascheroni;
dblRealMean = mean(matZ(intNofS,:),2)';


fGumbel = @(x) (1./dblBeta).*exp(-(((x(:)-dblMode)./dblBeta) + exp(-(x(:)-dblMode)./dblBeta)));
fGumbelCDF = @(x) exp(-exp(-(((x(:)-dblMode)./(dblBeta)))));

dblX = 2.2;
dblGumbelCDF = fGumbelCDF(dblX);
dblP = (1-dblGumbelCDF);
dblZ = norminv(1-dblP/2);

dblApproxVar = ((pi^2)/6)*(dblBeta.^2);
dblRealVar = var(matZ(intNofS,:),[],2)';


subplot(2,2,3)

dblStepBin = 0.1;
vecBinEdges = 0:dblStepBin:5;
vecBins = vecBinEdges(2:end) - dblStepBin/2;
objHist= histogram(matZ(intNofS,:),vecBinEdges);
vecDistroGumbel = (size(matZ,2)*dblStepBin)*fGumbel(vecBins);
vecDistroGumbelCDF = (size(matZ,2)*dblStepBin)*fGumbelCDF(vecBins);

hold on
plot(vecBins,vecDistroGumbel,'r--');
plot(vecBins,vecDistroGumbelCDF,'g--');
plot(vecBins,cumsum(vecDistroGumbel)/10,'b--');
hold off
xlabel('X_n');
ylabel('Frequency (count)');
title(sprintf('n=%d; r=Gumbel pdf, g=cdf, b=csum(pdf)',intS));


subplot(2,2,4)
hold on
plot([0 1],[0 1],'k--')
scatter(cumsum(objHist.Values)/intMaxIter,cumsum(vecDistroGumbel)/intMaxIter,[],lines(1))
hold off
xlabel('Cumulative fraction (empirical)')
ylabel('Cumulative fraction (Gumbel)')

%% second test
%given
dblVar = 5.2989e-06;
dblE= 0.0091;

%calc
dblAlpha=pi/8;
fMaxE = @(N) -norminv((1-dblAlpha)./((N*2)-2*dblAlpha+1));

intN_star = fminsearch(@(n) abs(dblE - fMaxE(n)),100)

intN = fMaxE(0.49:0.001:0.51)


dblBeta = (sqrt(6)*sqrt(dblVar))./(pi);

dblEulerMascheroni = vpa(eulergamma);
dblEulerMascheroni = 0.5772156649015328606065120900824;

dblMode = dblE - dblBeta.*dblEulerMascheroni;

dblApproxMean = dblMode + dblBeta*dblEulerMascheroni;

fGumbel = @(x) (1./dblBeta).*exp(-(((x(:)-dblMode)./dblBeta) + exp(-(x(:)-dblMode)./dblBeta)));
fGumbelCDF = @(x) exp(-exp(-(((x(:)-dblMode)./(dblBeta)))));

dblX=0:0.001:0.02;
vecG = fGumbel(dblX);

hold on
plot(dblX,vecG*3.5)
