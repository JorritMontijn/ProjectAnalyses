
function [mu, stddev, guessrate, lapserate, curve] = MOL_FitCumGauss(xdata,ydata)

% Function: F=@(g,l,u,v,x) g+(1-g-l)*0.5*(1+erf((x-u)/sqrt(2*v^2)));

% The psychometric function presented in Wichmann and Hill's paper
% is a cumulative gaussian function with 4 parameters:
% Mean (u): The mean value of the distribution representing subject bias.
% Standard deviation (v): The variation of the distribution representing
%   the subjects discrimination sensitivity.
% Guess rate (g) and lapse rate (l): Two additional parameters representing
%   the subjects fallibility (ie. potential inability to ever reach 100%
%   performance) at each end of the distribution/stimulus spectrum.

SPmu    = xdata(ceil((numel(xdata)+1)/2));

SPv     = mean(diff(xdata));
ULv     = SPv*10;
LLv     = SPv/10;

% SPv     = abs(mean(diff(ydata(2:end-1)))) * mean(diff(xdata));
% ULv     = SPv*100;
% LLv     = SPv;

% Set limits and starting positions for fit:
UL =       [0.5,    0.5,    inf,    ULv];     % Upper limits for g, l, u ,v
SP =       [0.2,    0.1,    SPmu,    SPv];      % Start points for g, l, u ,v
LL =       [0,      0,      -inf,   -LLv];        % Lower limits for  g, l, u ,v

%Check range of data:
if min(ydata)<0 || max(ydata)>1
    % Attempt to normalise data to range 0 to 1
    ydata = ydata/(mean(ydata)*2);
end

% Prepare fitting function
F=@(g,l,u,v,x) g+(1-g-l)*0.5*(1+erf((x-u)/sqrt(2*v^2)));

% SPs and limits specified, use while fitting
ffit=fit(xdata',ydata',F,'StartPoint',SP,'Upper',UL,'Lower',LL);

mu          = ffit.u;
stddev      = ffit.v;
guessrate   = ffit.g;
lapserate   = ffit.l;

% Create a new xAxis with higher resolution
fineX = linspace(min(xdata),max(xdata),numel(xdata)*50);
% Generate curve from fit
curve = feval(ffit, fineX);
curve = [fineX', curve];

%     case 'weibull'
%
%         disp('blabla')
%
%
%         figure;
% x = -6:0.1:6;
% alpha = 0;
% beta = 1;
% plot(x,1./(1+exp(-(x-alpha)./beta)), 'k-')
% title('The sigmoid function');
% xlabel('Stimulus value (e.g. intensity)');
% ylabel('Probability of stimulus detection');
% % Varying the alpha variable shifts the Sigmoid left and right. Varying the
% % beta variable scales the steepness of the slope.
%
% figure;
% pGuess.t =0.1;
% pGuess.b = 2;
%
% x = exp(linspace(log(min(results.intensity)),log(max(results.intensity)),101));
%
% %x = logspace(log(min(results.intensity)),log(max(results.intensity)),101);
% y= Weibull(pGuess,x);
%
%
%
% hold on
%
% plot(log(x),y*100,'r-','LineWidth',2);
%
%
% % Starting parameters
% pInit.t = .1;
% pInit.b = 2;
%
% [pBest,logLikelihoodBest] = fit('fitPsychometricFunction',pInit,{'b','t'},results,'Weibull')

end

