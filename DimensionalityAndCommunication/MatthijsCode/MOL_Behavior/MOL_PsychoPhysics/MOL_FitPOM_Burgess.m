function [c50, exponent, curve] = MOL_FitPOM_Burgess(xdata,ydata)


%% Example data:

% xdata = [0 0.01 0.02 0.05 0.1 0.25 0.5 1];
% ydata = [0 0 0.05 0.08 0.3 0.7 0.75 0.8]; 
% 
% xdata = 0:0.1:1;
% ydata = [0 0 0 0.1 0.3 0.5 0.7 0.9 0.95 1 1];
% 
% [c50, exponent, curve] = MOL_FitPOM_Burgess(xdata,ydata);

% zl = bl + sl * f(cl)


%%
% close all
% n = 5;
% chalf = 0.15;
% y = contrast.^n./(chalf^n + contrast.^n)
% figure; plot(contrast,y,'b','LineWidth',2); hold on; plot(contrast,resp,'k.','MarkerSize',20);

%% Fit


% F=@(h,n,x) x.^n./(h.^n + x.^n);
% 
% % Set limits and starting positions for fit:
% UL =       [0.9,    3];     % Upper limits for g, l, u ,v
% SP =       [0.5,    2];      % Start points for g, l, u ,v
% LL =       [0,      0];        % Lower limits for  g, l, u ,v
% 
% ffit=fit(xdata',ydata',F,'StartPoint',SP,'Upper',UL,'Lower',LL);
% 
% c50         = ffit.h;
% exponent    = ffit.n;

%% more complicated fit:

F=@(g,s,u,v,x) g + s * (x.^v./(u.^v + x.^v));

UL =       [0.5,    1,      2,      50];     % Upper limits for g, l, u ,v
SP =       [0.1,    0.5,    0.5,    2];      % Start points for g, l, u ,v
LL =       [0,      0,      0,      0];        % Lower limits for  g, l, u ,v

ffit=fit(xdata',ydata',F,'StartPoint',SP,'Upper',UL,'Lower',LL);

baseline    = ffit.g;
sensitivity = ffit.s;
c50         = ffit.u;
exponent    = ffit.v;

%%

% Create a new xAxis with higher resolution
fineX = linspace(min(xdata),max(xdata),numel(xdata)*50);
% Generate curve from fit
curve = feval(ffit, fineX);
curve = [fineX', curve];

figure;
plot(curve(:,1),curve(:,2),'b','LineWidth',2); 
hold on; 
plot(xdata,ydata,'k.','MarkerSize',20);


end
