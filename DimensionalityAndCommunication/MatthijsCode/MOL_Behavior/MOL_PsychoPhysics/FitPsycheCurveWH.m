%% http://matlaboratory.blogspot.co.uk/2015/05/fitting-better-psychometric-curve.html
function [ffit, curve] = FitPsycheCurveWH(xAxis, yData, varargin)

% Start points and limits
if numel(varargin{:})>1
    useLims=1;
    UL=varargin{1}(1,:);
    SP=varargin{1}(2,:);
    LM=varargin{1}(3,:);
else
    useLims=0;
end

% Transpose if necessary
if size(xAxis,1)<size(xAxis,2)
    xAxis=xAxis';
end
if size(yData,1)<size(yData,2)
    yData=yData';
end

% Check range of data
if min(yData)<0 || max(yData)>1  
     % Attempt to normalise data to range 0 to 1
     yData = yData/(mean(yData)*2);
end
    
% Prepare fitting function
F=@(g,l,u,v,x) g+(1-g-l)*0.5*(1+erf((x-u)/sqrt(2*v^2)));

% Fit using fit function from fit toolbox
if useLims==1
    % SPs and limits specified, use while fitting
    ffit=fit(xAxis,yData,F,'StartPoint',SP,'Upper',UL,'Lower',LM);
else
    % Fits not specified, don't use while fitting
    ffit=fit(xAxis,yData,F);
end

% Create a new xAxis with higher resolution
fineX = linspace(min(xAxis),max(xAxis),numel(xAxis)*50);
% Generate curve from fit
curve = feval(ffit, fineX);
curve = [fineX', curve];