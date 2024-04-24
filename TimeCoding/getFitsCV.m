function [sFits,vecCV,vecMean,vecSd] = getFitsCV(vecAllSpikeTime,dblTotDur,vecTimescales,boolFixedAsymptote)
	%getFitsCV Fits Lin, Exp decay and Root functions to CV/timescale curves
	%   [sFits,vecCV,vecMean,vecSd] = getFitsCV(vecAllSpikeTime,dblTotDur,vecTimescales,boolFixedAsymptote)
	
	
	%optional inputs
	if ~exist('boolFixedAsymptote','var') || isempty(boolFixedAsymptote)
		boolFixedAsymptote = true;
	end
	
	%% get lin fit for mean/sd
	intTimescaleNum = numel(vecTimescales);
	vecSlopes = nan(size(vecTimescales));
	vecR2 = nan(size(vecTimescales));
	vecMean = nan(size(vecTimescales));
	vecSd = nan(size(vecTimescales));
	intK=1;
	for intScale=1:intTimescaleNum
		vecBins = 0:vecTimescales(intScale):dblTotDur;
		vecCounts = histcounts( vecAllSpikeTime,vecBins);
		vecMean(intScale) = mean(vecCounts);
		vecSd(intScale) = std(vecCounts);
	end
	%fit linear model
	mdl = fitlm(vecMean,vecSd);
	dblSlope_SdMean = mdl.Coefficients.Estimate(2);
	vecFitY_SdMean =  mdl.Fitted;
	[dblR2_SdMean,dblSS_tot,dblSS_res,dblT,dblP,dblR2_adjusted,dblR2_SE] = getR2(vecSd,vecFitY_SdMean,intK);
	
	%% fit cv/timescale pop
	fExp = fittype('a+b*exp(-x/c)',...
		'dependent',{'y'},'independent',{'x'},...
		'coefficients',{'a','b','c'});
	
	vecX = vecTimescales';
	vecCV = (vecSd./vecMean)';
	vecStartCoeffs = [vecCV(end) vecCV(1)/vecTimescales(1) 1/2];
	vecUpper = [1e16 1e16 1];
	vecLower = [0 0 0];
	if boolFixedAsymptote
		fRoot = fittype('(1/((b*x)^c))',...
			'dependent',{'y'},'independent',{'x'},...
			'coefficients',{'b','c'});
		vecStartCoeffs = vecStartCoeffs(2:end);
		vecUpper = vecUpper(2:end);
		vecLower = vecLower(2:end);
	else
		fRoot = fittype('a+(1/((b*x)^c))',...
			'dependent',{'y'},'independent',{'x'},...
			'coefficients',{'a','b','c'});
	end
	
	intK=1;
	%fit linear model
	mdl = fitlm(vecX,vecCV);
	
	dblSlope_Lin = mdl.Coefficients.Estimate(2);
	vecFitY =  mdl.Fitted;
	[dblR2_Lin,dblSS_tot,dblSS_res,dblT,dblP,dblR2_adjusted,dblR2_SE] = getR2(vecCV,vecFitY,intK);
	
	%fit exp model
	vecStartCoeffsExp = [vecCV(end) vecCV(1)-vecCV(end) 0.1];
	[fitobject,gof] = fit(vecX,vecCV,fExp,'lower',[0 0 1e-6],'upper',[1e16 1e16 1e16],'startpoint',vecStartCoeffsExp);
	vecFitExp = fitobject(vecX);
	[dblR2_Exp,dblSS_tot,dblSS_res,dblT,dblP,dblR2_adjusted,dblR2_SE] = getR2(vecCV,vecFitExp,intK);
	dblHalfLife_Exp = fitobject.c*log(2);
	dblAsymptote_Exp = fitobject.a;
	dblScale_Exp = fitobject.b;
	
	%fit root model
	[fitobject,gof,output] = fit(vecX,vecCV,fRoot,'lower',vecLower,'upper',vecUpper,'startpoint',vecStartCoeffs,...
		'tolfun',1e-16,'tolx',1e-16,'maxfunevals',1e3,'maxiter',1e3);
	vecFitRoot = fitobject(vecX);
	[dblR2_Root,dblSS_tot,dblSS_res,dblT,dblP,dblR2_adjusted,dblR2_SE] = getR2(vecCV,vecFitRoot,intK);
	if boolFixedAsymptote
		dblAsymptote_Root=0;
	else
		dblAsymptote_Root = fitobject.a;
	end
	dblScale_Root = fitobject.b;
	dblExponent_Root = fitobject.c;
	
	%% save data
	sFits = struct;
	sFits.LinMuSd.R2=dblR2_SdMean;
	sFits.LinMuSd.Slope=dblSlope_SdMean;
	sFits.LinMuSd.FitY=vecFitY_SdMean;
	sFits.Lin.R2=dblR2_Lin;
	sFits.Lin.Slope=dblSlope_Lin;
	sFits.Lin.FitY = vecFitY;
	sFits.Exp.R2 = dblR2_Exp;
	sFits.Exp.HalfLife = dblHalfLife_Exp;
	sFits.Exp.Asymptote = dblAsymptote_Exp;
	sFits.Exp.Scale = dblScale_Exp;
	sFits.Exp.FitY = vecFitExp;
	sFits.Root.R2 = dblR2_Root;
	sFits.Root.Exponent = dblExponent_Root;
	sFits.Root.Asymptote = dblAsymptote_Root;
	sFits.Root.Scale = dblScale_Root;
	sFits.Root.FitY = vecFitRoot;
	
end

