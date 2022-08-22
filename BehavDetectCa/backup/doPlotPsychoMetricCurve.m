function doPlotPsychoMetricCurve(sExp,sPC,intPlotType,boolSavePlots,strSes)
	%get data location
	if isfield(sExp,'Contrast')
		sBehavior = sExp;
	elseif isfield(sExp.structStim,'Contrast')
		sBehavior = sExp.structStim;
	end
	
	%get total number of repetitions
	cellFields = cell(1);
	cellFields{1} = 'Contrast';
	sTypes = getStimulusTypes(sBehavior,cellFields);
	cellSelect = getSelectionVectors(sBehavior,sTypes);
	intNumContrasts = length(cellSelect);
	
	%get plot type
	if nargin < 3
		intPlotType = 1;
	end
	%get if to save
	if nargin < 4
		boolSavePlots = false;
	end
	if nargin < 5
		strSes = '';
	end
	
	%plot
	if intPlotType==1,figure;else subplot(2,2,3);end
	set(gca,'XScale','log','YScale','linear') 
	vecX = sTypes.matTypes*100;
	vecLabelsX = vecX;
	if vecX(1) == 0
		vecX(1) = 0.25;
		vecLabelsX(1) = 0;
		semilogx(vecX(2:end),sPC.vecCorrect(2:end),'-x')
		hold on
		semilogx(vecX(1),sPC.vecCorrect(1),'-x')
		hold off
	else
		semilogx(vecX,sPC.vecCorrect)
	end
	set(gca,'XScale','log','YScale','linear') 
	set(gca,'XTick',vecX,'XTickLabel',vecLabelsX)
	ylabel('Fraction correct')
	ylim([0 1])
	xlabel('Contrast (%)')
	grid on
	if boolSavePlots && intPlotType==1
		drawnow;
		strFig = sprintf('%s_behavioraldetection_raw',strSes);
		export_fig([strFig '.tif']);
		export_fig([strFig '.pdf']);
	end
	%test for significant difference
	%[table,chi2,p] = crosstab(sPC.cellResponses{1},sPC.cellResponses{end});
	%title(sprintf('Chi-square test lowest vs highest contrast; p=%.3f',p))
	
	
	if intPlotType==1,figure;else subplot(2,2,4);end
	
	matRTs = sPC.matRTs*1000;
	vecMean = nanmean(matRTs,1);
	vecErr = nanstd(matRTs,[],1)./sqrt(sum(~isnan(matRTs),1));
	vecX = sTypes.matTypes*100;
	vecLabelsX = vecX;
	if vecX(1) == 0
		vecX = vecX(2:end);
		vecMean = vecMean(2:end);
		vecErr = vecErr(2:end);
		vecLabelsX = vecLabelsX(2:end);
	end
	errorbar(vecX,vecMean,vecErr,'-x')
	set(gca,'XScale','log','YScale','linear') 
	ylabel('Mean reaction time (ms)')
	set(gca,'XTick',vecX,'XTickLabel',vecLabelsX)
	xlabel('Contrast (%)')
	grid on
	
	if boolSavePlots
		drawnow;
		strFig = sprintf('%s_behavioralRT_raw',strSes);
		export_fig([strFig '.tif']);
		export_fig([strFig '.pdf']);
	end
end

