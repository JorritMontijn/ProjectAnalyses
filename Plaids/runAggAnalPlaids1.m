%% load preprocessed data
clear all;
strDataPath = 'F:\Data\Processed\PlaidData\';
strFileOut = [strDataPath 'PreProAggregatePlaids.mat'];
load(strFileOut);



%% plot population tuning curves grouped by stimulus, preferred ori vs z-resp
vecHorzC = [0 0 0 0 25 25 25 nan 50 50 50 nan 100 nan nan nan];
vecVertC = [0 25 50 100 0 25 50 nan 0 25 50 nan 0 nan nan nan];

intPrefOriGroups = 12;
vecPrefOriEdges = linspace(dblMaxOriPlot-180,dblMaxOriPlot,intPrefOriGroups+1);
vecPrefOriCenters = vecPrefOriEdges(2:end)-mean(diff(vecPrefOriEdges))/2;

cellAggPlotData = cell(numel(vecHorzC),intPrefOriGroups);

%%
vecPV = [2 3 4 7 8];%example: 7

for intRec=vecPV%1:numel(cellPlaids)
	%sRec
	%sRec(intRec).vecStopFrame = vecStopFrame;
	%sRec(intRec).vecStartFrame = vecStartFrame;
	%sRec(intRec).vecContrast0 = vecContrast0;
	%sRec(intRec).vecContrast90 = vecContrast90;
	%sRec(intRec).vecPrecedingContrast0 = vecPrecedingContrast0;
	%sRec(intRec).vecPrecedingContrast90 = vecPrecedingContrast90;
	%sRec(intRec).matResp = matResp(neuron,trial);
	%sRec(intRec).matRespPV = matRespPV(neuron,trial);
	%sRec(intRec).vecPrefDeg = vecPrefDeg(neuron);
	
	%retrieve data
	vecInclude = sRec(intRec).vecInclude;
	vecContrast0 = sRec(intRec).vecContrast0;
	vecContrast90 = sRec(intRec).vecContrast90;
	vecPrecedingContrast0 = sRec(intRec).vecPrecedingContrast0;
	vecPrecedingContrast90 = sRec(intRec).vecPrecedingContrast90;
	vecPrefDeg = sRec(intRec).vecPrefDeg;
	matResp = sRec(intRec).matRespPV;
	matPreResp = sRec(intRec).matPreRespPV;
	indGrat = vecContrast0 == 100 | vecContrast90 == 100;
	indPlaid = ~indGrat;
	
	%calc vars
	boolRemPreceding = true;
	indInclPrec = true(size(vecContrast0));
	if boolRemPreceding
		indInclPrec = vecPrecedingContrast0 == 0 & vecPrecedingContrast90 == 0;
	end
	vecNeuronOriGroup = sum(bsxfun(@gt,vecPrefDeg,vecPrefOriEdges),2);
	
	%loop through combos
	for intCombo = 1:16
		if isnan(vecHorzC(intCombo)),continue;end
		if intCombo == 1 %pre-stim baseline
			intOri0Contrast = vecHorzC(intCombo);
			intOri90Contrast = vecVertC(intCombo);
			
			vecUseStims = vecPrecedingContrast0 == 0 & vecPrecedingContrast90 == 0;
			for intOriGroup = 1:intPrefOriGroups
				vecAct = mean(matPreResp(vecNeuronOriGroup==intOriGroup & vecInclude,vecUseStims),2);
				cellAggPlotData{intCombo,intOriGroup} = cat(1,cellAggPlotData{intCombo,intOriGroup},vecAct(:));
			end
		else
			intOri0Contrast = vecHorzC(intCombo);
			intOri90Contrast = vecVertC(intCombo);
			
			vecUseStims = vecContrast0 == intOri0Contrast & vecContrast90 == intOri90Contrast & indInclPrec;
			for intOriGroup = 1:intPrefOriGroups
				vecAct = mean(matResp(vecNeuronOriGroup==intOriGroup & vecInclude,vecUseStims),2);
				cellAggPlotData{intCombo,intOriGroup} = cat(1,cellAggPlotData{intCombo,intOriGroup},vecAct(:));
			end
		end
	end
end

%% plot
dblMeanBaseline = mean(cell2vec(cellAggPlotData(1,:)));
h=figure;
vecPlot = 1:16;
for intPlot = vecPlot
	if isnan(vecHorzC(intPlot)),continue;end
	subplot(4,4,intPlot)
	
	intOri0Contrast = vecHorzC(intPlot);
	intOri90Contrast = vecVertC(intPlot);
	
	cellData = cellAggPlotData(intPlot,:);
	%transform data to vector
	vecPrefOriRads = deg2rad(vecPrefOriCenters);
	vecNumCells = cellfun(@numel,cellData);
	cellPrefRad = cellfun(@(x,y) x*y,vec2cell(vecPrefOriRads),cellfun(@ones,mat2cell([ones(size(vecNumCells)); vecNumCells]',ones(1,numel(vecNumCells)),2),'uniformoutput',false)','uniformoutput',false);
	vecListX = cell2vec(cellPrefRad);
	vecListY = cell2vec(cellData);
	vecMeanR = cellfun(@nanmean,cellData);
	
	%fit
	kappa1 = 1;
	gain1 = 0.01;
	kappa2 = 1;
	gain2 = 0.01;
	baseline = 0.02;
	vecParams0 = [kappa1 gain1 kappa2 gain2 baseline];
	sOptions = curvefitoptimoptions('curvefitfun','MaxFunEvals',1000,'MaxIter',1000,'Display','off');
	%all points
	vecFittedP = lsqcurvefit(@vonMisesFixedAngles, vecParams0, 2*vecListX, vecListY,[eps eps eps eps -1],[10000 10000 10000 10000 1],sOptions);
	vecFitY = vonMisesFixedAngles(vecFittedP,vecPrefOriRads*2);
	vecFitY0 = vonMisesFixedAngles(vecParams0,vecPrefOriRads*2);
	
	
	%plot
	plot(vecPrefOriCenters,vecFitY);
	hold on
	errorbar(vecPrefOriCenters,vecMeanR,cellfun(@nanstd,cellAggPlotData(intPlot,:))./sqrt(cellfun(@numel,cellAggPlotData(intPlot,:))),'kx')
	hold off
	if intOri0Contrast == 100 || intOri90Contrast == 100
		ylim([0.01 0.08]);
	else
		ylim([0.01 0.045]);
	end
	
	%labels
	title(sprintf('Horz C%d, Vert C%d',intOri0Contrast,intOri90Contrast));
	if ismember(intPlot,[4 10 11 13])
	xlabel('Preferred ori (degs)');
	end
	if ismember(intPlot,[1:4:16])
	ylabel('Avg. dF/F0');
	end
	set(gca,'xtick',-45:45:135)
	fixfig;
end
drawnow;
maxfig([]);