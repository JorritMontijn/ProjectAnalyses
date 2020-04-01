%% analyze input strength dependency for dimensionality of pop responses

%% settings
clearvars;
boolSaveFigs = true;
vecRunSims = [21 24 31 34];

%run
for intLoadSim=vecRunSims
	clearvars -except vecRunSims intLoadSim bool*
	boolLoad = true;
	if intLoadSim == 11 && boolLoad
		strSimulation = 'xAreaDistributed_OriFull_2017-06-15';
	elseif intLoadSim == 12 && boolLoad
		strSimulation = 'xAreaDistributed_ContExcOnly_2017-06-08'; %contrast
	elseif intLoadSim == 13 && boolLoad
		strSimulation = 'xAreaDistributed_ContInhOnly_2017-06-12'; %contrast
	elseif intLoadSim == 14 && boolLoad
		strSimulation = 'xAreaDistributed_ContNone_2017-06-12'; %contrast
		
	elseif intLoadSim == 21 && boolLoad
		strSimulation = 'xAreaDistributed_LumFull_2017-06-26'; %luminance full
	elseif intLoadSim == 22 && boolLoad
		strSimulation = ''; %lum exc
	elseif intLoadSim == 23 && boolLoad
		strSimulation = ''; %lum inh
	elseif intLoadSim == 24 && boolLoad
		strSimulation = 'xAreaDistributed_LumNone_2017-06-26'; %lum none
		
		
	elseif intLoadSim == 31 && boolLoad
		strSimulation = 'xAreaDistributed_ContFull_2017-06-26'; %contrast full
	elseif intLoadSim == 32 && boolLoad
		strSimulation = ''; %contrast exc
	elseif intLoadSim == 33 && boolLoad
		strSimulation = ''; %contrast inh
	elseif intLoadSim == 34 && boolLoad
		strSimulation = 'xAreaDistributed_ContNone_2017-06-26'; %contrast none
		
	end
	
	%% RUN: #header
	strFigDir = 'D:\Data\Results\Block1\';
	cellIn = strsplit(strSimulation,'_');
	strFiller = cellIn{1};
	strType = cellIn{2};
	strDate = cellIn{3};
	strTag = [strType strDate];
	cellStrArea = {'V1','V2','V1-V2'};
	for intWithinArea=[1 2 3]
		%% load
		% load data
		strDataFile = ['D:\Data\Processed\V1_LIFmodel\Simulation_' strSimulation '_AB1_Area' num2str(intWithinArea) '.mat'];
		sSave = load(strDataFile);
		strArea = cellStrArea{intWithinArea};
		intParamNum = numel(sSave.cellI);
		mapC = redbluepurple(intParamNum);
		vecC = linspace(0,100,intParamNum);
		
		%% plot Fisher I
		if isfield(sSave,'matI') && ~all(isnan(sSave.matI(:)))
			vecGroupSizes = 2.^(0:(size(sSave.matI,2)-1));
			
			hFigI = figure;
			hAxI1 = subplot(2,2,1);ylabel('Bias-corrected Fisher information (d''^2)');xlabel('Population size (# of neurons)');
			title(sprintf('Unshuffled; %s; Area %s',strType,strArea));hold on;
			hAxI2 = subplot(2,2,2);ylabel('Bias-corrected Fisher information (d''^2)');xlabel('Population size (# of neurons)');
			title(sprintf('Shuffled; %s; Area %s',strType,strArea));hold on;
			%full screen
			drawnow;
			jFig = get(handle(gcf), 'JavaFrame');
			jFig.setMaximized(true);
			figure(gcf);
			drawnow;
			
			% plot information
			for intStimTypeIdx=1:intParamNum
				cellC{intStimTypeIdx} = sprintf('%03d',vecC(intStimTypeIdx));
				plot(hAxI1,vecGroupSizes,sSave.matI(intStimTypeIdx,:),'Color',mapC(intStimTypeIdx,:));
				plot(hAxI2,vecGroupSizes,sSave.matI_shuff(intStimTypeIdx,:),'Color',mapC(intStimTypeIdx,:));
			end
			drawnow;
			
			dblMaxY = max([max(get(hAxI1,'ylim')) max(get(hAxI2,'ylim'))]);
			ylim(hAxI1,[-1 dblMaxY]);fixfig(hAxI1);legend(cellC,'Location','northwest');
			ylim(hAxI2,[-1 dblMaxY]);fixfig(hAxI2);legend(cellC,'Location','northwest');
			
			%save figure
			if boolSaveFigs
				figure(hFigI);drawnow;
				export_fig([strFigDir  'AnalBlock1Info_Area' num2str(intWithinArea) '_' strTag '.tif']);
				export_fig([strFigDir  'AnalBlock1Info_Area' num2str(intWithinArea) '_' strTag '.pdf']);
			end
		end
		
		%% plot dim dep
		if isfield(sSave,'matPredDimDepR2')
			%prep
			hFigDD = figure;
			hAxD1 = subplot(2,3,1);
			xlabel('Number of predictive dimensions');
			ylabel('Predictive performance');
			title(sprintf('Semedo Fig 4; %s; Area %s',strType,strArea));
			hold on;
			
			hAxD2 = subplot(2,3,2);
			plot([0 10],[0 10],'k--');
			hold on;
			xlabel('Target population dimensionality')
			ylabel('Number of predictive dimensions')
			
			hAxD3 = subplot(2,3,3);
			xlabel('Number of dominant dimensions');
			ylabel('Predictive performance');
			hold on;
			%full screen
			drawnow;
			jFig = get(handle(gcf), 'JavaFrame');
			jFig.setMaximized(true);
			figure(gcf);
			drawnow;
			
			%% do plotting
			for intStimTypeIdx=1:intParamNum
				cellC{intStimTypeIdx} = sprintf('%03d',vecC(intStimTypeIdx));
				matPredictionsR2_SingleNeurons = sSave.cellPredictionsR2_SingleNeurons{intStimTypeIdx};
				matPredictionsR2 = sSave.cellPredictionsR2{intStimTypeIdx};
				matPredDimDepR2 = sSave.cellPredDimDepR2{intStimTypeIdx};
				matPredictiveDimensions = sSave.cellPredictiveDimensions{intStimTypeIdx};
				matTargetPopDimensionality = sSave.cellTargetPopDimensionality{intStimTypeIdx};
				matR2RemPredDim = sSave.cellR2RemPredDim{intStimTypeIdx};
				matR2UseDomDim = sSave.cellR2UseDomDim{intStimTypeIdx};
				
				%fig 4
				axes(hAxD1);
				vecPredDD = squeeze(nanmean(nanmean(matPredDimDepR2,2),1));
				plot(1:length(vecPredDD),vecPredDD,'Color',mapC(intStimTypeIdx,:));
				scatter(length(vecPredDD)-0.5,nanmean(matPredictionsR2(:)),80,mapC(intStimTypeIdx,:),'x');
				%if intContrastIdx==1,text(1,mean(matPredictionsR2(:))+0.02,'Full prediction','Color',mapC(intContrastIdx,:),'FontSize',14,'Rotation',45);end
				
				%fig 5
				axes(hAxD2);
				dblDomDim = mean(matTargetPopDimensionality(:));
				dblPredDim = mean(matPredictiveDimensions(:));
				scatter(dblDomDim,dblPredDim,80,mapC(intStimTypeIdx,:),'x');
				
				%fig 7
				axes(hAxD3);
				intMaxDimAnal = size(matR2UseDomDim,3);
				vecPlot = (1:intMaxDimAnal);
				vecDomDD =  squeeze(nanmean(nanmean(matR2UseDomDim,2),1));
				%plot(vecPlot,vecPredDD(vecPlot)','Color',mapC(intStimTypeIdx,:),'Marker','x');
				plot(vecPlot,vecDomDD(vecPlot)','Color',mapC(intStimTypeIdx,:));%,'Marker','o');
				drawnow;
			end
			
			%fix figs
			ylim(hAxD1,[0 max(get(hAxD1,'ylim'))]);fixfig(hAxD1);
			axes(hAxD2);fixfig;
			ylim(hAxD3,[0 max(get(hAxD3,'ylim'))]);fixfig(hAxD3);legend(cellC,'Location','southeast');
			
			%save figure
			if boolSaveFigs
				figure(hFigDD);drawnow;
				export_fig([strFigDir  'AnalBlock1Pred_Area' num2str(intWithinArea) '_' strTag '.tif']);
				export_fig([strFigDir  'AnalBlock1Pred_Area' num2str(intWithinArea) '_' strTag '.pdf']);
			end
		end
		
		
		%% pre-allocate variables mean
		if isfield(sSave,'vecPopMax')
			% get data 
				sSave.vecPopMax = vecPopMax;
				sSave.vecPopStDev =  vecPopStDev;
				sSave.vecPopMean = vecPopMean;
				sSave.vecPopHet = vecPopHet;
				
				sSave.matNoiseCorrsS1 = matNoiseCorrsS1;
				sSave.matNoiseCorrsS2 = matNoiseCorrsS2;
				
			%prepare figures
			hFigPM1 = figure;
			
			%full screen
			drawnow;
			jFig = get(handle(gcf), 'JavaFrame');
			jFig.setMaximized(true);
			figure(gcf);
			drawnow;
			
			%fig 2
			hFigPM = figure;
			hAxM1 = subplot(2,3,1);
			hold on;
			hAxM2 = subplot(2,3,2);
			hold on;
			hAxM3 = subplot(2,3,3);
			hold on;
			
			%full screen
			drawnow;
			jFig = get(handle(gcf), 'JavaFrame');
			jFig.setMaximized(true);
			figure(gcf);
			drawnow;
				
			%set props
			axes(hAxM1);
			ylabel('Normalized density');
			xlabel('Noise correlation (r)');fixfig;
			
			axes(hAxM2);
			plot(hAxM2,vecParamVals,vecNC_StDev,'k');
			hold on
			plot(hAxM2,vecParamVals,vecNC_Mean,'r');
			xlim([min(vecParamVals) max(vecParamVals)]);
			hold off
			legend('SD','Mean')
			xlabel('Stimulation Intensity')
			ylabel('Noise correlation (r)')
			fixfig;
			
			axes(hAxM3);
			plot(vecParamVals,vecPopMax);
			xlabel('Stimulation Intensity')
			ylabel('Maximum firing rate (Hz)')
			fixfig;
			
			subplot(2,3,4)
			plot(vecParamVals,vecPopStDev,'k');
			hold on
			plot(vecParamVals,vecPopMean,'r');
			hold off
			legend('SD','Mean')
			xlabel('Stimulation Intensity')
			ylabel('Population firing rate (Hz)')
			ylim([0 max(get(gca,'ylim'))]);
			fixfig;
			
			subplot(2,3,5)
			vecCV = vecPopStDev./vecPopMean;
			vecCV(abs(vecCV)>10)=nan;
			plot(vecParamVals,vecCV);
			xlabel('Stimulation Intensity')
			ylabel('Population CV')
			fixfig;
			
			subplot(2,3,6)
			plot(vecParamVals,vecPopHet);
			xlabel('Stimulation Intensity')
			ylabel('Population heterogeneity')
			fixfig;
			
			%drawnow;
			
			%save figure
			if boolSaveFigs
				figure(hFigPM);drawnow;
				%export_fig([strFigDir  'AnalBlock1PopM_Area' num2str(intWithinArea) '_' strTag '.tif']);
				%export_fig([strFigDir  'AnalBlock1PopM_Area' num2str(intWithinArea) '_' strTag '.pdf']);
			end
		end
		
		
	end
end
