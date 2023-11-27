%% load spiking data & plot tuning curves


%% set recording
%close all;
clear all;
cellUniqueAreas = {...
	'PoissonPeakHeight',...Area 1
	'PoissonDoublePeak',...Area 2
	};

if isfolder('F:\Drive\MontijnHeimel_TimeseriesZeta')
	strPath = 'F:\Drive\MontijnHeimel_TimeseriesZeta';
	strDataSourcePath = 'F:\Data\Processed\Neuropixels\';
else
	strPath = 'C:\Drive\MontijnHeimel_TimeseriesZeta';
	strDataSourcePath = 'E:\DataPreProcessed\';
end
strDataTargetPath = fullfile(strPath,'\Data\');
strFigPath = fullfile(strPath,'\Figs\');
intMakePlots =0; %0=none, 1=normal plot, 2=including raster
vecRandTypes = [1 2];%[1 2];%1=normal,2=rand
vecRestrictRange = [0 inf];
boolSave = true;
vecResamps = 500;%250;%10:10:90;%[10:10:100];
intResampIdx = 1;
intResampNum = vecResamps(intResampIdx);
intResamps= numel(vecResamps);
intUseGenN = 1000;
boolUnbalanced = true;
intArea = 2;
intNeurons = intUseGenN;
dblFracDiffSpikes = 0.25;%max 0.5
optLow = 2;
optHigh = 1e2;
dblFixedBinWidth = 50/1000;

%% pre-allocate output variables
cellNeuron = cell(intNeurons,2,intResamps);
matTtest2 = nan(intNeurons,2,intResamps);
matZeta2 = nan(intNeurons,2,intResamps);
matAnova2 = nan(intNeurons,2,intResamps);
matAnova2_unbalanced = nan(intNeurons,2,intResamps);
matAnova2_optimal = nan(intNeurons,2,intResamps);

strArea = cellUniqueAreas{intArea};
strDate = getDate();

%% generate data
hTicN = tic;
for intNeuron=1:intNeurons
	%% message
	if toc(hTicN) > 5
		fprintf('Processing neuron %d/%d [%s]\n',intNeuron,intNeurons,getTime);
		hTicN=tic;
	end
	
	%% tuning params
	dblBaseRate = exprnd(1);
	boolDoublePeaked = false; %orientation or direction tuned
	dblPrefOri = rand(1)*2*pi; %preferred orientation (rads)
	dblKappa = rand(1)*5+5; %von Mises concentration parameter
	dblPrefRate = dblBaseRate; %mean single-spike rate during stimulus (exponential ISI)
	dblJitter = 1; %in ms'
	intSU = intNeuron;
	
	%% stimulus data
	dblStimDur = 3;
	dblPreBaseDur = 3;
	dblPostBaseDur = 2;
	dblTrialDur = dblPreBaseDur + dblStimDur + dblPostBaseDur;
	dblUseMaxDur = dblStimDur;
	intOris = 24;
	dblStepDeg = 360/intOris;
	vecOris = linspace(0,360-dblStepDeg,intOris);
	intReps = 10;
	intTrials = intOris*intReps;
	vecTrialAngles1 = nan(1,intTrials);
	for intRep=1:intReps
		vecTrialAngles1(((intRep-1)*intOris+1):(intRep*intOris)) = vecOris(randperm(intOris));
	end
	vecTrialStart1 = 10 + (dblPreBaseDur:dblTrialDur:(dblTrialDur*intTrials+eps));
	matTrialT1 = cat(2,vecTrialStart1',vecTrialStart1'+dblStimDur);
	
	%balanced or unbalanced?
	if boolUnbalanced
		strUb = 'B0';
		vecTrialStart2 =  10 + (dblPreBaseDur:dblTrialDur:(dblTrialDur*(intTrials/2)+eps));
		matTrialT2 = cat(2,vecTrialStart2',vecTrialStart2'+dblStimDur);
		vecTrialAngles2 = vecTrialAngles1(1:size(matTrialT2,1));
	else
		strUb = 'B1';
		matTrialT2 = matTrialT1;
		vecTrialAngles2 = vecTrialAngles1;
	end
	
	%% generate data
	% generate peak
	dblStartDelay = 0.053;
	dblPeakDelay1 = dblStartDelay;
	intAddSpikes1 = round(intTrials*dblFracDiffSpikes);
	[vecSpikeTimes1,dblPrefOri] = getGeneratedSpikingDataWithPeak(vecTrialAngles1,matTrialT1,dblBaseRate,dblPrefRate,dblJitter,dblKappa,boolDoublePeaked,dblPrefOri,intAddSpikes1,dblStartDelay,dblPeakDelay1);
	intTrials2 = size(matTrialT2,1);
				
	%% real+rand
	for intRandType=vecRandTypes
		%randomize
		if intRandType ==2
			%generate n2, no diff
			intAddSpikes2 = round((intAddSpikes1/intTrials)*intTrials2);
			dblPeakDelay2 = dblPeakDelay1;
			[vecSpikeTimes2,dblPrefOri] = getGeneratedSpikingDataWithPeak(vecTrialAngles2,matTrialT2,dblBaseRate,dblPrefRate,dblJitter,dblKappa,boolDoublePeaked,dblPrefOri,intAddSpikes2,dblStartDelay,dblPeakDelay2);
		else
			if strcmp(strArea,'PoissonDoublePeak')
				%generate n2, diff in peak time
				dblPeakDelay2 = dblPeakDelay1+0.002;
				intAddSpikes2 = round(dblFracDiffSpikes*intTrials2);
				[vecSpikeTimes2,dblPrefOri] = getGeneratedSpikingDataWithPeak(vecTrialAngles2,matTrialT2,dblBaseRate,dblPrefRate,dblJitter,dblKappa,boolDoublePeaked,dblPrefOri,intAddSpikes2,dblStartDelay,dblPeakDelay2);
			elseif strcmp(strArea,'PoissonPeakHeight')
				%generate n2, diff in peak height
				dblPeakDelay2 = dblPeakDelay1;
				intAddSpikes2 = 2*round(dblFracDiffSpikes*intTrials2);
				[vecSpikeTimes2,dblPrefOri] = getGeneratedSpikingDataWithPeak(vecTrialAngles2,matTrialT2,dblBaseRate,dblPrefRate,dblJitter,dblKappa,boolDoublePeaked,dblPrefOri,intAddSpikes2,dblStartDelay,dblPeakDelay2);
			else
				error('not recognized');
			end
		end
		%plot if first
		if intNeuron == 1
			%f = @() getZeta(vecSpikeTimes,vecUseTrialStart,dblTrialDur,intResampNum,0);
			%t=timeit(f);
			%f2 = @() getZeta2(vecSpikeTimes,vecUseTrialStart,dblTrialDur,intResampNum,0);
			%t2=timeit(f2);
			intPlot = 0;
		else
			intPlot = 0;
		end
		
		%% run tests
		[dblZeta2P,sZETA] = zetatest2(vecSpikeTimes1,matTrialT1,vecSpikeTimes2,matTrialT2,dblUseMaxDur,intResampNum,intPlot);
		
		%% prep
		%if balanced
		hTic2 = tic;
		[vecTrialPerSpike1,vecTimePerSpike1] = getSpikesInTrial(vecSpikeTimes1,matTrialT1(:,1),dblUseMaxDur);
		[vecTrialPerSpike2,vecTimePerSpike2] = getSpikesInTrial(vecSpikeTimes2,matTrialT2(:,1),dblUseMaxDur);
		if numel(vecTimePerSpike1) < 3 && numel(vecTimePerSpike2) < 3,continue;end
		
		%% optimal ANOVA
		xComb = sort(cat(1,vecTimePerSpike1,vecTimePerSpike2));
		[optN, dblC, allN, allC] = opthist(xComb);
		if optN<optLow,optN=optLow;end %at least 2 bins
		if optN>optHigh,optN=optHigh;end %at least 2 bins
		
		intTrials1 = size(matTrialT1,1);
		intTrials2 = size(matTrialT2,1);
		dblBinWidth = dblUseMaxDur/optN;
		vecBins = 0:dblBinWidth:dblUseMaxDur;
		matPSTH1 = nan(intTrials1,optN);
		matLabelN1 = ones(size(matPSTH1));
		matLabelBin1 = repmat(1:optN,[intTrials1 1]);
		for intTrial=1:intTrials1
			matPSTH1(intTrial,:) = histcounts(vecTimePerSpike1(vecTrialPerSpike1==intTrial),vecBins);
		end
		matPSTH2 = nan(intTrials2,numel(vecBins)-1);
		matLabelN2 = 2*ones(size(matPSTH2));
		matLabelBin2 = repmat(1:optN,[intTrials2 1]);
		for intTrial=1:intTrials2
			matPSTH2(intTrial,:) = histcounts(vecTimePerSpike2(vecTrialPerSpike2==intTrial),vecBins);
		end
		
		%if not balanced
		y = cat(1,matPSTH1(:),matPSTH2(:));
		g1 = cat(1,matLabelN1(:),matLabelN2(:));
		g2 = cat(1,matLabelBin1(:),matLabelBin2(:));
		[vecP,tbl,stats] = anovan(y,{g1 g2},'model','interaction','display','off');
		[h crit_p adj_p]=fdr_bh(vecP([1 3]));
		dblAnova2P_optimal = min(adj_p);
		
		%% fixed bin ANOVA
		vecBins = 0:dblFixedBinWidth:dblUseMaxDur;
		optN = numel(vecBins)-1;
		matPSTH1 = nan(intTrials1,optN);
		matLabelN1 = ones(size(matPSTH1));
		matLabelBin1 = repmat(1:optN,[intTrials1 1]);
		for intTrial=1:intTrials1
			matPSTH1(intTrial,:) = histcounts(vecTimePerSpike1(vecTrialPerSpike1==intTrial),vecBins);
		end
		matPSTH2 = nan(intTrials2,numel(vecBins)-1);
		matLabelN2 = 2*ones(size(matPSTH2));
		matLabelBin2 = repmat(1:optN,[intTrials2 1]);
		for intTrial=1:intTrials2
			matPSTH2(intTrial,:) = histcounts(vecTimePerSpike2(vecTrialPerSpike2==intTrial),vecBins);
		end
		%if balanced
		if intTrials1==intTrials2
			matPSTH = matPSTH1 - matPSTH2;
			dblAnova2P=anova1(matPSTH,[],'off');
		else
			dblAnova2P=1;
		end
		
		y = cat(1,matPSTH1(:),matPSTH2(:));
		g1 = cat(1,matLabelN1(:),matLabelN2(:));
		g2 = cat(1,matLabelBin1(:),matLabelBin2(:));
		[vecP,tbl,stats] = anovan(y,{g1 g2},'model','interaction','display','off');
		[h crit_p adj_p]=fdr_bh(vecP([1 3]));
		dblAnova2P_unbalanced= min(adj_p);
		
		%% t-test
		%'vecTtestP','vecTtestTime'
		hTic3 = tic;
		vecRespBinsDur1 = sort(flat([matTrialT1(:,1) matTrialT1(:,2)]));
		vecR1 = histcounts(vecSpikeTimes1,vecRespBinsDur1);
		vecD1 = diff(vecRespBinsDur1)';
		vecMu1 = vecR1(1:2:end)./vecD1(1:2:end);
		
		vecRespBinsDur2 = sort(flat([matTrialT2(:,1) matTrialT2(:,2)]));
		vecR2 = histcounts(vecSpikeTimes2,vecRespBinsDur2);
		vecD2 = diff(vecRespBinsDur2)';
		vecMu2 = vecR2(1:2:end)./vecD2(1:2:end);
		
		%get metrics
		[h,dblTtest2P]=ttest2(vecMu1,vecMu2);
		
		%% save
		% assign data
		cellNeuron{intNeuron,intRandType,intResampIdx} = [strArea strDate 'N' num2str(intSU)];
		matTtest2(intNeuron,intRandType) = dblTtest2P;
		matZeta2(intNeuron,intRandType) = dblZeta2P;
		matAnova2(intNeuron,intRandType) = dblAnova2P;
		matAnova2_unbalanced(intNeuron,intRandType) = dblAnova2P_unbalanced;
		matAnova2_optimal(intNeuron,intRandType) = dblAnova2P_optimal;
		%vecInclusionP(intIdx) = dblTtestP;
	end
end

%% save
if boolSave
	save([strDataTargetPath 'Zeta2DataAnova' strArea strUb 'Resamp' num2str(intResampNum) '.mat' ],...
		'cellNeuron','matTtest2','matZeta2','matAnova2','matAnova2_unbalanced','matAnova2_optimal');
end