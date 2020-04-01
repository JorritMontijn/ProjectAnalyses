%% description
%{
Dat is goed nieuws! Hieronder mijn antwoorden op je vragen:

    -In the Kohn paper they take two different neuronal populations. So
    they look at simultaneous rate fluctuations. Our situation would be
    quite different because we look at rate fluctuations in the same
    neurons, but separated in time. What would be the consequences?

In principe maakt dit het alleen maar makkelijker. De eerste vraag die ik
zou willen beantwoorden, is of de populatieactiviteit in de 1e en 2e bump
in dezelfde space liggen of niet. Ik zou eerst met factor analysis (FA) of
PCA een ordering maken van de components die de meeste variance verklaren
in de 1e bump (B1), en een andere ordering voor de 2e bump (B2). Vervolgens
kun je een normale regression op B2 doen van de B1 projectie in de reduced
space voor de cumulatieve subspaces waarbij je steeds de volgende meest
belangrijke PCA/FA dimensie toevoegt. Dan krijg je dus een R^2 voor het
predicten van B2 door B1 door een funnel die je zelf hebt bepaald door de
belangrijkste principal components (PCs) in B1 space te kiezen. Vervolgens
je deze curve vergelijken met wat je krijgt als je een directe reduced-rank
regression doet van B1 naar B2. Deze laatste is dus eigenlijk vergelijkbaar
met een greedy decoder, en zal altijd sneller groeien dan je PC-reduced
regression. Als de curves perfect overlappen zijn de belangrijkste
predictive dimensions van B1->B2 ook degene die de meeste interne
variabiliteit in B1 laten zien.

Het mooie van dezelfde latent space voor B1 en B2 gebruiken is dat je nu
een symmetrische vergelijking kan doen, want je kan nu ook de B2-PCs
gebruiken om B1 in te projecteren, en de omgekeerde analyse doen. Als je
dezelfde R^2 curve krijgt uit B1(B1-PCs)->B2 als uit B1(B2-PCs)->B2, dan
weet je zeker dat B1 en B2 in identieke subspaces liggen, met als enige
caveat dat 1 van de 2 groter kan zijn, maar dat in dat geval iig alle
hogere dimensies een kleinere PC moeten hebben dan de dimensies die beide
subspaces delen. Als de curves niet identiek zijn, dan weet je dat de
spaces niet hetzelfde zijn, en kun je afleiden uit de groeisnelheid welke
space groter is.

Vervolgens kun je dit nog vergelijken met de B1(B1-PCs)->B1 met
B2(B2-PCs)->B2 predictions. Als deze interne predictions identiek zijn aan
de PC-matched reduction (dus B1(B1-PCs)->B2 en B2(B2-PCs)->B1), dan weet je
dat de spaces hetzelfde zijn. En als de interne predictions hoger zijn, dan
delen B1 en B2 een subspace, maar zijn ze niet identiek.

 Ter samenvatting, uiteindelijk hebben we de volgende curves, waarbij (Bx) de PC-space van x aangeeft.
1) B1(B1)->B1 (interne prediction B1: PCA/FA)
2) B1(B1)->B2 (hoeveel overlapt B2 met interne variance in B1?)
3) B1(RRR)->B2 (greedy regression: wat is de grootte van de B1->B2 subspace?
	[vergelijk met (2): zijn predictive en internal dimensions hetzelfde?])
4) B2(B2)->B2 (interne prediction B2: PCA/FA)
5) B1(B2)->B2 (inverse regression: zijn interne dimensies van B2 hetzelfde als B1
	[vergelijk met (2) en (4)]?)

Als bonus kun je B1 en B2 nog omdraaien in (2,3,5) om te controleren dat
het hetzelfde antwoord geeft. En je kan ook controleren of PCA en FA hetzelfde geven.

Wat heb je hier nu aan?
Het zou dus uiteindelijk mooi zijn als we bijvoorbeeld vinden dat de B1
voor 80% overlapt met B2, en dat B2 voor 50% overlapt met B1 (B2 is dan dus
groter). Dan hebben we drie subspaces: B1-unique, B1-B2-shared, en
B2-unique. We kunnen nu de populatieactiviteit op deze spaces projecteren
en kijken welke tot de beste stimulus decoding leidt, en welke het beste
correleert met de behavioural decision. Het geeft ook een aanknopingspunt
om het daarna met de andere gebieden te vergelijken. Ik zou het persoonlijk
wel interessant vinden om bijvoorbeeld te kijken of de B1-B2 overlap
verschilt tussen gebieden, en of de inter-areal communication specifiek via
1 van deze subspaces gaat. Uiteindelijk kun je deze inter-areal analyse ook
omdraaien en RRR tussen gebieden doen om te kijken in welk epoch je de
meeste activiteit vindt (zie onder).


    -In the Kohn paper they study ONLY the trial-to-trial fluctuations.
    This means that they only study the variability as proxy of
    communication between subpopulations. In our case we need to consider
    for a second how to go to the noise correlations given multiple
    orientations, multiple frequencies and multiple orientation/frequency
    change amplitudes.

Het feit dat ze alleen naar noise hebben gekeken klinkt inderdaad wat
vreemd, maar ze laten in de supp mat wel zien dat hun punt algemeen is voor
elke orientatie en er een smooth progression zit in de communication
subspace: hij draait dus als het ware mee met de stimulus manifold. Feature
conjunctions is wat dat betreft meer een open vraag, maar je had tijdens de
experimenten voornamelijk gefocust op near-threshold stimuli toch? Dus ik
weet niet hoeveel we daar nu mee kunnen doen. Wat we iig wel kunnen doen is
dezelfde stimulus nemen en de across-area RRR tijdens de 1st en 2nd bump
bekijken om te zien of de communication subspace onveranderlijk is of niet.
Als hij namelijk stationair is, is het aannemelijk dat de communication
subspace een direct anatomical substrate heeft. Als hij echter
non-stationair is over tijd (1st vs 2nd bump), maar wel vergelijkbaar over
trials, dan worden de signalen blijkbaar dynamisch doorgestuurd, maar ik
heb op het moment geen idee hoe je dat zou kunnen bereiken..

%}
%% clear
clear all;

%% set variables
%On which timestamp to align as t=0
strAlignOn = 'Change';
intFullIters = 10;
intResamplings = 10;
boolShuffle = false;
dblLambda = 500;
boolSavePlots = true;
strFigDir = 'D:\Data\Results\BumpsMatthijs\';

% parameters for window of interest:
dblStartBaseT = -0.4;
dblStopBaseT = -0.2;
dblStartEp1T = 0;
dblStopEp1T = 0.2;
dblStartEp2T = 0.2;
dblStopEp2T = 0.4;

%% load data
fprintf('Loading data... [%s]\n',getTime);
[Data] = MOL_GetData('D:\Data\Processed\MatthijsOudeLohuis\','CHDET',{'ChangeDetectionConflict'},{'2003' '2009' '2010' '2011' '2012' '2013'},[],{'sessionData' 'trialData' 'spikeData'});
sessionData     = Data.sessionData;
trialData       = Data.trialData;
spikeData       = Data.spikeData;

%% Remove last 20 trials:
trialData = MOL_RemoveLastnTrials(trialData,20);

%% Filter good neurons:
spikeData = MOL_filterNeurons(sessionData,trialData,spikeData);

%% split data by session
[cellSessions,ia,vecSessionIdx] = unique(spikeData.session_ID);

%% split by area
[cellAreas,ia,vecAreaIdx] = unique(spikeData.area);

%% Get significant first and second bump responses:
%[sessionData,trialData,spikeData] = MOL_calc_signresp(sessionData,trialData,spikeData);

%% make populations
[cellPopID,ia,vecPopulationIdx] = unique(vecSessionIdx*10 + vecAreaIdx);

%% define V1 populations
indV1=vecAreaIdx==find(strcmp(cellAreas,'V1'));
vecPopsV1 = unique(vecPopulationIdx(indV1));
intPops = numel(vecPopsV1);

%% pre-allocate output
vecCellsPerPop = nan(1,intPops);
cellM_xx2x = cell(intFullIters,intPops);
cellM_xx2y = cell(intFullIters,intPops);
cellM_xRRR2y = cell(intFullIters,intPops);
cellM_xy2y = cell(intFullIters,intPops);

cellM_yy2y = cell(intFullIters,intPops);
cellM_yy2x = cell(intFullIters,intPops);
cellM_yRRR2x = cell(intFullIters,intPops);
cellM_yx2x = cell(intFullIters,intPops);

%%  analyze
dblTimeConversion = 10^6;
for intPopIdx=1:intPops
	%% prepare, retrieve and transform data
	%retrieve cells
	intPopulation = vecPopsV1(intPopIdx);
	indUseCells = vecPopulationIdx==intPopulation;
	intCells = sum(indUseCells);
	intCutOffCellNr = 16;
	vecCellsPerPop(intPopIdx) = intCells;
	if intCells < intCutOffCellNr
		fprintf('Pop %d, Number of cells is %d, which is under %d, skipping... [%s]\n',intPopIdx,intCells,intCutOffCellNr,getTime);
		continue;
	end
	
	%get spikes and transform to seconds
	cellSpikes = cellfun(@rdivide,spikeData.ts(indUseCells),cellfill(dblTimeConversion,[intCells 1]),'uniformoutput',false);
	strSesID = unique(spikeData.session_ID(indUseCells));
	if numel(strSesID)~=1,error([mfilename ':SessionNotUnique'],'Something went wrong with population selection...');end
	
	%get trial data
	indSesTrials = strcmp(trialData.session_ID,strSesID);
	intSesTrials = sum(indSesTrials);
	
	vecTrialType = trialData.trialType(indSesTrials);
	vecTrialVisOriPre = trialData.visualOriPreChange(indSesTrials);
	vecTrialVisOriChange = trialData.visualOriChange(indSesTrials);
	vecTrialAudFreqPre = trialData.audioFreqPreChange(indSesTrials);
	vecTrialAudFreqChange = trialData.audioFreqChange(indSesTrials);
	vecTrialCorrect = trialData.correctResponse(indSesTrials);
	vecTrialRespSide = trialData.responseSide(indSesTrials);
	
	vecTrialLickSecs = cellfun(@rdivide,trialData.lickTime(indSesTrials),cellfill(dblTimeConversion,[intSesTrials 1]),'uniformoutput',false);
	
	vecTrialStartSecs = trialData.trialStart(indSesTrials)./dblTimeConversion;
	vecTrialChangeSecs = trialData.stimChange(indSesTrials)./dblTimeConversion;
	vecTrialEndSecs = trialData.trialEnd(indSesTrials)./dblTimeConversion;
	vecTrialStartITISecs = trialData.itiStart(indSesTrials)./dblTimeConversion;
	vecTrialEndITISecs = trialData.itiEnd(indSesTrials)./dblTimeConversion;
	vecTrialRespSecs = trialData.responseLatency(indSesTrials)./dblTimeConversion;
	
	%sub-select trials that are large change, visual, and correct
	indKeepTrials = abs(vecTrialVisOriChange) == 90;
	dblCutOff = 0.4;
	dblFracCorr = sum(vecTrialCorrect(indKeepTrials))/sum(indKeepTrials);
	if dblFracCorr < dblCutOff
		fprintf('Pop %d, Correct response proportion is %.3f, which is under %.3f, skipping... [%s]\n',intPopIdx,dblFracCorr,dblCutOff,getTime);
		continue;
	end
	%message
	fprintf('Analyzing pop %d/%d; %d neurons, %d trials [%s]\n',intPopIdx,intPops,intCells,sum(indKeepTrials),getTime);
	
	%define 1st vs 2nd bump windows
	if strcmp(strAlignOn,'Change')
		vecT0 = vecTrialChangeSecs(indKeepTrials);
	end
	
	%build vectors
	intTrials = numel(vecT0);
	vecStartBaseT = vecT0 + dblStartBaseT;
	vecStopBaseT = vecT0 + dblStopBaseT;
	vecStartEp1T = vecT0 + dblStartEp1T;
	vecStopEp1T = vecT0 + dblStopEp1T;
	vecStartEp2T = vecT0 + dblStartEp2T;
	vecStopEp2T = vecT0 + dblStopEp2T;
	vecEdges = sort(cat(1,vecStartBaseT,vecStopBaseT,vecStartEp1T,vecStopEp1T,vecStartEp2T,vecStopEp2T));
	vecBaseEdges = sort(cat(1,vecStartBaseT,vecStopBaseT));
	vecEp1Edges = sort(cat(1,vecStartEp1T,vecStopEp1T));
	vecEp2Edges = sort(cat(1,vecStartEp2T,vecStopEp2T));
	%calculate response matrix
	matRespBase = nan(intCells,intTrials);
	matRespEp1 = nan(intCells,intTrials);
	matRespEp2 = nan(intCells,intTrials);
	for intNeuron=1:intCells
		vecBaseR = histcounts(cellSpikes{intNeuron},vecBaseEdges);
		matRespBase(intNeuron,:) = vecBaseR(1:2:end);
		vecEp1R = histcounts(cellSpikes{intNeuron},vecEp1Edges);
		matRespEp1(intNeuron,:) = vecEp1R(1:2:end);
		vecEp2R = histcounts(cellSpikes{intNeuron},vecEp2Edges);
		matRespEp2(intNeuron,:) = vecEp2R(1:2:end);
	end
	
	%% prep data
	matData = cat(1,matRespEp1,matRespEp2);
	vecCellArea = cat(1,ones(intCells,1),ones(intCells,1)*2);
	vecTrialStimType = ones(1,intTrials);
	
	%% go through iters
	fprintf('    Iter (tot: %d) : \n',intFullIters);
	for intIter=1:intFullIters
		%% msg
		fprintf('\b[%d]\n',intIter);
		
		%% get data splits
		%set parameters
		dblCutOff = 0.9;
		sParams=struct;
		intUseSize = intCells;%16%floor(intCells/2);
		sParams.intSizeX = intUseSize;%floor(intCells/2);
		sParams.intSizeY = intUseSize;%floor(intCells/2);
		sParams.intResamplings = intResamplings;
		sParams.vecCellArea = vecCellArea;
		sParams.intWithinArea = 3;%0, across; 1, within 1; 2, within 2
		sParams.vecUseStimTypes = 1;
		cellStr = {'B1','B2','B1-B2'};
		strC = cellStr{sParams.intWithinArea};
		
		%get splits
		[cellMatX1,cellMatX2,cellNeuronsX,cellMatY1,cellMatY2,cellNeuronsY,cellTrials12,cellTrials1,cellTrials2] = ...
			doDimDataSplitsCV(matData,vecTrialStimType,sParams);
		%[cellMatX,	cellNeuronsX,cellMatY,cellNeuronsY,cellTrials] = doDimDataSplits(matData,vecTrialStimType,sParams);
		
		%% analyze
		intMaxDim = sParams.intSizeY;
		sOut = doDimCompareSpaceCV(cellMatX1,cellMatY1,cellMatX2,cellMatY2,intMaxDim,dblLambda,boolShuffle);
		%sOut = doDimCompareSpace(cellMatX,cellMatY,cellMatX,cellMatY,intMaxDim,dblLambda);
		
		%% save output
		vecM_xx2x = squeeze(nanmean(sOut.matR2_xx2x,1));
		vecM_xx2y = squeeze(nanmean(sOut.matR2_xx2y,1));
		vecM_xRRR2y = squeeze(nanmean(sOut.matR2_xRRR2y,1));
		vecM_xy2y = squeeze(nanmean(sOut.matR2_xy2y,1));
		
		vecM_yy2y = squeeze(nanmean(sOut.matR2_yy2y,1));
		vecM_yy2x = squeeze(nanmean(sOut.matR2_yy2x,1));
		vecM_yRRR2x = squeeze(nanmean(sOut.matR2_yRRR2x,1));
		vecM_yx2x = squeeze(nanmean(sOut.matR2_yx2x,1));
		
		
		cellM_xx2x{intIter,intPopIdx} = vecM_xx2x;
		cellM_xx2y{intIter,intPopIdx} = vecM_xx2y;
		cellM_xRRR2y{intIter,intPopIdx} = vecM_xRRR2y;
		cellM_xy2y{intIter,intPopIdx} = vecM_xy2y;
		
		cellM_yy2y{intIter,intPopIdx} = vecM_yy2y;
		cellM_yy2x{intIter,intPopIdx} = vecM_yy2x;
		cellM_yRRR2x{intIter,intPopIdx} = vecM_yRRR2x;
		cellM_yx2x{intIter,intPopIdx} = vecM_yx2x;
	end
end
warning('on','MATLAB:nearlySingularMatrix')

%% process data
cellM2_xx2x = cell(1,intPops);
cellM2_xx2y = cell(1,intPops);
cellM2_xRRR2y = cell(1,intPops);
cellM2_xy2y = cell(1,intPops);

cellM2_yy2y = cell(1,intPops);
cellM2_yy2x = cell(1,intPops);
cellM2_yRRR2x = cell(1,intPops);
cellM2_yx2x = cell(1,intPops);

cellS2_xx2x = cell(1,intPops);
cellS2_xx2y = cell(1,intPops);
cellS2_xRRR2y = cell(1,intPops);
cellS2_xy2y = cell(1,intPops);

cellS2_yy2y = cell(1,intPops);
cellS2_yy2x = cell(1,intPops);
cellS2_yRRR2x = cell(1,intPops);
cellS2_yx2x = cell(1,intPops);

% AUCs
matAUC_xx2x = nan(intFullIters,intPops);
matAUC_xx2y = nan(intFullIters,intPops);
matAUC_xRRR2y = nan(intFullIters,intPops);
matAUC_xy2y = nan(intFullIters,intPops);

matAUC_yy2y = nan(intFullIters,intPops);
matAUC_yy2x = nan(intFullIters,intPops);
matAUC_yRRR2x = nan(intFullIters,intPops);
matAUC_yx2x = nan(intFullIters,intPops);


for intPopIdx=1:intPops
	for intIter=1:intFullIters
		if isempty(cellM_xx2x{intIter,intPopIdx}),continue;end
		vecXX2X = cellM_xx2x{intIter,intPopIdx};
		matAUC_xx2x(intIter,intPopIdx) = sum((vecXX2X)/vecXX2X(end)) / numel(vecXX2X);
		
		vecXX2Y = cellM_xx2y{intIter,intPopIdx};
		matAUC_xx2y(intIter,intPopIdx) = sum((vecXX2Y)/vecXX2Y(end)) / numel(vecXX2Y);
		
		vecXR2Y = cellM_xRRR2y{intIter,intPopIdx};
		matAUC_xRRR2y(intIter,intPopIdx) = sum((vecXR2Y)/vecXR2Y(end)) / numel(vecXR2Y);
		
		vecXY2Y = cellM_xy2y{intIter,intPopIdx};
		matAUC_xy2y(intIter,intPopIdx) = sum((vecXY2Y)/vecXY2Y(end)) / numel(vecXY2Y);
		
		
		vecYY2Y = cellM_yy2y{intIter,intPopIdx};
		matAUC_yy2y(intIter,intPopIdx) = sum((vecYY2Y)/vecYY2Y(end)) / numel(vecYY2Y);
		
		vecYY2X = cellM_yy2x{intIter,intPopIdx};
		matAUC_yy2x(intIter,intPopIdx) = sum((vecYY2X)/vecYY2X(end)) / numel(vecYY2X);
		
		vecYR2X = cellM_yRRR2x{intIter,intPopIdx};
		matAUC_yRRR2x(intIter,intPopIdx) = sum((vecYR2X)/vecYR2X(end)) / numel(vecYR2X);
		
		vecYX2X = cellM_yx2x{intIter,intPopIdx};
		matAUC_yx2x(intIter,intPopIdx) = sum((vecYX2X)/vecYX2X(end)) / numel(vecYX2X);
	end
	%mean curves
	%process data
	cellM2_xx2x{intPopIdx} = mean(cell2mat(cellM_xx2x(:,intPopIdx)'),2);
	cellS2_xx2x{intPopIdx} = std(cell2mat(cellM_xx2x(:,intPopIdx)'),[],2)./sqrt(intFullIters);
	
	cellM2_xx2y{intPopIdx} = mean(cell2mat(cellM_xx2y(:,intPopIdx)'),2);
	cellS2_xx2y{intPopIdx} = std(cell2mat(cellM_xx2y(:,intPopIdx)'),[],2)./sqrt(intFullIters);
	
	cellM2_xRRR2y{intPopIdx} = mean(cell2mat(cellM_xRRR2y(:,intPopIdx)'),2);
	cellS2_xRRR2y{intPopIdx} = std(cell2mat(cellM_xRRR2y(:,intPopIdx)'),[],2)./sqrt(intFullIters);
	
	cellM2_xy2y{intPopIdx} = mean(cell2mat(cellM_xy2y(:,intPopIdx)'),2);
	cellS2_xy2y{intPopIdx} = std(cell2mat(cellM_xy2y(:,intPopIdx)'),[],2)./sqrt(intFullIters);
	
	cellM2_yy2y{intPopIdx} = mean(cell2mat(cellM_yy2y(:,intPopIdx)'),2);
	cellS2_yy2y{intPopIdx} = std(cell2mat(cellM_yy2y(:,intPopIdx)'),[],2)./sqrt(intFullIters);
	
	cellM2_yy2x{intPopIdx} = mean(cell2mat(cellM_yy2x(:,intPopIdx)'),2);
	cellS2_yy2x{intPopIdx} = std(cell2mat(cellM_yy2x(:,intPopIdx)'),[],2)./sqrt(intFullIters);
	
	cellM2_yRRR2x{intPopIdx} = mean(cell2mat(cellM_yRRR2x(:,intPopIdx)'),2);
	cellS2_yRRR2x{intPopIdx} = std(cell2mat(cellM_yRRR2x(:,intPopIdx)'),[],2)./sqrt(intFullIters);
	
	cellM2_yx2x{intPopIdx} = mean(cell2mat(cellM_yx2x(:,intPopIdx)'),2);
	cellS2_yx2x{intPopIdx} = std(cell2mat(cellM_yx2x(:,intPopIdx)'),[],2)./sqrt(intFullIters);
end

%% plot
%prep plot
%close all;
figure
vecAxPtr=nan(1,8);
for intPlot=1:8
	vecAxPtr(intPlot) = subplot(2,4,intPlot);
	hold on;
end
matColors = [0.1 0.1 0.8;...
	0.1 0.5 0.8;...
	0.5 0.1 0.8;...
	0.4 0.4 0.8;...
	0.8 0.1 0.1;...
	0.8 0.5 0.1;...
	0.8 0.1 0.5;...
	0.8 0.4 0.4];

%labels
cellLegend = {'A) B1(B1)-B1',...
	'B) B1(B1)-B2',...
	'C) B1(RRR)-B2',...
	'D) B1(B2)-B2'...
	'E) B2(B2)-B2',...
	'F) B2(B2)-B1',...
	'G) B2(RRR)-B1',...
	'H) B2(B1)-B1'};

cellLabelsX = {['B1 space, sorted by B1-' getGreek('lambda','lower')],...
	['B1 space, sorted by B1-' getGreek('lambda','lower')],...
	['B1 predictive dims'],...
	['B1 space, sorted by B2-' getGreek('lambda','lower')],...
	['B2 space, sorted by B2-' getGreek('lambda','lower')],...
	['B2 space, sorted by B2-' getGreek('lambda','lower')],...
	['B2 predictive dims'],...
	['B2 space, sorted by B1-' getGreek('lambda','lower')]};

cellLabelsY = {['Fraction of B1 activity (R^2)'],...
	['Fraction of B2 activity (R^2)'],...
	['Fraction of B2 activity (R^2)'],...
	['Fraction of B2 activity (R^2)'],...
	['Fraction of B2 activity (R^2)'],...
	['Fraction of B1 activity (R^2)'],...
	['Fraction of B1 activity (R^2)'],...
	['Fraction of B1 activity (R^2)']};

% plot with error bars
%{
for intPopIdx=1:intPops
	vecX = (1:numel(cellM2_xx2x{intPopIdx}))';%1:intMaxDim;%linspace(0,1,intMaxDim);
	errorbar(vecAxPtr(1),vecX,cellM2_xx2x{intPopIdx},cellS2_xx2x{intPopIdx},'Color',matColors(1,:));
	errorbar(vecAxPtr(2),vecX,cellM2_xx2y{intPopIdx},cellS2_xx2y{intPopIdx},'Color',matColors(2,:));
	errorbar(vecAxPtr(3),vecX,cellM2_xRRR2y{intPopIdx},cellS2_xRRR2y{intPopIdx},'Color',matColors(3,:));
	errorbar(vecAxPtr(4),vecX,cellM2_xy2y{intPopIdx},cellS2_xy2y{intPopIdx},'Color',matColors(4,:));
	
	errorbar(vecAxPtr(5),vecX,cellM2_yy2y{intPopIdx},cellS2_yy2y{intPopIdx},'Color',matColors(5,:));
	errorbar(vecAxPtr(6),vecX,cellM2_yy2x{intPopIdx},cellS2_yy2x{intPopIdx},'Color',matColors(6,:));
	errorbar(vecAxPtr(7),vecX,cellM2_yRRR2x{intPopIdx},cellS2_yRRR2x{intPopIdx},'Color',matColors(7,:));
	errorbar(vecAxPtr(8),vecX,cellM2_yx2x{intPopIdx},cellS2_yx2x{intPopIdx},'Color',matColors(8,:));
end
%}
% plot only means
for intPopIdx=1:intPops
	vecX = (1:numel(cellM2_xx2x{intPopIdx}))';%1:intMaxDim;%linspace(0,1,intMaxDim);
	plot(vecAxPtr(1),vecX,cellM2_xx2x{intPopIdx},'Color',matColors(1,:));
	plot(vecAxPtr(2),vecX,cellM2_xx2y{intPopIdx},'Color',matColors(2,:));
	plot(vecAxPtr(3),vecX,cellM2_xRRR2y{intPopIdx},'Color',matColors(3,:));
	plot(vecAxPtr(4),vecX,cellM2_xy2y{intPopIdx},'Color',matColors(4,:));
	
	plot(vecAxPtr(5),vecX,cellM2_yy2y{intPopIdx},'Color',matColors(5,:));
	plot(vecAxPtr(6),vecX,cellM2_yy2x{intPopIdx},'Color',matColors(6,:));
	plot(vecAxPtr(7),vecX,cellM2_yRRR2x{intPopIdx},'Color',matColors(7,:));
	plot(vecAxPtr(8),vecX,cellM2_yx2x{intPopIdx},'Color',matColors(8,:));
end

% finish plots
for intPlot=1:8
	axes(vecAxPtr(intPlot));
	hold off
	%ylim([0 max()]);
	title(cellLegend{intPlot});
	fixfig;
	xlabel(cellLabelsX{intPlot});
	ylabel(cellLabelsY{intPlot});
	if mod(intPlot,4) == 1
		ylim([0 0.3]);
	else
		ylim([0 0.05]);
	end
end

%save
maxfig();
if boolSavePlots
	strFigFile = ['DimCom1ST_MultiPredR2' getDate];
	export_fig([strFigDir strFigFile '.tif']);
	export_fig([strFigDir strFigFile '.pdf']);
end

% summary figure
figure
intUsedPops = sum(any(~isnan(matAUC_xx2x),1));
vecMeans = [nanmean(matAUC_xx2x(:)) ...
	nanmean(matAUC_xx2y(:)) ...
	nanmean(matAUC_xRRR2y(:)) ...
	nanmean(matAUC_xy2y(:)),...
	nanmean(matAUC_yy2y(:)) ...
	nanmean(matAUC_yy2x(:)) ...
	nanmean(matAUC_yRRR2x(:)) ...
	nanmean(matAUC_yx2x(:))];

vecSEMs = [nanstd(matAUC_xx2x(:)) ...
	nanstd(matAUC_xx2y(:)) ...
	nanstd(matAUC_xRRR2y(:)) ...
	nanstd(matAUC_xy2y(:)) ...
	nanstd(matAUC_yy2y(:)) ...
	nanstd(matAUC_yy2x(:)) ...
	nanstd(matAUC_yRRR2x(:)) ...
	nanstd(matAUC_yx2x(:)) ...
	]./sqrt(intUsedPops);

%ylim([0 0.7]);
ylim([0 0.5])
set(gca,'ytick',0:0.1:0.7);
xlim([0.5 8.5]);
hold on
for i=1:8
	errorbar(i,(1-vecMeans(i))*2,vecSEMs(i),'x','Color',matColors(i,:));
end
set(gca,'xtick',1:8,'xticklabel',cellLegend);
ylabel('Dimensionality index (2*(1-AUC))');
xtickangle(45)
fixfig;

if boolSavePlots
	strFigFile = ['DimCom1ST_DimensionalitySummary' getDate];
	export_fig([strFigDir strFigFile '.tif']);
	export_fig([strFigDir strFigFile '.pdf']);
end



