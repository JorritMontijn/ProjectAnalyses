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

% parameters for window of interest:
dblStartBaseT = -0.4;
dblStopBaseT = -0.2;
dblStartEp1T = 0;
dblStopEp1T = 0.2;
dblStartEp2T = 0.2;
dblStopEp2T = 0.4;

%% load data
[Data] = MOL_GetData('F:\Data\Processed\BumpsMatthijs\','CHDET',{'ChangeDetectionConflict'},{'2003' '2009' '2010' '2011' '2012' '2013'},[],{'sessionData' 'trialData' 'spikeData'});
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
[sessionData,trialData,spikeData] = MOL_calc_signresp(sessionData,trialData,spikeData);

%% make populations
[cellPopID,ia,vecPopulationIdx] = unique(vecSessionIdx*10 + vecAreaIdx);

%% define V1 populations
indV1=vecAreaIdx==find(strcmp(cellAreas,'V1'));
vecPopsV1 = unique(vecPopulationIdx(indV1));
intPops = numel(vecPopsV1);
intCutOffCellNr = 16;

%% analyze
dblTimeConversion = 10^6;
for intPopIdx=1:intPops
	%% prepare, retrieve and transform data
	%retrieve cells
	intPopulation = vecPopsV1(intPopIdx);
	indUseCells = vecPopulationIdx==intPopulation;
	intCells = sum(indUseCells);
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
	
	dblCutOff = 0.4;
	dblFracCorr = sum(vecTrialCorrect)/numel(vecTrialCorrect);
	if dblFracCorr < dblCutOff%|| ismember(intPopIdx,vecRemoveSessions)
		fprintf('Pop %d, Correct response proportion is %.3f, which is under %.3f, skipping... [%s]\n',intPopIdx,dblFracCorr,dblCutOff,getTime);
		continue;
	end
	
	%define 1st vs 2nd bump windows
	if strcmp(strAlignOn,'Change')
		vecT0 = vecTrialChangeSecs;
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
	intIters = 10;
	matData = cat(1,matRespEp1,matRespEp2);
	vecCellArea = cat(1,ones(intCells,1),ones(intCells,1)*2);
	vecTrialStimTypes = ones(1,intTrials);
	
	%% get data splits
	%set parameters
	dblLambda = 500;
	dblCutOff = 0.9;
	sDataParams=struct;
	sDataParams.intSizeX = floor(intCells/2);
	sDataParams.intSizeY = floor(intCells/2);
	sDataParams.intResamplings = intIters;
	sDataParams.vecCellArea = vecCellArea;
	sDataParams.intWithinArea = 1;%0, across; 1, within 1; 2, within 2
	sDataParams.vecUseStimTypes = 1;
	cellStr = {'B1','B2','B1-B2'};
	strC = cellStr{sDataParams.intWithinArea};
	
	%get splits
	[cellMatX,cellNeuronsX,cellMatY,cellNeuronsY,cellTrials] = doDimDataSplits(matData,vecTrialStimTypes,sDataParams);
	fprintf('  Created %dx%d data splits [%s]\n',size(cellMatX),getTime);
	warning off;
	
	%% run predictions
	% predict with single neurons
	matPredictionsR2_SingleNeurons = doDimPredSingle(cellMatX,cellMatY,dblLambda);
	fprintf('  %s; Mean single neuron predictions: <%s\b> [%s]\n',strC,sprintf('%.4f ',xmean(xmean(matPredictionsR2_SingleNeurons,1),3)),getTime);
	
	% predict with full pop
	matPredictionsR2 = doDimPredFull(cellMatX,cellMatY,dblLambda);
	fprintf('  %s; Mean full population predictions: <%s\b> [%s]\n',strC,sprintf('%.4f ',xmean(matPredictionsR2,1)),getTime);
	
	% predict with reduced rank regression, normal
	matPredDimDepR2 = doDimPredRR(cellMatX,cellMatY);
	vecMeanR2 = mean(squeeze(matPredDimDepR2),1);
	vecSdR2 = std(squeeze(matPredDimDepR2),[],1);
	matPredictiveDimensions = sum(bsxfun(@rdivide,matPredDimDepR2,matPredDimDepR2(:,:,end))<dblCutOff,3);
	fprintf('  %s; Mean number of predictive dimensions: <%s\b> [%s]\n',strC,sprintf('%.2f ',xmean(matPredictiveDimensions,1)),getTime);
	
	% predict with removing predictive dimensions
	[matR2RemPredDim,intRemPredDimMax] = doDimPredRemPred(cellMatX,cellMatY,matPredictiveDimensions,matPredictionsR2,dblLambda);
	intRemPredDimMax = sDataParams.intSizeY;
	fprintf('  %s; Prediction before removal of pred dims: <%s\b>; After removal (d=%d): <%s\b> [%s]\n',strC,sprintf('%.4f ',xmean(matR2RemPredDim(:,:,1),1)),intRemPredDimMax,sprintf('%.4f ',xmean(matR2RemPredDim(:,:,end),1)),getTime);
	
	% predict with dominant dimensions
	matR2UseDomDim = doDimPredRemDom(cellMatX,cellMatY,intRemPredDimMax,dblLambda);
	fprintf('  %s; Prediction using first dominant dim: <%s\b>; Using d=%d: <%s\b> [%s]\n',strC,sprintf('%.4f ',xmean(matR2UseDomDim(:,:,1),1)),intRemPredDimMax,sprintf('%.4f ',xmean(matR2UseDomDim(:,:,end),1)),getTime);
	
	% predict with principal dimensions
	matR2UsePrincDim = doDimPredRemPrinc(cellMatX,cellMatY,intRemPredDimMax,dblLambda);
	fprintf('  %s; Prediction using largest eigenvector: <%s\b>; Using d=%d: <%s\b> [%s]\n',strC,sprintf('%.4f ',xmean(matR2UsePrincDim(:,:,1),1)),intRemPredDimMax,sprintf('%.4f ',xmean(matR2UsePrincDim(:,:,end),1)),getTime);
	
	%% plot
	vecX = 1:sDataParams.intSizeY;
	vecEndX = (vecX(end)+0.5)*ones([sDataParams.intResamplings 1]);
	figure
	subplot(2,2,1)
	plot(vecX,squeeze(matPredDimDepR2)','Color',[0.1 0.1 0.8]);
	hold on
	%scatter(vecEndX,matPredictionsR2)
	hold off
	vecLimY = [0 max(get(gca,'ylim'))];
	ylim(vecLimY);
	fixfig;
	xlabel('Predictive dim #');
	ylabel(sprintf('R^2 of %s',strC));
	
	subplot(2,2,2)
	plot(vecX(1:intRemPredDimMax),squeeze(matR2UseDomDim)','Color',[0.8 0.1 0.1]);
	ylim(vecLimY);
	fixfig;
	xlabel('Dominant FA dim #');
	ylabel(sprintf('R^2 of %s',strC));
	
	subplot(2,2,3)
	plot(vecX(1:intRemPredDimMax),squeeze(matR2UsePrincDim)','Color',[0.6 0.1 0.6]);
	ylim(vecLimY);
	fixfig;
	xlabel('Principal component #');
	ylabel(sprintf('R^2 of %s',strC));
	
end

















