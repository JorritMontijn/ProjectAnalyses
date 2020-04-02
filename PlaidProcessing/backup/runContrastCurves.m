%{
Date		Type	BiDir			Plaid			Plane	Notes
20130612	SR101	xyt02			xyt01			1
20130625	PV		xyt02			xyt01			1
20131016	PV		xyt01/03		xyt02[/04:LQ]	1		xyt04 low qual
20131016	PV		xyt05							2		
20131016	PV		xyt06							3
20131022	PV		xyt01			xyt02			1
20131022	PV		xyt03							2
20140129	None	xyt01/03/04		xyt02			1		No SR101
20140129	None	xyt05							2		No SR101
20140314	SR101	xyt08			xyt09			1
20140423	Both	xyt01/03		xyt02			1		PV and SR101!
20140425	SR101	xyt09			xyt10			1
%}

strMasterPath = 'D:\Data\Processed\imagingdata\';
strDate = '20130612';
strRecGratings = 'xyt02';
strRecPlaids = 'xyt01';


%load gratings
strTargetDirGratings = [strMasterPath strDate filesep strRecGratings filesep];
strTargetFileGratings = [strDate strRecGratings '_ses.mat'];
sLoad = load([strTargetDirGratings strTargetFileGratings]);
sesGratings = sLoad.ses;

%load plaids
strTargetDirPlaids = [strMasterPath strDate filesep strRecPlaids filesep];
strTargetFilePlaids = [strDate strRecPlaids '_ses.mat'];
sLoad = load([strTargetDirPlaids strTargetFilePlaids]);
sesPlaids = sLoad.ses;

% switch to PV
%sesGratings.neuron = sesGratings.PV;
%sesPlaids.neuron = sesPlaids.PV;

%% transform gratings to resp mat
matRespG = getRespMat(sesGratings); %G=Gratings
vecOriG = getOriListFromTrials(sesGratings.structStim.Orientation);
sTypesG = getStimulusTypes(sesGratings.structStim);
cellSelectG = getSelectionVectors(sesGratings.structStim,sTypesG);
sTuningG = calcTuningRespMat(matRespG,cellSelectG,vecOriG);

%% transform resp mat to stimresp
matStimResponse = doMatRespTransform(matRespG,cellSelectG);
matMeanResp = squeeze(xmean(matStimResponse,2));

%% transform plaids to resp mat
matRespP = getRespMat(sesPlaids); %G=Gratings
vecOrisP = getOriListFromTrials(sesPlaids.structStim.Orientation);
vecContrastsP = [25 50];
sTypesP = getStimulusTypes(sesPlaids.structStim);
cellSelectP = getSelectionVectors(sesPlaids.structStim,sTypesP);

%% get baseline
structParams.intStartOffset = round(2*sesPlaids.samplingFreq); %remove 2 seconds after stim offset
indTrialStart = logical([1 diff(sesPlaids.structStim.TrialNumber)]);
matRespBaseP = getRespMat(sesPlaids,[],-2,structParams); %use pre-stim period
matRespBaseP = matRespBaseP(:,indTrialStart);

%% get trials
indOri1 = sesPlaids.structStim.Orientation==vecOrisP(1);
indOri2 = sesPlaids.structStim.Orientation==vecOrisP(2);
indCon1 = sesPlaids.structStim.Contrast==vecContrastsP(1);
indCon2 = sesPlaids.structStim.Contrast==vecContrastsP(2);


vecTotContrastAtOri1 = accumarray(sesPlaids.structStim.TrialNumber',(sesPlaids.structStim.Contrast.*indOri1)');
vecTotContrastAtOri2 = accumarray(sesPlaids.structStim.TrialNumber',(sesPlaids.structStim.Contrast.*indOri2)');

indGratingTrialsC100AtOri1 = cellSelectG{vecOriG==vecOrisP(1)};
indGratingTrialsC100AtOri2 = cellSelectG{vecOriG==vecOrisP(2)};


%% plot contrast response curve
vecContrasts = [0 25 50 100];
intNeurons = numel(sesPlaids.neuron);
matMeanConts1 = nan(intNeurons,4);
matStdsConts1 = nan(intNeurons,4);
matMeanConts2 = nan(intNeurons,4);
matStdsConts2 = nan(intNeurons,4);

for intNeuron=1:intNeurons
	%ori 1
vecRespCon00Ori1 = matRespP(intNeuron,vecTotContrastAtOri2==0);
vecRespCon25Ori1 = matRespP(intNeuron,vecTotContrastAtOri1==25&vecTotContrastAtOri2==0);
vecRespCon50Ori1 = matRespP(intNeuron,vecTotContrastAtOri1==50&vecTotContrastAtOri2==0);
vecRespCon100Ori1 = matRespG(intNeuron,indGratingTrialsC100AtOri1);
matMeanConts1(intNeuron,:) = [mean(vecRespCon00Ori1) mean(vecRespCon25Ori1) mean(vecRespCon50Ori1) mean(vecRespCon100Ori1)];
matStdsConts1(intNeuron,:) = [std(vecRespCon00Ori1) std(vecRespCon25Ori1) std(vecRespCon50Ori1) std(vecRespCon100Ori1)];

%ori 2
vecRespCon00Ori2 = matRespP(intNeuron,vecTotContrastAtOri1==0);
vecRespCon25Ori2 = matRespP(intNeuron,vecTotContrastAtOri2==25&vecTotContrastAtOri1==0);
vecRespCon50Ori2 = matRespP(intNeuron,vecTotContrastAtOri2==50&vecTotContrastAtOri1==0);
vecRespCon100Ori2 = matRespG(intNeuron,indGratingTrialsC100AtOri2);

matMeanConts2(intNeuron,:) = [mean(vecRespCon00Ori2) mean(vecRespCon25Ori2) mean(vecRespCon50Ori2) mean(vecRespCon100Ori2)];
matStdsConts2(intNeuron,:) = [std(vecRespCon00Ori2) std(vecRespCon25Ori2) std(vecRespCon50Ori2) std(vecRespCon100Ori2)];

%plot
subplot(2,4,1)
errorbar(vecContrasts,matMeanConts1(intNeuron,:),matStdsConts1(intNeuron,:)/sqrt(numel(vecRespCon25Ori2)));

subplot(2,4,2)
errorbar(vecContrasts,matMeanConts2(intNeuron,:),matStdsConts2(intNeuron,:)/sqrt(numel(vecRespCon25Ori2)));
pause
end

