%create carandini figure

%stimuli:
%{
		T1(3s)		T2(3)		T3(3s)		T4(3s)		Tot dur (+mask(1s)+ITI(5s))
	1)	Low on		High on		high off	low off			18s
	2)	High on		high off								12s
	3)	Low on								low off			18s
	4)	High on		low on		low off		high off		18s
	5)	low on		low off									12s
	6)	High on								high off		18s
	7)	H+L on		high off	low off						15s
	8)	H+L on		low off		high off					15s
	9)	Low on		low on		low off		low off			18s
	10)	High on		High on		high off	high off		18s
%}

%{
%%%% step 1

load data



%%%% step 2

%reformat data; output:
sTrialType(intTrialType).sPrimOri(intOriType).sEpoch(intEpoch).matFrames(intFrame,intRepetition)

intTrialType = [1:10]; trial types as above
intOriType = [1 2]; 1=horizontal; 2=vertical
intEpoch = [1:6]; 1-4: T1-4; 5=mask; 6=iti
intFrame = frame number (1 - max in epoch)
intRepetition = repetition number (1-4)

matData: matrix containing dF/F values
matFrames: matrix containing corresponding original frame numbers



%%%% step 3


%}

%% settings
dblPercentageCells = 0.5;


intPrefDirType = 3;
if intPrefDirType == 1
	maxOri = 180;
	intIncrement = 45/4;
	binEdges = (-intIncrement/2):intIncrement:(maxOri-(0.5*intIncrement));
	binCenters = 0:intIncrement:(maxOri-(0.5*intIncrement));
elseif intPrefDirType == 2
	maxOri = 180;
	intIncrement = 45;
	binEdges = (-intIncrement/2):intIncrement:(maxOri-(0.5*intIncrement));
	binCenters = 0:intIncrement:(maxOri-(0.5*intIncrement));
elseif intPrefDirType == 3
	maxOri = 360;
	intIncrement = 45;
	binEdges = (-intIncrement/2):intIncrement:(maxOri-(0.5*intIncrement));
	binCenters = 0:intIncrement:(maxOri-(0.5*intIncrement));
elseif intPrefDirType == 4
	
end

%% load data
if ~exist('sesPlaids','var')
	%load plaid data
	strPlaidSes = 'D:\Data\Processed\imagingdata\20130612\xyt01\20130612xyt01_ses.mat';
	sLoad = load(strPlaidSes);
	sesPlaids = sLoad.ses;
	clear sLoad;
end
if ~exist('sesGratings','var')
	%load ori tuning data
	strTuningSes = 'D:\Data\Processed\imagingdata\20130612\xyt02\20130612xyt02_ses.mat';
	sLoad = load(strTuningSes);
	sesGratings = sLoad.ses;
	clear sLoad;
end

%% pre-allocate
intNeurons = numel(sesGratings.neuron);
vecNormV100H0 = nan(1,intNeurons);
vecV50H25 = nan(1,intNeurons);
vecV50H0 = nan(1,intNeurons);
vecV25H50 = nan(1,intNeurons);
vecV25H0 = nan(1,intNeurons);
vecNormV0H100 = nan(1,intNeurons);
vecV0H50 = nan(1,intNeurons);
vecV0H25 = nan(1,intNeurons);
vecPrefOri = nan(1,intNeurons);
vecPrefDir = nan(1,intNeurons);
vecOSI = nan(1,intNeurons);

%% get tuning data
%get lookups
vecStimTypeLookup = round(getOriListFromTrials(sesGratings.structStim.Orientation));
vecStimTypeIndex = 1:length(vecStimTypeLookup);

for intNeuron=1:intNeurons
	%split activation by orientation
	[structActivity,matAct_TypeByRep,vecStimTypeLookup,sOut] = calcTuning(sesGratings,intNeuron);
	vecMean = mean(sOut.matStimAct_TypeByRep,2);
	
	%get activity at horizontal stimulus
	dblActH = mean([vecMean(vecStimTypeLookup == 0) vecMean(vecStimTypeLookup == 180)]);
	dblActH1 = vecMean(vecStimTypeLookup == 0);
	dblActH2 = dblActH1;%vecMean(vecStimTypeLookup == 180);
	
	%get activity at vertical stimulus
	dblActV = mean([vecMean(vecStimTypeLookup == 90) vecMean(vecStimTypeLookup == 270)]);
	dblActV1 = vecMean(vecStimTypeLookup == 90);
	dblActV2 = dblActV1;%vecMean(vecStimTypeLookup == 270);
	
	%normalize
	vecActNorm = imnorm([dblActH1 dblActV1 dblActH2 dblActV2]);
	
	%assign to vectors
	vecNormV100H0(intNeuron) = mean([vecActNorm(2) vecActNorm(4)]);
	vecNormV0H100(intNeuron) = mean([vecActNorm(1) vecActNorm(3)]);
end
sTuning = calcOriTuning(sesGratings);
vecPrefOri = sTuning.vecPrefAngle;
vecPrefDir = vecPrefOri;
vecOSI = sTuning.vecOSI;

%get neurons tuned strongest
[c,vecCells] = findmax(vecOSI,round(dblPercentageCells*intNeurons));

%% reformat
sTrialType = doReformatPlaidData(sesPlaids);

%% get applicable stimulus types and activity levels

%{
V\H		100				50				25				0
100		-				-				-				OT-V
50		-				-				7-V+8-H			2+4+6+10-V
25		-				7-H + 8-V		-				1+3+5+9-V
0		OT-H			2+4+6+10-H		1+3+5+9-H		-
%}

%1 hor, 2 vert
vecHorContrast = [-1 100 50 25 0];
vecVertContrast = [-1 100 50 25 0];
intSizeH = length(vecHorContrast);
intSizeV = length(vecVertContrast);
for intHorInd = 1:intSizeH
	for intVertInd = 1:intSizeV
		%get contrasts
		intHorCont = vecHorContrast(intHorInd);
		intVertCont = vecVertContrast(intVertInd);
		
		if intVertCont == 50 && intHorCont == 25
			for intNeuron=1:intNeurons
				%7-V+8-H
				
				%7-V
				vecFrames1 = sTrialType(7).sPrimOri(2).sEpoch(1).matFrames(:);
				vecAct1 = sesPlaids.neuron(intNeuron).dFoF(vecFrames1);
				
				%8-H
				vecFrames2 = sTrialType(8).sPrimOri(1).sEpoch(1).matFrames(:);
				vecAct2 = sesPlaids.neuron(intNeuron).dFoF(vecFrames2);
				
				%mean
				vecV50H25(intNeuron) = mean([vecAct1(:); vecAct2(:)]);
			end
		elseif intVertCont == 50 && intHorCont == 0
			
			for intNeuron=1:intNeurons
				%2+4+6+10-V
				
				%2-V
				vecFrames1 = sTrialType(2).sPrimOri(2).sEpoch(1).matFrames(:);
				vecAct1 = sesPlaids.neuron(intNeuron).dFoF(vecFrames1);
				
				%4-V
				vecFrames2 = sTrialType(4).sPrimOri(2).sEpoch(1).matFrames(:);
				vecAct2 = sesPlaids.neuron(intNeuron).dFoF(vecFrames2);
				
				%6-V
				vecFrames3 = sTrialType(6).sPrimOri(2).sEpoch(1).matFrames(:);
				vecAct3 = sesPlaids.neuron(intNeuron).dFoF(vecFrames3);
				
				%10-V
				vecFrames4 = sTrialType(10).sPrimOri(2).sEpoch(1).matFrames(:);
				vecAct4 = sesPlaids.neuron(intNeuron).dFoF(vecFrames4);
				
				%mean
				vecV50H0(intNeuron) = mean([vecAct1(:); vecAct2(:); vecAct3(:); vecAct4(:)]);
			end
		elseif intVertCont == 25 && intHorCont == 50
			%7-H + 8-V
			for intNeuron=1:intNeurons
				%7-H+8-V
				
				%7-H
				vecFrames1 = sTrialType(7).sPrimOri(1).sEpoch(1).matFrames(:);
				vecAct1 = sesPlaids.neuron(intNeuron).dFoF(vecFrames1);
				
				%8-V
				vecFrames2 = sTrialType(8).sPrimOri(2).sEpoch(1).matFrames(:);
				vecAct2 = sesPlaids.neuron(intNeuron).dFoF(vecFrames2);
				
				%mean
				vecV25H50(intNeuron) = mean([vecAct1(:); vecAct2(:)]);
			end
		elseif intVertCont == 25 && intHorCont == 0
			%1+3+5+9-V
			for intNeuron=1:intNeurons
				%1+3+5+9-V
				
				%1-V
				vecFrames1 = sTrialType(1).sPrimOri(2).sEpoch(1).matFrames(:);
				vecAct1 = sesPlaids.neuron(intNeuron).dFoF(vecFrames1);
				
				%3-V
				vecFrames2 = sTrialType(3).sPrimOri(2).sEpoch(1).matFrames(:);
				vecAct2 = sesPlaids.neuron(intNeuron).dFoF(vecFrames2);
				
				%5-V
				vecFrames3 = sTrialType(5).sPrimOri(2).sEpoch(1).matFrames(:);
				vecAct3 = sesPlaids.neuron(intNeuron).dFoF(vecFrames3);
				
				%9-V
				vecFrames4 = sTrialType(9).sPrimOri(2).sEpoch(1).matFrames(:);
				vecAct4 = sesPlaids.neuron(intNeuron).dFoF(vecFrames4);
				
				%mean
				vecV25H0(intNeuron) = mean([vecAct1(:); vecAct2(:); vecAct3(:); vecAct4(:)]);
			end
		elseif intVertCont == 0 && intHorCont == 50
			%2+4+6+10-H
			for intNeuron=1:intNeurons
				%2+4+6+10-H
				
				%2-H
				vecFrames1 = sTrialType(2).sPrimOri(1).sEpoch(1).matFrames(:);
				vecAct1 = sesPlaids.neuron(intNeuron).dFoF(vecFrames1);
				
				%4-H
				vecFrames2 = sTrialType(4).sPrimOri(1).sEpoch(1).matFrames(:);
				vecAct2 = sesPlaids.neuron(intNeuron).dFoF(vecFrames2);
				
				%6-H
				vecFrames3 = sTrialType(6).sPrimOri(1).sEpoch(1).matFrames(:);
				vecAct3 = sesPlaids.neuron(intNeuron).dFoF(vecFrames3);
				
				%10-H
				vecFrames4 = sTrialType(10).sPrimOri(1).sEpoch(1).matFrames(:);
				vecAct4 = sesPlaids.neuron(intNeuron).dFoF(vecFrames4);
				
				%mean
				vecV0H50(intNeuron) = mean([vecAct1(:); vecAct2(:); vecAct3(:); vecAct4(:)]);
			end
		elseif intVertCont == 0 && intHorCont == 25
			%1+3+5+9-H
			for intNeuron=1:intNeurons
				%1+3+5+9-H
				
				%1-H
				vecFrames1 = sTrialType(1).sPrimOri(1).sEpoch(1).matFrames(:);
				vecAct1 = sesPlaids.neuron(intNeuron).dFoF(vecFrames1);
				
				%3-H
				vecFrames2 = sTrialType(3).sPrimOri(1).sEpoch(1).matFrames(:);
				vecAct2 = sesPlaids.neuron(intNeuron).dFoF(vecFrames2);
				
				%5-H
				vecFrames3 = sTrialType(5).sPrimOri(1).sEpoch(1).matFrames(:);
				vecAct3 = sesPlaids.neuron(intNeuron).dFoF(vecFrames3);
				
				%9-H
				vecFrames4 = sTrialType(9).sPrimOri(1).sEpoch(1).matFrames(:);
				vecAct4 = sesPlaids.neuron(intNeuron).dFoF(vecFrames4);
				
				%mean
				vecV0H25(intNeuron) = mean([vecAct1(:); vecAct2(:); vecAct3(:); vecAct4(:)]);
			end
		end
	end
end

%% normalize activity per neuron
vecNormV50H25 = nan(1,intNeurons);
vecNormV50H0 = nan(1,intNeurons);
vecNormV25H50 = nan(1,intNeurons);
vecNormV25H0 = nan(1,intNeurons);
vecNormV0H50 = nan(1,intNeurons);
vecNormV0H25 = nan(1,intNeurons);

for intNeuron = 1:intNeurons
	vecRawAct = [vecV50H25(intNeuron) vecV50H0(intNeuron) vecV25H50(intNeuron) vecV25H0(intNeuron) vecV0H50(intNeuron) vecV0H25(intNeuron)];
	vecNormAct = imnorm(vecRawAct);
	
	vecNormV50H25(intNeuron) = vecNormAct(1);
	vecNormV50H0(intNeuron) = vecNormAct(2);
	vecNormV25H50(intNeuron) = vecNormAct(3);
	vecNormV25H0(intNeuron) = vecNormAct(4);
	vecNormV0H50(intNeuron) = vecNormAct(5);
	vecNormV0H25(intNeuron) = vecNormAct(6);
end



%% draw figure
h=figure;
set(h,'Color','w');
vecX = [binCenters binCenters(end)+intIncrement];
vecTickX = 0:45:max(binCenters);
vecLim = [min(vecX) max(vecX)];
%1 hor, 2 vert
for intHorInd = 1:intSizeH
	for intVertInd = 1:intSizeV
		%get plotnr
		intPlotNr = intHorInd + (intVertInd-1)*intSizeV;
		subplot(intSizeV,intSizeH,intPlotNr);
		
		%get contrasts
		intHorCont = vecHorContrast(intHorInd);
		intVertCont = vecVertContrast(intVertInd);
		
		%get contrast pics
		if intHorCont == -1 && intVertCont ~= -1
			[imTemp] = imread(sprintf('90_%d.png',intVertCont));
			subimage(imTemp);
			title(sprintf('Vert %d',intVertCont));
			axis off;
		elseif intVertCont == -1 && intHorCont ~= -1
			[imTemp] = imread(sprintf('0_%d.png',intHorCont));
			subimage(imTemp);
			title(sprintf('Hor %d',intHorCont));
			axis off;
		elseif intVertCont == 100 && intHorCont == 0
			%OT-V
			
			%[nVec,meanVec,stdVec] = makeBins(vecPrefOri(vecCells),vecNormV100H0(vecCells),binEdges);
			[nVec,meanVec,stdVec] = makeBins(vecPrefDir(vecCells),vecNormV100H0(vecCells),binEdges);
			
			plot(vecX,[meanVec meanVec(1)]);
			set(gca,'XTick',vecTickX);
			xlim(vecLim);
		elseif intVertCont == 50 && intHorCont == 25
			%7-V+8-H
			
			[nVec,meanVec,stdVec] = makeBins(vecPrefDir(vecCells),vecNormV50H25(vecCells),binEdges);
			
			plot(vecX,[meanVec meanVec(1)]);
			set(gca,'XTick',vecTickX);
			xlim(vecLim);
		elseif intVertCont == 50 && intHorCont == 0
			%2+4+6+10-V
			
			[nVec,meanVec,stdVec] = makeBins(vecPrefDir(vecCells),vecNormV50H0(vecCells),binEdges);
			
			plot(vecX,[meanVec meanVec(1)]);
			set(gca,'XTick',vecTickX);
			xlim(vecLim);
		elseif intVertCont == 25 && intHorCont == 50
			%7-H + 8-V
			
			[nVec,meanVec,stdVec] = makeBins(vecPrefDir(vecCells),vecNormV25H50(vecCells),binEdges);
			
			plot(vecX,[meanVec meanVec(1)]);
			set(gca,'XTick',vecTickX);
			xlim(vecLim);
		elseif intVertCont == 25 && intHorCont == 0
			%1+3+5+9-V
			
			[nVec,meanVec,stdVec] = makeBins(vecPrefDir(vecCells),vecNormV25H0(vecCells),binEdges);
			
			plot(vecX,[meanVec meanVec(1)]);
			set(gca,'XTick',vecTickX);
			xlim(vecLim);
		elseif intVertCont == 0 && intHorCont == 100
			%OT-H
			%[nVec,meanVec,stdVec] = makeBins(vecPrefOri(vecCells),vecNormV0H100(vecCells),binEdges);
			[nVec,meanVec,stdVec] = makeBins(vecPrefDir(vecCells),vecNormV0H100(vecCells),binEdges);
			
			plot(vecX,[meanVec meanVec(1)]);
			set(gca,'XTick',vecTickX);
			xlim(vecLim);
		elseif intVertCont == 0 && intHorCont == 50
			%2+4+6+10-H
			
			[nVec,meanVec,stdVec] = makeBins(vecPrefDir(vecCells),vecNormV0H50(vecCells),binEdges);
			
			plot(vecX,[meanVec meanVec(1)]);
			set(gca,'XTick',vecTickX);
			xlim(vecLim);
		elseif intVertCont == 0 && intHorCont == 25
			%1+3+5+9-H
			
			[nVec,meanVec,stdVec] = makeBins(vecPrefDir(vecCells),vecNormV0H25(vecCells),binEdges);
			
			plot(vecX,[meanVec meanVec(1)]);
			set(gca,'XTick',vecTickX);
			xlim(vecLim);
		else
			axis off;
		end
	end
end




