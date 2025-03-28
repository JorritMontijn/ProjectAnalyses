%% description
%Parameters:
%Each column of the data set corresponds to a Cell with its respective temporal components
%�	Animal: 19, 20, 21, 22, 23, 24, 26.
%�	Cell_id: cell id of the spatial components
%�	Cell_num: number of the cell into the data frame
%�	Session: number of the session
%�	sigraw_n: Set of temporal components. Saved as a structure with dimensions (1x1896 cell). Each column corresponds to a cell component denoised and deconvolved.
%�	Sigraw_F_dff_n: Set of DF/F normalized temporal components. Saved as a structure with dimensions (1x1896 cell). Each row corresponds to the DF/F fluorescence for the corresponding component.
%
%Condition file: US_stim_S1_Condition.mat
%
%Parameters:
%�	Animal: 19, 20, 21, 22, 23, 24, 26.
%�	Condition1 (High Intensity): Set of temporal onsets of the stimulation.  Saved as an array with dimension (# animal x # onsets (in frames)) Each row corresponds to an animal and each column corresponds to an onset of the stimulation (20 stim per condition)
%�	Condition2 (Low Intensity): Set of temporal onsets of the stimulation.  Saved as an array with dimension (# animal x # onsets (in frames)) Each row corresponds to an animal and each column corresponds to an onset of the stimulation (20 stim per condition)
%�	Condition3 (Inhibitory condition): Set of temporal onsets of the stimulation.  Saved as an array with dimension (# animal x # onsets (in frames)) Each row corresponds to an animal and each column corresponds to an onset of the stimulation (20 stim per condition)
%�	Condition4 (Control): Set of temporal onsets of the stimulation.  Saved as an array with dimension (# animal x # onsets (in frames)) Each row corresponds to an animal and each column corresponds to an onset of the stimulation (20 stim per condition)

%% set recording
close all;
clear all;
strDataTargetPath = 'D:\Data\Processed\DataCamiloValeria\';
vecRunTypes = [1 2];
intResampNum = 500;
boolSave = true;%true;
strFigPath = 'D:\Data\Results\TraceZeta\ultrasound\';
strName = 'UltrasoundMod2Still';

% check if prepro data exists
sDir = dir([strDataTargetPath 'AggRecUS2.mat']);
if numel(sDir) == 1
	sLoad = load(fullpath(sDir(1).folder,sDir(1).name));
	sRec = sLoad.sRec;
	clear sLoad;
	fprintf('Loaded prepro data from %s [%s]\n',strDataTargetPath,getTime);
end
if ~exist('sRec','var') || ~isfield(sRec(end),'matTsZetaP')
	%% load data
	sLoadCond = load([strDataTargetPath 'US_stim_S1_Conditions.mat']);
	sLoadCa = load([strDataTargetPath 'US_stim_S1_Calcium.mat']);
	sLoadMot = load([strDataTargetPath 'US_stim_S1_Motion.mat']);
	
	%% split data by recording
	sTypes = getStimulusTypes(sLoadCond,{'Animal','Session'});
	cellSelectStims = getSelectionVectors(sLoadCond,sTypes);
	sRec = struct;
	cellSelectNeurons = getSelectionVectors(sLoadCa,sTypes);
	
	for intRec=1:numel(cellSelectStims)
		cellCond = cell(1,4);
		cellCond{1} = sLoadCond.Condition1(cellSelectStims{intRec},:);
		cellCond{2} = sLoadCond.Condition2(cellSelectStims{intRec},:);
		cellCond{3} = sLoadCond.Condition3(cellSelectStims{intRec},:);
		cellCond{4} = sLoadCond.Condition4(cellSelectStims{intRec},:);
		intTrials = sum(cellfun(@numel,cellCond));
		cellCondIdx = cell(size(cellCond));
		for intStimType=1:4
			cellCondIdx{intStimType} = intStimType*ones(size(cellCond{intStimType}));
		end
		vecStimOn = cell2vec(cellCond);
		vecStimType = cell2vec(cellCondIdx);
		%sort
		[vecStimOn,vecReorder]=sort(vecStimOn);
		vecStimType = vecStimType(vecReorder);
		%build stim struct
		structStim = struct;
		structStim.vecStimOn = vecStimOn;
		structStim.vecStimType = vecStimType;
		structStim.strLabel = 'Stim types; 1=High Intensity, 2=Low Intensity, 3=Inhibitory condition, 4=Control';
		
		%add stim
		sRec(intRec).structStim = structStim;
		sRec(intRec).SampFreq = 15;
		
		%select neurons
		vecNeurons = find(cellSelectNeurons{intRec});
		intNumN = numel(vecNeurons);
		neuron = struct;
		for intN=1:intNumN
			%select neuron
			intIdxN = vecNeurons(intN);
			neuron(intN).idx = intN;
			neuron(intN).origIdx = intIdxN;
			neuron(intN).cell_num = sLoadCa.Cell_num(intIdxN);
			neuron(intN).cell_id = sLoadCa.Cell_Id(intIdxN);
			neuron(intN).dFoF = sLoadCa.sigraw_F_dff_n{intIdxN};
			neuron(intN).deconvF = sLoadCa.sigraw_n{intIdxN};
		end
		
		%add neurons
		sRec(intRec).neuron = neuron;
		
		%add running
		sRec(intRec).Velocity = sLoadMot.Velocity{intRec}; %10cm/s
	end
	
	%% run analysis
	intPlot  = 0;
	hTic = tic;
	for intRec=1:numel(sRec)
		[varDataOut,vecUnique,vecCounts,cellSelect,vecRepetition] = val2idx(sRec(intRec).structStim.vecStimType);
		intStimTypes = numel(vecUnique);
		intNeurons = numel(sRec(intRec).neuron);
		matTsZetaP = nan(intStimTypes,intNeurons);
		matAnovaP = nan(intStimTypes,intNeurons);
		matWilcoxP = nan(intStimTypes,intNeurons);
		
		%select still trials
		vecVelocity = sRec(intRec).Velocity;
		intStimDur = sRec(intRec).SampFreq*2;
		dblUseMaxDur = floor(min(diff(sRec(intRec).structStim.vecStimOn / sRec(intRec).SampFreq)));
		vecPreBaseStart = sRec(intRec).structStim.vecStimOn - intStimDur - 0.5;
		vecPreBaseEnd = sRec(intRec).structStim.vecStimOn - 0.5;
		vecPreBaseEdges = flat([vecPreBaseStart'; vecPreBaseEnd']);
		vecStimStart = sRec(intRec).structStim.vecStimOn - 0.5;
		vecStimEnd = sRec(intRec).structStim.vecStimOn + intStimDur - 0.5;
		vecStimEdges = flat([vecStimStart'; vecStimEnd']);
		%calc
		[vecCounts,vecMeans,vecSDs,cellVals,cellIDs] = ...
			makeBins(1:numel(vecVelocity),vecVelocity,vecPreBaseEdges);
		vecPreBaseVel = vecMeans(1:2:end);
		[vecCounts,vecMeans,vecSDs,cellVals,cellIDs] = ...
			makeBins(1:numel(vecVelocity),vecVelocity,vecStimEdges);
		vecStimVel = vecMeans(1:2:end);
		vecTotVel = (vecPreBaseVel + vecStimVel)/2;
		
		%start times
		vecStartT = sRec(intRec).structStim.vecStimOn ./ sRec(intRec).SampFreq;
		
		for intNeuron=1:intNeurons
			if toc(hTic) > 5
				fprintf('Rec %d/%d, neuron %d/%d [%s]\n',intRec,numel(sRec),intNeuron,intNeurons,getTime);
				hTic=tic;
			end
			vec_dFoF = sRec(intRec).neuron(intNeuron).dFoF;
			vecTimestamps = (1:numel(vec_dFoF)) / sRec(intRec).SampFreq;
			dblEndSecs = (numel(vec_dFoF) - sRec(intRec).structStim.vecStimOn(end)) / sRec(intRec).SampFreq;
			dblUseMaxDur = min(dblUseMaxDur,dblEndSecs/2);
			%transform to resp matrix
			[matTE,vecWindowBinCenters] = getRespMat(vecTimestamps,vec_dFoF,vecStartT,[0 dblUseMaxDur]);
			%calc pre & post act
			[vecCounts,vecMeans,vecSDs,cellVals,cellIDs] = ...
				makeBins(1:numel(vec_dFoF),vec_dFoF,vecPreBaseEdges);
			vecPreBaseAct = vecMeans(1:2:end);
			[vecCounts,vecMeans,vecSDs,cellVals,cellIDs] = ...
				makeBins(1:numel(vec_dFoF),vec_dFoF,vecStimEdges);
			vecStimAct = vecMeans(1:2:end);
			
			for intStimType=1:intStimTypes
				vecUseTrials = find(sRec(intRec).structStim.vecStimType==intStimType);
				indStill = vecTotVel(vecUseTrials) < median(vecTotVel(vecUseTrials));
				vecUseTrials = vecUseTrials(indStill);
				
				vecStimOnTimes = vecStartT(vecUseTrials);
				%zeta
				dblZetaP=zetatstest(vecTimestamps,vec_dFoF,vecStimOnTimes,dblUseMaxDur,intResampNum,intPlot);
				matTsZetaP(intStimType,intNeuron) = dblZetaP;
				%anova
				matSubResp = matTE(vecUseTrials,:);
				dblAnovaP = anova1(matSubResp,[],'off');
				matAnovaP(intStimType,intNeuron) = dblAnovaP;
				%wilcoxon ranksum
				dblWilcoxP = ranksum(vecPreBaseAct(vecUseTrials),vecStimAct(vecUseTrials));
				matWilcoxP(intStimType,intNeuron) = dblWilcoxP;
				if dblZetaP < -inf%0.01
					dblZetaP2=zetatstest(vecTimestamps,vec_dFoF,vecStimOnTimes,dblUseMaxDur,intResampNum,2);
					disp intNeuron
				end
			end
		end
		sRec(intRec).matTsZetaP = matTsZetaP;
		sRec(intRec).matAnovaP = matAnovaP;
		sRec(intRec).matWilcoxP = matWilcoxP;
	end
	%save
	save([strDataTargetPath 'AggRecUS2.mat'],'sRec');
end
%% plot
dblStep = 0.25;
vecBinE = 0:dblStep:10;
vecBinC = vecBinE(2:end) - dblStep/2;
cellTypes = {'High','Low','Inhibitory','Control'};
vecColA = [0.8 0 0];
vecColW = [0 0 0];
vecColZ = lines(1);

%plot all recordings
figure;maxfig;
for intRec=1:numel(sRec)
	subplot(2,4,intRec)
	vecP = sRec(intRec).matTsZetaP(1,:);
	dblTP = 100*(sum(vecP<0.05)/numel(vecP));
	vecP_Ctrl = sRec(intRec).matTsZetaP(4,:);
	dblFP = 100*(sum(vecP_Ctrl<0.05)/numel(vecP_Ctrl));
	dblPrecision = sum(vecP<0.05) / (sum(vecP<0.05) + sum(vecP_Ctrl<0.05));
	vecCounts = histcounts(-norminv(vecP./2),vecBinE);
	vecCounts_Ctrl = histcounts(-norminv(vecP_Ctrl./2),vecBinE);
	[h,dblDiffP] = ttest(-norminv(vecP./2),-norminv(vecP_Ctrl./2));
	hold on
	%plot(-norminv(0.05/2).*[1 1],[0 max(vecCounts)],'k--');
	plot(vecBinC,vecCounts,'color',lines(1));
	plot(vecBinC,vecCounts_Ctrl,'color',vecColA);
	hold off
	xlabel('Ts-ZETA (z)')
	ylabel('Number of cells (count)');
	title(sprintf('R%d: TP=%.1f%%; FP=%.1f%%; diff-p=%.1e',intRec,dblTP,dblFP,dblDiffP));
	fixfig;grid off;
end
export_fig(fullpath(strFigPath,[strName 'Pvalues.jpg']));
export_fig(fullpath(strFigPath,[strName 'Pvalues.pdf']));

%% plot ROCs
vecAggAucZ = [];
vecAggAucA = [];
vecAggAucW = [];
vecAggTPZ = [];
vecAggTPA = [];
vecAggTPW = [];
vecAggFPZ = [];
vecAggFPA = [];
vecAggFPW = [];
figure;maxfig;
for intRec=1:numel(sRec)
	%zeta
	vecZeta = sort(-norminv(sRec(intRec).matTsZetaP(1,:)./2));
	vecZeta_Ctrl = sort(-norminv(sRec(intRec).matTsZetaP(4,:)./2));
	vecBothData = cat(2,vecZeta,vecZeta_Ctrl);
	vecThresholds = sort(vecBothData,'descend');
	intCells = numel(vecZeta);
	vecTPZ = sum(vecZeta>=vecThresholds',2)/intCells;
	vecFPZ = sum(vecZeta_Ctrl>=vecThresholds',2)/intCells;
	dblThresh = -norminv(0.05/2);
	dblTPZ = sum(vecZeta>dblThresh)/intCells;
	dblFPZ = sum(vecZeta_Ctrl>dblThresh)/intCells;
	%anova
	vecAnova = -norminv(sRec(intRec).matAnovaP(1,:)./2);
	vecAnova_Ctrl = -norminv(sRec(intRec).matAnovaP(4,:)./2);
	vecBothData = cat(2,vecAnova,vecAnova_Ctrl);
	vecThresholds = sort(vecBothData,'descend');
	intCells = numel(vecAnova);
	vecTPA = sum(vecAnova>=vecThresholds',2)/intCells;
	vecFPA = sum(vecAnova_Ctrl>=vecThresholds',2)/intCells;
	dblTPA = sum(vecAnova>dblThresh)/intCells;
	dblFPA = sum(vecAnova_Ctrl>dblThresh)/intCells;
	%ranksum
	vecWilcox = -norminv(sRec(intRec).matWilcoxP(1,:)./2);
	vecWilcox_Ctrl = -norminv(sRec(intRec).matWilcoxP(4,:)./2);
	vecBothData = cat(2,vecWilcox,vecWilcox_Ctrl);
	vecThresholds = sort(vecBothData,'descend');
	intCells = numel(vecWilcox);
	vecTPW = sum(vecWilcox>=vecThresholds',2)/intCells;
	vecFPW = sum(vecWilcox_Ctrl>=vecThresholds',2)/intCells;
	dblTPW = sum(vecWilcox>dblThresh)/intCells;
	dblFPW = sum(vecWilcox_Ctrl>dblThresh)/intCells;
	
	%test
	[AUC_Z,AciZ,Ase_Z,pAuc_Z] = getAuc(vecZeta,vecZeta_Ctrl,0.05,'mann-whitney');
	[AUC_A,AciA,Ase_A,pAuc_A] = getAuc(vecAnova,vecAnova_Ctrl,0.05,'mann-whitney');
	[AUC_W,AciA,Ase_W,pAuc_W] = getAuc(vecWilcox,vecWilcox_Ctrl,0.05,'mann-whitney');
	
	%z vs a
	m0 = AUC_Z - AUC_A;
	s0 = (Ase_Z + Ase_A)/2;
	zZA = abs(m0/s0);
	AUC_pZA = normcdf(zZA,'upper')*2;
	
	%plot
	subplot(2,4,intRec)
	hold on
	plot(vecFPZ,vecTPZ,'color',lines(1));
	plot(vecFPA,vecTPA,'color',vecColA);
	plot(vecFPW,vecTPW,'color',vecColW);
	hold off
	xlabel('False positive (FP) rate');
	ylabel('True positive (TP) rate');
	h=legend({sprintf('Ts-ZETA, AUC=%.3f',AUC_Z),sprintf('ANOVA, AUC=%.3f',AUC_A),sprintf('Wilcox, AUC=%.3f',AUC_W)},'location','best');
	title(sprintf('FP; Z=%.3f,A=%.3f,W=%.3f',dblFPZ,dblFPA,dblFPW))
	if AUC_pZA < 0.01
		h.Title.String=sprintf('AUC Z-A, p=%.2e',AUC_pZA);
	else
		h.Title.String=sprintf('AUC Z-A, p=%.3f',AUC_pZA);
	end
	fixfig;
	
	%save if any is significant
	if (AUC_Z > 0.5 && pAuc_Z < 0.05) || (AUC_A > 0.5 && pAuc_A < 0.05) || (AUC_W > 0.5 && pAuc_W < 0.05)
		vecAggAucZ(end+1) = AUC_Z;
		vecAggAucA(end+1) = AUC_A;
		vecAggAucW(end+1) = AUC_W;
		vecAggTPZ(end+1) = dblTPZ;
		vecAggTPA(end+1) = dblTPA;
		vecAggTPW(end+1) = dblTPW;
		vecAggFPZ(end+1) = dblFPZ;
		vecAggFPA(end+1) = dblFPA;
		vecAggFPW(end+1) = dblFPW;
	end
end
export_fig(fullpath(strFigPath,[strName 'ROCs.jpg']));
export_fig(fullpath(strFigPath,[strName 'ROCs.pdf']));

%% plot aggregate AUC & FP
intSigNum = numel(vecAggAucZ);
figure,maxfig;
subplot(2,3,1)
hold on
errorbar(1,mean(vecAggAucZ),std(vecAggAucZ)./sqrt(intSigNum),'x','color',vecColZ);
errorbar(2,mean(vecAggAucA),std(vecAggAucA)./sqrt(intSigNum),'x','color',vecColA);
errorbar(3,mean(vecAggAucW),std(vecAggAucW)./sqrt(intSigNum),'x','color',vecColW);
hold off
xlim([0.5 3.5]);
ylim([0.5 0.8]);
ylabel('AUC (mean +/- st err)');
set(gca,'xtick',[1 2 3],'xticklabel',{'Ts-ZETA','ANOVA','Wilcox'});
xtickangle(30);
fixfig;
[h,pAZ]=ttest(vecAggAucZ,vecAggAucA);
[h,pWZ]=ttest(vecAggAucZ,vecAggAucW);
[h,pAW]=ttest(vecAggAucW,vecAggAucA);
title(sprintf('T-tests AUCs;AZ-p=%.3f,WZ-p=%.3f,WA-p=%.3f',pAZ,pWZ,pAW));

subplot(2,3,2)
hold on
plot([0.5 3.5],[0.05 0.05],'--','color',[0.3 0.3 0.3]);
text(1.5,0.05,'\alpha = 0.05','Fontsize',14,'horizontalalignment','center','verticalalignment','bottom','color',[0.3 0.3 0.3])
errorbar(1,mean(vecAggFPZ),std(vecAggFPZ)./sqrt(intSigNum),'x','color',vecColZ);
errorbar(2,mean(vecAggFPA),std(vecAggFPA)./sqrt(intSigNum),'x','color',vecColA);
errorbar(3,mean(vecAggFPW),std(vecAggFPW)./sqrt(intSigNum),'x','color',vecColW);
hold off
xlim([0.5 3.5]);
%ylim([0.5 0.8]);
ylabel('FP rate (mean +/- st err)');
set(gca,'xtick',[1 2 3],'xticklabel',{'Ts-ZETA','ANOVA','Wilcox'});
xtickangle(30);
%set(gca,'yscale','log');ylim([0.002 0.5]);
fixfig;
[h,pAlphaZ]=ttest(vecAggFPZ,0.05);
[h,pAlphaA]=ttest(vecAggFPA,0.05);
[h,pAlphaW]=ttest(vecAggFPW,0.05);
title(sprintf('T-tests FPR vs 0.05;Z-p=%.3f,A-p=%.2e,W-p=%.2e',pAlphaZ,pAlphaA,pAlphaW));

export_fig(fullpath(strFigPath,[strName 'Summary.jpg']));
export_fig(fullpath(strFigPath,[strName 'Summary.pdf']));
