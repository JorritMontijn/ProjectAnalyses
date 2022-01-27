clear all;
%close all;
strPath = 'D:\Data\Processed\TraceZeta\';
strFigPath = 'D:\Data\Results\TraceZeta\';

%% prep
intResamps = 250;
strFileSearch = ['TsZetaFeatureConjunctionsResamp' num2str(intResamps) '.mat'];
sDir = dir(fullpath(strPath,strFileSearch));
strFile = sDir(1).name;
sLoad = load(fullpath(sDir(1).folder,strFile));
matTsZeta = sLoad.matTsZeta;
matZeta = sLoad.matZeta;
matTtest = sLoad.matTtest;
dblTau0 = 63/1000;%s
vecSampFreqs = sLoad.vecSampFreqs;
vecTau = sLoad.vecTau*dblTau0;
vecTrialDur = sLoad.vecTrialDur;
strRec = 'FeatureConjunctions';
intTauNum = numel(vecTau);

%% plot
vecZetaAuc = nan(1,numel(vecTrialDur));
vecTtestAuc = nan(1,numel(vecTrialDur));
matZetaReal = squeeze(matZeta(:,1,:));
matZetaRand = squeeze(matZeta(:,2,:));
matTtestReal = squeeze(matTtest(:,1,:));
matTtestRand = squeeze(matTtest(:,2,:));
figure
for intDurIdx=1:numel(vecTrialDur)
	vecTP = matZetaReal(:,intDurIdx);
	vecFP = matZetaRand(:,intDurIdx);
	vecZetaAuc(intDurIdx) = getAuc(vecTP,vecFP);
	
	vecTTP = matTtestReal(:,intDurIdx);
	vecTFP = matTtestRand(:,intDurIdx);
	vecTtestAuc(intDurIdx) = getAuc(vecTTP,vecTFP);
end
matTsAuc3 = nan(numel(vecTrialDur),numel(vecTau),numel(vecSampFreqs));
for intSampFreqIdx=1:numel(vecSampFreqs)
	dblSampFreq = vecSampFreqs(intSampFreqIdx);
	%matTsZeta(intNeuron,intRunType,intTrialDurIdx,intTauIdx,intSampFreqIdx)
	matTsZetaReal = squeeze(matTsZeta(:,1,:,:,intSampFreqIdx)); %neuron x trial dur x tau
	matTsZetaRand = squeeze(matTsZeta(:,2,:,:,intSampFreqIdx));
	
	for intDurIdx=1:numel(vecTrialDur)
		for intTauIdx=1:numel(vecTau)
			vecTP = matTsZetaReal(:,intDurIdx,intTauIdx);
			vecFP = matTsZetaRand(:,intDurIdx,intTauIdx);
			
			matTsAuc3(intDurIdx,intTauIdx,intSampFreqIdx) = getAuc(vecTP,vecFP);
		end
	end
	subplot(2,4,intSampFreqIdx)
	imagesc(matTsAuc3(:,:,intSampFreqIdx),[0.6 0.95]);
	title(sprintf('Sampling frequency %.1f Hz',dblSampFreq))
	h=colorbar;
	title(h,'AUC');
	set(gca,'xtick',1:numel(vecTau),'xticklabel',vecTau,'ytick',1:numel(vecTrialDur),'yticklabel',roundi(vecTrialDur,2));
	axis xy
	xtickangle(45);
	xlabel('Ca Indicator lifetime tau (s)');
	ylabel('Trial duration (s)');
end
