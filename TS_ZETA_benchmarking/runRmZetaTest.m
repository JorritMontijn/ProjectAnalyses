vecT = (1:3:20)./20;
vecT = vecT*pi;
intN1 = 10;
intN2 = 15;
vecSNR = 1./exp(3:-0.5:0);
intReps = 100;
matZetaP = nan(numel(vecSNR),intReps,2);
matRmAnovaP = nan(numel(vecSNR),intReps,2);
hTic=tic;
for intSNR=1:numel(vecSNR)
	dblSNR = vecSNR(intSNR);
	for intRep=1:intReps
		if toc(hTic)>5
			fprintf('SNR %d/%d, rep %d/%d [%s]\n',intSNR,numel(vecSNR),intRep,intReps,getTime);
			hTic=tic;
		end
		for intRand=[2 1]
			matCond1 = (1/dblSNR)*rand(intN1,numel(vecT)) + (1*rand(intN1,1))*sin(vecT);
			
			if intRand == 1
				matCond2 = (1/dblSNR)*rand(intN2,numel(vecT)) + (1*rand(intN2,1))*cos(vecT);
			else
				matCond2 = (1/dblSNR)*rand(intN2,numel(vecT)) + (1*rand(intN2,1))*sin(vecT);
			end
			intPlot = 0;
			intResampNum = 250;
			dblZetaP = zetarmtest2(vecT,matCond1,vecT,matCond2,intResampNum,intPlot);
			
			t=table;
			t.Cond = cat(1,false(size(matCond1,1),1),true(size(matCond2,1),1));
			for intT=1:numel(vecT)
				t.(['t' num2str(intT)]) = cat(1,matCond1(:,intT),matCond2(:,intT));
			end
			rm = fitrm(t,'t1-t7 ~ Cond','WithinDesign',vecT','WithinModel','orthogonalcontrasts');
			ranovatbl = ranova(rm);
			dblRmAnova = ranovatbl.pValue(2);
			
			matZetaP(intSNR,intRep,intRand) = dblZetaP;
			matRmAnovaP(intSNR,intRep,intRand) = dblRmAnova;
		end
	end
end

%% plot
dblAlpha = 0.05;
[vecTPRZ,matTPRZ_ci] = binofit(sum(matZetaP(:,:,1)<dblAlpha,2),intReps);
[vecFPRZ,matFPRZ_ci] = binofit(sum(matZetaP(:,:,2)<dblAlpha,2),intReps);

[vecTPRA,matTPRA_ci] = binofit(sum(matRmAnovaP(:,:,1)<dblAlpha,2),intReps);
[vecFPRA,matFPRA_ci] = binofit(sum(matRmAnovaP(:,:,2)<dblAlpha,2),intReps);
figure
subplot(2,3,1)
hold on
errorbar(vecT,mean(matCond1,1),std(matCond1,[],1),'xk-');
errorbar(vecT,mean(matCond2,1),std(matCond2,[],1),'xb-');
title('example data');
legend({'condition1','condition2'},'location','best');
xlabel('time');
ylabel('signal');

subplot(2,3,2)
set(gca,'xscale','log');
hold on
errorbar(vecSNR,vecTPRZ,matTPRZ_ci(:,1)-vecTPRZ,matTPRZ_ci(:,2)-vecTPRZ,'color',lines(1));
errorbar(vecSNR,vecTPRA,matTPRA_ci(:,1)-vecTPRA,matTPRA_ci(:,2)-vecTPRA,'color',[0.8 0 0]);
xlabel('SNR')
ylabel('True positive rate');
title(sprintf('TPR, N=%d reps, alpha=%.3f',intReps,dblAlpha));
legend({'R-ZETA2','R-ANOVA'},'location','best');

subplot(2,3,3)
set(gca,'xscale','log');
hold on
errorbar(vecSNR,vecFPRZ,matFPRZ_ci(:,1)-vecFPRZ,matFPRZ_ci(:,2)-vecFPRZ,'color',lines(1));
errorbar(vecSNR,vecFPRA,matFPRA_ci(:,1)-vecFPRA,matFPRA_ci(:,2)-vecFPRA,'color',[0.8 0 0]);
xlabel('SNR')
ylabel('False positive rate');
title(sprintf('FPR, N=%d reps, alpha=%.3f',intReps,dblAlpha));
legend({'R-ZETA2','R-ANOVA'},'location','best');
maxfig;fixfig;

