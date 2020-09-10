Contrastmin1 = 80;
Contrastmax1 = 100;
Contrastmin2 = 10;
Contrastmax2 = 170;

cont.c0_25 = GetContrastMean(Contrastmin1, Contrastmax1, orientall, plaidsactall, condition.condition1);
cont.c0_50 = GetContrastMean(Contrastmin1, Contrastmax1, orientall, plaidsactall, condition.condition2);
cont.c0_100 = GetContrastMean(Contrastmin1, Contrastmax1, orientall, orientactall, condition.condition10);
cont.c25_0 = GetContrastMean2(Contrastmin2, Contrastmax2, orientall, plaidsactall, condition.condition6);
cont.c50_0 = GetContrastMean2(Contrastmin2, Contrastmax2, orientall, plaidsactall, condition.condition7);
cont.c100_0 = GetContrastMean2(Contrastmin2, Contrastmax2, orientall, orientactall, condition.condition9);

contac.mat.c0_25 = vertcat(cont.c0_25(1:8).mat);
contac.mat.c0_50 = vertcat(cont.c0_50(1:8).mat);
contac.mat.c0_100 = vertcat(cont.c0_100(1:8).mat);
contac.mat.c25_0 = vertcat(cont.c25_0(1:8).mat);
contac.mat.c50_0 = vertcat(cont.c50_0(1:8).mat);
contac.mat.c100_0 = vertcat(cont.c100_0(1:8).mat);

contac.neur.c0_25 = mean(contac.mat.c0_25, 2);
contac.neur.c0_50 = mean(contac.mat.c0_50, 2);
contac.neur.c0_100 = mean(contac.mat.c0_100, 2);
contac.neur.c25_0 = mean(contac.mat.c25_0, 2);
contac.neur.c50_0 = mean(contac.mat.c50_0, 2);
contac.neur.c100_0 = mean(contac.mat.c100_0, 2);

contac.neur.c25 = vertcat(contac.neur.c0_25, contac.neur.c25_0);
contac.neur.c50 = vertcat(contac.neur.c0_50, contac.neur.c50_0);
contac.neur.c100 = vertcat(contac.neur.c0_100, contac.neur.c100_0);

%% contrast 0-0 code
for i = 1:8
    datco(ContrastVect(i)).ses.neuron = dat(ContrastVect(i)).ses.neuron(orientall(i).significantindex);
    datco(OrientVect(i)).ses.neuron = dat(OrientVect(i)).ses.neuron(orientall(i).significantindex);
end

activity90.peak = max(vertcat(plaidneurons(10).tuning.activity));
activity0.peak = max(vertcat(plaidneurons(9).tuning.activity));

for j = 1:8
    for i = 1:length(datco(ContrastVect(j)).ses.neuron)
        for k = 1:length(plaid(j).contrindex)
            for h = 1:length(dat(OrientVect(j)).ses.structStim.FrameOn)
            baseline.means(ContrastVect(j)).neuron(i).timeframe(k) = mean(datco(ContrastVect(j)).ses.neuron(i).dFoF((plaid(j).contrindex(k).FrameOn-25):plaid(j).contrindex(k).FrameOn));
            baseline.means(OrientVect(j)).neuron(i).timeframe(h) = mean(datco(OrientVect(j)).ses.neuron(i).dFoF((dat(OrientVect(j)).ses.structStim.FrameOn(h)-25):dat(OrientVect(j)).ses.structStim.FrameOn(h)));
            end
        end
    end
end

for j = 1:8
    for i = 1:length(datco(ContrastVect(j)).ses.neuron)
        baseline.combined(j).neuron(i).timeframe = horzcat(baseline.means(ContrastVect(j)).neuron(i).timeframe, baseline.means(OrientVect(j)).neuron(i).timeframe);
    end
end

for j = 1:8
    for i = 1:length(baseline.combined(j).neuron)
        baseline.combined(j).neuron(i).gm = mean(baseline.combined(j).neuron(i).timeframe);
        baselinegms(j).neuron = vertcat(baseline.combined(j).neuron.gm);
    end
end

CO.means = transpose(vertcat(baselinegms(1:8).neuron));
CO.sd = std(vertcat(baselinegms(1:8).neuron));
CO.se = CO.sd/sqrt(length(vertcat(baselinegms(1:8).neuron)));
CO.gm = mean(CO.means)
%%

contac.mean.c25 = mean(contac.neur.c25);
contac.mean.c50 = mean(contac.neur.c50);
contac.mean.c100 = mean(contac.neur.c100);
contac.sd.c25 = std(contac.neur.c25);
contac.sd.c50 = std(contac.neur.c50);
contac.sd.c100 = std(contac.neur.c100);
contac.se.c25 = contac.sd.c25/sqrt(length(contac.neur.c25));
contac.se.c50 = contac.sd.c50/sqrt(length(contac.neur.c50));
contac.se.c100 = contac.sd.c100/sqrt(length(contac.neur.c100));

contacplotdata.contrast_ttest_bar = bar(1:4, [CO.gm, contac.mean.c25, contac.mean.c50, contac.mean.c100]);
hold on
errorbar(1:4,[CO.gm, contac.mean.c25, contac.mean.c50, contac.mean.c100],[CO.se, contac.se.c25, contac.se.c50, contac.se.c100],'.', 'Color', [0 0 0]);
grid on
box off
contacplotdata.contrast_ttest_bar.FaceColor = [1 1 1];
xlabels = ["Contrast 0"; "Contrast 25"; "Contrast 50"; "Contrast 100"];
set(gca, 'xticklabels', xlabels)
ylabel ('Mean neuron activity (\DeltaF/F)');
hold off
print(gcf, '-bestfit','ContrastTuning','-dpdf')

[H, P, CI, STATS] = ttest2(CO.means, contac.neur.c25)
[H, P, CI, STATS] = ttest2(contac.neur.c25, contac.neur.c50)
[H, P, CI, STATS] = ttest2(contac.neur.c50, contac.neur.c100)
%%

[p,tbl,stats] = anova1([CO.means, contac.neur.c25, contac.neur.c50, contac.neur.c100]);