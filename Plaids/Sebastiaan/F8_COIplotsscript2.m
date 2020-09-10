COIcorrelations2.data(1) = COIfunc2(dat, plaidsPVactall, plaidsact, plaidsactall, orientactall, orientpvactall, orientall, conditions, 2);
%COIcorrelations.data(2) = COIfunc2(dat, plaidsPVactall, plaidsact, plaidsactall, orientactall, orientpvactall, orientall, conditions, 3);
COIcorrelations2.data(2) = COIfunc2(dat, plaidsPVactall, plaidsact, plaidsactall, orientactall, orientpvactall, orientall, conditions, 4);
COIcorrelations2.data(3) = COIfunc2(dat, plaidsPVactall, plaidsact, plaidsactall, orientactall, orientpvactall, orientall, conditions, 7);
COIcorrelations2.data(4) = COIfunc2(dat, plaidsPVactall, plaidsact, plaidsactall, orientactall, orientpvactall, orientall, conditions, 8);
for i = 1:length(COIcorrelations2.data)
    corrvec2_25_25(i).vec = reshape(COIcorrelations2.data(i).vals.val_25_25, [1, numel(COIcorrelations2.data(i).vals.val_25_25)]);
    corrvec2_25_50(i).vec = reshape(COIcorrelations2.data(i).vals.val_25_50, [1, numel(COIcorrelations2.data(i).vals.val_25_50)]);
    corrvec2_50_25(i).vec = reshape(COIcorrelations2.data(i).vals.val_50_25, [1, numel(COIcorrelations2.data(i).vals.val_50_25)]);
    corrvec2_50_50(i).vec = reshape(COIcorrelations2.data(i).vals.val_50_50, [1, numel(COIcorrelations2.data(i).vals.val_50_50)]);
    corrvec2_0_25(i).vec = reshape(COIcorrelations2.data(i).vals.val_0_25, [1, numel(COIcorrelations2.data(i).vals.val_0_25)]);
    corrvec2_25_0(i).vec = reshape(COIcorrelations2.data(i).vals.val_25_0, [1, numel(COIcorrelations2.data(i).vals.val_25_0)]);
    corrvec2_0_50(i).vec = reshape(COIcorrelations2.data(i).vals.val_0_50, [1, numel(COIcorrelations2.data(i).vals.val_0_50)]);
    corrvec2_50_0(i).vec = reshape(COIcorrelations2.data(i).vals.val_50_0, [1, numel(COIcorrelations2.data(i).vals.val_50_0)]);
    corrvec2_0_100(i).vec = reshape(COIcorrelations2.data(i).vals.val_0_100, [1, numel(COIcorrelations2.data(i).vals.val_0_100)]);
    corrvec2_100_0(i).vec = reshape(COIcorrelations2.data(i).vals.val_100_0, [1, numel(COIcorrelations2.data(i).vals.val_100_0)]);
end

COIcorrelations2.datafull.vals.val_25_25 = horzcat(corrvec2_25_25(1:4).vec);
COIcorrelations2.datafull.vals.val_25_50 = horzcat(corrvec2_25_50(1:4).vec);
COIcorrelations2.datafull.vals.val_50_25 = horzcat(corrvec2_50_25(1:4).vec);
COIcorrelations2.datafull.vals.val_50_50 = horzcat(corrvec2_50_50(1:4).vec);
COIcorrelations2.datafull.vals.val_0_25 = horzcat(corrvec2_0_25(1:4).vec);
COIcorrelations2.datafull.vals.val_25_0 = horzcat(corrvec2_25_0(1:4).vec);
COIcorrelations2.datafull.vals.val_0_50 = horzcat(corrvec2_0_50(1:4).vec);
COIcorrelations2.datafull.vals.val_50_0 = horzcat(corrvec2_50_0(1:4).vec);
COIcorrelations2.datafull.vals.val_0_100 = horzcat(corrvec2_0_100(1:4).vec);
COIcorrelations2.datafull.vals.val_100_0 = horzcat(corrvec2_100_0(1:4).vec);


for i = 1:8
    datco(ContrastVect(i)).ses.neuron = dat(ContrastVect(i)).ses.neuron(orientall(i).significantindex);
    datco(OrientVect(i)).ses.neuron = dat(OrientVect(i)).ses.neuron(orientall(i).significantindex);
end

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
for j = 1:5
    for i = 1:length(dat(ContrastPVVect(j)).ses.PV)
        for k = 1:length(plaid(j).contrindex)
            for h = 1:length(dat(OrientPVVect(j)).ses.structStim.FrameOn)
                baseline.means(ContrastPVVect(j)).PV(i).timeframe(k) = mean(dat(ContrastPVVect(j)).ses.PV(i).dFoF((plaid(find(ismember(ContrastVect, ContrastPVVect(j)))).contrindex(k).FrameOn-25):plaid(find(ismember(ContrastVect, ContrastPVVect(j)))).contrindex(k).FrameOn));
                baseline.means(OrientPVVect(j)).PV(i).timeframe(h) = mean(dat(OrientPVVect(j)).ses.PV(i).dFoF((dat(OrientVect(j)).ses.structStim.FrameOn(h)-25):dat(OrientVect(j)).ses.structStim.FrameOn(h)));
            end
        end
    end
end

ActiveSets = [3,4,7,8,13,14,15,16];
for j = 1:8
    activity.baseline(j).neuron = baseline.means(ActiveSets(j)).neuron;
    activity.baseline(j).PV = baseline.means(ActiveSets(j)).PV;
end

for j = 1:8
    for i = 1:length(activity.baseline(j).neuron)
        for k = 1:length(activity.baseline(j).neuron(i).timeframe)
            activitymat(j).neuron(i, k) = activity.baseline(j).neuron(i).timeframe(k);
        end
    end
end

for j = 1:8
    for i = 1:length(activity.baseline(j).PV)
        for k = 1:length(activity.baseline(j).PV(i).timeframe)
            activitymat(j).PV(i, k) = activity.baseline(j).PV(i).timeframe(k);
            activitymat(j).PVmean(k) = mean(activitymat(j).PV(:,k));
        end
    end
end

for j = 1:8
    for i = 1:length(activity.baseline(j).neuron)
             ActivityCorr2(j).neuron(i) = corr(activitymat(j).neuron(i,:)', activitymat(j).PVmean');
             ActivityCorrvec2(j).dataset = reshape(ActivityCorr2(j).neuron, [1, numel(ActivityCorr2(j).neuron)]);
    end
end
ActivityCorrAll2 = horzcat(ActivityCorrvec2(1:8).dataset);

%% Subplots method 2
figure(2)

subplot(4,4,1)
histogram(ActivityCorrAll2, 'FaceColor', 'none', 'NumBins', 15)
hold on
xlim ([-1, 1])
box off
title (" ")
ylabel (" ")
xlabel (" ")

subplot(4,4,2)
histogram(COIcorrelations2.datafull.vals.val_0_25, 'FaceColor', 'none', 'NumBins', 15)
xlim ([-1, 1])
box off
title (" ")
ylabel (" ")
xlabel (" ")

subplot(4,4,3)
histogram(COIcorrelations2.datafull.vals.val_0_50, 'FaceColor', 'none', 'NumBins', 15)
xlim ([-1, 1])
box off
title (" ")
ylabel(" ")
xlabel(" ")

subplot(4,4,4)
histogram(COIcorrelations2.datafull.vals.val_0_100, 'Facecolor', 'none', 'NumBins', 15)
xlim ([-1, 1])
box off
title (" ")
ylabel (" ")
xlabel(" ")

subplot(4,4,5)
histogram(COIcorrelations2.datafull.vals.val_25_0, 'FaceColor', 'none', 'NumBins', 15)
xlim ([-1, 1])
box off
title (" ")
ylabel (" ")
xlabel (" ")

subplot(4,4,6)
histogram(COIcorrelations2.datafull.vals.val_25_25, 'FaceColor', 'none', 'NumBins', 15)
xlim ([-1, 1])
box off
title (" ")
ylabel (" ")
xlabel (" ")

subplot(4,4,7)
histogram(COIcorrelations2.datafull.vals.val_25_50, 'FaceColor', 'none', 'NumBins', 15)
xlim ([-1, 1])
box off
title (" ")
ylabel (" ")
xlabel (" ")

subplot(4,4,9)
histogram(COIcorrelations2.datafull.vals.val_50_0, 'FaceColor', 'none', 'NumBins', 15)
xlim ([-1, 1])
box off
title (" ")
ylabel (" ")
xlabel(" ")

subplot(4,4,10)
histogram(COIcorrelations2.datafull.vals.val_50_25, 'FaceColor', 'none', 'NumBins', 15)
xlim ([-1, 1])
box off
title (" ")
ylabel (" ")
xlabel (" ")

subplot(4,4,11)
histogram(COIcorrelations2.datafull.vals.val_50_50, 'FaceColor', 'none', 'NumBins', 15)
xlim ([-1, 1])
box off
title (" ")
ylabel (" ")
xlabel (" ")

subplot(4,4,13)
histogram(COIcorrelations2.datafull.vals.val_100_0, 'FaceColor', 'none', 'NumBins', 15)
xlim ([-1, 1])
box off
title (" ")
ylabel (" ")
xlabel (" ")

normaxes(gcf, 'y', 1:11)
hold off

print(gcf, '-bestfit','COI_singlePV','-dpdf')
