for j = 1:8
    meanactivity(j).a0_25 = mean(mean(plaidsactall(j).activ(:,conditions(j).condition1)));
    meanactivity(j).a0_50 = mean(mean(plaidsactall(j).activ(:,conditions(j).condition2)));
    meanactivity(j).a25_25 = mean(mean(plaidsactall(j).activ(:,conditions(j).condition3)));
    meanactivity(j).a25_50 = mean(mean(plaidsactall(j).activ(:,conditions(j).condition4)));
    meanactivity(j).a50_50 = mean(mean(plaidsactall(j).activ(:,conditions(j).condition5)));
    meanactivity(j).a25_0 = mean(mean(plaidsactall(j).activ(:,conditions(j).condition6))); 
    meanactivity(j).a50_0 = mean(mean(plaidsactall(j).activ(:,conditions(j).condition7)));
    meanactivity(j).a50_25 = mean(mean(plaidsactall(j).activ(:,conditions(j).condition8)));
    meanactivity(j).a100_0 = mean(mean(orientactall(j).activ(:,conditions(j).condition9)));
    meanactivity(j).a0_100 = mean(mean(orientactall(j).activ(:,conditions(j).condition10)));
end

for j = [2, 3, 4, 7, 8]
    PVmeanactivity(j).a0_25 = mean(mean(plaidsPVactall(j).activ(:,conditions(j).condition1)));
    PVmeanactivity(j).a0_50 = mean(mean(plaidsPVactall(j).activ(:,conditions(j).condition2)));
    PVmeanactivity(j).a25_25 = mean(mean(plaidsPVactall(j).activ(:,conditions(j).condition3)));
    PVmeanactivity(j).a25_50 = mean(mean(plaidsPVactall(j).activ(:,conditions(j).condition4)));
    PVmeanactivity(j).a50_50 = mean(mean(plaidsPVactall(j).activ(:,conditions(j).condition5)));
    PVmeanactivity(j).a25_0 = mean(mean(plaidsPVactall(j).activ(:,conditions(j).condition6))); 
    PVmeanactivity(j).a50_0 = mean(mean(plaidsPVactall(j).activ(:,conditions(j).condition7)));
    PVmeanactivity(j).a50_25 = mean(mean(plaidsPVactall(j).activ(:,conditions(j).condition8)));
    PVmeanactivity(j).a100_0 = mean(mean(orientpvactall(j).activ(:,conditions(j).condition9)));
    PVmeanactivity(j).a0_100 = mean(mean(orientpvactall(j).activ(:,conditions(j).condition10)));
end
PVgrandmeanactivity.a0_25 = mean(vertcat(PVmeanactivity(1:8).a0_25));
PVgrandmeanactivity.a0_50 = mean(vertcat(PVmeanactivity(1:8).a0_50));
PVgrandmeanactivity.a25_25 = mean(vertcat(PVmeanactivity(1:8).a25_25));
PVgrandmeanactivity.a25_50 = mean(vertcat(PVmeanactivity(1:8).a25_50));
PVgrandmeanactivity.a50_50 = mean(vertcat(PVmeanactivity(1:8).a50_50));
PVgrandmeanactivity.a25_0 = mean(vertcat(PVmeanactivity(1:8).a25_0));
PVgrandmeanactivity.a50_0 = mean(vertcat(PVmeanactivity(1:8).a50_0));
PVgrandmeanactivity.a50_25 = mean(vertcat(PVmeanactivity(1:8).a50_25));
PVgrandmeanactivity.a100_0 = mean(vertcat(PVmeanactivity(1:8).a100_0));
PVgrandmeanactivity.a0_100 = mean(vertcat(PVmeanactivity(1:8).a0_100));

corrmatrix(1:146, 1) = COIcorrelations2.datafull.vals.val_0_25;
corrmatrix(1:146, 2) = COIcorrelations2.datafull.vals.val_0_50;
corrmatrix(1:146, 3) = COIcorrelations2.datafull.vals.val_0_100;
corrmatrix(1:146, 4) = COIcorrelations2.datafull.vals.val_25_0;
corrmatrix(1:146, 5) = COIcorrelations2.datafull.vals.val_25_25;
corrmatrix(1:146, 6) = COIcorrelations2.datafull.vals.val_25_50;
corrmatrix(1:146, 7) = COIcorrelations2.datafull.vals.val_50_0;
corrmatrix(1:146, 8) = COIcorrelations2.datafull.vals.val_50_25;
corrmatrix(1:146, 9) = COIcorrelations2.datafull.vals.val_50_50;
corrmatrix(1:146, 10) = COIcorrelations2.datafull.vals.val_100_0;

corrmatrix2(1:146, 1) = COIcorrelations2.datafull.vals.val_25_25;
corrmatrix2(1:146, 2) = COIcorrelations2.datafull.vals.val_50_50;
corrmatrix2(1:146, 3) = COIcorrelations2.datafull.vals.val_25_50;
corrmatrix2(147:292, 3) = COIcorrelations2.datafull.vals.val_50_25;

corrmatrix2(corrmatrix2==0) = NaN;
%%
figure (1)
[p,tbl,stats] = kruskalwallis(corrmatrix)
hold on
xlabs  =  {'Contrasts 0-25', 'Contrasts 0-50', 'Contrasts 0-100', 'Contrasts 25-0', 'Contrasts 25-25', 'Contrasts 25-50', 'Contrasts 50-0', 'Contrasts 50-25', 'Contrasts 50-50', 'Contrasts 100-0'};
xticklabels (xlabs)

box off
grid on
hold off
print(gcf, '-bestfit','corrs_bplot_all','-dpdf')
[p,h,stats] = signrank(corrmatrix2(:,1), corrmatrix2(:,3))
[p,h,stats] = signrank(corrmatrix2(:,2), corrmatrix2(:,3))
[p,h,stats] = signrank(corrmatrix2(:,1), corrmatrix2(:,2))
%%
figure (2)
[p,tbl,stats] = kruskalwallis(corrmatrix2)
hold on
xlabs  =  {'Contrasts 25-25', 'Contrasts 50-50', 'Contrasts 25-50'};
xticklabels (xlabs)

box off
grid on
hold off
print(gcf, '-bestfit','corrs_bplot','-dpdf')
[p,h,stats] = signrank(corrmatrix2(:,1), corrmatrix2(:,3))
[p,h,stats] = signrank(corrmatrix2(:,2), corrmatrix2(:,3))
[p,h,stats] = signrank(corrmatrix2(:,1), corrmatrix2(:,2))