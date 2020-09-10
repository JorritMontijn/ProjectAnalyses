%%Preallocation
forplots = struct('Zero', cell(1,8), 'Twenty', cell(1,8), 'Forty', cell(1,8), 'Sixty', cell(1,8), 'Ninety', cell(1,8), 'HundredTen', cell(1,8), 'HundredThirty', cell(1,8), 'HundredFifty', cell(1,8)); 
%%finding the scores corresponding to the different stimulus orientations
%%to index for plots.  
for j = 1:8
    for i = 1:length(orientall(j).significant)
        forplots(j).Zero(i)             =                  mean(orientactall(j).activZG(orientall(j).significantindex(i),dat(OrientVect(j)).ses.structStim.Orientation == 0));
        forplots(j).Twenty(i)           =                  mean(orientactall(j).activZG(orientall(j).significantindex(i),dat(OrientVect(j)).ses.structStim.Orientation == 22.5));
        forplots(j).Forty(i)            =                  mean(orientactall(j).activZG(orientall(j).significantindex(i),dat(OrientVect(j)).ses.structStim.Orientation == 45));
        forplots(j).Sixty(i)            =                  mean(orientactall(j).activZG(orientall(j).significantindex(i),dat(OrientVect(j)).ses.structStim.Orientation == 67.5));
        forplots(j).Ninety(i)           =                  mean(orientactall(j).activZG(orientall(j).significantindex(i),dat(OrientVect(j)).ses.structStim.Orientation == 90));
        forplots(j).HundredTen(i)       =                  mean(orientactall(j).activZG(orientall(j).significantindex(i),dat(OrientVect(j)).ses.structStim.Orientation == 112.5));
        forplots(j).HundredThirty(i)    =                  mean(orientactall(j).activZG(orientall(j).significantindex(i),dat(OrientVect(j)).ses.structStim.Orientation == 135));
        forplots(j).HundredFifty(i)     =                  mean(orientactall(j).activZG(orientall(j).significantindex(i),dat(OrientVect(j)).ses.structStim.Orientation == 157.5));
    end
end


subplot(2,4,1);
degrees_0 = scatter(orientall(1).significant, forplots(1).Zero,10,'black', "filled");
hold on
%scatter(orientall(2).significant, forplots(2).Zero,10,'yellow', "filled"); niet een significant aantal significante tuning curves --> eruit laten
%scatter(orientall(3).significant, forplots(3).Zero,10,'black', "filled"); niet een significant aantal significante tuning curves --> eruit laten
scatter(orientall(4).significant, forplots(4).Zero,10,'black', "filled");
scatter(orientall(5).significant, forplots(5).Zero,10,'black', "filled");
scatter(orientall(6).significant, forplots(6).Zero,10,'black', "filled");
scatter(orientall(7).significant, forplots(7).Zero,10,'black', "filled");
scatter(orientall(8).significant, forplots(8).Zero,10,'black', "filled");
title("0° Stimulus")
ylabel("Activity (Z-scored)")
xlabel("Preferred orientation")
%legend ("dataset 1", "dataset 4", "dataset 5", "dataset 6", "dataset 7", "dataset 8")
hold off

subplot(2,4,2);
degrees_22_5 = scatter(orientall(1).significant, forplots(1).Twenty,10,'black', "filled");
hold on
%scatter(orientall(2).significant, forplots(2).Twenty,10,'yellow', "filled"); niet een significant aantal significante tuning curves --> eruit laten
%scatter(orientall(3).significant, forplots(3).Twenty,10,'black', "filled"); niet een significant aantal significante tuning curves --> eruit laten
scatter(orientall(4).significant, forplots(4).Twenty,10,'black', "filled");
scatter(orientall(5).significant, forplots(5).Twenty,10,'black', "filled");
scatter(orientall(6).significant, forplots(6).Twenty,10,'black', "filled");
scatter(orientall(7).significant, forplots(7).Twenty,10,'black', "filled");
scatter(orientall(8).significant, forplots(8).Twenty,10,'black', "filled");
title("22.5° Stimulus")
ylabel("Activity (Z-scored)")
xlabel("Preferred orientation")
%legend ("dataset 1", "dataset 4", "dataset 5", "dataset 6", "dataset 7", "dataset 8")
hold off

subplot(2,4,3);
degrees_45 = scatter(orientall(1).significant, forplots(1).Forty,10,'black', "filled");
hold on
%scatter(orientall(2).significant, forplots(2).Forty,10,'yellow', "filled"); niet een significant aantal significante tuning curves --> eruit laten
%scatter(orientall(3).significant, forplots(3).Forty,10,'black', "filled"); niet een significant aantal significante tuning curves --> eruit laten
scatter(orientall(4).significant, forplots(4).Forty,10,'black', "filled");
scatter(orientall(5).significant, forplots(5).Forty,10,'black', "filled");
scatter(orientall(6).significant, forplots(6).Forty,10,'black', "filled");
scatter(orientall(7).significant, forplots(7).Forty,10,'black', "filled");
scatter(orientall(8).significant, forplots(8).Forty,10,'black', "filled");
title("45° Stimulus")
ylabel("Activity (Z-scored)")
xlabel("Preferred orientation")
%legend ("dataset 1", "dataset 4", "dataset 5", "dataset 6", "dataset 7", "dataset 8")
hold off

subplot(2,4,4);
degrees_67_5 = scatter(orientall(1).significant, forplots(1).Sixty,10,'black', "filled");
hold on
%scatter(orientall(2).significant, forplots(2).Sixty,10,'yellow', "filled"); niet een significant aantal significante tuning curves --> eruit laten
%scatter(orientall(3).significant, forplots(3).Sixty,10,'black', "filled"); niet een significant aantal significante tuning curves --> eruit laten
scatter(orientall(4).significant, forplots(4).Sixty,10,'black', "filled");
scatter(orientall(5).significant, forplots(5).Sixty,10,'black', "filled");
scatter(orientall(6).significant, forplots(6).Sixty,10,'black', "filled");
scatter(orientall(7).significant, forplots(7).Sixty,10,'black', "filled");
scatter(orientall(8).significant, forplots(8).Sixty,10,'black', "filled");
title("67.5° Stimulus")
ylabel("Activity (Z-scored)")
xlabel("Preferred orientation")
%legend ("dataset 1", "dataset 4", "dataset 5", "dataset 6", "dataset 7", "dataset 8")
hold off

subplot(2,4,5);
degrees_90 = scatter(orientall(1).significant, forplots(1).Ninety,10,'black', "filled");
hold on
%scatter(orientall(2).significant, forplots(2).Ninety,10,'yellow', "filled"); niet een significant aantal significante tuning curves --> eruit laten
%scatter(orientall(3).significant, forplots(3).Ninety,10,'black', "filled"); niet een significant aantal significante tuning curves --> eruit laten
scatter(orientall(4).significant, forplots(4).Ninety,10,'black', "filled");
scatter(orientall(5).significant, forplots(5).Ninety,10,'black', "filled");
scatter(orientall(6).significant, forplots(6).Ninety,10,'black', "filled");
scatter(orientall(7).significant, forplots(7).Ninety,10,'black', "filled");
scatter(orientall(8).significant, forplots(8).Ninety,10,'black', "filled");
title("90° Stimulus")
ylabel("Activity (Z-scored)")
xlabel("Preferred orientation")
%legend ("dataset 1", "dataset 4", "dataset 5", "dataset 6", "dataset 7", "dataset 8")
hold off

subplot(2,4,6);
degrees_112_5 = scatter(orientall(1).significant, forplots(1).HundredTen,10,'black', "filled");
hold on
%scatter(orientall(2).significant, forplots(2).HundredTen,10,'yellow', "filled"); niet een significant aantal significante tuning curves --> eruit laten
%scatter(orientall(3).significant, forplots(3).HundredTen,10,'black', "filled"); niet een significant aantal significante tuning curves --> eruit laten
scatter(orientall(4).significant, forplots(4).HundredTen,10,'black', "filled");
scatter(orientall(5).significant, forplots(5).HundredTen,10,'black', "filled");
scatter(orientall(6).significant, forplots(6).HundredTen,10,'black', "filled");
scatter(orientall(7).significant, forplots(7).HundredTen,10,'black', "filled");
scatter(orientall(8).significant, forplots(8).HundredTen,10,'black', "filled");
title("112.5° Stimulus")
ylabel("Activity (Z-scored)")
xlabel("Preferred orientation")
%legend ("dataset 1", "dataset 4", "dataset 5", "dataset 6", "dataset 7", "dataset 8")
hold off

subplot(2,4,7);
degrees_135 = scatter(orientall(1).significant, forplots(1).HundredThirty,10,'black', "filled");
hold on
%scatter(orientall(2).significant, forplots(2).HundredThirty,10,'yellow', "filled"); niet een significant aantal significante tuning curves --> eruit laten
%scatter(orientall(3).significant, forplots(3).HundredThirty,10,'black', "filled"); niet een significant aantal significante tuning curves --> eruit laten
scatter(orientall(4).significant, forplots(4).HundredThirty,10,'black', "filled");
scatter(orientall(5).significant, forplots(5).HundredThirty,10,'black', "filled");
scatter(orientall(6).significant, forplots(6).HundredThirty,10,'black', "filled");
scatter(orientall(7).significant, forplots(7).HundredThirty,10,'black', "filled");
scatter(orientall(8).significant, forplots(8).HundredThirty,10,'black', "filled");
title("135° Stimulus")
ylabel("Activity (Z-scored)")
xlabel("Preferred orientation")
%legend ("dataset 1", "dataset 4", "dataset 5", "dataset 6", "dataset 7", "dataset 8")
hold off

subplot(2,4,8)
degrees_157_5 = scatter(orientall(1).significant, forplots(1).HundredFifty,10,'black', "filled");
hold on
%scatter(orientall(2).significant, forplots(2).HundredFifty,10,'yellow', "filled"); niet een significant aantal significante tuning curves --> eruit laten
%scatter(orientall(3).significant, forplots(3).HundredFifty,10,'black', "filled"); niet een significant aantal significante tuning curves --> eruit laten
scatter(orientall(4).significant, forplots(4).HundredFifty,10,'black', "filled");
scatter(orientall(5).significant, forplots(5).HundredFifty,10,'black', "filled");
scatter(orientall(6).significant, forplots(6).HundredFifty,10,'black', "filled");
scatter(orientall(7).significant, forplots(7).HundredFifty,10,'black', "filled");
scatter(orientall(8).significant, forplots(8).HundredFifty,10,'black', "filled");
title("157.5° Stimulus")
ylabel("Activity (Z-scored)")
xlabel("Preferred orientation")
%legend ("dataset 1", "dataset 4", "dataset 5", "dataset 6", "dataset 7", "dataset 8")
hold off

print(gcf, '-bestfit','tuning_scatterplots','-dpdf')