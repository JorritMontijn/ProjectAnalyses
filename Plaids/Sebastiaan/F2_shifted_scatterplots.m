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
%%For every orientation I shift the Tuning Curves to 90 degrees:
%0 degrees
for j = 1:length(orientall)
    for i = 1:length(orientall(j).significant)
        if orientall(j).significant(i) <= 90
            shifttest_0(j).significant(i) = orientall(j).significant(i)+90;
        else
            shifttest_0(j).significant(i) = orientall(j).significant(i)- 90;
        end
    end
end
%22.5 degrees
for j = 1:length(orientall)
    for i = 1:length(orientall(j).significant)
        if orientall(j).significant(i) <= 90
            shifttest_22(j).significant(i) = orientall(j).significant(i)+67.5;
        else
            shifttest_22(j).significant(i) = orientall(j).significant(i)- 67.5;
        end
    end
end
%45 degrees
for j = 1:length(orientall)
    for i = 1:length(orientall(j).significant)
        if orientall(j).significant(i) <= 90
            shifttest_45(j).significant(i) = orientall(j).significant(i)+45;
        else
            shifttest_45(j).significant(i) = orientall(j).significant(i)- 45;
        end
    end
end
%67.5 degrees
for j = 1:length(orientall)
    for i = 1:length(orientall(j).significant)
        if orientall(j).significant(i) <= 90
            shifttest_67(j).significant(i) = orientall(j).significant(i)+22.5;
        else
            shifttest_67(j).significant(i) = orientall(j).significant(i)- 22.5;
        end
    end
end
%90 degrees
for j = 1:length(orientall)
    for i = 1:length(orientall(j).significant)
        if orientall(j).significant(i) <= 90
            shifttest_90(j).significant(i) = orientall(j).significant(i);
        else
            shifttest_90(j).significant(i) = orientall(j).significant(i);
        end
    end
end
%112.5 degrees
for j = 1:length(orientall)
    for i = 1:length(orientall(j).significant)
        if orientall(j).significant(i) >= 22.5
            shifttest_112(j).significant(i) = orientall(j).significant(i)-22.5;
        else
            shifttest_112(j).significant(i) = orientall(j).significant(i) +157.5;
        end
    end
end
%135 degrees
for j = 1:length(orientall)
    for i = 1:length(orientall(j).significant)
        if orientall(j).significant(i) >= 45
            shifttest_135(j).significant(i) = orientall(j).significant(i)-45;
        else
            shifttest_135(j).significant(i) = orientall(j).significant(i)+135;
        end
    end
end
%157.5 degrees
for j = 1:length(orientall)
    for i = 1:length(orientall(j).significant)
        if orientall(j).significant(i) >= 67.5
            shifttest_157(j).significant(i) = orientall(j).significant(i)-67.5;
        else
            shifttest_157(j).significant(i) = orientall(j).significant(i)+112.5;
        end
    end
end


scatter(shifttest_0(1).significant, forplots(1).Zero,10,'black', "filled");
hold on
%scatter(shifttest_0(2).significant, forplots(2).Zero,10,'yellow', "filled"); niet een significant aantal significante tuning curves --> eruit laten
%scatter(shifttest_0(3).significant, forplots(3).Zero,10,'black', "filled"); niet een significant aantal significante tuning curves --> eruit laten
scatter(shifttest_0(4).significant, forplots(4).Zero,10,'black', "filled");
scatter(shifttest_0(5).significant, forplots(5).Zero,10,'black', "filled");
scatter(shifttest_0(6).significant, forplots(6).Zero,10,'black', "filled");
scatter(shifttest_0(7).significant, forplots(7).Zero,10,'black', "filled");
scatter(shifttest_0(8).significant, forplots(8).Zero,10,'black', "filled");

scatter(shifttest_22(1).significant, forplots(1).Twenty,10,'black', "filled");
%scatter(shifttest_22(2).significant, forplots(2).Twenty,10,'yellow', "filled"); niet een significant aantal significante tuning curves --> eruit laten
%scatter(shifttest_22(3).significant, forplots(3).Twenty,10,'black', "filled"); niet een significant aantal significante tuning curves --> eruit laten
scatter(shifttest_22(4).significant, forplots(4).Twenty,10,'black', "filled");
scatter(shifttest_22(5).significant, forplots(5).Twenty,10,'black', "filled");
scatter(shifttest_22(6).significant, forplots(6).Twenty,10,'black', "filled");
scatter(shifttest_22(7).significant, forplots(7).Twenty,10,'black', "filled");
scatter(shifttest_22(8).significant, forplots(8).Twenty,10,'black', "filled");

scatter(shifttest_45(1).significant, forplots(1).Forty,10,'black', "filled");
%scatter(shifttest_45(2).significant, forplots(2).Forty,10,'yellow', "filled"); niet een significant aantal significante tuning curves --> eruit laten
%scatter(shifttest_45(3).significant, forplots(3).Forty,10,'black', "filled"); niet een significant aantal significante tuning curves --> eruit laten
scatter(shifttest_45(4).significant, forplots(4).Forty,10,'black', "filled");
scatter(shifttest_45(5).significant, forplots(5).Forty,10,'black', "filled");
scatter(shifttest_45(6).significant, forplots(6).Forty,10,'black', "filled");
scatter(shifttest_45(7).significant, forplots(7).Forty,10,'black', "filled");
scatter(shifttest_45(8).significant, forplots(8).Forty,10,'black', "filled");

scatter(shifttest_67(1).significant, forplots(1).Sixty,10,'black', "filled");
%scatter(shifttest_67(2).significant, forplots(2).Sixty,10,'yellow', "filled"); niet een significant aantal significante tuning curves --> eruit laten
%scatter(shifttest_67(3).significant, forplots(3).Sixty,10,'black', "filled"); niet een significant aantal significante tuning curves --> eruit laten
scatter(shifttest_67(4).significant, forplots(4).Sixty,10,'black', "filled");
scatter(shifttest_67(5).significant, forplots(5).Sixty,10,'black', "filled");
scatter(shifttest_67(6).significant, forplots(6).Sixty,10,'black', "filled");
scatter(shifttest_67(7).significant, forplots(7).Sixty,10,'black', "filled");
scatter(shifttest_67(8).significant, forplots(8).Sixty,10,'black', "filled");

scatter(shifttest_90(1).significant, forplots(1).Ninety,10,'black', "filled");
%scatter(shifttest_90(2).significant, forplots(2).Ninety,10,'yellow', "filled"); niet een significant aantal significante tuning curves --> eruit laten
%scatter(shifttest_90(3).significant, forplots(3).Ninety,10,'black', "filled"); niet een significant aantal significante tuning curves --> eruit laten
scatter(shifttest_90(4).significant, forplots(4).Ninety,10,'black', "filled");
scatter(shifttest_90(5).significant, forplots(5).Ninety,10,'black', "filled");
scatter(shifttest_90(6).significant, forplots(6).Ninety,10,'black', "filled");
scatter(shifttest_90(7).significant, forplots(7).Ninety,10,'black', "filled");
scatter(shifttest_90(8).significant, forplots(8).Ninety,10,'black', "filled");

scatter(shifttest_112(1).significant, forplots(1).HundredTen,10,'black', "filled");
%scatter(shifttest_112(2).significant, forplots(2).HundredTen,10,'yellow', "filled"); niet een significant aantal significante tuning curves --> eruit laten
%scatter(shifttest_112(3).significant, forplots(3).HundredTen,10,'black', "filled"); niet een significant aantal significante tuning curves --> eruit laten
scatter(shifttest_112(4).significant, forplots(4).HundredTen,10,'black', "filled");
scatter(shifttest_112(5).significant, forplots(5).HundredTen,10,'black', "filled");
scatter(shifttest_112(6).significant, forplots(6).HundredTen,10,'black', "filled");
scatter(shifttest_112(7).significant, forplots(7).HundredTen,10,'black', "filled");
scatter(shifttest_112(8).significant, forplots(8).HundredTen,10,'black', "filled");

scatter(shifttest_135(1).significant, forplots(1).HundredThirty,10,'black', "filled");
%scatter(shifttest_135(2).significant, forplots(2).HundredThirty,10,'yellow', "filled"); niet een significant aantal significante tuning curves --> eruit laten
%scatter(shifttest_135(3).significant, forplots(3).HundredThirty,10,'black', "filled"); niet een significant aantal significante tuning curves --> eruit laten
scatter(shifttest_135(4).significant, forplots(4).HundredThirty,10,'black', "filled");
scatter(shifttest_135(5).significant, forplots(5).HundredThirty,10,'black', "filled");
scatter(shifttest_135(6).significant, forplots(6).HundredThirty,10,'black', "filled");
scatter(shifttest_135(7).significant, forplots(7).HundredThirty,10,'black', "filled");
scatter(shifttest_135(8).significant, forplots(8).HundredThirty,10,'black', "filled");

scatter(shifttest_157(1).significant, forplots(1).HundredFifty,10,'black', "filled");
%scatter(shifttest_157(2).significant, forplots(2).HundredFifty,10,'yellow', "filled"); niet een significant aantal significante tuning curves --> eruit laten
%scatter(shifttest_157(3).significant, forplots(3).HundredFifty,10,'black', "filled"); niet een significant aantal significante tuning curves --> eruit laten
scatter(shifttest_157(4).significant, forplots(4).HundredFifty,10,'black', "filled");
scatter(shifttest_157(5).significant, forplots(5).HundredFifty,10,'black', "filled");
scatter(shifttest_157(6).significant, forplots(6).HundredFifty,10,'black', "filled");
scatter(shifttest_157(7).significant, forplots(7).HundredFifty,10,'black', "filled");
scatter(shifttest_157(8).significant, forplots(8).HundredFifty,10,'black', "filled");
xticklabels(linspace(-90, 90, 10))
title("All stimuli - Circ-Shifted")
ylabel("Activity (Z-scored)")
xlabel("Preferred orientation")
hold off

print(gcf, '-bestfit','circshift_final_black','-dpdf')


