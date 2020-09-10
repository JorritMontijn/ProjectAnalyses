PV_orient_places = find(ismember(OrientVect, OrientPVVect));
for j = 1:length(OrientPVVect)
    [nRows, nCols] = size(orientpvactall(PV_orient_places(j)).activ);
    for i = 1:nRows
        forplotsPV(j).Zero(i)             =                  mean(orientpvactall(PV_orient_places(j)).activZ(i,find(dat(OrientPVVect(j)).ses.structStim.Orientation == 0)));
        forplotsPV(j).Twenty(i)           =                  mean(orientpvactall(PV_orient_places(j)).activZ(i,find(dat(OrientPVVect(j)).ses.structStim.Orientation == 22.5)));
        forplotsPV(j).Forty(i)            =                  mean(orientpvactall(PV_orient_places(j)).activZ(i,find(dat(OrientPVVect(j)).ses.structStim.Orientation == 45)));
        forplotsPV(j).Sixty(i)            =                  mean(orientpvactall(PV_orient_places(j)).activZ(i,find(dat(OrientPVVect(j)).ses.structStim.Orientation == 67.5)));
        forplotsPV(j).Ninety(i)           =                  mean(orientpvactall(PV_orient_places(j)).activZ(i,find(dat(OrientPVVect(j)).ses.structStim.Orientation == 90)));
        forplotsPV(j).HundredTen(i)       =                  mean(orientpvactall(PV_orient_places(j)).activZ(i,find(dat(OrientPVVect(j)).ses.structStim.Orientation == 112.5)));
        forplotsPV(j).HundredThirty(i)    =                  mean(orientpvactall(PV_orient_places(j)).activZ(i,find(dat(OrientPVVect(j)).ses.structStim.Orientation == 135)));
        forplotsPV(j).HundredFifty(i)     =                  mean(orientpvactall(PV_orient_places(j)).activZ(i,find(dat(OrientPVVect(j)).ses.structStim.Orientation == 157.5)));
    end
end

for j = 1:length(orientPVall)
    for i = 1:length(orientPVall(j).degrees)
        if orientPVall(j).degrees(i) <= 90
            shifttest_0(j).degrees(i) = orientPVall(j).degrees(i)+90;
        else
            shifttest_0(j).degrees(i) = orientPVall(j).degrees(i)- 90;
        end
    end
end
for j = 1:length(orientPVall)
    for i = 1:length(orientPVall(j).degrees)
        if orientPVall(j).degrees(i) <= 90
            shifttest_22(j).degrees(i) = orientPVall(j).degrees(i)+67.5;
        else
            shifttest_22(j).degrees(i) = orientPVall(j).degrees(i)- 67.5;
        end
    end
end
for j = 1:length(orientPVall)
    for i = 1:length(orientPVall(j).degrees)
        if orientPVall(j).degrees(i) <= 90
            shifttest_45(j).degrees(i) = orientPVall(j).degrees(i)+45;
        else
            shifttest_45(j).degrees(i) = orientPVall(j).degrees(i)- 45;
        end
    end
end

for j = 1:length(orientPVall)
    for i = 1:length(orientPVall(j).degrees)
        if orientPVall(j).degrees(i) <= 90
            shifttest_67(j).degrees(i) = orientPVall(j).degrees(i)+22.5;
        else
            shifttest_67(j).degrees(i) = orientPVall(j).degrees(i)- 22.5;
        end
    end
end

for j = 1:length(orientPVall)
    for i = 1:length(orientPVall(j).degrees)
        if orientPVall(j).degrees(i) <= 90
            shifttest_90(j).degrees(i) = orientPVall(j).degrees(i);
        else
            shifttest_90(j).degrees(i) = orientPVall(j).degrees(i);
        end
    end
end

for j = 1:length(orientPVall)
    for i = 1:length(orientPVall(j).degrees)
        if orientPVall(j).degrees(i) >= 22.5
            shifttest_112(j).degrees(i) = orientPVall(j).degrees(i)-22.5;
        else
            shifttest_112(j).degrees(i) = orientPVall(j).degrees(i) +157.5;
        end
    end
end

for j = 1:length(orientPVall)
    for i = 1:length(orientPVall(j).degrees)
        if orientPVall(j).degrees(i) >= 45
            shifttest_135(j).degrees(i) = orientPVall(j).degrees(i)-45;
        else
            shifttest_135(j).degrees(i) = orientPVall(j).degrees(i)+135;
        end
    end
end

for j = 1:length(orientPVall)
    for i = 1:length(orientPVall(j).degrees)
        if orientPVall(j).degrees(i) >= 67.5
            shifttest_157(j).degrees(i) = orientPVall(j).degrees(i)-67.5;
        else
            shifttest_157(j).degrees(i) = orientPVall(j).degrees(i)+112.5;
        end
    end
end

for j = 1:5
    scatter(shifttest_0(j).degrees, forplotsPV(j).Zero,10,'black', "filled");
    hold on
    scatter(shifttest_22(j).degrees, forplotsPV(j).Twenty,10,'black', "filled");
    scatter(shifttest_45(j).degrees, forplotsPV(j).Forty,10,'black', "filled");
    scatter(shifttest_67(j).degrees, forplotsPV(j).Sixty,10,'black', "filled");
    scatter(shifttest_90(j).degrees, forplotsPV(j).Ninety,10,'black', "filled");
    scatter(shifttest_112(j).degrees, forplotsPV(j).HundredTen,10,'black', "filled");
    scatter(shifttest_135(j).degrees, forplotsPV(j).HundredThirty,10,'black', "filled");
    scatter(shifttest_157(j).degrees, forplotsPV(j).HundredFifty,10,'black', "filled");
end


xticklabels(linspace(-90, 90, 10))
title("PV-cell activity response - circular-shifted stimuli")
ylabel("Activity (Z-scored)")
xlabel("Prefered orientation")
hold off

print(gcf, '-bestfit','PVshift_final_black','-dpdf')
