clear all;
intP = 1000;
vecT = linspace(0,1,intP)*1.5-0.3;
dblStep=(2*pi)/intP;
vecTheta = dblStep:dblStep:(2*pi);
matCol = circcol(intP);
intShiftCol = round(intP*(0.3/1.5));
matCol = circshift(matCol,[-(intShiftCol-1) 0]);
%h=cline(cos(vecTheta),sin(vecTheta),matCol,true);
%arrayfun(@(x) set(x,'LineWidth',50),h)

%patches
figure;
hAx=axes;
colormap(hAx,matCol);colorbar;

vE=[1:intP];
vX = [cos(vecTheta(vE))';];
vY = [sin(vecTheta(vE))';];
mC = [matCol(vE,:);];


%%
dbl0 = (intShiftCol/intP)*2*pi;
make_colorwheel(matCol)
%hold on
%polarplot(dbl0,1.9,'xb');


%% save plot in selected format
print(['PLOT_colwheel_' getDate],'-dpdf')