vecNoise = 10.^[-5:-2];
vecNoise = [0 0.00003 0.0001];%*10
fSat = @(x,c) x./(1+c.*x);

vecPopSizes = [1:10000];
matR = fSat(vecPopSizes,vecNoise')';
dblMaxVal = max(matR(:));
matC = redbluepurple(7);

clf;
hold on;
plot(vecPopSizes/1000,matR(:,1)./dblMaxVal,'Color',matC(6,:));
plot(vecPopSizes/1000,matR(:,2)./dblMaxVal,'Color',matC(4,:));
plot(vecPopSizes/1000,matR(:,3)./dblMaxVal,'Color',matC(2,:));
plot([vecPopSizes(1) vecPopSizes(end)]./1000,0.4.*[1 1],'k--');
text((vecPopSizes(end)/1000)*0.7,0.9,'unsat','FontSize',18,'Color',matC(6,:))
text((vecPopSizes(end)/1000)*0.8,0.6,'sat, high','FontSize',18,'Color',matC(4,:))
text((vecPopSizes(end)/1000)*0.6,0.33,'sat, low','FontSize',18,'Color',matC(2,:))
text((vecPopSizes(end)/1000)*0.05,0.46,'Behaviour','FontSize',18,'Color',[0 0 0])

hold off;
xlabel('# of neurons in analysis (x 1000)');
ylabel('Total V1 information (a.u.)');
fixfig;grid off;

export_fig('Fig1_raw.tif');
export_fig('Fig1_raw.pdf');