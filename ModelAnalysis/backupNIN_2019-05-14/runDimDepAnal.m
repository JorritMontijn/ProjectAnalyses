strOutputDir = 'D:\Data\Results\DimDep\';
load([strOutputDir 'DimDepOriDec9_2017__1_16.mat']);
load([strOutputDir 'DimDepOriDecShuff9_2017__1_16.mat']);


figure
subplot(2,2,1)
plot([0 numel(vecDimensionalityDepLogis)],(1/12)*[1 1],'k--');
hold on
plot(vecDimensionalityDepLogis,'b')
plot(vecDimensionalityDepLogisShuff,'b--')
hold off
title('logistic')
ylim([0 1])

subplot(2,2,2)
plot([0 numel(vecDimensionalityDepLogis)],(1/12)*[1 1],'k--');
hold on
plot(vecDimensionalityDepMahal,'r')
plot(vecDimensionalityDepMahalShuff,'r--')
hold off
title('mahal')
ylim([0 1])

%%
figure
subplot(2,2,1)
plot(vecDimensionalityDepLogis(1:10),'b')
title('logistic')
xlabel('Dimensionality')
ylabel('Decoding accuracy')
fixfig

subplot(2,2,2)
plot(vecDimensionalityDepMahal(1:10),'r')
title('mahalanobis')
xlabel('Dimensionality')
ylabel('Decoding accuracy')
fixfig
