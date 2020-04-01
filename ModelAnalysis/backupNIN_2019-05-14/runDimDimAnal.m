strOutputDir = 'D:\Data\Results\DimDep\';
load([strOutputDir 'DimDepOriDec9_2017__1_16.mat']);
load([strOutputDir 'DimDepOriDecShuff9_2017__1_16.mat']);


figure
subplot(2,2,1)
plot(vecDimensionalityDepLogis,'b')
hold on
plot(vecDimensionalityDepLogisShuff,'b--')
hold off
title('logistic')

subplot(2,2,2)
plot(vecDimensionalityDepMahal,'r')
hold on
plot(vecDimensionalityDepMahalShuff,'r--')
hold off
title('mahal')
