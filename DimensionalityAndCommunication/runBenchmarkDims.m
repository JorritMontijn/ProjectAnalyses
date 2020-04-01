clear all;
vecSubspXY = 1:5;
vecSubspYX = 6:10;

intN = 20;
intT = 50;
intSharedN = max([vecSubspYX vecSubspXY]);
intPrivSepN = 5;
intPrivShN = intN-intSharedN-intPrivSepN;
matX1 = cat(2,rand(intT,intSharedN),0.001*rand(intT,intPrivSepN),rand(intT,intPrivShN));
matX2 = cat(2,rand(intT,intSharedN),0.001*rand(intT,intPrivSepN),rand(intT,intPrivShN));
matY1 = cat(2,rand(intT,intSharedN),rand(intT,intPrivShN),0.001*rand(intT,intPrivSepN));
matY2 = cat(2,rand(intT,intSharedN),rand(intT,intPrivShN),0.001*rand(intT,intPrivSepN));

matTransferXY = zeros(1,intN);
matTransferXY(vecSubspXY) = 1;
matTransferXY = diag(matTransferXY);
matTXY1 = matX1 * matTransferXY;
matTXY2 = matX2 * matTransferXY;
matY1 = matY1 + matTXY1;
matY2 = matY2 + matTXY2;

matTransferYX = zeros(1,intN);
matTransferYX(vecSubspYX) = 1;
matTransferYX = diag(matTransferYX);
matTYX1 = matY1 * matTransferYX;
matTYX2 = matY2 * matTransferYX;
matX1 = matX1 + matTYX1;
matX2 = matX2 + matTYX2;

dblLambda = 1;

[vecR2XY,vecR2YX,matSubspace_Shared,matSubspace_PrivateX,matSubspace_PrivateY,matSeparator,vecR2PrivX,vecR2PrivY] = ...
	getSpacesSharedPrivateCV(matX1,matY1,matX2,matY2,dblLambda);

%% plot
figure
subplot(4,6,1)
matCovX = cov(matX2);
imagesc(matCovX,max(abs(matCovX(:)))*[-1 1]);
title('A1) Cov(X)')
fixfig;grid off;

subplot(4,6,2)
imagesc(matTransferXY,max(abs(matTransferXY(:)))*[-1 1]);
title('A2) Transfer(X->Y)')
fixfig;grid off;

subplot(4,6,7)
matCovY = cov(matY2);
imagesc(matCovY,max(abs(matCovY(:)))*[-1 1]);
title('A3) Cov(Y)')
fixfig;grid off;

subplot(4,6,8)
imagesc(matTransferYX,max(abs(matTransferYX(:)))*[-1 1]);
title('A4) Transfer(Y->X)')
fixfig;grid off;

%check magnitudes
vecXS = mean(abs(matX2 * matSubspace_Shared),1)
vecXPX = mean(abs(matX2 * matSubspace_PrivateX),1)
vecXPY = mean(abs(matX2 * matSubspace_PrivateY),1)

vecYS = mean(abs(matY2 * matSubspace_Shared),1)
vecYPX = mean(abs(matY2 * matSubspace_PrivateX),1)
vecYPY = mean(abs(matY2 * matSubspace_PrivateY),1)

subplot(2,3,2)
plot(vecXS,'m')
hold on
plot(vecXPX,'r')
plot(vecXPY,'b')
hold off
xlabel('Dimension (neuron #)');
title('B) X proj in subsp')
legend({'Shared','X-priv','Y-priv',},'Location','Best')
fixfig;

subplot(2,3,3)
plot(vecYS,'m')
hold on
plot(vecYPX,'r')
plot(vecYPY,'b')
hold off
xlabel('Dimension (neuron #)');
title('C) Y proj in subsp')
legend({'Shared','X-priv','Y-priv',},'Location','Best')
fixfig;

subplot(2,3,4)
imagesc(matSubspace_Shared,max(abs(matSubspace_Shared(:)))*[-1 1])
title('D) Shared subspace');
xlabel('Dimension (neuron #)');
ylabel('Dimension (neuron #)');
colormap(redblue)
fixfig;grid off;

subplot(2,3,5)
imagesc(matSubspace_PrivateX,max(abs(matSubspace_PrivateX(:)))*[-1 1])
title('E) Private X subspace');
fixfig;grid off;
xlabel('Dimension (neuron #)');
ylabel('Dimension (neuron #)');

subplot(2,3,6)
imagesc(matSubspace_PrivateY,max(abs(matSubspace_PrivateY(:)))*[-1 1])
title('F) Private Y subspace');
fixfig;grid off;
xlabel('Dimension (neuron #)');
ylabel('Dimension (neuron #)');

maxfig;
%export_fig('Example_SharedPrivate.tif')
%export_fig('Example_SharedPrivate.pdf')

%check magnitudes
dblXS = mean(vecXS(:))
dblXPX = mean(vecXPX(:))
dblXPY = mean(vecXPY(:))

dblYS = mean(vecYS(:))
dblYPX = mean(vecYPX(:))
dblYPY = mean(vecYPY(:))
