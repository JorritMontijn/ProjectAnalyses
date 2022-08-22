close all
%% plot 3D example figures
%% H0: uniform around axis
%subplot(2,2,1)
figure
vecX = linspace(0,1,10);
vecY = vecX;
vecZ = vecY;
vecR = 0.02;
vecC = vecZ;

intPointsPerRadius = 27;


[handle,x,y,z,c] = plotTube3D(vecX,vecY,vecZ,vecR,vecC,intPointsPerRadius);
colormap(hot);
hold on
surf(x,y,0*z,c);

%calculate locations for circle around axis
vecV = [1 1 1];
vecV = vecV/norm(vecV,2);
vecA = [-0.5*sqrt(2) -0.5*sqrt(2) sqrt(2)];
vecA = vecA/norm(vecA,2);
vecB = cross(vecA,vecV);

intPoints = 100;
dblStep = ((pi/intPoints)*2);
vecTheta= -pi:dblStep:(pi+dblStep);
vecTheta = [vecTheta vecTheta(1)];

%near-ring
dblRadius = 0.1;
dblDist = 0.5;
dblWidth = 0.05;
vecCircleX = dblDist+dblRadius*cos(vecTheta)*vecA(1)+dblRadius*sin(vecTheta)*vecB(1);
vecCircleY = dblDist+dblRadius*cos(vecTheta)*vecA(2)+dblRadius*sin(vecTheta)*vecB(2);
vecCircleZ = dblDist+dblRadius*cos(vecTheta)*vecA(3)+dblRadius*sin(vecTheta)*vecB(3);
[handle,x,y,z,c] = tubeplot([vecCircleX;vecCircleY;vecCircleZ],dblWidth,0.5,27,0.0);
surf(x,y,0*z,c);

%far-ring
dblRadius = 0.3;
dblDist = 0.5;
dblWidth = 0.05;
vecCircleX = dblDist+dblRadius*cos(vecTheta)*vecA(1)+dblRadius*sin(vecTheta)*vecB(1);
vecCircleY = dblDist+dblRadius*cos(vecTheta)*vecA(2)+dblRadius*sin(vecTheta)*vecB(2);
vecCircleZ = dblDist+dblRadius*cos(vecTheta)*vecA(3)+dblRadius*sin(vecTheta)*vecB(3);
[handle,x,y,z,c] = tubeplot([vecCircleX;vecCircleY;vecCircleZ],dblWidth,0.5,27,0.0);
surf(x,y,0*z,c);

colormap(hot);
xlabel('Norm. act. neuron 1')
ylabel('Norm. act. neuron 2')
zlabel('Norm. act. neuron 3')

view([-26,14]);
axis tight
shading interp;
lighting gouraud


%% H1: non-uniform around axis
%subplot(2,2,2)
figure
[handle,matX_final,matY_final,matZ_final,matC] = plotTube3D(vecX,vecY,vecZ,vecR,vecC,intPointsPerRadius);
colormap(hot);
hold on
surf(matX_final,matY_final,0*matZ_final,matC);

%calculate locations for circle around axis
vecV = [1 1 1];
vecV = vecV/norm(vecV,2);
vecA = [-0.5*sqrt(2) -0.5*sqrt(2) sqrt(2)];
vecA = vecA/norm(vecA,2);
vecB = cross(vecA,vecV);


dblStep = ((pi/intPoints)*2);
vecTheta= -pi:dblStep:(pi+dblStep);
vecTheta = [vecTheta vecTheta(1)];

%near blob
dblRadius = 0.1;
dblDist = 0.5;
dblWidth = 0.07;
intPoints = 100;
[x,y,z] = sphere(intPoints);
x=x*dblWidth+dblDist;
y=y*dblWidth+dblDist;
z=z*dblWidth+dblDist;

x=x+vecA(3)*dblRadius;
y=y+vecA(1)*dblRadius;
z=z+vecA(2)*dblRadius;
c = 0.5*ones(size(x));
surf(x,y,z,c);
surf(x,y,zeros(size(x)),c);


%far blob
dblRadius = 0.3;
dblDist = 0.5;
dblWidth = 0.07;
intPoints = 100;
[x,y,z] = sphere(intPoints);
x=x*dblWidth+dblDist;
y=y*dblWidth+dblDist;
z=z*dblWidth+dblDist;

x=x+vecB(1)*dblRadius;
y=y+vecB(2)*dblRadius;
z=z+vecB(3)*dblRadius;
c = 0.5*ones(size(x));
surf(x,y,z,c)
surf(x,y,zeros(size(x)),c);

xlabel('Norm. act. neuron 1')
ylabel('Norm. act. neuron 2')
zlabel('Norm. act. neuron 3')
colormap(hot);
axis xy

%camlight(-170,-56)
%camlight(-170,-56)
%camlight(-215,70)
%view([-10,22]);
%view([-150,16]);
%view([-60,26]);
view([-26,14]);
axis tight
shading interp;
lighting gouraud

axis off


%create background axes
figure
vecLim = [-0.0163 1.0163];
xlim(vecLim);
ylim(vecLim);
zlim(vecLim);
view([-26,14]);
set(gca,'xtick',0:0.2:1)
set(gca,'ytick',0:0.2:1)
set(gca,'ztick',0:0.2:1)

xlabel('Activity neuron 1')
ylabel('Activity neuron 2')
zlabel('Activity neuron 3')