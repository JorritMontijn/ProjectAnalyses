strRunStim = 'DG';
runHeaderPopTimeCoding;
%%
intRec=1;%1
runRecPrepNpx;
intNumN = intRespN;
%%
% make plot
dblStartT = 2000;%vecStimOnTime(1);
dblPlotDur = 500;%vecStimOffTime(end) - vecStimOnTime(1);
dblBinW = 0.1;
vecEdgesT = dblStartT:dblBinW:(dblStartT+dblPlotDur);
vecCenterT = vecEdgesT(2:end)-dblBinW/2;
intNumT = numel(vecCenterT);
matPlot = nan(intNumN,intNumT);
for intN=1:intNumN
	matPlot(intN,:) = (0+histcounts(cellSpikeTimes{intN},vecEdgesT));
end
vecMeanR = mean(matPlot,2)/dblBinW;
indRem = vecMeanR<0;
intPlotN = sum(~indRem);
matPlot(indRem,:)=[];
% plot

% filter
dblScale = 500;
dblFactor = (dblPlotDur/dblBinW)/dblScale;
vecFilt = normpdf((-2*dblFactor):(2*dblFactor),0,dblFactor);
matFilt = vecFilt'*vecFilt;
matFilt = matFilt./sum(matFilt(:));
matPlotFilt = imfilt(matPlot,vecFilt,0);
vecMinR =  min(matPlotFilt,[],2);
vecMaxR =  max(matPlotFilt,[],2);
matPlotFilt = (matPlotFilt-vecMinR)./(vecMaxR-vecMinR);

%sort
vecNewOrder = clustsort(matPlot');

figure

imagesc(vecCenterT,1:intPlotN,matPlotFilt(vecNewOrder,:));
colormap(flipud(bone));
xlabel('Time (s)')
ylabel('Neuron #');
title('Norm. act (Hz)');
colorbar
fixfig;

%% make lognormal plot
vecOut = logmvnrnd(4,4,10000);
vecE = 0:0.5:15;
vecC = histcounts(vecOut,vecE);
stairs(vecE(2:end)-diff(vecE(1:2))/2,vecC);

%% plot constant sd/mean with overlap
dblCV = 0.8;

intNumQ = 5;
vecGains = 1:intNumQ;


vecMu = vecGains;
vecSd = vecMu*dblCV;

dblRange = ceil(max(vecMu) + max(vecSd)*2);
dblStep = 0.001;
vecBinE = -dblRange:dblStep:dblRange;
vecBinC = vecBinE(2:end)-diff(vecBinE(1:2))/2;
figure
for i=1:intNumQ
	subplot(intNumQ,1,intNumQ-i+1); hold on
	
	%s1
	vecP1 = normpdf(vecBinC,vecMu(i),vecSd(i));
	vecP1 = vecP1./sum(vecP1(:));
	plot(vecBinC,vecP1,'r');
	
	%s2
	vecP2 = normpdf(vecBinC,-vecMu(i),vecSd(i));
	vecP2 = vecP2./sum(vecP2(:));
	plot(vecBinC,vecP2,'k');
	
	%find decision boundary
	normcdf(0,vecMu(i),vecSd(i))
	
	%calc overlap
	dblOverlap = sum(vecP2.*vecP1) ./ sum(vecP2.^2);
	ylim([0 5e-4]);
end

%% test
inum=80;
vecP1stP2 = false(1,inum);
vecP1 = normrnd(vecMu(i),vecSd(i),[1,inum]);
vecP2 = normrnd(-vecMu(i),vecSd(i),[1,inum]);
dblError_Emp = (sum(vecP2>=0)+sum(vecP1<=0))/(2*inum)
dblError_MuSd = normcdf(-vecMu(i)/vecSd(i))
dblError_CV = normcdf(-1/dblCV,0,1)

vecX = 0:inum;
vecP = binopdf(vecX,inum,dblError_MuSd);
plot(vecX/inum,vecP);

%% plot ring manifold
vecP = linspace(0,2*pi,50);
dblFactor = 0.3;
dblConeWidth = 1/8;
figure;
hold on
vecFactor = [0.1:0.1:0.5];
dblPlaneHeight = 0;

%build base cone
vecX = dblFactor*cos(vecP);
vecY = dblFactor*sin(vecP);
vecZ = dblFactor/sin(dblConeWidth*pi)+zeros(size(vecP));%dblFactor*2+2*dblFactor*sin(vecP);%dblFactor*sin(vecP)+dblFactor;

%rotate 45 degrees along x axis and y axis
%get rotation
dblYawDV = 0; %yaw
dblPitchML = 45; %pitch
dblRollAP = -45; %roll

% build rotation matrix in yaw pitch roll
a = deg2rad(dblYawDV);
c = deg2rad(dblPitchML);
b = deg2rad(dblRollAP);
matR = [...
	cos(a)*cos(b)	cos(a)*sin(b)*sin(c)-sin(a)*cos(c)	cos(a)*sin(b)*cos(c)+sin(a)*sin(c);...
	sin(a)*cos(b)	sin(a)*sin(b)*sin(c)+cos(a)*cos(c)	sin(a)*sin(b)*cos(c)-cos(a)*sin(c);...
	-sin(b)			cos(b)*sin(c)						cos(b)*cos(c)]';

%transform to points, rotate, and transform back
matP = cat(2,vecX(:),vecY(:),vecZ(:))';
matRotP = matR*matP;
vecX = reshape(matRotP(1,:)',size(vecX));
vecY = reshape(matRotP(2,:)',size(vecY));
vecZ = reshape(matRotP(3,:)',size(vecZ));

%move slice to atlas space
%vecX = round(vecX + sSlice.Center(1)); %ML
%vecY = round(vecY + sSlice.Center(2)); %AP
%vecZ = round(vecZ + sSlice.Center(3)); %DV

[pnts,conn,line1,line2] = generate_open_tube(0,0.03,9,vecX,vecY,vecZ);
matC = circcol(length(conn));
for k=1:length(conn)
	patch('Faces',     conn{k}, ...
		'Vertices',  pnts', ...
		'FaceAlpha', 0.8, ...
		'FaceColor', matC(k,:), ...
		'EdgeColor', 'none' ) ;
end
colormap(matC)
%xlim([-1.1 1.1]);
%ylim([-1.1 1.1]);
%zlim([-1.1 1.1]);
set(gcf,'Renderer','zbuffer')
drawnow;
%%{
%patch([-1.1 -1.1 1.1 1.1],[-1.1 1.1 1.1 -1.1],dblPlaneHeight*[1 1 1 1],[0.9 0.9 0.9],...
%	'FaceAlpha',0.2,'FaceColor',[0.5 0.5 0.5],'EdgeColor','none');
h= cline(vecX,vecY,dblPlaneHeight*ones(size(vecX)),1:length(vecX),true);
set(h,'LineWidth',4,'EdgeColor','interp')
%}
%%{
%lighting and stuff

hold off
grid on

light('Position',[5 -5 7],'Style','infinite');
lighting gouraud

xlabel('Norm. act. neuron 1')
ylabel('Norm. act. neuron 2')
zlabel('Norm. act. neuron 3')
xlabel(get(get(gca,'xlabel'), 'String'),'FontSize',24); %set x-label and change font size
ylabel(get(get(gca,'ylabel'), 'String'),'FontSize',24);%set y-label and change font size
zlabel(get(get(gca,'zlabel'), 'String'),'FontSize',24); %set x-label and change font size
dblFontSize = 18;
set(gca,'FontSize',dblFontSize,'Linewidth',2); %set grid line width and change font size of x/y ticks

vecCamPos = [7.5210   -3.5373    7.0120];
vecCamPos2 = [7.8771   -3.9503    6.3114];

set(gca,'cameraposition',vecCamPos2)
axis equal;
vecLim = [0 max([get(gca,'xlim') get(gca,'ylim') get(gca,'zlim')])];
vecLim = [0 1.5];
xlim(vecLim);
ylim(vecLim);
zlim(vecLim);
drawnow;

%% plot cone manifold
intNumP = 50;
vecP = linspace(0,2*pi,intNumP);
dblFactor = 0.5;
dblConeWidth = 1/8;
figure;
hold on
matCol = redbluepurple(5);
vecFactor = [0.1:0.1:0.5];
dblPlaneHeight = 0;

%set rotation
dblYawDV = 0; %yaw
dblPitchML = 45; %pitch
dblRollAP = -45; %roll

% build rotation matrix in yaw pitch roll
a = deg2rad(dblYawDV);
c = deg2rad(dblPitchML);
b = deg2rad(dblRollAP);
matR = [...
	cos(a)*cos(b)	cos(a)*sin(b)*sin(c)-sin(a)*cos(c)	cos(a)*sin(b)*cos(c)+sin(a)*sin(c);...
	sin(a)*cos(b)	sin(a)*sin(b)*sin(c)+cos(a)*cos(c)	sin(a)*sin(b)*cos(c)-cos(a)*sin(c);...
	-sin(b)			cos(b)*sin(c)						cos(b)*cos(c)]';
	
%build base cone
r = linspace(0,0.55,intNumP) ;
th = linspace(0,2*pi,intNumP) ;
[R,T] = meshgrid(r,th) ;
X = R.*cos(T) ;
Y = R.*sin(T) ;
Z = (R./sin(dblConeWidth*pi));

%rotate cone
matP = cat(2,X(:),Y(:),Z(:))';
matRotP = matR*matP;
X = reshape(matRotP(1,:)',size(X));
Y = reshape(matRotP(2,:)',size(Y));
Z = reshape(matRotP(3,:)',size(Z));

colormap(redbluepurple(intNumP));
C = repmat(r,[intNumP 1]);
%colormap(circcol(intNumP));
%C = repmat(th',[1 intNumP]);
%h=surf(X,Y,Z,C,'edgecolor','none','facealpha',0.5);

for i=[]%3%1:5%[1 2 4]%3 5];%1:numel(vecFactor)
	
	
	%build base ring
	dblFactor = vecFactor(i);
	vecX = dblFactor*cos(vecP);
	vecY = dblFactor*sin(vecP);
	vecZ = dblFactor/sin(dblConeWidth*pi)+zeros(size(vecP));%dblFactor*2+2*dblFactor*sin(vecP);%dblFactor*sin(vecP)+dblFactor;
	
	%transform to points, rotate, and transform back
	matP = cat(2,vecX(:),vecY(:),vecZ(:))';
	matRotP = matR*matP;
	vecX = reshape(matRotP(1,:)',size(vecX));
	vecY = reshape(matRotP(2,:)',size(vecY));
	vecZ = reshape(matRotP(3,:)',size(vecZ));
	
	%move slice to atlas space
	%vecX = round(vecX + sSlice.Center(1)); %ML
	%vecY = round(vecY + sSlice.Center(2)); %AP
	%vecZ = round(vecZ + sSlice.Center(3)); %DV
	dblWidthFactor = ((dblFactor/0.3)-1)/2+1;
	[pnts,conn,line1,line2] = generate_open_tube(0,0.03*dblWidthFactor,9,vecX,vecY,vecZ);
	matC = circcol(length(conn));
	for k=1:length(conn)
		patch('Faces',     conn{k}, ...
			'Vertices',  pnts', ...
			'FaceAlpha', 0.8, ...
			...'FaceColor',matCol(i,:), ...
			'FaceColor',matC(k,:), ...
			'EdgeColor', 'none' ) ;
	end
end
%colormap(matC)
%xlim([-1.1 1.1]);
%ylim([-1.1 1.1]);
%zlim([-1.1 1.1]);
set(gcf,'Renderer','zbuffer')
drawnow;
%patch([-1.1 -1.1 1.1 1.1],[-1.1 1.1 1.1 -1.1],dblPlaneHeight*[1 1 1 1],[0.9 0.9 0.9],...
%	'FaceAlpha',0.2,'FaceColor',[0.5 0.5 0.5],'EdgeColor','none');
%h= cline(vecX,vecY,dblPlaneHeight*ones(size(vecX)),1:length(vecX),true);
h=surf(X,Y,0*Z,C,'edgecolor','none','facealpha',0.5);

set(h,'LineWidth',4,'EdgeColor','interp')

%lighting and stuff

hold off
grid on

light('Position',[5 -5 7],'Style','infinite');
lighting gouraud

xlabel('Norm. act. neuron 1')
ylabel('Norm. act. neuron 2')
zlabel('Norm. act. neuron 3')
xlabel(get(get(gca,'xlabel'), 'String'),'FontSize',24); %set x-label and change font size
ylabel(get(get(gca,'ylabel'), 'String'),'FontSize',24);%set y-label and change font size
zlabel(get(get(gca,'zlabel'), 'String'),'FontSize',24); %set x-label and change font size
dblFontSize = 18;
set(gca,'FontSize',dblFontSize,'Linewidth',2); %set grid line width and change font size of x/y ticks

vecCamPos = [7.5210   -3.5373    7.0120];
vecCamPos2 = [7.8771   -3.9503    6.3114];

set(gca,'cameraposition',vecCamPos2)
axis equal;
vecLim = [0 max([get(gca,'xlim') get(gca,'ylim') get(gca,'zlim')])];
vecLim = [0 1.5];
xlim(vecLim);
ylim(vecLim);
zlim(vecLim);
drawnow;
