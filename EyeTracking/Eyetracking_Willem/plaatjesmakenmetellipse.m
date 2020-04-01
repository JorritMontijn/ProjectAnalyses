close all

%	strFile = '20140207M42xyt01.avi_video_heel';

%Load data
if ~exist('video_heel','var')
	strDir = 'D:\Data\Processed\video\';
	%strFile = '20140425M50_xyt01.avi_video_heel.mat';
	strFile = '20140430xyt08_video_heel.mat';
	%strFile = '20140314xyt01_video_heel.mat';
	load([strDir strFile]);
end

%Lopen van frame f1 tot f2 met stapjes van fx

%Detects first frame in the video of a non-black image
run eyetrackdetectstart
%f1 = 650;
f2 = length(video_heel.frames);
fx = 3;

%updaten van eyepick (kleur voor oog) om de 'fex' frames genomen over een bereik van
%'ex'frames.
fex = 1;
ex = 15;

%updaten van centerofmass (zwaartepunt ellipse) om de 'cmex' frames genomen over een bereik van
%'cmx'frames.
cmex = 1;
cmx = 3;

%preallocate vectors
frame = f1:fx:f2;
ecclipse_area = zeros(1,length(f1:f2));
dblCenterX = zeros(1,length(f1:f2));
dblCenterY = zeros(1,length(f1:f2));
eyeposxovertime= zeros(1,length(f1:f2));
eyeposyovertime= zeros(1,length(f1:f2));
semiminorovertime= zeros(1,length(f1:f2));
semimajorovertime= zeros(1,length(f1:f2));
eyekleurovertime = zeros(1,length(f1:f2));

[r, c, a] =size(video_heel.frames(1,1).cdata);

%User moet middenpunt van oog aangeven.
image(video_heel.frames(1,f1).cdata)
[x0,y0] = ginput(1)

x0 = round(x0);
y0 = round(y0);
xbase = x0;
ybase = y0;
contin = 0;

clearvars -except video_heel eyekleurpick rcvaluesovertime f1 f2 fx ex fex x0 y0 ecclipse_area xbase ybase contin strFile frame

close all
figure;
scale = 1;

uicontrol('Style', 'pushbutton', 'String', 'Continue to blink analysis',...
	'Position', [220 10 160 20],...
	'Callback', 'contin = 1;');

for f = f1:fx:f2
	
	if contin == 1
		break
	end
	
	if f == f1
		eyekleurpick_first = pickeyecolor((f+ex),ex,x0,y0,video_heel,f1);
		eyekleurpick_next = eyekleurpick_first;
		eyekleurpick = eyekleurpick_next;
		
	elseif f/fex == round(f/fex) && f>f1
		
		if isnan(x0) == 1
			x0 = xbase;
			y0 = ybase;
			
		elseif x0<2 && mean(dblCenterX)>1
			xtemp = dblCenterX;
			ytemp = dblCenterY;
			xtemp(xtemp<2)=[];
			ytemp(ytemp<2)=[];
			x0 = round(mean(xtemp));
			y0 = round(mean(ytemp));
		end
		
		eyekleurpick_next = pickeyecolor(f,ex,x0,y0,video_heel,f1);
		eyekleurovertime(f)=eyekleurpick_next;
		
	end
	
	if abs(eyekleurpick_next-eyekleurpick) < 0.25
		eyekleurpick = eyekleurpick_next;
	else
		eyekleurpick = ((eyekleurpick_first*2)+eyekleurpick+eyekleurpick_next)/4
	end
	
	tic
	
	img1 = im2double(video_heel.frames(1,scale*f).cdata);
	t = video_heel.times(1,scale*f);
	[r, c, a] =size(img1);
	
	%Average pixels with surrounding, new matrix called img2
	img2 = zeros(r,c,a);
	for m = 2:r-1
		for n = 2:c-1
			img2(m,n,:)=((img1(m,n,:)+img1(m,n+1,:)+img1(m,n-1,:)+img1(m+1,n,:)+img1(m+1,n+1,:)+img1(m+1,n-1,:)+img1(m-1,n,:)+img1(m-1,n+1,:)+img1(m-1,n-1,:))/9);
		end
	end
	
	%Average again, new matrix called img3
	img3=zeros(r,c,a);
	for o=3:r-2
		for p=3:c-3
			img3(o,p,:)=((img2(o,p,:)+img2(o,p+1,:)+img2(o,p-1,:)+img2(o+1,p,:)+img2(o+1,p+1,:)+img2(o+1,p-1,:)+img2(o-1,p,:)+img2(o-1,p+1,:)+img2(o-1,p-1,:))/9);
		end
	end
	
	img4=zeros(r,c,a);
	for o=3:r-2
		for p=3:c-3
			img4(o,p,:)=((img3(o,p,:)+img3(o,p+1,:)+img3(o,p-1,:)+img3(o+1,p,:)+img3(o+1,p+1,:)+img3(o+1,p-1,:)+img3(o-1,p,:)+img3(o-1,p+1,:)+img3(o-1,p-1,:))/9);
		end
	end
	
	img5(:,:,1) = imadjust(img4(:,:,1));
	img5(:,:,2) = imadjust(img4(:,:,1));
	img5(:,:,3) = imadjust(img4(:,:,1));
	
	eye = zeros(r,c,a);
	
	%neemt een kleinere ROI
	
	q1=round(x0-0.25*r);
	if q1 < 1
		q1 = 1;
	end
	q2=round(x0+0.25*r);
	if q2 > r
		q2 = r;
	end
	s1=round(y0-0.25*r);
	if s1 < 1
		s1 = 1;
	end
	s2=round(y0+0.25*r);
	if s2 > c
		s2 = c;
	end
	
	%vind pixels met soortgelijke intensiteit van aangegeven oog
	
	for q = q1:q2
		for s = s1:s2
			
			if img5(q,s,:)>= 0.65*eyekleurpick & img5(q,s,:)<= 1.75*eyekleurpick;
				eye(q,s,:)=1;
			else eye(q,s,:)=0;
			end
		end
	end
	
	
	%beeldbewerking van berekende stukken oog
	se = strel('ball',20,20);
	erodedI = imerode(eye,se);
	erodedII = imerode(erodedI,se);
	
	eroded = mat2gray(erodedII);
	eroded = ~bwareaopen(~eroded,2000);
	eroded = bwareaopen(eroded,900);
	
	eroded = eroded(:,:,1);
	eroded = bwmorph(eroded,'skel',3);
	eroded = bwmorph(eroded,'majority',15);
	
	eroded = bwareaopen(eroded,150);
	
	eroded(:,:,1) = eroded(:,:,1);
	eroded(:,:,2) = eroded(:,:,1);
	eroded(:,:,3) = eroded(:,:,1);
	
	
	%vind eventuele cirkels in gebied
	[centers,radii,metric] = imfindcircles(eroded(:,:,1),[10,30]);
	
	%Calling countblackspots: gives a matrix of [r c] points where intensity of
	%color is infinitive. (eye=1) rcvalues = countblackspots(eye);
	rcvalues = countblackspots(eroded(:,:,1));
	
	%create dummy matrix same size as eroded matrix on which initial
	%circular fit was done
	
	matX = eroded(:,:,1);
	matY = false(size(matX));
	for intVal=1:size(rcvalues,2);
		matY(rcvalues(2,intVal),rcvalues(1,intVal)) = true;
	end
	
	if  isempty(rcvalues) == 0
		
		%If there ARE rcvalues found, but no circles fit within, then take mean
		%ellipse center coordinates as new x and y coordinates for the eye.
		
		if isempty(radii) == 1
			%no circles
			
			ecclipse_area(f) = 1;
			semimajorovertime(f) = 1;
			semiminorovertime(f) = 1;
			dblCenterY(f) = 1;
			dblCenterX(f) = 1;
			
			if isnan(x0) == 1
				
				x0 = xbase;
				y0 = ybase;
				
			end
			
			if sum(dblCenterX) > length(f1:fx:f)
				
				xtemp = dblCenterX;
				ytemp = dblCenterY;
				xtemp(xtemp<2)=[];
				ytemp(ytemp<2)=[];
				
				x0 = round(mean(xtemp));
				y0 = round(mean(ytemp));
				
			end
		end
		
		rcdist = (std(rcvalues(1,:)))+(std(rcvalues(2,:))) / 2;
		
		del1 = find(rcvalues(1,:)<x0-25);
		del2 = find(rcvalues(1,:)>x0+25);
		del3 = find(rcvalues(2,:)<y0-25);
		del4 = find(rcvalues(2,:)>y0+25);
		del=unique([del1,del2,del3,del4]);
		rcvalues(:,del) = [];
		
		if  isempty(centers) == 0
			
			del1 = find(rcvalues(1,:)<centers(1)-1.5*radii(1));
			del2 = find(rcvalues(1,:)>centers(1)+1.5*radii(1));
			del3 = find(rcvalues(2,:)<centers(2)-1.5*radii(1));
			del4 = find(rcvalues(2,:)>centers(2)+1.5*radii(1));
			del=unique([del1,del2,del3,del4]);
			rcvalues(:,del) = [];
		else
		end
		
		if rcvalues > 300
			
			del1 = find(rcvalues(1,:)<mean(rcvalues(1,:))-2*std(rcvalues(1,:)));
			del2 = find(rcvalues(1,:)>mean(rcvalues(1,:))+2*std(rcvalues(1,:)));
			del3 = find(rcvalues(2,:)<mean(rcvalues(2,:))-2*std(rcvalues(2,:)));
			del4 = find(rcvalues(2,:)>mean(rcvalues(2,:))+2*std(rcvalues(2,:)));
			del=unique([del1,del2,del3,del4]);
			rcvalues(:,del) = [];
			
		else
		end
		
		steps = 36;
		eyeposxovertime(f) = x0;
		eyeposyovertime(f) = y0;
		
		if  isempty(radii) == 0
			%RCvalues found AND circle fit
			
			%set initial parameters
			vecParamsInit = [y0 x0 -(1/8)*pi radii(1) radii(1)];
			
			%perform fitting
			vecParamsFit = MLFit('ellipseArea', vecParamsInit, matX, matY, 1); % no mse
			
			%get parameters
			
			dblCenterX(f) = vecParamsFit(1);
			dblCenterY(f) = vecParamsFit(2);
			
			dblMajorAxisAngle = vecParamsFit(3);
			dblMajorRadius = vecParamsFit(4);
			dblMinorRadius = vecParamsFit(5);
			
			%get fitted circle
			matFitY = ellipseArea(matX,vecParamsFit);
			
			%put in output
			semimajorovertime(f) = dblMajorRadius;
			%semimajorovertime(semimajorovertime==0)=[];
			
			semiminorovertime(f) = dblMinorRadius;
			%semiminorovertime(semiminorovertime==0)=[];
			
			ecclipse_area(f) = sum(matFitY(:));
			%ecclipse_area(ecclipse_area==0)=[];
			
			x0 = dblCenterX(f);
			y0 = dblCenterY(f);
			
		end
		
	else
		%If there are no rcvalues found, then
		eyeposxovertime(f) = xbase;
		eyeposyovertime(f) = ybase;
		x0 = round(( 2*eyeposxovertime(f) + eyeposxovertime(f-(1*fx) ) + eyeposxovertime(f-(2*fx) ) ) / 4 );
		y0 = round(( 2*eyeposyovertime(f) + eyeposyovertime(f-(1*fx) ) + eyeposyovertime(f-(2*fx) ) ) / 4 );
		
	end
	
	%Region of interest for plotted elipse is taken as mean of last 3 frames.
	if f>=f1+(3*fx)
		if eyeposxovertime(f-(1*fx))>0 && eyeposxovertime(f-(2*fx))>0
			x0 = round(( 2*eyeposxovertime(f) + eyeposxovertime(f-(1*fx) ) + eyeposxovertime(f-(2*fx) ) ) / 4 );
			y0 = round(( 2*eyeposyovertime(f) + eyeposyovertime(f-(1*fx) ) + eyeposyovertime(f-(2*fx) ) ) / 4 );
		end
	end
	
	if f<f1+(10*fx)
		
		if f>f1
			if exist('texti','var')
				delete(texti)
			end
		else
		end
		
		image1 = image(img1);
		hold on
		
		if isempty(radii) == 0
			[row,col] = find(matY);
			plot(row,col,'ro');
			%   plot(abs(x0),abs(y0),'*')
			if isempty(rcvalues) == 0
				[row,col] = find(matFitY);
				plot(row,col,'bx');
			else
			end
			axis([1 c 1 r])
			set(gca,'YDir','reverse');
			
		else
			
		end
		
		texti = text(3*pi/4,-1,['time = ',int2str(t) ' frame = ' int2str(f)],'FontSize',16);
		drawnow
		
	else
	end
	
	eyekleurpick
	toc
	
end

clearvars -except video_heel eyekleurpick rcvaluesovertime ecclipse_area xbase ybase strFile eyeposxovertime eyeposyovertime frame dblCenterX dblCenterY
eyetracking.ecclipse_area = ecclipse_area;
eyetracking.frame = frame;
eyetracking.dblCenterX = dblCenterX;
eyetracking.dblCenterY = dblCenterY;
savefile = [strFile(1:(length(strFile)-11)),'_eyetracking.mat'];
current = pwd;
cd('D:\Data\Processed\video\');
save(savefile, 'eyetracking');
disp('Eyetracking saved.');
cd(current);
waitforbuttonpress
run blinkdet