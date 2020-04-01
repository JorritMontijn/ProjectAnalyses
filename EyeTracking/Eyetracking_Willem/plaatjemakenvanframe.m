%Lopen van frame f1 tot f2 met stapjes van fx
f1 = 29110;
f2 = 29140;
fx = 1;
%updaten van eyepick om de 'fex' frames genomen over een bereik van
%'ex'frames.
fex = 1;
ex = 1000;

%close all
figure
scale = 1;

%for f = 1:(length(video_heel.frames)/scale)
for f = f1:fx:f2
    tic
    
    img1 = im2double(video_heel.frames(1,scale*f).cdata);
    t = video_heel.times(1,scale*f);
    [r c a] =size(img1);
    
    if f>f1
        delete(texti)
    else
    end
    
    image1 = image(img1);
    hold on

    texti = text(3*pi/4,-1,['time = ',int2str(t) ' frame = ' int2str(f)],'FontSize',16);
    drawnow
    
    %pause(1.5)
    toc
    
end