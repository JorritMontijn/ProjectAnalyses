function [iw,ahp,pr,ppd,slope,had]=waveforms_features(spike,fs)
% Function for finding features of spike waveforms
%
%Inputs:
%spike: average spike waveform for one specific cluster/neuron
%fs: sampling frequency
%N.B. Spikes shall be aligned by the positive (highest) peak, the opposite
%of both [Niell & Stryker 2008] and [Buzsaki 2004], but nonetheless the way
%patch clampers (and positive-thinking people) like spikes
%
%Extracted features are Initial Width (iw) and AfterHyperPolarization
%(ahp), following [Bruno & Simons 2002] in us
%Also, following [Niell & Stryker 2008], we extract peak ratio (pr), peak-to-peak delay in ms (ppd)
%and slope 500 ms after the first peak (slope)
%Moreover, we extract half-amplitude duration (had), following [Buzsaki 2004]
%
%--------------------------------------------------------------------------

%Upsampling
interp_factor=10;
spike=interp(spike,interp_factor);
[value, max_point]=max(spike(1:(ceil(length(spike)*2/3))));
fs=fs*interp_factor;

iw=0;
ahp=0;
pr=0;
ppd=0;
slope=0;
had=0;

%Look for iw before peak
found=0;
i=max_point;
while ~found && i>0
    if spike(i)>0
        iw=iw+1;
        i=i-1;
    else
        found=1;
    end
end

%Look for iw after peak
found=0;
i=max_point+1;
while ~found && i<(length(spike)+1)
    if spike(i)>0
        iw=iw+1;
        i=i+1;
    else
        found=1;
    end
end

%Look for ahp after peak
found=0;
i=i;
while ~found && i<(length(spike)+1)
    if spike(i)<0
        ahp=ahp+1;
        i=i+1;
    else
        found=1;
    end
end

iw=iw/fs*1e6;
ahp=ahp/fs*1e6;

%Look for hwa
hw=spike(max_point)/2; %Half width
[app pre_]=min(abs(spike(1:max_point)-hw));
pre_=max_point-pre_;
[app post_]=min(abs(spike(max_point+1:end)-hw));
had=(post_+pre_)/fs*1000;

%Look for pr, ppd and slope
%[v_max, t_max]=max(spike);
t_max=max_point;
v_max=spike(max_point);

[v_min, t_min]=min(spike(max_point+1:end));%(1:max_point));
%[v_min, t_min]=min(spike);

pr=abs(v_min/v_max);
%ppd=abs(t_min-t_max)/fs*1000;
ppd=t_min/fs*1000;
spike=spike/v_max;
app=gradient(spike);
% if t_min>t_max
    slope=app(round(max_point+0.5e-3*fs));
% else
%     slope=app(round(max_point-0.5e-3*fs));
% end