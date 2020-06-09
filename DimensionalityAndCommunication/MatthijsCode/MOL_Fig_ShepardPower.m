
Par.strAuStimType               = 'shepardtone';
Par.intSamplingRate             = 192000;                    %Preferred sampling rate
Par.FreqSpacing                 = 1/256;                    %Spacing between consecutive frequencies in octaves
Par.allOctaves                  = 13:Par.FreqSpacing:14-Par.FreqSpacing; %Init all octaves
Par.vecLoadFreqs                = 2.^Par.allOctaves;        %Init all frequencies
Par.vecLoadFreqs                = round(Par.vecLoadFreqs);  %Round them off to closest integer
% Par.vecLoadFreqs                = 8000:10:16000-10;  %Round them off to closest integer

Par.nFreqs                      = length(Par.vecLoadFreqs); %store number of loaded frequencies for module circularity
Par.ShepardTones                = 7;                        %Number of shepard tones (partial tones above and below)
Par.ShepardWeights              = gausswin(length(Par.vecLoadFreqs) * Par.ShepardTones) / 2.5; %Generate partial tone weights divide for sound range 0-1
Par.vecFreqChange               = 0.5;                      %Vector with possible delta frequency in octave
Par.minSampleDur                = 0.02;                     %Minimum duration of loaded tone that is repeated (too low and the schedule might get empty)
Par.maxSampleDur                = 0.1;                      %Maximum duration of loaded tone that is repeated (this determines resolution for changing ad hoc)

intThisTrial                    = 1;
Par.Stim.vecAuStimInt(intThisTrial) = 1;
nSamples                        = 100000;

psdx_all        = NaN(Par.nFreqs,nSamples/2+1);
SPL             = NaN(Par.nFreqs,1);
for iFreq = 1:Par.nFreqs
    
    Par.Stim.vecFreq(intThisTrial)    = Par.vecLoadFreqs(iFreq);
    varaudio            = load_auditory(Par,intThisTrial);
    sound               = repmat(varaudio.vecSound(1,:),1,5000);
    sound               = sound(1:nSamples);
%     sound               = filterA(sound,Par.intSamplingRate); 
    SPL(iFreq)          = spl(sound,'air');
    x                   = sound + randn(size(sound))/10000; %Add noise to facilitate log10 imagesc
    N                   = length(x);
    xdft                = fft(x);
    xdft                = xdft(1:N/2+1);
    psdx                = (1/(Par.intSamplingRate*N)) * abs(xdft).^2;
    psdx(2:end-1)       = 2*psdx(2:end-1);
%     psdx                = smooth(psdx,10,'sgolay',2)';
    psdx                = smooth(psdx,10)';
    psdx_all(iFreq,:)   = psdx;
end

freq                    = 0:Par.intSamplingRate/nSamples:Par.intSamplingRate/2;

figure; set(gcf,'units','normalized','Position',[0.2 0.1 0.5 0.7],'color','w')
imagesc(10*log10(psdx_all)); hold on;
% imagesc(psdx_all); hold on;

[~,lowIdx]              = min(abs(freq-Par.vecLoadFreqs(1)));
[~,highIdx]             = min(abs(freq-Par.vecLoadFreqs(end)));
plot([lowIdx highIdx highIdx lowIdx lowIdx],[Par.nFreqs Par.nFreqs 1 1 Par.nFreqs],'k','LineWidth',2);

title('Power Spectrum - 8-16kHz with 10 Hz changes, 48 kHz samplerate, 9 shepards','FontSize',15)
title('Power Spectrum - 8.192-16.384kHz with 1/256 octave changes, 192 kHz samplerate, 5 shepards','FontSize',15)

set(gca, 'YTick', 1:Par.nFreqs/10:Par.nFreqs, 'YTickLabels',Par.vecLoadFreqs(1:Par.nFreqs/10:Par.nFreqs)/1000,'FontSize', 15)
set(gca, 'XTick', 1:nSamples/25:length(freq), 'XTickLabels',freq(1:nSamples/25:length(freq))/1000,'FontSize', 12)
xlabel('Frequency (kHz)')
ylabel('Ground Tone (kHz)')

figure; set(gcf,'units','normalized','Position',[0.2 0.1 0.5 0.7],'color','w')
plot(Par.vecLoadFreqs,SPL,'b','LineWidth',3);
xlabel('Ground frequency (kHz)')

