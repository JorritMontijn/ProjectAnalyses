function [vecFrequency,vecAmplitude] = getFFT(vecSignal,dblSamplingFreq)
	%perform FFT
	intL = length(vecSignal);
	NFFT = 2^nextpow2(intL); % Next power of 2 from length of signal
	vecAmplitude = fft(vecSignal,NFFT)/intL;
	vecFrequency = dblSamplingFreq/2*linspace(0,1,NFFT/2+1);
	
	%{
% Plot single-sided amplitude spectrum.
plot(f,2*abs(Y(1:NFFT/2+1)))
title('Single-Sided Amplitude Spectrum of y(t)')
xlabel('Frequency (Hz)')
ylabel('|Y(f)|')
	%}
end