function [x] = TSNRWiener(y,fs,LT)
%	Input y  :   noisy speech
%   Input fs :   sampling frequency
%   Input LT :   listening time
%   Output x :   filtered output speech signal

%% Initialization
l = length(y);
wt = 25; %window time in ms
wl = fix((wt/1000)*fs); %length of each segment in samples
NFFT = 2*wl; %length of FFT
win = hanning(wl); 

%% Computation of FFT noise magnitude
N = zeros(NFFT,1);
for i = 0:(LT-wl)
    Nwindowed = y(i+1:i+wl).*win;
    N = N + abs(fft(Nwindowed,NFFT)).^2;
end
N = N / (LT - wl + 1); %estimated spectral noise magnitude
%% Variable declarations
oldXmag = zeros(NFFT,1);
newX = zeros(l,1);
overlap = fix(0.5*wl); %50% overlap for windowing
offset = wl - overlap;
max_it = fix((l-NFFT)/offset);

beta = 0.98;
epsilon1 = 0.15;
epsilon2 = 0.26;

%% Main algorithm
for i = 0:max_it
    itbegin = i*offset+1;
    itend = i*offset+wl;
    ffty = fft(win.*y(itbegin:itend),NFFT); %FFT of windowed speech segment
    ymag = abs(ffty); %spectral magnitude of speech
    yphase = unwrap(angle(ffty)); %spectral phase of speech
    
    %Decision-Directed step
    SNRpost = max(((ymag.^2)./N) - 1,epsilon1); %a posteriori SNR
    SNRddprior = beta*(oldXmag.^2)./N + (1-beta)*SNRpost; %DD a priori SNR
    Gdd = SNRddprior./(1+SNRddprior); %DD Gain
    newXmag = Gdd .* ymag; %Spectral magnitude estimate of clean speech
    
    %TSNR step 2
    SNRtsnrprior = (newXmag.^2)./N; %TSNR step 2 a priori SNR
    Gtsnr = max(SNRtsnrprior./(1+SNRtsnrprior),epsilon2); %TSNR step 2 Gain
    newXmag = Gtsnr.*ymag; 
    ffty = newXmag.*exp(1i*yphase); %Updated spectral estimate of speech signal
    oldXmag = abs(newXmag);
    
    %Add the estimated segment to the past values
    newX(itbegin:itbegin+NFFT-1) = newX(itbegin:itbegin+NFFT-1) + 0.5*real(ifft(ffty,NFFT));
end
x = newX; %TSNR estimate of time domain speech signal
x = x * max(abs(y))/max(abs(x));
end
