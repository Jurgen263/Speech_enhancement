function [s] = kalman_speechV2(x,Fs,clean,Qv)
%UNTITLED7 Summary of this function goes here
%   Input x  : noisy speech
%   Input Fs : sampling rate 
%   Output s : speech after kalman filtering

%---------------------settings----------------------
wt = 25; %window time in ms
order = 10;


wl = fix((wt/1000)*Fs); %length of each section in samples

y2 = buffer(x,wl);
y = buffer(x,wl,wl*0.5,'nodelay'); %split input series into sections of the window length
y_clean = buffer(clean,wl,wl*0.5,'nodelay');


%Qv = var(y(:,1));                   %Estimate Qv as the variance of the first section
x_filtered = nan(size(y));
for i = 1:size(y,2)
    section = y(:,i).*sqrt(hanning(wl));  %take section
    clean_section = y_clean(:,i).*sqrt(hanning(wl));
    corr = xcorr(clean_section);              %determine correlation sequence
    corr(1:length(clean_section)-1) = [];  %make sure, sequence starts at index 0
    [coefs,Qw] = levinson(corr,order);  %determine coefficients and Qw
    A = [zeros(order-1,1) eye(order-1);fliplr(-coefs(2:end))];  %create A matrix
    C=[zeros(1,order-1),1];                                     %create C vector
                                                                %determine Qv
    
    x_filtered(:,i) = kalman_filter(section,A,C,Qw,Qv);
end

s = zeros(1,size(y2,1)*size(y2,2)); %bit hardcoded for overlap 0.5
s(1:wl) = x_filtered(:,1);
for i = 1:size(y,2)-1
    s((0.5+0.5*(i-1))*wl + 1:(1.5+0.5*(i-1))*wl) = s((0.5+0.5*(i-1))*wl + 1:(1.5+0.5*(i-1))*wl) + x_filtered(:,i+1)';    
end
s = s(1:length(x));
end

