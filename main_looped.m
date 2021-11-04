close all
clear all
rng(4599136)

start_path = fullfile(matlabroot, '\toolbox');
topLevelFolder = 'C:\Users\Wervers\Desktop\TIMIT\Speech'; %Change to current directory
filePattern = sprintf('%s/**/*.wav', topLevelFolder);
allFileInfo = dir(filePattern);

%Extract the speech samples from the dataset
dataset = zeros(floor(length(allFileInfo)/20),121140);
for j = 1:floor(length(allFileInfo)/20)
    i = 2*j                            %Due to weird format in dataset, we need to skip odd numbers           
    part1 = allFileInfo(i).folder;
    part2 = allFileInfo(i).name;
    target = append(part1,'\', part2);
    addvec = audioread(target)';
    
    addvec = sqrt(length(addvec)) * addvec/ sqrt(sum(addvec.^2));
    dataset(j,1:length(addvec)) = addvec;
    lengths(j) = length(addvec);
end


SNR_target = -3; %The desired SNR after noise is added

%Load in the babble noise
[babble, fs] = audioread('babble noise.wav');

%Create empty vectors for the results
SNR_TSNR = nan(1,size(dataset,1));
SNR_block = nan(1,size(dataset,1));
SNR_sure = nan(1,size(dataset,1));
SNR_kalman = nan(1,size(dataset,1));
SNR_TSNR_babble = nan(1,size(dataset,1));
SNR_block_babble = nan(1,size(dataset,1));
SNR_sure_babble = nan(1,size(dataset,1));
SNR_kalman_babble = nan(1,size(dataset,1));
%Loop through the main for each of the selected speech samples
parfor i = 1:size(dataset,1)
i
%Get data and contaminate with white noise    
noisetarget = randn(size(dataset(i,1:lengths(i))))*std(dataset(i,1:lengths(i)))/db2mag(SNR_target);
noisysample = dataset(i,1:lengths(i)) + noisetarget;

%Scale babble noise for the desired SNR
babbletarget = (babble(1:lengths(i))/std(babble(1:lengths(i))))*std(dataset(i,1:lengths(i)))/db2mag(SNR_target);
babblesample = dataset(i,1:lengths(i)) + babbletarget';


%perform enhancement algorithms
%TSNR
outTSNR= TSNRWiener(noisysample',16000,1000);
outTSNR =  (outTSNR/sqrt(sum(outTSNR.^2))) * sqrt(length(outTSNR)); %scale the power of the signal the same way as the input signal
outTSNR_babble = TSNRWiener(babblesample',16000,1000);
outTSNR_babble =  (outTSNR_babble/sqrt(sum(outTSNR_babble.^2))) * sqrt(length(outTSNR_babble)); %scale the power of the signal the same way as the input signal

%Block James-Stein shrinkage
out_block = wdenoise(noisysample,'DenoisingMethod','BlockJS');
out_block =  (out_block/sqrt(sum(out_block.^2))) * sqrt(length(out_block)); %scale the power of the signal the same way as the input signal
out_block_babble = wdenoise(babblesample,'DenoisingMethod','BlockJS');
out_block_babble =  (out_block_babble/sqrt(sum(out_block_babble.^2))) * sqrt(length(out_block_babble)); %scale the power of the signal the same way as the input signal

%SURE wavelet shrinkage
out_sure = wdenoise(noisysample,'DenoisingMethod','SURE', 'ThresholdRule','soft');
out_sure = (out_sure/sqrt(sum(out_sure.^2))) * sqrt(length(out_sure)); %scale the power of the signal the same way as the input signal
out_sure_babble = wdenoise(babblesample,'DenoisingMethod','SURE', 'ThresholdRule','soft');
out_sure_babble = (out_sure_babble/sqrt(sum(out_sure_babble.^2))) * sqrt(length(out_sure_babble)); %scale the power of the signal the same way as the input signal

%Kalman filtering
outkalman = kalman_speechV2(noisysample,16000,dataset(i,1:lengths(i)), std(dataset(i,1:lengths(i)))/db2mag(SNR_target));
outkalman =  (outkalman/sqrt(sum(outkalman.^2))) * sqrt(length(outkalman)); %scale the power of the signal the same way as the input signal
outkalman_babble = kalman_speechV2(babblesample,16000,dataset(i,1:lengths(i)), std(dataset(i,1:lengths(i)))/db2mag(SNR_target));
outkalman_babble =  (outkalman_babble/sqrt(sum(outkalman_babble.^2))) * sqrt(length(outkalman_babble)); %scale the power of the signal the same way as the input signal


%get SNRs
%SNR_babble(i) = snr(dataset(i,1:lengths(i)), babbletarget');
SNR_TSNR(i) = snr(dataset(i,1:lengths(i)),dataset(i,1:lengths(i)) - outTSNR');
SNR_block(i) = snr(dataset(i,1:lengths(i)),dataset(i,1:lengths(i)) - out_block);
SNR_sure(i) = snr(dataset(i,1:lengths(i)),dataset(i,1:lengths(i)) - out_sure);
SNR_kalman(i) = snr(dataset(i,1:lengths(i)),dataset(i,1:lengths(i)) - outkalman);

SNR_TSNR_babble(i) = snr(dataset(i,1:lengths(i)),dataset(i,1:lengths(i)) - outTSNR_babble');
SNR_block_babble(i) = snr(dataset(i,1:lengths(i)),dataset(i,1:lengths(i)) - out_block_babble);
SNR_sure_babble(i) = snr(dataset(i,1:lengths(i)),dataset(i,1:lengths(i)) - out_sure_babble);
SNR_kalman_babble(i) = snr(dataset(i,1:lengths(i)),dataset(i,1:lengths(i)) - outkalman_babble);

end

%Compute the avarages
avg_white_TSNR = mean(SNR_TSNR);
avg_white_kalman = mean(SNR_kalman);
avg_white_block = mean(SNR_block);
avg_white_sure = mean(SNR_sure);

avg_babble_TSNR = mean(SNR_TSNR_babble);
avg_babble_kalman = mean(SNR_kalman_babble);
avg_babble_block = mean(SNR_block_babble);
avg_babble_sure = mean(SNR_sure_babble);
