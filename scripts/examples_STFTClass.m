clear all
% Include path to STFTClass.m
addpath(genpath('library'))

% Let say wee have a audio signal 2s long
fs = 16000;
x = rand(fs * 2, 1);

%% STFT Paramters
winlen = uint32(256); % it means 256 samples. 
% winlen = 256; % It means 256 ms

hop = 0.25;  % 75% overlap. Default is 50%, or 0.5

nfft = 512; % Default is same length as winlen


%% Get STFTClass Object
stftObj = STFTClass(fs, winlen, hop, nfft);

% Asking NOT to display information while processing. It's optional 
stftObj.setVerbose(false);


%% Perform STFT on x
X = stftObj.stft(x);

% We can also get first T(user-defined) time frames using
% X = stftObj.stft(x, T);

% And we can also plot the spectrogram with
% X = stftObj.stft(x, T, true);


%% Get STFT output parameters

% Get time (in sec) for each time frames
TArr = stftObj.getTimeBins();

% Get center frequency for each frequency bin
FArr = stftObj.getFrequencyBins();

% Get frequencies for certain STFT bins
freqBinMatlab = [4, 20, 49];
freq = stftObj.getFrequency(freqBinMatlab);

% Get frequency bin number for certain frequencies
freq_arr = [1200, 3900, 5000];
freqBins = stftObj.getFreqBinMatlab(freq_arr);

% Get STFT setup parameters
[windowsize, hopsize, nfft] = stftObj.getParams();

% Print STFT setup parameters
stftObj.printParams();

% Plot spectrogram
isBinned = false; % True, if you want to display bin number instead of time and frequency
stftObj.plotSpectrogram(X, isBinned)



%% Get back time domain signal from STFT
x1 = stftObj.istft(X);

% We can keep signal length to sigL(user-defined) by
% x = stftObj.istft(X, sigL);
