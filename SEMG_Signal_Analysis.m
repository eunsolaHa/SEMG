clc;close all;clear;

%% SEMG signal Time domain and frequency domain analysis
    % Time domain : RMS,IEMG and Zero-crossing rate
    % Frequency domain : Mean frequency and Median frequency 
    % Time-Frequency domain : Mean frequency and Median frequency 
    
%% load data form one file
data = load('°ûÈñ¼· ÃÖ´ë±Ù·Â.txt');
totalSamples = size(data,1); % size : returns the number of rows of data
%totalSmples meams total number of rows of stored data
 
xin = data(2000:totalSamples-2000,1); % exclude data at the beginning and end
fs = 1/0.001; %(=1000Hz) indicates 1000 samples per second

%visualization of Raw data
plot(xin)
xlabel('Time (ms)'); 
ylabel('amplitude')
title('SEMG Raw data');

%% 4 order band-pass filter
fmin = 20;
fmax = 450;
N=4;
[b a] = butter(N, [fmin, fmax]/(fs/2), 'bandpass'); % 
[H W] = freqz(b,a,1000);  
BBF_xin = filter(b,a, xin);

% Calculate frequency axis in Hz
w = (0:999) * (fs / 2) / 1000;

% visualization of bandpass filter
figure;
plot(w,abs(H));
axis([0 fs/2 0 1.2]); %Set x-axis to span from 0 to fs/2
xlabel('frequency');
ylabel('amplitude')
title('4th order bandpass Filter');

figure;
p=plot(BBF_xin,'b');
set(p,'LineWidth',0.1);
xlabel('frequency');
ylabel('amplitude')
title('Filtered EMG');

%% rolling window : Divide data by time to identify characteristics in small sections

l = length(BBF_xin);
window = 500;
r = floor(l/window);
z = zeros(window, r);

for ii=1:r
    start = int32(500*(ii-1)+1); % 
    finish = int32(500*ii);
    z(:,ii)=BBF_xin(start:finish,1);
end

%% Time domain
% RMS(Root Mean Square) : Indicates the size of force through changes in amplitude

rms_c = zeros(r, 1);
for iii = 1:r
    rms_c(iii, 1) = sqrt(mean(z(:, iii).^2)); 
    % rms_c(iii,1) = rms(z(:,iii))
end

% RMS visualization
figure; plot(rms_c)
xlabel('Time(=rolling window)');
ylabel('Amplitude');
title('Root Mean Square(RMS)');

% integrated EMG : Muscle contractility and motor unit recruitment ability

iemg = zeros(r, 1);
for iiii = 1:r
    rectified_signal = abs(z(:, iiii));
    integrated_signal = cumsum(rectified_signal); % z(:, 1) 

    iemg(iiii,1) = integrated_signal(end); 
end

% IEMG visualization 
figure; 
plot(iemg)
xlabel('Time');
ylabel('Amplitude');
title('integrated EMG(IEMG)');

% zero crossing : The ratio between positive and negative numbers within a signal 
                % and represents frequency characteristics.

zc_r = zeros(r, 1);
zc_c = zeros(r, 1);

% use function of zerocrossrate
%for iiiii = 1:r
    %[rate, count] = zerocrossrate(z(:,iiiii));
    %zc_r(iiiii, 1) = rate;
    %zc_c(iiiii, 1) = count;
%end

% calcultation of zerocrossrate
for iiiii = 1:r
    zero_crossings = 0;
    for i = 2:length(z(:, iiiii))
        if (z(i, iiiii) >= 0 && z(i - 1, iiiii) < 0) || (z(i, iiiii) < 0 && z(i - 1, iiiii) >= 0)
            zero_crossings = zero_crossings + 1;
        end
    end
    zc_r(iiiii) = zero_crossings / (length(z(:, iiiii)) - 1);
end

% zero-crossing rate visualization 
figure;
plot(zc_r);
xlabel('Rolling Window(n=500)');
ylabel('Zero-Crossing Rate');
title('Zero-Crossing Rate');

%% Frequency domain 

% Fast Fourier trasnsform(fft)
emg1 = z(:,1);
A = fft(emg1); 
n = length(emg1);
p = abs(A); 
p = p(1:n/2);
p(1)=0; 

% FFT visualization 
figure;
plot(p);
xlabel('Frequency');
ylabel('Amplitude');
title('FFT');

% When drawing two graphs
%subplot(2,1,1); % First subplot of 2x1 grid
%plot(emg1);
%xlabel('Time');
%ylabel('Amplitude');
%title('Raw data');

%subplot(2,1,2); % Second subplot of 2x1 grid
%plot(p);
%xlabel('Frequency');
%ylabel('Amplitude');
%title('Frequency Domain (FFT)');


% Mean frequency(MNF) and Median frequency(MDF) through FFT 

mnf = zeros(r, 1);
mdf = zeros(r, 1);
for iiiiii = 1:r
    emg_s = z(:,iiiiii);
    emg_l = length(emg_s);
    emg_fft = abs(fft(emg_s));
    mnf(iiiiii, 1) = meanfreq(emg_fft, 500);
    
    mdf(iiiiii, 1) = medfreq(emg_fft, 500);
end

figure; plot(mnf)
xlabel('Time (s)');
ylabel('Frequency(Hz)');
title('MNF Analysis in Frequency');

figure; plot(mdf)
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('MDF Analysis in Frequency');

%% Time-Frequency domain
% STFT(Short Time fourier transform)
fs = 1000;
windowLength = 500; 
Noverlap = 0; % no overlap
gaussian = gausswin(windowLength); % Improved frequency resolution
[s, f, t] = spectrogram(BBF_xin, gaussian, Noverlap, windowLength, fs);
           %spectrogram : Indicates frequency change over time

% STFT visualization 
figure;
imagesc(t, f, abs(s));
axis xy;
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('STFT');

% MNF visualization
mnf_f = meanfreq(abs(s),f);
figure;plot(mnf_f)
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('MNF Analysis in Time-Frequency');

% MDF visualization
mdf_f = medfreq(abs(s), f);
figure;plot(mdf_f)
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('MDF Analysis in Time-Frequency');