% Script to perform super-resolution of spikes from low frequency
addpath ./matlab2tikz-matlab2tikz-816f875/src
clear all;
close all;

% Physical constants
responsivity = 0.7;
electron_charge = 1.6e-19;

% Receiving power
% ary = linspace(0, 2.5, 11);
% ary = 1.7;
received_photons = linspace(0, 41, 20);
% received_photons = 50;
arraylength = length(received_photons);
% received_photons = 10.^ary;

% FFT frame creation
repetition = 10000;
frequency = 15.5;
bandwidth = 100;
sampling_rate = 2*bandwidth;
sampling_period = 1/sampling_rate;
initial_phase = 0;
number_of_points = 200;
observation_time = sampling_period * (number_of_points);
frequency_resolution = 1/observation_time;

% Create input sample sinusoidal tone
index = linspace(0, number_of_points-1, number_of_points);
t = sampling_period * index;
f = frequency_resolution * index;
Ps = electron_charge * received_photons / (responsivity * observation_time);
amplitude = 2*responsivity*sqrt(Ps);
x = amplitude' * sin(2*pi*frequency*t + initial_phase);
x = repmat(x,repetition,1);

% Additive Gaussian noise creation
noise_floor = 2 * electron_charge * responsivity;
vari = (sampling_rate / 2) * noise_floor;
sigma = sqrt(vari);
noise = sigma * randn(repetition*arraylength,number_of_points);
SNR = 20*log10((amplitude/sqrt(2))/sigma);

% Taking FFT
fx = zeros(arraylength, repetition, number_of_points);
for i=1:repetition
    for j=1:arraylength
        fx(j,i,:) = fft(x((i-1)*arraylength+j,:)+noise((i-1)*arraylength+j,:), number_of_points);
    end
end

% Take MUSIC
result_music = zeros(arraylength, repetition);
for i=1:repetition
    for j=1:arraylength
        [s, w] = pmusic(x((i-1)*arraylength+j,:)+noise((i-1)*arraylength+j,:),2,4096);
        [pwr, idx] = max(s);
        result_music(j,i) = w(idx)*100/pi;
    end
end

result_fx = 2*abs(fx(:,:,1:end/2)/(number_of_points)).^2;
[maxpower, detection_result] = max(result_fx,[],3);
detection_result = detection_result - 1;

RMSE_fft = rms(detection_result' - frequency);
RMSE_music = rms(result_music' - frequency);

single_spectrum = result_fx(:,1,:);
averaged_spectrum = mean(result_fx(:,:,:),2);

sgl = zeros(1,bandwidth);
avge = zeros(1,bandwidth);

for i=1:bandwidth
    sgl(i) = single_spectrum(1,1,i);
    avge(i) = averaged_spectrum(1,1,i);
end
% 
% figure(1);
% plot(f(1:end/2),db(sgl,'power'));
% set(gca,'yscale','log');
% 
% figure(2);
% plot(f(1:end/2),db(avge,'power'));
% set(gca,'yscale','log');
% 
% figure(3);
% % scatter( ones(1,length(detection_result)) ,detection_result);
% histogram(detection_result, 'BinWidth', 1, 'BinLimits', [-0.5, 99.5]);
% 
% figure(4);
% % scatter( ones(1,length(detection_result)) ,detection_result);
% histogram(result_music, 'BinWidth', 1, 'BinLimits', [-0.5, 99.5]);


figure(1);
plot(received_photons, RMSE_fft,'LineWidth',2);
hold on;
plot(received_photons, RMSE_music,'LineWidth',2);
title('Detection Error for DFT/MUSIC (Low SNR)','FontSize',18,'FontWeight','normal');
xlabel('SNR','FontSize',15);
ylabel('Error (Hz)','FontSize',15);
legend('DFT','MUSIC');
grid on;
set(gca,'xlim',[0 41])
set(gca,'xtick',[0:10:40])
set(gca,'ylim',[0 50])
set(gca,'ytick',[0:10:50])
matlab2tikz('./final_report/figure2a.tex');
