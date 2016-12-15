addpath ./matlab2tikz-matlab2tikz-816f875/src

clear all;
close all;

received_photons = linspace(1, 41, 21);

data_fft = load('./saved_data/fft.mat');
data_music = load('./saved_data/music.mat');
data_l1 = load('./saved_data/l1.mat');
data_l2 = load('./saved_data/l2.mat');

RMSE_fft    = data_fft.RMSE_fft;
RMSE_music  = data_music.RMSE_music;
RMSE_l1     = data_l1.RMSE_l1;
RMSE_l2     = data_l2.RMSE_l2;

figure(1);
plot(received_photons, RMSE_fft,'LineWidth',2);
hold on;
plot(received_photons, RMSE_music,'LineWidth',2);
hold on;
plot(received_photons, RMSE_l1,'LineWidth',2);
hold on;
plot(received_photons, RMSE_l2,'LineWidth',2);
title('Detection Error for DFT/MUSIC/L1/L2 (Low SNR)','FontSize',18,'FontWeight','normal');
xlabel('SNR','FontSize',15);
ylabel('Error (Hz)','FontSize',15);
legend('DFT','MUSIC','L1','L2');
grid on;
set(gca,'xlim',[0 42])
set(gca,'xtick',[0:10:40])
set(gca,'ylim',[0 50])
set(gca,'ytick',[0:10:50])
matlab2tikz('./plots/comparison.tex');
