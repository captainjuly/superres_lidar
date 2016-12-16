    %% Script to perform super-resolution of spikes from low frequency
addpath ./matlab2tikz-matlab2tikz-816f875/src
clear all;
close all;
%

%% Physical constants
responsivity = 0.7;
electron_charge = 1.6e-19;
iterations = 41;
%

%% Receiving power
% ary = linspace(0, 2.5, 11);
% ary = 1.7;
received_photons = linspace(1, 40, iterations);
%  received_photons = [50];
arraylength = length(received_photons);
% received_photons = 10.^ary;
%

%% Problem settings
SRF = 10;
%

%% FFT frame creation
repetition = 10000;
frequency = 15.5;
bandwidth = 50;
sampling_rate = 2*bandwidth;
sampling_period = 1/sampling_rate;
initial_phase = 0;
number_of_points = 100;
observation_time = sampling_period * (number_of_points);
frequency_resolution = 1/observation_time;
%

%% Create input sample sinusoidal tone
index = linspace(0, number_of_points-1, number_of_points);
index_SR = linspace(0, number_of_points*SRF-1, number_of_points*SRF);
t = sampling_period * index;
f = frequency_resolution * index;
t_SR = sampling_period * index_SR; 
f_SR = (frequency_resolution / SRF) * index_SR;
Ps = electron_charge * received_photons / (responsivity * observation_time);
% Ps_ref = Ps(1);
amplitude = 2*responsivity*sqrt(Ps);
amplitude = amplitude/amplitude(end);
% amplitude = ;
x = amplitude' * sin(2*pi*frequency*t + initial_phase);
x = repmat(x,repetition,1);
%

%% Additive Gaussian noise creation
noise_floor = 2 * electron_charge * responsivity;
vari = (sampling_rate / 2) * noise_floor;
% sigma = sqrt(vari);
sigma = sqrt(vari)/(2*responsivity*sqrt(Ps(end)));
noise = sigma * randn(repetition*arraylength,number_of_points);
SNR = 20*log10((amplitude/sqrt(2))/sigma);
delta_noise = sqrt(number_of_points + sqrt(2*number_of_points))*sigma; % good approximation of norm(z)

y = x + noise;
%%
% 
% %% l1-min for super-res search
% asize = 101;
% a = linspace(0.5,1.5,asize);
% solutionl1 = zeros(repetition,arraylength); 
% solutionl2 = zeros(repetition,arraylength); 
% RMSE_l1 = zeros(1,arraylength);
% RMSE_l2 = zeros(1,arraylength);
% 
% % for k=1:arraylength
% %     for l=1:repetition
% %         error = zeros(length(f_SR)/2,asize);
% % 
% %         for i=1:length(f_SR)/2
% %             temp = a'*sin(2*pi*f_SR(i)*t+initial_phase);
% %             for j=1:asize
% %                 error(i,j) = norm(temp(j,:)-y((l-1)*arraylength+k,:));
% %             end
% %         end
% %         % figure();
% %         % surf(a, f_SR(1:end/2), error);
% %         % xlabel('amplitude','FontSize',15);
% %         % ylabel('frequency','FontSize',15);
% % 
% %         [ymin,idxmin] = min(min(error,[],2));
% %         solution(l,k) = f_SR(idxmin);
% %     end
% %     RMSE_l1(k) = rms(solution(:,k) - frequency);
% % end
% 
% 
% for k=1:arraylength
%     parfor l=1:repetition
% %         error = zeros(length(f_SR)/2,asize);
%         errorl1 = zeros(1, asize);
%         errorl2 = zeros(1, asize);
%         avectorl1 = zeros(1, length(f_SR)/2);
%         avectorl2 = zeros(1, length(f_SR)/2);
% %         avector(1) = inf;
% 
%         for i=1:length(f_SR)/2
%             temp = a'*sin(2*pi*f_SR(i)*t+initial_phase);
%             for j=1:asize
% %                 error(i,j) = norm(temp(j,:)-y((l-1)*arraylength+k,:),1);
%                 errorl1(j) = norm(temp(j,:)-y((l-1)*arraylength+k,:),1);
%                 errorl2(j) = norm(temp(j,:)-y((l-1)*arraylength+k,:),2);
% %                 err = norm(temp(j,:)-y((l-1)*arraylength+k,:));
% %                 if err <= delta_noise
% %                     tempa = [tempa a(j)];
% %                 end
%             end
%             avectorl1(i) = min(errorl1);
%             avectorl2(i) = min(errorl2);
% %             if length(tempa) ~= 0
% %                 avector(i) = min(tempa);
% %             else
% %                 avector(i) = inf;
% %             end
%         end
%         % figure();
%         % surf(a, f_SR(1:end/2), error);
%         % xlabel('amplitude','FontSize',15);
%         % ylabel('frequency','FontSize',15);
% 
% %         [ymin,idxmin] = min(min(error,[],2));
%         [ymin,idxmin] = min(avectorl1);
%         solutionl1(l,k) = f_SR(idxmin);
%         [ymin,idxmin] = min(avectorl2);
%         solutionl2(l,k) = f_SR(idxmin);
%         fprintf('array: %d, repetition: %d\n',k,l)
%     end
%     RMSE_l1(k) = rms(solutionl1(:,k) - frequency);
%     RMSE_l2(k) = rms(solutionl2(:,k) - frequency);
% end

%% theoretical DFT error

SER = zeros(1,iterations);
% NN = number_of_points/2;
NN = bandwidth;

for j=1:iterations
    for i=1:NN-1
        temp = ((-1)^(i+1))*nchoosek(NN-1,i)*(1/(i+1))*exp(-(i*received_photons(j))/(i+1));
        SER(j) = SER(j)+temp;
    end
end

P1 = ones(1,iterations)-SER;
P2 = SER/(NN-1);

RMSE = zeros(1,iterations);
for j=1:iterations
    for i=1:NN
        if i ~= 15
            RMSE(j) = RMSE(j) + (i-frequency)^2 * P2(j);
        else
            RMSE(j) = RMSE(j) + (i-frequency)^2 * P1(j);
        end
    end
end
RMSE = sqrt(RMSE);

%


%% Taking FFT
fx = zeros(arraylength, repetition, number_of_points);
for i=1:repetition
    for j=1:arraylength
        fx(j,i,:) = fft(y((i-1)*arraylength+j,:), number_of_points);
    end
end
%

% %% Take MUSIC
% result_music = zeros(arraylength, repetition);
% for i=1:repetition
%     for j=1:arraylength
%         [s, w] = pmusic(y((i-1)*arraylength+j,:),2,4096);
%         [pwr, idx] = max(s);
%         result_music(j,i) = w(idx)*100/pi;
%     end
% end

%% Summarize the result
result_fx = 2*abs(fx(:,:,1:end/2)/(number_of_points)).^2;
[maxpower, detection_result] = max(result_fx,[],3);
detection_result = detection_result - 1;

RMSE_fft = rms(detection_result' - frequency,1);
% RMSE_music = rms(result_music' - frequency,1);

single_spectrum = result_fx(:,1,:);
averaged_spectrum = mean(result_fx(:,:,:),2);

sgl = zeros(1,bandwidth);
avge = zeros(1,bandwidth);

for i=1:bandwidth
    sgl(i) = single_spectrum(1,1,i);
    avge(i) = averaged_spectrum(1,1,i);
end
%

%% Plot the result
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

% 
figure(1);
plot(received_photons, RMSE_fft,'LineWidth',2);
hold on;
plot(received_photons, RMSE,'LineWidth',2);
% hold on;
% plot(received_photons, RMSE_l1,'LineWidth',2);
% hold on;
% plot(received_photons, RMSE_l2,'LineWidth',2);
title('Detection Error for DFT (Low SNR)','FontSize',18,'FontWeight','normal');
xlabel('SNR','FontSize',15);
ylabel('Error (Hz)','FontSize',15);
% legend('DFT','MUSIC','L1','L2');
legend('DFT','DFT Theory');
grid on;
set(gca,'xlim',[0 41])
set(gca,'xtick',[0:5:40])
set(gca,'ylim',[0 20])
set(gca,'ytick',[0:5:20])
% matlab2tikz('./final_report/figure2a.tex');
%

%% Save data
% 
% savefile_fft = './saved_data/fft.mat';
% savefile_music = './saved_data/music.mat';
% savefile_l1 = './saved_data/l1.mat';
% savefile_l2 = './saved_data/l2.mat';
% 
% save(savefile_fft,'RMSE_fft');
% save(savefile_music,'RMSE_music');
% save(savefile_l1,'RMSE_l1');
% save(savefile_l2,'RMSE_l2');
% 
