%% Script to perform super-resolution of spikes from low frequency
addpath ./matlab2tikz-matlab2tikz-816f875/src
clear all;
close all;
%

%% Physical constants
responsivity = 0.7;
electron_charge = 1.6e-19;
%

%% Receiving power
% ary = linspace(0, 2.5, 11);
% ary = 1.7;
% received_photons = linspace(1, 41, 21);
received_photons = [40];
arraylength = length(received_photons);
% received_photons = 10.^ary;
%

%% Problem settings
SRF = 10;
%

%% FFT frame creation
repetition = 1;
frequency = 15.5;
bandwidth = 100;
sampling_rate = 2*bandwidth;
sampling_period = 1/sampling_rate;
initial_phase = 0;
number_of_points = 200;
number_of_points_SR = number_of_points * SRF;
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
x = amplitude' * cos(2*pi*frequency*t + initial_phase);
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
y_SR = zeros(repetition*arraylength,number_of_points_SR);
y_SR(:,1:length(y)) = y;

y_SR_freq = abs(fft(y_SR,[],2));
%

%% l1-min for super-res search
asize = 21;
a = linspace(0.6,1.1,asize);
solution = zeros(repetition,arraylength); 
solution_new = zeros(repetition,arraylength); 
RMSE_l1 = zeros(1,arraylength);
RMSE_l1_new = zeros(1,arraylength);
% margin = 0.5*(length(f_SR)-length(y));

for k=1:arraylength
    for l=1:repetition
        l1norm = zeros(length(f_SR)/2,asize);
        error = zeros(length(f_SR)/2,asize);
        l1norm_new = zeros(length(f_SR)/2,asize);
        error_new = zeros(length(f_SR)/2,asize);
%         error_th = zeros(length(f_SR)/2,asize);
%         errorvector = zeros(1, length(f_SR)/2);

        for i=1:length(f_SR)/2
            temp = a'*cos(2*pi*f_SR(i)*t+initial_phase);
%             temp = zeros(asize,length(f_SR));
            
            temp_new = zeros(asize,length(f_SR));
            temp_new(:,1:length(y)) = a'*cos(2*pi*f_SR(i)*t+initial_phase);
           
            
            temperror = [];
            for j=1:asize
                temp_lp = fft(temp_new(j,:));
                l1norm_new(i,j) = norm(temp_lp-y_SR_freq((l-1)*arraylength+k,:),1);
                error_new(i,j) = norm(temp_lp-y_SR_freq((l-1)*arraylength+k,:),1);
                 
                l1norm(i,j) = norm(temp(j,:)-y((l-1)*arraylength+k,:),1);
                error(i,j) = norm(temp(j,:)-y((l-1)*arraylength+k,:),1);
%                 err = norm(temp(j,:)-y((l-1)*arraylength+k,:));
%                 if error(i,j) <= delta_noise
% %                     temperror = [temperror err];
%                     error_th(i,j) = error(i,j);
%                 else
%                     error_th(i,j) = inf;
%                 end
            end
%             if length(temperror) ~= 0
%                 errorvector(i) = min(temperror);
%             else
%                 errorvector(i) = inf;
%             end
        end
        figure(1);
        surf(a, f_SR(1:end/2), error);
        xlabel('Amplitude','FontSize',15);
        ylabel('Frequency (Hz)','FontSize',15);
        zlabel('L1 Error','FontSize',15);
%         matlab2tikz('./final_report/figure4b.tex');
        
        figure(2);
        surf(a, f_SR(1:end/2), error_new);
        xlabel('Amplitude','FontSize',15);
        ylabel('Frequency (Hz)','FontSize',15);
        zlabel('L1 Error','FontSize',15);
%         matlab2tikz('./final_report/figure4b.tex');

        [ymin,idxmin] = min(min(error,[],2));
%         [ymin,idxmin] = min(errorvector);
        solution(l,k) = f_SR(idxmin);
        [ymin,idxmin] = min(min(error_new,[],2));
        solution_new(l,k) = f_SR(idxmin);
    end
    RMSE_l1(k) = rms(solution(:,k) - frequency);
    RMSE_l1_new(k) = rms(solution_new(:,k) - frequency);
end



%

%%

% 
% %% Taking FFT
% fx = zeros(arraylength, repetition, number_of_points);
% for i=1:repetition
%     for j=1:arraylength
%         fx(j,i,:) = fft(y((i-1)*arraylength+j,:), number_of_points);
%     end
% end
% %
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
% 
% %% Summarize the result
% result_fx = 2*abs(fx(:,:,1:end/2)/(number_of_points)).^2;
% [maxpower, detection_result] = max(result_fx,[],3);
% detection_result = detection_result - 1;
% 
% RMSE_fft = rms(detection_result' - frequency,1);
% RMSE_music = rms(result_music' - frequency,1);
% 
% single_spectrum = result_fx(:,1,:);
% averaged_spectrum = mean(result_fx(:,:,:),2);
% 
% sgl = zeros(1,bandwidth);
% avge = zeros(1,bandwidth);
% 
% for i=1:bandwidth
%     sgl(i) = single_spectrum(1,1,i);
%     avge(i) = averaged_spectrum(1,1,i);
% end
% %
% 
% %% Plot the result
% % 
% % figure(1);
% % plot(f(1:end/2),db(sgl,'power'));
% % set(gca,'yscale','log');
% % 
% % figure(2);
% % plot(f(1:end/2),db(avge,'power'));
% % set(gca,'yscale','log');
% % 
% % figure(3);
% % % scatter( ones(1,length(detection_result)) ,detection_result);
% % histogram(detection_result, 'BinWidth', 1, 'BinLimits', [-0.5, 99.5]);
% % 
% % figure(4);
% % % scatter( ones(1,length(detection_result)) ,detection_result);
% % histogram(result_music, 'BinWidth', 1, 'BinLimits', [-0.5, 99.5]);
% 
% % 
% figure(1);
% plot(received_photons(end-10:end), RMSE_fft(end-10:end),'LineWidth',2);
% hold on;
% plot(received_photons(end-10:end), RMSE_music(end-10:end),'LineWidth',2);
% hold on;
% plot(received_photons(end-10:end), RMSE_l1(end-10:end),'LineWidth',2);
% title('Detection Error for DFT/MUSIC/L1 (Low SNR)','FontSize',18,'FontWeight','normal');
% xlabel('SNR','FontSize',15);
% ylabel('Error (Hz)','FontSize',15);
% legend('DFT','MUSIC','L1');
% grid on;
% set(gca,'xlim',[0 41])
% set(gca,'xtick',[0:10:40])
% set(gca,'ylim',[0 50])
% set(gca,'ytick',[0:10:50])
% % matlab2tikz('./final_report/figure2a.tex');
% %
