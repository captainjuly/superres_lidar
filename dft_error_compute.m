addpath ./matlab2tikz-matlab2tikz-816f875/src

% Compute theoretical RMSE
iterbegin = 1;
iterations = 21;
% received_photons = linspace(iterbegin,iterbegin+iterations-1,iterations);
received_photons = linspace(1, 20, iterations);

SER = zeros(1,iterations);
% NN = number_of_points/2;
NN = 50;
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

figure(1);
plot(received_photons, SER,'LineWidth',2);
title('Symbol Error Rate in 50-ary FSK Receiver','FontSize',18,'FontWeight','normal');
xlabel('SNR','FontSize',15);
ylabel('Symbol Error Rate','FontSize',15);
grid on;
set(gca,'xtick',[0:5:20])
set(gca,'ylim',[-0.1 1])
set(gca,'ytick',[0:0.2:1])
matlab2tikz('./final_report/figure1a.tex');

figure(2);
plot(received_photons, RMSE,'LineWidth',2);
title('RMS Error for Frequency Estimation','FontSize',18,'FontWeight','normal');
xlabel('SNR','FontSize',15);
ylabel('Error (Hz)','FontSize',15);
grid on;
set(gca,'xtick',[0:5:20])
set(gca,'ylim',[0 20])
set(gca,'ytick',[0:5:20])
matlab2tikz('./final_report/figure1b.tex');
% set(gca,'xscale','log');
% set(gca,'yscale','log');
% figure(2);
% plot(received_photons, RMSE);
% title('RMS Error for Frequency Estimation');
% xlabel('SNR');
% ylabel('Error (Hz)');
% grid on;
% set(gca,'xscale','log');
% set(gca,'yscale','log');