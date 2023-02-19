start_yr = 5;
end_yr = 21;
dt = 0.1;
crop_i = start_yr/dt:(end_yr/dt);
md_grid_mid = squeeze(md_grid(mid_i,:,:));
md_grid_mid = md_grid_mid(:,crop_i);
LTs = zeros(size(md_grid_mid,1), length(crop_i));
for i = 1:size(md_grid_mid,1)
    [LT, ST] = trenddecomp(md_grid_mid(i,:));
%     % inspect the period from the ST array using FFT
%     % stop when the period is less than 3 years
%     j = 1; % only inspect the first one in ST array
%     data = ST(:,j);
%     dt = 0.1; % year
%     Fs = 1/dt; % sampling frequency: n of samples per year (or per s, based on what time unit you use)
%     T = 1/Fs; % sampling period: inverse of freq
%     total_t = length(data)/Fs; % number of year
%     L = total_t/T;
%     data_fft = fft(data);
%     % construct power spectrum
%     P2 = abs(data_fft/L);
%     P1 = P2(1:L/2+1);
%     P1(2:end-1) = 2*P1(2:end-1);
%     % get the prominent
%     f = Fs*(0:(L/2))/L;
%     [~, max_i] = max(P1);
%     period = 1/f(max_i);
%     if period >= 3 && period < inf
%         % add to long term
%         LT = LT + data;
%         disp(['Found one LT hiding in ST! The period is ',num2str(period)])
%     end
%     LT = LT-LT(1); % displaced to match delta H = 0 at t=0
    LTs(i,:) = LT;
end
STs = md_grid_mid - LTs;

%% plot timeseries
t = linspace(0,26,260);
figure;
subplot(1,2,1)
plot(t, LTs);
subplot(1,2,2)
plot(t, STs);

%% Plot center flow line evolution
t = linspace(5,21,160);
t_axis = 50:209;
x = linspace(0,60000,1200);
window_size = 50;
figure
for i = 1:length(t)
    ti = t_axis(i);
    plot(x,movmean(LTs(:,ti),window_size),'-b'); hold on
    pause(0.05)
    ylim([-5,2])
end

%% Plot center flow line evolution (non-detrended ts)
t = linspace(0,26,260);
x = linspace(0,60000,1200);
figure
for i = 1:length(t)
    plot(x,md_grid_mid(:,i),'-b'); hold on
    pause(0.05)
    ylim([-5,2])
end

%% Find prominent period of the sine wave
dt = 0.1; % year
Fs = 1/dt; % sampling frequency: n of samples per year (or per s, based on what time unit you use)
T = 1/Fs; % sampling period: inverse of freq
total_t = length(data)/Fs; % number of year
L = total_t/T;
data_fft = fft(data);
% construct power spectrum
P2 = abs(data_fft/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
% get the prominent 
f = Fs*(0:(L/2))/L;
[~, max_i] = max(P1);
period = 1/f(max_i);
