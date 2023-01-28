phis = zeros(size(perturb_xi));
count = 0;

for i = perturb_xi
    count = count + 1;
    signal = STs_fine(:,i);
    
    T = 0.01; % sampling period, think of it as the smallest unit of t. Can be simply dt in our data
    Fs = 1/T; % sampling frequency
    last_t = 15.89; % last year, counting from t = 0
    L = (last_t/T)+1; % length of the signal (total number of sampling periods)
    t = (0:L-1)*T; % make the time axis
    
    y = fft(signal);
    z = fftshift(y);
    ly = length(y);
    f_new = (-ly/2:ly/2-1)/ly*Fs;
    % plot the double-sided frequency spectrum
    %figure; stem(f_new, abs(z))
    
    % now if you don't like double-sided, you can do the following, which keeps
    % just one half
    P2 = abs(y/L);
    P1 = P2(1:L/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    f = Fs*(0:(L/2))/L;
    figure;plot(f,P1)
    
    tol = 1e-6;
    y(abs(y)<tol) = 0;
    theta = angle(y);
    theta2 = theta(1:L/2+1);
    
    % find the prominent freq
    [~,idx] = max(P1);
    phi = theta2(idx)/pi;
    phis(count) = phi;
end