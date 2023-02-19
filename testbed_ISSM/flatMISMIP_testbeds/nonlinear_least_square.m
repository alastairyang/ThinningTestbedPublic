%% decompose the sine wave and linear trend
x = 0:0.01:10; 
p_t = [0.05, 0.02, 10, 0.01];
y = p_t(1)*x + p_t(2)*x.*sin(p_t(3)*x + p_t(4));
fun1 = @(p, xdata) p(1)*xdata + p(2)*xdata.*sin(p(3)*x + p(4));
p0 = [0.1, 0.1, 7, 0.03];
options = optimset('TolFun', 1e-10);
p_inv = lsqcurvefit(fun1, p0, x, y,[],[],options);
figure; plot(x, y, '-r', x, fun1(p_inv, x), '-.b')

%% decompose the multiples of exponential decay
x = 0:0.01:10;
p_t = [3, 10, 0.3]; % true parameters
y = p_t(1)*x.*exp(-p_t(2)*rem(x,2)) + p_t(3)*x;
fun2 =  @(p,xdata) p(1)*x.*exp(-p(2)*rem(xdata,2)) + p(3)*xdata;
p0 = [1, 8, 0.5];
options = optimset('TolFun', 1e-10);
p_inv = lsqcurvefit(fun2, p0, x, y,[],[],options);
figure; plot(x, y, '-r', x, fun2(p_inv, x), '-.b')