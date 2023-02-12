clear; close all;

xn = readmatrix('../../data/transmitter.txt'); %% 读取gen_input.m生成的信号数据
Nx = length(xn);
n = (0: Nx-1);
xn_f = single(xn); %% 将数据转为单精度浮点数，后续fft()函数对单精度数据会自动切换为单精度FFT
xn_x_lim = [0, Nx-1];

time = hat(); %% 计时器计时开始
Xk = fft(xn_f);
diff = hat() - time; %% 计时器计时结束，计算计时数值差，单位s，精确到us;
duration = diff * 1e6;
Nk = length(Xk);
k = (0: Nk-1);
Xk_x_lim = [0, Nk-1]; 
Xk_amax = max(abs(Xk)); Xk_amin = min(abs(Xk)); 
Xk_arange = Xk_amax - Xk_amin; Xk_ay_lim = [Xk_amin-Xk_arange*0.1, Xk_amax+Xk_arange*0.1]; %% 计算幅度最值，用于确定绘图坐标轴范围
Xk_pmax = max(angle(Xk)); Xk_pmin = min(angle(Xk));
Xk_prange = Xk_pmax - Xk_pmin; Xk_py_lim = [Xk_pmin-Xk_prange*0.1, Xk_pmax+Xk_prange*0.1];
xn_k = ifft(Xk); %% IDFT恢复原信号

Xkc = readmatrix('../../data/output.txt'); %% 读取C语言FFT算法的输出数据
Nkc = length(Xkc);
kc = (0: Nkc-1);
Xkc_f = single(Xkc);
Xkc_x_lim = [0, Nkc-1];
Xkc_amax = max(abs(Xkc)); Xkc_amin = min(abs(Xkc));
Xkc_arange = Xkc_amax - Xkc_amin; Xkc_ay_lim = [Xkc_amin-Xkc_arange*0.1, Xkc_amax+Xkc_arange*0.1]; %% 计算幅度最值，用于确定绘图坐标轴范围
Xkc_pmax = max(angle(Xkc)); Xkc_pmin = min(angle(Xkc));
Xkc_prange = Xkc_pmax - Xkc_pmin; Xkc_py_lim = [Xkc_pmin-Xkc_prange*0.1, Xkc_pmax+Xkc_prange*0.1]; %% 计算幅度最值，用于确定绘图坐标轴范围
xnc = ifft(Xkc);
Nxc = length(xnc);
nc = (0: Nxc-1);
xnc_x_lim = [0, Nxc-1];

fs = 1e6;
res = fs/Nx; %% omg_res = 2*pi*fs/Nx;
resc = fs/Nkc; %% omg_resc = 2*pi*fs/Nkc;
xt_x_lim = xn_x_lim ./ fs;  Xomg_x_lim = [0, fs/2]; Xomg_py_lim = [-pi, pi];
xtc_x_lim = xnc_x_lim ./ fs; Xomgc_x_lim = [0, fs/2]; Xomgc_py_lim = [-pi, pi];

figure1 = figure('name', '基于C语言的单精度FFT算法实现验证');
if ~isreal(xn)
    xn_amax = max(abs(xn)); xn_amin = min(abs(xn)); %% 计算幅度最值，用于确定绘图坐标轴范围
    xn_arange = xn_amax - xn_amin; xn_ay_lim = [xn_amin-xn_arange*0.1, xn_amax+xn_arange*0.1]; %% 计算幅度最值，用于确定绘图坐标轴范围
    %% 若原信号为复信号，则图像用幅时图像表示
    %% 注意y轴范围为xn_amin - xn_amin * 0.1 
    subplot(3, 3, [1, 2, 3]); stem(n, abs(xn_f)); axis([xn_x_lim, xn_ay_lim]); title('原信号x[n]时域图像');
    subplot(3, 3, 6); stem(n, abs(xn_k)); axis([xn_x_lim, xn_ay_lim]); title('matlab IFFT恢复信号'); 
else
    %% 若原信号为是信号，则无需计算幅度
    %% 注意y轴范围为xn_amin + xn_amin * 0.1
    xn_amax = max(xn); xn_amin = min(xn); 
    xn_arange = xn_amax - xn_amin; xn_ay_lim = [xn_amin-xn_arange*0.1, xn_amax+xn_arange*0.1]; %% 计算幅度最值，用于确定绘图坐标轴范围
    subplot(3, 3, [1, 2, 3]); stem(n, xn_f); axis([xn_x_lim, xn_ay_lim]); title('原信号x[n]时域图像');
    subplot(3, 3, 6); stem(n, xn_k); axis([xn_x_lim, xn_ay_lim]); title('matlab IFFT恢复信号'); 
end
if ~isreal(xnc)
    xnc_amax = max(abs(xnc)); xnc_amin = min(abs(xnc));
    xnc_arange = xnc_amax - xnc_amin; xnc_ay_lim = [xnc_amin-xnc_arange*0.1, xnc_amax+xnc_arange*0.1]; %% 计算幅度最值，用于确定绘图坐标轴范围
    subplot(3, 3, 9); stem(nc, abs(xnc)); axis([xnc_x_lim, xnc_ay_lim]); title('C IFFT恢复信号'); 
else
    xnc_amax = max(xnc); xnc_amin = min(xnc);
    xnc_arange = xnc_amax - xnc_amin; xnc_ay_lim = [xnc_amin-xnc_arange*0.1, xnc_amax+xnc_arange*0.1]; %% 计算幅度最值，用于确定绘图坐标轴范围
    subplot(3, 3, 9); stem(nc, xnc); axis([xnc_x_lim, xnc_ay_lim]); title('C IFFT恢复信号'); 
end
    subplot(3, 3, 4); stem(k, abs(Xk));  axis([Xk_x_lim, Xk_ay_lim]); title('matlab FFT幅频特性');
    subplot(3, 3, 5); stem(k, angle(Xk));  axis([Xk_x_lim, Xk_py_lim]); title('matlab FFT相频特性');
    subplot(3, 3, 7); stem(kc, abs(Xkc_f));  axis([Xkc_x_lim, Xkc_ay_lim]); title('C FFT幅频特性');
    subplot(3, 3, 8); stem(kc, angle(Xkc_f));  axis([Xkc_x_lim, Xkc_py_lim]); title('C FFT相频特性');

figure2 = figure('name', 'FFT的频域分析');
if ~isreal(xn)
    subplot(3, 3, [1, 2, 3]); plot(n/fs, abs(xn_f)); axis([xt_x_lim, xn_ay_lim]); title('原信号x[n]时域图像');
    subplot(3, 3, 6); plot(n/fs, abs(xn_k)); axis([xt_x_lim, xn_ay_lim]); title('matlab IFFT恢复信号'); 
else
    subplot(3, 3, [1, 2, 3]); plot(n/fs, xn_f); axis([xt_x_lim, xn_ay_lim]); title('原信号x[n]时域图像');
    subplot(3, 3, 6); plot(n/fs, xn_k); axis([xt_x_lim, xn_ay_lim]); title('matlab IFFT恢复信号'); 
end
if ~isreal(xnc)
    subplot(3, 3, 9); plot(nc/fs, abs(xnc)); axis([xtc_x_lim, xnc_ay_lim]); title('C IFFT恢复信号'); 
else
    subplot(3, 3, 9); plot(nc/fs, xnc); axis([xtc_x_lim, xnc_ay_lim]); title('C IFFT恢复信号'); 
end
    subplot(3, 3, 4); plot(k*res, abs(Xk));  axis([Xomg_x_lim, Xk_ay_lim]); title('matlab FFT幅频特性');
    subplot(3, 3, 5); plot(k*res, angle(Xk));  axis([Xomg_x_lim, Xomg_py_lim]); title('matlab FFT相频特性');
    subplot(3, 3, 7); plot(kc*resc, abs(Xkc_f));  axis([Xomgc_x_lim, Xkc_ay_lim]); title('C FFT幅频特性');
    subplot(3, 3, 8); plot(kc*resc, angle(Xkc_f));  axis([Xomgc_x_lim, Xomgc_py_lim]); title('C FFT相频特性');

fprintf("duration is %fus\n", duration);
