clear; close all;

xn = readmatrix('../../data/transmitter.txt'); %% 读取gen_input.m生成的信号数据
Nx = length(xn);
n = (0: Nx-1);
xn_f = single(xn); %% 将数据转为单精度浮点数，后续fft()函数对单精度数据会自动切换为单精度FFT
xn_amax = max(abs(xn)); xn_amin = min(abs(xn)); %% 计算幅度最值，用于确定绘图坐标轴范围

time = hat(); %% 计时器计时开始
Xk = fft(xn_f);
diff = hat() - time; %% 计时器计时结束，计算计时数值差，单位ms
Nk = length(Xk);
k = (0: Nk-1);
Xk_amax = max(abs(Xk)); Xk_amin = min(abs(Xk));
Xk_pmax = max(angle(Xk)); Xk_pmin = min(angle(Xk));
xn_k = ifft(Xk); %% IDFT恢复原信号

Xkc = readmatrix('../../data/output.txt'); %% 读取C语言FFT算法的输出数据
Nkc = length(Xkc);
kc = (0: Nkc-1);
Xkc_f = single(Xkc);
Xkc_amax = max(abs(Xkc)); Xkc_amin = min(abs(Xkc));
Xkc_pmax = max(angle(Xkc)); Xkc_pmin = min(angle(Xkc));
xnc = ifft(Xkc);
Nxc = length(xnc);
nc = (0: Nxc-1);
xnc_amax = max(abs(xnc)); xnc_amin = min(abs(xnc));

figure('name', '基于C语言的单精度FFT算法实现验证');
subplot(3, 3, [1, 2, 3]); stem(n, abs(xn_f)); axis([0, Nx-1, xn_amin - xn_amin * 0.1, xn_amax + xn_amax * 0.1]); title('原信号x[n]时域图像');
subplot(3, 3, 4); stem(k, abs(Xk));  axis([0, Nk-1, Xk_amin - Xk_amin * 0.1, Xk_amax + Xk_amax * 0.1]); title('matlab FFT幅频特性');
subplot(3, 3, 5); stem(k, angle(Xk));  axis([0, Nk-1, Xk_pmin - Xk_pmin * 0.1, Xk_pmax + Xk_pmax * 0.1]); title('matlab FFT相频特性');
subplot(3, 3, 6); stem(n, abs(xn_k)); axis([0, Nx-1, xn_amin - xn_amin * 0.1, xn_amax + xn_amax * 0.1]); title('matlab IFFT恢复信号'); 
subplot(3, 3, 7); stem(kc, abs(Xkc_f));  axis([0, Nkc-1, Xkc_amin - Xkc_amin * 0.1, Xkc_amax + Xkc_amax * 0.1]); title('C FFT幅频特性');
subplot(3, 3, 8); stem(kc, angle(Xkc_f));  axis([0, Nkc-1, Xkc_pmin - Xkc_pmin * 0.1, Xkc_pmax + Xkc_pmax * 0.1]); title('C FFT相频特性');
subplot(3, 3, 9); stem(nc, abs(xnc)); axis([0, Nxc-1, xnc_amin - xnc_amin * 0.1, xnc_amax + xnc_amax * 0.1]); title('C IFFT恢复信号'); 
fprintf('FFT cost %fms\n', diff);