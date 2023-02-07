clear; close all;
N1 = 2048; N2 = 4096; %% 实序列和复序列的长度
f = 100;
%% 测试信号
%% n = (0: N1-1); xn = cos(f*n) + 1j * sin(f/2*n); %% 复序列
%% n = (0: N2-1); xn = cos(1/f*n);  %% 实序列
%% n = (0: N2-1); xn = cos(f*n);    %% 补零实序列对照
%% n = (0: N2/2-1); xn = cos(f*n);  %% 补零实序列对照
%% n = (0: N1/2-1); xn = cos(f*n) + 1j * sin(f/2*n); %% 补零复序列
%% n = (0: N2 + N1); xn = cos(1/f*n);  %% 截断实序列

%% 常见信号
%% n = (0: N2-1); xn = stepfun(n, N2/2-1); %% 阶跃信号
%% n = (0: N2-1); xn = rectpuls(n-N1, N2/4-1);    %% 方波脉冲信号
%% n = (0: N2-1); xn = tripuls(n-N1, N2/4-1, 0);      %% 三角波脉冲信号
n = (0: N2-1); xn = square(50*n, 50);      %% 50%占空比的周期方波信号
%% n = (0: N2-1); xn = dirac(n);   %% 单位脉冲信号，不支持

writematrix(xn, "../../data/transmitter"); %% 传输给FFT_figure.m的数据文件
file_dir = fopen('../../data/test.txt', 'w'); %% 输出文件路径
if ~isreal(xn)
    fprintf(file_dir, '%f %f\n', [real(xn); imag(xn)]); %% 若为复序列则以空格分隔，输出实部和虚部
else
    fprintf(file_dir, '%f\n', real(xn)); %% 若为实序列则直接输出
end
fclose(file_dir);