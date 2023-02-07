clear; close all;
N1 = 2048; N2 = 4096; %% 实序列和复序列的长度
f = 100;
%% 测试信号
%% n = (0: N1-1); xn = cos(f*n) + 1j * sin(f/2*n);
%% n = (0: N2-1); xn = cos(f*n);
n = (0: N2/2-1); xn = cos(f*n);
writematrix(xn, "../../data/transmitter"); %% 传输给FFT_figure.m的数据文件
file_dir = fopen('../../data/test.txt', 'w'); %% 输出文件路径
if ~isreal(xn)
    fprintf(file_dir, '%f %f\n', [real(xn); imag(xn)]); %% 若为复序列则以空格分隔，输出实部和虚部
else
    fprintf(file_dir, '%f\n', real(xn)); %% 若为实序列则直接输出
end
fclose(file_dir);