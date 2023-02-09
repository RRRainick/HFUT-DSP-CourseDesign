%% 生成复序列（长度2048）的对应旋转因子查找表
clear;
N = 2048; %% 2048 = 2^{11}
loop = 2.^(log2(N)-1: -1:0); %% 长度为2048的序列共有11个stage
n = (0: N-1); %% 每个stage有2048个旋转因子
WN = exp(-1j*2*pi./N.*n'.*loop); %% 生成旋转因子
n_WN = reshape(WN, [16, 1408]); %% reshape矩阵使其长度限制为16，方便在编辑器中查看
file_dir = fopen('../../data/complex_WN.txt', 'w'); %% 打开输出文件
for i = (1: 1408)
    real_WN = real(n_WN(:, i));
    imag_WN = imag(n_WN(:, i));
    fprintf(file_dir, '%f+%f*I,\t', [real_WN';imag_WN']); %% 由于printf()函数输出矩阵按列输出，因此需要转换为列向量
    fprintf(file_dir, '\n');
end
fclose(file_dir); %% 关闭文件

%{
%% 生成实序列（长度4096）的对应旋转因子查找表
clear;
N = 4096; %% 2^{12}
loop = 2.^(log2(N)-1: -1: 0);
n = (0: N-1);
WN = exp(-1j*2*pi./N.*n'.*loop);
n_WN = reshape(WN, [16, 3072]); 
file_dir = fopen('../../data/complex_WN.txt', 'w');
for i = (1: 3072)
    real_WN = real(n_WN(:, i));
    imag_WN = imag(n_WN(:, i));
    fprintf(file_dir, '%f+%f*I,\t', [real_WN';imag_WN']);
    fprintf(file_dir, '\n');
end
fclose(file_dir);
%}

%% 生成实序列（长度4096）的对应旋转因子查找表
clear;
N = 4096; %% 2^{12}
n = (0: N-1);
WN = exp(-1j*2*pi./N.*n');
n_WN = reshape(WN, [16, N/16]); 
file_dir = fopen('../../data/real_WN.txt', 'w');
for i = (1: N/16)
    real_WN = real(n_WN(:, i));
    imag_WN = imag(n_WN(:, i));
    fprintf(file_dir, '%f+%f*I,\t', [real_WN';imag_WN']);
    fprintf(file_dir, '\n');
end
fclose(file_dir);