%% 生成复序列（长度2048）的对应bitreverse序列查找表
clear;
N = 2048;
n = (0: N-1);
x = bitrevorder(n); %% 生成bitreverse序列
xm = reshape(x, [16, 128]); %% reshape矩阵使其长度限制为16，方便在编辑器中查看
file_dir = fopen('../../data/complex_br_table.txt', 'w'); %% 打开输出文件
for i = (1: 128)
    fprintf(file_dir, '%d,\t', xm(:, i));
    fprintf(file_dir, '\n');
end
fclose(file_dir);

%% 生成实序列（长度4096）的对应bitreverse序列查找表
clear;
N = 4096;
n = (0: N-1);
x = bitrevorder(n);
xm = reshape(x, [16, 256]);
file_dir = fopen('../../data/real_br_table.txt', 'w');
for i = (1: 256)
    fprintf(file_dir, '%d,\t', xm(:, i));
    fprintf(file_dir, '\n');
end
fclose(file_dir);