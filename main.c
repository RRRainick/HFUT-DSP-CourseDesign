#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <windows.h>
#include "br_lookup_table.h"
#include "WN_lookup_table.h"

#ifndef LENGTH
    #define LENGTH 2048
#endif
#ifndef WIDTH
    #define WIDTH 2
#endif
#define epsilon 1e-7

LARGE_INTEGER freq;
LARGE_INTEGER start, end;
double rollover_us, reverse_us, btf_us, duration_us, real_us;

_Bool read_file(float *ptr_seq, int argc, char *argv[]);
void rollover_seq(float const (*ptr_input)[WIDTH], float complex *complex_xn);
void bitreorder_seq(float complex *xn);
void compute_butterfly(float complex *const r_xn);
void process_real_seq(_Bool complex_sign, float complex *Xk, float complex *Real_Xk);
void write_file(_Bool complex_sign, float complex *Xk, float complex *real_Xk);

int main(int argc, char *argv[])
{
    LARGE_INTEGER main_freq;
    LARGE_INTEGER main_start, main_end;
    float input_seq[LENGTH][WIDTH] = {0}; // 默认零矩阵，实现序列自动补零
    _Bool complex_sign = read_file(input_seq, argc, argv);
    float const (*ptr_input)[WIDTH] = input_seq;
    float complex xn[LENGTH]; // 原序列数组；FFT的结果由于原位运算也储存于这个数组
    float complex real_Xk[LENGTH * WIDTH];  // 实序列FFT数组，其长度与复序列FFT数组长度不同
    rollover_seq(ptr_input, xn);    // 重排输入序列为float complex数组
   
    QueryPerformanceFrequency(&main_freq);   // 获取CPU时钟周期
    QueryPerformanceCounter(&main_start);    // 获取当前经过的CPU时钟次数

    bitreorder_seq(xn); // 序列反比特排序

    compute_butterfly(xn);  // 蝶形运算
    
    process_real_seq(complex_sign, xn, real_Xk);    // 实序列修正

    QueryPerformanceCounter(&main_end);  // 获取当前经过的CPU时钟次数
    duration_us = (double)(main_end.QuadPart - main_start.QuadPart) * 1e6 / (double)(main_freq.QuadPart);    // 时钟次数作差并除以时钟频率得到经过的时间，以us为单位的

    write_file(complex_sign, xn, real_Xk);
    printf("rollover cost %fus, reverse cost %fus, compute butterfly cost %fus, real process cost %fus.\n", rollover_us, reverse_us, btf_us, real_us);
    printf("duration is %fus.\n", duration_us);
    return 0;
}

/*
 * 用于从matlab的输出文件中读取数据，将其由字符转为浮点数后储存于相应数组
 * ptr_seq：用于储存输入数据的数组指针
 * argc：命令参数个数
 * argv[]：命令参数字符数组
 * 返回值：复数标识符
 */
_Bool read_file(float *ptr_seq, int argc, char *argv[])
{
    unsigned char tmp_array[32]; // 用于储存str2float的数字，float类型在x64上占32bit
    char ch;
    unsigned int digit = 0; // 表示数字位数
    unsigned int count = LENGTH * WIDTH; // 下标计数器，避免数组越界
    double i_decp = 0.0;    // index_decimal_point，表示左起小数点的位置
    double tmp_num = 0.0;   // 临时储存数字和
    _Bool plus_sign = 1;    // 正数指示变量（主要用于判断负数）
    _Bool complex_sign = 0; // 复数指示变量
    FILE *fp;
    /*
     * 打开储存输入序列的txt文件
     * 对不满足要求的参数报错并退出
     */
    if(argc != 2){
        printf("usage of %s is incorrect.\n", argv[0]);
        exit(EXIT_FAILURE);
    }
    if((fp = fopen(argv[1], "r")) == NULL){
        printf("can not open file %s.\n", argv[1]);
        exit(EXIT_FAILURE);
    }
    /*
     * 按字符读取txt文件
     * 将char转为float
     * 数据储存到ptr_sequence指向的数组
     */
    while((ch = getc(fp)) != EOF && count >= 1){    //当文件结束或数组即将越界时终止循环, 同时截断长度超过2048的复序列和超过4096的实序列
        switch(ch){
            case '.':
                i_decp = digit;
                break;
            case '-':
                plus_sign = 0;
                break;
            case 'i':   // 读取到i即判断为复序列
                complex_sign = 1;
                break;
            case ' ':
            case '\n':
                if(i_decp == 0){    // 输入为整数的情况
                    for(int i = 0; i < digit; i++){
                            tmp_num += tmp_array[i] * pow(10.0, digit-1.0-i);
                        }
                }
                else{   // 输入为浮点数的情况
                    for(int i = 0; i < digit; i++){
                        tmp_num += tmp_array[i] * pow(10.0, i_decp-1.0-i);
                    }
                }
                if(!plus_sign)
                    tmp_num = (-1.0)*tmp_num;
                *ptr_seq = tmp_num; // 计算值写入数组
                ptr_seq++; count--;
                digit = 0; tmp_num = 0.0; plus_sign = 1; i_decp = 0.0;  // 变量复位
                break;
            default:
                tmp_array[digit++] = ch - '0'; // char2short
                break;
        }
    }
    fclose(fp);

    return complex_sign;
}

/*
 * 用于根据储存数组的实数和复数类型，将其进行数组降维和数据类型转换（float-->float complex)
 * ptr_input：二维数组指针
 * xn：转存的复序列指针
 * 返回值：对数组元素处理的对应数组指针
 */
void rollover_seq(float const (*ptr_input)[WIDTH], float complex *xn)
{
    QueryPerformanceFrequency(&freq);   // 获取CPU时钟周期
    QueryPerformanceCounter(&start);    // 获取当前经过的CPU时钟次数
    float complex *output_xn;
    /*
     * 对实序列：偶序列作为实部，奇序列作为虚部
     * 对复序列：实部对应实部，虚部对应虚部
     */
    for(int i = 0; i < LENGTH; i++){
        xn[i] = (*ptr_input)[0] + (*ptr_input)[1] * I;
        ptr_input++;
    }

    QueryPerformanceCounter(&end);  // 获取当前经过的CPU时钟次数
    rollover_us = (double)(end.QuadPart - start.QuadPart) * 1e6 / (double)(freq.QuadPart);    // 时钟次数作差并处以时钟频率得到经过的时间，以us为单位的
    
    return;
} 

/*
 * 将信号序列进行bitreverse排序
 * xn：信号序列的数组指针
*/
void bitreorder_seq(float complex *xn)
{
    QueryPerformanceFrequency(&freq);
    QueryPerformanceCounter(&start);
    // log_{2}^{2048}=11, log_{2}^{4096}=12
    unsigned int count = LENGTH; // 排序次数计数器
    unsigned int const *ptr_br_table = br_table;   // 查找表指针
    _Bool i_arr_sign[LENGTH * WIDTH] = {0}; // inverse_array_sign, 标识已经交换过的元素，防止重复交换
    
    for(unsigned int k = 0; k < count; k++){
        unsigned int index = ptr_br_table[k];
        if(i_arr_sign[k] == 1 || i_arr_sign[index] == 1)
            continue;  
        else{
            float complex tmp = xn[k];
            xn[k] = xn[index]; xn[index] = tmp; 
            i_arr_sign[k] = 1; i_arr_sign[index] = 1;
        }
    }

    QueryPerformanceCounter(&end);
    reverse_us = (double)(end.QuadPart - start.QuadPart) * 1e6 / (double)(freq.QuadPart);

    return;    
}

/*
 * 对bitreverse后的序列进行蝶形运算
 * r_xn：reverse_xn，已bitreverse后的序列的数组指针
 */
void compute_butterfly(float complex * const r_xn)
{
    QueryPerformanceFrequency(&freq);
    QueryPerformanceCounter(&start);
    unsigned int length = LENGTH;
    float complex const *r_ptr_WN_table = complex_WN_table;    // reference_ptr_WN_table，指向每个stage最初的旋转因子
    float complex const *ptr_WN_table = complex_WN_table;      // 旋转因子查找表指针
    float complex *ptr_r_xn = r_xn; 

    unsigned int max_stage = (unsigned int)log2((double)length);    //计算所需的stage个数
    unsigned int j_max = length / 2;    // 每个stage内循环次数
    unsigned int k_max = 1;     // 每个循环内的蝶形运算次数
    for(unsigned int i = 1; i <= max_stage; i++){
        for(unsigned int j = 1; j <= j_max; j++){
            for(unsigned int k = 1; k <= k_max; k++){
                float complex tmp1, tmp2;
                tmp1 = ptr_r_xn[0] + ptr_WN_table[0] * ptr_r_xn[k_max];
                tmp2 = ptr_r_xn[0] + ptr_WN_table[k_max] * ptr_r_xn[k_max];
                ptr_r_xn[0] = tmp1; ptr_r_xn[k_max] = tmp2;     // in-place computation
                ptr_r_xn++;
                ptr_WN_table++;
            }
            ptr_r_xn = r_xn + j * k_max * 2;    // 每轮j循环结束前，需要将序列和查找表指针指向相应元素
            ptr_WN_table = r_ptr_WN_table + j * k_max * 2;
        }
        ptr_r_xn = r_xn;    // j循环结束后，序列指针复位
        r_ptr_WN_table += length;   // 旋转因子查找表指针跳过4096个元素，原因取决于旋转因子查找表结构
        ptr_WN_table = r_ptr_WN_table;
        j_max /= 2;
        k_max *= 2;
    }

    QueryPerformanceCounter(&end);
    btf_us = (double)(end.QuadPart - start.QuadPart) * 1e6 / (double)(freq.QuadPart);

    return;
}

/*
 * 若输入序列为实序列，需要通过DFT的共轭对称性进行修正
 * complex_sign：复序列标识符
 * Xk：经过复序列FFT后的数组指针
 * real_Xk：用于储存实序列FFT后的数组指针*/
void process_real_seq(_Bool complex_sign, float complex *Xk, float complex *real_Xk)
{
    QueryPerformanceFrequency(&freq);
    QueryPerformanceCounter(&start);
    unsigned int length = LENGTH;
    float complex const *ptr_WN_table = real_WN_table;  // 旋转因子查找表指针
    float complex Xek[LENGTH]; float complex Xok[LENGTH];

    if(~complex_sign){
        for(unsigned int i = 1; i <= LENGTH - 1; i++){
            //Xek[i] = (Xk[i] + creal(Xk[length - i]) - I * cimag(Xk[length - i])) / 2;   // Xe[k] = 1/2 * (Y[k] + Y^{*}[N-k])    0<=k<=N-1
            //Xok[i] = (Xk[i] - (creal(Xk[length - i]) - I * cimag(Xk[length - i]))) / (2 * I);  // Xo[k] = /2i * (Y[k] - Y^{*}[N-k])    0<=k<=N-1
            Xek[i] = (Xk[i] + conjf(Xk[length - i])) / 2;   // Xe[k] = 1/2 * (Y[k] + Y^{*}[N-k])    0<=k<=N-1
            Xok[i] = (Xk[i] - conjf(Xk[length - i])) / (2 * I);  // Xo[k] = /2i * (Y[k] - Y^{*}[N-k])    0<=k<=N-1
        }
        //Xek[0] = (Xk[0] + creal(Xk[0]) - I * cimag(Xk[0])) / 2;   // 由于DFT的周期性，X[0]=X[N]，因此n=0时需要单独考虑
        //Xok[0] = (Xk[0] - (creal(Xk[0]) - I * cimag(Xk[0]))) / (2 * I); 
        Xek[0] = (Xk[0] + conjf(Xk[0])) / 2;   // 由于DFT的周期性，X[0]=X[N]，因此n=0时需要单独考虑
        Xok[0] = (Xk[0] - conjf(Xk[0])) / (2 * I); 
        for(unsigned int j = 1; j <= LENGTH - 1; j++){
            real_Xk[j] = Xek[j] + Xok[j] * ptr_WN_table[j];
            real_Xk[LENGTH * WIDTH - j] = conjf(real_Xk[j]);    // X[N-k] = X^{*}[k]    0<=k<=N-1
        }
        real_Xk[0] = Xek[0] + Xok[0] * ptr_WN_table[0]; //对应的real_Xk[LENGTH * WIDTH]为下一个周期的元素，因此需要单独考虑
        real_Xk[LENGTH] = Xek[0] + Xok[0] * ptr_WN_table[LENGTH];
    }
    else{}
    QueryPerformanceCounter(&end);
    real_us = (double)(end.QuadPart - start.QuadPart) * 1e6 / (double)(freq.QuadPart);
    return;

}

/*
 * 将DFT后的数据写入output.txt文件
 * complex_sign：复数标识符
 * Xk：DFT后序列的数组指针
 */
void write_file(_Bool complex_sign, float complex *Xk, float complex *real_Xk)
{
    float complex const *str_output = Xk;   // 定义const指针，避免数据修改
    unsigned int count;
    FILE *fp;
    if(complex_sign){
        count = LENGTH;
        str_output = Xk;
    }
    else{
        count = LENGTH * WIDTH;
        str_output = real_Xk;
    }
    
    if((fp = fopen(".\\data\\output.txt", "w")) == NULL){
        fprintf(stdout, "can not open file output.txt.\n");
        exit(EXIT_FAILURE);
    }
    while(count > 0){
        float real = creal(*str_output);
        float imag = cimag(*str_output);
        if(imag < epsilon)  // 浮点数判断正负必须设置精度边界
            fprintf(fp, "%f%fi\n", real, imag); // 虚部<0时需要去掉加号，以适配matlab中readmatrix()的格式要求
        else
            fprintf(fp, "%f+%fi\n", real, imag);
        str_output++;
        count--;
    }
    fclose(fp);
    
    return;
}

   
   
   
   