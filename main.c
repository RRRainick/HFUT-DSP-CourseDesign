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

LARGE_INTEGER freq;
LARGE_INTEGER start, end;
double resort_ms, reverse_ms, btf_ms, duration_ms;

_Bool read_file(float *ptr_seq, int argc, char *argv[]);
float complex *rollover_seq(_Bool complex_sign, float const (*ptr_input)[WIDTH], float complex *real_xn, float complex *complex_xn);
void bitreorder_seq(_Bool complex_sign, float complex *xn);
void compute_butterfly(_Bool complex_sign, float complex *const r_xn);
void write_file(_Bool complex_sign, float complex *Xk);

int main(int argc, char *argv[])
{
    float input_seq[LENGTH][WIDTH] = {0}; // 默认零矩阵，实现序列自动补零
    _Bool complex_sign = read_file(input_seq, argc, argv);
    float const (*ptr_input)[WIDTH] = input_seq;
    float complex real_xn[LENGTH*WIDTH]; float complex complex_xn[LENGTH]; // 由于C语言不支持函数重载，因此将real_xn和complex_xn均定义为float complex类型
    float complex *xn = rollover_seq(complex_sign, ptr_input, real_xn, complex_xn);
   
    QueryPerformanceFrequency(&freq);   // 获取CPU时钟周期
    QueryPerformanceCounter(&start);    // 获取当前经过的CPU时钟次数

    bitreorder_seq(complex_sign, xn);

    compute_butterfly(complex_sign, xn);
    
    QueryPerformanceCounter(&end);  // 获取当前经过的CPU时钟次数
    duration_ms = (double)(end.QuadPart - start.QuadPart) * 1e3 / (double)(freq.QuadPart);    // 时钟次数作差并除以时钟频率得到经过的时间，以ms为单位的

    write_file(complex_sign, xn);
    printf("resort cost %fms, reverse cost %fms, compute butterfly cost %fms.\nduration is %fms.\n", resort_ms, reverse_ms, btf_ms, duration_ms);

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
    double i_decp = 0.0;    // 表示左起小数点的位置
    double tmp_num = 0.0;   // 临时储存数字和
    _Bool plus_sign = 1;    // 正数指示变量（主要用于判断负数）
    _Bool complex_sign = 0; // 复数指示变量
    FILE *fp;
    /*
     * 打开储存输入序列的txt文件
     * 对不满足要求的参数报错并退出
     */
    if(argc != 3){
        printf("usage of %s is incorrect.\n", argv[0]);
        exit(EXIT_FAILURE);
    }
    if((fp = fopen(argv[1], "r")) == NULL){
        printf("can not open file %s.\n", argv[1]);
        exit(EXIT_FAILURE);
    }
    /*
     * 根据命令参数判别实序列和复序列
     */
    if(*(argv[2]) == 'c' || *(argv[2]) == 'C'){
        complex_sign = 1;
    }
    else if(*(argv[2]) == 'r' || *(argv[2]) == 'R'){
        complex_sign = 0;
    }
    else{
        printf("please indicate sequence format in %s.\n", argv[1]);
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
            case ' ':
            case '\n':
                if(i_decp == 0){
                    for(int i = 0; i < digit; i++){
                            tmp_num += tmp_array[i] * pow(10.0, digit-1.0-i);
                        }
                }
                else{
                    for(int i = 0; i < digit; i++){
                        tmp_num += tmp_array[i] * pow(10.0, i_decp-1.0-i);
                    }
                }
                if(!plus_sign)
                    tmp_num = (-1.0)*tmp_num;
                *ptr_seq = tmp_num;
                ptr_seq++; count--;
                // printf("%f\n", tmp_num);
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
 * complex_sign：复数标识符
 * ptr_input：二维数组指针
 * real_xn：转存的实数数组指针
 * complex_xn：转存的复数数组指针
 * 返回值：对数组元素处理的对应数组指针
 */
float complex *rollover_seq(_Bool complex_sign, float const (*ptr_input)[WIDTH], float complex *real_xn, float complex *complex_xn)
{
    QueryPerformanceFrequency(&freq);   // 获取CPU时钟周期
    QueryPerformanceCounter(&start);    // 获取当前经过的CPU时钟次数
    float complex *output_xn;
    if(complex_sign){
        for(int i = 0; i < LENGTH; i++){
            complex_xn[i] = (*ptr_input)[0] + (*ptr_input)[1] * I;
            ptr_input++;
        }
        output_xn = complex_xn;
    }
    else{
        for(int j = 0; j < LENGTH * WIDTH; j = j + 2){
           real_xn[j] = (*ptr_input)[0];
           real_xn[j+1] = (*ptr_input)[1];
           ptr_input++; 
        }
        output_xn = real_xn;
    }

    QueryPerformanceCounter(&end);  // 获取当前经过的CPU时钟次数
    resort_ms = (double)(end.QuadPart - start.QuadPart) * 1e3 / (double)(freq.QuadPart);    // 时钟次数作差并处以时钟频率得到经过的时间，以ms为单位的
    
    return output_xn;
} 

/*
 * 将信号序列进行bitreverse排序
 * complex_sign：复数标识符
 * xn：信号序列的数组指针
*/
void bitreorder_seq(_Bool complex_sign, float complex *xn)
{
    QueryPerformanceFrequency(&freq);
    QueryPerformanceCounter(&start);
    // log_{2}^{2048}=11, log_{2}^{4096}=12
    unsigned int count; // 排序次数计数器
    unsigned int const *ptr_br_table;   // 查找表指针
    _Bool i_arr_sign[LENGTH * WIDTH] = {0}; // inverse_array_sign, 标识已经交换过的元素，防止重复交换
    if(complex_sign){
        count = LENGTH;
        ptr_br_table = complex_br_table;
    }
    else{
        count = LENGTH * WIDTH;
        ptr_br_table = real_br_table;
    }
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
    reverse_ms = (double)(end.QuadPart - start.QuadPart) * 1e3 / (double)(freq.QuadPart);

    return;    
}

/*
 * 对bitreverse后的序列进行蝶形运算
 * complex_sign：复数标识符
 * r_xn：reverse_xn，已bitreverse后的序列的数组指针
 */
void compute_butterfly(_Bool complex_sign, float complex * const r_xn)
{
    QueryPerformanceFrequency(&freq);
    QueryPerformanceCounter(&start);
    unsigned int length;
    float complex const *r_ptr_WN_table;    // reference_ptr_WN_table，指向每个stage最初的旋转因子
    float complex const *ptr_WN_table;      // 旋转因子查找表指针
    float complex *ptr_r_xn = r_xn; 

    if(complex_sign){
        length = LENGTH;
        r_ptr_WN_table = complex_WN_table;
        ptr_WN_table = complex_WN_table;
    }
    else{
        length = LENGTH * WIDTH;
        r_ptr_WN_table = real_WN_table;
        ptr_WN_table = real_WN_table;
    }

    unsigned int max_stage = (unsigned int)log2((double)length);    //计算所需的stage个数
    unsigned int j_max = length / 2;    // 每个stage内循环次数
    //unsigned int j_max = (unsigned int)pow(2.0, (double)(max_stage-1)); 
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
    btf_ms = (double)(end.QuadPart - start.QuadPart) * 1e3 / (double)(freq.QuadPart);

    return;
}

/*
 * 将DFT后的数据写入output.txt文件
 * complex_sign：复数标识符
 * Xk：DFT后序列的数组指针
 */
void write_file(_Bool complex_sign, float complex *Xk)
{
    float complex const *str_output = Xk;   // 定义const指针，避免数据修改
    unsigned int count;
    FILE *fp;
    if(complex_sign)
        count = LENGTH;
    else
        count = LENGTH * WIDTH;
    
    if((fp = fopen(".\\data\\output.txt", "w")) == NULL){
        fprintf(stdout, "can not open file output.txt.\n");
        exit(EXIT_FAILURE);
    }
    while(count > 0){
        float real = creal(*str_output);
        float imag = cimag(*str_output);
        if(imag < 0)
            fprintf(fp, "%f%fi\n", real, imag); // 虚部<0时需要去掉加号，以适配matlab中readmatrix()的格式要求
        else
            fprintf(fp, "%f+%fi\n", real, imag);
        str_output++;
        count--;
    }
    fclose(fp);
    
    return;
}

   
   
   
   