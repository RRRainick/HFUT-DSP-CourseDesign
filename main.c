#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <windows.h>
#include "br_lookup_table.h"
#include "WN_lookup_table.h"

#ifdef LENGTH
#else
    #define LENGTH 2048
    #define WIDTH 2
#endif

LARGE_INTEGER freq;
LARGE_INTEGER start, end;
double resort_ms, reverse_ms, btf_ms;

_Bool read_file(float *ptr_seq, int argc, char *argv[]);
float complex *resort_seq(_Bool complex_sign, float const const (*str_input)[WIDTH], void *xn, float complex *real_xn, float complex *complex_xn);
void reverse_seq(_Bool complex_sign, float complex *xn);
void compute_butterfly(_Bool complex_sign, float complex *const r_xn);
void write_file(_Bool complex_sign, float complex *Xk);

int main(int argc, char *argv[])
{
    float input_seq[LENGTH][WIDTH]; // 默认零矩阵，满足序列自动补零
    _Bool complex_sign = read_file(input_seq, argc, argv);
    float const (*str_input)[WIDTH] = input_seq;
    void *xn;   // 输入序列指针
    float complex real_xn[LENGTH*WIDTH]; float complex complex_xn[LENGTH]; // 由于C语言不支持函数重载，因此将real_xn和complex_xn均定义为float complex类型
    xn = resort_seq(complex_sign, str_input, xn, real_xn, complex_xn);
   
    reverse_seq(complex_sign, xn);

    compute_butterfly(complex_sign, xn);
    
    write_file(complex_sign, xn);

    printf("resort cost %fms, reverse cost %fms, compute butterfly cost %fms.\nduration is %fms.\n", resort_ms, reverse_ms, btf_ms, resort_ms+reverse_ms+btf_ms);
    return 0;
}


_Bool read_file(float *ptr_seq, int argc, char *argv[])
{
    unsigned char tmp_array[32]; // 用于储存str2float的数字，float类型在x64上占64bit
    char ch;
    unsigned int digit = 0;
    unsigned int count = LENGTH * WIDTH; // 下标计数器，避免数组越界
    double i_decp = 0.0;
    double tmp_num = 0.0;
    _Bool plus_sign = 1;
    _Bool complex_sign = 0;
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
    if(*(argv[2]) == 'c'){
        complex_sign = 1;
    }
    else if(*(argv[2]) == 'r'){
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
    while((ch = getc(fp)) != EOF){
        if(count == 0)  //避免数组越界
            break;
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

float complex *resort_seq(_Bool complex_sign, float const (*str_input)[WIDTH], void *xn, float complex *real_xn, float complex *complex_xn)
{
    QueryPerformanceFrequency(&freq);
    QueryPerformanceCounter(&start);
    float complex *output_xn;
    if(complex_sign){
        for(int i = 0; i < LENGTH; i++){
            complex_xn[i] = (*str_input)[0] + (*str_input)[1] * I;
            str_input++;
        }
        output_xn = complex_xn;
    }
    else{
        for(int j = 0; j < LENGTH * WIDTH; j = j + 2){
           real_xn[j] = (*str_input)[0];
           real_xn[j+1] = (*str_input)[1];
           str_input++; 
        }
        output_xn = real_xn;
    }
    QueryPerformanceCounter(&end);
    resort_ms = (double)(end.QuadPart - start.QuadPart) * 1e6 / (double)(freq.QuadPart);
    return output_xn;
} 

void reverse_seq(_Bool complex_sign, float complex *xn)
{
    QueryPerformanceFrequency(&freq);
    QueryPerformanceCounter(&start);
    // log_{2}^{2048}=11, log_{2}^{4096}=12
    unsigned int count;
    unsigned int const *ptr_br_table;
    _Bool i_arr_sign[LENGTH * WIDTH] = {0}; // 防止重复交换
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
    reverse_ms = (double)(end.QuadPart - start.QuadPart) * 1e6 / (double)(freq.QuadPart);
    
}

void compute_butterfly(_Bool complex_sign, float complex * const r_xn)
{
    QueryPerformanceFrequency(&freq);
    QueryPerformanceCounter(&start);
    unsigned int length;
    float complex const *r_ptr_WN_table;
    float complex const *o_ptr_WN_table;
    float complex const *ptr_WN_table;

    if(complex_sign){
        length = LENGTH;
        r_ptr_WN_table = complex_WN_table;
        o_ptr_WN_table = complex_WN_table;
        ptr_WN_table = complex_WN_table;
    }
    else{
        length = LENGTH * WIDTH;
        r_ptr_WN_table = real_WN_table;
        o_ptr_WN_table = real_WN_table;
        ptr_WN_table = real_WN_table;
    }

    unsigned int max_stage = (unsigned int)log2((double)length);
    unsigned int j_max = length / 2; //unsigned int j_max = (unsigned int)pow(2.0, (double)(max_stage-1)); 
    unsigned int k_max = 1;
    float complex *ptr_r_xn = r_xn;
    //unsigned int count = 0;
    int gap;
    int gap_WN;
    for(unsigned int i = 1; i <= max_stage; i++){
        for(unsigned int j = 1; j <= j_max; j++){
            for(unsigned int k = 1; k <= k_max; k++){
                float complex tmp1, tmp2;
                tmp1 = ptr_r_xn[0] + ptr_WN_table[0] * ptr_r_xn[k_max];
                tmp2 = ptr_r_xn[0] + ptr_WN_table[k_max] * ptr_r_xn[k_max];
                ptr_r_xn[0] = tmp1; ptr_r_xn[k_max] = tmp2;
                ptr_r_xn++;
                gap = ptr_r_xn - r_xn;
                ptr_WN_table++;
                gap_WN = ptr_WN_table - r_ptr_WN_table;
                //count++;
            }
            ptr_r_xn = r_xn + j * k_max * 2;
            ptr_WN_table = r_ptr_WN_table + j * k_max * 2;
            gap = ptr_r_xn - r_xn;
            gap_WN = ptr_WN_table - r_ptr_WN_table;
        }
        ptr_r_xn = r_xn;
        gap = ptr_r_xn - r_xn;
        r_ptr_WN_table += length;
        ptr_WN_table = r_ptr_WN_table;
        gap_WN = ptr_WN_table - o_ptr_WN_table;
        j_max /= 2;
        k_max *= 2;
    }
    QueryPerformanceCounter(&end);
    btf_ms = (double)(end.QuadPart - start.QuadPart) * 1e6 / (double)(freq.QuadPart);
}

void write_file(_Bool complex_sign, float complex *Xk)
{
    float complex const *str_output = Xk;
    unsigned int count;
    FILE *fp;
    if(complex_sign)
        count = LENGTH;
    else
        count = LENGTH * WIDTH;
    
    if((fp = fopen("D:\\Rainick\\Code_Repo\\Matlab\\output.txt", "w")) == NULL){
        fprintf(stdout, "can not open file D:\\Rainick\\Code_Repo\\Matlab\\output.txt.\n");
        exit(EXIT_FAILURE);
    }
    while(count > 0){
        float real = creal(*str_output);
        float imag = cimag(*str_output);
        if(imag < 0)
            fprintf(fp, "%f%fi\n", real, imag);
        else
            fprintf(fp, "%f+%fi\n", real, imag);
        str_output++;
        count--;
    }
    fclose(fp);
}

   
   
   
   