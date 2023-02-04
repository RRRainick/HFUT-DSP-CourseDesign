#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<complex.h>
#include"br_lookup_table.h"

#ifdef LENGTH
#else
    LENGTH = 2048;
    WIDTH = 2;
#endif

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
                printf("%f\n", tmp_num);
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


void reverse_seq(_Bool complex_sign, float complex *xn){
    // log_{2}^{2048}=11, log_{2}^{4096}=12
    unsigned int count;
    unsigned int const *ptr_br_table;
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
        complex tmp = xn[k];
        xn[k] = xn[index]; xn[index] = tmp; 
    }


}

int main(int argc, char *argv[])
{
    float input_seq[LENGTH][WIDTH]; // 默认零矩阵，满足序列自动补零
    _Bool complex_sign = read_file(input_seq, argc, argv);
    float const (*str_input)[WIDTH] = input_seq;
    void *xn;   // 输入序列指针
    void *Xk;
    float complex real_xn[LENGTH*WIDTH]; float complex complex_xn[LENGTH]; // 由于C语言不支持函数重载，因此将real_xn和complex_xn均定义为float complex类型
    float real_Xk[LENGTH*WIDTH]; float complex complex_Xk[LENGTH];

    if(complex_sign){
        for(int i = 0; i < LENGTH; i++){
            complex_xn[i] = (*str_input)[0] + (*str_input)[1] * I;
            str_input++;
        }
        xn = complex_xn;
    }
    else{
        for(int j = 0; j < LENGTH * WIDTH; j = j + 2){
           real_xn[j] = (*str_input)[0];
           real_xn[j+1] = (*str_input)[1];
           str_input++; 
        }
        xn = real_xn;
    }
    
    reverse_seq(complex_sign, xn);
    return 0;
}

   
   
   
   