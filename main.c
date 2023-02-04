#include<stdio.h>
#include<stdlib.h>
#include<math.h>
// #include<complex.h>
#define LENGTH 2048
#define WIDTH 2

_Bool read_file(float *ptr_sequence, int argc, char *argv[])
{
    unsigned char tmp_array[32]; // 用于储存str2float的数字，float类型在x64上占64bit
    char ch;
    int digit = 0;
    int count = LENGTH * WIDTH; // 下标计数器，避免数组越界
    double index_decp = 0.0;
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
                index_decp = digit;
                break;
            case '-':
                plus_sign = 0;
                break;
            case ' ':
            case '\n':
                if(index_decp == 0){
                    for(int i = 0; i < digit; i++){
                            tmp_num += tmp_array[i] * pow(10.0, digit-1.0-i);
                        }
                }
                else{
                    for(int i = 0; i < digit; i++){
                        tmp_num += tmp_array[i] * pow(10.0, index_decp-1.0-i);
                    }
                }
                if(!plus_sign)
                    tmp_num = (-1.0)*tmp_num;
                *ptr_sequence = tmp_num;
                ptr_sequence++; count--;
                printf("%f\n", tmp_num);
                digit = 0; tmp_num = 0.0; plus_sign = 1; index_decp = 0.0;  // 变量复位
                break;
            default:
                tmp_array[digit++] = ch - '0'; // char2short
                break;
        }
    }
    fclose(fp);
    return complex_sign;

}
void reverse_bit(){

}

int main(int argc, char *argv[])
{
    float input_sequence[LENGTH][WIDTH]; // 默认零矩阵，满足序列自动补零
    _Bool complex_sign = read_file(input_sequence, argc, argv);
    /*
    float (*str_input)[WIDTH] = input_sequence;
    float complex c_sequence[LENGTH]; float r_sequence[LENGTH];
    if(complex_sign){
        for(int i = 0; i < LENGTH; i++){
            c_sequence[i] = (*str_input)[0] + (*str_input)[1] * I;
            str_input++;
        }
    } 
    else{
        for(int j = 0; j < LENGTH * WIDTH; j = j + 2){
           r_sequence[j] = (*str_input)[0];
           r_sequence[j+1] = (*str_input)[1];
           str_input++; 
        }
    }
    */
    
    return 0;
}

   
   
   
   