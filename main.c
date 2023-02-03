#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#define LENGTH 2048
#define WIDTH 2

int main(int argc, char *argv[])
{
    float input_sequence[LENGTH][WIDTH]; // 默认零矩阵，满足序列自动补零
    unsigned char tmp_array[32]; // float类型在x64上占64bit
    char ch;
    int digit = 0;
    double index_decp = 0.0;
    int i_input = 0; int j_input = 0;
    _Bool plus_sign = 1;
    double tmp_num = 0.0;
    FILE *fp;
    /*
     * open input sequence file
     * store in array input_sequence[LENGTH][WIDTH]
     */
    if(argc != 2){
        printf("usage of %s is incorrect.\n", argv[0]);
        exit(EXIT_FAILURE);
    }
    if((fp = fopen(argv[1], "r")) == NULL){
        printf("can not open file %s.\n", argv[1]);
        exit(EXIT_FAILURE);
    }
    while((ch = getc(fp)) != EOF){
        // sign check
        if(ch == '-'){
            plus_sign = 0;
            continue;
        }
        else if(ch == '.'){
            index_decp = digit;
            continue;
        }
        else if(ch == ' ' || ch == '\n'){
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
            if(j_input == 2){
                i_input++;
                j_input = 0;
            }
            input_sequence[i_input][j_input] = tmp_num;
            printf("%f\n", tmp_num);
            digit = 0; tmp_num = 0.0; plus_sign = 1; index_decp = 0.0; j_input++; // variable reset
        }
        else{
           tmp_array[digit++] = ch - '0'; // char2short
        }
    }
    fclose(fp);


    return 0;
}

   
   
   
   