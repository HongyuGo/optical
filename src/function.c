#include "function.h"
#include <math.h>
#include "matrix.h"
#include "commom.h"

double h_response(double t, double TL, double S) { return erf(t / S / TL); }

double readback(double t, double jitter, Matrix* d, double S, double T, double TL) {
    int len = d->column;
    double rs = 0;
    double tmp = 0;
    for (int i = 0; i < len; i++) {
        tmp = d->data[i] * h_response(t - (i + 1) * T + jitter, TL, S);
        rs = rs + tmp;
    }
    return rs;
}

void gen_firtaps_v2(Matrix *random, Matrix *sample, MATRIX_TYPE* gpr_coeff, int fir_len, char constraint,char method)
{
    int gpr_len = 5;
    // int i,j;
    //int datalength = random->column;
    Matrix *I_matrix = M_Zeros(5,1);
    switch(constraint){
        case '1':
                break;
        case 'c':
                M_value_one(I_matrix,gpr_len/2 + 1,1,1.0); 
                break;
        case '2':
                break;
        default:
                break;
    }
    int K = (fir_len - 1) / 2;
    Matrix * R_matrix = auto_corr(sample,-K,K);
    R_matrix = M_Inverse(R_matrix);
    //M_print(R_matrix,"R_matrix");
    // for(i = 1; i < K; i++){
    //     for(j = 0; j < K; j++){
    //         M_value_one(R_matrix,i + 1,j + 1,M_get_one(R_matrix,1,abs(j - i) + 1));
    //     }
    // }
    int L = gpr_len;
    Matrix * A_matrix = auto_corr(random,0,L - 1);
    //A_matrix = M_Inverse(A_matrix);
    //M_print(A_matrix,"A_matrix");


    free(A_matrix);
    free(I_matrix);
    free(R_matrix);
}
Matrix * auto_corr(Matrix * x, int bot , int top){
    int L = top - bot + 1;
    int i, j,k,end = x->column;
    MATRIX_TYPE sum = 0.0,prod = 0.0;
    Matrix * ac = M_Zeros(L,L);
    for(i = bot ; i <= top; i++){
        for(j = bot; j <= top; j++){
           if(end + bot < top + 1)
                break;
           else{
            for(k = top + 1; k <= end + bot; k++){
                sum += x->data[k - i - 1] * x->data[k - j - 1];
            }
            prod = sum / (end + bot - i - (top + 1 -i) + 1);
            sum = 0.0;
            M_value_one(ac,i-bot+1,j-bot+1,prod);
           } 
        }
    }
    return ac;
}
#if 0
void auto_corr(Matrix * A, Matrix * B,int shift , Matrix * output, int row, int col)
{
    int len_x = A->column;
    int len_y = B->column;
    int len_output = output->column;
    int start , finish ,k;
    if(shift >= 0)
        start = shift;
    else
        start = 0;
    
    if(len_x - 1 <= shift + len_y - 1)
        finish = len_x  - 1;
    else
        finish = shift + len_y - 1;
    if(start > finish)return;
    int location = (row - 1) * len_output + (col - 1); 
    for(k = start; k <= finish; k++){
        output->data[location] += A->data[k] * B->data[k - shift];
    }
}
#endif