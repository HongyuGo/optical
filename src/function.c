#include "function.h"
#include "commom.h"
#include "matrix.h"
#include "params.h"
#include <math.h>

double h_response(double t, double TLL, double SS) { return erf(t / SS / TLL); }

double readback(double t, double jitter, Matrix *d, double SS, double TT,
                double TLL) {
    int len = d->column;
    double rs = 0;
    double tmp = 0;
    for (int i = 0; i < len; i++) {
        tmp = d->data[0][i] * h_response(t - (i + 1) * TT + jitter, TLL, SS);
        rs = rs + tmp;
    }
    return rs;
}

Matrix **gen_firtaps_v2(Matrix *random, Matrix *sampled,
                        MATRIX_TYPE *gpr_coeff_data, int fir_len,
                        char constraint, char method) {
    int gpr_len = 5;
    // int i,j;
    // int datalength = random->column;
    Matrix *I_matrix = M_Zeros(gpr_len, 1);
    switch (constraint) {
    case '1':
        break;
    case 'c':
        M_value_one(I_matrix, gpr_len / 2 + 1, 1, 1.0);
        // M_print(I_matrix,"I_matrix");
        break;
    case '2':
        break;
    default:
        break;
    }
    //-----------R_matrix.. Calculate the autocorrelation matrix of the sampled
    // data---------------
    // Size:(2K+1)x(2K+1)
    int K = (fir_len - 1) / 2;
    // printf("%d",K);
    Matrix *R_matrix = auto_corr(sampled, -K, K);
    // R_matrix = M_Inverse(R_matrix);
    // M_print(R_matrix,"R_matrix");
    // for(i = 1; i < K; i++){
    //     for(j = 0; j < K; j++){
    //         M_value_one(R_matrix,i + 1,j + 1,M_get_one(R_matrix,1,abs(j - i)
    //         + 1));
    //     }
    // }

    //-----------A_matrix.. Caculate the autocorrelation matrix of the random
    // data----------------- Size: LxL
    int L = gpr_len;
    Matrix *A_matrix = auto_corr(random, 0, L - 1);
    // A_matrix = M_Inverse(A_matrix);
    // M_print(A_matrix,"A_matrix");

    //-----------T_matrix.. Caculate the cross-correlation matrix of the sampled
    // and the random----------- data. Size: (2K+1)xL
    Matrix *T_matrix = cross_corr(sampled, random, fir_len, gpr_len);
    // M_print(T_matrix,"T_matrix");

    MATRIX_TYPE lagrange =
        Caculate_lagrange(R_matrix, A_matrix, T_matrix, I_matrix);
    // printf("lagrange = %lf\n",lagrange);
    // M_print(M_lagrange,"data");

    Matrix *gpr_coeff =
        Caculate_gpr_coeff(R_matrix, A_matrix, T_matrix, I_matrix, lagrange);
    // M_print(gpr_coeff,"gpr_coeff");

    Matrix *fir_coeff = Caculate_fir_coeff(R_matrix, T_matrix, gpr_coeff);
    // M_print(fir_coeff,"fir_coeff");

    Matrix **return_back = NULL;
    return_back = (Matrix **)malloc(sizeof(Matrix *) * 2);
    return_back[1] = gpr_coeff;
    return_back[0] = fir_coeff;

    M_free(T_matrix);
    M_free(A_matrix);
    M_free(I_matrix);
    M_free(R_matrix);
    return return_back;
}
Matrix *auto_corr(Matrix *x, int bot, int top) {
    int L = top - bot + 1;
    int i, j, k, end = x->column;
    MATRIX_TYPE sum = 0.0, prod = 0.0;
    Matrix *ac = M_Zeros(L, L);
    for (i = bot; i <= top; i++) {
        for (j = bot; j <= top; j++) {
            if (end + bot < top + 1)
                break;
            else {
                for (k = top + 1; k <= end + bot; k++) {
                    sum += x->data[0][k - i - 1] * x->data[0][k - j - 1];
                }
                prod = sum / (end + bot - (top + 1) + 1);
                sum = 0.0;
                M_value_one(ac, i - bot + 1, j - bot + 1, prod);
            }
        }
    }
    return ac;
}
Matrix *cross_corr(Matrix *x, Matrix *y, int K, int L) {
    int hK = (K - 1) / 2;
    int i, j, k;
    MATRIX_TYPE prod = 0.0, sum = 0.0;
    int end = x->column;
    Matrix *xc = M_Zeros(K, L);
    int tmp = hK > L - 1 ? hK : L - 1;
    for (i = -hK; i <= hK; i++) {
        for (j = 0; j <= L - 1; j++) {
            if (1 + tmp > end - tmp)
                break;
            else {
                for (k = 1 + tmp; k <= end - tmp; k++) {
                    sum += x->data[0][k - i - 1] * y->data[0][k - j - 1];
                }
                prod = sum / (end - tmp - (1 + tmp) + 1);
                sum = 0.0;
                M_value_one(xc, i + hK + 1, j + 1, prod);
            }
        }
    }
    return xc;
}
MATRIX_TYPE Caculate_lagrange(Matrix *R_matrix, Matrix *A_matrix,
                              Matrix *T_matrix, Matrix *I_matrix) {
    // M_print(R_matrix,"R_matrix");
    Matrix *R_inv = M_Inverse(R_matrix);
    // M_print(R_inv,"R_inv");

    Matrix *T_T = M_T(T_matrix);
    Matrix *TR = M_mul(T_T, R_inv);
    M_free(T_T);
    M_free(R_inv);

    Matrix *TRT = M_mul(TR, T_matrix);
    Matrix *A_TRT = M_add_sub(1.0, A_matrix, 1.0, TRT);
    M_free(TRT);
    M_free(TR);

    Matrix *A_TRT_inv = M_Inverse(A_TRT);
    Matrix *I_T = M_T(I_matrix);
    Matrix *I_A_TRT = M_mul(I_T, A_TRT_inv);
    M_free(A_TRT_inv);
    M_free(A_TRT);
    M_free(I_T);

    Matrix *M_lagrange = M_mul(I_A_TRT, I_matrix);
    M_free(I_A_TRT);

    MATRIX_TYPE lagrange = 1.0 / M_lagrange->data[0][0];
    M_free(M_lagrange);
    return lagrange;
}
Matrix *Caculate_gpr_coeff(Matrix *R_matrix, Matrix *A_matrix, Matrix *T_matrix,
                           Matrix *I_matrix, MATRIX_TYPE lagrange) {
    Matrix *R_inv = M_Inverse(R_matrix);
    // M_print(R_inv,"R_inv");
    Matrix *T_T = M_T(T_matrix);
    Matrix *TR = M_mul(T_T, R_inv);
    M_free(T_T);
    M_free(R_inv);

    Matrix *TRT = M_mul(TR, T_matrix);
    Matrix *A_TRT = M_add_sub(1.0, A_matrix, 1.0, TRT);
    // M_print(A_TRT,"A_TRT");
    M_free(TRT);
    M_free(TR);

    Matrix *A_TRT_inv = M_Inverse(A_TRT);
    // M_print(A_TRT_inv,"A_TRT_inv");
    Matrix *gpr_coeff = M_mul(A_TRT_inv, I_matrix);
    M_free(A_TRT);
    M_free(A_TRT_inv);

    gpr_coeff = M_nummul(gpr_coeff, lagrange);
    return gpr_coeff;
}
Matrix *Caculate_fir_coeff(Matrix *R_matrix, Matrix *T_matrix,
                           Matrix *gpr_coeff) {
    Matrix *R_inv = M_Inverse(R_matrix);
    Matrix *TR = M_mul(R_inv, T_matrix);
    M_free(R_inv);
    Matrix *TRgpr = M_mul(TR, gpr_coeff);
    M_free(TR);
    return TRgpr;
}

int compare(int _dataget_col, int _array_col,MATRIX_TYPE (*dataget)[_dataget_col],MATRIX_TYPE (*array)[_array_col],int row){
    int i;
    for(i = 1; i <= 4; i++){
        if(dataget[0][0] == array[row - 1][i * 2 -2] && dataget[0][1] == array[row - 1][i * 2 - 1]){
            return i;
        }
    }
    return -1;
}

void encode_17pp(int _src_col, int _des_col, MATRIX_TYPE (*_src)[_src_col], MATRIX_TYPE (*_des)[_des_col]){
    static MATRIX_TYPE code1[4] = {0,0,0,2};
    static MATRIX_TYPE code2[4] = {0,0,1,2};
    static MATRIX_TYPE code3[4] = {0,1,0,2};
    static MATRIX_TYPE code4[7] = {0,1,0,1,0,0,2};
    static MATRIX_TYPE code5[7] = {0,1,0,0,0,0,2};
    static MATRIX_TYPE code6[7] = {0,0,0,1,0,0,2};
    static MATRIX_TYPE code7[10] = {0,0,0,1,0,0,1,0,0,2};
    static MATRIX_TYPE code8[10] = {0,0,0,1,0,0,0,0,0,2};
    static MATRIX_TYPE code9[10] = {0,1,0,1,0,0,1,0,0,2};
    static MATRIX_TYPE code10[10] = {0,1,0,1,0,0,0,0,0,2};
    static MATRIX_TYPE code11[10] = {0,0,1,0,0,0,0,0,0,2};
    static MATRIX_TYPE code12[13] = {0,0,0,1,0,0,1,0,0,1,0,0,2};
    static MATRIX_TYPE code13[13] = {0,1,0,1,0,0,1,0,0,1,0,0,2};
    static MATRIX_TYPE code14[4] ={1,0,1,2};

    MATRIX_TYPE *out_array[15];
    out_array[1] = code1;
    out_array[2] = code2;
    out_array[3] = code3;
    out_array[4] = code4;
    out_array[5] = code5;
    out_array[6] = code6;
    out_array[7] = code7;
    out_array[8] = code8;
    out_array[9] = code9;
    out_array[10] = code10;
    out_array[11] = code11;
    out_array[12] = code12;
    out_array[13] = code13;
    out_array[14] = code14;


    int r = 0;
    int NRZIUserLen = _src_col;
    int NRZIUserLen_keep = _src_col;
    int cur = 0;
    MATRIX_TYPE array[4][8] = {{0,0,1,1,1,0,0,1},
                               {0,0,1,1,1,0,0,1},
                               {0,0,1,1,1,0,0,1},
                               {0,0,1,1,1,0,0,1}};
    int out_index = 0;
    int last_channel_bit = 0;
    int vssf = 0;
    int special_cur = 0;
    int special_index = 0;
    int flag = 1;
    int index_new_turn = 0;
    int row = 1;
    int flag_terminate = 0;
    MATRIX_TYPE dataget[1][2];
    int index = 0;
    int des_cur = 0;
    int very_special_state_index = 0;
    while(NRZIUserLen > 0){
        flag = 1;
        index_new_turn = 0;
        row = 1;
        flag_terminate = 0;
        while(r < 3 && flag == 1 && flag_terminate == 0){
            r = r + 1;
            dataget[0][0] = _src[0][cur];
            dataget[0][1] = _src[0][cur + 1];
            cur = cur + 2;
            NRZIUserLen = NRZIUserLen - 2;
            index = compare(2,8,dataget,array,row);
            switch(index){
                case 1:
                    vssf = vssf != 4? 0 : 4;
                    index_new_turn = index_new_turn + 3;
                    if(r == 3){
                        out_index = 10;
                        if(NRZIUserLen_keep - cur >= 2){
                            dataget[0][0] = _src[0][cur];
                            dataget[0][1] = _src[0][cur + 1];
                            index = compare(2,8,dataget,array,4);
                            if(index == 1){
                                out_index = 13;
                                cur = cur + 2;
                                NRZIUserLen = NRZIUserLen - 2;
                            }
                        }
                    }
                    row = row + 1;
                    break;

                case 2: //11
                    out_index  = 1 + index_new_turn;
                    if(vssf == 0 || vssf == 2){
                        vssf = vssf + 1;
                    }else{
                        if(vssf != 4){
                            vssf = 0;
                        }
                    }
                    if(r == 1 && last_channel_bit == 0){
                        out_index = 14;
                    }
                    flag = 0;
                    break;

                case 3: //10
                    out_index = 2 + index_new_turn;
                    vssf = vssf != 4? 0: 4;
                    if(r == 3 && NRZIUserLen_keep - cur  >= 2){
                        special_cur = cur;
                        dataget[0][0] = _src[0][special_cur];
                        dataget[0][1] = _src[0][special_cur + 1];
                        special_index = compare(2,8,dataget,array,4);
                        if(special_index == 1){
                            out_index = 12;
                            cur = special_cur + 2;
                            NRZIUserLen = NRZIUserLen - 2;
                        }
                    }
                    flag = 0;
                    break;
                
                case 4: //01
                    out_index = 3 + index_new_turn;
                    if(vssf == 1){
                        vssf = vssf + 1;
                    }else{
                        if(vssf != 4){
                            vssf = 0;
                        }
                    }
                    flag = 0;
                    break;
            }


        }
        MATRIX_TYPE * value_cur = out_array[out_index];
        while(*value_cur != 2){
            _des[0][des_cur++] = *value_cur;
            value_cur++;
        }
        if(vssf == 4){
            if(out_array[out_index][0] == 0 && out_array[out_index][1] == 1 && out_array[out_index][2] == 0){
                MATRIX_TYPE *value_cur2 = out_array[11];
                while(*value_cur2 != 2){
                    _des[0][very_special_state_index++] = *value_cur2;
                    value_cur2++;
                }
            }   
            vssf = 0;
        }
        if(vssf == 3){
            very_special_state_index = des_cur - 9;
            vssf = 4;
        }
        last_channel_bit = _des[0][des_cur - 1];
        out_index = 0;
        r = 0;
    } 
}