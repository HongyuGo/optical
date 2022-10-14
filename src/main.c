#include "main.h"
#include <math.h>
#include "commom.h"

int main() {
    MATRIX_TYPE gpr_target[5] ={1, 2,2,2, 1};
    //int gpr_length = sizeof(gpr_target)/sizeof(MATRIX_TYPE);
    int fir_length = 13;
    // int pr2_target[5] = {1,2,2,2,1};
    // int pr2_length = sizeof(pr2_target)/sizeof(int);
    int SNR_len = sizeof(SNR)/sizeof(int);
    Matrix* BER1 = M_Zeros(SNR_len,1);
    Matrix* BER2 = M_Zeros(SNR_len, 1);
    Matrix* BER3 = M_Zeros(SNR_len, 1);
    
    for(int isnr = 1; isnr <= SNR_len; isnr++)
    {
        int numBits = 0;
        // int numErrs = 0;
        // int numERRs1 = 0;
        // int numErrs2 = 0;
        int numErrs3 = 0;
        if(numErrs3 < maxErrs && numBits < maxBits)
        {
            #if 0
            //------------PRBS data---------------------
            srand(3);
            MATRIX_TYPE rand_data[SectorLength + 1];
            for(int i = 0; i < SectorLength; i++)
                rand_data[i] = rand()%2;
            Matrix *ChannelBits = Matrix_gen(1,SectorLength,rand_data);
            M_print(ChannelBits,"ChannelBits");
            //M_print(ChannelBits);
            #endif
            #define Test_len 30
            MATRIX_TYPE test_data[Test_len] = {0,1,0,0,1,0,1,1,1,0,0,0,1,1,0,0,1,1,1,1,0,1,0,1,1,1,0,1,1,0};
            //MATRIX_TYPE test_data[15] = {0,1,0,0,1,0,1,1,1,0,0,0,1,0,1};
            Matrix *ChannelBits = Matrix_gen(1,Test_len,test_data);
            //M_print(ChannelBits,"CH");
            Matrix *tempInputPad = M_full(ChannelBits,0,0,1,0,0);
            //M_print(tempInputPad);

            //RLL encoder 
            Matrix * codedwords = Matrix_copy(tempInputPad);
            int codedlen = codedwords->column;
            int CodedBitsLength = codedlen - 1;

            Matrix *ak = M_full(codedwords,0,0,KWinLen,0,0);
            ak = M_numsub(M_numul(ak,2.0),1.0);
            //M_print(ak,"ak");
            Matrix *dk = M_Transition(ak);
            dk = M_numul(dk,0.5);
            //M_print(dk,"dk");
            Matrix* rk = M_Zeros(1,CodedBitsLength + 1 + KWinLen);
            for(int i = 0; i < CodedBitsLength + 1 + KWinLen; i++){
                double jitter = sigma_jitter;
                double stdpos_d = (i + 1) * T;
                rk->data[i] = readback(stdpos_d,jitter,dk,S,T,TL);
            }
            //-------------Normalization----------------
            MATRIX_TYPE min_rk = M_Min_value(rk->data,rk->column);
            MATRIX_TYPE max_rk = M_Max_value(rk->data,rk->column);
            for(int Norma_i = 0; Norma_i < rk->column; Norma_i++)
            {
                rk->data[Norma_i] = 2.0 * (rk->data[Norma_i] - min_rk)/(max_rk - min_rk) - 1.0;
            }
            Matrix *rk_normarlized = Matrix_copy(rk);
            //------------Equqlization and Detection of Original Signal------------------------
            // PR Equalizer-ML
            //M_print(ak,"ak");
            //M_print(rk_normarlized,"rk_normarlized");
            gen_firtaps_v2(ak,rk_normarlized,
                        gpr_target,fir_length,constraint,method);
            

            M_free(rk_normarlized);
            M_free(rk);
            M_free(dk);
            M_free(ak);
            M_free(codedwords);
            M_free(ChannelBits);
            M_free(tempInputPad);
        }
    }
    M_free(BER1);
    M_free(BER2);
    M_free(BER3);
}
