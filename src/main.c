#include "main.h"

int main() {
    // int gpr_target[5] ={1, 2,2,2, 1};
    // int gpr_length = sizeof(gpr_target)/sizeof(int);
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
            //------------PRBS data---------------------
            srand(1);
            int rand_data[SectorLength + 1];
            for(int i = 0; i < SectorLength; i++)
                rand_data[i] = rand()%2;
            Matrix *ChannelBits = Matrix_gen(1,SectorLength,rand_data);
            //M_print(ChannelBits);
            Matrix *tempInputPad = M_full(ChannelBits,0,0,1,0,0);
            //M_print(tempInputPad);

            //RLL encoder 
            Matrix * codedwords = Matrix_copy(tempInputPad);
            // int codedlen = codedwords->column;
            // int CodedBitsLength = codedlen - 1;

            Matrix *ak = M_full(codedwords,0,0,KWinLen,0,0);
            int mul = 2;
            ak = M_numsub(M_numul(ak,&mul,'i'),1);
            M_print(ak);
            Matrix *dk = Matrix_Transition(ak);
            M_print(dk);
            double mul2 = 0.5;
            dk = M_numul(dk,&mul2,'d');
            M_print(dk);

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
