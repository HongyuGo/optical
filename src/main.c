#include "main.h"
#include "function.h"
#include "matrix.h"
#include "params.h"

int main() {
    MATRIX_TYPE gpr_target[5] = {1, 2, 2, 2, 1};
    // int gpr_length = sizeof(gpr_target) / sizeof(MATRIX_TYPE);
    int fir_length = 13;
    int SNR_len = sizeof(SNR) / sizeof(int);
    Matrix* BER1 = M_Zeros(SNR_len, 1);
    Matrix* BER2 = M_Zeros(SNR_len, 1);
    Matrix* BER3 = M_Zeros(SNR_len, 1);
    Matrix* gpr_targets = Matrix_gen(1, 5, gpr_target);
    for (int isnr = 1; isnr <= SNR_len; isnr++) {
        int numBits = 0;
        int numErrs = 0;
        if (numErrs < maxErrs && numBits < maxBits) {
            //------------Generate random
            // data---------------------------------------
            srand(1);
            MATRIX_TYPE rand_data[SectorLength];
            for (int i = 0; i < SectorLength; i++) rand_data[i] = rand() % 2;
            Matrix* ChannelBits = Matrix_gen(1, SectorLength, rand_data);

            // MATRIX_TYPE test[10]={0,0,1,0,1,1,0,0,1,0};
            // ChannelBits=Matrix_gen(1, 10, test);
            // M_print(ChannelBits);

            //-------------------------------------------------------------------------

            //------------Read data through file
            // operation-----------------------------
            // FILE *fp = NULL;
            // fp = fopen("output.txt","r");
            // MATRIX_TYPE test_data[Test_len];
            // for(int i = 0; i < Test_len; i++){
            //     fscanf(fp,"%lf",&test_data[i]);
            // }
            // for(int i = 0; i < Test_len; i++){
            //     printf("%lf ",test_data[i]);
            // }
            // printf("\n");
            // fclose(fp);
            //--------------------------------------------------------------------------
            Matrix* tempInputPad = Matrix_copy(ChannelBits);
            // M_print(tempInputPad);
            MATRIX_TYPE src[1][SectorLength];
            MATRIX_TYPE des[1][SectorLength * 3 / 2];
            for (int i = 0; i < SectorLength; i++) {
                src[0][i] = tempInputPad->data[0][i];
            }
            Matrix* codedwords = M_Zeros(1, SectorLength * 3 / 2);
            encode_17pp(SectorLength, SectorLength * 3 / 2, src, des);
            for (int i = 0; i < SectorLength * 3 / 2; i++) {
                codedwords->data[0][i] = des[0][i];
            }
            // RLL encoder
            int codedlen = codedwords->column;
            int CodedBitsLength = codedlen;
            Matrix* ak = M_Zeros(1, codedlen + 1);
            for (int i = 1; i < codedlen + 1; i++) {
                ak->data[0][i] = Xor(codedwords->data[0][i - 1], ak->data[0][i - 1]);
            }
            // M_print(codedwords,"co");

            ak = M_numsub(M_nummul(ak, 2.0), 1.0);

            Matrix* rk = M_Zeros(1, CodedBitsLength + 1 + KWinLen);
            for (int i = 0; i < CodedBitsLength + 1 + KWinLen; i++) {
                // double jitter = sigma_jitter;
                double jitter = 0;
                double stdpos_d = (i + 1) * TL;
                rk->data[0][i] = readback(stdpos_d, jitter, ak, S, T, TL);
            }

            //-------------Normalization----------------
            Matrix* rk_normarlized = mapminmax(rk);
            // M_print(rk_normarlized, "rk_norma");
            // MATRIX_TYPE
            // fw[16]={-0.227109519750194,-0.0842158158970587,0.382645484641129,0.864440958767909,1.03843688740715,0.893878516295095,0.561414192286542,0.262979919285228,0.0477087341605870,-0.0165421190636973,-0.144869870150452,-0.426865942195521,-0.819682351563388,-1.06331073638984,-0.997975102313725,-0.725613028137702};
            // rk_normarlized=Matrix_gen(1, 16, fw);
            //------------Equqlization and Detection of Original
            // Signal------------------------
            // PR Equalizer-ML
            // M_print(ak,"ak");
            // M_print(rk_normarlized, "rk_normarlized");
            // Matrix** return_back = gen_firtaps_v2(ak, rk_normarlized,
            // gpr_target, fir_length, constraint, method); Matrix* fir_taps1 =
            // return_back[0]; M_print(fir_taps1, "fir_taps1"); Matrix*
            // gpr_coeff = return_back[1]; M_print(gpr_coeff, "gpr_coeff");
            // Write_fir_gpr(fir_taps1,gpr_coeff);
            MATRIX_TYPE filter_temp[13] = {-0.0025, -0.0090, -0.0145, 0.0177,  0.1166,  0.2417, 0.3001,
                                           0.2417,  0.1166,  0.0177,  -0.0145, -0.0090, -0.0025};
            Matrix* fir_taps_filter = Matrix_gen(1, 13, filter_temp);
            Matrix* fk_filter = M_Conv(rk_normarlized, fir_taps_filter);
            fk_filter = mapminmax(fk_filter);
            // M_print(fk_filter, "fk_filter");
            Matrix* temp_output = M_Conv(ak, gpr_targets);
            Matrix* ideal_output = M_Cut(temp_output, 1, 1, 2, temp_output->column - 3);
            ideal_output = mapminmax(ideal_output);
            // Matrix* un = M_Cut(fk_filter, 1, 1, 7, ideal_output->column + 6);
            // Matrix* un_T = M_T(un);
            // Matrix* un_mul = M_mul(un, un_T);
            // MATRIX_TYPE fe = M_eigen_struct(un_mul);
            // MATRIX_TYPE mu = 2 * (1 / fe);
            MATRIX_TYPE mu = 0.0023;
            Matrix* fk_filter_cut = M_Cut(fk_filter, 1, 1, 7, ideal_output->column + 6);
            // M_print(ideal_output, "ideal");
            Matrix* fir_taps_lms = LMS(fk_filter_cut, ideal_output, mu * 13, fir_length);
            // M_print(fir_taps_lms, "taps_LMS");
            // MATRIX_TYPE
            // temp[13]={2.23493217334475,2.71260601587680,2.21307172365564,1.05168680655958,-0.106425207969233,-0.716184441966804,-0.671748722283417,-0.262117936800589,0.122480498842974,0.258782110383948,0.161146210156126,-0.0188990342407967,-0.136380600668192};
            // Matrix* fir_taps_lms =Matrix_gen(1,13,temp);
            Matrix* fk_lms = M_Conv(fk_filter_cut, fir_taps_lms);
            fk_lms = mapminmax(fk_lms);
            fk_lms = M_numsub(fk_lms, -1);
            fk_lms = M_nummul(fk_lms, 4);
            // M_print(gpr_targets,"gpr_targets");
            // M_print(fir_taps_filter,"fir");
            // M_print(fk_lms, "fk_lms");
            Matrix* detected = viterbi_mlse(gpr_targets->column, fk_lms, M_T(gpr_targets));
            Matrix* uhat = M_Cut(detected, 1, 1, 1, codedlen);
            ak = M_nummul(M_numsub(ak, -1), 0.5);
            ak = M_Cut(ak, 1, 1, 2, ak->column);
            int cur_err = M_Compare(uhat, ak);
            printf("err:%d\n", cur_err);
            // M_print(fk1,"fk1");
            // Matrix* detected = viterbi_mlse(gpr_length, fk1, gpr_coeff);
            // M_print(detected, "detected");

            // free(return_back);
            // M_free(detected);
            // M_free(fk1);
            // M_free(fir_taps1);
            // M_free(gpr_coeff);
            M_free(rk_normarlized);
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
