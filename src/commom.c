#include "commom.h"
void Write_fir_gpr(Matrix *fir_taps1, Matrix *gpr_coeff){
    FILE *fp;
    int i,j;
    const char path[] = {"./output/output_fir_gpr.txt"};
    const char fir_taps1_str[] = {"fir_taps1\n"};
    const char gpr_coeff_str[] = {"gpr_coeff\n"};
    fp = fopen(path,"w");
    fputs(fir_taps1_str,fp);
    for(i = 0; i < fir_taps1->row; i++){
        for(j = 0; j < fir_taps1->column; j++){
            fprintf(fp,"%lf ",fir_taps1->data[i][j]);
        }
        fputc('\n',fp);
    }
    fputs(gpr_coeff_str,fp);
    for(i = 0; i < gpr_coeff->row; i++){
        for(j = 0; j < gpr_coeff->column; j++){
            fprintf(fp,"%lf ",gpr_coeff->data[i][j]);
        }
        fputc('\n',fp);
    }
    fclose(fp);
}