#include "function.h"
#include <math.h>
#include "matrix.h"
#include "commom.h"
#include "params.h"

double h_response(double t, double TLL, double SS) { return erf(t / SS / TLL); }

double readback(double t, double jitter, Matrix* d, double SS, double TT, double TLL) {
    int len = d->column;
    double rs = 0;
    double tmp = 0;
    for (int i = 0; i < len; i++) {
        tmp = d->data[0][i] * h_response(t - (i + 1) * TT + jitter, TLL, SS);
        rs = rs + tmp;
    }
    return rs;
}

Matrix** gen_firtaps_v2(Matrix *random, Matrix *sampled, MATRIX_TYPE* gpr_coeff_data, int fir_len, char constraint,char method){
    int gpr_len = 5;
    // int i,j;
    //int datalength = random->column;
    Matrix *I_matrix = M_Zeros(gpr_len,1);
    switch(constraint){
        case '1':
                break;
        case 'c':
                M_value_one(I_matrix,gpr_len/2 + 1,1,1.0);
                //M_print(I_matrix,"I_matrix");
                break;
        case '2':
                break;
        default:
                break;
    }
    //-----------R_matrix.. Calculate the autocorrelation matrix of the sampled data---------------
    // Size:(2K+1)x(2K+1)
    int K = (fir_len - 1) / 2;
    //printf("%d",K);
    Matrix * R_matrix = auto_corr(sampled,-K,K);
    // R_matrix = M_Inverse(R_matrix);
    // M_print(R_matrix,"R_matrix");
    // for(i = 1; i < K; i++){
    //     for(j = 0; j < K; j++){
    //         M_value_one(R_matrix,i + 1,j + 1,M_get_one(R_matrix,1,abs(j - i) + 1));
    //     }
    // }

    //-----------A_matrix.. Caculate the autocorrelation matrix of the random data-----------------
    //Size: LxL
    int L = gpr_len;
    Matrix * A_matrix = auto_corr(random,0,L - 1);
    //A_matrix = M_Inverse(A_matrix);
    //M_print(A_matrix,"A_matrix");
    
    //-----------T_matrix.. Caculate the cross-correlation matrix of the sampled and the random-----------
    //data. 
    //Size: (2K+1)xL
    Matrix *T_matrix = cross_corr(sampled, random, fir_len, gpr_len);
    //M_print(T_matrix,"T_matrix");


    MATRIX_TYPE lagrange = Caculate_lagrange(R_matrix,A_matrix,T_matrix,I_matrix);
    //printf("lagrange = %lf\n",lagrange);
    //M_print(M_lagrange,"data");

    Matrix *gpr_coeff = Caculate_gpr_coeff(R_matrix,A_matrix,T_matrix,I_matrix,lagrange);
    //M_print(gpr_coeff,"gpr_coeff");

    Matrix *fir_coeff = Caculate_fir_coeff(R_matrix,T_matrix,gpr_coeff);
    //M_print(fir_coeff,"fir_coeff");

    Matrix **return_back = NULL;
    return_back = (Matrix**)malloc(sizeof(Matrix*) * 2);
    return_back[1] = gpr_coeff;
    return_back[0] = fir_coeff;

    M_free(T_matrix);
    M_free(A_matrix);
    M_free(I_matrix);
    M_free(R_matrix);
    return return_back;
    
}
Matrix * auto_corr(Matrix * x, int bot , int top){
    int L = top - bot + 1;
    int i, j,k,end = x->column;
    MATRIX_TYPE sum = 0.0,prod = 0.0;
    Matrix * ac = M_Zeros(L,L);
    for(i = bot ; i <= top ; i++){
        for(j = bot; j <= top; j++){
           if(end + bot < top + 1)
                break;
           else{
            for(k = top + 1; k <= end + bot; k++){
                sum += x->data[0][k - i - 1] * x->data[0][k - j - 1];
            }
            prod = sum / (end + bot - (top + 1) + 1);
            sum = 0.0;
            M_value_one(ac,i-bot+1,j-bot+1,prod);
           } 
        }
    }
    return ac;
}
Matrix *cross_corr(Matrix * x, Matrix * y, int K, int L){
    int hK = (K - 1) / 2;
    int i, j, k;
    MATRIX_TYPE prod = 0.0, sum = 0.0;
    int end = x->column;
    Matrix *xc = M_Zeros(K,L);
    int tmp = hK > L-1 ? hK: L - 1;
    for(i = -hK; i <= hK; i++){
        for(j = 0 ; j <= L - 1 ; j++){
            if(1 + tmp > end - tmp)
                break;
            else{
                for(k = 1 + tmp; k <= end - tmp; k++){
                    sum += x->data[0][k - i - 1] * y->data[0][k - j - 1];
                }
                prod = sum / (end - tmp - (1 + tmp) + 1);
                sum = 0.0;
                M_value_one(xc,i + hK + 1,j + 1,prod);
            }
        }
    }
    return xc;
}
MATRIX_TYPE Caculate_lagrange(Matrix *R_matrix, Matrix *A_matrix, Matrix *T_matrix, Matrix *I_matrix){
    //M_print(R_matrix,"R_matrix");
    Matrix *R_inv = M_Inverse(R_matrix);
    //M_print(R_inv,"R_inv");
    
    Matrix *T_T = M_T(T_matrix);
    Matrix *TR = M_mul(T_T,R_inv);
    M_free(T_T); 
    M_free(R_inv);

    Matrix *TRT = M_mul(TR,T_matrix);
    Matrix *A_TRT = M_add_sub(1.0,A_matrix,1.0,TRT);
    M_free(TRT);
    M_free(TR);

    Matrix *A_TRT_inv = M_Inverse(A_TRT);
    Matrix *I_T = M_T(I_matrix);
    Matrix *I_A_TRT = M_mul(I_T,A_TRT_inv);
    M_free(A_TRT_inv);
    M_free(A_TRT);
    M_free(I_T);

    Matrix *M_lagrange = M_mul(I_A_TRT,I_matrix);
    M_free(I_A_TRT);
    
    MATRIX_TYPE lagrange = 1.0/M_lagrange->data[0][0];
    M_free(M_lagrange);
    return lagrange;
}
Matrix *Caculate_gpr_coeff(Matrix *R_matrix, Matrix *A_matrix, Matrix *T_matrix, Matrix *I_matrix, MATRIX_TYPE lagrange){
    Matrix *R_inv = M_Inverse(R_matrix);
    //M_print(R_inv,"R_inv");
    Matrix *T_T = M_T(T_matrix);
    Matrix *TR = M_mul(T_T,R_inv);
    M_free(T_T); 
    M_free(R_inv);

    Matrix *TRT = M_mul(TR,T_matrix);
    Matrix *A_TRT = M_add_sub(1.0,A_matrix,1.0,TRT);
    //M_print(A_TRT,"A_TRT");
    M_free(TRT);
    M_free(TR);

    Matrix *A_TRT_inv = M_Inverse(A_TRT);
    //M_print(A_TRT_inv,"A_TRT_inv");
    Matrix *gpr_coeff = M_mul(A_TRT_inv,I_matrix);
    M_free(A_TRT);
    M_free(A_TRT_inv);

    gpr_coeff = M_nummul(gpr_coeff,lagrange);
    return gpr_coeff;
}
Matrix *Caculate_fir_coeff(Matrix *R_matrix, Matrix *T_matrix, Matrix *gpr_coeff){
    Matrix *R_inv = M_Inverse(R_matrix);
    Matrix *TR = M_mul(R_inv,T_matrix);
    M_free(R_inv);
    Matrix *TRgpr = M_mul(TR,gpr_coeff);
    M_free(TR);
    return TRgpr;
}
Matrix *viterbi_mlse(int gpr_len,Matrix *fk1, Matrix *gpr_coeff){
    /*
    mat_detected_output, the detected signal sequence through viterbi equalizer
    the length of mat_detected_output is equal to fk1
    */
    Matrix *mat_detected_output = NULL;
    mat_detected_output = (Matrix*)malloc(sizeof(Matrix));
    mat_detected_output->row = 1;
    if(fk1->row==1) mat_detected_output->column = fk1->column;
    else mat_detected_output->column = fk1->row;
    mat_detected_output->data = GetMemory(mat_detected_output->row,mat_detected_output->column); 
    int stateSize = gpr_len -1;
    int numOfStates = pow(stateSize, 2);
    Trellis_nst **trellis_nst = NULL;
    Trellis_pst **trellis_pst = NULL;
    trellis_nst = (Trellis_nst **)malloc(sizeof(Trellis_nst*) * numOfStates * 2);
    trellis_pst = (Trellis_pst **)malloc(sizeof(Trellis_pst*) * numOfStates * 2);
    /*initial trellis_nst*/
    for (int i = 0; i < numOfStates*2; ++i) {
        trellis_nst[i] = (Trellis_nst*)malloc(sizeof(Trellis_nst));   
    }
    /*initial trellis_pst*/
    for(int i=0;i< numOfStates;++i){
        trellis_pst[i] = (Trellis_pst*)malloc(sizeof(Trellis_pst));
        trellis_pst[i]->counter = 0;
    }
    //printf("test:%d",(2*(0&1)-1));
    for (int s = 0; s < numOfStates; ++s) {
        /*b is the branch value*/
        for(int b=1;b<=2;++b){
            int cur = s + numOfStates * (b-1); /*can understand cur = (s,b) or (nextstate,b)*/
            trellis_nst[cur]->input = 2*b - 3;
            trellis_nst[cur]->output = (2*b - 3) * gpr_coeff->data[0][0];
            for(int i=1,tmp = s;i<=stateSize; ++i){
                /*int tmp;put the initialization here is incorrect*/
                trellis_nst[cur]->output += gpr_coeff->data[i][0] * (2*(tmp&1)-1);
                //printf("tmp:%d,gpr_coeff->data[i][0]:%lf,(2*(tmp&1)-1):%d  ",tmp,gpr_coeff->data[i][0],(2*(tmp&1)-1));
                tmp = tmp >> 1;
                //printf("tmp:%d  ",tmp);
            }
            trellis_nst[cur]->next = ((s<<1|1)&((1<<stateSize)+b-3));/*update s*/
            printf("nst - curstate:%d, input:%d, output:%lf, nextstate:%d\n",s, trellis_nst[cur]->input, trellis_nst[cur]->output,trellis_nst[cur]->next);
            int nextstate = trellis_nst[cur]->next ;
            trellis_pst[nextstate]->counter +=1;
            trellis_pst[nextstate]->input[trellis_pst[nextstate]->counter-1] = trellis_nst[cur]->input;
            trellis_pst[nextstate]->output[trellis_pst[nextstate]->counter-1] = trellis_nst[cur]->output;
            trellis_pst[nextstate]->pre[trellis_pst[nextstate]->counter-1] = s;
            printf("trellis.pst(%d,%d).prestate= %d\n",nextstate,trellis_pst[nextstate]->counter,s);
        }
    }

    /*the part of viterbi_mlse*/
    int depth = mat_detected_output->column + 1;
    Matrix *mat_state_metric = NULL;
    mat_state_metric = (Matrix*)malloc(sizeof(Matrix));
    mat_state_metric->row = numOfStates;
    mat_state_metric->column = depth;
    mat_state_metric->data = GetMemory(mat_state_metric->row, mat_state_metric->column);

    int  **survivor_path = NULL;
    survivor_path = (int **)malloc(numOfStates *sizeof(int *));
    for(int i =0;i<numOfStates;++i){
        survivor_path[i] = (int *)malloc((depth-1) * sizeof(int));
    }

    for(int i=0;i<mat_state_metric->row;++i){
        for(int j=0;j<mat_state_metric->column;++j){
            mat_state_metric->data[i][j] = 1000;
        }
    }
    mat_state_metric->data[0][0] = 0; /*assume that the initial state is 0*/
    double *state_value = (double *)malloc(2*sizeof(double));
    double branch_metric;
    for(int i = 1;i < depth;++i){
        for(int j=0;j < numOfStates;++j){
            
            for(int b=1;b<=2;++b){
                double out = trellis_pst[j]->output[b-1];
                double f = fk1->data[0][i-1];
                branch_metric = (f-out)*(f-out);
                int c = trellis_pst[j]->pre[b-1];
                c = c*1;
                state_value[b-1] = mat_state_metric->data[trellis_pst[j]->pre[b-1]][i-1] + branch_metric;
            }
            //printf("v0:%lf,v1:%lf->",state_value[0],state_value[1]);
            if(state_value[0]>state_value[1]){
                mat_state_metric->data[j][i] = state_value[1];
                survivor_path[j][i-1] = 2;
            } 
            else if (state_value[0]<state_value[1]){
                mat_state_metric->data[j][i] = state_value[0];
                survivor_path[j][i-1] = 1;
            } 
            else{
                mat_state_metric->data[j][i] = state_value[0];
                int rand_num = rand();
                if(rand_num%2 == 1) survivor_path[j][i-1] = 1;//for branch 1
                else survivor_path[j][i-1] = 2;//for branch 2 
            }
            //printf("j%d,survivor_path%d ",j,survivor_path[j][i-1]);
            printf("%d  ",survivor_path[j][i-1]);
            //printf("%.4lf  ",mat_state_metric->data[j][i-1]);
        }
        printf("\n");
    }
    int state = 0;
    for(int i = depth-1;i>=1;--i){
        mat_detected_output->data[0][i-1] = trellis_pst[state]->input[survivor_path[state][i-1] - 1];
        printf("state:%d, survivor_path:%d\n",state, survivor_path[state][i-1]);
        //printf("%.4lf",mat_detected_output->data[0][i-1]);
        state = trellis_pst[state]->pre[-1+survivor_path[state][i-1]];
    }
    /*map output (-1,1)->(0,1)*/
    mat_detected_output = M_nummul(M_numsub(mat_detected_output,-1.0),0.5);/*(m+1)*2*/
    /*free the space of avoid Memory leak*/
    for (int i = 0; i < numOfStates * 2; ++i) {
        free(trellis_nst[i]);  
    }
    for (int i = 0; i < numOfStates ; ++i) {
        free(trellis_pst[i]);   
    }
    for(int i=0;i<numOfStates;++i){
        free(survivor_path[i]);
    }
    free(trellis_nst);
    free(trellis_pst);
    free(state_value);
    free(survivor_path);
    return mat_detected_output;
}
