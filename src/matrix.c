#include "matrix.h"
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "params.h"
             
MATRIX_TYPE** GetMemory(int row, int col){
    MATRIX_TYPE ** data = NULL;
    data = (MATRIX_TYPE **)malloc(sizeof(MATRIX_TYPE*) * row);
    for(int i = 0; i < row; i++){
        data[i] = (MATRIX_TYPE*)malloc(sizeof(MATRIX_TYPE) * col);
    }
    return data;
}
/* Generate Matrix Struct */
Matrix* Matrix_gen(int row, int column, MATRIX_TYPE* data) {
    Matrix* _mat = (Matrix*)malloc(sizeof(Matrix));
    if (_mat == NULL) return 0;
    int i,j;
    _mat->row = row;
    _mat->column = column;
    _mat->data = GetMemory(_mat->row,_mat->column); 
    for (i = 0; i < row; i++) {
        for(j = 0; j < column; j++){
            _mat->data[i][j] = data[i * row + j];
        }
    }
    return _mat;
}
/* Copy Mtrix(gen new one)*/
Matrix* Matrix_copy(Matrix* _mat_sourse) {
    int row = _mat_sourse->row;
    int col = _mat_sourse->column;
    int size = row * col;
    MATRIX_TYPE *data = (MATRIX_TYPE*)malloc(sizeof(MATRIX_TYPE) * size);
    for(int i = 0; i < row; i++){
        for(int j = 0; j < col; j++){
            data[i * row + j] = _mat_sourse->data[i][j];
        }
    }
    Matrix* _mat_copy = Matrix_gen(_mat_sourse->row, _mat_sourse->column, data);
    free(data);
    return _mat_copy;
}
/* Free Memory*/
int M_free(Matrix* _mat) {
    for(int i = 0; i < _mat->row; i++){
        free(_mat->data[i]);
    }
    free(_mat->data);
    //(_DETAILED_ >= 3) ? printf(">>Matrix_%x has been freed.\n", _mat) : 0;
    free(_mat);
    return 0;
}

/* Add & Sub*/
Matrix* M_add_sub(MATRIX_TYPE scale_mat_subed, Matrix* _mat_subed, MATRIX_TYPE scale_mat_minus, Matrix* _mat_minus) {
    Matrix* _mat_result = NULL;
    if ((_mat_subed->column == _mat_minus->column) && (_mat_subed->row == _mat_minus->row)) {
        _mat_result = Matrix_copy(_mat_subed);
        for(int i = 0; i < _mat_subed->row; i++){
            for(int j = 0; j < _mat_subed->column; j++){
                _mat_result->data[i][j] = (_mat_result->data[i][j]) * scale_mat_subed - (_mat_minus->data[i][j]) * scale_mat_minus;
            }
        }
    } 
    else {
        printf(M_add_sub_003);
    }
    return _mat_result;
}
/*Cut_out_part_of_Matrix*/
Matrix* M_Cut(Matrix* _mat, int row_head, int row_tail, int column_head, int column_tail) {
    Matrix* mat_result = NULL;
    if (row_tail < 0) {
        if (row_tail == _END_) {
            row_tail = _mat->row;
        } else {
            printf(M_Cut_007);
            system("pause");
        }
    }

    if (row_head < 0) {
        if (row_head == _END_) {
            row_head = _mat->row;
        } else {
            printf(M_Cut_007);
            system("pause");
        }
    }

    if (column_tail < 0) {
        if (column_tail == _END_) {
            column_tail = _mat->column;
        } else {
            printf(M_Cut_007);
            system("pause");
        }
    }

    if (column_head < 0) {
        if (column_head == _END_) {
            column_head = _mat->column;
        } else {
            printf(M_Cut_007);
            system("pause");
        }
    }

    if ((row_tail > _mat->row) || (column_tail > _mat->column)) {
        printf(M_Cut_005);
        system("pause");
    } else {
        if ((row_head > row_tail) || (column_head > column_tail)) {
            printf(M_Cut_006);
            system("pause");
        } else {
            row_head = row_head - 1;
            column_head = column_head - 1;
            mat_result = (Matrix*)malloc(sizeof(Matrix));
            mat_result->row = row_tail - row_head;
            mat_result->column = column_tail - column_head;
            mat_result->data = GetMemory(mat_result->row,mat_result->column);
            int i, j;
            for (i = 0; i < (row_tail - row_head); i++) {
                for (j = 0; j < (column_tail - column_head); j++) {
                    // mat_result->data[i * (mat_result->column) + j] =
                    //     _mat->data[(i + row_head) * (_mat->column) + (j + column_head)];
                    mat_result->data[i][j] = _mat->data[i+row_head][j+column_head];
                }
            }
        }
    }
    return mat_result;
}
/*Full*/
Matrix* M_full(Matrix* _mat, int row_up, int row_down, int column_left, int column_right, MATRIX_TYPE full_data) {
    Matrix* mat_result = NULL;
    mat_result = (Matrix*)malloc(sizeof(Matrix));
    mat_result->row = (_mat->row + row_up + row_down);
    mat_result->column = (_mat->column + column_left + column_right);
    mat_result->data = GetMemory(mat_result->row,mat_result->column);
    int i, j;
    for (i = 0; i < mat_result->row; i++) {
        for (j = 0; j < mat_result->column; j++) {
            if ((i >= row_up) && (i < (row_up + _mat->row))) { /*这里的双判断，可以优化*/
                if ((j >= column_left) && (j < (column_left + _mat->column))) {
                    // mat_result->data[i * (mat_result->column) + j] =
                    //     _mat->data[(_mat->column) * (i - row_up) + (j - column_left)];
                    mat_result->data[i][j] = _mat->data[i -row_up][j - column_left];
                } else {
                    // mat_result->data[i * (mat_result->column) + j] = full_data;
                    mat_result->data[i][j] = full_data;
                }
            } else {
                // mat_result->data[i * (mat_result->column) + j] = full_data;
                mat_result->data[i][j] = full_data;
            }
        }
    }
    return mat_result;
}
/*Generate Zeros _matrix*/
Matrix* M_Zeros(int row, int column) {
    Matrix* Zero_mat = (Matrix*)malloc(sizeof(Matrix));
    int i,j;
    Zero_mat->column = column;
    Zero_mat->row = row;
    MATRIX_TYPE** data = GetMemory(Zero_mat->row,Zero_mat->column);
    for (i = 0; i < Zero_mat->row; i++) {
        for(j = 0; j < Zero_mat->column; j++){
            data[i][j] = 0;
        }
    }
    Zero_mat->data = data;
    return Zero_mat;
}
/*Print Matrix*/
int M_print(Matrix* _mat, const char *name) {
    int i, j;
    printf("%s row:%d col:%d\n",name, _mat->row, _mat->column);
    for (i = 0; i < _mat->row; i++) {
        for (j = 0; j < _mat->column; j++) {
            printf(PRECISION, _mat->data[i][j]);
        }
        printf("\n");
    }
    return 0;
}
/*Matrix Multiply*/
Matrix* M_nummul(Matrix* _mat, double _num) {
    MATRIX_TYPE** data = _mat->data;
    int i, j;
    for (i = 0; i < _mat->row; i++) {
        for(j = 0; j < _mat->column; j++){
            data[i][j] = data[i][j] * _num;
        }
    }
    return _mat;
}
/*Matrix sub*/
Matrix* M_numsub(Matrix* _mat, MATRIX_TYPE _num) {
    MATRIX_TYPE** data = _mat->data;
    int i, j; 
    for (i = 0; i < _mat->row; i++) {
        for(j = 0; j < _mat->column; j++){
            data[i][j] = data[i][j] - _num;
        }
    }
    return _mat;
}
/*Generation of transition matrix*/
Matrix* M_Transition(Matrix* _mat) {
    if (_mat == NULL) {
        printf(M_Transition_001);
        return NULL;
    }
    if (_mat->column < 2) {
        printf(M_Transition_002);
        return NULL;
    }
    Matrix* _mat_result = NULL;
    int row = _mat->row,col = _mat->column,i,j;
    MATRIX_TYPE** _data = GetMemory(row,col);
    for (i = 0; i < row; i++) {
        for(j = 0; j < col - 1; j++){
            _data[i][j] = _mat->data[i][j + 1] - _mat->data[i][j];
        }
    }
    _mat_result = (Matrix*)malloc(sizeof(Matrix));
    _mat_result->row = row;
    _mat_result->column = col - 1;
    _mat_result->data = _data;
    return _mat_result;
}
/*Find min value in a MATRIX_TYPE[*]*/
MATRIX_TYPE M_Min_value(MATRIX_TYPE *data, int size) {
    MATRIX_TYPE Val_min = data[size - 1];
    for (int i = size - 2; i >= 0; i--) {
        if (data[i] <= Val_min) {
            Val_min = data[i];
        }
    }
    return Val_min;
}
/*Find max value in a MATRIX_TYPE[*]*/
MATRIX_TYPE M_Max_value(MATRIX_TYPE *data, int size) {
    MATRIX_TYPE Val_max = data[size - 1];
    for (int i = size - 2; i >= 0; i--) {
        if (data[i] >= Val_max) {
            Val_max = data[i];
        }
    }
    return Val_max;
}
/*Assign a value to a position of the matrix*/
void M_value_one(Matrix *_mat, int row , int col,MATRIX_TYPE value){
    int _m_row = _mat->row;
    int _m_column = _mat->column;
    if(row <= 0 || row > _m_row || col <= 0 || col > _m_column){
        printf("%s\n",M_value_one_001);
        return;
    }
    _mat->data[row - 1][col - 1] = value;
}
/*Get a value to a position of the matrix*/
MATRIX_TYPE M_get_one(Matrix *_mat, int row, int col){
    int _m_row = _mat->row;
    int _m_column = _mat->column;
    if(row <= 0 || row > _m_row || col <= 0 || col > _m_column){
        printf("%s\n",M_get_one_001);
        return 0.0;
    }
    return _mat->data[row - 1][col - 1];
}
/*Transpose*/
Matrix *M_T(Matrix *_mat_source) {
    Matrix *_mat = (Matrix *) malloc(sizeof(Matrix));
    _mat->column = _mat_source->row;
    _mat->row = _mat_source->column;
    MATRIX_TYPE **data = GetMemory(_mat->row,_mat->column);
    _mat->data = data;
    int i, j;
    for (i = 0; i < (_mat->row); i++) {
        for (j = 0; j < _mat->column; j++) {
            data[i][j] = _mat_source->data[j][i];
        }
    }
    return _mat;
}
/*Matrix Multiply*/
Matrix *M_mul(Matrix *_mat_left, Matrix *_mat_right) {
	/*_mat_result = _mat_left*_mat_right */
    //(_DETAILED_>=3)?printf(">>Matrix_%x * Matrix_%x =\n", _mat_left, _mat_right):0;
    /*Determine_Matrix_Dimensions*/
    Matrix *_mat_result = NULL;
    if (_mat_left->column != _mat_right->row) {
        printf(M_mul_001);
    } else {
        _mat_result = (Matrix *) malloc(sizeof(Matrix));
        int row = _mat_left->row;
        int mid = _mat_left->column;
        int column = _mat_right->column;
        int i, j, k;
        MATRIX_TYPE **_data = GetMemory(row,column);
        MATRIX_TYPE temp = 0;
        /*Ergodic*/
        for (i = 0; i < row; i++) {
            for (j = 0; j < column; j++) {
                /*Caculate Element*/
                temp = 0;
                for (k = 0; k < mid; k++) {
                    temp += (_mat_left->data[i][k]) * (_mat_right->data[k][j]);
                }
                _data[i][j] = temp;
            }
        }
        _mat_result->row = row;
        _mat_result->column = column;
        _mat_result->data = _data;
    }
    //(_DETAILED_>=3)?printf("\tMatrix_%x\n", _mat_result):0;
    return _mat_result;
}
/*Matrix inverse*/
Matrix *M_Inverse(Matrix *_mat){
    MATRIX_TYPE **TempA = _mat->data;
    Matrix *_mat_result = (Matrix*)malloc(sizeof(Matrix));
    MATRIX_TYPE **B = GetMemory(_mat->row,_mat->column);
    for(int i = 0; i < _mat->row; i++){
        for(int j = 0; j < _mat->column; j++){
            B[i][j] = 0.0;
        }
    }
	int c,i,j,m = _mat->row;
	double temp;
	double **A;

	/* Copying the input matrix TempA into matrix A */
	A = (double**)malloc(m*sizeof(double));
	for(i=0;i<m;i++)
	{
		A[i] = (double*)malloc(m*sizeof(double));
	}

	for(i=0;i<m;i++)
	{
		for(j=0;j<m;j++)
		{
			A[i][j] = TempA[i][j];
		}
	}

	/* Creating the identity matrix B */
	/* Memory for Matrix B should already be allocated and initialized to 0 */
	for(i=0;i<m;i++)
	{
			B[i][i] = 1;
	} 

	/*** Performing Gaussian elimination on A ***/

	for(c=0;c<m;c++)
	{/* loop over all columns of A */

		/*
		If diagnol element of column c is 0, then interchange row c with any other row
		so that diagnol element of column c is non-zero.
		*/
		if(A[c][c] == 0)
		{
			for(i=0;i<m;i++)
			{
				if(A[i][c] != 0)
				{
					/** Interchange row i and row c **/
					for(j=0;j<m;j++)
					{
						temp = A[i][j];
						A[i][j] = A[c][j];
						A[c][j] = temp;
						/* Do the same operation for B */
						temp = B[i][j];
						B[i][j] = B[c][j];
						B[c][j] = temp;
					}
					break;
				}
			}
		}

		/*
		Making the diagnol element of column c, A[c][c] 1 by dividing
		all the elements of row c, i.e. A[c][0 1 ... m-1] by A[c][c].
		*/
		if(A[c][c] != 1)
		{
			temp = A[c][c];
			for(j=0;j<m;j++)
			{
				A[c][j] = A[c][j]/temp;
				/* Do the same operation on the matrix B */
				B[c][j] = B[c][j]/temp;
			}
		}

		for(i=0;i<m;i++)
		{/* loop over all rows of A */

			/*
			Making all elements of column c as zero, except for the diagnol element,
			which is 1 at this point.
			*/
			if(i!=c)
			{
				temp = A[i][c];
				for(j=0;j<m;j++)
				{
					A[i][j] = A[i][j] - temp*A[c][j];
					B[i][j] = B[i][j] - temp*B[c][j];
				}
			}
		}/* loop over all rows of A */

	}/* loop over all columns of A */
	
	for(i=0;i<m;i++)
	{
		free(A[i]);
	}
	free(A);
    _mat_result->row = _mat->row;
    _mat_result->column = _mat->column;
    _mat_result->data = B;
    return _mat_result;
}

/*Matrix Convolution*//*this is for the matrix only like 1*n or n*1*/
Matrix *M_Conv(Matrix *_mat1, Matrix *_mat2){
    int flag1,flag2,len1,len2,min,max;
    if(_mat1->row > 1){
        flag1 = 0;
        len1 = _mat1->row;
    } 
    else{
        flag1 = 1;
        len1 = _mat1->column;
    } 
    if(_mat2->row > 1){
        flag2 = 0;
        len2 = _mat2->row;
    } 
    else{
        flag2 = 1;
        len2 = _mat2->column;
    }
    int len = len1+len2-1;
    Matrix *_mat_result = NULL;
    _mat_result = (Matrix*)malloc(sizeof(Matrix));
    _mat_result->row = 1;
    _mat_result->column = len;
     _mat_result->data = GetMemory(_mat_result->row,_mat_result->column); 
    for(int k=0;k<len;++k){
        max = 0 > (k+1-len2)?0:(k+1-len2);
        min = k<(len1-1)?k:(len1-1);
        if(flag1&&flag2){
            for(int i = max;i<=min;i++){
                _mat_result->data[0][k] += _mat1->data[0][i] * _mat2->data[0][k-i]; 
            }
        }
        else if(flag1 == 1 && flag2 == 0){
            for(int i = max;i<=min;i++){
                _mat_result->data[0][k] += _mat1->data[0][i] * _mat2->data[k-i][0]; 
            }
        }
        else if(flag1 == 0 && flag2 == 1){
            for(int i = max;i<=min;i++){
                _mat_result->data[0][k] += _mat1->data[i][0] * _mat2->data[0][k-i]; 
            }
        }
        else{
            for(int i = max;i<=min;i++){
                _mat_result->data[0][k] += _mat1->data[i][0] * _mat2->data[k-i][0]; 
            }
        }
        
    }
    return _mat_result;
}
/**/
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