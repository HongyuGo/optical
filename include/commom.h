#ifndef _COMMOM_H_
#define _COMMOM_H_

#define _END_ -1
#define _HEAD_ 1
#define _DETAILED_ 2
#define MATRIX_TYPE double
#define PRECISION "%d\t"

#define Matrix_Transition_001 "@ERROR: 输入矩阵不能为空"
#define Matrix_Transition_002 "@ERROR: 输入矩阵column必须大于1，才能生成跃迁矩阵"
#define M_add_sub_003 "@ERROR: Matrix_Dimensions Wrong!\n\tDetails:(M_add_sub_003)_mat_subed != _mat_minus\n"
#define M_Cut_005 "@ERROR: Matrix_Cut Over!\n\tDetails:(M_Cut_005)_Cut_tail over_the limited\n"
#define M_Cut_006 "@ERROR: Matrix_Cut Wrong!\n\tDetails:(M_Cut_006)_Head_>_Tail\n"
#define M_Cut_007 "@ERROR: Matrix_Cut Wrong!\n\tDetails:(M_Cut_007)_Range_can't_be_negative!'\n"
#endif
