#ifndef _PARAMS_H_
#define _PARAMS_H_

#define maxErrs (1000)
#define maxBits (1e6)
#define lamda (0.405)
#define NA (0.85)
#define t0 (0.86 * (lamda) / (NA))
#define T (1)
#define TL (0.05587)
#define S ((t0) / (TL))
#define SectorLength (4096)
// #define delta (0.1)
#define edge_width (((fir_length)-1) / 2)
#define sigma_jitter ((0.0) * (TL))
#define KWinLen 0
#define max_metric 1000000
#define pi (3.1415926)

#define _END_ -1
#define _HEAD_ 1
#define _DETAILED_ 2
#define _ROW_ 1
#define _COLUMN_ 0
#define MATRIX_TYPE double
#define TRANS_TYPE double
#define PRECISION "%.6lf "

#define M_mul_001                                     \
    "@ERROR: Matrix_Dimensions "                      \
    "Wrong!\n\tDetails:(M_mul_001)_mat_left->column " \
    "!= _mat_right->row\n"
#define M_value_one_001 "@ERROR: M_value_one_001:The row or column entered is out of range"
#define M_get_one_001 "@ERROR: M_get_one_001:The row or column entered is out of range"
#define M_Transition_001 "@ERROR: The input matrix cannot be empty"
#define M_Transition_002                                                 \
    "@ERROR: The input matrix column must greater than 1 to generate a " \
    "transition matrix"
#define M_add_sub_003                                                        \
    "@ERROR: Matrix_Dimensions Wrong!\n\tDetails:(M_add_sub_003)_mat_subed " \
    "!= _mat_minus\n"
#define M_Cut_005                                                        \
    "@ERROR: Matrix_Cut Over!\n\tDetails:(M_Cut_005)_Cut_tail over_the " \
    "limited\n"
#define M_Cut_006 "@ERROR: Matrix_Cut Wrong!\n\tDetails:(M_Cut_006)_Head_>_Tail\n"
#define M_Cut_007         \
    "@ERROR: Matrix_Cut " \
    "Wrong!\n\tDetails:(M_Cut_007)_Range_can't_be_negative!'\n"
#define M_swap_004                                                           \
    "@ERROR: Matrix_Swap_Line Over!\n\tDetails:(M_swap_004)_Swap_line over " \
    "the limited\n"
#define M_Dia_Inv_002            \
    "@ERROR: Matrix_Dimensions " \
    "Wrong!\n\tDetails:(M_Dia_Inv_002)_mat_left->column != _mat_right->row\n"
#define M_Dia_Inv_023                                                      \
    "@ERROR: Matrix is not invertible!\n\t Details:(M_Dia_Inv_023)Please " \
    "Check: Inverse element of Dia == 0! \n"
#define M_eigen_max_021                                                      \
    "@ERROR: Matrix_Dimensions "                                             \
    "Wrong!\n\tDetails:(M_eigen_max_021)Mat->column != Mat->row!\n\t\t(For " \
    "eigen, the "                                                            \
    "Matrix must be a square matrix!)\n"
extern int SNR;
extern char constraint;
extern char method;
#endif
