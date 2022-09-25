#ifndef _PARAMS_H_
#define _PARAMS_H_

#define maxErrs 100
#define maxBits 1e10
#define lamda 0.405
#define NA 0.65
#define t0 0.86 * lamda / NA
#define T 1
#define TL 0.102
#define S t0 / TL
#define SectorLength 4096
#define delta 0.1
#define fir_length 13
#define edge_width (fir_length - 1) / 2
#define sigma_jitter 0.01 * TL
#define KWinLen 1

#endif
