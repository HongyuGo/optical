#include "main.h"

int main() {
	int* ptr = (int*)malloc(sizeof(int));
	*ptr = 7;
    free(ptr);
	// int maxErrs = 100;
	// int maxBits = 1e6;
    // int lamda = 0.405;
    // double NA = 0.65;
    // double t0 = 0.86 * lamda / NA;
    // int T = 1;
    // double TL = 0.102;
    // double S = t0 / TL;

#if 0

    int gpr_target = (1, 2, 2, 2, 1);
    int gpr_length = len(gpr_target);
    int fir_length = 13; // Must be an odd number.
    int edge_width = (fir_length - 1) / 2;


    string constraint = "centre";
    char method = '1';
    int rate = 3 / 5;
    sigma_jitter = 0.01 * TL;
    % jitter noise
            approx = 'f';
    % channel approximation mode : 'f', '1', '2'.

                                             SNR = 29;
    BER = zeros(length(SNR), 1);

    numBits = 0;
    numErrs = 0;
    while (numErrs1 < maxErrs && numBits < maxBits)
    {
        /*----------PRBS data--------------*/
        ChannelBits = randi([0 1], 1, SectorLength);

        ak = [ zeros(1, KWinLen), SectorLength ];
        ak = ak * 2 - 1; // BPSK
        dk = 1 / 2 * (ak(2
                         : end) -
                      ak(1
                         : end - 1)); // (1-D) precoder
        /*----------Readback signal--------*/
        rk = zeros(1, SectorLength + KWinLen);
        for (i = 1; i = i + 1; i <= SectorLength + 1 + KWinLen)
        {
            jitter = normrnd(0, sigma_jitter);
            stdpos_d = i * T;
            rk(i) = readback(stdpos_d, jitter, dk, S, T, TL);
        }
    }
#endif

}
