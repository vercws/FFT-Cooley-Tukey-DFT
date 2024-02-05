#include "FFT.h"

int main()
{
    const int MAX_ORDER = 13;
    const bool PRINT_COEFS = false;

    for (int o = 1; o <= MAX_ORDER; o++)
    {
        const int N = 1 << o;
        std::cout << "N: " << N << std::endl;

        double* f = new double[N];
        for (int n = 0; n < N; n++)
        {
            f[n] = n / (double)N;
        }

        clock_t t1 = clock();
        std::complex<double>* cDFT = dft(f, N);
        clock_t t2 = clock();
        double dft_time = (t2 - t1) / (double)CLOCKS_PER_SEC * 1000.0;
        std::cout << "DFT time[ms]: " << dft_time << std::endl;
        if (PRINT_COEFS == true)
        {
            for (int i = 0; i < N; i++)
            {
                std::cout << std::real(cDFT[i]) << " " << std::imag(cDFT[i]) << std::endl;
            }
        }

        t1 = clock();
        std::complex<double>* cFFT = fft(f, N);
        t2 = clock();
        double fft_time = (t2 - t1) / (double)CLOCKS_PER_SEC * 1000.0;
        std::cout << "FFT time[ms]: " << fft_time << std::endl;
        if (PRINT_COEFS == true)
        {
            for (int i = 0; i < N; i++)
            {
                std::cout << cFFT[i].real() << " " << cFFT[i].imag() << std::endl;
            }
        }

    }
}