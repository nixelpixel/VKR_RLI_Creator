/*
 * CDSPCPU.cpp
 *
 *  Created on: Nov 10, 2019
 *      Author: user
 */

#ifdef __x86_64__
#define X86MATH
#endif

#include "CDSPFFTW.hpp"
#include <omp.h>

#ifdef X86MATH
#include "x86math.hpp"
#endif

CDSP_FFTW::CDSP_FFTW(CLog *log) /*: CDSP(log)*/ {
    Log = log;
    FFTPlan = NULL;

#ifdef X86MATH
    x86CreateSinCosTable();
#endif
}

CDSP_FFTW::CDSP_FFTW() {
    // TODO Auto-generated constructor stub
    FFTPlan = NULL;
}

CDSP_FFTW::~CDSP_FFTW() {
    // TODO Auto-generated destructor stub
    if (FFTPlan) {
        fftwf_destroy_plan(FFTPlan);
    }

#ifdef X86MATH
    x86FreeSinCosTable();
#endif
}

/*!
   \brief Вывод информации о версии библиотеки, частоте процессора, количестве
   ядер и возможности процессора

    \return      Функция не возвращает значения
*/

//============================
// fft()
//============================
int CDSP_FFTW::fft(Ipp32f *pImageData, int N, int M) {
    int order = log2(N);
    assert(N == pow(2, order)); // N должно быть степенью двойки

    if (!FFTPlan) {
        FFTPlan = fftwf_plan_dft_r2c_1d(N, NULL, NULL, FFTW_ESTIMATE);
    }

#pragma omp parallel for
    for (int j = 0; j < M; j++) {
        float *p = &pImageData[j * N];
        p[0] = 0.0f;
        p[1] = 0.0f;
        fftwf_execute_dft_r2c(FFTPlan, p, (fftwf_complex *)p);
        p[0] = 0;
        p[1] = 0; // Убираем постоянную составляющую
    }

    return -1;
}

int CDSP_FFTW::fft(Ipp32fc *pImageData, int N, int M) {
    int order = log2(N);
    assert(N == pow(2, order)); // N должно быть степенью двойки

    int nbthreads = omp_get_max_threads();
    fftwf_init_threads();
    fftwf_plan_with_nthreads(nbthreads);
    omp_set_num_threads(nbthreads);

    fftwf_plan p;
    int rank = 1;
    int n[] = {(int)N};

#pragma omp critical
    p = fftwf_plan_many_dft(rank, n, M, (fftwf_complex *)pImageData, NULL, 1, N,
                            (fftwf_complex *)pImageData, NULL, 1, N,
                            FFTW_FORWARD, FFTW_ESTIMATE);
    fftwf_execute(p);

#pragma omp critical
    fftwf_destroy_plan(p);

    /*  //#pragma omp parallel for
        for(int j=0; j<M; j++){
            fftwf_plan p1;
            fftwf_complex *p = (fftwf_complex*)(&pImageData[j*N]);

            //#pragma omp critical (make_plan)
            p1 = fftwf_plan_dft_1d(N,  p, (fftwf_complex*)p, FFTW_FORWARD,
       FFTW_ESTIMATE );

            fftwf_execute(p1);
            fftwf_destroy_plan(p1);

            fftwf_cleanup_threads();
        }
    */
    return -1;
}

int CDSP_FFTW::ifft(Ipp32fc *pImageData, int N, int M) {
    int order = log2(N);
    assert(N == pow(2, order)); // N должно быть степенью двойки

    int nbthreads = omp_get_max_threads();
    fftwf_init_threads();
    fftwf_plan_with_nthreads(nbthreads);
    omp_set_num_threads(nbthreads);

    fftwf_plan p;
    int rank = 1;
    int n[] = {(int)N};

#pragma omp critical
    p = fftwf_plan_many_dft(rank, n, M, (fftwf_complex *)pImageData, NULL, 1, N,
                            (fftwf_complex *)pImageData, NULL, 1, N,
                            FFTW_BACKWARD, FFTW_ESTIMATE);
    fftwf_execute(p);

#pragma omp critical
    fftwf_destroy_plan(p);

    // #pragma omp parallel for
    //  for(int j=0; j<M; j++){
    //      fftwf_plan p1;
    //      fftwf_complex *p = (fftwf_complex*)(&pImageData[j*N]);

    // #pragma omp critical (make_plan)
    //      p1 = fftwf_plan_dft_1d(N,  p, (fftwf_complex*)p, FFTW_BACKWARD,
    //      FFTW_ESTIMATE );

    //      fftwf_execute(p1);
    //      fftwf_destroy_plan(p1);
    //  }

    return -1;
}

int CDSP_FFTW::bp(CImage *image, CRadar *radar, float *outdata, float *wnd,
                  float *V, float x0, float y0, float dx, float dy, float *h,
                  int nx, int ny, bool bCorr, Ipp32fc *pc) {
    int cnt = 0;

    progress = 0;
#pragma omp parallel for
    for (int i = 0; i < ny; i++) {
        float y = y0 + i * dy;
        for (int j = 0; j < nx; j++) {
            float x = x0 + j * dx;
            outdata[i * nx + j] =
                bp_pixel(image, radar, wnd, V, x, y, h, bCorr, pc);
        }
        cnt++;
        progress = cnt * 100 / ny;
        printf("<%3d%%>\e[6D", progress);
        fflush(stdout);
    }

    return 0;
}

int CDSP_FFTW::bp_progressive(CImage *image, CRadar *radar, float *outdata,
                              float *wnd, float *V, float x0, float y0,
                              float dx, float dy, float *h, int nx, int ny,
                              bool bCorr, Ipp32fc *pc) {
    int cnt = 0;

    progress = 0;

    int npmax = image->ImageH0;
    float k = (radar->GetLambda() * radar->Fs) /
              (2 * mean(V, npmax) * dx * radar->Ns);

    for (int j = 0; j < nx; j++) {
        float x = x0 + j * dx;
        int np = x * k;

        if (image->ImageH0 != np) {
            HammingWindow(wnd, np);

            if (np > npmax) {
                image->ImageH0 = npmax;
            } else {
                image->ImageH0 = np;
            }
        }

#pragma omp parallel for
        for (int i = 0; i < ny; i++) {
            float y = y0 + i * dy;
            outdata[i * nx + j] =
                bp_pixel(image, radar, wnd, V, x, y, h, bCorr, pc);
        }
        cnt++;
        progress = cnt * 100 / nx;
        printf("<%3d%%>\e[6D", progress);
        fflush(stdout);
    }

    return 0;
}

int CDSP_FFTW::bpm(CImage *image, CRadar *radar, float *outdata, float *wnd,
                   float *V, float x0, float y0, float dx, float dy, float *h,
                   int nx, int ny, bool bCorr, float angle, float div,
                   Ipp32fc *pc) {
    int cnt = 0;

    float div2 = div / 2;

    float cosAngle = cos(angle);
    float sinAngle = sin(angle);

#pragma omp parallel for
    for (int j = 0; j < nx; j++) {
        float r = x0 + j * dx;

        int b = (cone_edge(r, -div2) - y0) / dy;
        int e = (cone_edge(r, div2) - y0) / dy;
        if (b < 0) {
            b = 0;
        }
        if (e > ny) {
            e = ny;
        }

        for (int i = b; i < e; i++) {

            float x = (x0 + j * dx) * cosAngle - (y0 + i * dy) * sinAngle;
            float y = (x0 + j * dx) * sinAngle + (y0 + i * dy) * cosAngle;

            outdata[i * nx + j] =
                bp_pixel(image, radar, wnd, V, x, y, h, bCorr, pc);
        }
        cnt++;
        printf("<%3d%%>\e[6D", cnt * 100 / nx);
        fflush(stdout);
    }

    return 0;
}

int CDSP_FFTW::bpm_progressive(CImage *image, CRadar *radar, float *outdata,
                               float *wnd, float *V, float x0, float y0,
                               float dx, float dy, float *h, int nx, int ny,
                               bool bCorr, float angle, float div,
                               Ipp32fc *pc) {
    int cnt = 0;

    float div2 = div / 2;

    float cosAngle = cos(angle);
    float sinAngle = sin(angle);

    int npmax = image->ImageH0;
    float k = (radar->GetLambda() * radar->Fs) /
              (2 * mean(V, npmax) * dx * radar->Ns);

    for (int j = 0; j < nx; j++) {
        float r = x0 + j * dx;

        int np = r * k;

        if (image->ImageH0 != np) {

            if (np > npmax) {
                image->ImageH0 = npmax;
            } else {
                image->ImageH0 = np;
            }

            HammingWindow(wnd, image->ImageH0);
        }

        int b = (cone_edge(r, -div2) - y0) / dy;
        int e = (cone_edge(r, div2) - y0) / dy;
        if (b < 0) {
            b = 0;
        }
        if (e > ny) {
            e = ny;
        }

#pragma omp parallel for
        for (int i = b; i < e; i++) {

            float x = (x0 + j * dx) * cosAngle - (y0 + i * dy) * sinAngle;
            float y = (x0 + j * dx) * sinAngle + (y0 + i * dy) * cosAngle;

            outdata[i * nx + j] =
                bp_pixel(image, radar, wnd, V, x, y, h, bCorr, pc);
        }
        cnt++;
        printf("<%3d%%>\e[6D", cnt * 100 / nx);
        fflush(stdout);
    }

    return 0;
}

int CDSP_FFTW::bp_strip(CImage *in, CImage *out, CRadar *radar, float *window,
                        float *V, float x0, float y0, float dx, float *h,
                        int nx, int line, bool bCorr, Ipp32fc *pc,
                        float angle) {

    float *outdata = (float *)out->pImageData + nx * line;
    // memset(out->pImageData + nx*line, 0, nx*sizeof(float));
    // Ipp32f *outdata = (Ipp32f*)out->pImageData + nx*line;

    float cosAngle = cos(angle);
    float tanAngle = tan(angle);

#pragma omp parallel for
    for (int j = 0; j < nx; j++) {
        float x = (x0 + j * dx) * cosAngle;
        float y = x * tanAngle;

        outdata[j] = bp_pixel(in, radar, window, V, x, y, h, bCorr, pc);
    }

    // free(window);
    return 0;
}

float CDSP_FFTW::bp_pixel(CImage *image, CRadar *radar, float *wnd, float *V,
                          float x, float y, float *h, bool bCorr, Ipp32fc *pc) {
    Ipp32fc *data = (Ipp32fc *)image->pImageData;
    int Np = image->ImageH0;
    int Ns = image->ImageWidth;

    float Tp = radar->GetTp();
    float Mu = radar->GetMu();
    float F0 = radar->F0;
    float Fs = radar->Fs;

    float dt = Tp;
    float c1 = -4 * M_PI / SPEED_OF_LIGHT * F0;
    float c2 = -M_PI * Tp;
    float c5 = 2 * Mu / SPEED_OF_LIGHT;
    float dfr = Fs / (2 * (Ns - 1));
    float dr = Fs * SPEED_OF_LIGHT / (4 * Mu * (Ns - 1));

    Ipp32fc *signal = (Ipp32fc *)malloc(Np * 2 * sizeof(float));
    Ipp32fc *sop = (Ipp32fc *)malloc(Np * 2 * sizeof(float));
    Ipp32fc *sk = (Ipp32fc *)malloc(Np * 2 * sizeof(float));

    // #pragma omp parallel for
    // #pragma block_loop
    for (int k = 0; k < Np; k++) {
        float t = k * dt;
        float c4 = powf(x, 2) + powf(h[k], 2);
        float r = sqrtf(c4 + powf(y + V[k] * t, 2));

#ifdef X86MATH
        int ind = x86Round(r / dr);
#else
        int ind = roundf(r / dr);
#endif

        //signal[k] = data[k * Ns + ind];
        signal[k] = data[(k * Ns + ind + image->actualIndex / image->GetPixelSize(image->ImageType))%(image->ImageHeight*image->ImageWidth)];
    
        float phase = c1 * r;

#ifdef X86MATH
        x86SinCos(phase, &sop[k].im, &sop[k].re);
#else
        sop[k].re = cosf(phase);
        sop[k].im = sinf(phase);
#endif

        // Коррекция скачков фазы
        if (bCorr) {
            float fp = (ind - 1) * dfr;
            // float fr = range2fr(r);
            float fr = r * c5;
            float phi = c2 * (fr - fp);

#ifdef X86MATH
            x86SinCos(phi, &sk[k].im, &sk[k].re);
#else
            sk[k].re = cosf(phi);
            sk[k].im = sinf(phi);
#endif
        }
    }

    if (bCorr) {
        complexMul(signal, sk, Np);
    }

    free(sk);

    complexMul(signal, wnd, Np); // Оконное взвешивание
    complexMul(signal, sop, Np);
    free(sop);

    if (pc) {
        complexMul(signal, pc, Np);
    }

    Ipp32fc Sum;
    Ipp32f fSum = 0;
    complexSum(signal, &Sum, Np);
    fSum = sqrtf(Sum.im * Sum.im + Sum.re * Sum.re);

    free(signal);

    return fSum;
}
/*
int CDSP_FFTW::bp_strip(CLinkedImage *in, CImage *out, CRadar *radar,
                        float *window, float *V, float x0, float y0, float dx,
                        float *h, int nx, int line, bool bCorr, Ipp32fc *pc,
                        float angle) {

    float *outdata = (float *)out->pImageData + nx * line;

    float cosAngle = cos(angle);
    float tanAngle = tan(angle);

    //    #pragma omp parallel for
    for (int j = 0; j < nx; j++) {
        float x = (x0 + j * dx) * cosAngle;
        float y = x * tanAngle;

        outdata[j] = bp_pixel(in, radar, window, V, x, y, h, bCorr, pc);
    }

    return 0;
}

float CDSP_FFTW::bp_pixel(CLinkedImage *image, CRadar *radar, float *wnd,
                          float *V, float x, float y, float *h, bool bCorr,
                          Ipp32fc *pc) {
    //  Ipp32fc* data = (Ipp32fc*)image->pImageData;
    int Np = image->ImageH0;
    int Ns = image->ImageWidth;

    float Tp = radar->GetTp();
    float Mu = radar->GetMu();
    float F0 = radar->F0;
    float Fs = radar->Fs;

    float dt = Tp;
    float c1 = -4 * M_PI / SPEED_OF_LIGHT * F0;
    float c2 = -M_PI * Tp;
    float c5 = 2 * Mu / SPEED_OF_LIGHT;
    float dfr = Fs / (2 * (Ns - 1));
    float dr = Fs * SPEED_OF_LIGHT / (4 * Mu * (Ns - 1));

    Ipp32fc *signal = (Ipp32fc *)malloc(Np * 2 * sizeof(float));
    Ipp32fc *sop = (Ipp32fc *)malloc(Np * 2 * sizeof(float));
    Ipp32fc *sk = (Ipp32fc *)malloc(Np * 2 * sizeof(float));

    // #pragma omp parallel for
    // #pragma block_loop
    auto iter = image->data.into_iter();

    printf("x, y = %f, %f\n", x, y);
    // printf("image->data = %d\n", image->data.size());
    for (int k = 0; k < Np; k++) {
        float t = k * dt;
        float c4 = powf(x, 2) + powf(h[k], 2);
        float r = sqrtf(c4 + powf(y + V[k] * t, 2));

#ifdef X86MATH
        int ind = x86Round(r / dr);
#else
        int ind = roundf(r / dr);
#endif

        Ipp32fc *data = iter.peek();
        iter.next();
        printf("Next ok %p\n", data);
        signal[k] = data[Ns + ind];
        printf("signal ok\n");
        float phase = c1 * r;

#ifdef X86MATH
        x86SinCos(phase, &sop[k].im, &sop[k].re);
#else
        sop[k].re = cosf(phase);
        sop[k].im = sinf(phase);
#endif

        // Коррекция скачков фазы
        if (bCorr) {
            float fp = (ind - 1) * dfr;
            // float fr = range2fr(r);
            float fr = r * c5;
            float phi = c2 * (fr - fp);

#ifdef X86MATH
            x86SinCos(phi, &sk[k].im, &sk[k].re);
#else
            sk[k].re = cosf(phi);
            sk[k].im = sinf(phi);
#endif
        }
    }

    if (bCorr) {
        complexMul(signal, sk, Np);
    }

    free(sk);

    complexMul(signal, wnd, Np); // Оконное взвешивание
    complexMul(signal, sop, Np);
    free(sop);

    if (pc) {
        complexMul(signal, pc, Np);
    }

    Ipp32fc Sum;
    Ipp32f fSum = 0;
    complexSum(signal, &Sum, Np);
    fSum = sqrtf(Sum.im * Sum.im + Sum.re * Sum.re);

    free(signal);

    return fSum;
}
*/
int CDSP_FFTW::dspProgress() { return progress; }
