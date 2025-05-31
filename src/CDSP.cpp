/*
 * CDSP.cpp
 *
 *  Created on: 14 мая 2017 г.
 *      Author: user
 */

#include "CDSP.hpp"
#include "helpers.hpp"

CDSP::CDSP(CLog *log) : Log(log) {}

CDSP::CDSP() { Log = NULL; }

CDSP::~CDSP() {
    // TODO Auto-generated destructor stub
}

int CDSP::pointTargetSignal(void *data, int Ns, int Fs, int Np, float *H,
                            float *V, float X, float Y, float A, float F0,
                            float DF, float Tp, float Mu, bool complex,
                            bool symmetric_ramp) {
    float c1 = 1.0f / Fs;

    float f0 = F0;
    float mu = Mu;

    // #pragma omp parallel for  // Почему-то сигнал рваный получается
    for (int i = 0; i < Np; i++) {
        if (symmetric_ramp && (i & 1)) {
            f0 = F0 + DF;
            Mu = -Mu;
        } else {
            f0 = F0;
            mu = Mu;
        }

        float t, R, p1, p2, p3, phase;
        float c2 = pow(X, 2) + pow(H[i], 2);
        for (int j = 0; j < Ns; j++) {
            // Вычисление фазы опорной функции для каждого периода
            t = (float)j * c1; // Отсчеты времени в пределах периода
            R = sqrt(pow((Y + V[i] * (t + i * Tp)), 2) + c2);
            p1 = f0;
            p2 = mu * t;
            p3 = -mu * R / SPEED_OF_LIGHT;
            phase = 4 * M_PI * R / SPEED_OF_LIGHT * (p1 + p2 + p3);
            if (!complex)
                ((Ipp32f *)data)[i * Ns + j] += (A * cos(phase) + 1) / 0.5f;
            else {
                (((Ipp32fc *)data)[i * Ns + j]).re +=
                    (A * cos(phase) + 1) / 0.5f;
                (((Ipp32fc *)data)[i * Ns + j]).im +=
                    (A * sin(phase) + 1) / 0.5f;
            }
        }
    }
    return 0;
}

int CDSP::printDSPInfo(void) {
    Log->LogPrintf((char *)"============\n");
    
#ifdef X86MATH
    Log->LogPrintf("CPU x86_64\n");
#endif

    std::vector<std::string> tags;
    std::vector<std::string> tags_rus;


    tags.push_back("model name");
    tags.push_back("siblings");
    
    tags_rus.push_back("Модель процессора");
    tags_rus.push_back("Количество потоков");

    tags = sys_get_cpu_info(tags);

    for (int i = 0; i < tags.size(); i++) {
        // Log->LogPrintf("%s:\t%s\n", tags_rus[i].c_str(), tags[i].c_str());
        printf("%s:\t%s\n", tags_rus[i].c_str(), tags[i].c_str());
    }
    
    return 0;

}

float CDSP::mean(Ipp32f *pImageData, int N) {
    float m = 0;
    for (int i = 0; i < N; i++) {
        m += pImageData[i];
    }
    return m / N;
}

int CDSP::removeDC(Ipp32fc *pImageData, int W0, int W, int H) {

    Ipp32fc summ;

    for (int x = 0; x < W; x++) {

        summ.im = 0;
        summ.re = 0;
        for (int y = 0; y < H; y++) {

            summ.im += pImageData[W0 * y + x].im;
            summ.re += pImageData[W0 * y + x].re;
        }

        summ.im /= H;
        summ.re /= H;

        for (int y = 0; y < H; y++) {

            pImageData[W0 * y + x].im -= summ.im;
            pImageData[W0 * y + x].re -= summ.re;
        }
    }

    return 0;
}

//============================
// fftshift()
//============================
int CDSP::fftshift(Ipp32fc *pImageData, int N, int M) {
#pragma omp parallel for
    for (int j = 0; j < M; j++) {
        // fftshift
        Ipp32fc *pCopyBuffer = (Ipp32fc *)malloc(N / 2 * sizeof(Ipp32fc));
        memcpy(pCopyBuffer, &pImageData[j * N], N / 2 * sizeof(Ipp32fc));
        memcpy(&pImageData[j * N], &pImageData[j * N + N / 2],
               N / 2 * sizeof(Ipp32fc));
        memcpy(&pImageData[j * N + N / 2], pCopyBuffer,
               N / 2 * sizeof(Ipp32fc));
        free(pCopyBuffer);
    }

    return 0;
}

//==================================
// rmc()
//==================================
int CDSP::rmc(Ipp32fc *pData, int W, int H, int Hb, float Lambda, float Tp,
              float V) {
    float Fdmax = 1 / (2.0f * Tp);
    float dFd = 2 * Fdmax / W;
    float a = Lambda / (2.0f * V);

#pragma omp parallel for
    for (int w = 0; w < W; w++) {
        float fd = -Fdmax + w * dFd;
        float k = sqrt(1.0 - pow(a * fd, 2));

        if (k > 1)
            k = 1;

        for (int h = 0; h < H; h++) {
            // int index = round(h/k);
            int index = round((h + Hb) / k - Hb);
            if (index >= H)
                continue;

            pData[h * W + w] = pData[index * W + w];
        }
    }

    return 0;
}

int CDSP::hamming(float *pImageData, int W, int H, int W0, int H0) {
    float *window = (float *)malloc(W0 * sizeof(float));

#pragma omp parallel for
    for (int i = 0; i < W0; i++) {
        window[i] = 0.54f - 0.46f * cos(2 * M_PI * i / W0);
    }

#pragma omp parallel for
    for (int i = 0; i < H0; i++) {
        for (int j = 0; j < W0; j++) {
            pImageData[i * W + j] *= window[j];
        }
    }

    free(window);
    return 0;
}

int CDSP::hamming(Ipp32fc *pImageData, int W, int H, int W0, int H0) {
    float *window = (float *)malloc(W0 * sizeof(float));

#pragma omp parallel for
    for (int i = 0; i < W0; i++) {
        window[i] = 0.54f - 0.46f * cos(2 * M_PI * i / W0);
    }

#pragma omp parallel for
    for (int i = 0; i < H0; i++) {
        for (int j = 0; j < W0; j++) {
            pImageData[i * W + j].im *= window[j];
            pImageData[i * W + j].re *= window[j];
        }
    }

    free(window);
    return 0;
}

int CDSP::complexMul(Ipp32fc *SrcDst1, Ipp32fc *Src2, int len) {

    for (int i = 0; i < len; i++) {
        float re = SrcDst1[i].re * Src2[i].re - SrcDst1[i].im * Src2[i].im;
        float im = SrcDst1[i].re * Src2[i].im + SrcDst1[i].im * Src2[i].re;
        SrcDst1[i].re = re;
        SrcDst1[i].im = im;
    }

    return 0;
}

int CDSP::complexMul(Ipp32fc *SrcDst1, float *Src2, int len) {
    for (int i = 0; i < len; i++) {
        SrcDst1[i].re = SrcDst1[i].re * Src2[i];
        SrcDst1[i].im = SrcDst1[i].im * Src2[i];
    }

    return 0;
}

int CDSP::complexSum(Ipp32fc *Src, Ipp32fc *Dst, int len) {

    Dst[0].re = 0;
    Dst[0].im = 0;
    for (int i = 0; i < len; i++) {
        Dst[0].re += Src[i].re;
        Dst[0].im += Src[i].im;
    }

    return 0;
}

int CDSP::HammingWindow(float *pImageData, int W0) {

#pragma omp parallel for
    for (int i = 0; i < W0; i++) {
        pImageData[i] = 0.54f - 0.46f * cosf(2 * M_PI * i / W0);
    }

    return 0;
}

int CDSP::BlackmanWindow(float *pImageData, int W0) {

    const float a0 = 0.42f;
    const float a1 = 0.5f;
    const float a2 = 0.08f;

#pragma omp parallel for
    for (int i = 0; i < W0; i++) {
        pImageData[i] =
            a0 - a1 * cosf(2 * M_PI * i / W0) + a2 * cosf(4 * M_PI * i / W0);
    }

    return 0;
}

int CDSP::NuttallWindow(float *pImageData, int W0) {

    const float a0 = 0.355768f;
    const float a1 = 0.487396f;
    const float a2 = 0.134232f;
    const float a3 = 0.012604f;

#pragma omp parallel for
    for (int i = 0; i < W0; i++) {
        pImageData[i] = a0 - a1 * cos(2 * M_PI * i / W0) +
                        a2 * cosf(4 * M_PI * i / W0) -
                        a3 * cosf(6 * M_PI * i / W0);
    }

    return 0;
}

int CDSP::windowWeighting(float *pImageData, float *window, int W, int W0,
                          int H0) {

#pragma omp parallel for
    for (int i = 0; i < H0; i++) {
        for (int j = 0; j < W0; j++) {
            pImageData[i * W + j] *= window[j];
        }
    }
    return 0;
}

int CDSP::windowWeighting(Ipp32fc *pImageData, float *window, int W, int W0,
                          int H0) {
#pragma omp parallel for
    for (int i = 0; i < H0; i++) {
        for (int j = 0; j < W0; j++) {
            pImageData[i * W + j].im *= window[j];
            pImageData[i * W + j].re *= window[j];
        }
    }
    return 0;
}

//==================================
// entropyf() Вычисляет энергию вектора
//==================================
float CDSP::entropyf(float *pData, int W, int H) {
    float P = 0;

    // #pragma omp parallel for
    for (int i = 0; i < W * H; i++) {
        if (pData[i]) {
            P += pData[i] * pData[i];
        }
    }

    float E = 0;
    // #pragma omp parallel for
    for (int i = 0; i < W * H; i++) {
        if (pData[i]) {
            float p = pData[i] * pData[i] / P;
            E += p * log(p);
        }
    }

    return -E;
}

float CDSP::entropy(Ipp32fc *pData, int W, int H) {
    float P = 0;
    // #pragma omp parallel for
    for (int i = 0; i < W * H; i++) {
        if (pData[i].re && pData[i].im) {
            float m =
                sqrt(pData[i].re * pData[i].re + pData[i].im * pData[i].im);
            P += m * m;
        }
    }

    float E = 0;
    // #pragma omp parallel for
    for (int i = 0; i < W * H; i++) {
        if (pData[i].re && pData[i].im) {
            float m =
                sqrt(pData[i].re * pData[i].re + pData[i].im * pData[i].im);
            float p = m * m / P;
            E += p * log(p);
        }
    }

    return -E;
}

//==================================
// transpose()
//==================================
int CDSP::transpose(Ipp32fc *pImageData, int N, int M) {
    int buflen = N * M * sizeof(Ipp32fc);
    Ipp32fc *buff = (Ipp32fc *)malloc(buflen);

    size_t block = 16;

#pragma omp parallel for
    for (size_t i = 0; i < M; i += block) {
        for (size_t j = 0; j < N; ++j) {
            for (size_t b = 0; b < block && i + b < M; ++b) {
                buff[j * M + i + b] = ((Ipp32fc *)pImageData)[(i + b) * N + j];
            }
        }
    }

    memcpy(pImageData, buff, buflen);
    free(buff);

    return 0;
}

//==================================
// phocus()
//==================================
int CDSP::phocus(Ipp32fc *pData, int W, int H, int Hmax, float Rb, float Tp,
                 float lambda, float V, float Rmax, int Np) {
#pragma omp parallel for
    for (int h = 0; h < H; h++) {
        float r = Rb + Rmax / Hmax * h;
        if (r == 0)
            r = Rb + Rmax / Hmax;
        float c1 = 2 * M_PI / lambda / r;
        // Формирование опорной функции
        Ipp32fc *S = (Ipp32fc *)aligned_alloc(64, W * sizeof(Ipp32fc));
        memset(S, 0, W * sizeof(Ipp32fc)); // Это обязательно
        for (int w = 0; w < Np; w++) {
            float t = w - Np / 2.0f;
            float phase = -c1 * pow(V * t * Tp, 2);

#ifdef X86MATH
            x86SinCos(phase, &S[w].im, &S[w].re);
#else
            S[w].im = sinf(phase);
            S[w].re = cosf(phase);
#endif
        }

        complexMul(&pData[h * W], S, W);
        free(S);
    }

    return 0;
}

float CDSP::bp_pixel(CImage *image, CRadar *radar, float *wnd, float *V,
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

        signal[k] = data[k * Ns + ind];
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

int CDSP::cone_edge(float x, float phi) {

    float r = x / cos(phi);
    float y = r * sin(phi);

    return (int)(y + 0.5);
}

int CDSP::bp(CImage *image, CRadar *radar, float *outdata, float *wnd, float *V,
             float x0, float y0, float dx, float dy, float *h, int nx, int ny,
             bool bCorr, Ipp32fc *pc) {
    Log->LogPrintf((char *)"\n\t*Использована реализация %s по-умолчанию!\n",
                   __func__);
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

int CDSP::bp_progressive(CImage *image, CRadar *radar, float *outdata,
                         float *wnd, float *V, float x0, float y0, float dx,
                         float dy, float *h, int nx, int ny, bool bCorr,
                         Ipp32fc *pc) {
    Log->LogPrintf((char *)"\n\t*Использована реализация %s по-умолчанию!\n",
                   __func__);
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

int CDSP::bpm(CImage *image, CRadar *radar, float *outdata, float *wnd,
              float *V, float x0, float y0, float dx, float dy, float *h,
              int nx, int ny, bool bCorr, float angle, float div, Ipp32fc *pc) {
    Log->LogPrintf((char *)"\n\t*Использована реализация %s по-умолчанию!\n",
                   __func__);
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

int CDSP::bpm_progressive(CImage *image, CRadar *radar, float *outdata,
                          float *wnd, float *V, float x0, float y0, float dx,
                          float dy, float *h, int nx, int ny, bool bCorr,
                          float angle, float div, Ipp32fc *pc) {
    Log->LogPrintf((char *)"\n\t*Использована реализация %s по-умолчанию!\n",
                   __func__);
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

int CDSP::bp_strip(CImage *in, CImage *out, CRadar *radar, float *window,
                   float *V, float x0, float y0, float dx, float *h, int nx,
                   int line, bool bCorr, Ipp32fc *pc, float angle) {

    Log->LogPrintf((char *)"\n\t*Использована реализация %s по-умолчанию!\n",
                   __func__);
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
