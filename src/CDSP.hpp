/*
 * CDSP.hpp
 *
 *  Created on: 14 мая 2017 г.
 *      Author: user
 */

#ifndef CDSP_H_
#define CDSP_H_

#include <math.h>
#include <string.h>
#ifdef QUASAR_IPP
#include <xmmintrin.h> //SSE
#endif
#include "CImage.hpp"
#include "CLinkedImage.hpp"
#include "CLog.hpp"
#include "CRadar.hpp"
#include "helpers.hpp"
#include "ipp/ipp.h"
#include <assert.h>

// #include <arrayfire.h>

#define SPEED_OF_LIGHT 300000000.0f

#define check_sts(st)                                                          \
    if ((st) != ippStsNoErr)                                                   \
        goto exitLine; /* Go to Exit if IPP function returned status different \
                          from ippStsNoErr */

void transpose_scalar_block(Ipp32fc *A, Ipp32fc *B, const int lda,
                            const int ldb, const int block_size);
void transpose4x4_SSE(Ipp32fc *A, Ipp32fc *B, const int lda, const int ldb);

/**
 \brief Виртуальный класс, предоставляющий функции системы ЦОС.
 */

class CDSP {
  public:
    CDSP();
    CDSP(CLog *log);
    virtual ~CDSP();

    // Переменные
    CLog *Log;

    // Функции
#ifdef SUPPORT_AF
    virtual af::array af_mean(const af::array &X) { return 0; };
    virtual int af_fft_r2c(af::array &X) { return 0; };
    virtual int af_fft(af::array &X) { return 0; };
    virtual int af_hamming(af::array &X, int W0) { return 0; };
    virtual int af_transpose(af::array &X) { return 0; };
    virtual int af_fftshift(af::array &X) { return 0; };
    virtual int af_rmc(af::array &X, int hb, float lambda, float tp, float v) {
        return 0;
    };
    virtual int af_ifft(af::array &X) { return 0; };
    virtual int af_phocus(af::array &X, int Hmax, float Rb, float Tp,
                          float lambda, float V, float Rmax, int Np) {
        return 0;
    };
    virtual int af_complexMul(af::array &X, af::array &S, af::seq i) const {
        return 0;
    };
    virtual float af_entropy(af::array &X) { return 0; };
    virtual int af_compensateShift(af::array &X, float mp) { return 0; };
#endif

    virtual int printDSPInfo(void); // Вывод информации о параметрах системы
    virtual float mean(Ipp32f *pImageData, int N);
    virtual int removeDC(Ipp32fc *pImageData, int W, int W0, int H);
    virtual int hamming(Ipp32f *pImageData, int W, int H, int W0, int H0);
    virtual int hamming(Ipp32fc *pImageData, int W, int H, int W0, int H0);
    virtual int fft(Ipp32f *pImageData, int N, int M) { return 0; };
    virtual int fft(Ipp32fc *pImageData, int N, int M) { return 0; };
    virtual int ifft(Ipp32fc *pImageData, int N, int M) { return 0; };
    virtual int idft(Ipp32fc *pImageData, int N, int M) { return 0; };

    virtual int fftshift(Ipp32fc *pImageData, int N, int M);
    virtual int transpose(Ipp32fc *pImageData, int N, int M);
    virtual int transpose(Ipp32f *pImageData, int N, int M) { return 0; };
    virtual int compensateShift(Ipp32fc *pImageData, int W, int H, float mp,
                                int PixelSize) {
        return 0;
    };
    virtual int magnitude(Ipp32fc *pInData, Ipp32f *pOutData, int N, int M) {
        return 0;
    };
    //  virtual int     af_magnitude(af::array* afImageData, af::array magn, int
    //  N, int M){return 0;};
    virtual int rmc(Ipp32fc *pData, int w, int h, int hb, float lambda,
                    float tp, float v);
    virtual int phocus(Ipp32fc *pData, int W, int H, int Hmax, float Rb,
                       float Tp, float lambda, float V, float Rmax, int Np);
    virtual float entropy(Ipp32fc *pData, int W, int H);
    virtual float entropyf(Ipp32f *pData, int W, int H);
    virtual int complexMul(Ipp32fc *SrcDst1, Ipp32fc *Src2, int len);
    virtual int complexMul(Ipp32fc *SrcDst1, float *Src2, int len);
    virtual int complexSum(Ipp32fc *Src, Ipp32fc *Dst, int len);

    virtual int hilbert(Ipp32f *pInData, Ipp32fc *pOutData, int W, int H) {
        return 0;
    };
    virtual float bp_slow_pixel(Ipp32fc *data, int W, int H, float V, float Rx,
                                float Ry, float F0, float Tp, float Mu) {
        return 0;
    };
    virtual int pointTargetSignal(void *data, int Ns, int Fs, int Np, float *H,
                                  float *V, float X, float Y, float A, float F0,
                                  float DF, float Tp, float Mu, bool complex,
                                  bool symmetric_ramp);
    virtual int fastconv(Ipp32fc *srcdst, Ipp32fc *src, int W) { return 0; };

    virtual int HammingWindow(float *pImageData, int W0);
    virtual int BlackmanWindow(float *pImageData, int W0);
    virtual int NuttallWindow(float *pImageData, int W0);

    virtual int windowWeighting(float *pImageData, float *window, int W, int W0,
                                int H0);

    virtual int windowWeighting(Ipp32fc *pImageData, float *window, int W,
                                int W0, int H0);

    virtual int bp(CImage *image, CRadar *radar, float *outdata, float *wnd,
                   float *V, float x0, float y0, float dx, float dy, float *h,
                   int nx, int ny, bool bCorr, Ipp32fc *pc);
    virtual int bp_progressive(CImage *image, CRadar *radar, float *outdata,
                               float *wnd, float *V, float x0, float y0,
                               float dx, float dy, float *h, int nx, int ny,
                               bool bCorr, Ipp32fc *pc);
    virtual int bpm(CImage *image, CRadar *radar, float *outdata, float *wnd,
                    float *V, float x0, float y0, float dx, float dy, float *h,
                    int nx, int ny, bool bCorr, float angle, float div,
                    Ipp32fc *pc);
    virtual int bpm_progressive(CImage *image, CRadar *radar, float *outdata,
                                float *wnd, float *V, float x0, float y0,
                                float dx, float dy, float *h, int nx, int ny,
                                bool bCorr, float angle, float div,
                                Ipp32fc *pc);
    virtual int bp_strip(CImage *in, CImage *out, CRadar *radar, float *window,
                         float *V, float x0, float y0, float dx, float *h,
                         int nx, int line, bool bCorr, Ipp32fc *pc,
                         float angle);

    /*
    virtual int bp_strip(CLinkedImage *in, CImage *out, CRadar *radar,
                         float *window, float *V, float x0, float y0, float dx,
                         float *h, int nx, int line, bool bCorr, Ipp32fc *pc,
                         float angle);
*/
    virtual int dspProgress() { return 0; };

    int progress = 0;

    int cone_edge(float x, float phi);

  private:
    float bp_pixel(CImage *image, CRadar *radar, float *wnd, float *V, float x,
                   float y, float *h, bool bCorr, Ipp32fc *pc);
};

#endif /* CDSP_H_ */
