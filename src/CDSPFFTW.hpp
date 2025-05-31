/*
 * CDSPCPU.hpp
 *
 *  Created on: Nov 10, 2019
 *      Author: user
 */

#ifndef CDSPCPU_H_
#define CDSPCPU_H_

#include "CDSP.hpp"
#include "CLog.hpp"
#include "fftw3.h"

class CDSP_FFTW : public CDSP {
  public:
    CDSP_FFTW();
    CDSP_FFTW(CLog *log);
    virtual ~CDSP_FFTW();

    // int fft_r2c(float* inData, fftwf_complex* outData, int N, int M);
    int fft(Ipp32f *pImageData, int N, int M);
    int fft(Ipp32fc *pImageData, int N, int M);
    int ifft(Ipp32fc *pImageData, int N, int M);

    int bp(CImage *image, CRadar *radar, float *outdata, float *wnd, float *V,
           float x0, float y0, float dx, float dy, float *h, int nx, int ny,
           bool bCorr, Ipp32fc *pc);
    int bp_progressive(CImage *image, CRadar *radar, float *outdata, float *wnd,
                       float *V, float x0, float y0, float dx, float dy,
                       float *h, int nx, int ny, bool bCorr, Ipp32fc *pc);
    int bpm(CImage *image, CRadar *radar, float *outdata, float *wnd, float *V,
            float x0, float y0, float dx, float dy, float *h, int nx, int ny,
            bool bCorr, float angle, float div, Ipp32fc *pc);
    int bpm_progressive(CImage *image, CRadar *radar, float *outdata,
                        float *wnd, float *V, float x0, float y0, float dx,
                        float dy, float *h, int nx, int ny, bool bCorr,
                        float angle, float div, Ipp32fc *pc);
    int bp_strip(CImage *in, CImage *out, CRadar *radar, float *window,
                 float *V, float x0, float y0, float dx, float *h, int nx,
                 int line, bool bCorr, Ipp32fc *pc, float angle);
    /*
    int bp_strip(CLinkedImage *in, CImage *out, CRadar *radar, float *window,
                 float *V, float x0, float y0, float dx, float *h, int nx,
                 int line, bool bCorr, Ipp32fc *pc, float angle);

    float bp_pixel(CLinkedImage *image, CRadar *radar, float *wnd, float *V,
                   float x, float y, float *h, bool bCorr, Ipp32fc *pc);
*/
    int dspProgress();

  private:
    fftwf_plan FFTPlan;
    float bp_pixel(CImage *image, CRadar *radar, float *wnd, float *V, float x,
                   float y, float *h, bool bCorr, Ipp32fc *pc);
};

#endif /* CDSPCPU_H_ */
