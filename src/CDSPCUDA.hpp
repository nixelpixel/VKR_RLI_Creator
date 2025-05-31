/*
 * CDSPCUDA.hpp
 *
 *  Created on: 01 янв. 2020 г.
 *      Author: user
 */

#ifndef CDSPCUDA_H_
#define CDSPCUDA_H_

#include "CDSP.hpp"
#include "CLog.hpp"
#include "fftw3.h"
#ifdef QUASAR_CUDA
#include <cuda.h>
#include <cuda_runtime.h>
#include <cufft.h>
// #include <cuda_runtime_api.h>
// #include <driver_types.h>
// #define __CUDA_RUNTIME_H__
// #define   __DRIVER_TYPES_H__
// #include <CUDA_DRIVER_API.h>
#include "./cuda/helper_cuda.h"
#endif
// #include <omp.h>

#define MAX_BUF_SIZE 512

class CDSP_CUDA : public CDSP {
  public:
    CDSP_CUDA(CLog *log);
    CDSP_CUDA();
    virtual ~CDSP_CUDA();
#ifdef QUASAR_CUDA
    int printDSPInfo(void); // Вывод информации о параметрах системы
    int fft(Ipp32f *pImageData, int N, int M);
    int fft(Ipp32fc *pImageData, int N, int M);
    int ifft(Ipp32fc *pImageData, int N, int M);

    int bp(CImage *image, CRadar *radar, float *outdata, float *wnd, float *V,
           float x0, float y0, float dx, float dy, float *h, int nx, int ny,
           bool bCorr, Ipp32fc *pc);
    int bpm(CImage *image, CRadar *radar, float *outdata, float *wnd, float *V,
            float x0, float y0, float dx, float dy, float *h, int nx, int ny,
            bool bCorr, float angle, float div, Ipp32fc *pc);
    int bp_strip(CImage *in, CImage *out, CRadar *radar, float *window,
                 float *V, float x0, float y0, float dx, float *h, int nx,
                 int line, bool bCorr, Ipp32fc *pc, float angle);

  private:
    int square(Ipp32f *pInData, Ipp32f *pOutData, int N, int M);
    int log10(Ipp32f *pInData, Ipp32f *pOutData, int N, int M);

    int devID;
    CUdevice cuDevice;
    cudaDeviceProp prop;

    // float bp_pixel(CImage* image, CRadar *radar, float *wnd, float* V, float
    // x, float y, float* h, bool bCorr);

  private:
    void printCPUFeatures(void); // Вывод информации о параметрах процессора
#endif
};

#endif /* CDSPCUDA_H_ */
