#ifndef CBP_H_
#define CBP_H_

#include "CDSP.hpp"
#include "CDSPCUDA.hpp"
#include "CDSPFFTW.hpp"
#include "CDSPIPP.hpp"
#include "CDSPNEON.hpp"
#include "CLog.hpp"
#include "CRadar.hpp"
#include "ipp/ipp.h"
#include <math.h>
#include <pybind11/pybind11.h>
#include <stdio.h>
#include <string.h>
// #include "CSim.hpp"
#include "CImage.hpp"
#include "CSAR.hpp"
#include "base.hpp"
#include "helpers.hpp"

class CBackProjection {
  public:
    CBackProjection();
    CBackProjection(CSAR *sar, CImage *image);
    virtual ~CBackProjection();

    void setPixelSize(float pixSize);
    void setPixelSize(float _dx, float _dy);
    int setVH(float v, float h, int H0);
    int setVH(pybind11::list &v, pybind11::list &h, int H0);
    int setVH(pybind11::list &v, pybind11::list &h);
    int setWindowType(CSAR::WindowType w_type);
    int setProgressive(bool p);

    int telescopic(float x0, float lx, float y0, float ly);
    int telescopicAngle(float x0, float lx, float angle, float div);

    int phaseCorrection(float x, float y);

    int stripInit(float lx, float Tstrip, float Ts, pybind11::list &v,
                  pybind11::list &h);
    int strip(float x0, float y0, int n, float angle);
    CImage *getStripImage();

    float autofocus(float v0, float v1, float point_x, float point_y, float Ts,
                    float ang = 0.0, float ls = 40.0, int n = 20,
                    char *report = nullptr);

    // int setConsts(float x0_, float y0_, float dx_, int nx_, bool bCorr_);
    // int bp_strip(CImage* AccImage, CImage* out, int n);
    // int createWindow(int size, CSAR::WindowType w_type);

  private:
    CImage *Image = nullptr;

    float *V;
    float *H;
    float *Window;

    Ipp32fc *pc = NULL;

    float dx = 1.0;
    float dy = 1.0;
    CSAR::WindowType W_type = CSAR::WindowType::HAMMING;
    bool bCorr = true;
    bool progressive = false;

    CSAR *SAR;

    // strip
    CImage *AccImage = nullptr;
    CImage *StripImage = nullptr;

    int snx;
};

#endif /* CBP_H_ */
