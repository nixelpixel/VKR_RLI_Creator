/*
 * CImage.hpp
 *
 *  Created on: 19 июн. 2017 г.
 *      Author: user
 */

#ifndef CIMAGE_H_
#define CIMAGE_H_

#include <stdint.h>
// #include "ipp/ipp.h"
#include "CLog.hpp"
#include "CPlot.hpp"
// #include "CDSPIPP.hpp"
#include "helpers.hpp"
#include <assert.h>
#include <png.h>
#include <stdlib.h>
#include <string.h>
// #include <arrayfire.h>
#include <pybind11/pybind11.h>


#include "vector_to_nparray.hpp"

#include <jpeglib.h>

#define DEF_PLOT_WIDTH 1024
#define DEF_PLOT_HEIGHT 600

/**
 \brief Класс, описывает параметры Радиолокационного изображения и предоставляет
 функции для работы с ним.
 */

class CImage {
  public:
    CImage();
    CImage(CLog *log);
    virtual ~CImage();

    enum eImageType {
        IT_UNDEFINED,
        IT_INT16,
        IT_FLOAT,
        IT_FCOMPLEX
    }; ///< Формат РЛИ
    enum ColorMap { COLORMAP_GRAY = 0, COLORMAP_JET = 1 }; ///< Цветовая схема
    enum CompressType { STANDART, PROGRESSIVE };

    struct metadata_t {
        double lat = 0;
        double lon = 0;
        float dx = 0;
        float dy = 0;
        float x0 = 0;
        float y0 = 0;
        float ang = 0;
        float driftAngle = 0;
        float lx = 0;
        float ly = 0;
        float div = 0;

        float v = 0;
        float h = 0;
        float kr = 0;
        float t = 0;
        float ts = 0;

        float r1 = 0;
        float r2 = 0;

        uint8_t mode = 0;
        uint8_t type = 0;

        uint32_t r4 = 0;
        uint32_t r5 = 0;

        uint16_t crc = 0;
    } __attribute__((packed));

    struct metadata_t metadata;

    CPlot Plot; // Построение графиков

    void *pImageData; ///< Указатель на массив данных РЛИ
#ifdef SUPPORT_AF
    af::array afImageData; ///< Массив данных РЛИ в формате Array Fire
#endif
    int ImageWidth;  ///< Шиина РЛИ
    int ImageHeight; ///< Высота РЛИ
    eImageType ImageType; ///< Тип РЛИ (амплитудное, комплексное и т.д.)
    int ImageW0; ///< Количесво отсчетов (реальное) по ширине изображения без
                 ///< дополнения до степени двойки
    int ImageH0; ///< Количесво отсчетов (реальное) по высоте изображения без
                 ///< дополнения до степени двойки
    int ImageNearRangeIndex; ///< Индекс ближней границы вырезанного изображения
                             ///< на исходном
    int actualIndex;
    // enum      ImageType {IT_INT16 = sizeof(int16_t), IT_F = sizeof(Ipp32f),
    // IT_FC = sizeof(Ipp32fc)}; ///< Формат РЛ

  private:
    CLog *Log;
    //  CDSP    *DSP;

    std::vector<double> pyValues;

  public:
    // int       SaveSignalToPng8(char* filename, int ns);
    // int       SaveSignalToPng8MagNorm(char* filename, int ns, float
    // brightness);
    int SaveToPng8_(char *fileName, int Width, int Height, uint8_t *pData,
                    ColorMap colormap);
    int SaveToJpg8_(char *fileName, int Width, int Height, uint8_t *pData,
                    int quality);
    // int       SaveToPng8MagNorm(char* fileName, float brightness, ColorMap
    // colormap);
    int SaveToPng8(char *fileName, float brightness, ColorMap colormap);

    static void JpegCompress(unsigned char *image, int width, int height,
                             int quality, unsigned long *jpegSize,
                             unsigned char **jpegBuf,
                             CompressType compress_type = PROGRESSIVE);
    int SaveToJpg8(char *fileName, float brightness, int quality);
    int SaveToJpg8(char *fileName, float brightness, int quality, char *label);
    static int SaveToJpg8(char *fileName, int Width, int Height, uint8_t *pData,
                          int quality);

    uint8_t *Prepocessing(uint8_t *savebuff, float brightness, int quality,
                          char *label);
    int SaveToBin(char *fileName);

    int LoadBin(char *fileName, int width, int height, eImageType type);

    png_color GetJetColor(float v);
    int Create(int width, int height, eImageType type);
    int SetPixel(int w, int h, float value);
    void Free(void);
    void CloneFrom(CImage *image);
    void CutH(int hb, int he);
    void CutW(int wb, int we);
    int magnitude(Ipp32fc *pInData, float *pOutData, int N, int M);
    int Roll(int n);
    float AutoBrigtnessFull(float *rawimg);
    float AutoBrigtness(float *rawimg);
    int RangeBrightness(float x0, float dx);

    void PrepareForAF(void);
    void SetImageWidth(int width);
    void SetImageHeight(int height);
    void SetImageType(eImageType type);
    static int GetPixelSize(eImageType type);
    int GetWidth(void);
    int GetHeight(void);
    int CopyImageTop(CImage *srcimg);

    int Norm();

    pybind11::object GetData(int begin, int cnt);
    pybind11::object GetImage(void);

    int ReadMetadata(char *fileName);

    void AddLabel(uint8_t *img, char *label);

    void VerticalMirror();
    void Rotate(float angle);

    void PlotSignal(char *fileName);
};

#endif /* CIMAGE_H_ */
