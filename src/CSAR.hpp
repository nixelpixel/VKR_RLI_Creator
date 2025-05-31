/*
 * CSAR.h
 *
 *  Created on: 3 янв. 2017 г.
 *      Author: User
 */

#ifndef CSAR_H_
#define CSAR_H_

#include "CDSP.hpp"
#include "CDSPCUDA.hpp"
#include "CDSPFFTW.hpp"
#include "CDSPIPP.hpp"
#include "CDSPNEON.hpp"
#include "CLog.hpp"
#include "CRadar.hpp"
#include "ipp/ipp.h"
#include <math.h>
#include <stdio.h>
#include <string.h>
// #include "CSim.hpp"
#include "CImage.hpp"
#include "base.hpp"
#include "helpers.hpp"
// #include "AltMath.hpp"
// #include "CBackProjection.hpp"

/**
 \brief Класс, реализующий обработку радиоголограммы и синтез радиолокационноых
 изображений.
 */

class CSAR {
  public:
    CSAR();
    CSAR(CLog *log, CRadar *radar, DSPType dsp_type);
    virtual ~CSAR();

    /**
     \brief Структура, описывающая параметры радиолокационного изображения
     */

    enum WindowType { HAMMING, BLACKMAN, NUTTALL, NONE };

    void printDSPInfo(void); // Вывод информации о системе ЦОС
    int removeDC(CImage *image); // Удаление постоянной составляющей из сигнала
    int rangeCompress(CImage *image);   // Сжатие по дальности
    int transpose(CImage *image);       // Транспонирование
    int azimuthCompress(CImage *image); // Сжатие по азимуту
    int
    compensateShift(CImage *image); // Коррекция смещения отметок по дальности
                                    // вследствие дробного отношения количества
                                    // отсчетов АЦП за период модуляции
    int rmc(CImage *image, float V); // Коррекция миграции
    int fsrmc(CImage *image,
              float V); // Коррекция миграции методом частотного масштабирования
    int phocus(CImage *image, float V); // Фокусировка изображения
    float enthropy(CImage *image); // Вычисление энтропии
    int hamming(CImage *image);    // Оконное взвешивание
    int windowWeighting(CImage *image, uint64_t window_addr);
    int interpolate(CImage *image, int kr); // Частотная интерполяция
    int hilbert(CImage *image); // Преобразование Гильберта

    int dspProgress();

    int acc_dcm(CImage *image, float Vmax,
                int Np); // Формирование двумерной корреляционной матрицы

    float fr2range(
        float fr); // Пересчет значения разностной частоты к наклонной дальности
    float range2fr(float r);
    float index2range(
        int ind); // Пересчет значения индекса изображения к наклонной дальности
    inline int range2index(
        float range); // Пересчет наклонной дальности к индексу изображения
    inline float
    get_range(float x, float y, float h, float V,
              float t); // Вычисление наклонной дальности по координатам

    int downrange(CImage *image, float decim, float rb, float re,
                  float H);                                 // RDA
    int angleCorrection(CImage *image, float x0, float lx); // RDA

    uint64_t createWindow(int size, WindowType w_type);
    void freeWindow(uint64_t window_addr);

    int resize(CImage *image, int newW, int newH);

    float driftAngle(CImage *image, float V, char *FnPlot = nullptr);

  private:
    CLog *Log;
    CRadar *Radar;
    CDSP *DSP;
    DSPType dspType;

    friend class CBackProjection;
    friend class CStripStream;
};

inline int CSAR::range2index(float range) {
    float Rmax = fr2range(Radar->Fs / 2);
    return round(range * (float)(Radar->Ns2 / 2) / Rmax);
}

inline float CSAR::get_range(float x, float y, float h, float V, float t) {
    float c1 = y + V * t;
    return sqrt(x * x + c1 * c1 + h * h);
}

#endif /* CSAR_H_ */
