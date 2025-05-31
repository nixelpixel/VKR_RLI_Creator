/*
 * CRadar.h
 *
 *  Created on: 31 мая 2017 г.
 *      Author: Рязанцев Леонид
 */

#ifndef CRADAR_H_
#define CRADAR_H_

/**
 \brief Класс, описывающий параметры аппаратуры РЛС и формируемых аппаратурой
 сигналов. Предоставляет данные и реализует функции для работы с ними.
 */

#include "CImage.hpp"
#include "CLog.hpp"
#include <math.h>
#include <pybind11/numpy.h>
#include <stdlib.h>
#include <string.h>

#define SPEED_OF_LIGHT 300000000.0f
#define HEADER_KEY 0xE3F0D0E3

// using namespace pybind11::python;
// using namespace pybind11::literals;

class CRadar {
  public:
    /// Набор возможных значений, определяющих количество бит на один отсчет
    enum SampleSize {
        SAMPLE_SIZE_UNDEFINED = 0, ///< Значение не определено
        SAMPLE_SIZE_1BIT = 1,      ///< 1 бит на отсчет
        SAMPLE_SIZE_2BIT = 2,      ///< 2 бита на отсчет
        SAMPLE_SIZE_4BIT = 4,      ///< 4 бита на отсчет
        SAMPLE_SIZE_8BIT = 8,      ///< 8 бит на отсчет
        SAMPLE_SIZE_16BIT = 16,    ///< 16 бит на отсчет
    };

    /// Набор возможных значений, определяющих закон ленейной частотной
    /// модуляции зондирующего сигнала
    enum RampType {
        RAMP_UNDEFINED = 0, ///< Закон не определен
        RAMP_SAW = 1,       ///< Несимметричная модуляция
        RAMP_TRI = 2, ///< Симметричная (треугольная) модуляция
    };

  public:
    float F0; ///< Начальная частота ЛЧМ зондирующего сигнала АЦП, Гц
    float DF; ///< Девиация ЛЧМ зондирующего сигнала, Гц
    int Fs; ///< Частота дискретизации АЦП, Гц
    int Ns; ///< Количество отсчетов за период зондирования
    int Ns2; ///< Количество отсчетов за период зондирования (дополненное до
             ///< степени двойки)
    int Nsync; ///< Количество отсчетов в синхропаузе
    int Nppb; ///< Количество отсчетов переходных процессов в начале периода
    int Nppe; ///< Количество отсчетов переходных процессов в конце периода
    float
        Mp; ///< Коэффициент коррекции смещения отметок по дальности вследствие
            ///< дробного отношения количества отсчетов АЦП за период модуляции
    SampleSize Sample_Size; ///< Количество бит, отводимых на один отсчет
    RampType
        Ramp_Type; ///< Тип модуляции сигнала (симметричный, несимметричный)
    float Tm;
    int NSeek; ///  Текущая позиция указателя для блочного считывания

    struct Hologram {
        uint8_t *pHologram;
        int BlockSize;
        int Np2;
        int Np;
        int pImageDataSize;
    } HologramData;

    struct header_t {
        uint32_t key;
        uint32_t size;
        uint32_t Ns;
        uint32_t Fs;
        uint16_t Nsync;
        uint16_t Nppb;
        uint16_t Nppe;
        float Tm;
        float F0;
        float DF;
        float Mp;
        uint8_t RampType;
        uint8_t SampleSize;
    } __attribute__((packed));

  protected:
    CLog *Log;

  public:
    CRadar();
    CRadar(CLog *log);
    CRadar(CLog *log, const char *address, const char *port);
    virtual ~CRadar();

    virtual int ReadDataFromFile(char *filename, CImage *image, int Np,
                                 float tshift, int syncpos, int kr, int comp2,
                                 int rmdist) {
        return 0;
    };
    virtual int ReadBlockFromFile(char *filename, CImage *image, int Np,
                                  int kr) {
        return 0;
    };
    virtual int ReadBlockSetPos(int syncpos, float tshift) { return 0; };
    virtual int ReadBlockResetPos(void) { return 0; };

    virtual int ReadDataFromSocket(CImage *image) { return 0; };
    virtual int FindSyncPos(CImage *image, int *sync_begin, int *sync_width) {
        return 0;
    };
    virtual pybind11::tuple FindSyncPosPy(CImage *image) {
        return pybind11::make_tuple(0, 0);
    };
    virtual float GetFileDuration(char *filename) { return 0; };
    virtual int ReadBlock(CImage *image, int block) { return 0; };

    virtual int RemoveTransients(CImage *image) { return 0; };

    void SetParams(float f0, float df, int fs, int ns, int nsync, int nppb,
                   int nppe, float mp, float tm, SampleSize sample_size,
                   RampType ramp_type);
    int SetParamsFromFile(char *path);
    float GetLambda(void);
    float GetTp(void);
    float GetTm(void);
    float GetMu(void);
    int Time2Np(float Ts);
};

#endif /* CRADAR_H_ */
